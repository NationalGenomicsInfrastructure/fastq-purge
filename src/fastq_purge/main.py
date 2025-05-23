import argparse
import logging
import os
import pathlib
import re
import sys
from hashlib import sha256
from importlib.metadata import version
from itertools import repeat
from pickle import dumps

import dnaio
import multiprocess as mp
import psutil
from loky import get_reusable_executor
from multiprocess import Process, Semaphore
from multiprocess.managers import BaseManager
from multiprocess.pool import Pool
from rbloom import Bloom
from rich.logging import RichHandler
from rich.progress import track
from rich_argparse import ArgumentDefaultsRichHelpFormatter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="%H:%M:%S",
    force=True,
    handlers=[
        RichHandler(rich_tracebacks=True, tracebacks_show_locals=True, markup=True)
    ],
)
# Create a local logger for the module
_logger = logging.getLogger(__name__.split(".")[0])

# Define a global variable to store the undetermined set
# This is needed to avoid passing the object to the child process
undetermined_set = set()


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog=f"{sys.argv[0].split('/')[-1]}",
        usage="%(prog)s [options]",
        description="Package to purge fastq files from unwanted reads",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "--undetermined-path",
        type=pathlib.Path,
        help="""Path to the target fastq file(s) to purge. It can be gzipped or not.
        It can be a single file or a directory. If a directory is provided,
        all files in the directory that match the pattern '*.fq*' or '*.fastq*'
        will be used as target files. The search is non-recursive.""",
        required=True,
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=None,
        help="Path to the output fastq file",
        required=False,
    )
    parser.add_argument(
        "--assigned-path",
        type=pathlib.Path,
        default=None,
        help="""Path to the bloom filter sources file(s). They can be gzipped or not.
        It can be a single file or multiple files separated by spaces, or a directory.
        If a directory is provided, all files in the directory that match the pattern
        '*.fq*' or '*.fastq*' will be used as bloom filter sources. The search is recursive.""",
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="exact",
        choices=["exact", "approx"],
        help="""Method to use for filtering. The 'exact' method uses a python set to store 
        the undetermined reads. The 'approx' method uses a bloom filter to store the undetermined
        reads. The 'exact' method is more resource intensive, but it is more accurate. The
        'approx' method is less resource intensive, but it has a false positive rate associated with it.
        """,
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=1000000000,
        help="""The estimated number of items in the bloom filter. It is used together with the
        false positive rate to calculate the size of the bloom filter. The default is 1 billion.
        This is ignored if the method is 'exact'.""",
    )
    parser.add_argument(
        "--fpr",
        type=float,
        default=0.0001,
        help="False positive rate for the bloom filter. The default is 0.0001.",
    )
    parser.add_argument(
        "--threading-method",
        type=str,
        default="mp_pool",
        choices=["loky", "mp_pool", "mp_manager"],
        help="""Threading library and method to use for processing (used for testing purposes).
        'loky' uses the loky library, which is a robust, cross-platform and cross-version
        implementation of the ProcessPoolExecutor class of concurrent.futures. 'mp_pool' uses
        the multiprocess library, which is a fork of the multiprocessing library with enhanced
        serialization using the dill library. 'mp_manager' uses the multiprocess library with a custom
        manager to share data between processes. The 'mp_pool' and 'mp_manager' methods are
        equivalent, but the 'mp_pool' method less resource intensive.""",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for processing",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('fastq-purge')}",
        help="Show version and exit",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate command line arguments.

    Args:
        args (argparse.Namespace): Parsed command line arguments.
    Returns:
        argparse.Namespace: Validated arguments.
    """

    def explode_path(path: pathlib.Path, recursive: bool = False) -> list[pathlib.Path]:
        """
        Explode a path into a list of files or directories.
        If the path is a directory, it will return all files in the directory that match
        the pattern '*.fq*' or '*.fastq*'. If the path is a file, it will return the file itself.

        Args:
            path (pathlib.Path): Path to the target file or directory.
            recursive (bool): Whether to search recursively in the directory.
        Returns:
            list[pathlib.Path]: List of files or directories.
        """
        patterns = ["**/*.fq*", "**/*.fastq*"] if recursive else ["*.fq*", "*.fastq*"]
        if path.is_dir():
            files = []
            for pattern in patterns:
                files += list(pathlib.Path(path).glob(pattern))
            return files
        elif path.is_file():
            return [path]
        else:
            _logger.error(f"Path '{path}' is neither a file nor a directory")
            exit(1)

    def create_paired_dict(path: list[pathlib.Path]) -> dict[str, list[pathlib.Path]]:
        """
        Create a dictionary of paired files from the list of paths.
        The keys are the lane numbers and the values are tuples of paired files.

        Args:
            path (list[pathlib.Path]): List of paths to the target files.
        Returns:
            dict[str, list[pathlib.Path]]: Dictionary of paired files.
        """
        # Create the list of all target files
        path = [
            x for x in explode_path(path, recursive=False) if "purged" not in x.name
        ]
        patterns = sorted(
            list(
                {
                    tuple(([x.parent, re.sub(r"_R[12]_", "_R[12]_", x.name)]))
                    for x in path
                }
            )
        )
        path = {
            pattern.split("_")[2]: tuple(pathlib.Path(parent).glob(f"{pattern}"))
            for parent, pattern in patterns
        }
        return {
            key: tuple(sorted(list(values)))
            if len(values) == 2
            else tuple([values[0], None])
            for key, values in path.items()
        }

    args.undetermined_path = create_paired_dict(args.undetermined_path)

    _logger.info("Undetermined files:")
    for lane, targets in args.undetermined_path.items():
        for target in targets:
            _logger.info(f"    {lane}: '{target}'")

    # Check whether the output path exists, if not, create it
    args.output_path = (
        pathlib.Path(args.output_path)
        if args.output_path
        else pathlib.Path(
            args.undetermined_path[list(args.undetermined_path.keys())[0]][0].parent
        )
    )
    if not args.output_path.is_dir():
        _logger.debug("The output path does not exist! It will be created...")
        args.output_path.mkdir(parents=True, exist_ok=True)
    output_dict = dict()
    for key, values in args.undetermined_path.items():
        output_list = []
        for input_file in values:
            if input_file is None:
                output_list.append(None)
            else:
                if input_file.suffix == ".gz":
                    fq_suffix = input_file.with_suffix("").suffix
                    output_list.append(
                        args.output_path.joinpath(
                            input_file.with_suffix("")
                            .with_suffix(f".purged{fq_suffix}.gz")
                            .name
                        )
                    )
                else:
                    fq_suffix = input_file.suffix
                    output_list.append(
                        args.output_path.joinpath(
                            input_file.with_suffix(f".purged{fq_suffix}").name
                        )
                    )
            output_dict[key] = tuple(output_list)
    args.output_path = output_dict

    _logger.info("Output files:")
    for lane, targets in args.output_path.items():
        for target in targets:
            _logger.info(f"    {lane}: '{target}'")

    # Create the list of all bloom sources
    args.assigned_path = list(
        set(
            [
                bs
                for assigned_file in args.assigned_path
                for bs in explode_path(assigned_file, recursive=True)
            ]
        )
    )

    # Valiedate the bloom sources, removing unnecessary files (e.g. retain only one of the paired reads)
    sources_set = set()
    clean_assigned_files = []
    for assigned_file in sorted(args.assigned_path):
        source_basename = re.sub(r"_[IR][0-9]_", "_", assigned_file.name)
        if source_basename in sources_set:
            _logger.warning(
                f"Ignoring '{assigned_file}' as it is either a paired or index file"
            )
            continue
        sources_set.add(source_basename)
        clean_assigned_files.append(assigned_file)
    args.assigned_path = clean_assigned_files

    assigned_dict = dict()
    for assigned_file in args.assigned_path:
        try:
            key = re.search(r"_L[0-9]{3}_", assigned_file.name).group().strip("_")
        except AttributeError:
            _logger.warning(
                f"File '{assigned_file}' does not match the expected pattern 'L[0-9]{3}'"
            )
            continue
        if key in assigned_dict:
            assigned_dict[key].append(assigned_file)
        else:
            assigned_dict[key] = [assigned_file]
    args.assigned_path = assigned_dict
    _logger.info("Assigned files:")
    for key, value in args.assigned_path.items():
        for val in value:
            _logger.info(f"    {lane}: '{val}'")

    return args


def memory_usage(logger, pid: int = None) -> None:
    """
    Print the memory usage of the current process.

    Args:
        logger (logging.Logger): Logger to use for logging the memory usage.
        pid (int, optional): Process ID to check memory usage for. If None, uses the current process.
    """
    vmem = psutil.virtual_memory()
    smem = psutil.swap_memory()
    prefix = f"PID {pid} - " if pid else ""
    logger.debug(f"{prefix}Memory used: {vmem.used / 1024**3:.2f} GB ({vmem.percent}%)")
    logger.debug(f"{prefix}Swap used: {smem.used / 1024**3:.2f} GB ({smem.percent}%)")


def build_bloom_filter(
    bloom_sources: list[pathlib.Path], max_items: int, fpr: float
) -> Bloom:
    """
    Build a bloom filter from the given sources.

    Args:
        bloom_sources (list[pathlib.Path]): List of paths to the bloom filter source files.
        max_items (int): Maximum number of items in the bloom filter.
        fpr (float): False positive rate for the bloom filter.

    Returns:
        Bloom: The constructed bloom filter.
    """

    def hash_func(obj):
        """
        Hash function to use for the bloom filter. It uses SHA256 to hash the object and
        returns the first 16 bytes as an integer.
        This is a simple hash function, but it can be replaced with a more complex one if needed.

        Args:
            obj: The object to hash.
        Returns:
            int: The hash value. Hash function to be used for the bloom filter.
        """
        # Use the first 16 bytes of the SHA256 hash as the hash value
        # This is a simple hash function, but it can be replaced with a more complex one if needed
        # The hash value is converted to an integer
        # using the system byte order and signed=True to allow negative values
        h = sha256(dumps(obj)).digest()
        return int.from_bytes(h[:16], sys.byteorder, signed=True)

    # Create a bloom filter with the given parameters
    bf = Bloom(max_items, fpr, hash_func)
    for bloom_source in bloom_sources:
        _logger.debug(f"Adding '{bloom_source.name}' to bloom filter")
        with dnaio.open(bloom_source) as reader:
            for record in track(
                reader,
                description=f"Parsing '{bloom_source.name}'...",
                total=None,
                transient=True,
            ):
                bf.add(record.name)
    return bf


def process_target_file(
    target: pathlib.Path,
    output_path: pathlib.Path,
    bloom_filter: Bloom,
    threads: int = 1,
) -> None:
    """
    Process the target file and remove reads that are in the bloom filter.

    Args:
        target (pathlib.Path): Path to the target fastq file.
        bloom_filter (Bloom): The bloom filter to use for filtering.
    """

    # Read the target fastq file and write the purged reads to the output file
    # Building the index allows access to the read.raw property, which contains the full read name
    # If the full read name is not needed, it is possible drop the index generation
    total = 0
    valid = 0
    if target.suffix == ".gz":
        fq_suffix = target.with_suffix("").suffix
        output_file = target.with_suffix("").with_suffix(f".purged{fq_suffix}.gz")
    else:
        fq_suffix = target.suffix
        output_file = target.with_suffix(f".purged{fq_suffix}")

    # If provided, use the output path to save the purged fastq file
    if output_path:
        output_file = output_path.joinpath(output_file.name)

    _logger.info(f"Output file: '{output_file}'")
    with dnaio.open(
        output_file, mode="w", compression_level=7, open_threads=(threads + 1) // 2
    ) as writer:
        _logger.info("Reading target fastq file...")
        with dnaio.open(target, open_threads=(threads + 1) // 2) as reader:
            for record in reader:
                if record.name not in bloom_filter:
                    writer.write(record)
                    valid += 1
                total += 1
                if total % 1000000 == 0:
                    _logger.info(f"    Processed {total} reads")

    logging.info(f"    Removed {total - valid} reads")
    logging.info(f"    Kept {valid} reads")
    logging.info(f"    Purged fastq file saved to '{output_file}'")


def build_undetermined_set(path: pathlib.Path) -> set:
    """
    Build a set of undetermined reads from the given path.

    Args:
        path (pathlib.Path): Path to the undetermined fastq file.
    Returns:
        set: Set of undetermined reads.
    """
    _logger.debug(f"Adding '{path.name}' to set")
    with dnaio.open(path) as reader:
        for record in track(
            reader,
            description=f"Parsing '{path.name}'...",
            total=None,
            transient=True,
        ):
            undetermined_set.add(record.name.split(" ")[0])
    return undetermined_set


def get_logger(name: str, log_level: str = "INFO") -> logging.Logger:
    """
    Get a logger with the given name and log level. Used in child processes.

    Args:
        name (str): Name of the logger.
        log_level (str): Log level to set for the logger.
    Returns:
        logging.Logger: The logger with the given name and log level.
    """
    # Get the logger for the given name
    logger = logging.getLogger(name)
    # Set the logging level
    logger.setLevel(log_level)
    return logger


def process_target_set(
    path: pathlib.Path,
    log_level: str = "INFO",
    buffer_size: int = 1000000,
) -> set:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        path (pathlib.Path): Path to the target fastq file.
        undetermined_set (set): Set of undetermined reads.
        buffer_size (int): Size of the buffer to use for processing.
    Returns:
        tuple: Tuple containing the process ID, target file name, and set of already assigned reads.
    """
    # Use the global variable to access the undetermined set
    undetermined_set = os.undetermined_set

    # Create a logger for the process
    logger = get_logger("loky", log_level)
    pid = os.getpid()

    logger.info(f"Processing target file '{path.name}' with pid {pid}")
    already_assigned = set()
    tmp_set = set()
    memory_usage(logger, pid)
    with dnaio.open(path) as reader:
        for i, record in enumerate(reader):
            tmp_set.add(record.name.split(" ")[0])
            if i % buffer_size == 0 and i > 0:
                logger.debug(f"    {pid}: Processed {i} reads")
                tmp_set.intersection_update(undetermined_set)
                already_assigned.update(tmp_set)
                tmp_set.clear()
        logger.debug(f"    {pid}: Processed {i} reads")
        tmp_set.intersection_update(undetermined_set)
        already_assigned.update(tmp_set)
    logger.info(
        f"Process {pid} ('{path.name}') finished yielding {len(already_assigned)} assigned reads"
    )
    return (pid, path.name, already_assigned)


def loky_process_target_set(
    assigned_dict: dict,
    threads: int = 1,
    log_level: str = "INFO",
) -> set:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        assigned_dict (dict): Dictionary of target files to process.
        threads (int): Number of threads to use for processing.
        log_level (str): Log level to set for the logger.
    Returns:
        set: Set of already assigned reads.
    """

    # Hackish trick to pass a global variable to the child process
    # by mutating a global variable from a module (such as the os module)
    def set_value(value):
        """
        Set the value of the global variable in the child process.

        Args:
            value (set): The value to set for the global variable.
        """
        _logger.info(f"Setting value in child process {os.getpid()}")
        os.undetermined_set = value

    _logger.info("Initializing loky executor...")
    executor = get_reusable_executor(
        max_workers=threads,
        timeout=30,
        kill_workers=True,
        initializer=set_value,
        initargs=(undetermined_set,),
    )
    _logger.info("Purging already assigned reads...")
    already_assigned = set()
    for as_lane, targets in assigned_dict.items():
        results = executor.map(
            process_target_set,
            targets,
            repeat(log_level),
        )
        already_assigned.update(*[res[2] for res in results])
    return already_assigned


class CustomManager(BaseManager):
    """Custom manager to share data between processes."""

    # nothing
    pass


class MultiprocessingCustom:
    """Custom class to be used with the multiprocessing manager."""

    # constructor
    def __init__(self, data):
        # store the data in the instance
        self.undetermined = data
        self.buffer_size = 1000000
        self.assigned = set()
        self.log = dict()

    def intersect_and_update(self, subset) -> int:
        """
        Intersect the subset with the undetermined set and update the assigned set.

        Args:
            subset (set): The subset to intersect with the undetermined set.
        Returns:
            int: The number of reads that were already assigned.
        """
        subset.intersection_update(self.undetermined)
        count = len(subset)
        self.assigned.update(subset)
        subset.clear()
        return count

    # perform the main task
    def task(self, path):
        """
        Perform the main task of the custom class.
        Args:
            path (pathlib.Path): Path to the target fastq file.
        Returns:
            tuple: Tuple containing the process ID, target file name, and number of already assigned reads.
        """
        tmp_set = set()
        counts = 0
        with dnaio.open(path) as reader:
            for i, record in enumerate(reader):
                tmp_set.add(record.name.split(" ")[0])
                if i % self.buffer_size == 0 and i > 0:
                    counts += self.intersect_and_update(tmp_set)
        counts += self.intersect_and_update(tmp_set)
        _logger.debug(
            f"Process ('{path.name}') finished yielding {counts} already assigned reads"
        )
        self.log[path.name] = counts
        return (path.name, counts)

    # view the log
    def view_log(self):
        """View the log of the custom class."""
        for key, value in self.log.items():
            _logger.info(f"{key}: {value}")

    # get all stored values
    def get_assigned(self):
        """Get all stored values."""
        return self.assigned

    # remove all stored values
    def clear_assigned(self):
        """Remove all stored values."""
        return self.assigned.clear()


def manager_work(shared_custom, path, semaphore, log_level="INFO"):
    """
    Custom function to be executed in a child process.

    Args:
        shared_custom (MultiprocessingCustom): Shared custom class instance.
        path (pathlib.Path): Path to the target fastq file.
        semaphore (Semaphore): Semaphore to limit the number of concurrent processes.
        log_level (str): Log level to set for the logger.
    """
    # acquire the semaphore
    semaphore.acquire()
    # create a logger for the process
    logger = get_logger("multiprocessing", log_level)
    pid = os.getpid()
    logger.debug(f"Process {pid} started: {path.name} with {mp.get_start_method()}")
    # call the function on the shared custom instance
    name, number = shared_custom.task(path)
    # return the result
    logger.info(f"Process {pid} finished: {name} {number}")
    # release the semaphore
    semaphore.release()


def pool_task(path, log_level="INFO") -> tuple:
    """
    Custom function to be executed in a child process.

    Args:
        path (pathlib.Path): Path to the target fastq file.
        log_level (str): Log level to set for the logger.
    Returns:
        tuple: Tuple containing the process ID, target file name, and set of already assigned reads.
    """
    logger = get_logger("multiprocessing", log_level)
    pid = os.getpid()
    logger.debug(f"Process {pid} started: {path.name} with {mp.get_start_method()}")
    global undetermined_set

    memory_usage(logger, pid)

    tmp_set = set()
    assigned_set = set()
    with dnaio.open(path) as reader:
        for i, record in enumerate(reader):
            tmp_set.add(record.name.split(" ")[0])
            if i % 1000000 == 0 and i > 0:
                tmp_set.intersection_update(undetermined_set)
                assigned_set.update(tmp_set)
                tmp_set.clear()
        tmp_set.intersection_update(undetermined_set)
        assigned_set.update(tmp_set)
    _logger.info(
        f"Process {pid} ('{path.name}') finished yielding {len(assigned_set)} already assigned reads"
    )
    return (path.name, assigned_set)


def multiprocessing_process_target_set(
    assigned_dict: dict,
    threads: int = 1,
    log_level: str = "INFO",
    method: str = "pool",  # either "manager" or "pool"
) -> None:
    """
    Process the target file and remove reads that are in the undetermined set.

    Args:
        assigned_dict (dict): Dictionary of target files to process.
        threads (int): Number of threads to use for processing.
        log_level (str): Log level to set for the logger.
        method (str): Method to use for multiprocessing. Either "manager" or "pool".
    Returns:
        set: Set of already assigned reads.
    """
    mp.set_start_method("fork", force=True)

    already_assigned = set()
    memory_usage(_logger)

    if method == "manager":
        _logger.info("Using multiprocessing custom manager...")
        # Method 1: Using a custom manager
        CustomManager.register("MultiprocessingCustom", MultiprocessingCustom)
        with CustomManager() as manager:
            # create a shared custom class instance
            shared_custom = manager.MultiprocessingCustom(undetermined_set)
            memory_usage(_logger)
            semaphore = Semaphore(threads)
            for as_lane, targets in assigned_dict.items():
                _logger.debug(
                    f"Creating {len(targets)} child processes for lane '{as_lane}'"
                )
                processes = [
                    Process(
                        target=manager_work,
                        args=(shared_custom, target, semaphore, log_level),
                    )
                    for target in targets
                ]
                _logger.debug(f"Starting {len(processes)} child processes")
                for process in processes:
                    process.start()
                _logger.debug(f"Waiting for {len(processes)} child processes to finish")
                for process in processes:
                    process.join()
                _logger.debug("Child process finished")
                already_assigned.update(shared_custom.get_assigned())
                shared_custom.clear_assigned()
    else:
        _logger.info("Using multiprocessing pool...")
        # Method 2: Using a process pool
        _logger.info("Purging already assigned reads...")
        for as_lane, targets in assigned_dict.items():
            _logger.debug(f"Creating {len(targets)} child processes for lane {as_lane}")
            # create and configure the process pool
            with Pool(processes=threads) as pool:
                # issue tasks to the process pool
                results = pool.map(pool_task, targets)
            for res in results:
                _logger.info(
                    f"Process finished: {res[0]} with {len(res[1])} already assigned reads"
                )
            for res in results:
                already_assigned.update(res[1])
    return already_assigned


def write_purged_fastq(
    path_undetermined: pathlib.Path,
    path_purged: pathlib.Path,
    assigned_set: set,
    threads: int = 1,
) -> None:
    """
    Write the purged fastq file.

    Args:
        path_undetermined (pathlib.Path): Path to the undetermined fastq file.
        path_purged (pathlib.Path): Path to the purged fastq file.
        assigned_set (set): Set of already assigned reads.
        threads (int): Number of threads to use for processing.
    """
    if path_undetermined[1] is None or path_purged[1] is None:
        _logger.info(f"Writing purged fastq file '{path_purged[0].name}'")
        with dnaio.open(
            path_purged[0],
            mode="w",
            compression_level=7,
            open_threads=(threads + 1) // 2,
        ) as writer:
            with dnaio.open(
                path_undetermined[0], open_threads=(threads + 1) // 2
            ) as reader:
                for record in reader:
                    if record.name.split(" ")[0] not in assigned_set:
                        writer.write(record)
    else:
        _logger.info(
            f"Writing purged fastq files '{path_purged[0].name}' and '{path_purged[1].name}'"
        )
        with dnaio.open(
            path_undetermined[0], path_undetermined[1], open_threads=(threads + 1) // 2
        ) as reader:
            with dnaio.open(
                path_purged[0],
                path_purged[1],
                mode="w",
                compression_level=7,
                open_threads=(threads + 1) // 2,
            ) as writer:
                for r1, r2 in reader:
                    if r1.name.split(" ")[0] not in assigned_set:
                        writer.write(r1, r2)


def main() -> None:
    """Main function to run the script."""

    args = parse_args()

    # Set the logging level based on the command line argument
    _logger.setLevel(args.log_level)

    # Validate the command line arguments
    args = validate_args(args)

    if args.method == "approx":
        _logger.info("Using bloom filter method")
        _logger.info("Building bloom filter...")
        # Create a bloom filter with the given parameters
        bf = build_bloom_filter(args.undetermined_path, args.max_items, args.fpr)
        # Log the bloom filter parameters
        _logger.debug(f"Bloom filter size: {bf.size_in_bits} bits")
        # _logger.debug(f"Hash functions: {bf.hash_func}")
        _logger.debug(f"Number of items: {bf.approx_items:.1f}")

        # Process the target files and remove reads that are in the bloom filter
        for target in args.assigned_path:
            _logger.info(f"Processing target file '{target}'")
            process_target_file(target, args.output_path, bf, args.threads)
    else:
        _logger.info("Using exact method")
        memory_usage(_logger)
        _logger.info("Building python set from undetermined reads...")
        for un_lane, (un_file_1, un_file_2) in args.undetermined_path.items():
            undetermined_set = build_undetermined_set(un_file_1)
            undetermined_count = len(undetermined_set)
            set_size_mb = sys.getsizeof(undetermined_set) / 1024**2
            _logger.info(
                f"Number of undetermined reads: {undetermined_count} ({set_size_mb:.2f} MB)"
            )

            memory_usage(_logger)
            if args.threading_method == "loky":
                # Loky executor
                already_assigned = loky_process_target_set(
                    args.assigned_path,
                    threads=args.threads,
                    log_level=args.log_level,
                )
            else:
                # Multiprocess
                method = "pool" if args.threading_method == "mp_pool" else "manager"
                already_assigned = multiprocessing_process_target_set(
                    args.assigned_path,
                    threads=args.threads,
                    log_level=args.log_level,
                    method=method,
                )

            _logger.info(
                f"Found {len(already_assigned)} duplicates in the undetermined file."
            )
            if already_assigned:
                write_purged_fastq(
                    (un_file_1, un_file_2),
                    args.output_path[un_lane],
                    already_assigned,
                    args.threads,
                )

            else:
                _logger.info(
                    "Found no duplicates in the undetermined reads. Exiting..."
                )

    _logger.info("Done!")


if __name__ == "__main__":
    main()
