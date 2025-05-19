import argparse
import logging
import pathlib
import re
import sys
from hashlib import sha256
from importlib.metadata import version
from pickle import dumps

import pgzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from needletail import (
    parse_fastx_file,
)
from rbloom import Bloom
from rich.logging import RichHandler
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
# Create a logger
_logger = logging.getLogger(__name__.split(".")[0])


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
        "--raw-names",
        action="store_true",
        help="Retain the full read name, including the read index information",
        default=False,
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
        "--max-items",
        type=int,
        default=1000000000,
        help="Maximum number of items in the bloom filter",
    )
    parser.add_argument(
        "--fpr",
        type=float,
        default=0.01,
        help="False positive rate for the bloom filter",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="exact",
        choices=["exact", "bloom"],
        help="Method to use for filtering",
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
        If the path is a directory, it will return all files in the directory that match the pattern '*.fq*' or '*.fastq*'.
        If the path is a file, it will return the file itself.
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
        Create a dictionary of paired files from the undetermined path.
        The keys are the lane numbers and the values are lists of paired files.
        """
        # Create the list of all target files
        path = [
            x for x in explode_path(path, recursive=False) if "purged" not in x.name
        ]
        patterns = sorted(
            list(
                {tuple([x.parent, re.sub(r"_R[12]_", "_R[12]_", x.name)]) for x in path}
            )
        )
        path = {
            pattern.split("_")[2]: tuple(pathlib.Path(parent).glob(f"{pattern}"))
            for parent, pattern in patterns
        }
        return {
            key: values if len(values) == 2 else tuple([values[0], None])
            for key, values in path.items()
        }

    args.undetermined_path = create_paired_dict(args.undetermined_path)

    _logger.info("Undetermined files:")
    for lane, targets in args.undetermined_path.items():
        _logger.info(f"    {lane}: '{targets}'")

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
        _logger.info(f"    {lane}: '{targets}'")

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
        _logger.info(f"    {key}: '{value}'")

    return args


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
        Hash function to be used for the bloom filter.
        It uses SHA256 to hash the object and returns the first 16 bytes as an integer.
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
        for record in parse_fastx_file(bloom_source):
            bf.add(record.name)
        # # Alternatively, it is possible to retrieve all names at once, however this might
        # # be memory intensive for large files, though it has not been tested yet
        # fq = pyfastx.Fq(bloom_source, build_index=False)
        # bf.update(fq.keys())
    return bf


def process_target_file(
    target: pathlib.Path,
    output_path: pathlib.Path,
    bloom_filter: Bloom,
    keep_raw_names: bool = False,
) -> None:
    """
    Process the target file and remove reads that are in the bloom filter.
    Args:
        target (pathlib.Path): Path to the target fastq file.
        bloom_filter (Bloom): The bloom filter to use for filtering.
    """

    # # Read the target fastq file and write the purged reads to the output file
    # # Building the index allows access to the read.raw property, which contains the full read name
    # # If the full read name is not needed, it is possible drop the index generation
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
    with (
        pgzip.open(output_file, "wt", thread=4)
        if output_file.suffix == ".gz"
        else open(output_file, "wt") as output_fastq
    ):
        if keep_raw_names:
            # Although the biopython implementation might be slower, it is allows to
            # retain the full read name, including the read index information
            _logger.info("Reading target fastq file using biopython...")
            with (
                pgzip.open(target, "rt", thread=8)
                if target.suffix == ".gz"
                else open(target, "rt") as input_fastq
            ):
                for title, seq, qual in FastqGeneralIterator(input_fastq):
                    # The read name is the first part of the title
                    read_name = title.split(" ")[0]
                    if read_name not in bloom_filter:
                        output_fastq.write(f"@{title}\n{seq}\n+\n{qual}\n")
                        valid += 1
                    total += 1
                    if total % 1000000 == 0:
                        _logger.info(f"    Processed {total} reads")

        else:
            # If the full read name is not needed, it is possible to use pyfastx
            # to read the fastq file, which is faster than biopython
            _logger.info("Reading target fastq file using pyfastx...")
            for record in parse_fastx_file(str(target)):
                if record.name not in bloom_filter:
                    output_fastq.write(
                        f"@{record.name}\n{record.seq}\n+\n{record.qual}\n"
                    )
                    valid += 1
                total += 1
                if total % 1000000 == 0:
                    _logger.info(f"    Processed {total} reads")

    logging.info(f"    Removed {total - valid} reads")
    logging.info(f"    Kept {valid} reads")
    logging.info(f"    Purged fastq file saved to '{output_file}'")


def main() -> None:
    """Main function to run the script."""

    args = parse_args()

    # Set the logging level based on the command line argument
    _logger.setLevel(args.log_level)

    # Validate the command line arguments
    args = validate_args(args)

    if args.method == "bloom":
        _logger.info("Using bloom filter method")
        _logger.info("Building bloom filter...")
        bf = build_bloom_filter(args.undetermined_path, args.max_items, args.fpr)
        # Log the bloom filter parameters
        _logger.debug(f"Bloom filter size: {bf.size_in_bits} bits")
        # _logger.debug(f"Hash functions: {bf.hash_func}")
        _logger.debug(f"Number of items: {bf.approx_items:.1f}")

        # Process the target files and remove reads that are in the bloom filter
        for target in args.assigned_path:
            _logger.info(f"Processing target file '{target}'")
            process_target_file(target, args.output_path, bf, args.raw_names)
    else:
        _logger.info("Using exact method")
        _logger.info("Building python set from undetermined reads...")
        undetermined_set = set()
        output_file = None
        for source_file in args.undetermined_path:
            _logger.debug(f"Adding '{source_file.name}' to set")
            for i, record in enumerate(parse_fastx_file(str(source_file))):
                undetermined_set.add(record.name.split(" ")[0])
                if i % 10000000 == 0:
                    _logger.info(f"    Processed {i} reads")
            if source_file.suffix == ".gz":
                fq_suffix = source_file.with_suffix("").suffix
                output_file = source_file.with_suffix("").with_suffix(
                    f".purged{fq_suffix}.gz"
                )
            else:
                fq_suffix = source_file.suffix
                output_file = source_file.with_suffix(f".purged{fq_suffix}")

            # If provided, use the output path to save the purged fastq file
            if args.output_path:
                output_file = args.output_path.joinpath(output_file.name)

        undetermined_count = len(undetermined_set)
        _logger.info(f"Number of undetermined reads: {undetermined_count}")
        _logger.info("Purging already assigned reads...")
        for target in args.assigned_path:
            _logger.info(f"Processing target file '{target}'")
            already_assigned = set()
            subset = set()
            for i, record in enumerate(parse_fastx_file(str(target))):
                subset.add(record.name.split(" ")[0])
                if i % 100000 == 0:
                    _logger.info(f"    Processed {i} reads")
                    subset.intersection_update(undetermined_set)
                    _logger.info(f"Number of already assigned reads: {len(subset)}")
                    already_assigned.update(subset)
                    subset.clear()
            _logger.info(f"    Processed {i} reads")
            subset.intersection_update(undetermined_set)
            _logger.info(f"Number of already assigned reads: {len(subset)}")
            already_assigned.update(subset)
            _logger.info(f"Found {len(already_assigned)} already assigned reads")
            undetermined_set.difference_update(already_assigned)
        _logger.info(
            f"Found {undetermined_count - len(undetermined_set)} duplicates in the undetermined file."
        )

        if undetermined_count > len(undetermined_set):
            if output_file:
                _logger.info(f"Writing purged reads to '{output_file}'...")
                with (
                    pgzip.open(output_file, "wt", thread=4)
                    if output_file.suffix == ".gz"
                    else open(output_file, "wt") as output_fastq
                ):
                    for record in parse_fastx_file(str(args.undetermined_path[0])):
                        if record.name.split(" ")[0] in undetermined_set:
                            output_fastq.write(
                                f"@{record.name}\n{record.seq}\n+\n{record.qual}\n"
                            )
            else:
                _logger.warning(
                    "Output file not specified! Purged reads will not be saved!"
                )
        else:
            _logger.info("Found no duplicates in the undetermined reads. Exiting...")

    _logger.info("Done!")


if __name__ == "__main__":
    main()
