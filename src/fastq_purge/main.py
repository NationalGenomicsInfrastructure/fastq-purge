import argparse
import gzip
import logging
import pathlib
from importlib.metadata import version
import pyfastx
import re
import sys
from hashlib import sha256
from pickle import dumps

from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
        prog="%(prog)s",
        usage="%(prog)s [options]",
        description="Package to purge fastq files from unwanted reads",
        formatter_class=ArgumentDefaultsRichHelpFormatter,
    )
    parser.add_argument(
        "--target-path",
        type=pathlib.Path,
        help= ("Path to the target fastq file(s) to purge. It can be gzipped or not. ",
             "It can be a single file or a directory. If a directory is provided, ",
             "all files in the directory that match the pattern '*.fq*' or '*.fastq*' ",
             "will be used as target files. The search is non-recursive."),
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
        "--bloom-sources",
        type=pathlib.Path,
        default=None,
        help=("Path to the bloom filter sources file(s). They can be gzipped or not. ",
              "It can be a single file or multiple files separated by spaces, or a directory. ",
              "If a directory is provided, all files in the directory that match the pattern ",
              "'*.fq*' or '*.fastq*' will be used as bloom filter sources. The search is recursive."),
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
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('fastq_purge')}",
        help="Show version and exit",
    )
    return parser.parse_args()


def validate_args(
    args: argparse.Namespace
) -> argparse.Namespace:
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
        patterns = ["**/*.fq*", "**/*.fastq*"] if recursive else ["*.fq", "*.fastq"]
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

    # Create the list of all target files
    args.target_path = [ x for x in explode_path(args.target_path, recursive=False) if "purged" not in x.name ]
    _logger.info("Target files:")
    for i, target in enumerate(args.target_path):
        _logger.info(f"    {i + 1:02d}: '{target}'")

    # Check whether the output path exists, if not, create it
    if args.output_path:
        args.output_path = pathlib.Path(args.output_path)
        if not args.output_path.is_dir():
            _logger.debug(
                f"The output path does not exist! It will be created..."
            )
            args.output_path.mkdir(parents=True, exist_ok=True)

    # Create the list of all bloom sources
    args.bloom_sources = list(set([ bs for bloom_source in args.bloom_sources for bs in explode_path(bloom_source, recursive=True) ]))

    # Valiedate the bloom sources, removing unnecessary files (e.g. retain only one of the paired reads)
    sources_set = set()
    clean_bloom_sources = []
    for bloom_source in sorted(args.bloom_sources):
        source_basename = re.sub(r"_R[12]_", "_", bloom_source.name)
        if source_basename in sources_set:
            _logger.warning(
                f"Ignoring '{bloom_source}' as it is the paired read another source file"
            )
            continue
        sources_set.add(source_basename)
        clean_bloom_sources.append(bloom_source)
    args.bloom_sources = clean_bloom_sources
    _logger.info("Bloom filter sources:")
    for i, bloom_source in enumerate(args.bloom_sources):
        _logger.info(f"    {i + 1:02d}: '{bloom_source}'")
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
        for name,_,_ in pyfastx.Fastq(str(bloom_source), build_index=False):
            bf.add(name)
        # # Alternatively, it is possible to retrieve all names at once, however this might
        # # be memory intensive for large files, though it has not been tested yet
        # fq = pyfastx.Fq(bloom_source, build_index=False)
        # bf.update(fq.keys())
    return bf

def process_target_file(
    target: pathlib.Path, output_path: pathlib.Path, bloom_filter: Bloom, keep_raw_names: bool = False) -> None:
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
        output_file = target.with_suffix("").with_suffix(
            f".purged{fq_suffix}.gz"
        )
    else:
        fq_suffix = target.suffix
        output_file = target.with_suffix(f".purged{fq_suffix}")

    # If provided, use the output path to save the purged fastq file
    if output_path:
        output_file = output_path.joinpath(output_file.name)

    _logger.info(f"Output file: '{output_file}'")
    with (
        gzip.open(output_file, "wt")
        if output_file.suffix == ".gz"
        else open(output_file, "wt") as output_fastq
    ):
        if keep_raw_names:
            # Although the biopython implementation might be slower, it is allows to
            # retain the full read name, including the read index information
            _logger.info(f"Reading target fastq file using biopython...")
            with gzip.open(target, "rt") if target.suffix == ".gz" else open(target, "rt") as input_fastq:
                for title, seq, qual in FastqGeneralIterator(input_fastq):
                    # The read name is the first part of the title
                    read_name = title.split(" ")[0]
                    if not read_name in bloom_filter:
                        output_fastq.write(f"@{title}\n{seq}\n+\n{qual}\n")
                        valid += 1
                    total += 1
                    if total % 1000000 == 0:
                        _logger.info(f"    Processed {total} reads")
                        

        else:
            # If the full read name is not needed, it is possible to use pyfastx
            # to read the fastq file, which is faster than biopython
            _logger.info(f"Reading target fastq file using pyfastx...")
            for name, seq, qual  in pyfastx.Fastq(str(target), build_index=False):
                if not name in bloom_filter:
                    output_fastq.write(f"@{name}\n{seq}\n+\n{qual}\n")
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

    _logger.info(f"Building bloom filter...")
    bf = build_bloom_filter(
        args.bloom_sources, args.max_items, args.fpr
    )
    # Log the bloom filter parameters
    _logger.debug(f"Bloom filter size: {bf.size_in_bits} bits")
    # _logger.debug(f"Hash functions: {bf.hash_func}")
    _logger.debug(f"Number of items: {bf.approx_items:.1f}")

    # Process the target files and remove reads that are in the bloom filter
    for target in args.target_path:
        _logger.info(f"Processing target file '{target}'")
        process_target_file(target, args.output_path, bf, args.raw_names)
    _logger.info("Done!")
    

if __name__ == "__main__":
    main()
