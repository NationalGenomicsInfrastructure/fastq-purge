# FastQ-purge

Python package for removing from a set of FASTQ files reads that are present in another set of FASTQ files.

## Motivation

When demultiplexing a sequencing run that contains multiple sequencing designs, it is common to have reads belonging to one design ending up in the _undetermined_ FASTQ files of another design. In such cases, to avoid including duplicates, the _undetermined_ files for one or more lanes are not preserved, and their metrics will not available in GenStat. To resolve this, the _undetermined_ FASTQ files need to be purged of reads that are have been already assigned in either of the samples run on the same lane/flowcell. This package aims at providing a performance-optimized solution to this problem.

## Requirements

The package requires Python 3.11 or higher, and the [uv](https://docs.astral.sh/uv/) package and project manager to be installed. The package is compatible with both Linux and MacOS, but it has not been tested on Windows machines.

## Installation

Although it is recommended to use `uv` to run the package, it is not strictly required. The package can be built and installed using the standard Python packaging tools. To build the package, run the following command:

```bash
uv build
```

This will create a `dist` directory containing the built package. The package can then be installed using pip (or any other package manager that supports wheels). To install the package with pip, run the following command:

```bash
pip install dist/fastq_purge-0.1.0-py3-none-any.whl
```

## Usage

The package can be used as a command line tool, either by running the `fastq-purge` (if installed with pip) or `uv run fastq-purge` (if using uv) command. The minimum required arguments are the paths to the _undetermined_ and _assigned_ FASTQ files, which can be provided as follows:

```bash
fastq-purge --undetermined-path <undetermined_path> --assigned-path <assigned_path>
```

### Main arguments

Both `--undetermined-path` and `--assigned-path` can be provided as a path to a single FASTQ file or a directory containing multiple FASTQ files. When providing a directory, the main difference is that the _undetermined_ FASTQ files will be searched for only in the specified directory, while the _assigned_ FASTQ files will be searched for in the specified directory and all its subdirectories. This is because the _undetermined_ FASTQ files are usually stored in the same directory, while the _assigned_ FASTQ files are usually stored in a subdirectory named after the project and the sample name.

The package will automatically detect the file format (gzipped or not), and the output will generated accordingly. By default, the output will be written to the same directory as the _undetermined_ FASTQ files, with the suffix `.purged` added to the filename (e.g. `<sample_name>.fastq.gz` -> `<sample_name>.purged.fastq.gz`). Alternatively, the output can be written to a different directory by using the `--output-dir` argument. The output directory will be created if it does not exist.

By default, the package will use a single thread to process the files. However, this can be changed by using the `--threads` argument. In additon, there are different libraries and methods that have been implemented for parallelization, and the user can choose which one to use by using the `--threading-method` argument. The available methods are `mp_pool` (default), `mp_manager`, and `loky`. The `mp_pool` method uses the `multiprocessing.Pool` class to create a pool of worker processes, while the `mp_manager` method uses the `multiprocessing.Manager` class to create a manager process that manages a pool of worker processes. The `loky` method uses the `loky` library, which is a robust and efficient library for parallel processing in Python. The `loky` and `mp_manager` methods are recommended for small files (less than 100M reads in each _undetermined_ FASTQ file), while the `mp_pool` method is recommended for larger files.

### Optional arguments

The package also provides several optional arguments that can be used to customize the behavior of the package. These arguments are:

- `--method`: The method to use for purging the reads. The available methods are `exact` (default) and `approx`. The `exact` method uses the `set` data structure to store the read names and check for duplicates, while the `approx` method uses a Bloom filter to store the read names and check for duplicates. The `approx` method is recommended when the _undetermined_ read names cannot fit in memory, while the `exact` method is recommended in all other cases. The `approx` method is also faster than the `exact` method, but it may produce false positives (i.e. reads that are not actually duplicates).
- `--max-items`: The maximum number of items expected to be stored in the Bloom filter. This argument is only used when the `--method` argument is set to `bloom`. The default value is 1e9 (1 billion) items. Together with the `--fpr` argument, this value determines the final size of the Bloom filter.
- `--fpr`: The false positive rate of the Bloom filter. This argument is only used when the `--method` argument is set to `bloom`. The default value is 0.0001 (0.01%). Together with the `--max-items` argument, this value determines the final size of the Bloom filter.
- `--log-level`: The log level to use. The available log levels are `DEBUG`, `INFO`, `WARNING`, `ERROR`, and `CRITICAL`. The default value is `INFO`. The log level can be set to `DEBUG` to enable debug logging, which will provide more detailed information about the progress of the package.
