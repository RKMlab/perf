# PERF
## Introduction
PERF is a Python package developed for fast and accurate identification of microsatellites from DNA sequences. Microsatellites or Simple Sequence Repeats (SSRs) are short tandem repeats of 1-6nt motifs. They are present in all genomes, and have a wide range of uses and functional roles. The existing tools for SSR identification have one or more caveats in terms of speed, comprehensiveness, accuracy, ease-of-use, flexibility and memory usage. PERF was designed to address all these problems.

PERF is a recursive acronym that stands for "PERF is an Exhaustive Repeat Finder". It is compatible with both Python 2 (tested on Python 2.7) and 3 (tested on Python 3.5). Its key features are:
  - Fast run time, despite being a single-threaded application. As an example, identification of all SSRs from the entire human genome takes less than 10 minutes. The speed can be approved 5-fold using PyPy (has a runtime of less than 2min using PyPy 5.8)
  - Linear time and space complexity (O(n))
  - 100% accurate and comprehensive - Does not miss any repeats or does not pick any incorrect ones
  - Easy to use - The only required argument is the input DNA sequence in FASTA format
  - Flexible. While having sensible defaults, most of the parameters are customizable by the user at runtime.
  - Repeat cutoffs can be specified either in terms of the total repeat length or in terms of number of repeating units.
  - TSV and HTML report. The default output is an easily parseable and exportable tab-separated format. Optionally, PERF also generates an interactive HTML report that depicts trends in repeat data as concise charts and tables.

## Installation
PERF can be directly installed using pip with the package name `ssr-perf`. 
```bash
$ pip install ssr-perf
```

This name was chosen for the package so as not to clash with the existing `perf` command of Linux.

Alternatively, it can also be installed from the source code:
```bash
# Download the git repo
$ git clone https://github.com/RKMlab/perf.git

# Install
$ cd perf
$ python setup.py install
```
Both of the methods add a console command `ssr-perf`, which can be executed from any directory.

## Usage instructions
The help message and available options can be accessed using
```bash
$ ssr-perf -h # Short option
$ ssr-perf --help # Long option
```
which gives the following output
```
usage: ssr-perf [-h] -i <FILE> [-o <FILE>] [-a] [-rep <FILE>] [-m <INT>]
                [-M <INT>] [--min-length <INT> | --min-units INT or FILE]
                [--no-prefix] [--no-suffix] [--version]

Required arguments:
  -i <FILE>, --input <FILE>
                        Input file in FASTA format

Optional arguments:
  -o <FILE>, --output <FILE>
                        Output file name. Default: Input file name + _perf.tsv
  -a, --analyse         Generate a summary HTML report.
  -rep <FILE>, --repeats <FILE>
                        File with list of repeats (Not allowed with -m and/or
                        -M)
  -m <INT>, --min-motif-size <INT>
                        Minimum size of a repeat motif in bp (Not allowed with
                        -rep)
  -M <INT>, --max-motif-size <INT>
                        Maximum size of a repeat motif in bp (Not allowed with
                        -rep)
  --min-length <INT>    Minimum length cutoff of repeat
  --min-units INT or FILE
                        Minimum number of repeating units to be considered.
                        Can be an integer or a file specifying cutoffs for
                        different motif sizes.
  --no-prefix           Avoid optional prefixes. Only applicable with --min-
                        units)
  --no-suffix           Avoid optional suffixes. Only applicable with --min-
                        units)
  --version             show program's version number and exit
```
The details of each option are given below:

### `-i or --input`
This is the only required argument for the program. The input file must be a valid FASTA file. PERF uses [Biopython's](http://biopython.org/wiki/SeqIO) FASTA parser to read the input files. It accepts both single-line and multi-line sequences. Files with multiple sequences are also valid. To see more details about the FASTA format, see [this page](http://bioperl.org/formats/sequence_formats/FASTA_sequence_format).

### `-o or --output`
## HTML report
In addition to the default tab-separated output, PERF can also generate a fully interactive HTML report for easy downstream analysis of the repeat data. An example HTML report can be accessed [here](https://raw.githubusercontent.com/RKMlab/perf/html-report/test_data/test_input_perf.html)