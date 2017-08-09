# PERF
## Introduction
PERF is a Python package developed for fast and accurate identification of microsatellites from DNA sequences. Microsatellites or Simple Sequence Repeats (SSRs) are short tandem repeats of 1-6nt motifs. They are present in all genomes, and have a wide range of uses and functional roles. The existing tools for SSR identification have one or more caveats in terms of speed, comprehensiveness, accuracy, ease-of-use, flexibility and memory usage. PERF was designed to address all these problems.

PERF is a recursive acronym that stands for "PERF is an Exhaustive Repeat Finder". It is compatible with both Python 2 (test on Python 2.7) and 3 (test on Python 3.5). Its key features are:
  - Fast run time, despite being a single-threaded application. As an example, all SSRs from the entire human genome takes less than 10 minutes. The speed can be approved 5-fold using PyPy (tested on PyPy 5.8)
  - Linear time and space complexity (O(n))
  - 100% accurate and comprehensive - Does not miss any repeats or does not pick any incorrect ones
  - Easy to use - The only required argument is the input DNA sequence in FASTA format
  - Flexible. While having sensible defaults, most of the parameters are customizable by the user at runtime.
  - Repeat cutoffs can be specified either in terms of the total repeat length or in terms of number of repeating units.
  - TSV and HTML report. The default output is an easily parseable and exportable tab-separated format. Optionally, PERF also generates an interactive HTML report that depicts trends in repeat data as concise charts and tables.

## Installation
PERF can be directly installed using pip with the package name `ssr-perf`. 
```bash
pip install ssr-perf
```

This name was chosen for the package so as not to clash with the existing `perf` command of Linux.

Alternatively, it can also be installed from the source code:
```bash
# Download the git repo
git clone https://github.com/RKMlab/perf.git

# Install
cd perf
python setup.py install
```
Both of the methods add a console command `ssr-perf`, which can be executed from any directory.

## Usage instructions

## HTML report
In addition to the default tab-separated output, PERF can also generate a fully interactive HTML report for easy downstream analysis of the repeat data. An example HTML report can be accessed [here](https://raw.githubusercontent.com/RKMlab/perf/html-report/test_data/test_input_perf.html)