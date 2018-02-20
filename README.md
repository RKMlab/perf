# PERF
[![Build](https://img.shields.io/badge/Build-passing-brightgreen.svg)]()
[![PyPI](https://img.shields.io/badge/PyPI-v0.2.5-blue.svg)]()
[![License](https://img.shields.io/badge/Licence-MIT-blue.svg)]()
## Introduction
PERF is a Python package developed for fast and accurate identification of microsatellites from DNA sequences. Microsatellites or Simple Sequence Repeats (SSRs) are short tandem repeats of 1-6nt motifs. They are present in all genomes, and have a wide range of uses and functional roles. The existing tools for SSR identification have one or more caveats in terms of speed, comprehensiveness, accuracy, ease-of-use, flexibility and memory usage. PERF was designed to address all these problems.

PERF is a recursive acronym that stands for "PERF is an Exhaustive Repeat Finder". It is compatible with both Python 2 (tested on Python 2.7) and 3 (tested on Python 3.5). Its key features are:
  - Fast run time, despite being a single-threaded application. As an example, identification of all SSRs from the entire human genome takes less than 7 minutes. The speed can be further improved ~3 to 4 fold using [PyPy](https://pypy.org/) (human genome finishes in less than 2 minutes using PyPy v5.8.0)
  - Linear time and space complexity (O(n))
  - Identifies perfect SSRs
  - 100% accurate and comprehensive - Does not miss any repeats or does not pick any incorrect ones
  - Easy to use - The only required argument is the input DNA sequence in FASTA format
  - Flexible - Most of the parameters are customizable by the user at runtime
  - Repeat cutoffs can be specified either in terms of the total repeat length or in terms of number of repeating units
  - TSV output and HTML report. The default output is an easily parseable and exportable tab-separated format. Optionally, PERF also generates an interactive HTML report that depicts trends in repeat data as concise charts and tables

## Installation
PERF can be directly installed using pip with the package name `perf_ssr`. 
```bash
$ pip install perf_ssr
```

This name was chosen for the package so as not to clash with the existing `perf` package.

Alternatively, it can be installed from the source code:
```bash
# Download the git repo
$ git clone https://github.com/RKMlab/perf.git

# Install
$ cd perf
$ python setup.py install
```
Both of the methods add a console command `PERF`, which can be executed from any directory. It can also be used without installation by running the `core.py` file in the `PERF` subfolder:

```bash
$ git clone https://github.com/RKMlab/perf.git
$ cd perf/PERF
$ python core.py -h # Print the help message of PERF (see below)
```

## Usage
The help message and available options can be accessed using
```bash
$ PERF -h # Short option
$ PERF --help # Long option
```
which gives the following output
```
usage: PERF [-h] -i <FILE> [-o <FILE>] [-a] [-l <INT> | -u INT or FILE]
            [-rep <FILE>] [-m <INT>] [-M <INT>] [-s <INT>] [-S <FLOAT>]
            [-f <FILE> | -F <FILE>] [--version]

Required arguments:
  -i <FILE>, --input <FILE>
                        Input file in FASTA format

Optional arguments:
  -o <FILE>, --output <FILE>
                        Output file name. Default: Input file name + _perf.tsv
  -a, --analyse         Generate a summary HTML report.
  -l <INT>, --min-length <INT>
                        Minimum length cutoff of repeat
  -u INT or FILE, --min-units INT or FILE
                        Minimum number of repeating units to be considered.
                        Can be an integer or a file specifying cutoffs for
                        different motif sizes.
  -rep <FILE>, --repeats <FILE>
                        File with list of repeats (Not allowed with -m and/or
                        -M)
  -m <INT>, --min-motif-size <INT>
                        Minimum size of a repeat motif in bp (Not allowed with
                        -rep)
  -M <INT>, --max-motif-size <INT>
                        Maximum size of a repeat motif in bp (Not allowed with
                        -rep)
  -s <INT>, --min-seq-length <INT>
                        Minimum size of sequence length for consideration (in
                        bp)
  -S <FLOAT>, --max-seq-length <FLOAT>
                        Maximum size of sequence length for consideration (in
                        bp)
  -f <FILE>, --filter-seq-ids <FILE>
  -F <FILE>, --target-seq-ids <FILE>
  --version             show program's version number and exit
```
The details of each option are given below:

### `-i or --input`
**Expects:** *FILE*<br>
**Default:** *None*<br>
This is the only required argument for the program. The input file must be a valid FASTA file. PERF uses [Biopython's](http://biopython.org/wiki/SeqIO) FASTA parser to read the input files. It accepts both single-line and multi-line sequences. Files with multiple sequences are also valid. To see more details about the FASTA format, see [this page](http://bioperl.org/formats/sequence_formats/FASTA_sequence_format).

### `-o or --output`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *Input Filename + _perf.tsv (see below)*<br>
The output is a tab-delimited file, with one SSR record per line. If this option is not provided, the default output filename will be the same as the input filename, with its extension replaced with '_perf.tsv'. For example, if the input filename is `my_seq.fa`, the default output filename will be `my_seq_perf.tsv`. If the input filename does not have any extension, `_perf.tsv` will be appended to the filename. Please note that even in the case of no identified SSRs, the output file is still created (therefore overwriting any previous file of the same name) but with no content in the file.

The output columns follow the [BED](https://genome.ucsc.edu/FAQ/FAQformat.html) format. The details of the columns are given below:

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 1 | Chromosome | Chromosome or Sequence Name as specified by the first word in the FASTA header |
| 2 | Repeat Start | 0-based start position of SSR in the Chromosome |
| 3 | Repeat Stop | End position of SSR in the Chromosome |
| 4 | Repeat Class | Class of repeat as grouped by their cyclical variations |
| 5 | Repeat Length | Total length of identified repeat in nt |
| 6 | Repeat Strand | Strand of SSR based on their cyclical variation |
| 7 | Motif Number | Number of times the base motif is repeated |
| 8 | Actual Repeat | Starting sequence of the SSR irrespective of Repeat class and strand|

An example output showing some of the largest repeats from *Drosophila melanogaster* is given below
```
X       22012826  22014795  ACTGGG  1969    -       328     TCCCAG
2RHet   591337    591966    AATACT  629     -       104     ATTAGT
4       1042143   1042690   AAATAT  547     +       91      AAATAT
2RHet   598244    598789    AATACT  545     -       90      AGTATT
XHet    122       663       AGAT    541     +       135     GATA
X       22422335  22422827  AGAT    492     +       123     GATA
3R      975265    975710    AAAT    445     -       111     TTAT
X       15442288  15442724  ACAGAT  436     +       72      ACAGAT
2L      22086818  22087152  AATACT  334     -       55      TATTAG
YHet    137144    137466    AAGAC   322     -       64      CTTGT
```

### `-a or --analyze`
**Expects:** *None*<br>
**Default:** *False*<br>
In addition to the default tab-separated output, PERF can also generate a fully interactive HTML report for easy downstream analysis of the repeat data. The filename will be the same prefix as that of the main output. For example, if the input filename was `my_seq.fa`, the analysis report will be  `my_seq_perf.html`. An example HTML report, generated from the repeat data of *Homo sapiens* (build hg19), can be accessed [here](https://raw.githubusercontent.com/RKMlab/perf/html-report/test_data/Human_hg19_perf.html) (Right click -> Save As).

### `-l or --min-length`
**Expects:** *INTEGER*<br>
**Default:** *12*<br>
Minimum length cut-off to be considered when finding an SSR. The same cut-off will apply for SSRs of all motif lengths, even if the motif length is not a divisor of this value. In such cases, SSRs that end with a partial motif are also picked if they pass the length cut-off.

### `-u or --min-units`
**Expects:** *INTEGER* OR *FILE*<br>
**Default:** *None*<br>
This option finds SSRs with a minimum number of repeating motifs. The argument accepts either an integer or file. If an integer is specified, the same value is used for all motif lengths. Instead, a specific value for each motif length using a two-column tab-separated file as demonstrated below:

```bash
$ cat repeat_units.txt
1	10
2	6
3	4
4	3
5	2
6	2
```

The first column specifies the motif length, and the second column specifies the minimum number of times the motif should be repeated to be considered an SSR. This file can be used to identify repeats with different number of repeating motifs: monomers repeated at least 10 times, dimers repeated at least 6 times etc., using the following command
``` bash
$ PERF -i my_seq.fa -m 1 -M 6 -u repeat_units.txt
```

### `-rep or --repeats`
**Expects:** *FILE*<br>
**Default:** *None*<br>
PERF provides an option to limit the search to specific repeat motifs. The repeats of interest should be specified via a file containing 4 tab-separated columns, as shown below:

```bash
$ cat my_repeats.txt
A   A   1   +                                                                
T   A   1   -
AG  AG  2   +
CT  AG  2   -
GA  AG  2   +
TC  AG  2   -
$ PERF -i my_seq.fa -rep my_repeats.txt # Find all A and AG repeats from my_seq.fa
```

**Note:** This option is not allowed when `-m` or `-M` options are used.
### `-m or --min-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *1*<br>
Minimum length of motifs to be considered. By default, PERF ignores redundant motifs. For example, a stretch of 12 A's is considered a monomer repeat of 12 A's rather than a dimer repeat of 6 AA's. However, this is only true if `-m` is set to 1. If for example, `-m` is set to 2, then stretches of 12 A's are reported as dimer AA repeats. If this behavior isn't desired, we suggest using the `-rep` option (see above) to specify the motifs that should/shouldn't be included.

### `-M or --max-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *6*<br>
Maximum length of motifs to be considered. Setting a large value of `-M` has a non-trivial effect on both the runtime and memory usage of PERF. This is noticeable with `-M` values above 10.

### `-s or --min-seq-length`
**Expects:** *INTEGER*<br>
**Default:** *0*<br>
Minimum length of the input sequence to be searched for SSRs in bp. All sequences in the input file that are smaller than this length will be ignored.

### `-S or --max-seq-length`
**Expects:** *INTEGER*<br>
**Default:** *Infinity*<br>
Maximum length of the input sequence to be searched for SSRs in bp. All sequences in the input file that are larger than this length will be ignored.

### `-f or --filter-seq-ids`
**Expects:** *FILE*<br>
**Default:** *None*<br>
This option accepts a file with a list of sequence IDs in the input file that should be ignored. Useful for ignoring contigs, scaffold, or other poor quality sequences. The IDs can be FASTA headers (starting with '>' symbol) or just the names without the '>' symbol.

### `-F or --target-seq-ids`
**Expects:** *FILE*<br>
**Default:** *None*<br>
This option accepts a file with a list of sequence IDs in the input file that should be analyzed. All other sequences will be ignored. Useful for analyzing specific chromosomes from a large input file. The IDs can be FASTA headers (starting with '>' symbol) or just the names without the '>' symbol.

### `--version`
Prints the version info of PERF.

## Examples

The following examples assume that the file with input sequence in FASTA format is named `my_seq.fa`.

#### Basic Usage
``` bash
# Find all monomer to hexamer repeats of >=12nt length
$ PERF -i my_seq.fa
# Specify output filename
$ PERF -i my_seq.fa -o PERF_output.tsv
```

#### Generate Analysis Report
``` bash
# Find all monomer to hexamer repeats of >=12nt length and generate an HTML report
$ PERF -i my_seq.fa -a
# Specify output filename
$ PERF -i my_seq.fa -o PERF_out.tsv -a # HTML file is called PERF_out.html
```

#### Set Cut-off Criteria
```bash
# Find all monomer to hexamer repeats of >=15nt length
$ PERF -i my_seq.fa -l 15
# Find SSRs with at least 6 repeating motifs (for all motif lengths)
$ PERF -i my_seq.fa -u 6
```

#### Identify Specific Repeats
``` bash
$ cat my_repeats.txt
AG  AG  2   +
CT  AG  2   -
GA  AG  2   +
TC  AG  2   -
# Find all AG repeats and generate an HTML report
$ PERF -i my_seq.fa -rep my_repeats.txt -a
```

#### Change Motif Length Cut-offs
```bash
# Ignore monomer and dimer repeats, and repeats with <4 repeating units
$ PERF -i my_seq.fa -m 3 -u 4
# Report only tetramer repeats of >=16nt length, and generate HTML report
$ PERF -i my_seq.fa -m 4 -M 4 -l 16 -a

```

In all the above examples, the output of PERF is saved to `my_seq_perf.tsv` and the HTML report is saved to `my_seq_perf.html` unless `-o` is specified.


## Citation

If you find PERF useful for your research, please cite it as follows:

PERF: an exhaustive algorithm for ultra-fast and efficient identification of microsatellites from large DNA sequences<br>
*Akshay Kumar Avvaru, Divya Tej Sowpati, Rakesh Kumar Mishra*<br>
Bioinformatics, , btx721<br>
[doi](https://doi.org/10.1093/bioinformatics/btx721): 10.1093/bioinformatics/btx721

## Contact
For queries or suggestions, please contact:

Divya Tej Sowpati - <tej@ccmb.res.in><br>
Akshay Kumar Avvaru - <avvaru@ccmb.res.in>

