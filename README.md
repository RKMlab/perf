# PERF
[![Build](https://img.shields.io/badge/Build-passing-brightgreen.svg)]()
[![PyPI](https://img.shields.io/badge/PyPI-v0.4.5-blue.svg)]()
[![License](https://img.shields.io/badge/Licence-MIT-blue.svg)]()
## Introduction
PERF is a Python package developed for fast and accurate identification of microsatellites from DNA sequences. Microsatellites or Simple Sequence Repeats (SSRs) are short tandem repeats of 1-6nt motifs. They are present in all genomes, and have a wide range of uses and functional roles. The existing tools for SSR identification have one or more caveats in terms of speed, comprehensiveness, accuracy, ease-of-use, flexibility and memory usage. PERF was designed to address all these problems.

PERF is a recursive acronym that stands for "PERF is an Exhaustive Repeat Finder". It is compatible with both Python 2 (tested on Python 2.7) and 3 (tested on Python 3.5). Its key features are:
  - Fast run time. As an example, identification of all SSRs from the entire human genome takes less than 7 minutes. The speed can be further improved ~3 to 4 fold using [PyPy](https://pypy.org/) (human genome finishes in less than 2 minutes using PyPy v5.8.0)
  - Linear time and space complexity (O(n))
  - Identifies perfect SSRs
  - 100% accurate and comprehensive - Does not miss any repeats or does not pick any incorrect ones
  - Easy to use - The only required argument is the input DNA sequence in FASTA format
  - Flexible - Most of the parameters are customizable by the user at runtime
  - Repeat cutoffs can be specified either in terms of the total repeat length or in terms of number of repeating units
  - TSV output and HTML report. The default output is an easily parseable and exportable tab-separated format. Optionally, PERF also generates an interactive HTML report that depicts trends in repeat data as concise charts and tables

## Change log 

## [0.4.6] - 2021-04-22
### Fixes
 - Fixed usage of unit options file input for fastq input.
 - Fixed usage of repeats input file.

## [0.4.5] - 2020-05-07
### Added
 - Annotation of repeats w.r.t to genomic context using a GFF or GTF file. (option -g).
 - Multi-threading. Parallel identification of repeats in different sequences.
 - Identification of perfect repeats in fastq files.
 - Analysis report for repeats in fastq files.
 - Option to identify atomic repeats.

### Changed
 - Analysis report rebuilt with Semantic ui and Apex Charts.
 - Visualisation of repeat annotation data in analysis report.

### Fixes 
 - Python2 compatability fixed.
 - Bug fixes for PyPi compatability.
 - Import error issues.

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
usage: core.py [-h] -i <FILE> [-o <FILE>] [--format <STR>] [--version]
               [-rep <FILE>] [-m <INT>] [-M <INT>] [-s <INT>] [-S <FLOAT>]
               [--include-atomic] [-l <INT> | -u INT or FILE] [-a] [--info]
               [-g <FILE>] [--anno-format <STR>] [--gene-key <STR>]
               [--up-promoter <INT>] [--down-promoter <INT>]
               [-f <FILE> | -F <FILE>] [-t <INT>]

Required arguments:
  -i <FILE>, --input <FILE>
                        Input sequence file.

Optional arguments:
  -o <FILE>, --output <FILE>
                        Output file name. Default: Input file name + _perf.tsv
  --format <STR>        Input file format. Default: fasta, Permissible: fasta,
                        fastq
  --version             show program's version number and exit
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
  --include-atomic      An option to include factor atomic repeats for minimum
                        motif sizes greater than 1.
  -l <INT>, --min-length <INT>
                        Minimum length cutoff of repeat
  -u INT or FILE, --min-units INT or FILE
                        Minimum number of repeating units to be considered.
                        Can be an integer or a file specifying cutoffs for
                        different motif sizes.
  -a, --analyse         Generate a summary HTML report.
  --info                Sequence file info recorded in the output.
  -f <FILE>, --filter-seq-ids <FILE>
                        List of sequence ids in fasta file which will be
                        ignored.
  -F <FILE>, --target-seq-ids <FILE>
                        List of sequence ids in fasta file which will be used.
  -t <INT>, --threads <INT>
                        Number of threads to run the process on. Default is 1.

Annotation arguments:
  -g <FILE>, --annotate <FILE>
                        Genic annotation input file for annotation, Both GFF
                        and GTF can be processed. Use --anno-format to specify
                        format.
  --anno-format <STR>   Format of genic annotation file. Valid inputs: GFF,
                        GTF. Default: GFF
  --gene-key <STR>      Attribute key for geneId. The default identifier is
                        "gene". Please check the annotation file and pick a
                        robust gene identifier from the attribute column.
  --up-promoter <INT>   Upstream distance(bp) from TSS to be considered as
                        promoter region. Default 1000
  --down-promoter <INT>
                        Downstream distance(bp) from TSS to be considered as
                        promoter region. Default 1000
```
The details of each option are given below:

### `-i or --input`
**Expects:** *FILE*<br>
**Default:** *None*<br>
This is the only required argument for the program. The input file must be a valid FASTA/FASTQ file. PERF uses [Biopython's](http://biopython.org/wiki/SeqIO) FASTA parser to read the input fasta files. It accepts both single-line and multi-line sequences. Files with multiple sequences are also valid. To see more details about the FASTA format, see [this page](http://bioperl.org/formats/sequence_formats/FASTA_sequence_format).

### `-o or --output`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *Input Filename + _perf.tsv (see below)*<br>
If this option is not provided, the default output filename will be the same as the input filename, with its extension replaced with '_perf.tsv'. For example, if the input filename is `my_seq.fa`, the default output filename will be `my_seq_perf.tsv`. If the input filename does not have any extension, `_perf.tsv` will be appended to the filename. Please note that even in the case of no identified SSRs, the output file is still created (therefore overwriting any previous file of the same name) but with no content in the file.
#### Output for fasta
The output is a tab-delimited file, with one SSR record per line. 
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

#### Output for fastq
The output is a tab-delimited file, with data on each repeat class per line.
| S.No | Column | Description |
|:----:| ------ | ----------- |
| 1 | Repeat Class | Class of repeat as grouped by their cyclical variations |
| 2 | Number of reads | Number of reads having an instance of the repeat |
| 3 | Frequency | Total number of instances of the repeat  |
| 4 | Bases | Total number of bases covered by the repeat |
| 5 | Repeat reads per million reads | Number of  |
| 6 | Instances per million reads | Strand of SSR based on their cyclical variation |
| 7 | Repeat Bases per MB of sequence | Number of times the base motif is repeated |
| 8 | Length distribution | Starting sequence of the SSR irrespective of Repeat class and strand|
| 9 | Motif distribution | Starting sequence of the SSR irrespective of Repeat class and strand|


### `--format`
**Expects:** *STRING (specifying format of the file)*<br>
**Default:** *fasta*<br>
PERF was originally developed to identify repeats in FASTA files. In version 4.0.0 PERF can identify repeats in FASTQ sequence files as well. The default format the program expects is fasta. Specify input format as 'fasta' for FASTA files and 'fastq' for FASTQ files.

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

### `--include-atomic`
**Expects:** *None*<br>
**Default:** *False*<br>
Searches for atomic repeats when set to *True*. For example, when minimum motif size is set to 2bp, PERF ignores monomer repeats. When include atomic repeats is set to *True*, PERF identifies AA, CC, GG and TT as dimer repeats.

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

### `--info`
**Expects:** *None*<br>
**Default:** *False*<br>
This option when set to *True*, includes information about the input sequence files and repeat summary data in the output file.

```bash
$ tail -5 test_input_perf.tsv
gi|514486271|gb|KE346361.1|	2667759	2667775	ATC	16	+	5	CAT
#File_name: test_input.fa
#Total_sequences: 2
#Total_bases: 6462134
#GC: 53.970000
```

### `-g or --annotate`
**Expects:** *FILE*<br>
**Default:** *None*<br>
Input a genomic feature file to annotate the repeats in the genomic context. PERF accepts both GFF and GTF format genomic feature files. Each repeat is annotated w.r.t the closest gene and classified either as Genic, Exonic, Intronic and Intergenic according to the position of the repeat. Besides this, the repeat is also checked if it falls in the promoter region of the gene. Annotation adds 7 columns to the default perf output which already consist 8 columns.

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 9 | Gene name | Name of the closest gene |
| 10 | Gene Start | Start position of gene in the Chromosome |
| 11 | Gene Stop | End position of gene in the Chromosome |
| 12 | Strand | The strand orientation of the gene |
| 13 | Genomic annotation | Annotation of the repeat w.r.t to the gene. Possible annotations are {Genic, Exonic, Intronic, Intergenic} |
| 14 | Promoter annotation | If repeat falls in the promoter region of the closest gene. The default promoter region is 1Kb upstream and downstream of TSS. |
| 15 | Distance from TSS | Distance of the repeat from the TSS of the gene. |

### `--anno-format`
**Expects:** *STRING*<br>
**Default:** *GFF*<br>
Option to specify the format of the input genomic feature file. Accepted  inputs are GFF or GTF. More details about the GFF and GTF formats can be found [here](https://asia.ensembl.org/info/website/upload/gff.html).

### `--gene-key`
**Expects:** *STRING*<br>
**Default:** *gene*<br>
The attribute key used for the name of the gene in the GFF/GTF file. In the below example GFF file, we have the location of a gene and it's mRNA and exon locations. The last column of the file specifies attributes associated with each feature, like ID, Parent, gene etc. PERF uses on of the attribute to identify the gene and also it's exons. In th below example the key "gene" can be used to identify gene and the exons of the gene as they have the same gene name. Please check your GFF/GTF file for a robust attribute key which can identify all genes and their corresponding exons. We are actively working on better annotation where we can identify genes and their exons based on the ID and Parent.

```
# Sample GFF
NC_004354.4	RefSeq	gene	124370	126714	.	-	.	ID=gene1;Name=CG17636;gbkey=Gene;gene=CG17636;gene_biotype=protein_coding;gene_synonym=DmelCG17636,EG:23E12.1;
NC_004354.4	RefSeq	mRNA	124370	126714	.	-	.	ID=rna1;Parent=gene1;Name=NM_001103384.3;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	126626	126714	.	-	.	ID=id13;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	125495	126259	.	-	.	ID=id14;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
```

### `--up-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Upstream distance(bp) from the TSS of the gene to be considered as promoter region. Default 1000.

### `--down-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Downstream distance(bp) from the TSS of the gene to be considered as promoter region. Default 1000.

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
# Specify fastq format
$ PERF -i my_seq.fastq --format fastq
```

#### Generate Analysis Report
``` bash
# Find all monomer to hexamer repeats of >=12nt length and generate an HTML report
$ PERF -i my_seq.fa -a
# Specify output filename
$ PERF -i my_seq.fa -o PERF_out.tsv -a # HTML file is called PERF_out.html
```

#### Annotate Repeats
``` bash
# Find all monomer to hexamer repeats of >=12nt length and annotate them
$ PERF -i my_seq.fa -g my_seq.gff
# Specify feature file format and set downstream promoter region to 500bp
$ PERF -i my_seq.fa -g my_seq.gtf --anno-format gtf --down-promoter 500 # HTML file is called PERF_out.html
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

