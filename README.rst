PERF
====

Introduction
------------

PERF is a Python package developed for fast and accurate identification of microsatellites from DNA sequences. Microsatellites or Simple Sequence Repeats (SSRs) are short tandem repeats of 1-6nt motifs. They are present in all genomes, and have a wide range of uses and functional roles. The existing tools for SSR identification have one or more caveats in terms of speed, comprehensiveness, 
accuracy, ease-of-use, flexibility and memory usage. PERF was designed to address all these problems.

PERF is a recursive acronym that stands for “PERF is an Exhaustive Repeat Finder”. It is compatible with both Python 2 (tested on Python 2.7) and 3 (tested on Python 3.5). Its key features are:

 - Fast run time, despite being a single-threaded application. As an example, identification of all SSRs from the entire human genome takes less than 7 minutes. The speed can be further improved ~3 to 4 fold using (human genome finishes in less than 2 minutes using PyPy v5.8.0) 
 - Linear time and space complexity (O(n))
 - Identifies perfect SSRs
 - 100% accurate and comprehensive
 - Does not miss any repeats or does not pick any incorrect ones
 - Easy to use
 - The only required argument is the input DNA sequence in FASTA format
 - Flexible
 - Most of the parameters are customizable by the user at runtime
 - Repeat cutoffs can be specified either in terms of the total repeat length or in terms of number of repeating units
 - TSV output and HTML report

The default output is an easily parseable and exportable tab-separated format. Optionally, PERF also generates an interactive HTML report that depicts trends in repeat data as concise charts and tables.

Installation
------------

PERF can be directly installed using pip with the package name
``perf_ssr``.

.. code:: bash

    $ pip install perf_ssr

This name was chosen for the package so as not to clash with the existing ``perf`` package.

Alternatively, it can also be installed from the source code:

.. code:: bash

    # Download the git repo
    $ git clone https://github.com/RKMlab/perf.git

    # Install
    $ cd perf
    $ python setup.py install

Both of the methods add a console command ``PERF``, which can be executed from any directory. It can also be used without installation by running the ``core.py`` file in the ``PERF`` subfolder:

.. code:: bash

    $ git clone https://github.com/RKMlab/perf.git
    $ cd perf/PERF
    $ python core.py -h # Print the help message of PERF (see below)

Usage
-----

The help message and available options can be accessed using

.. code:: bash

    $ PERF -h # Short option
    $ PERF --help # Long option

which gives the following output

::

    usage: PERF [-h] -i <FILE> [-o <FILE>] [-a] [-l <INT> | -u INT or FILE]
                [-rep <FILE>] [-m <INT>] [-M <INT>] [--version]

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
      --version             show program's version number and exit