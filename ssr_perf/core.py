#! /usr/bin/env python

# pylint: disable=C0103

from __future__ import print_function
import sys
from os.path import splitext
import argparse
from tqdm import tqdm
from Bio import SeqIO
from utils import generate_repeats, get_ssrs, build_rep_set
from analyse import analyse

def getArgs():
    """
    Parses command line arguments and returns them to the caller
    """
    __version__ = 'v0.2.0'
    parser = argparse.ArgumentParser()
    helper = parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input file in FASTA format')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<FILE>', default=sys.stdout, help='Output file name. Default: Input file name + _perf.tsv')
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')
    optional.add_argument('-rep', '--repeats', type=argparse.FileType('r'), metavar='<FILE>', help='File with list of repeats (Not allowed with -m and/or -M)')
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', default=1, help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', default=6, help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    cutoff_group = optional.add_mutually_exclusive_group()
    cutoff_group.add_argument('--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    cutoff_group.add_argument('--min-units', metavar='INT or FILE', help="Minimum number of repeating units to be considered. Can be an integer or a file specifying cutoffs for different motif sizes.")
    optional.add_argument('--no-prefix', action='store_true', help='Avoid optional prefixes. Only applicable with --min-units)')
    optional.add_argument('--no-suffix', action='store_true', help='Avoid optional suffixes. Only applicable with --min-units)')
    # optional.add_argument('-h', '--help', action='help', help="Show this help message and exit")
    optional.add_argument('--version', action='version', version='ssr-perf ' + __version__)

    # repeats_input = parser.add_mutually_exclusive_group()
    # motif_size_group = repeats_input.add_argument_group('motif_size')


    args = parser.parse_args()
    if args.repeats and (args.min_motif_size or args.max_motif_size):
        parser.error("-rep is not allowed with -m/-M")
    if args.output.name == "<stdout>":
        args.output = open(splitext(args.input)[0] + '_perf.tsv', 'w')
    return args

def getSSRNative(args):
    """
    Identifies microsatellites using native string matching.
    As the entire sequence is scanned in a single iteration, the speed is vastly improved
    """
    length_cutoff = args.min_length
    repeat_file = args.repeats
    seq_file = args.input
    out_file = args.output
    repeats_info = build_rep_set(repeat_file, length_cutoff=length_cutoff)
    repeat_set = set(repeats_info.keys())
    print('Using length cutoff of %d' % (length_cutoff), file=sys.stderr)

    num_records = 0
    # with open(seq_file, "rt") as handle:
    #     records = SeqIO.parse(handle, 'fasta')
    #     for r in records:
    #         num_records += 1

    with open(seq_file, "rt") as handle:
        records = SeqIO.parse(handle, 'fasta')
        for record in tqdm(records, total=num_records):
            get_ssrs(record, repeats_info, repeat_set, out_file)
    out_file.close()


def getSSR_units(args, unit_cutoff):
    """
    Identifies microsatellites using native string matching.
    The repeat length cutoffs vary for different motif sizes.
    """
    repeat_file = args.repeats
    seq_file = args.input
    out_file = args.output
    
    repeats_info = build_rep_set(repeat_file, unit_cutoff=unit_cutoff)
    repeat_set = set(repeats_info.keys())
    
    print('Using unit cutoff of ', unit_cutoff, file=sys.stderr)

    num_records = 0
    # with open(seq_file, "rt") as handle:
    #     records = SeqIO.parse(handle, 'fasta')
    #     for r in records:
    #         num_records += 1

    with open(seq_file, "rt") as handle:
        records = SeqIO.parse(handle, 'fasta')
        for record in tqdm(records, total=num_records):
            get_ssrs(record, repeats_info, repeat_set, out_file)
    out_file.close()


def main():
    """
    Main function of the script
    """
    args = getArgs()

    # User specifies motif size range instead of giving a repeats file
    if args.repeats is None:
        min_motif_size = args.min_motif_size
        max_motif_size = args.max_motif_size
        args.repeats = generate_repeats(min_motif_size, max_motif_size)


    # User specifies minimum length
    if args.min_length:
        getSSRNative(args)

    # User specific minimum number of units
    elif args.min_units:
        unit_cutoff = dict()
        try:
            args.min_units = int(args.min_units)
            unit_cutoff[0] = args.min_units
        except ValueError:
            try:
                with open(args.min_units, 'r') as unitsIn:
                    for line in unitsIn:
                        L = line.strip().split()
                        try:
                            L[0] = int(L[0])
                            L[1] = int(L[1])
                            if L[1] == 1:
                                print('Warning: Repeat unit of 1 used for size %d.' % (L[0]), file=sys.stderr)
                            unit_cutoff[L[0]] = L[1]
                        except ValueError:
                            sys.exit('Invalid file format given for minimum units. Refer to help for more details')
            except FileNotFoundError:
                sys.exit('Units file specified is not found. Please provide a valid file')
        getSSR_units(args, unit_cutoff)

    # Default settings
    elif args.min_length is None and args.min_units is None:
        args.min_length = 12
        getSSRNative(args)

    # Specifies to generate a HTML report
    if args.analyse:
        analyse(args)


if __name__ == '__main__':
    main()