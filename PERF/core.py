#! /usr/bin/env python

# pylint: disable=C0103, C0301

from __future__ import print_function, division
import sys
from os.path import splitext
import argparse
from tqdm import tqdm
from Bio import SeqIO
from collections import Counter
import gzip

if sys.version_info.major == 2:
    from utils import generate_repeats, get_ssrs, build_rep_set, univset, rawcharCount
    from analyse import analyse
elif sys.version_info.major == 3:
    from .utils import generate_repeats, get_ssrs, build_rep_set, univset, rawcharCount
    from .analyse import analyse

inf = float('inf')

def getArgs():
    """
    Parses command line arguments and returns them to the caller
    """
    __version__ = 'v0.2.5'
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input file in FASTA format')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<FILE>', default=sys.stdout, help='Output file name. Default: Input file name + _perf.tsv')
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')
    cutoff_group = optional.add_mutually_exclusive_group()
    cutoff_group.add_argument('-l', '--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    cutoff_group.add_argument('-u', '--min-units', metavar='INT or FILE', help="Minimum number of repeating units to be considered. Can be an integer or a file specifying cutoffs for different motif sizes.")
    optional.add_argument('-rep', '--repeats', type=argparse.FileType('r'), metavar='<FILE>', help='File with list of repeats (Not allowed with -m and/or -M)')
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-s', '--min-seq-length', type=int, metavar = '<INT>', default=0, help='Minimum size of sequence length for consideration (in bp)')
    optional.add_argument('-S', '--max-seq-length', type=float, metavar='<FLOAT>', default=inf, help='Maximum size of sequence length for consideration (in bp)')
    seqid_group = optional.add_mutually_exclusive_group()
    seqid_group.add_argument('-f', '--filter-seq-ids', metavar='<FILE>')
    seqid_group.add_argument('-F', '--target-seq-ids', metavar='<FILE>')
    optional.add_argument('--version', action='version', version='PERF ' + __version__)

    args = parser.parse_args()
    if args.repeats and (args.min_motif_size or args.max_motif_size):
        parser.error("-rep is not allowed with -m/-M")
    if args.repeats is None:
        if args.min_motif_size is None:
            args.min_motif_size = 1
        if args.max_motif_size is None:
            args.max_motif_size = 6
    if args.output.name == "<stdout>":
        args.output = open(splitext(args.input)[0] + '_perf.tsv', 'w')

    return args

def get_targetids(filter_seq_ids, target_seq_ids):
    target_ids = univset()
    if filter_seq_ids:
        target_ids = univset()
        filter_ids = []
        with open(filter_seq_ids) as fh:
            for line in fh:
                line = line.strip()
                filter_ids.append(line)
        target_ids = target_ids - set(filter_ids)
    elif target_seq_ids:
        target_ids = []
        with open(target_seq_ids) as fh:
            for line in fh:
                line = line.strip()
                target_ids.append(line)
        target_ids = set(target_ids)

    return target_ids

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
    min_seq_length = args.min_seq_length
    max_seq_length = args.max_seq_length
    target_ids = get_targetids(args.filter_seq_ids, args.target_seq_ids)
    print('Using length cutoff of %d' % (length_cutoff), file=sys.stderr)

    num_records = rawcharCount(seq_file, '>')
    if seq_file.endswith('gz'):
        handle = gzip.open(seq_file, 'rt')
    else:
        handle = open(seq_file, 'r')

    records = SeqIO.parse(handle, 'fasta')
    records = tqdm(records, total=num_records)
    for record in records:
        records.set_description("Processing %s" %(record.id))
        if  min_seq_length <= len(record.seq) <= max_seq_length and record.id in target_ids:
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
    min_seq_length = args.min_seq_length
    max_seq_length = args.max_seq_length
    target_ids = get_targetids(args.filter_seq_ids, args.target_seq_ids)
    print('Using unit cutoff of ', unit_cutoff, file=sys.stderr)

    num_records = rawcharCount(seq_file, '>')
    if seq_file.endswith('gz'):
        handle = gzip.open(seq_file, 'rt')
    else:
        handle = open(seq_file, 'r')

    records = SeqIO.parse(handle, 'fasta')
    records = tqdm(records, total=num_records)
    for record in records:
        records.set_description("Processing %s" %(record.id))
        if  (min_seq_length <= len(record.seq) <= max_seq_length) and record.id in target_ids:
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
