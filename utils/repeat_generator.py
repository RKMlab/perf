#! /usr/bin/env python

# pylint: disable=C0103
from __future__ import print_function
import sys
import argparse
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', default=1, help='Minimum size of a repeat motif in bp')
parser.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', default=6, help='Maximum size of a repeat motif in bp')
parser.add_argument('-fo', '--out', type=argparse.FileType('w'), metavar='<FILE>', default=sys.stdout, help='Output file')
args = parser.parse_args()

def rev_comp(string):
    complement = string.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]


def expand_repeat(string, size):
    return_string = ''
    i = 0
    while len(return_string) < size:
        return_string += string[i]
        i += 1
        if i >= len(string):
            i = 0
    return return_string


def get_cycles(string):
    cycles = []
    for i in range(len(string)):
        cycles.append(string[i:] + string[:i])
    return cycles


def generate_repeats(min_size, max_size, output_file):
    alphabet = ['A', 'C', 'G', 'T']
    expanded_set = set()
    repeat_set = set()
    for i in range(min_size, max_size+1):
        for combination in product(alphabet, repeat=i):
            repeat = ''.join(combination)
            repeat_revcomp = rev_comp(repeat)
            expanded = expand_repeat(repeat, max_size)
            if expanded in expanded_set:
                continue
            else:
                repeat_cycles = get_cycles(repeat)
                for cycle in repeat_cycles:
                    strand = '+'
                    string = expand_repeat(cycle, max_size)
                    expanded_set.add(string)
                    if cycle not in repeat_set:
                        repeat_set.add(cycle)
                        print(cycle, repeat, str(len(cycle)), strand, sep='\t', file=output_file)
                if repeat_revcomp == repeat:
                    continue
                repeat_cycles = get_cycles(repeat_revcomp)
                for cycle in repeat_cycles:
                    strand = '-'
                    string = expand_repeat(cycle, max_size)
                    expanded_set.add(string)
                    if cycle not in repeat_set:
                        repeat_set.add(cycle)
                        print(cycle, repeat, str(len(cycle)), strand, sep='\t', file=output_file)

min_motif_size = args.min_motif_size
max_motif_size = args.max_motif_size
generate_repeats(min_motif_size, max_motif_size, args.out)
