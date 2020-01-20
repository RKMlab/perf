#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
from itertools import product
from utils import rev_comp


def expand_repeat(string, size):
    """Expands a motif to highest motif size, used for checking duplicates"""
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


def generate_repeats(min_size, max_size, atomic):
    """Generates all possible motifs for repeats in a given length range"""
    generated_repeats = []
    alphabet = ['A', 'C', 'G', 'T']
    expanded_set = set()
    repeat_set = set()
    init_size = 1
    if atomic:
        init_size = min_size
    for i in range(init_size, max_size+1):
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
                        if len(cycle) >= min_size:
                            generated_repeats.append('\t'.join([cycle, repeat, str(len(cycle)), strand]))
                if repeat_revcomp == repeat:
                    continue
                repeat_cycles = get_cycles(repeat_revcomp)
                for cycle in repeat_cycles:
                    strand = '-'
                    string = expand_repeat(cycle, max_size)
                    expanded_set.add(string)
                    if cycle not in repeat_set:
                        repeat_set.add(cycle)
                        if len(cycle) >= min_size:
                            generated_repeats.append('\t'.join([cycle, repeat, str(len(cycle)), strand]))
    return generated_repeats


def build_rep_set(repeat_file, length_cutoff=None, unit_cutoff=None):
    """
        Outputs the repeats info dictionary used by the get_ssrs function.
        Takes list of repeat motifs from repeats file(output by generate_repeats function) as input.
        Creates a dictionary with expanded repeat as the key and (class, motif_length, strand) as values.
        Works either by "length_cutoff=" or by "unit_cutoff=" arguments.
    """
    repeats_out = dict()
    motif_fallback = dict()
    repeat_lengths = set()
    if length_cutoff is not None:
        motif_lengths = set()
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            motif_lengths.add(motif_length)
            motif = expand_repeat(motif, length_cutoff)
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        for m in motif_lengths:
            motif_fallback[m] = max(motif_lengths)
        repeats_out['fallback'] = motif_fallback
        repeats_out['rep_lengths'] = [length_cutoff]

    elif unit_cutoff is not None:
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            try:
                motif = motif*unit_cutoff[motif_length]
            except KeyError:
                motif = motif*unit_cutoff[0]
            repeat_lengths.add(len(motif))
            motif_fallback[motif_length] = len(motif) - 1
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        repeat_lengths = sorted(list(repeat_lengths))
        repeats_out['rep_lengths'] = repeat_lengths
        repeats_out['fallback'] = motif_fallback

    return repeats_out



def get_ssrs(seq_record, repeats_info, out_file):
    """Native function that identifies repeats in fasta files."""
    
    repeat_lengths = repeats_info['rep_lengths'] # All possible length cutoffs
    input_seq = str(seq_record.seq).upper()
    input_seq_length = len(input_seq)
    for length_cutoff in repeat_lengths:
        fallback = length_cutoff - 1
        sub_start = 0  # substring start
        sub_stop = sub_start + repeat_lengths[-1]  # substring stop
        while sub_stop <= input_seq_length:
            sub_stop = sub_start + length_cutoff
            sub_seq = input_seq[sub_start:sub_stop]
            if sub_seq in repeats_info:
                match = True
                motif_length = repeats_info[sub_seq]['motif_length']
                offset = length_cutoff % motif_length
                repeat_seq = input_seq[sub_start+offset:sub_start+offset+motif_length]
                i = 0
                while match:
                    j = sub_stop
                    if sub_stop >= input_seq_length:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, repeats_info[sub_seq]['class'], match_length, repeats_info[sub_seq]['strand'], num_units, sub_seq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
                    elif input_seq[j] == repeat_seq[i]:
                        sub_stop += 1
                        i += 1
                        if i >= motif_length:
                            i = 0
                    else:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, repeats_info[sub_seq]['class'], match_length, repeats_info[sub_seq]['strand'], num_units, sub_seq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1
