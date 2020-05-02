#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
from itertools import product
from Bio import SeqIO
from tqdm import tqdm
import gzip, os
from os import remove as del_file
import multiprocessing as multi

from utils import rev_comp, rawcharCount, getGC, get_targetids



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
    rep_minlengths = dict()
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
            rep_minlengths[m] = length_cutoff
        repeats_out['fallback'] = motif_fallback
        repeats_out['rep_lengths'] = [length_cutoff]
        repeats_out['rep_minlengths'] = rep_minlengths

    elif unit_cutoff is not None:
        rep_minlengths = dict()
        for motif_length in unit_cutoff:
            rep_minlengths[motif_length] = motif_length*unit_cutoff[motif_length]
        min_length_cutoff = min(list(rep_minlengths.values()))
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            # try:
            #     motif = motif*unit_cutoff[motif_length]
            # except KeyError:
            #     motif = motif*unit_cutoff[0]
            motif = expand_repeat(motif, min_length_cutoff)
            repeat_lengths.add(len(motif))
            motif_fallback[motif_length] = len(motif) - 1
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        # repeat_lengths = sorted(list(repeat_lengths))
        repeats_out['rep_lengths'] = list(repeat_lengths)
        repeats_out['fallback'] = motif_fallback
        repeats_out['rep_minlengths'] = rep_minlengths

    return repeats_out



def get_ssrs(seq_record, repeats_info, out):
    """Native function that identifies repeats in fasta files."""
    if type(out) == str:
        out_file = open(out, 'w')
    else:
        out_file = out
    repeat_lengths = repeats_info['rep_lengths'] # All possible length cutoffs
    input_seq = str(seq_record.seq).upper()
    input_seq_length = len(input_seq)
    rep_minlengths = repeats_info['rep_minlengths']
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
                        if match_length >= rep_minlengths[motif_length]:
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
                        if match_length >= rep_minlengths[motif_length]:
                            print(seq_record.id, sub_start, sub_stop, repeats_info[sub_seq]['class'], match_length, repeats_info[sub_seq]['strand'], num_units, sub_seq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1
    if type(out) == str:
        out_file.close()
    

def fasta_ssrs(args, repeats_info):
    
    if args.input.endswith('gz'):
        handle = gzip.open(args.input, 'rt')
    else:
        handle = open(args.input, 'r')

    seq_nucleotide_info = dict()
    num_records = rawcharCount(args.input, '>')
    records = SeqIO.parse(handle, 'fasta')
    target_ids = get_targetids(args.filter_seq_ids, args.target_seq_ids)
    
    if args.threads > 1:
        i = 0
        pool = multi.Pool(processes=args.threads)
        for record in records:
            out_name = './temp_%s.tsv' %(i)
            i += 1
            if (args.info or args.analyse)==True:
                for a in record.seq.upper():
                    try: seq_nucleotide_info[a] += 1
                    except KeyError: seq_nucleotide_info[a] = 1
            if  args.min_seq_length <= len(record.seq) <= args.max_seq_length and record.id in target_ids:
                pool.apply_async(get_ssrs, (record, repeats_info, out_name,))
    
        pool.close() 
        pool.join()

        # Concat all the output files into one.
        temp_outs = tqdm(range(num_records), total=num_records)
        for o in temp_outs:
            name = './temp_%s.tsv' %(o)
            temp_outs.set_description("Concatenating file: %d " %(o))
            with open(name, 'r') as fh:
                for line in fh:
                    print(line.strip(), file=args.output)
            del_file(name)
    
    elif args.threads == 1:
        records = tqdm(records, total=num_records)
        for record in records:
            records.set_description("Processing %s" %(record.id))
            if (args.info or args.analyse)==True:
                for a in record.seq.upper():
                    try: seq_nucleotide_info[a] += 1
                    except KeyError: seq_nucleotide_info[a] = 1
            if  args.min_seq_length <= len(record.seq) <= args.max_seq_length and record.id in target_ids:
                get_ssrs(record, repeats_info, args.output)

    if (args.info or args.analyse)==True:
        line = "#File_name: %s\n#Total_sequences: %d\n#Total_bases: %d\n#GC: %f"\
        %(os.path.basename(args.input), num_records, sum(seq_nucleotide_info.values()),\
        round(getGC(seq_nucleotide_info), 2))
        print(line, file=args.output)
    
    args.output.close()