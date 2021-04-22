#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
from itertools import product
from Bio import SeqIO
from tqdm import tqdm
import sys, gzip, os
from os import remove as del_file
import multiprocessing as multi

if sys.version_info.major == 2:
    from utils import rev_comp, rawcharCount, getGC, get_targetids
    from analyse import analyse_fasta
    from annotation import annotate
elif sys.version_info.major == 3:
    from .utils import rev_comp, rawcharCount, getGC, get_targetids
    from .analyse import analyse_fasta
    from .annotation import annotate

def num_factors(num):
    factors = []
    for i in range(1,num):
        if num%i == 0: factors.append(i)
    return factors

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


def generate_repeats(sizes, atomic):
    """Generates all possible motifs for repeats in a given length range"""
    generated_repeats = []
    alphabet = ['A', 'C', 'G', 'T']
    expanded_set = set()
    repeat_set = set()
    sizes.sort()
    min_size = sizes[0]
    max_size = sizes[-1]
    non_atomic_repeats = dict()
    for s in range(1, max_size):
        if s not in sizes:
            non_atomic_repeats[s] = set()
            if atomic:
                for combination in product(alphabet, repeat=s):
                    repeat = ''.join(combination)
                    expanded = expand_repeat(repeat, max_size)
                    non_atomic_repeats[s].add(expanded)
    for i in sizes:
        factors = num_factors(i)
        for combination in product(alphabet, repeat=i):
            repeat = ''.join(combination)
            repeat_revcomp = rev_comp(repeat)
            expanded = expand_repeat(repeat, max_size)
            atomic_check = False
            if atomic:
                for factor in factors:
                    if factor not in sizes and expanded in non_atomic_repeats[factor]:
                        atomic_check = True
            if expanded in expanded_set:
                continue
            elif atomic and atomic_check:
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
    motif_cutoff = dict()
    repeat_lengths = set()
    if length_cutoff is not None:
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            motif = expand_repeat(motif, length_cutoff)
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        repeats_out['cutoff'] = [length_cutoff]

    elif unit_cutoff is not None:
        cutoffs = set()
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            motif = motif*unit_cutoff[motif_length]
            cutoffs.add(len(motif))
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        repeats_out['cutoff'] = sorted(list(cutoffs))

    return repeats_out



def get_ssrs(seq_record, repeats_info, out):
    """Native function that identifies repeats in fasta files."""
    if type(out) == str:
        out_file = open(out, 'w')
    else:
        out_file = out
    length_cutoffs = repeats_info['cutoff']
    input_seq = str(seq_record.seq).upper()
    input_seq_length = len(input_seq)
    for length_cutoff in length_cutoffs:
        fallback = length_cutoff - 1
        sub_start = 0  # substring start
        sub_stop = sub_start + length_cutoff  # substring stop
        while sub_stop <= input_seq_length:
            sub_stop = sub_start + length_cutoff
            sub_seq = input_seq[sub_start:sub_stop]
            if sub_seq in repeats_info:
                match = True
                repeat_data = repeats_info[sub_seq]
                motif_length = repeat_data['motif_length']
                rep_class = repeat_data['class']
                strand = repeat_data['strand']
                offset = length_cutoff % motif_length
                repeat_seq = input_seq[sub_start+offset:sub_start+offset+motif_length]
                i = 0
                while match:
                    j = sub_stop
                    if sub_stop >= input_seq_length:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, rep_class, match_length, strand, num_units, sub_seq[:motif_length], sep="\t", file=out_file)
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
                        print(seq_record.id, sub_start, sub_stop, rep_class, match_length, strand, num_units, sub_seq[:motif_length], sep="\t", file=out_file)
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
           
    if args.annotate is not None:
        annotate(args)

    # Specifies to generate a HTML report
    if args.analyse:
        analyse_fasta(args)
    