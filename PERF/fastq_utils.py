#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
from datetime import datetime
from itertools import islice
from collections import Counter, defaultdict
import sys
import multiprocessing as multi

from utils import dotDict


def process_fastq(handle, min_seq_length, max_seq_length, repeats_info):
    """Processes fastq files for identification of repeats."""
    n,b = [0,0]
    total_repeats = 0
    reads_with_repeats = 0
    total_repeat_bases = 0
    fastq_repeat_info = dict()
    readlen_freq = Counter()
    start_time = datetime.now()
    while 1:
        lines_gen = list(islice(handle, 4))
        if len(lines_gen) == 0:
            break
        read_id = lines_gen[0].strip()
        read_seq = lines_gen[1].strip()
        n += 1 # total reads
        b += len(read_seq) # total bases
        readlen_freq.update([len(read_seq)]) #updating the read length frequencies
        if n%50000 == 0:
            time_diff = datetime.now() - start_time
            print('Processed reads: %d | Time elapsed: %s | Rate: %d iters/s\r' %(n, time_diff, n/(time_diff.seconds+1)), end = '')
            sys.stdout.flush()
        record = dotDict({'id': read_id, 'seq': read_seq})
        # read length should be greater than minimum repeat length
        if  min_seq_length <= len(record.seq) <= max_seq_length:
            rep_identified = get_ssrs_fastq(record, repeats_info)
            for rep in rep_identified:
                try:
                    fastq_repeat_info[rep]['lengths'].update(rep_identified[rep])
                    fastq_repeat_info[rep]['instances'] += len(rep_identified[rep])
                    fastq_repeat_info[rep]['reads'] += 1
                    fastq_repeat_info[rep]['bases'] += sum(rep_identified[rep])
                except KeyError:
                    fastq_repeat_info[rep] = {'lengths': Counter(rep_identified[rep]), 
                                            'instances': len(rep_identified[rep]),
                                            'reads': 1,
                                            'bases': sum(rep_identified[rep])}
                reads_with_repeats += 1
                total_repeats += len(rep_identified[rep])
    print('') #A line for proper printing of the output
    min_readlen = min(readlen_freq.keys())
    max_readlen = max(readlen_freq.keys())
    if min_readlen == max_readlen:
        readlen_range = '%d bp' %(min_readlen)
    else:
        readlen_range = '%d-%d (bp)' %(min_readlen, max_readlen)
    return {'info': { 'readInfo': {'total_reads': n, 'total_bases': b, 'reads_with_repeats': reads_with_repeats, 
            'total_repeats': total_repeats, 'total_repeat_bases': total_repeat_bases, 'readlen_freq': readlen_freq, 
            'readlen_range': readlen_range}, 'repInfo': fastq_repeat_info}}


def get_ssrs_fastq(seq_record, repeats_info):
    """Native function to identify repeats in fastq files"""

    rep_identified = defaultdict(list)
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
                        rep_identified[repeats_info[sub_seq]['class']].append(match_length)
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
                        rep_identified[repeats_info[sub_seq]['class']].append(match_length)
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1
    return rep_identified


def ssr_fastq_output(fastq_out, out_file):
    """PERF OUTPUT for fastq files."""
    rep_fastq_info = fastq_out['info']['repInfo']
    total_repeats = fastq_out['info']['readInfo']['total_repeats']
    reads_with_repeats = fastq_out['info']['readInfo']['reads_with_repeats']
    n = fastq_out['info']['readInfo']['total_reads']
    b = fastq_out['info']['readInfo']['total_bases']
    readlen_freq = fastq_out['info']['readInfo']['readlen_freq']
    total_repeat_classes = 0
    
    print('#Total_reads: %d\n#Total_bases: %d\n#Total_repeat_instances:\
        %d\n#Total_reads_with_repeats: %d\n#Total_repeats_per_million_reads: %d'\
    %(n, b, total_repeats, reads_with_repeats, round((total_repeats/n)*1000000, 3)),
    file=out_file)
    print('#Read_length_distribution: ', readlen_freq.most_common(), file=out_file)
    print('repeatClass', 'reads', 'instances', 'bases', 'reads_norm', 'instances_norm',
        'bases_norm', 'length_distribution', sep='\t', file=out_file)
    
    for rep in sorted(rep_fastq_info, key= lambda k: (len(k), k)):
        try:
            total_repeat_classes += 1
            print(rep, rep_fastq_info[rep]['reads'], rep_fastq_info[rep]['instances'], rep_fastq_info[rep]['bases'],
                round((rep_fastq_info[rep]['reads']/n)*1000000, 3), round((rep_fastq_info[rep]['instances']/n)*1000000, 3), 
                round((rep_fastq_info[rep]['bases']/b)*1000000, 3), '-'.join([':'.join([str(y) for y in x]) for x in sorted(rep_fastq_info[rep]['lengths'].items())]),
                sep='\t', file=out_file)
        except TypeError:
            print(rep, '\t'.join([str(0) for i in range(6)]), '-', sep="\t", file=out_file)

