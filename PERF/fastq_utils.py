#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
from datetime import datetime
from itertools import islice
from collections import Counter, defaultdict
import sys, gzip
import multiprocessing as multi

if sys.version_info.major == 2:
    from utils import dotDict, build_cycVariations
    from analyse import analyse_fastq
elif sys.version_info.major == 3:
    from .utils import dotDict, build_cycVariations
    from .analyse import analyse_fastq


def get_ssrs_fastq(seq_record, repeats_info):
    """Native function to identify repeats in fastq files"""

    repeats = defaultdict(list)
    num_repeats = 0
    read = 0
    
    length_cutoffs = repeats_info['cutoff'] # All possible length cutoffs
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
                num_repeats += 1
                read = 1
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
                        repeats[rep_class].append('%s-%d'%(sub_seq[:motif_length], match_length))
                        sub_start = sub_stop - fallback
                    elif input_seq[j] == repeat_seq[i]:
                        sub_stop += 1
                        i += 1
                        if i >= motif_length:
                            i = 0
                    else:
                        match = False
                        match_length = sub_stop - sub_start
                        repeats[rep_class].append('%s-%d'%(sub_seq[:motif_length], match_length))
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1
    return {'read': read, 'num_repeats': num_repeats, 'repeats': repeats}


def process_fastq(args, repeats_info):
    """Processes fastq files for identification of repeats."""

    print('\nProcessing fastq file...')

    if args.input.endswith('gz'):
        handle = gzip.open(args.input, 'rt')
    else:
        handle = open(args.input, 'r')

    n,b = [0,0]
    readlen_freq = Counter()
    start_time = datetime.now()
    
    total_repeats = 0
    reads_with_repeats = 0
    total_repeat_bases = 0
    fastq_repeats = dict()

    fastq_repeat_info = dict()
    repeat_class_info = dict()
    lenFreq_data = dict()
    
    while 1:
        lines_gen = list(islice(handle, 4))
        if len(lines_gen) == 0:
            for rep in fastq_repeats:
                cycles = build_cycVariations(rep)
                for instance in fastq_repeats[rep]:
                    f = fastq_repeats[rep][instance]
                    m = instance.split('-')[0]
                    l = int(instance.split('-')[1])
                    if l not in lenFreq_data[rep]:
                        lenFreq_data[rep][l] = [0]*len(cycles)
                    lenFreq_data[rep][l][cycles.index(m)] += f
                    repeat_class_info[rep]['lengths'].update([l]*f)
                    repeat_class_info[rep]['motifs'].update([m]*f)
                    repeat_class_info[rep]['instances'] += f
                    repeat_class_info[rep]['bases'] += l*f
                    total_repeat_bases += l*f
                fastq_repeats[rep] = Counter()
            time_diff = datetime.now() - start_time
            print('Processed reads: %d | Time elapsed: %s | Rate: %d iters/s\r' %(n, time_diff, n/(time_diff.seconds+1)), end = '')
            sys.stdout.flush()
            break            
        
        read_id = lines_gen[0].strip()
        read_seq = lines_gen[1].strip()
        read_len = len(read_seq)
        n += 1 # total reads
        b += len(read_seq) # total bases
        readlen_freq.update([read_len]) #updating the read length frequencies

        if n%50000 == 0:
            for rep in fastq_repeats:
                cycles = build_cycVariations(rep)
                for instance in fastq_repeats[rep]:
                    f = fastq_repeats[rep][instance]
                    m = instance.split('-')[0]
                    l = int(instance.split('-')[1])
                    if l not in lenFreq_data[rep]:
                        lenFreq_data[rep][l] = [0]*len(cycles)
                    lenFreq_data[rep][l][cycles.index(m)] += f
                    repeat_class_info[rep]['lengths'].update([l]*f)
                    repeat_class_info[rep]['motifs'].update([m]*f)
                    repeat_class_info[rep]['instances'] += f
                    repeat_class_info[rep]['bases'] += l*f
                    total_repeat_bases += l*f
                fastq_repeats[rep] = Counter()
            time_diff = datetime.now() - start_time
            print('Processed reads: %d | Time elapsed: %s | Rate: %d iters/s\r' %(n, time_diff, n/(time_diff.seconds+1)), end = '')
            sys.stdout.flush()
        
        record = dotDict({'id': read_id, 'seq': read_seq})
        # read length should be greater than minimum repeat length
        if  args.min_seq_length <= read_len <= args.max_seq_length:            
            rep_identified = get_ssrs_fastq(record, repeats_info)
            for rep in rep_identified['repeats']:
                try:
                    fastq_repeats[rep].update(rep_identified['repeats'][rep])
                    repeat_class_info[rep]['reads'] += 1
                except KeyError:
                    fastq_repeats[rep] = Counter(rep_identified['repeats'][rep])
                    repeat_class_info[rep] = {'lengths': Counter(), 'motifs': Counter(), 'instances': 0, 'reads': 1, 'bases': 0}
                    lenFreq_data[rep] = {}
            total_repeats += rep_identified['num_repeats']
            reads_with_repeats += rep_identified['read']
    
    fastq_repeat_info['totalRepReads'] = reads_with_repeats
    fastq_repeat_info['totalRepFreq'] = total_repeats
    fastq_repeat_info['totalRepBases'] = total_repeat_bases
    fastq_repeat_info['numRepClasses'] = len(lenFreq_data.keys())
    fastq_repeat_info['lenFrequency'] = lenFreq_data
    fastq_repeat_info['repClassInfo'] = repeat_class_info

    print('') #A line for proper printing of the output
    min_readlen = min(readlen_freq.keys())
    max_readlen = max(readlen_freq.keys())
    if min_readlen == max_readlen:
        readlen_range = '%d bp' %(min_readlen)
    else:
        readlen_range = '%d-%d (bp)' %(min_readlen, max_readlen)
    return {'info': { 'seqInfo': {'Total_reads': n, 'Total_bases': b, 'Readlen_freq': readlen_freq, 
            'Readlen_range': readlen_range}, 'repInfo': fastq_repeat_info}}


def ssr_fastq_output(fastq_out, out_file):
    """PERF OUTPUT for fastq files."""
    
    
    n = fastq_out['info']['seqInfo']['Total_reads']
    b = fastq_out['info']['seqInfo']['Total_bases']
    readlen_freq = fastq_out['info']['seqInfo']['Readlen_freq']
    readlen_range = fastq_out['info']['seqInfo']['Readlen_range']

    fastq_repeat_info = fastq_out['info']['repInfo']
    reads_with_repeats = fastq_repeat_info['totalRepReads']
    total_repeats = fastq_repeat_info['totalRepFreq']
    total_repeat_classes = fastq_repeat_info['numRepClasses']
    repeat_class_info = fastq_repeat_info['repClassInfo']
    
    print('#Total_reads: %d'%(n), file=out_file)
    print('#Total_bases: %d' %(b), file=out_file)
    print('#Total_repeat_instances: %d' %(total_repeats), file=out_file)
    print('#Total_reads_with_repeats: %d' %(reads_with_repeats), file=out_file)
    print('#Total_repeats_per_million_reads: %f' %(round((total_repeats/n)*1000000, 3)), file=out_file)
    print('#Read_length_distribution: ', readlen_freq.most_common(), file=out_file)

    print('repeatClass', 'reads', 'instances', 'bases', 'reads_per_million', 'instances_per_million',
        'bases_norm', 'length_distribution', 'motif_distribution', sep='\t', file=out_file)
    
    for rep in sorted(repeat_class_info, key= lambda k: (len(k), k)):
        rep_info = repeat_class_info[rep]
        print(  
                rep, int(rep_info['reads']), int(rep_info['instances']), int(rep_info['bases']),
                round((rep_info['reads']/n)*1000000, 3),
                round((rep_info['instances']/n)*1000000, 3), 
                round((rep_info['bases']/b)*1000000, 3),
                ';'.join(['-'.join([str(y) for y in x]) for x in sorted(rep_info['lengths'].items())]),
                ';'.join(['-'.join([str(y) for y in x]) for x in sorted(rep_info['motifs'].items())]),
                sep='\t', file=out_file  
            )


def fastq_ssrs(args, repeats_info):

    fastq_out = process_fastq(args, repeats_info)
    ssr_fastq_output(fastq_out, args.output)
    if args.analyse:
        analyse_fastq(args, fastq_out)
    args.output.close()