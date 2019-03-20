#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
import sys
from itertools import product, takewhile, repeat, islice
from tqdm import tqdm
import gzip
from collections import Counter, defaultdict
from datetime import datetime


def getGC(basesCounter):
    totalBases = sum(basesCounter.values())
    try:
        GC = (float(basesCounter["G"] + basesCounter["C"])/(totalBases-basesCounter["N"]))*100
    except KeyError:
        GC = (float(basesCounter["G"] + basesCounter["C"])/totalBases)*100
    return GC


def rev_comp(string):
    """Outputs reverse complement of a nucleotide sequence"""
    if sys.version_info.major == 2:
        import string as st
        complement = string.translate(st.maketrans('ACGT', 'TGCA'))
    else:
        complement = string.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]


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


def rawcharCount(filename, char):
    if filename.endswith('gz'):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(char.encode('ASCII')) for buf in bufgen if buf )


def generate_repeats(min_size, max_size):
    """Generates all possible motifs for repeats in a given length range"""
    generated_repeats = []
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
                        generated_repeats.append('\t'.join([cycle, repeat, str(len(cycle)), strand]))
    return generated_repeats


def build_rep_set(repeat_file, **kwargs):
    """
        Outputs the repeats info dictionary used by the get_ssrs function.
        Takes list of lines from repeats file(output by generate_repeats function) as an input.
        Works either by "length_cutoff=" or by "unit_cutoff=" arguments.
    """
    repeats_out = dict()
    motif_fallback = dict()
    repeat_lengths = []
    length_cutoff = kwargs.get('length_cutoff', None)
    unit_cutoff = kwargs.get('unit_cutoff', None)
    if length_cutoff is not None:
        motif_lengths = []
        for line in repeat_file:
            motif_dict = dict()
            L = line.strip().split('\t')
            motif = L[0]
            motif_length = int(L[2])
            if motif_length not in motif_lengths:
                motif_lengths.append(motif_length)
            i = 0
            input_seq_length = motif_length
            while i < motif_length and input_seq_length < length_cutoff:
                motif += motif[i]
                i += 1
                if i >= motif_length:
                    i = 0
                input_seq_length = len(motif)
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
            if len(motif) not in repeat_lengths:
                repeat_lengths.append(len(motif))
            motif_fallback[motif_length] = len(motif) - 1
            motif_dict['class'] = L[1]
            motif_dict['motif_length'] = motif_length
            motif_dict['strand'] = L[3]
            repeats_out[motif] = motif_dict
        repeat_lengths.sort()
        repeats_out['rep_lengths'] = repeat_lengths
        repeats_out['fallback'] = motif_fallback

    return repeats_out


def get_targetids(filter_seq_ids, target_seq_ids):
    """
        The function returns the set of desired sequence ids 
        across which repeats will be identified.
    """
    target_ids = univset()
    if filter_seq_ids:
        target_ids = univset()
        filter_ids = []
        with open(filter_seq_ids) as fh:
            for line in fh:
                line = line.strip()
                line = line.lstrip('>')
                filter_ids.append(line)
        target_ids = target_ids - set(filter_ids)
    
    elif target_seq_ids:
        target_ids = []
        with open(target_seq_ids) as fh:
            for line in fh:
                line = line.strip()
                line = line.lstrip('>')
                target_ids.append(line)
        target_ids = set(target_ids)

    return target_ids


def get_ssrs(seq_record, repeats_info, repeats, out_file):
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
            subseq = input_seq[sub_start:sub_stop]
            if subseq in repeats:
                match = True
                motif_length = repeats_info[subseq]['motif_length']
                offset = length_cutoff % motif_length
                repeat_seq = input_seq[sub_start+offset:sub_start+offset+motif_length]
                i = 0
                while match:
                    j = sub_stop
                    if sub_stop >= input_seq_length:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, repeats_info[subseq]['class'], match_length, repeats_info[subseq]['strand'], num_units, subseq[:motif_length], sep="\t", file=out_file)
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
                        print(seq_record.id, sub_start, sub_stop, repeats_info[subseq]['class'], match_length, repeats_info[subseq]['strand'], num_units, subseq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1


def process_fastq(handle, min_seq_length, max_seq_length, repeats_info, repeat_set):
    """Processes fastq files and identifies repeats."""
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
        readlen_freq.update([len(read_seq)]) #updatind the read length frequencies
        if n%50000 == 0:
            time_diff = datetime.now() - start_time
            print('Processed reads: %d | Time elapsed: %s | Rate: %d iters/s\r' %(n, time_diff, n/(time_diff.seconds+1)), end = '')
            sys.stdout.flush()
        record = dotDict({'id': read_id, 'seq': read_seq})
        # read length should be greater than minimum repeat length
        if  min_seq_length <= len(record.seq) <= max_seq_length:
            rep_identified = get_ssrs_fastq(record, repeats_info, repeat_set)
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


def get_ssrs_fastq(seq_record, repeats_info, repeats):
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
            subseq = input_seq[sub_start:sub_stop]
            if subseq in repeats:
                match = True
                motif_length = repeats_info[subseq]['motif_length']
                offset = length_cutoff % motif_length
                repeat_seq = input_seq[sub_start+offset:sub_start+offset+motif_length]
                i = 0
                while match:
                    j = sub_stop
                    if sub_stop >= input_seq_length:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        rep_identified[repeats_info[subseq]['class']].append(match_length)
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
                        rep_identified[repeats_info[subseq]['class']].append(match_length)
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
    print('#Read_length_distribution: ', readlen_freq.most_common(), end="", file=out_file)
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



class univset(object):
    def __init__(self):
        self._diff = set()
 
    def __sub__(self, other):
        S = univset()
        if type(other) == set:
            S._diff = self._diff | other
            return S
        else:
            S._diff = self._diff | other._diff
            return S
 
    def __rsub__(self, other):
        return other &amp; self._diff
 
    def __contains__(self, obj):
        return not obj in self._diff
 
    def __and__(self, other):
        return other - self._diff
 
    def __rand__(self, other):
        return other - self._diff
 
    def __repr__(self):
        if self._diff == set():
            return "ANY"
        else:
            return "ANY - %s"%self._diff
 
    def __or__(self, other):
        S = univset()
        S._diff = self._diff - other
        return S
 
    def __xor__(self, other):
        return (self - other) | (other - self)
 
    def add(self, elem):
        if elem in self._diff:
            self._diff.remove(elem)
 
    def update(self, elem):
        self._diff = self._diff - other
 
    def __ror__(self, other):
        return self.__or__(other)
 
    def union(self, other):
        return self.__or__(other)
 
    def difference(self, other):
        return self.__sub__(other)
 
    def intersection(self, other):
        return self.__and__(other)
 
    def symmetric_difference(self, other):
        return self.__xor__(other)
 
    def __lt__(self, other):
        return self.issubset(other)
 
    def __eq__(self, other):
        if type(other) == set:
            return False
        try:
            return self._diff == other._diff
        except AttributeError:
            return False
 
    def __ne__(self, other):
        return not self.__eq__(other)
 
    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)
 
    def __gt__(self, other):
        return self.issuperset(other)
 
    def __gt__(self, other):
        return self.issuperset(other) or self == other


class dotDict(dict):
    """
    Example:
    m = dotDict({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, *args, **kwargs):
        super(dotDict, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(dotDict, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(dotDict, self).__delitem__(key)
        del self.__dict__[key]