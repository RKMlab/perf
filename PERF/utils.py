#! /usr/bin/env python
# pylint: disable=C0111, C0301

from __future__ import print_function, division
import sys, gzip
from itertools import takewhile, repeat, islice
from tqdm import tqdm
from collections import Counter, defaultdict


kmers = { 
    1: 'Monomer', 2: 'Dimer', 3: 'Trimer', 4: 'Tetramer', 5: 'Pentamer',
    6: 'Hexamer', 7: 'Heptamer', 8: 'Octamer', 9: 'Nonamer', 10: 'Decamer',
    11: 'Undecamer', 12: 'Dodecamer', 13: 'Tridecamer', 14: 'Tetradecamer', 15: 'Pentadecamer',
    16: 'Hexadecamer', 17: 'Heptadecamer', 18: 'Octadecamer', 19: 'Nonadecamer', 20: 'Icosamer',
    21: 'Uncosamer', 22: 'Docosamer', 23: 'Tricosamer', 24: 'Tetracosamer', 25: 'Pentacosamer',
    26: 'Hexacosamer', 27: 'Heptacosamer', 28: 'Octacosamer', 29: 'Nonacosamer', 30: 'Triacontamer',
    31: 'Untriacontamer', 32: 'Dotriacontamer', 33: 'Tritriacontamer', 34: 'Tetratriacontamer', 35: 'Pentatriacontamer',
    36: 'Hexatriacontamer', 37: 'Heptatriacontamer', 38: 'Octatriacontamer', 39: 'Nonatriacontamer', 40: 'Tetracontamer',
    41: 'Untetracontamer', 42: 'Dotetracontamer', 43: 'Tritetracontamer', 44: 'Tetratetracontamer', 45: 'Pentatetracontamer',
    46: 'Hexatetracontamer', 47: 'Heptatetracontamer', 48: 'Octatetracontamer', 49: 'Nonatetracontamer', 50: 'Pentacontamer',
}


def get_cycles(string):
    cycles = set()
    for i in range(len(string)):
        cycles.add(string[i:] + string[:i])
    cycles = sorted(list(cycles))
    return cycles


def build_cycVariations(string):
    cycles = get_cycles(string)
    rev_cycles = get_cycles(rev_comp(string))
    for r in rev_cycles:
        if r not in cycles: cycles.append(r)
    return cycles


def getGC(basesCounter):
    totalBases = sum(basesCounter.values())
    try:
        GC = (float(basesCounter['G'] + basesCounter['C'])/(totalBases-basesCounter['N']))*100
    except KeyError:
        GC = (float(basesCounter['G'] + basesCounter['C'])/totalBases)*100
    return GC


def rev_comp(string):
    """Outputs reverse complement of a nucleotide sequence"""
    if sys.version_info.major == 2:
        import string as st
        complement = string.translate(st.maketrans('ACGT', 'TGCA'))
    else:
        complement = string.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]


def rawcharCount(filename, char):
    if filename.endswith('gz'):
        f = gzip.open(filename, 'rb')
    else:
        f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(char.encode('ASCII')) for buf in bufgen if buf )


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
                line = line.split(' ')[0]
                filter_ids.append(line)
        target_ids = target_ids - set(filter_ids)
    
    elif target_seq_ids:
        target_ids = []
        with open(target_seq_ids) as fh:
            for line in fh:
                line = line.strip()
                line = line.lstrip('>')
                line = line.split(' ')[0]
                target_ids.append(line)
        target_ids = set(target_ids)

    return target_ids


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