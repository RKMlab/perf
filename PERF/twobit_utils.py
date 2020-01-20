"""
twobitreader

Licensed under Perl Artistic License 2.0
No warranty is provided, express or implied
"""
from array import array
from bisect import bisect_right
from errno import ENOENT, EACCES
from os import R_OK, access
try:
    from os import strerror
except ImportError:
    strerror = lambda x: 'strerror not supported'
from bitarray import bitarray
ba = bitarray()
from os.path import exists
from itertools import izip
import logging
import textwrap
import sys

def true_long_type():
    """
    OS X uses an 8-byte long, so make sure L (long) is the right size
    and switch to I (int) if needed
    """
    for type_ in ['L', 'I']:
        test_array = array(type_, [0])
        long_size = test_array.itemsize
        if long_size == 4:
            return type_
    raise ImportError("Couldn't determine a valid 4-byte long type to use \
as equivalent to LONG")

LONG = true_long_type()

def byte_to_bases(x):
    """convert one byte to the four bases it encodes"""
    c = (x >> 4) & 0xf
    f = x & 0xf
    cc = (c >> 2) & 0x3
    cf = c & 0x3
    fc = (f >> 2) & 0x3
    ff = f & 0x3
    return map(bits_to_base, (cc, cf, fc, ff))

def bits_to_base(x):
    """convert integer representation of two bits to correct base"""
    if x is 0: return 'T'
    if x is 1: return 'C'
    if x is 2: return 'A'
    if x is 3: return 'G'

def base_to_bin(x):
    """
    provided for user convenience
    convert a nucleotide to its bit representation 
    """
    if x == 'T': return '00'
    if x == 'C': return '01'
    if x == 'A': return '10'
    if x == 'G': return '11'

def create_byte_table():
    """create BYTE_TABLE"""
    d = {}
    for x in xrange(2**8):
        d[x] = byte_to_bases(x)
    return d

def split16(x):
    """
    split a 16-bit number into integer representation
    of its course and fine parts in binary representation
    """
    c = (x >> 8) & 0xff
    f = x & 0xff
    return c, f

def create_twobyte_table():
    """create TWOBYTE_TABLE"""
    d = {}
    for x in xrange(2**16):
        c, f = split16(x)
        d[x] = byte_to_bases(c) + byte_to_bases(f)
    return d

BYTE_TABLE = create_byte_table()
TWOBYTE_TABLE = create_twobyte_table()

def longs_to_char_array(longs, first_base_offset, last_base_offset, array_size):
    """
    takes in a iterable of longs and converts them to bases in a char array
    returns a ctypes string buffer
    """
    longs_len = len(longs)

    # dna = ctypes.create_string_buffer(array_size)
    dna = array('c', 'N' * longs_len)
    # translate from 32-bit blocks to bytes
    # this method ensures correct endianess (byteswap as neeed)
    bytes = array('B')
    bytes.fromstring(longs.tostring())

    # first block
    first_block = ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(4)])
    i = 16 - first_base_offset
    if array_size < i: i = array_size
    dna[0:i] = array('c', first_block[first_base_offset:first_base_offset + i])
    if longs_len == 1: return dna
    # middle blocks (implicitly skipped if they don't exist)
    for byte in bytes[4:-4]:
        dna[i:i + 4] = array('c', BYTE_TABLE[byte])
        i += 4
    # last block
    last_block = array('c', ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(-4,0)]))
    dna[i:i + last_base_offset] = last_block[0:last_base_offset]
    print(dna)
    return dna

class TwoBitFile(dict):
    """
        python-level reader for .2bit files (i.e., from UCSC genome browser)
        (note: no writing support)

        TwoBitFile inherits from dict
        You may access sequences by name, e.g.
        >>> genome = TwoBitFile('hg18.2bit')
        >>> chr20 = genome['chr20']

        Sequences are returned as TwoBitSequence objects
        You may access intervals by slicing or using str() to dump the entire entry
        e.g.
        >>> chr20[100100:100200]
        'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggtat
        ttgcatgttaagtttttttcc'
        >>> whole_chr20 = str(chr20)

        Fair warning: dumping the entire chromosome requires a lot of memory

        See TwoBitSequence for more info
    """
    
    def __init__(self, foo):
        super(TwoBitFile, self).__init__()
        if not exists(foo):
            raise IOError(ENOENT, strerror(ENOENT), foo)
        if not access(foo, R_OK):
            raise IOError(EACCES, strerror(EACCES), foo)
        self._filename = foo
        self._file_handle = open(foo, 'rb')
        self._load_header()
        self._load_index()
        for name, offset in self._offset_dict.iteritems():
            self[name] = TwoBitSequence(self._file_handle, offset,
                                        self._byteswapped)
        return        
        
    def _load_header(self):
        file_handle = self._file_handle
        header = array(LONG)
        header.fromfile(file_handle, 4)
        # check signature -- must be 0x1A412743
        # if not, swap bytes
        byteswapped = False
        (signature, version, sequence_count, reserved) = header
        if not signature == 0x1A412743:
            byteswapped = True
            header.byteswap()
            (signature2, version, sequence_count, reserved) = header
            if not signature2 == 0x1A412743:
                raise TwoBitFileError('Signature in header should be 0x1A412743'
                                    + ', instead found 0x%X' % signature)
        if not version == 0: 
            raise TwoBitFileError('File version in header should be 0.')
        if not reserved == 0:
            raise TwoBitFileError('Reserved field in header should be 0.')
        self._byteswapped = byteswapped
        self._sequence_count = sequence_count
        
    def _load_index(self):
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        remaining = self._sequence_count
        sequence_offsets = []
        file_handle.seek(16)
        while True:
            if remaining == 0: break
            name_size = array('B')
            name_size.fromfile(file_handle, 1)
            if byteswapped: name_size.byteswap()
            name = array('c')
            if byteswapped: name.byteswap()
            name.fromfile(file_handle, name_size[0])
            offset = array(LONG)
            offset.fromfile(file_handle, 1)
            if byteswapped: offset.byteswap()
            sequence_offsets.append((name.tostring(), offset[0]))
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._offset_dict = dict(sequence_offsets)

    def sequence_sizes(self):
        """returns a dictionary with the sizes of each sequence"""
        d = {}
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        for name, offset in self._offset_dict.iteritems():
            file_handle.seek(offset)
            dna_size = array(LONG)
            dna_size.fromfile(file_handle, 1)
            if byteswapped: dna_size.byteswap()
            d[name] = dna_size[0]
        return d