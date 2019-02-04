#! /usr/bin/env python

# pylint: disable=C0103, C0301

from __future__ import print_function, division
import sys
from os.path import splitext
import argparse
from tqdm import tqdm
from Bio import SeqIO
from collections import Counter, defaultdict
import gzip
from itertools import islice
from datetime import datetime

if sys.version_info.major == 2:
    from utils import generate_repeats, get_ssrs, get_ssrs_fastq, build_rep_set, univset, rawcharCount, dotDict
    from analyse import analyse
    from annotation import annotate
elif sys.version_info.major == 3:
    from utils import generate_repeats, get_ssrs, get_ssrs_fastq, build_rep_set, univset, rawcharCount, dotDict
    from analyse import analyse
    from annotation import annotate

inf = float('inf')

def getArgs():
    """
    Parses command line arguments and returns them to the caller
    """
    __version__ = 'v0.3.2'
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input file in FASTA format')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--format', metavar='<STR>', default='fasta', help='Input file format. Default: fasta, Permissible: fasta, fastq')
    # optional.add_argument('--fastq-report', type=argparse.FileType('w'), metavar='<FILE>', default=sys.stdout, help='Fastq report file name. Only usable if input file format is fastq. Default: Input file name + _perfFastqReport.tsv')
    optional.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='<FILE>', default=sys.stdout, help='Output file name. Default: Input file name + _perf.tsv')
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')
    optional.add_argument('-rep', '--repeats', type=argparse.FileType('r'), metavar='<FILE>', help='File with list of repeats (Not allowed with -m and/or -M)')
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-s', '--min-seq-length', type=int, metavar = '<INT>', default=0, help='Minimum size of sequence length for consideration (in bp)')
    optional.add_argument('-S', '--max-seq-length', type=float, metavar='<FLOAT>', default=inf, help='Maximum size of sequence length for consideration (in bp)')
    optional.add_argument('-g', '--annotate', metavar='<FILE>', help='Genic annotation input file for annotation, Both GFF and GTF can be processed. Use --anno-format to specify format.')
    optional.add_argument('--anno-format', type=str, default='GFF', help='Format of genic annotation file. Valid inputs: GFF, GTF. Default: GFF')
    optional.add_argument("--gene-attribute", metavar="<STR>", default="gene", type=str, help='Attribute key for geneId')
    optional.add_argument('--up-promoter', metavar="<INT>", type=int, default=1000, help='Upstream distance(bp) from TSS to be considered as promoter region. Default 1000')
    optional.add_argument('--down-promoter', metavar="<INT>", type=int, default=1000, help='Downstream distance(bp) from TSS to be considered as promoter region. Default 1000')    
    optional.add_argument('--version', action='version', version='PERF ' + __version__)
    cutoff_group = optional.add_mutually_exclusive_group()
    cutoff_group.add_argument('-l', '--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    cutoff_group.add_argument('-u', '--min-units', metavar='INT or FILE', help="Minimum number of repeating units to be considered. Can be an integer or a file specifying cutoffs for different motif sizes.")
    seqid_group = optional.add_mutually_exclusive_group()
    seqid_group.add_argument('-f', '--filter-seq-ids', metavar='<FILE>')
    seqid_group.add_argument('-F', '--target-seq-ids', metavar='<FILE>')

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
    # if args.format == 'fastq' and args.fastq_report.name == "<stdout>":
    #     args.fastq_report = open(splitext(args.input)[0] + '_perfFastqReport.tsv', 'w')

    if args.anno_format == 'GTF' and args.gene_attribute == 'gene':
        args.gene_attribute = 'gene_id'

    return args

def process_fastq(handle, min_seq_length, max_seq_length, repeats_info, repeat_set, out_file):
    n = 0
    b = 0
    reads_with_repeats = 0
    total_repeats = 0
    readlen_range = Counter()
    repFastq_info = {}
    startTime = datetime.now()
    while 1:
        lines_gen = list(islice(handle, 4))
        if len(lines_gen) == 0:
            break
        read_id = lines_gen[0].strip()
        read_seq = lines_gen[1].strip()
        n += 1 # total reads
        b += len(read_seq) # total bases
        readlen_range.update([len(read_seq)])
        # read length should be greater than minimum repeat length
        if n%50000 == 0:
            timeDiff = datetime.now() - startTime
            print('Processed reads: %d | Time elapsed: %s | Rate: %d iters/s\r' %(n, timeDiff, n/(timeDiff.seconds+1)), end = '')
            sys.stdout.flush()
        record = dotDict({'id': read_id, 'seq': read_seq})
        if  min_seq_length <= len(record.seq) <= max_seq_length:
            rep_identified = get_ssrs_fastq(record, repeats_info, repeat_set, out_file)
            for rep in rep_identified:
                try:
                    repFastq_info[rep]['lengths'].update(rep_identified[rep])
                    repFastq_info[rep]['instances'] += len(rep_identified[rep])
                    repFastq_info[rep]['reads'] += 1
                    repFastq_info[rep]['bases'] += sum(rep_identified[rep])
                except KeyError:
                    repFastq_info[rep] = {'lengths': Counter(rep_identified[rep]), 
                                            'instances': len(rep_identified[rep]),
                                            'reads': 1,
                                            'bases': sum(rep_identified[rep])}
                reads_with_repeats += 1
                total_repeats += len(rep_identified[rep])

    return {'info': {'total_reads': n, 'total_bases': b, 'reads_with_repeats': reads_with_repeats, 
            'total_repeats': total_repeats, 'readlen_range': readlen_range}, 'repInfo': repFastq_info}

def get_targetids(filter_seq_ids, target_seq_ids):
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

def getSSRNative(args):
    """
    Identifies microsatellites using native string matching.
    As the entire sequence is scanned in a single iteration, the speed is vastly improved
    """

    length_cutoff = args.min_length
    repeat_file = args.repeats
    seq_file = args.input
    out_file = args.output
    # fastq_report = args.fastq_report
    repeats_info = build_rep_set(repeat_file, length_cutoff=length_cutoff)
    repeat_set = set(repeats_info.keys())
    min_seq_length = args.min_seq_length
    max_seq_length = args.max_seq_length
    target_ids = get_targetids(args.filter_seq_ids, args.target_seq_ids)
    inFormat = args.format
    print('Using length cutoff of %d' % (length_cutoff), file=sys.stderr)

    if seq_file.endswith('gz'):
        handle = gzip.open(seq_file, 'rt')
    else:
        handle = open(seq_file, 'r')
    if inFormat == 'fasta':
        num_records = rawcharCount(seq_file, '>')
        records = SeqIO.parse(handle, 'fasta')
        records = tqdm(records, total=num_records)
        for record in records:
            records.set_description("Processing %s" %(record.id))
            if  min_seq_length <= len(record.seq) <= max_seq_length and record.id in target_ids:
                get_ssrs(record, repeats_info, repeat_set, out_file)
    elif inFormat == 'fastq':
        fastq_out = process_fastq(handle, min_seq_length, max_seq_length, repeats_info, repeat_set, out_file)
        repFastq_info = fastq_out['repInfo']
        total_repeats = fastq_out['info']['total_repeats']
        reads_with_repeats = fastq_out['info']['reads_with_repeats']
        n = fastq_out['info']['total_reads']
        b = fastq_out['info']['total_bases']
        readlen_range = fastq_out['info']['readlen_range']
        print('#Total reads: %d\n#Total bases: %d\n#Total repeat instances: %d\n#Total reads with repeats: %d\n#Total repeats per million reads: %d' \
        %(n, b, total_repeats, reads_with_repeats, round((total_repeats/n)*1000000, 3)), file=out_file)
        print('#Read lenghth distribution: ', end="", file=out_file)
        print(readlen_range.most_common(), file=out_file)
        print('repeatClass', 'reads', 'instances', 'bases', 'reads_norm', 'instances_norm', 'bases_norm', 'length_distribution', sep='\t', file=out_file)
        for rep in sorted(repFastq_info, key= lambda k: (len(k), k)):
            try:
                print(rep, repFastq_info[rep]['reads'], repFastq_info[rep]['instances'], repFastq_info[rep]['bases'],
                    round((repFastq_info[rep]['reads']/n)*1000000, 3), round((repFastq_info[rep]['instances']/n)*1000000, 3), 
                    round((repFastq_info[rep]['bases']/b)*1000000, 3), '-'.join([':'.join([str(y) for y in x]) for x in sorted(repFastq_info[rep]['lengths'].items())]),
                    sep='\t', file=out_file)
            except TypeError:
                print(rep, '\t'.join([str(0) for i in range(6)]), '-', sep="\t", file=out_file)
        print('')
    out_file.close()


def getSSR_units(args, unit_cutoff):
    """
    Identifies microsatellites using native string matching.
    The repeat length cutoffs vary for different motif sizes.
    """

    repeat_file = args.repeats
    seq_file = args.input
    inFormat = args.format
    out_file = args.output
    # fastq_report = args.fastq_report
    repeats_info = build_rep_set(repeat_file, unit_cutoff=unit_cutoff)
    repeat_set = set(repeats_info.keys())
    min_seq_length = args.min_seq_length
    max_seq_length = args.max_seq_length
    target_ids = get_targetids(args.filter_seq_ids, args.target_seq_ids)
    num_records = rawcharCount(seq_file, '>')
    inFormat = args.format

    print('Using unit cutoff of ', unit_cutoff, file=sys.stderr)
    if seq_file.endswith('gz'):
        handle = gzip.open(seq_file, 'rt')
    else:
        handle = open(seq_file, 'r')
    if inFormat == 'fasta':
        records = SeqIO.parse(handle, 'fasta')
        records = tqdm(records, total=num_records)
        for record in records:
            records.set_description("Processing %s" %(record.id))
            if  min_seq_length <= len(record.seq) <= max_seq_length and record.id in target_ids:
                get_ssrs(record, repeats_info, repeat_set, out_file)
    elif inFormat == 'fastq':
        fastq_out = process_fastq(handle, min_seq_length, max_seq_length, repeats_info, repeat_set, out_file)
        repFastq_info = fastq_out['repInfo']
        total_repeats = fastq_out['info']['total_repeats']
        reads_with_repeats = fastq_out['info']['reads_with_repeats']
        n = fastq_out['info']['total_reads']
        b = fastq_out['info']['total_bases']
        readlen_range = fastq_out['info']['readlen_range']
        print('#Total reads: %d\n#Total bases: %d\n#Total repeat instances: %d\n#Total reads with repeats: %d\n#Total repeats per million reads: %d' \
        %(n, b, total_repeats, reads_with_repeats, round((total_repeats/n)*1000000, 3)), file=out_file)
        print('#Read lenghth distribution: ', end="", file=out_file)
        print(readlen_range.most_common(), file=out_file)
        print('repeatClass', 'reads', 'instances', 'bases', 'reads_norm', 'instances_norm', 'bases_norm', 'length_distribution', sep='\t', file=out_file)
        for rep in sorted(repFastq_info, key= lambda k: (len(k), k)):
            try:
                print(rep, repFastq_info[rep]['reads'], repFastq_info[rep]['instances'], repFastq_info[rep]['bases'],
                    round((repFastq_info[rep]['reads']/n)*1000000, 3), round((repFastq_info[rep]['instances']/n)*1000000, 3), 
                    round((repFastq_info[rep]['bases']/b)*1000000, 3), '-'.join([':'.join([str(y) for y in x]) for x in sorted(repFastq_info[rep]['lengths'].items())]),
                    sep='\t', file=out_file)
            except TypeError:
                print(rep, '\t'.join([str(0) for i in range(6)]), '-', sep="\t", file=out_file)
        print('')
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
           
    if args.annotate is not None:
        annotate(args)

    # Specifies to generate a HTML report
    if args.analyse:
        analyse(args)

if __name__ == '__main__':
    main()
