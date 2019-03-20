#! /usr/bin/env python

from collections import defaultdict

defaultInfo = { 'info': { 'readsInfo': {}, 'repInfo': defaultdict() } }

def analyse_fastq(args):
    out_file = args.output
    with open(out_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                defaultInfo['info']['readsInfo'][fields[0]] = fields[1]
            else:
                fields = line.strip().split('\t')
                
                repClass = line[0]
                repReads = line[1]
                repFreq = line[2]
                repBases = line[3]
                repReadsNorm = line[4]
                repFreqNorm = line[5]
                repBasesNorm = line[6]
                repLengthDist = {}
                
                for i in line[7].split('-'):
                    entry = i.split(':')
                    repLengthDist[entry[0]] = entry[1]

                defaultInfo['info']['repInfo'][repClass] = {}
                defaultInfo['info']['repInfo'][repClass]['reads'] = repReads
                defaultInfo['info']['repInfo'][repClass]['freq'] = repFreq
                defaultInfo['info']['repInfo'][repClass]['bases'] = repBases
                defaultInfo['info']['repInfo'][repClass]['read_norm'] = repReadsNorm
                defaultInfo['info']['repInfo'][repClass]['freq_norm'] = repFreqNorm
                defaultInfo['info']['repInfo'][repClass]['bases_norm'] = repBasesNorm
                defaultInfo['info']['repInfo'][repClass]['length_distribution'] = repLengthDist

