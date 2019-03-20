#! /usr/bin/env python

from __future__ import print_function, division
import sys
import os
import json
from collections import Counter, defaultdict
from Bio import SeqIO
import gzip
import numpy as np


def writetoHTML(html_file, defaultInfo):
    html_handle = open(html_file, 'w')
    current_dir = os.path.dirname(__file__)
    with open(current_dir + '/lib/template.html') as report:
        for line in report:
            line = line.strip()
            print(line, file=html_handle)
            try:
                start_index = line.index("^^")
                stop_index = line.index("$$")
                if (line[start_index+2: stop_index] == 'defaultInfo'):
                    print(defaultInfo, file=html_handle)
                else:
                    file_path = current_dir + '/lib' + line[start_index+2: stop_index]
                    with open(file_path) as fh:
                        for subline in fh:
                            subline = subline.strip()
                            print(subline, file=html_handle)
            except ValueError:
                pass
    html_handle.close()
    print("HTML report successfully saved to " + html_file)


def get_parameters(args):
    runCommand = 'PERF' + ' '.join(sys.argv)


def analyse(args):
    outFile = open('./analyse_data.js', 'w')
    repeatsOutFile = args.output.name
    current_dir = os.path.dirname(__file__)
    html_report = os.path.splitext(repeatsOutFile)[0] + '.html'
    inFormat = args.format
    print("Generating HTML report. This may take a while..")
    

    all_repeat_classes = list(map(lambda x: x.split('\t')[1], args.repeats))
    temp = []
    for a in all_repeat_classes:
        if a not in temp:
            temp.append(a)
    all_repeat_classes = temp
    del temp

    inf = float('inf')
    defaultInfo = {}
    defaultInfo['info'] = {'genomeInfo': {}, 'repInfo': {}, 'plotInfo': {}}
    
    if args.annotate: #if annotation is on the data is taken from t
        repeatsOutFile = os.path.splitext(repeatsOutFile)[0] + '_annotation.tsv'
        promUp = args.up_promoter
        promDown = args.down_promoter
        defaultInfo['info']['annoInfo'] = {'promUp': promUp, 'promDown': promDown}
        repAnno = {}
        TSS_dist = {}
        annoKeyDict = {}
    
    totalRepBases = 0
    totalRepFreq = 0
    longestLengths = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
    mostUnits = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
    minLength = inf
    minUnits = inf

    plot_data = defaultdict(defaultdict)
    with open(repeatsOutFile, 'r') as repFile:
        for line in repFile:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                defaultInfo['info']['genomeInfo'][fields[0]] = fields[1]
            else:
                fields = line.split('\t')
                fields[1] = int(fields[1])
                fields[2] = int(fields[2])
                fields[4] = int(fields[4])
                fields[6] = int(fields[6])

                seq = fields[0]
                start = fields[1]
                end = fields[2]
                repClass = fields[3]
                repLength = fields[4]
                repOri = fields[5]
                repUnit = fields[6]
                actualRepeat = fields[7]
                if args.annotate:
                    if repClass not in repAnno:
                        repAnno[repClass] = {'EP': 0, 'GP': 0, 'GN': 0, 'IP': 0, 'DP': 0, 'EN': 0, 'IN': 0, 'DN': 0, 'UU': 0}
                        TSS_dist[repClass] = []
                    genicKey = fields[12]
                    promKey = fields[13]
                    try:
                        tssD = int(fields[-1])
                        if -5000 <= tssD <= 5000:
                            TSS_dist[repClass].append(tssD)
                    except:
                        pass
                    if genicKey == 'Intergenic':
                        genicKey = 'Distal Intergenic'
                    elif genicKey == '-':
                        genicKey = 'Unannotated'
                        promKey = 'Unannotated'
                    annoKey = genicKey[0]+promKey[0]
                    annoKeyDict[annoKey] = genicKey + '+' + promKey
                    repAnno[repClass][annoKey] += 1 

                totalRepBases += repLength
                totalRepFreq += 1
                if repLength in plot_data[repClass]:
                    plot_data[repClass][repLength] += 1
                else:
                    plot_data[repClass][repLength] = 1
                
                if minUnits > repUnit:
                    minUnits = repUnit
                if minLength > repLength:
                    minLength = repLength

                if longestLengths[-1][4] < repLength:
                    longestLengths[-1] = fields
                elif longestLengths[-1][4] == repLength:
                    if repClass < longestLengths[-1][3]:
                        longestLengths[-1] = fields
                longestLengths.sort(key=lambda x: x[4])
                longestLengths.reverse()
                if mostUnits[-1][6] < repUnit:
                    mostUnits[-1] = fields
                elif mostUnits[-1][6] == repUnit:
                    if repClass < longestLengths[-1][3]:
                        longestLengths[-1] = fields
                mostUnits.sort(key=lambda x: x[6])
                mostUnits.reverse()
    totalBases = int(defaultInfo['info']['genomeInfo']['Total_bases'])
    defaultInfo['info']['plotInfo'] = plot_data
    defaultInfo['info']['repInfo']['numRepClasses'] = len(plot_data.keys())
    defaultInfo['info']['repInfo']['totalRepBases'] = totalRepBases
    defaultInfo['info']['repInfo']['totalRepFreq'] = totalRepFreq
    defaultInfo['info']['repInfo']['percentGenomeCovered'] = str(round((totalRepBases/totalBases)*100, 2)) + "%"
    defaultInfo['info']['repInfo']['repDensityByFreq'] = round((totalRepFreq/totalBases)*1000000, 2)
    defaultInfo['info']['repInfo']['repDensityByBases'] = round((totalRepBases/totalBases)*1000000, 2)
    defaultInfo['info']['repInfo']['minLength'] = minLength
    defaultInfo['info']['repInfo']['minUnits'] = minUnits
    defaultInfo['info']['repInfo']['longestRepeats'] = []
    defaultInfo['info']['repInfo']['mostRepeatUnits'] = []
    defaultInfo['info']['repInfo']['allRepClasses'] = all_repeat_classes
    if args.annotate:
        for r in TSS_dist:
            hist_values = np.histogram(TSS_dist[r], bins=200, range=(-5000,5000))
            TSS_dist[r] = list(hist_values[0])
            TSS_dist[r] = list(map(lambda x: int(x), TSS_dist[r]))
        defaultInfo['info']['annoInfo']['TSS_histBinEdges'] = list(map(lambda x: int(x), hist_values[1]))
        defaultInfo['info']['annoInfo']['repAnno'] = repAnno
        defaultInfo['info']['annoInfo']['TSS_dist'] = TSS_dist
        defaultInfo['info']['annoInfo']['annoKeyObj'] = annoKeyDict
    for a in longestLengths:
        testDict = {'seq': a[0], 'start': a[1], 'end': a[2], 'repClass': a[3], 'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]}
        defaultInfo['info']['repInfo']['longestRepeats'].append(testDict)
    for a in mostUnits:
        testDict = {'seq': a[0], 'start': a[1], 'end': a[2], 'repClass': a[3], 'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]}
        defaultInfo['info']['repInfo']['mostRepeatUnits'].append(testDict)
    defaultInfo = 'const data =' + json.dumps(defaultInfo)
    print(defaultInfo, file=outFile)
    outFile.close()

def analyse_fastq(args, fastq_out):
    """Generates HTML report for fastq files."""
    fastq_out['info']['readInfo']['file_name'] = args.input.split('/')[-1]
    total_repeats = fastq_out['info']['readInfo']['total_repeats']
    reads_with_repeats = fastq_out['info']['readInfo']['reads_with_repeats']
    total_repeat_bases = fastq_out['info']['readInfo']['total_repeat_bases']
    all_repeat_classes = list(map(lambda x: x.split('\t')[1], args.repeats))
    total_repeat_classes = len(fastq_out['info']['readInfo'].keys())
    n = fastq_out['info']['readInfo']['total_reads']
    b = fastq_out['info']['readInfo']['total_bases']
    temp = []
    for a in all_repeat_classes:
        if a not in temp:
            temp.append(a)
    all_repeat_classes = temp
    del temp
    fastq_out['info']['readInfo']['all_rep_classes'] = all_repeat_classes
    fastq_out['info']['readInfo']['repeatsnorm'] = round((total_repeats/n)*1000000, 2)
    fastq_out['info']['readInfo']['reads_with_repeats_norm'] = str(round((reads_with_repeats/n)*100, 2)) + '%'
    fastq_out['info']['readInfo']['percent_repeat_bases'] = str(round((total_repeat_bases/b)*100, 2)) + '%'
    fastq_out['info']['readInfo']['total_repeat_classes'] = str(total_repeat_classes) + '/501'
    info_out = open('../../analysis_new/fastq/src/analyse_fastqdata.js', 'w')
    default_info = 'const data =' + json.dumps(fastq_out)
    print(default_info, file=info_out)
    info_out.close()