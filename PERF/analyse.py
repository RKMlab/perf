#! /usr/bin/env python

from __future__ import print_function, division
import sys
import os
import json
from collections import Counter, defaultdict
from Bio import SeqIO
from datetime import datetime
import gzip


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

def analyse(args):
    seq_file = args.input
    repeatsOutFile = args.output.name
    current_dir = os.path.dirname(__file__)
    html_report = os.path.splitext(repeatsOutFile)[0] + '.html'
    print("Generating HTML report. This may take a while..")
    inf = float('inf')
    defaultInfo = {}
    defaultInfo['info'] = {}
    defaultInfo['info']['name'] = seq_file
    defaultInfo['info']['genomeSize'] = 0
    defaultInfo['info']['numSeq'] = 0
    defaultInfo['info']['seqInfo'] = []
    totalSeq = 0
    totalBases = 0
    basesCounter = Counter()
    seqSizes = {}
    if seq_file.endswith('gz'):
        fastaFile = gzip.open(seq_file, 'rt')
    else:
        fastaFile = open(seq_file, 'r')
    for record in SeqIO.parse(fastaFile, 'fasta'):
        totalSeq += 1
        seq = str(record.seq).upper()
        totalBases += len(seq)
        basesCounter.update(seq)
        seqSizes[record.id] = len(seq)
    try:
        GC = (float(basesCounter["G"] + basesCounter["C"])/(totalBases-basesCounter["N"]))*100
    except KeyError:
        GC = (float(basesCounter["G"] + basesCounter["C"])/totalBases)*100
    defaultInfo['info']['genomeSize'] = totalBases
    defaultInfo['info']['GC'] = round(GC, 2)
    defaultInfo['info']['numSeq'] = totalSeq
    totalRepBases = 0
    totalRepFreq = 0
    repFreqByClass = []
    repBasesByClass = []
    chrFreq = {}
    chrBases = {}
    plotData = {'replen': {}, 'repunit': {}}
    plotInfo = {'len': {}, 'unit': {}}
    longestLengths = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
    mostUnits = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
    minLength = inf
    minUnits = inf
    starttime = datetime.now()
    with open(repeatsOutFile, 'r') as repFile:
        for line in repFile:
            line = line.strip()
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
            totalRepBases += repLength
            totalRepFreq += 1

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

            if repClass not in plotData['replen']:
                plotData['replen'][repClass] = {}
                plotData['replen'][repClass][repLength] = 1
                plotData['repunit'][repClass] = {}
                plotData['repunit'][repClass][repUnit] = 1
            elif repClass in plotData['replen']:
                if repLength not in plotData['replen'][repClass]:
                    plotData['replen'][repClass][repLength] = 1
                elif repLength in plotData['replen'][repClass]:
                    plotData['replen'][repClass][repLength] += 1
                if repUnit not in plotData['repunit'][repClass]:
                    plotData['repunit'][repClass][repUnit] = 1
                elif repUnit in plotData['repunit'][repClass]:
                    plotData['repunit'][repClass][repUnit] += 1

    for rep in plotData['replen']:
        freqs = list(plotData['replen'][rep].values())
        repFreqByClass.append({ 'name': rep, 'value': sum(freqs) })
        repBasesByClass.append({ 'name': rep, 'value': sum(list(map(lambda x: x*(minLength + freqs.index(x)), freqs))) })
        lenfreqs = []
        unitfreqs = []
        lengths = sorted(list(plotData['replen'][rep].keys()))
        units = sorted(list(plotData['repunit'][rep].keys()))
        for i in range(minLength, max(lengths) + 1):
            try:
                lenfreqs.append(plotData['replen'][rep][i])
            except KeyError:
                lenfreqs.append(0)
        for i in range(minUnits, max(units) + 1):
            try:
                unitfreqs.append(plotData['repunit'][rep][i])
            except KeyError:
                unitfreqs.append(0)
        plotInfo['len'][rep] = lenfreqs
        plotInfo['unit'][rep] = unitfreqs

    defaultInfo['info']['plotInfo'] = plotInfo
    defaultInfo['info']['numRepClass'] = len(repFreqByClass)
    defaultInfo['info']['totalRepBases'] = totalRepBases
    defaultInfo['info']['totalRepFreq'] = totalRepFreq
    defaultInfo['info']['repFreqByClass'] = repFreqByClass
    defaultInfo['info']['repBasesByClass'] = repBasesByClass
    defaultInfo['info']['percentGenomeCovered'] = str(round((totalRepBases/totalBases)*100, 2)) + "%"
    defaultInfo['info']['repDensityByFreq'] = round((totalRepFreq/totalBases)*1000000, 2)
    defaultInfo['info']['repDensityByBases'] = round((totalRepBases/totalBases)*1000000, 2)
    defaultInfo['info']['minLength'] = minLength
    defaultInfo['info']['minUnits'] = minUnits
    defaultInfo['info']['longestRepeats'] = []
    defaultInfo['info']['mostRepeatUnits'] = []
    for a in longestLengths:
        testDict = {'seq': a[0], 'start': a[1], 'end': a[2], 'repClass': a[3], 'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]}
        defaultInfo['info']['longestRepeats'].append(testDict)
    for a in mostUnits:
        testDict = {'seq': a[0], 'start': a[1], 'end': a[2], 'repClass': a[3], 'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]}
        defaultInfo['info']['mostRepeatUnits'].append(testDict)
    seqSizes = Counter(seqSizes)
    for a in seqSizes.most_common(100):
        seqInfo = {}
        seqInfo['name'] = a[0]
        seqInfo['size'] = a[1]
        try:
            seqInfo['freq'] = chrFreq[a[0]]
            seqInfo['bp'] = chrBases[a[0]]
            seqInfo['freqDens'] = str(round((chrFreq[a[0]]/int(a[1]))*1000000, 2))
            seqInfo['bpDens'] = str(round((chrBases[a[0]]/int(a[1]))*1000000, 2))
        except KeyError:
            seqInfo['freq'] = 0
            seqInfo['bp'] = 0
            seqInfo['freqDens'] = "0.0"
            seqInfo['bpDens'] = "0.0"
        defaultInfo['info']['seqInfo'].append(seqInfo)
    defaultInfo = 'const data =' + json.dumps(defaultInfo)
    writetoHTML(html_report, defaultInfo)
