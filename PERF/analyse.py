#! /usr/bin/env python

from __future__ import print_function, division
import sys
import os
import json
from collections import OrderedDict as OD
from collections import Counter
from Bio import SeqIO


def writetoHTML(html_file):
    html_handle = open(html_file, 'w')
    current_dir = os.path.dirname(__file__)
    with open(current_dir + '/lib/template.html') as report:
        for line in report:
            line = line.strip()
            print(line, file=html_handle)
            try:
                start_index = line.index("^^")
                stop_index = line.index("$$")
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
    analyseDataOUT = open(current_dir + '/lib/src/data.js', 'w')
    html_report = os.path.splitext(repeatsOutFile)[0] + '.html'
    print("Generating HTML report. This may take a while..")

    defaultInfo = OD()
    defaultInfo['info'] = OD()
    defaultInfo['info']['name'] = seq_file
    defaultInfo['info']['genomeSize'] = 0
    defaultInfo['info']['numSeq'] = 0
    defaultInfo['info']['seqInfo'] = []
    totalSeq = 0
    totalBases = 0
    basesCounter = Counter()
    seqSizes = {}
    with open(seq_file, "rt") as fastaFile:
        for record in SeqIO.parse(fastaFile, 'fasta'):
            totalSeq += 1
            # print("Processing %s" % (record.id), file=sys.stderr)
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
    repFreqByClass = {}
    repBasesByClass = {}
    chrFreq = {}
    chrBases = {}
    plotData = {}
    plotInfo = {}
    longestLengths = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
    mostUnits = [['seq', 'start', 'stop', 'repClass', 0, '+', 0, 'actualrep']]*100
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

            if repClass not in plotData:
                plotData[repClass] = {}
                plotData[repClass][repLength] = 1
            elif repClass in plotData:
                if repLength not in plotData[repClass]:
                    plotData[repClass][repLength] = 1
                elif repLength in plotData[repClass]:
                    plotData[repClass][repLength] += 1

            totalRepBases += repLength
            totalRepFreq += 1
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
            try:
                repFreqByClass[repClass] += 1
                repBasesByClass[repClass] += repLength
            except KeyError:
                repFreqByClass[repClass] = 1
                repBasesByClass[repClass] = repLength
            try:
                chrFreq[seq] += 1
                chrBases[seq] += repLength
            except KeyError:
                chrFreq[seq] = 1
                chrBases[seq] = repLength

    temp = []
    for a in repFreqByClass:
        Obj = {'name': a, 'value': repFreqByClass[a]}
        temp.append(Obj)
    repFreqByClass = temp
    temp = []
    for a in repBasesByClass:
        Obj = {'name': a, 'value': repBasesByClass[a]}
        temp.append(Obj)
    repBasesByClass = temp
    temp = []
    for o in plotData:
        freqs = []
        lengths = sorted(list(plotData[o].keys()))
        for i in range(min(lengths), max(lengths) + 1):
            try:
                freqs.append(plotData[o][i])
            except KeyError:
                freqs.append(0)
        plotInfo[o] = freqs

    defaultInfo['info']['plotInfo'] = plotInfo
    defaultInfo['info']['numRepClass'] = len(repFreqByClass)
    defaultInfo['info']['totalRepBases'] = totalRepBases
    defaultInfo['info']['totalRepFreq'] = totalRepFreq
    defaultInfo['info']['repFreqByClass'] = repFreqByClass
    defaultInfo['info']['repBasesByClass'] = repBasesByClass
    defaultInfo['info']['percentGenomeCovered'] = str(round((totalRepBases/totalBases)*100, 2)) + "%"
    defaultInfo['info']['repDensityByFreq'] = round((totalRepFreq/totalBases)*1000000, 2)
    defaultInfo['info']['repDensityByBases'] = round((totalRepBases/totalBases)*1000000, 2)
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
    print('const data =', json.dumps(defaultInfo), file=analyseDataOUT)
    analyseDataOUT.close()
    writetoHTML(html_report)
