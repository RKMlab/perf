#! /usr/bin/env python

from __future__ import print_function, division
import sys, os, json
from collections import Counter, defaultdict
import numpy as np
from pprint import pprint

if sys.version_info.major == 2:
    from utils import rev_comp, kmers, get_cycles, build_cycVariations
elif sys.version_info.major == 3:
    from .utils import rev_comp, kmers, get_cycles, build_cycVariations


def writetoHTML(html_file, defaultInfo, repeat_options, input_format):
    html_handle = open(html_file, 'w')
    current_dir = os.path.dirname(__file__)

    template = open('%s/lib/template_%s.html' %(current_dir, input_format), 'r').read()

    fontawesome_js = open('%s/lib/src/all.js' %(current_dir), 'r').read()
    semantic_css = open('%s/lib/styles/semantic.min.css' %(current_dir), 'r').read()
    multiselect_css = open('%s/lib/styles/multi-select.min.css' %(current_dir), 'r').read()
    apexcharts_css = open('%s/lib/styles/apexcharts.min.css' %(current_dir), 'r').read()
    main_css = open('%s/lib/styles/main.css' %(current_dir), 'r').read()

    jquery_js = open("%s/lib/src/jquery-3.5.0.min.js" %(current_dir), "r").read()
    semantic_js = open("%s/lib/src/semantic.min.js" %(current_dir), "r").read()
    multiselect_js = open('%s/lib/src/jquery.multi-select.min.js' %(current_dir), 'r').read()
    apexcharts_js = open('%s/lib/src/apexcharts.min.js' %(current_dir), 'r').read()
    lodash_js = open('%s/lib/src/lodash.min.js' %(current_dir), 'r').read()
    main_js = open('%s/lib/src/main_%s.js' %(current_dir, input_format), 'r').read()
    tables_js = open('%s/lib/src/tables_%s.js' %(current_dir, input_format), 'r').read()
    annocharts_js = ''
    if input_format == 'fasta':
        annocharts_js = open('%s/lib/src/anno_charts.js' %(current_dir), 'r').read()

    template = template.format(
        fontawesome_js = fontawesome_js, 
        semantic_css = semantic_css, 
        multiselect_css = multiselect_css, 
        apexcharts_css = apexcharts_css, 
        main_css = main_css, 
        jquery_js = jquery_js, 
        semantic_js = semantic_js, 
        multiselect_js = multiselect_js, 
        apexcharts_js = apexcharts_js, 
        lodash_js = lodash_js, 
        analyse_data_js = defaultInfo, 
        main_js = main_js, 
        tables_js = tables_js, 
        annocharts_js = annocharts_js,
        repeat_options = repeat_options,
    )

    print(template, file=html_handle)
    html_handle.close()
    print("HTML report successfully saved to " + html_file)


def get_parameters(args):
    runCommand = 'PERF' + ' '.join(sys.argv)


def analyse_fasta(args):
    repeatsOutFile = args.output.name
    html_report = os.path.splitext(repeatsOutFile)[0] + '.html'
    print("\nGenerating HTML report. This may take a while..", end="\n\n")
    

    all_repeat_classes = []
    kmer_classes = defaultdict(list)
    cyclical_variations = dict()
    for r in args.repeats:
        r = r.split('\t')[1]
        if r not in all_repeat_classes:
            all_repeat_classes.append(r)
            cyclical_variations[r] = build_cycVariations(r)

    inf = float('inf')
    defaultInfo = {}
    defaultInfo['info'] = { 'seqInfo': {}, 'repInfo': {} }
    
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

    plot_data = dict()
    with open(repeatsOutFile, 'r') as repFile:
        for line in repFile:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                defaultInfo['info']['seqInfo'][fields[0]] = fields[1]
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

                if repClass not in plot_data:
                    plot_data[repClass] = dict()
                    plot_data[repClass][repLength] = [0]*len(cyclical_variations[repClass])
                if repLength not in plot_data[repClass]: 
                    plot_data[repClass][repLength] = [0]*len(cyclical_variations[repClass])
                plot_data[repClass][repLength][cyclical_variations[repClass].index(actualRepeat)] += 1

                if minUnits > repUnit: minUnits = repUnit
                if minLength > repLength: minLength = repLength

                if (longestLengths[-1][4] < repLength) or (longestLengths[-1][4] == repLength and repClass < longestLengths[-1][3]):
                    longestLengths[-1] = fields
                    longestLengths.sort(key=lambda x: x[4])
                    longestLengths.reverse()
                if (mostUnits[-1][6] < repUnit) or (mostUnits[-1][6] == repUnit and repClass < longestLengths[-1][3]):
                    mostUnits[-1] = fields
                    mostUnits.sort(key=lambda x: x[6])
                    mostUnits.reverse()
    for r in all_repeat_classes:
        kmer_classes[kmers[len(r)]].append(r)
        if r not in plot_data:
            plot_data[r] = 0

    repeat_options = ""
    for kmer in kmer_classes:
        repeat_options += '<optgroup label="%s">' %(kmer)
        for r in kmer_classes[kmer]:
            repeat_options += '<option value="%s">%s</option>' %(r, r)
        repeat_options += '</optgroup>'
    
    totalBases = int(defaultInfo['info']['seqInfo']['Total_bases'])
    defaultInfo['info']['repInfo']['lenFrequency'] = plot_data
    defaultInfo['info']['repInfo']['numRepClasses'] = str(len(plot_data.keys())) + '/' + str(len(all_repeat_classes))
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
    writetoHTML(html_report, defaultInfo, repeat_options, 'fasta')

def analyse_fastq(args, fastq_out):

    """Generates HTML report for fastq files."""
    html_report = os.path.splitext(args.output.name)[0] + '.html'
    print("\nGenerating HTML report. This may take a while..", end="\n\n")
    
    fastq_out['info']['seqInfo']['File_name'] = args.input.split('/')[-1]
    n = fastq_out['info']['seqInfo']['Total_reads']
    b = fastq_out['info']['seqInfo']['Total_bases']
    total_repeats = fastq_out['info']['repInfo']['totalRepFreq']
    reads_with_repeats = fastq_out['info']['repInfo']['totalRepReads']
    total_repeat_bases = fastq_out['info']['repInfo']['totalRepBases']
    all_repeat_classes = list(map(lambda x: x.split('\t')[1], args.repeats))
    temp = []
    for a in all_repeat_classes:
        if a not in temp:
            temp.append(a)
    all_repeat_classes = temp
    del temp

    kmer_classes = defaultdict(list)
    for r in all_repeat_classes:
        kmer_classes[kmers[len(r)]].append(r)
    repeat_options = ""
    for kmer in kmer_classes:
        repeat_options += '<optgroup label="%s">' %(kmer)
        for r in kmer_classes[kmer]:
            repeat_options += '<option value="%s">%s</option>' %(r, r)
        repeat_options += '</optgroup>'

    fastq_out['info']['repInfo']['numRepClasses'] = str(fastq_out['info']['repInfo']['numRepClasses']) + '/' + str(len(all_repeat_classes))
    fastq_out['info']['repInfo']['allRepClasses'] = all_repeat_classes
    fastq_out['info']['repInfo']['totalRepFreqNorm'] = round((total_repeats/n)*1000000, 2)
    fastq_out['info']['repInfo']['totalRepReadsNorm'] = str(round((reads_with_repeats/n)*100, 2)) + '%'
    fastq_out['info']['repInfo']['percentRepBases'] = str(round((total_repeat_bases/b)*100, 2)) + '%'

    rep_fastq_info = fastq_out['info']['repInfo']['repClassInfo']
    for rep in sorted(rep_fastq_info, key= lambda k: (len(k), k)):
        fastq_out['info']['repInfo']['repClassInfo'][rep]['reads_norm'] = round((rep_fastq_info[rep]['reads']/n)*1000000, 3)
        fastq_out['info']['repInfo']['repClassInfo'][rep]['instances_norm'] = round((rep_fastq_info[rep]['instances']/n)*1000000, 3)
        fastq_out['info']['repInfo']['repClassInfo'][rep]['bases_norm'] = round((rep_fastq_info[rep]['bases']/b)*1000000, 3)
    
    defaultInfo = 'const data =' + json.dumps(fastq_out)
    writetoHTML(html_report, defaultInfo, repeat_options, 'fastq')
