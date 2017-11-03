#!usr/bin/python
from __future__ import print_function, division
from operator import itemgetter
import argparse
from tqdm import tqdm
import os
import gzip
from utils import rawcharCount

"""

    CAUTION: Works currently for only sorted bed files.

    Preferential order of assigning annotation:
        Promoter >> Overlapping >> Intergenic

    Defaults:
        > Promoter distance is 1kb upstream and downstream of TSS.
        > Gene id considered is "gene".
"""

def selectAnnotation(List):
    if 'Exon' in List:
        return 'Exon'
    elif 'Intron' in List:
        return 'Intron'
    elif 'Genic' in List:
        return 'Genic'
    elif 'Intergenic' in List:
        return 'Intergenic'


def promoter(check):
    if check == 1:
        return 'Promoter'
    else:
        return 'Non-Promoter'


# Need to be updated for better parsing of the attributes
def processAttrs(attribute):
    attrObj = {}
    attributes = attribute.split(";")
    for a in attributes:
        attr = a.split("=")
        attrName = attr[0]
        attrObj[attrName] = attr[1]
    return attrObj


def processGFF(GFF, geneId):
    #Order of columns seqname, source, feature, start, end, score, strand, frame, attribute
    geneObj = {}
    subGeneObj = {}
    if (GFF.endswith('gz')):
        gff = gzip.open(GFF)
    else:
        gff = open(GFF)
    for line in gff:
        line = line.decode('utf-8')
        line = line.strip()
        if line.startswith('#'):
            pass
        else:
            fields = line.split('\t')
            seqname = fields[0]
            source = fields[1]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            frame = fields[7]
            attribute = fields[8]
            attrObj = processAttrs(attribute)
            if feature == "gene":
                try:
                    geneObj[seqname].append([attrObj[geneId], start, end, strand])
                except KeyError:
                    geneObj[seqname] = [[attrObj[geneId], start, end, strand]]
            elif feature == 'exon':
                try:
                    subGeneObj[attrObj[geneId]][feature].append([start, end, strand])
                except KeyError:
                    try:
                        subGeneObj[attrObj[geneId]][feature] = [[start, end, strand]]
                    except KeyError:
                        subGeneObj[attrObj[geneId]] = {feature: [[start, end, strand]]}
    for i in geneObj:
        geneObj[i] = sorted(geneObj[i], key=itemgetter(1))
    for a in subGeneObj:
        for b in subGeneObj[a]:
            subGeneObj[a][b] = sorted(subGeneObj[a][b], key=itemgetter(0))

    return {'gene': geneObj, 'subgene': subGeneObj}


def annotate(args):
    rep_file = args.output.name
    gff_file = args.gff_input
    output_file = open(os.path.splitext(rep_file)[0] + '_annotation.tsv', 'w')
    geneId = args.gene_attribute

    promUp = args.up_promoter
    promDown = args.down_promoter

    gffObject = processGFF(gff_file, geneId)
    geneObj = gffObject['gene']
    subGeneObj = gffObject['subgene']

    print('', end='\n')
    # Counting the number of lines in bed -------------------------------------
    num_records = rawcharCount(rep_file, '\n')
    with open(rep_file) as bed:
        prevSeqName = "Initialise" # Initialise for checking the prevSeqName
        minStartIndex = 0
        for line in tqdm(bed, total=num_records):
            # Object for the output entries to be appended --------------------
            Annotations = {'Genic': [], 'Exon': [], 'Intron': []}
            line = line.strip()
            fields = line.split('\t')
            seqname = fields[0]
            """
                If the seqname is not same the previous seq name the check
                starts from the first gene on the sequence.
            """
            if seqname != prevSeqName:
                minStartIndex = 0
            prevSeqName = seqname
            S1 = int(fields[1])
            E1 = int(fields[2])
            leastDist = float('inf')
            breakCheck = 0
            promoterCheck = 0
            try:
                for i, a in enumerate(geneObj[seqname][minStartIndex:]):
                    annotation = ''
                    geneName = a[0]
                    try:
                        subgeneElements = subGeneObj[geneName]
                    except KeyError:
                        subgeneElements = {}
                    S2 = a[1]
                    E2 = a[2]
                    Ori = a[3]
                    # Transcription Start site
                    TSS = S2
                    if Ori == '-':
                        TSS = E2
                    """
                        Calculating leastStart distance -
                          > Storing the index of the gene which has least distance,
                            greater than 1000 , upstream from the start of the repeat.
                    """
                    if i == 0:
                        leastStart = S1 - E2
                        minIndex = i
                    else:
                        if S1 - E2 > 1000:
                            if leastStart > (S1 - E2):
                                leastStart = S1 - E2
                                minIndex = i
                    if breakCheck == 1:
                        break
                    if (S2 - E1 > 1000):
                        breakCheck = 1
                    # Checking if region comes in promoter --------------------
                    # For positive strand orientation -------------------------
                    if Ori == '-' and (TSS-promDown <= S1 <= TSS+promUp or TSS-promDown <= E1 <= TSS+promUp):
                        promoterCheck = 1
                    elif Ori == '+' and (TSS-promUp <= S1 <= TSS+promDown or TSS-promUp <= E1 <= TSS+promDown):
                        promoterCheck = 1
                    # If no Promoter found ------------------------------------
                    # Checking if it overlaps ---------------------------------
                    if (E2 - S1 >=0 and S2 - S1 <=0) or (E2 - E1 >=0 and S2 - E1 <= 0):
                        annotation = 'Genic'
                        # Removes the Intergenic entries ------------------
                        Annotations['Intergenic'] = []
                        TSS = S2
                        diffSS = S2 - S1
                        diffES = E2 - S1
                        diffSE = S2 - E1
                        diffEE = E2 - E1
                        if abs(diffSS) < abs(diffSE):
                            TSSdist = diffSS
                        else:
                            TSSdist = diffSE
                        distance = TSSdist
                        # Checking overlap with subgene
                        if 'exon' in subgeneElements:
                            for site in subgeneElements['exon']:
                                S3 = site[0]
                                E3 = site[1]
                                if (E3 - S1 >=0 and S3 - S1 <=0) or (E3 - E1 >=0 and S3 - E1 <= 0):
                                    annotation = "Exon"
                                    break
                                else:
                                    annotation = "Intron"
                    elif len(Annotations['Exon']) == 0 and len(Annotations['Intron']) == 0 and len(Annotations['Genic']) == 0:
                        TSS = S2
                        diffSS = S2 - S1
                        diffES = E2 - S1
                        diffSE = S2 - E1
                        diffEE = E2 - E1
                        if abs(diffSS) < abs(diffSE):
                            TSSdist = diffSS
                        else:
                            TSSdist = diffSE
                        minDistance = min([abs(diffSS), abs(diffEE), abs(diffSE), abs(diffES)])
                        annotation = 'Intergenic'
                        distance = TSSdist
                        if minDistance < leastDist:
                            leastDist = minDistance
                            Annotations[annotation] = [line + '\t' + '\t'.join(str(b) for b in a) + '\t' + annotation + '\t' + promoter(promoterCheck) + '\t' + str(distance)]

                    if (annotation == "Exon" or annotation == "Intron" or annotation == "Genic"):
                        Annotations[annotation].append(line + '\t' + '\t'.join(str(b) for b in a) + '\t' + annotation + '\t' + promoter(promoterCheck) + '\t' + str(distance))

                minStartIndex += minIndex
                if minStartIndex > 0:
                    minStartIndex = minStartIndex - 1
            #If sequence is not found, reports as annotation not available
            except KeyError:
                Annotations = {'Genic': [], 'Exon': [], 'Intron': []}
                print(line + '\t' + '\t'.join(['-']*7), file = output_file)
            exclude = []
            for anno in list(Annotations.keys()):
                if len(Annotations[anno]) == 0:
                    del Annotations[anno]
            for anno in Annotations:
                leastEntryDist = float('inf')
                leastEntry = ""
                for entry in Annotations[anno]:
                    EntryDist = int(entry.split('\t')[-1])
                    if EntryDist < leastEntryDist:
                        leastEntryDist = EntryDist
                        leastEntry = entry
                if leastEntry != "":
                    Annotations[anno] = leastEntry
            if len(Annotations) > 1:
                annoSelected = selectAnnotation(list(Annotations.keys()))
                print(Annotations[annoSelected], file = output_file)
            else:
                for anno in Annotations:
                    print(Annotations[anno], file = output_file)
    output_file.close()