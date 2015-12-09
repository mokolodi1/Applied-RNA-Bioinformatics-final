import numpy as np
from scipy.stats import ttest_ind
from scipy.special import stdtr

# omg python pls
from collections import namedtuple # for JS-like objects

# NOTE: aligner is either 0 or 1
ElementTuple = namedtuple("ElementTuple", "aligner value sample")

interestingGenes = []
with open("data/interesting_genes.txt") as f:
    for line in f:
        # remove trailing "\n"
        interestingGenes.append(line[:-1])

print(interestingGenes)

# geneMapping should be a dictionary from ENSEMBL IDs to HUGO gene names
geneMapping = {}
headerLine = True
with open("data/ensembl_mapping.tsv") as f:
    for line in f:
        line = line[:-1] # remove trailing "\n"
        split = line.split("\t")
        if headerLine:
            headerLine = False
        else:
            geneLabel = split[0]

            # only create mapping if it's in in interestingGenes
            if geneLabel in interestingGenes:
                geneMapping[split[3]] = geneLabel

print("done creating ENSEMBL to gene name mapping")

sampleArray = False
with open("data/star.tsv") as starf, open("data/tophat.tsv") as tophatf:
    for star, tophat in zip(starf, tophatf):
        if not sampleArray:
            assert star == tophat
            # remove trailing "\n", remove first cell ("Gene ID")
            sampleArray = star[:-1].split("\t")[1:]
        else:
            star = star.strip().split("\t")
            tophat = tophat.strip().split("\t")
            assert star[0] == tophat[0]

            # "ENSG00000000003.10" ==> "ENSG00000000003"
            geneLabel = star[0].split(".")[0]
            if geneLabel not in geneMapping:
                continue

            # map from ENSEMBL to gene name
            # "ENSG00000000003" ==> "TSPAN6"
            geneLabel = geneMapping[geneLabel]

            # convert strings to values, remove first column (gene name)
            star = [float(x) for x in star[1:]]
            tophat = [float(x) for x in tophat[1:]]

            # create array of all values (combining star, tophat)
            chartValues = []
            for i in range(len(star)):
                chartValues.append(ElementTuple(0, star[i], sampleArray[i]))
                chartValues.append(ElementTuple(1, tophat[i], sampleArray[i]))

            # sort based on value
            chartValues = sorted(chartValues, key=lambda element: element.value)

            # create a dictionary with sample as key
            samplesDict = {}
            for index, element in enumerate(chartValues):
                if element.sample not in samplesDict:
                    samplesDict[element.sample] = [0, 0]
                # map from value to rank
                samplesDict[element.sample][element.aligner] = index

            totalRankDifference = 0
            for sample in sampleArray:
                dictEntry = samplesDict[sample]
                totalRankDifference += dictEntry[1] - dictEntry[0]

            print(geneLabel, totalRankDifference)


            # # use scipy.stats.ttest_ind.
            # t, p = ttest_ind(star, tophat, equal_var=False)
            # print(p)
