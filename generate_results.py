from statistics import pstdev
from scipy.stats import norm
from scipy.stats import wilcoxon

from random import randint

# omg python pls
from collections import namedtuple # for JS-like objects

# NOTE: aligner is either 0 or 1
ElementTuple = namedtuple("ElementTuple", "aligner value sample")

interestingGenes = []
with open("data/interesting_genes.txt") as f:
    for line in f:
        # remove trailing "\n"
        interestingGenes.append(line[:-1])

print("interestingGenes length: ", len(interestingGenes))

# geneMapping should be a dictionary from ENSEMBL IDs to HUGO gene names
geneMapping = {}
headerLine = True
multiplyBy = 319;
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

print("geneMapping length: ", len(geneMapping))

print()
print("gene\taverage rank difference\tp-value")

line = 0
genesTested = 0
sampleArray = False
with open("data/star.tsv") as starf, open("data/tophat.tsv") as tophatf:
    for star, tophat in zip(starf, tophatf):
        line += 1
        # if line % 1000 == 0:
        #     print("line ", line)
        if not sampleArray:

            # remove trailing "\n", remove first cell ("Gene ID")
            starArray = star[:-1].split("\t")[1:]
            tophatArray = star[:-1].split("\t")[1:]
            for first, second in zip(starArray, tophatArray):
                if first != second:
                    print(starArray, tophatArray)
            assert starArray == tophatArray

            sampleArray = starArray
        else:
            star = star.strip().split("\t")
            tophat = tophat.strip().split("\t")
            assert star[0] == tophat[0]

            # "ENSG00000000003.10" ==> "ENSG00000000003"
            ensembl = star[0]
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
            # lastValue = False
            samplesDict = {}
            for index, element in enumerate(chartValues):
                # if lastValue:
                #     if element.value == lastValue:
                #         print(element.value, " duplicated")
                # lastValue = element.value

                if element.sample not in samplesDict:
                    samplesDict[element.sample] = [0, 0]
                # map from value to rank
                samplesDict[element.sample][element.aligner] = index

            differences = []
            for sample in sampleArray:
                dictEntry = samplesDict[sample]
                diff = dictEntry[1] - dictEntry[0]
                differences.append(diff)
            totalRankDifference = sum(differences)

            statistic, pvalue = wilcoxon(star, tophat)
            print(ensembl, "\t", geneLabel, "\t", totalRankDifference / len(star), "\t", pvalue * multiplyBy)

            genesTested += 1

assert multiplyBy == genesTested
print("genes tested: ", genesTested)
