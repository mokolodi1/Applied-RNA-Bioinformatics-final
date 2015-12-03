import numpy as np
from scipy.stats import ttest_ind
from scipy.special import stdtr

# geneMapping should be a dictionary from ENSEMBL IDs to HUGO gene names
geneMapping = {};
headerLine = True;
with open("ensembl_mapping.tsv") as f:
    for line in f:
        line = line[:-1] # remove trailing "\n"
        split = line.split("\t")
        if headerLine:
            headerLine = False
        else:
            geneMapping[split[3]] = split[0]

        if split[3] == "ENSG00000000003.10":
            print("asdfasdf")

headerLine = True;
with open("star_topleft.tsv") as starf, open("tophat_topleft.tsv") as tophatf:
    for star, tophat in zip(starf, tophatf):
        if headerLine:
            headerLine = False
        else:
            star = star.strip().split("\t")
            tophat = tophat.strip().split("\t")
            assert star[0] == tophat[0]

            geneLabel = geneMapping[star[0]];
            star = [float(x) for x in star[1:]]
            tophat = [float(x) for x in tophat[1:]]

            print(star)
            print(tophat)

            # use scipy.stats.ttest_ind.
            t, p = ttest_ind(star, tophat, equal_var=False)
            print(p)
