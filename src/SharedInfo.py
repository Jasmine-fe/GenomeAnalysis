#!/usr/bin/env python
# coding: utf-8
import re
import os
from DataStructure import RepeatFragNInfo


# Dataset
humanDataset = "hg38.fa"
human_chr1 = "chr1.fa"
human_chr17 = "chr17.fa"
humanChrKey = [i for i in range(1, 23)] + ["X", "Y"]
dmChrX = "dm6/chrX_sequence.fasta"
currDataset = dmChrX
currDatasetName = re.search(r"([^.]+)", currDataset).group()
datasetPath = os.path.join(os.getcwd()) + "/../dataset/" + currDataset

# match pattern
chrPattern = "^chr(?:\d*|[A-Z]*)$"
noPattern = "*"
matchPattern = noPattern

# shared var
fragmentN = 1

cutterA = "GATC"
cutterB = "AAGCTT"

# Plot
colorA = "#DF6F31"
colorB = "#658EA9"
