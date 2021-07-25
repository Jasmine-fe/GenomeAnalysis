#!/usr/bin/env python
# coding: utf-8
import re
import os
from DataStructure import RepeatFragNInfo


# In[1]:


# dataset
humanDataset = "hg38.fa"
human_chr1 = "chr1.fa"
human_chr17 = "chr17.fa"
humanChrKey = [i for i in range(1, 23)] + ["X", "Y"]
dmChrX = "dm6/chrX_sequence.fasta"
# match pattern
chrPattern = "^chr(?:\d*|[A-Z]*)$"
noPattern = "*"

# shared var
currDataset = dmChrX
currDatasetName = re.search(r"([^.]+)", currDataset).group()

datasetPath = os.path.join(os.getcwd()) + "/../dataset/" + currDataset
matchPattern = noPattern
fragmentN = 1

cutterA = "GATC"
cutterB = "AAGCTT"
cutter = cutterB  # current

commonCount = 5


# %%
