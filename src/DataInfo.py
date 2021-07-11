#!/usr/bin/env python
# coding: utf-8
import re
import os
from DataStructure import RepeatFragNInfo


# In[1]:


# dataset
mouseDataset = "mm39.fa"
danioRerioDataset = "DanioRerio.fasta"
humanDataset = "hg38.fa"
human_chr1 = "chr1.fa"
human_chrY = "chrY.fa"
human_chr17 = "chr17.fa"
human_chr19 = "chr19.fa"

# match pattern
chrPattern = "^chr(?:\d*|[A-Z]*)$"
noPattern = "*"
# shared var
currDataset = human_chr19
currDatasetName = re.search(r"([^.]+)", currDataset).group()

datasetPath = os.path.join(os.getcwd()) + "/../dataset/" + currDataset
matchPattern = noPattern
fragmentN = 1
cutter = "GATC"
cutterLen = len(cutter)
commonCount = 5


# %%
