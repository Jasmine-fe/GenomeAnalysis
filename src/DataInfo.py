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
# match pattern
humanPattern = "^chr(?:\d*|[A-Z]*)$"
noPattern = "*"
# shared var
currDataset = human_chr1
currDatasetName = re.search(r"([^.]+)", currDataset).group()

datasetPath = os.path.join(os.getcwd()) + "/../dataset/" + currDataset
matchPattern = noPattern
fragmentN = 5
cutter = "GATC"
cutterLen = len(cutter)
commonCount = 5


# %%
