#!/usr/bin/env python
# coding: utf-8
import os
from DataStructure import SeqRepeatInfo


# In[1]:


# dataset
mouseDataset = 'mm39.fa'
danioRerioDataset = 'DanioRerio.fasta'
humanDataset = 'hg38.fa'
human_chr1 = 'chr1.fa'

# match pattern
humanPattern = "^chr(?:\d*|[A-Z]*)$"
noPattern = "*"
# shared var
currDataset = human_chr1

datasetPath = os.path.join(os.getcwd()) + '/../dataset/' + currDataset
matchPattern = noPattern
fragmentN = 5
cutter = 'GATC'
cutterLen = len(cutter)
commonCount = 5


# %%
