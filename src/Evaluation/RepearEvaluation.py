#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import math
import numpy as np 
import pandas as pd 
from collections import namedtuple


# In[162]:


def getFamilyAccIds():
    dfamHitDataset = pd.read_csv('./Evaluation/hg38_dfam_nrph_chr1.hits', sep="\t", header=0)
    # family accession   
    familyAccIds = list(set(dfamHitDataset['family_acc']))
    familyAccIds = [ re.search("DF\d*", i).group() for i in familyAccIds]
    familyAccIds = sorted(familyAccIds)
    return familyAccIds


# In[165]:


# Fix first row id and length problem
# Fix too long seq Problem

def getDfamSearchData():
    # init DfamRefSeqInfo value
    dfId = ''
    seqLineCount = 1
    currLineCount = 0
    seqLength = 0
    seq = ''
    seqFlag = 0
    DfamRefDataset = open('./Evaluation/Dfam.embl.txt', 'r')
    lines = DfamRefDataset.readlines()
    DfamRefSeqList = []
    for eachLine in lines:
        if re.search("^AC", eachLine):
            dfId = re.search("DF\d*", eachLine).group()
        elif re.search("^SQ", eachLine):
            seqFlag = 1
            currLineCount = 0
            seqLineCount = math.ceil(seqLength/60)
            seqLength = int(re.search("\d+", eachLine).group())
        elif seqFlag == 1 and (currLineCount <= seqLineCount):
            matchList = re.findall("[a-zA-Z]+", eachLine)
            segmentSeq = ''.join(matchList)
            seq += segmentSeq
            currLineCount += 1
            if currLineCount >= seqLineCount:
                DfamRefSeqList.append(DfamRefSeqInfo(id = dfId, length = seqLength, seq = seq))
                seq = ''
                seqFlag = 0
    return DfamRefSeqList

