#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import math
import numpy as np
import pandas as pd
from collections import namedtuple
from DataStructure import DfamConSeqInfo

# In[162]:


def getHitFamilyACIds():
    dfamHitDataset = pd.read_csv(
        "Evaluation/Source/hg38_dfam_nrph_chr1.hits", sep="\t", header=0
    )
    # family accession
    familyAccIds = list(set(dfamHitDataset["family_acc"]))
    familyAccIds = [re.search("DF\d*", i).group() for i in familyAccIds]
    familyAccIds = sorted(familyAccIds)
    return familyAccIds


# In[165]:


# Fix first row id and length problem
# Fix too long seq Problem


def getDfamConSeqData():
    # init DfamConSeqInfo value
    dfId = ""
    seqLineCount = 1
    currLineCount = 0
    seqLength = 0
    seq = ""
    seqFlag = 0

    DfamRefDataset = open("Evaluation/Source/Dfam.embl.txt", "r")
    lines = DfamRefDataset.readlines()
    DfamRefSeqList = []
    for eachLine in lines:
        if re.search("^AC", eachLine):
            dfId = re.search("DF\d*", eachLine).group()
        elif re.search("^SQ", eachLine):
            seqFlag = 1
            currLineCount = 0
            seqLineCount = math.ceil(seqLength / 60)
            seqLength = int(re.search("\d+", eachLine).group())
        elif seqFlag == 1 and (currLineCount <= seqLineCount):
            matchList = re.findall("[a-zA-Z]+", eachLine)
            segmentSeq = "".join(matchList)
            seq += segmentSeq
            currLineCount += 1
            if currLineCount >= seqLineCount:
                DfamRefSeqList.append(
                    DfamConSeqInfo(id=dfId, length=seqLength, seq=seq)
                )
                seq = ""
                seqFlag = 0
    return DfamRefSeqList


# For export ConsensusSeqs file

# consensusSeqInfo = getDfamConSeqData()
# consensusSeqList = [Seq(i.seq) for i in consensusSeqInfo]
# idList = [i.id for i in consensusSeqInfo]
# conSeqfragLenList, conSeqfragSeqList = parseSeqByCutter(consensusSeqList)
# filteredLenList = filter(lambda x: len(x) >= 3, conSeqfragSeqList)

# filteredIntegratedConSeqInfo = []
# for id, lengthList, seqList in zip(idList, conSeqfragLenList, conSeqfragSeqList):
#     if len(lengthList) >= 3:
#         dataDic = {"id": id, "lengthList": lengthList[1:-1], "seqList": seqList[1:-1]}
#         filteredIntegratedConSeqInfo.append(dataDic)

# filteredIntegratedConSeqDf = pd.DataFrame(
#     filteredIntegratedConSeqInfo, columns=["id", "lengthList", "seqList"]
# )

# filteredIntegratedConSeqDf.to_csv("./Evaluation/Source/ConsensusSeqs", index=False)
