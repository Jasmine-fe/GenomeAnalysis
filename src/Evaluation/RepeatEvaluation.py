#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import math
import numpy as np
import pandas as pd
from collections import namedtuple
from DataStructure import DfamConSeqInfo, OutputMatchPositionInfo
from DataInfo import (
    cutterLen,
    fragmentN,
)

# In[162]:


def getConsensusSeqIds():
    dfamHitDataset = pd.read_csv(
        "Evaluation/Source/hg38_dfam_nrph_chr1.hits", sep="\t", header=0
    )

    familyAcIds = list(set(dfamHitDataset["family_acc"]))
    familyAcIds = [re.search("DF\d*", i).group() for i in familyAcIds]
    familyAcIds = sorted(familyAcIds)
    familyAcIdsDf = pd.DataFrame(familyAcIds, columns=["id"])
    return familyAcIdsDf


# In[165]:
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


# In[ ]:
def getOutputMatchPosition(repeatInfoList):
    matchList = []
    for repeatN in repeatInfoList:
        repeatFragNLen = sum(repeatN.fragmentLenList) + (cutterLen * fragmentN)
        for fragposition in repeatN.position:
            matchList.append(
                OutputMatchPositionInfo(
                    fragposition.baseIdx, fragposition.baseIdx + repeatFragNLen
                )
            )
    matchList.sort(key=lambda x: x[0])
    return matchList


# In[ ]:
"""
Params
outputPositionList: [ OutputMatchPositionInfo(startIdx, endIdx), ...]
"""


def outputPositionBucketClassifier(outputPositionList, bucketNum=10):
    outputPositionLookupDic = {}
    eachBucketNum = int(len(outputPositionList) / bucketNum)
    for i in range(bucketNum):
        firstStartIdx, lastStartIdx = (
            outputPositionList[eachBucketNum * i].startIdx,
            outputPositionList[(eachBucketNum * (i + 1) - 1)].startIdx,
        )
        outputPositionLookupDic[i] = (firstStartIdx, lastStartIdx)
    outputPositionLookupDic[bucketNum - 1] = (
        outputPositionLookupDic[bucketNum - 1][0],
        outputPositionList[-1].startIdx,
    )  # handle eachBucketNum remainder
    return outputPositionLookupDic


# In[ ]:
def checkDfamMatch(
    outputPositionLookupDic,
    OutputMatchPositionList,
    DfamFileName="hg38_dfam_nrph_chr1",
    bucketNum=10,
):
    dfamDataset = pd.read_csv(
        f"Evaluation/Source/{DfamFileName}.hits",
        sep="\t",
        header=0,
        usecols=["familyAcc", "familyName", "aliSt", "aliEn", "envSt", "envEn"],
    )
    dfamDatasetMatchList = [False] * len(dfamDataset)
    eachBucketNum = int(len(OutputMatchPositionList) / bucketNum)
    for idx in range(len(dfamDataset)):
        row = dfamDataset.iloc[idx]
        li = sorted([row.aliSt, row.aliEn, row.envSt, row.envEn])
        start = li[0]
        end = li[-1]
        length = end - start
        flag = 0
        for bucketIdx in range(bucketNum):
            if (
                start >= outputPositionLookupDic[bucketIdx][0]
                and start <= outputPositionLookupDic[bucketIdx][1]
            ):
                for outputRow in OutputMatchPositionList[
                    eachBucketNum * bucketIdx : eachBucketNum * (bucketIdx + 1) - 1
                ]:
                    if (
                        (outputRow.startIdx >= start and outputRow.startIdx <= end)
                        or (outputRow.endIdx <= start and outputRow.endIdx >= end)
                        or (
                            (outputRow.startIdx + outputRow.endIdx) / 2 >= start
                            and (outputRow.startIdx + outputRow.endIdx) / 2 <= end
                        )
                    ):
                        flag = 1
                        dfamDatasetMatchList[idx] = True
                        break
                if flag:
                    break
    return dfamDatasetMatchList
