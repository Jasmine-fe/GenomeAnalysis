#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import os
import math
import numpy as np
import pandas as pd
from collections import namedtuple
from DataStructure import DfamConSeqInfo, PositionInfo
from DataInfo import (
    cutterLen,
    fragmentN,
)

"""
Params
positionList: [(startIdx, endIdx), ...]
"""


def positionBucketClassifier(positionList, bucketNum=10):
    positionLookupDic = {}
    eachBucketNum = int(len(positionList) / bucketNum)
    for i in range(bucketNum):
        firstStartIdx, lastStartIdx = (
            positionList[eachBucketNum * i].startIdx,
            positionList[(eachBucketNum * (i + 1) - 1)].startIdx,
        )
        positionLookupDic[i] = (firstStartIdx, lastStartIdx)
    # handle eachBucketNum remainder
    positionLookupDic[bucketNum - 1] = (
        positionLookupDic[bucketNum - 1][0],
        positionList[-1].startIdx,
    )
    return positionLookupDic


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
def getRepeatPositionList(repeatInfoList):
    matchList = []
    for repeatN in repeatInfoList:
        repeatFragNLen = sum(repeatN.fragmentLenList) + (cutterLen * fragmentN)
        for fragposition in repeatN.position:
            matchList.append(
                PositionInfo(
                    fragposition.baseIdx, fragposition.baseIdx + repeatFragNLen
                )
            )
    matchList.sort(key=lambda x: x[0])
    return matchList


def generateDfamPositionData(
    readFileName,
    outputFileName,
):
    dfamDataset = pd.read_csv(
        f"Evaluation/Source/{readFileName}.hits",
        sep="\t",
        header=0,
        usecols=["familyAcc", "familyName", "aliSt", "aliEn", "envSt", "envEn"],
    )
    ouputFilePath = f"Evaluation/Source/{outputFileName}.txt"
    outputCols = ["familyAcc", "familyName", "startIdx", "endIdx"]
    with open(ouputFilePath, "w") as outputFile:
        outputFile.write("\t".join(outputCols) + "\n")
        for idx in range(len(dfamDataset)):
            row = dfamDataset.iloc[idx]
            li = sorted([row.aliSt, row.aliEn, row.envSt, row.envEn])
            start, end = li[0], li[-1]
            output = f"{row.familyAcc}\t{row.familyName}\t{start}\t{end}\n"
            outputFile.write(output)


def getDfamPositionList(readFileName):
    dfamDataset = pd.read_csv(
        f"Evaluation/Source/{readFileName}.txt",
        sep="\t",
        header=0,
        usecols=["startIdx", "endIdx"],
    )
    positionList = []
    for idx in range(len(dfamDataset)):
        row = dfamDataset.iloc[idx]
        positionList.append(PositionInfo(row.startIdx, row.endIdx))
    positionList.sort(key=lambda tup: tup[0])
    return positionList


# In[ ]:
def checkDfamMatch(
    outputPositionLookupDic, outputMatchPositionList, readFileName, bucketNum=10
):
    dfamDataset = pd.read_csv(
        f"Evaluation/Source/{readFileName}.txt",
        sep="\t",
        header=0,
        usecols=["familyAcc", "familyName", "startIdx", "endIdx"],
    )
    dfamDatasetMatchList = [False] * len(dfamDataset)
    eachBucketNum = int(len(outputMatchPositionList) / bucketNum)
    for idx in range(len(dfamDataset)):
        row = dfamDataset.iloc[idx]
        start, end = row.start, row.end
        flag = 0
        for bucketIdx in range(bucketNum):
            if (
                start >= outputPositionLookupDic[bucketIdx][0]
                and start <= outputPositionLookupDic[bucketIdx][1]
            ):
                for outputRow in outputMatchPositionList[
                    eachBucketNum * bucketIdx : eachBucketNum * (bucketIdx + 1) - 1
                ]:
                    if outputRow.startIdx in range(
                        start, end
                    ) or outputRow.endIdx in range(start, end):
                        flag = 1
                        dfamDatasetMatchList[idx] = True
                        break
            if flag:
                break
    return dfamDatasetMatchList


# In[ ]:
# check output position with TE in Dfam
def checkOutputMatch(
    outputMatchPositionList, dfamPositionList, dfamPositionLookupDic, bucketNum=10
):
    outputMatchList = [False] * len(outputMatchPositionList)
    eachBucketNum = int(len(dfamPositionList) / bucketNum)
    for idx, outputPosition in enumerate(outputMatchPositionList):
        start, end = outputPosition[0], outputPosition[1]
        flag = 0
        for bucketIdx in range(bucketNum):
            if (
                start >= dfamPositionLookupDic[bucketIdx][0]
                and start <= dfamPositionLookupDic[bucketIdx][1]
            ):
                for dfamPosition in dfamPositionList[
                    eachBucketNum * bucketIdx : eachBucketNum * (bucketIdx + 1) - 1
                ]:
                    if dfamPosition.startIdx in range(
                        start, end
                    ) or dfamPosition.endIdx in range(start, end):
                        flag = 1
                        outputMatchList[idx] = True
                        break
            if flag == 1:
                break
    return outputMatchList
