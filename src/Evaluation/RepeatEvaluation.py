#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import os
import math
import numpy as np
import pandas as pd
from Bio import pairwise2
from collections import namedtuple
from DataStructure import DfamConSeqInfo, PositionInfo
from DataInfo import (
    cutterLen,
    fragmentN,
)


class RepeatEvaluation:
    def __init__(self, repeatInfoList):
        self.repeatInfoList = repeatInfoList
        self.repeatPositionList = []
        self.repeatPositionLookupDic = dict()
        self.bucketAmount = 10
        self.eachBucketNum = 0

    def getRepeatPositionList(self):
        for repeatN in self.repeatInfoList:
            repeatFragNLen = sum(repeatN.fragmentLenList) + cutterLen * fragmentN
            for fragposition in repeatN.position:
                self.repeatPositionList.append(
                    PositionInfo(
                        fragposition.baseIdx, fragposition.baseIdx + repeatFragNLen
                    )
                )
        self.repeatPositionList.sort(key=lambda x: x[0])
        return self.repeatPositionList

    """
    Params
    positionList: [(startIdx, endIdx), ...]
    """

    def positionBucketClassifier(self):
        self.eachBucketNum = int(len(self.repeatPositionList) / self.bucketAmount)
        for i in range(self.bucketAmount):
            firstStartIdx, lastStartIdx = (
                self.repeatPositionList[self.eachBucketNum * i].startIdx,
                self.repeatPositionList[(self.eachBucketNum * (i + 1) - 1)].startIdx,
            )
            self.repeatPositionLookupDic[i] = (firstStartIdx, lastStartIdx)
        # handle eachBucketNum remainder
        self.repeatPositionLookupDic[self.bucketAmount - 1] = (
            self.repeatPositionLookupDic[self.bucketAmount - 1][0],
            self.repeatPositionList[-1].startIdx,
        )
        return self.repeatPositionLookupDic

    # check output position with TE in Dfam
    # def checkOutputMatch(self, dfamPositionList, dfamPositionLookupDic, bucketNum=10):
    #     outputMatchList = [False] * len(self.repeatPositionList)
    #     matchedFamilyAccList, matchedFamilyNameList = [], []
    #     for idx, repeatPosition in enumerate(self.repeatPositionList):
    #         start, end = repeatPosition[0], repeatPosition[1]
    #         flag = 0
    #         for bucketIdx in range(bucketNum):
    #             if (
    #                 dfamPositionLookupDic[bucketIdx][0]
    #                 <= start
    #                 <= dfamPositionLookupDic[bucketIdx][1]
    #             ):
    #                 for dfamPosition in dfamPositionList[
    #                     (self.eachBucketNum * bucketIdx) : (
    #                         self.eachBucketNum * (bucketIdx + 1) - 1
    #                     )
    #                 ]:
    #                     if dfamPosition.startIdx in range(
    #                         start, end
    #                     ) or dfamPosition.endIdx in range(start, end):
    #                         matchedFamilyAccList.append(dfamPosition.familyAcc)
    #                         matchedFamilyNameList.append(dfamPosition.familyName)
    #                         flag = 1
    #                         print("Match")
    #                         outputMatchList[idx] = True
    #                         break
    #             if flag == 1:
    #                 break
    #     return outputMatchList, matchedFamilyAccList, matchedFamilyNameList


# ---------------------------- For single result ---------------------------------


def printPairSeq(repeatInfoListEle):
    print("FragmentN: ", fragmentN)
    print("FragmentLenList: ", repeatInfoListEle.fragmentLenList)
    print("Count: ", repeatInfoListEle.count)
    print("Position: ")
    for idx, row in enumerate(repeatInfoListEle.position):
        print(f"{idx} - baseIdx:{row.baseIdx} , seq:{row.seq}")


def getAlignResult(repeatInfoListEle):
    alignments = pairwise2.align.localxx(
        repeatInfoListEle.position[0].seq, repeatInfoListEle.position[1].seq
    )
    print(
        "Alignment Result\t",
        f"length:{(alignments[0].end - alignments[0].start)}\t",
        "score:",
        alignments[0].score,
    )
    print("seqA:", alignments[0].seqA)
    print("seqA:", alignments[0].seqB)