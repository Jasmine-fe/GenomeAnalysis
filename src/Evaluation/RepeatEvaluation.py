#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import os
import math
import numpy as np
import pandas as pd
from collections import namedtuple
from DataStructure import DfamConSeqInfo, PositionInfo, DfamPositionInfo
from DataInfo import (
    cutterLen,
    fragmentN,
)


class RepeatEvaluation:
    def __init__(self, repeatInfoList):
        self.repeatInfoList = repeatInfoList
        self.repeatPositionList = []
        self.repeatPositionLookupDic = dict()

    def getRepeatPositionList(self):
        for repeatN in self.repeatInfoList:
            repeatFragNLen = sum(repeatN.fragmentLenList) + cutterLen * (fragmentN - 1)
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

    def positionBucketClassifier(self, bucketNum=10):
        eachBucketNum = int(len(self.repeatPositionList) / bucketNum)
        for i in range(bucketNum):
            firstStartIdx, lastStartIdx = (
                self.repeatPositionList[eachBucketNum * i].startIdx,
                self.repeatPositionList[(eachBucketNum * (i + 1) - 1)].startIdx,
            )
            self.repeatPositionLookupDic[i] = (firstStartIdx, lastStartIdx)
        # handle eachBucketNum remainder
        self.repeatPositionLookupDic[bucketNum - 1] = (
            self.repeatPositionLookupDic[bucketNum - 1][0],
            self.repeatPositionList[-1].startIdx,
        )
        return self.repeatPositionLookupDic

    # check output position with TE in Dfam
    def checkOutputMatch(self, dfamPositionList, dfamPositionLookupDic, bucketNum=10):
        outputMatchList = [False] * len(self.repeatPositionList)
        eachBucketNum = int(len(dfamPositionList) / bucketNum)
        matchedFamilyAccList, matchedFamilyNameList = [], []
        for idx, repeatPosition in enumerate(self.repeatPositionList):
            start, end = repeatPosition[0], repeatPosition[1]
            flag = 0
            for bucketIdx in range(bucketNum):
                if (
                    dfamPositionLookupDic[bucketIdx][0]
                    <= start
                    <= dfamPositionLookupDic[bucketIdx][1]
                ):
                    for dfamPosition in dfamPositionList[
                        eachBucketNum * bucketIdx : eachBucketNum * (bucketIdx + 1) - 1
                    ]:
                        if dfamPosition.startIdx in range(
                            start, end
                        ) or dfamPosition.endIdx in range(start, end + 100):
                            matchedFamilyAccList.append(dfamPosition.familyAcc)
                            matchedFamilyNameList.append(dfamPosition.familyName)
                            flag = 1
                            outputMatchList[idx] = True
                            break
                if flag == 1:
                    break
        return outputMatchList, matchedFamilyAccList, matchedFamilyNameList
