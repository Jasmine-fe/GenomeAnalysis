#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import os
import math
import numpy as np
import pandas as pd
from Bio import pairwise2
from collections import namedtuple, Counter
from DataStructure import DfamConSeqInfo, PositionInfo, RepeatFragNInfo
from SharedInfo import fragmentN


class RepeatEvaluation:
    def __init__(self, repeatPositionList):
        self.repeatPositionList = repeatPositionList
        self.repeatPositionLookupDic = dict()
        self.bucketAmount = 10
        self.eachBucketNum = 0

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