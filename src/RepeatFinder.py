#!/usr/bin/env python
# coding: utf-8
# In[ ]:
import os
import re
import sys
import time
import copy
import timeit
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statistics

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaIterator
from itertools import combinations
from collections import Counter, namedtuple
from prettytable import PrettyTable
import importlib

# importlib.reload(sys.modules["DataStructure"])

import import_ipynb
from Util.SeqUtil import parseSeq
from DataInfo import (
    currDataset,
    datasetPath,
    matchPattern,
    cutter,
    cutterLen,
    fragmentN,
    commonCount,
)
from DataStructure import (
    RepeatFragNInfo,
    IRSPositionInfo,
    TRSPositionInfo,
    RepeatEvaInfo,
)

# In[ ]:
"""
Return
repeatFragNLenList: repeat seq (if fragNLenCounter > 1)
repeatFragNPositionDict: {(46, 221, 139, 482, 101): [(0, 1), (0, 176)],...} // (a,b) chr(a+1), start from (b+1) th fragments
"""


def findRepeatSeqs(fragmentsLenList):
    print(f"... start finding repeat seq ...")
    start = time.time()
    repeatFragNPositionDict = {}
    fragNLenCounter = Counter()
    for index, i in enumerate(fragmentsLenList):
        for j in range(len(i) - fragmentN + 1):
            seqCombination = tuple(x for x in i[j : j + fragmentN])
            fragNLenCounter[seqCombination] += 1
            if seqCombination in repeatFragNPositionDict:
                repeatFragNPositionDict.get(seqCombination).append(tuple([index, j]))
            else:
                repeatFragNPositionDict[seqCombination] = [tuple([index, j])]
    end = time.time()
    print(f"...cost{ end-start } sec to finding repeat seq  ...")
    repeatFragNLenList = [(k, v) for k, v in fragNLenCounter.items() if v > 1]
    return repeatFragNLenList, repeatFragNPositionDict


# In[ ]:
def commonRepeatFragLenTable(commonRepeatFragLenCounter, repeatFragNPositionDict):
    # compareTable = PrettyTable(['FragmentLen', 'Count', 'Position'])
    table = PrettyTable(["FragLen", "Count"])
    for i in range(commonCount):
        fragmentLenValue = commonRepeatFragLenCounter[i][0]
        count = commonRepeatFragLenCounter[i][1]
        # position = repeatFragNPositionDict[fragmentLenValue]
        # table.add_row([fragmentLenValue, count, position])
        table.add_row([fragmentLenValue, count])
    print(table)


# In[ ]:
"""
Params
repeatType: 1-Tandem Repeat Sequence, 2-Interspersed Repeat Sequence
Return 
repeatInfoList: RepeatFragNInfo [(fragmentLenList, count, position: TRSPositionInfo / IRSPositionInfo (chrIdx, fragmentIdx, baseIdx, seq / seqList)), ... ]  
"""


def integrateRepeatInfo(
    fragmentsSeqList,
    fragmentsLenList,
    repeatFragNLenList,
    repeatFragNPositionDict,
    repeatType=1,
):
    repeatInfoList = []
    listLen = len(repeatFragNLenList)
    for i in range(listLen):
        positionList = []
        lenData = repeatFragNLenList[i][0]
        countData = repeatFragNLenList[i][1]
        targetPositionDic = repeatFragNPositionDict[lenData]
        for j in range(len(targetPositionDic)):
            chrIdx = targetPositionDic[j][0]
            fragmentIdx = targetPositionDic[j][1]
            baseIdx = sum(fragmentsLenList[chrIdx][:fragmentIdx]) + (
                cutterLen * fragmentIdx
            )
            # repeat type, TRS or IRS
            if repeatType == 1:
                seq = fragmentsSeqList[chrIdx][fragmentIdx : fragmentIdx + fragmentN]
                positionList.append(TRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seq))
            elif repeatType == 2:
                seqList = Seq("").join(
                    fragmentsSeqList[chrIdx][fragmentIdx : fragmentIdx + fragmentN]
                )
                positionList.append(
                    IRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seqList)
                )
        repeatInfoList.append(RepeatFragNInfo(lenData, countData, positionList))
    return repeatInfoList


# In[ ]:
def getIRComb(repeatInfoList):
    seqComb = []
    for i in range(len(repeatInfoList)):
        perm = combinations(repeatInfoList[i].position, 2)
        seqComb.append(list(perm))
    return seqComb


# In[ ]:
def evaluateRepeat(seq1, seq2, match=1, mismatch=-0.5):
    score = 0
    mismatchCount = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 == base2 and base1 != "-":
            score += match
        elif base1 != base2 and base1 != "-" and base2 != "-":
            score += mismatch
            mismatchCount += 1
    seqLength = max(len(seq1), len(seq2))
    mismatchRatio = round((mismatchCount / seqLength), 2)
    return RepeatEvaInfo(score, seqLength, mismatchRatio)


# In[ ]:
#
def longestRepeatLenInN(fragNLenList):
    fragmentsLenCounter = Counter(fragNLenList)
    mostCommonFragment = fragmentsLenCounter.most_common(1)
    mostCommonLen, mostCommonCount = mostCommonFragment[0][0], mostCommonFragment[0][1]
    repeatFragmentLen = [mostCommonLen]
    matchIndices = [
        idx for idx, value in enumerate(fragNLenList) if value == repeatFragmentLen[0]
    ]
    repeatFragmentIndices = copy.copy(matchIndices)
    # Udate repeatFragmentLen List
    for i in range(int(len(fragNLenList) / 2)):
        if matchIndices[-1] + 1 < len(
            fragNLenList
        ):  # check not the last element in the list
            flag = 1
            value = fragNLenList[matchIndices[0] + 1]
            for i in matchIndices:
                if fragNLenList[i + 1] != value:
                    flag = -1
            if flag == 1:
                matchIndices = [i + 1 for i in matchIndices]
                repeatFragmentLen = repeatFragmentLen + [value]
            else:
                break
    return repeatFragmentLen, repeatFragmentIndices


# In[ ]:
def checkTandemRepeatExist(repeatSeq):
    fragNLenList = repeatSeq[0]
    print("fragNLenList", len(fragNLenList))
    return len(fragNLenList) != len(set(fragNLenList))


# In[ ]:
def findRepeatInFragmentN(fragmentsLenList):
    fragmentsLenCounter = Counter(fragmentsLenList)
    filteredLenCounter = {
        length: count for length, count in fragmentsLenCounter.items() if count > 1
    }
    repeatFragmentDict = dict()
    for repeatLength, count in filteredLenCounter.items():
        matchIndices = [
            idx for idx, value in enumerate(fragmentsLenList) if repeatLength == value
        ]
        repeatFragmentDict[repeatLength] = matchIndices
    return repeatFragmentDict


# %%
def generateFragmentOutputFile(
    tandemRepeatInfoList, outputFileName="outputFragment", matchRatioOfSum=0.6
):
    filePath = os.path.join(os.getcwd()) + f"/../outputFile/{outputFileName}.txt"
    with open(filePath, "w") as outputFile:
        for i in range(len(tandemRepeatInfoList)):
            repeatFragmentDict = findRepeatInFragmentN(
                tandemRepeatInfoList[i].fragmentLenList
            )
            for repeatFragmentLen, repeatFragmentIndices in repeatFragmentDict.items():
                fragmentNum = 1
                idxComb = list(combinations(repeatFragmentIndices, 2))
                for seqPosition in tandemRepeatInfoList[i].position:
                    for k in idxComb:
                        seq1 = Seq("").join(seqPosition.seq[k[0] : k[0] + fragmentNum])
                        seq2 = Seq("").join(seqPosition.seq[k[1] : k[1] + fragmentNum])
                        seq1BaseSum = seqPosition.baseIdx + sum(
                            [len(i) for i in seqPosition.seq[: k[0]]]
                        )
                        seq2BaseSum = seqPosition.baseIdx + sum(
                            [len(i) for i in seqPosition.seq[: k[1]]]
                        )
                        repeatEvaInfo = evaluateRepeat(seq1, seq2)
                        if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                            output = f"""score:{repeatEvaInfo.score}, length:{repeatFragmentLen}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seqPosition.chrIdx}, {seq1BaseSum}) {seq1}\nSeq2:({seqPosition.chrIdx}, {seq2BaseSum}) {seq2}\n\n"""
                            outputFile.write(output)


# In[2]:
def generateTROutputFile(
    tandemRepeatInfoList, outputFileName="outputTR", matchRatioOfSum=0.6
):
    filePath = os.path.join(os.getcwd()) + f"/../outputFile/{outputFileName}.txt"
    with open(filePath, "w") as outputFile:
        for i in range(len(tandemRepeatInfoList)):
            repeatFragmentLen, repeatFragmentIndices = longestRepeatLenInN(
                tandemRepeatInfoList[i].fragmentLenList
            )
            fragmentNum = len(repeatFragmentLen)
            idxComb = list(combinations(repeatFragmentIndices, 2))
            for seqPosition in tandemRepeatInfoList[i].position:
                for k in idxComb:
                    seq1 = Seq("").join(seqPosition.seqList[k[0] : k[0] + fragmentNum])
                    seq2 = Seq("").join(seqPosition.seqList[k[1] : k[1] + fragmentNum])
                    seq1BaseSum = seqPosition.baseIdx + sum(
                        [len(i) for i in seqPosition.seqList[: k[0]]]
                    )
                    seq2BaseSum = seqPosition.baseIdx + sum(
                        [len(i) for i in seqPosition.seqList[: k[1]]]
                    )
                    repeatEvaInfo = evaluateRepeat(seq1, seq2)
                    if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                        output = f"""score:{repeatEvaInfo.score}, length:{sum(repeatFragmentLen)}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seqPosition.chrIdx}, {seq1BaseSum}) {seq1}\nSeq2:({seqPosition.chrIdx}, {seq2BaseSum}) {seq2}\n\n"""
                        outputFile.write(output)


# In[ ]:


def generateIROutputFile(
    seqPermutation, outputFileName="outputIR", matchRatioOfSum=0.6
):
    filePath = os.path.join(os.getcwd()) + f"/../outputFile/{outputFileName}.txt"
    with open(filePath, "w") as outputFile:
        for i in range(len(seqPermutation)):
            for j in range(len(seqPermutation[i])):
                seq1 = seqPermutation[i][j][0]
                seq2 = seqPermutation[i][j][1]
                if len(seq1.seq) != 0:
                    repeatEvaInfo = evaluateRepeat(seq1.seq, seq2.seq)
                    if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                        output = f"""score:{repeatEvaInfo.score}, length:{repeatEvaInfo.length}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seq1.chrIdx}, {seq1.baseIdx}) {seq1.seq}\nSeq2:({seq2.chrIdx}, {seq2.baseIdx}) {seq2.seq}\n\n"""
                        outputFile.write(output)
