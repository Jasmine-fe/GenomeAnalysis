#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import os
import re
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
# importlib.reload(sys.modules['DataStructure'])

import import_ipynb
from Util.SeqUtil import parseSeq
from DataInfo import currDataset, datasetPath, matchPattern, cutter, cutterLen, fragmentN, commonCount
from DataStructure import SeqRepeatInfo, IRSPositionInfo, TRSPositionInfo, RepeatEvaInfo
# In[ ]:
'''
Return
allRepeatSeqCount: repeat seq (if repeatCount > 1)
commonRepeatSeq: top fragmentN repeat seq
Parameter
allRepeatSeqType: 0 dict, 1 list of tuple
'''
def findRepeatSeqs(fragmentsLenList, cutter , fragmentN, commonCount, allRepeatSeqType = 1):
    cutterLen = len(cutter)
    repeatCount, positionDic = findRepeatSeq(fragmentsLenList, fragmentN , cutterLen)
    allRepeatSeqCount = allRepeatSeqType and [(k, v) for k, v in repeatCount.items() if v > 1] or {k: v for k, v in repeatCount.items() if v > 1}
    commonRepeatSeq = repeatCount.most_common(commonCount)
    return allRepeatSeqCount, commonRepeatSeq, repeatCount, positionDic


# In[ ]:
def commonRepeatTable(commonRepeatCounter, positionDic, commonCount):
    # compareTable = PrettyTable(['FragmentLen', 'Count', 'Position'])
    compareTable = PrettyTable(['FragmentLen', 'Count'])
    for i in range(commonCount):
        fragmentLenValue = commonRepeatCounter[i][0]
        count = commonRepeatCounter[i][1]
        position = positionDic[fragmentLenValue]
    #     compareTable.add_row([fragmentLenValue, count, position])
        compareTable.add_row([fragmentLenValue, count])
    print(compareTable)


# In[ ]:
def findRepeatSeq(lenList, N, cutterLen):
    print(f'... start finding repeat seq ...')
    start = time.time()
    positionDic = {}
    repeatCount = Counter()
    preIndex = 0
    currentIndex = 0
    for index, i in enumerate(lenList):
        for j in range(len(i) - N + 1):
            seqCombination = tuple(x for x in i[j:j+N])
            repeatCount[seqCombination] += 1
            if seqCombination in positionDic:  
                positionDic.get(seqCombination).append(tuple([index, j]))   
            else:
                positionDic[seqCombination] = [tuple([index, j])]
            
    end = time.time() 
    print(f'...cost{ end-start } sec to finding repeat seq  ...')
    return repeatCount, positionDic


# In[ ]:
def commonRepeatSeqsPosition(fragmentsLenList, commonRepeatInfo, positionDic, commonCount):
    positionList = []
    for i in range(commonCount):
        targetFragmentLenList = commonRepeatInfo[i][0]
        targetFragmentDic = positionDic[targetFragmentLenList]
        for j in range(len(targetFragmentDic)):
            chrIndex = targetFragmentDic[j][0]
            fragmentIndex = targetFragmentDic[j][1]
            positionValue = sum(fragmentsLenList[chrIndex][:fragmentIndex])+(cutterLen*fragmentIndex)
            positionList.append(tuple([chrIndex, positionValue]))
    positionList.sort(key=lambda tup: tup[0])
    return positionList


# In[ ]:
# Repeat Position {chr i: startIndex1, startIndex2, startIndex3 ...}
def commonPositionListToDic(commonPositionList):
    commonPositionDic = dict()
    [ commonPositionDic [t[0]].append(t [1]) if t [0] in list(commonPositionDic.keys()) else commonPositionDic.update({t [0]: [t [1]]}) for t in commonPositionList ]
    return commonPositionDic


# In[ ]:
def printCommonPositionDic(commonPositionDic):
    keys = list(commonPositionDic.keys())
    for i in range(len(keys)):
        commonPositionDic[keys[i]].sort()
        print(f"\n{keys[i]}:\n{','.join( str(v) for v in commonPositionDic[keys[i]] )}\n")


# In[ ]:
def topRepeatSeqCounter(fragmentsSeqList, commonRepeatInfo, positionDic, fragmentN):
    topRepeatValue = commonRepeatInfo[0][0]
    topRepeatCount = commonRepeatInfo[0][1]
    print(f"the top repeat sequence: {topRepeatValue},{topRepeatCount}")
    topRepeatPosition = positionDic[topRepeatValue]
    topOriginalRepeatSeqs = []
    topCombinedRepeatSeqs = []
    for i in range(len(topRepeatPosition)):
        chrIndex, fragmentIndex = topRepeatPosition[i][0], topRepeatPosition[i][1]
        originalSeq = fragmentsSeqList[chrIndex][fragmentIndex]
        combinedSeq = Seq('').join(fragmentsSeqList[chrIndex][fragmentIndex:fragmentIndex+fragmentN])
        topOriginalRepeatSeqs.append(originalSeq)
        topCombinedRepeatSeqs.append(combinedSeq)
    singleRepeatSeqCount = Counter(topOriginalRepeatSeqs)
    comRepeatSeqCount = Counter(topCombinedRepeatSeqs)
    return singleRepeatSeqCount, comRepeatSeqCount


# In[ ]:
def printTopSeqPercentageTable(repeatSeqCount, showPartialSeq=True):
    repeatTable = PrettyTable(['Seq', 'Count', 'Percentage(%)'])
    sumCount = sum(repeatSeqCount.values())
    commonSeqs = repeatSeqCount.most_common(commonCount)
    for i in range(commonCount):
        currentSeq = showPartialSeq and commonSeqs[i][0][:10] or commonSeqs[i][0]
        currentCount = commonSeqs[i][1]
        currentPro = currentCount*100 / sumCount
        repeatTable.add_row([currentSeq, currentCount, currentPro])
    print(repeatTable)

# In[ ]:
'''
Return 
repeatInfoList: SeqRepeatInfo [(fragmentLenList, count, position: TRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seqList)), ... ]  
'''
def integrateTandemRepeatInfo(fragmentsSeqList, fragmentsLenList, allRepeatSeqCount, positionDic):
    repeatInfoList = []
    listLen = len(allRepeatSeqCount)
    for i in range(listLen):
        positionList = []
        lenData = allRepeatSeqCount[i][0]
        countData = allRepeatSeqCount[i][1]
        targetPositionDic = positionDic[lenData]
        for j in range(len(targetPositionDic)):
            chrIdx = targetPositionDic[j][0]
            fragmentIdx = targetPositionDic[j][1]
            baseIdx = sum(fragmentsLenList[chrIdx][:fragmentIdx])+(cutterLen*fragmentIdx)
            seq = fragmentsSeqList[chrIdx][fragmentIdx:fragmentIdx+fragmentN]
            positionList.append(TRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seq))
        repeatInfoList.append(SeqRepeatInfo(lenData, countData, positionList))
    return repeatInfoList
# In[ ]:
'''
Return 
repeatInfoList: SeqRepeatInfo [(fragmentLenList, count, position: IRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seq)), ... ]  
'''
def integrateRepeatInfo(fragmentsSeqList, fragmentsLenList, allRepeatSeqCount, positionDic):
    repeatInfoList = []
    listLen = len(allRepeatSeqCount)
    for i in range(listLen):
        positionList = []
        lenData = allRepeatSeqCount[i][0]
        countData = allRepeatSeqCount[i][1]
        targetPositionDic = positionDic[lenData]
        for j in range(len(targetPositionDic)):
            chrIdx = targetPositionDic[j][0]
            fragmentIdx = targetPositionDic[j][1]
            baseIdx = sum(fragmentsLenList[chrIdx][:fragmentIdx])+(cutterLen*fragmentIdx)
            seqList = Seq('').join(fragmentsSeqList[chrIdx][fragmentIdx:fragmentIdx+fragmentN])
            positionList.append(IRSPositionInfo(chrIdx, fragmentIdx, baseIdx, seqList))
        repeatInfoList.append(SeqRepeatInfo(lenData, countData, positionList))
    return repeatInfoList

# In[ ]:
def getIRComb(repeatInfoList):
    seqComb = []
    for i in range (len(repeatInfoList)):
        perm = combinations(repeatInfoList[i].position, 2) 
        seqComb.append(list(perm))
    return seqComb


# In[ ]:


def evaluateRepeat(seq1, seq2, match=1, mismatch=-0.5):
    score = 0
    mismatchCount = 0
    for base1, base2 in zip(seq1, seq2):
            if base1 == base2 and base1 != '-':
                score += match
            elif base1 != base2 and base1 != '-' and base2 != '-':
                score += mismatch
                mismatchCount += 1
    seqLength = max(len(seq1), len(seq2))
    mismatchRatio = round((mismatchCount/seqLength), 2)
    return RepeatEvaInfo(score ,seqLength, mismatchRatio) 


# In[ ]:


def generateIROutputFile(seqPermutation, matchRatioOfSum = 0.6):
    filePath = os.path.join(os.getcwd()) + '/../outputFile/outputIR.txt'
    with open( filePath, 'w') as outputFile:
        for i in range(len(seqPermutation)):
            for j in range(len(seqPermutation[i])):
                seq1 = seqPermutation[i][j][0] 
                seq2 = seqPermutation[i][j][1] 
                repeatEvaInfo = evaluateRepeat(seq1.seq , seq2.seq)
                if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                    output = f"""score:{repeatEvaInfo.score}, length:{repeatEvaInfo.length}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seq1.chrIdx}, {seq1.baseIdx}) {seq1.seq}\nSeq2:({seq2.chrIdx}, {seq2.baseIdx}) {seq2.seq}\n\n"""
                    outputFile.write(output)


# In[ ]:
def longestCommonLength(fragmentsLenList):
    fragmentsLenCounter = Counter(fragmentsLenList)
    mostCommonFragment = fragmentsLenCounter.most_common(1)
    mostCommonLen , mostCommonCount = mostCommonFragment[0][0], mostCommonFragment[0][1]
    repeatFragmentLen = [mostCommonLen]
    matchIndices = [ idx for idx, value in enumerate(fragmentsLenList) if value == repeatFragmentLen[0]]
    repeatFragmentIndices = copy.copy(matchIndices)
    # Udate repeatFragmentLen List
    for i in range(int(len(fragmentsLenList)/2)):
        if(matchIndices[-1]+1 < len(fragmentsLenList)): #check not the last element in the list
            flag = 1 
            value = fragmentsLenList[matchIndices[0]+1]
            for i in matchIndices:
                if fragmentsLenList[i+1] != value:
                    flag = -1
            if flag == 1:
                matchIndices = [i+1 for i in matchIndices]
                repeatFragmentLen = repeatFragmentLen + [value]
            else:
                break
    return repeatFragmentLen, repeatFragmentIndices


# In[ ]:
def checkTandemRepeatExist(repeatSeq):
    fragmentlengths = repeatSeq[0]
    count = repeatSeq[1]
    return len(fragmentlengths) != len(set(fragmentlengths))


# In[2]:

def generateTROutputFile(tandemRepeatInfoList, matchRatioOfSum= 0.6):
    filePath = os.path.join(os.getcwd()) + '/../outputFile/outputTR.txt'
    with open( filePath, 'w') as outputFile:
        for i in range(len(tandemRepeatInfoList)):
            repeatFragmentLen, repeatFragmentIndices = longestCommonLength(tandemRepeatInfoList[i].fragmentLenList)
            fragmentNum = len(repeatFragmentLen)
            idxComb = list(combinations(repeatFragmentIndices, 2))
            for seqPosition in tandemRepeatInfoList[i].position:
                for k in idxComb:
                    seq1 = Seq('').join(seqPosition.seqList[k[0]:k[0]+fragmentNum])
                    seq2 = Seq('').join(seqPosition.seqList[k[1]:k[1]+fragmentNum])
                    seq1BaseSum = seqPosition.baseIdx + sum([len(i) for i in seqPosition.seqList[:k[0]]])
                    seq2BaseSum = seqPosition.baseIdx + sum([len(i) for i in seqPosition.seqList[:k[1]]])
                    repeatEvaInfo = evaluateRepeat(seq1 , seq2)
                    if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                        output = f"""score:{repeatEvaInfo.score}, length:{sum(repeatFragmentLen)}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seqPosition.chrIdx}, {seq1BaseSum}) {seq1}\nSeq2:({seqPosition.chrIdx}, {seq2BaseSum}) {seq2}\n\n"""
                        outputFile.write(output)

# In[ ]:
def findRepeatInFragmentN(fragmentsLenList):
    fragmentsLenCounter = Counter(fragmentsLenList)
    filteredLenCounter =  { length: count for length, count in fragmentsLenCounter.items() if count > 1 }
    repeatFragmentDict = dict()
    for repeatLength, count in filteredLenCounter.items():
        matchIndices = [ idx for idx, value in enumerate(fragmentsLenList) if repeatLength == value ]
        repeatFragmentDict[repeatLength] = matchIndices
    return repeatFragmentDict 

# %%

def generateFragmentOutputFile(tandemRepeatInfoList, matchRatioOfSum= 0.6):
    filePath = os.path.join(os.getcwd()) + '/../outputFile/outputFragment.txt'
    with open(filePath, 'w') as outputFile:
        for i in range(len(tandemRepeatInfoList)):
            repeatFragmentDict = findRepeatInFragmentN(tandemRepeatInfoList[i].fragmentLenList)
            for repeatFragmentLen, repeatFragmentIndices in repeatFragmentDict.items(): 
                fragmentNum = 1
                idxComb = list(combinations(repeatFragmentIndices, 2))
                for seqPosition in tandemRepeatInfoList[i].position:
                    for k in idxComb:
                        seq1 = Seq('').join(seqPosition.seq[k[0]:k[0]+fragmentNum])
                        seq2 = Seq('').join(seqPosition.seq[k[1]:k[1]+fragmentNum])
                        seq1BaseSum = seqPosition.baseIdx + sum([len(i) for i in seqPosition.seq[:k[0]]])
                        seq2BaseSum = seqPosition.baseIdx + sum([len(i) for i in seqPosition.seq[:k[1]]])
                        repeatEvaInfo = evaluateRepeat(seq1 , seq2)
                        if repeatEvaInfo.score > repeatEvaInfo.length * matchRatioOfSum:
                            output = f"""score:{repeatEvaInfo.score}, length:{repeatFragmentLen}, mismatch ratio:{repeatEvaInfo.mismatchRatio}\nSeq1:({seqPosition.chrIdx}, {seq1BaseSum}) {seq1}\nSeq2:({seqPosition.chrIdx}, {seq2BaseSum}) {seq2}\n\n"""
                            outputFile.write(output)