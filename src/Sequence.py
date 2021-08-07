import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, pairwise2
from prettytable import PrettyTable
from collections import Counter, namedtuple

import importlib
from DataStructure import PositionInfo
from Util.SeqUtil import seqInfo, parseFasta, parseSeqByCutter
from Evaluation.DfamEvaluation import DfamEvaluation

from SharedInfo import (
    currDataset,
    datasetPath,
    matchPattern,
    fragmentN,
)
from RepeatFinder import (
    findRepeatSeqs,
    integrateRepeatInfo,
)

from DataStructure import RepeatFragNInfo

# importlib.reload(sys.modules['Evaluation'])


class Sequence:
    def __init__(self, cutter):
        self.fragmentLenList = []
        self.fragmentSeqList = []
        self.repeatInfoList = []
        self.filterRepeatInfoList = []
        self.repeatPositionList = []
        self.cutter = cutter
        self.parseFastaSeqs = None

    def parseFasta(self):
        self.parseFastaSeqs = parseFasta(
            currDataset, datasetPath, matchPattern, matchMode=False
        )
        return self.parseFastaSeqs

    def parseSeqByCutter(self):
        self.fragmentLenList, self.fragmentSeqList = parseSeqByCutter(
            self.parseFastaSeqs, cutter=self.cutter
        )
        return self.fragmentLenList, self.fragmentSeqList

    def findRepeatSeqs(self, lengthLimit=False):
        repeatFragNLenList, repeatFragNPositionDict = findRepeatSeqs(
            self.fragmentLenList
        )
        self.repeatInfoList = integrateRepeatInfo(
            self.cutter,
            self.fragmentSeqList,
            self.fragmentLenList,
            repeatFragNLenList,
            repeatFragNPositionDict,
            repeatType=2,
        )
        if lengthLimit:
            self.repeatInfoList = self.seqLengthLimit()
        return self.repeatInfoList

    def seqLengthLimit(self):
        lowerBound, upperBound = 23, 1500
        lengthLimitList = []
        for i in self.repeatInfoList:
            seqLen = i.fragmentLenList[0]  # For N = 1
            if seqLen >= lowerBound and seqLen <= upperBound:
                lengthLimitList.append(i)
        self.repeatInfoList = lengthLimitList
        return self.repeatInfoList

    def getRepeatPositionList(self, filter=True):
        """
        Return
        repeatPositionList: [(startIdx, endIdx), ...]
        """
        seqList = self.filterRepeatInfoList if filter else self.repeatInfoList
        cutterLen = len(self.cutter)
        for repeatN in seqList:
            repeatFragNLen = sum(repeatN.fragmentLenList) + cutterLen * fragmentN
            for fragposition in repeatN.position:
                self.repeatPositionList.append(
                    PositionInfo(
                        fragposition.baseIdx, fragposition.baseIdx + repeatFragNLen
                    )
                )
        self.repeatPositionList.sort(key=lambda x: x[0])
        return self.repeatPositionList

    def filterRepeatInfo(self):
        for i in self.repeatInfoList:
            filterPosition = self.filterSeqPosition(i.position)
            if len(filterPosition) > 1:
                positionList = filterPosition
                count = len(filterPosition)
                self.filterRepeatInfoList.append(
                    RepeatFragNInfo(i.fragmentLenList, count, positionList)
                )
        return self.filterRepeatInfoList

    def filterSeqPosition(self, positionList):
        filterPositionList = []
        seqCounter = Counter([i.seq for i in positionList])
        repeatSeqs = [x for x, count in seqCounter.items() if count > 1]
        for i in positionList:
            if i.seq in repeatSeqs:
                filterPositionList.append(i)
        return filterPositionList
