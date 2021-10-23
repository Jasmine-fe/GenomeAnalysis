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
        self.repeatPositionTable = None
        self.cutter = cutter
        self.parseFastaSeqs = None
        self.seqStateList = None

    def parseFasta(self):
        self.parseFastaSeqs = parseFasta(
            currDataset, datasetPath, matchPattern, matchMode=False
        )
        self.initSeqStateList()
        return self.parseFastaSeqs

    def initSeqStateList(self):
        self.seqStateList = [0] * len(self.parseFastaSeqs[0])

    def parseSeqByCutter(self):
        self.fragmentLenList, self.fragmentSeqList = parseSeqByCutter(
            self.parseFastaSeqs, cutter=self.cutter
        )
        return self.fragmentLenList, self.fragmentSeqList

    def findRepeatSeqs(self, lengthLimit=False):
        repeatFragNLenList, repeatFragNPositionDict = findRepeatSeqs(
            self.fragmentLenList
        )
        temrepeatInfoList = integrateRepeatInfo(
            self.cutter,
            self.fragmentSeqList,
            self.fragmentLenList,
            repeatFragNLenList,
            repeatFragNPositionDict,
            repeatType=2,
        )
        self.repeatInfoList = (
            self.seqLengthLimit(temrepeatInfoList) if lengthLimit else temrepeatInfoList
        )
        return self.repeatInfoList

    def seqLengthLimit(self, temrepeatInfoList):
        lowerBound, upperBound = 23, 2000
        lengthLimitList = []
        for i in temrepeatInfoList:
            seqLen = i.fragmentLenList[0]  # For N = 1
            if seqLen <= upperBound:
                lengthLimitList.append(i)
        return lengthLimitList

    def getRepeatPositionList(self, filter=True):  # [ modify filter ]
        """
        Return
        repeatPositionList: [(startIdx, endIdx), ...]
        """
        seqList = self.filterRepeatInfoList if filter else self.repeatInfoList
        cutterLen = len(self.cutter)
        for repeatN in seqList:
            repeatFragNLen = sum(repeatN.fragmentLenList) + cutterLen * (fragmentN + 1)
            for fragposition in repeatN.position:
                self.repeatPositionList.append(
                    PositionInfo(
                        repeatFragNLen,
                        fragposition.baseIdx,
                        fragposition.baseIdx + repeatFragNLen + 1,
                    )
                )
        self.repeatPositionList.sort(key=lambda x: x[0])
        self.repeatPositionTable = pd.DataFrame(self.repeatPositionList)
        self.getSeqFromPositionIdx()
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

    def getSeqFromPositionIdx(self):
        self.repeatPositionTable["seq"] = [
            str(self.parseFastaSeqs[0][row["startIdx"] : row["endIdx"]])
            for idx, row in self.repeatPositionTable.iterrows()
        ]
        return self.repeatPositionTable

    def getRepeatFragentInfo(self):
        frgmentCount = len(self.repeatPositionTable)
        lengthCount = self.repeatPositionTable["length"].nunique()
        return lengthCount, frgmentCount

    def generateRepeatFragentFile(self, filePath="tem.txt"):
        lengthCount, frgmentCount = self.getRepeatFragentInfo()
        uniqueLengthList = set(self.repeatPositionTable["length"])
        with open(filePath, "w") as outputFile:
            outputFile.write(
                f"Info:\nlengthCount:{lengthCount}, frgmentCount:{frgmentCount}\n"
            )
            for length in uniqueLengthList:
                output = f"{length}\n"
                targetPositionTable = self.repeatPositionTable[
                    self.repeatPositionTable["length"] == length
                ]
                for key, row in targetPositionTable.iterrows():
                    output += (
                        f'({row["startIdx"]}, {row["endIdx"]})'
                        + "\n"
                        + row["seq"]
                        + "\n"
                    )
                outputFile.write(output + "\n")

    def calculateDiversityRatio(self, repeatPositionTable):
        """
        repeatPositionTable (length , seq)
        """
        fragmentGroupByLenCountDf = (
            repeatPositionTable.groupby("length", as_index=False)
            .agg({"seq": "count"})
            .rename(columns={"seq": "lengthSeqCount"})
        )
        fragmentGroupBySeqAndLenKeys = (
            repeatPositionTable.groupby(["length", "seq"], as_index=False)
            .agg({"seq": "count"})
            .rename(columns={"seq": "seqTypeCount"})
        )
        mergedFragmentGroupByDf = fragmentGroupBySeqAndLenKeys.merge(
            fragmentGroupByLenCountDf, how="left", on="length"
        )
        mergedFragmentGroupByDf["diversityRatio"] = (
            mergedFragmentGroupByDf["seqTypeCount"]
            * mergedFragmentGroupByDf["seqTypeCount"]
        ) / (
            mergedFragmentGroupByDf["lengthSeqCount"]
            * mergedFragmentGroupByDf["lengthSeqCount"]
        )
        diversityRatioDf = mergedFragmentGroupByDf.groupby(
            "length", as_index=False
        ).agg({"diversityRatio": "sum"})
        return diversityRatioDf

    def seqStateGenerator(self):
        """
        state - 0 unmatch, 1 union, 2 intersection
        """
        for position in self.repeatPositionList:
            for base in range(position.startIdx, position.endIdx):
                # for base in range(
                #     position.startIdx - len(self.cutter), position.endIdx + len(self.cutter)
                # ):
                self.seqStateList[base] += 1
        self.seqStateList = list(map(lambda x: 1 if x > 1 else x, self.seqStateList))
        return self.seqStateList