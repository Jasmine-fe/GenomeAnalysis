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

# importlib.reload(sys.modules['Evaluation'])


class Sequence:
    def __init__(self, cutter):
        self.fragmentLenList = []
        self.fragmentSeqList = []
        self.repeatInfoList = []
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

    def findRepeatSeqs(self):
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
        return self.repeatInfoList

    def getRepeatPositionList(self):
        """
        Return
        repeatPositionList: [(startIdx, endIdx), ...]
        """
        cutterLen = len(self.cutter)
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
