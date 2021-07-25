import re
import os
import sys
import numpy as np
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt
from Bio import SeqIO, pairwise2
from prettytable import PrettyTable
from collections import Counter, OrderedDict, namedtuple

import importlib
from Util.FragmentUtil import fragmentLenPlot
from Util.SeqUtil import seqInfo, parseFasta, parseSeq, parseSeqByCutter
from DataStructure import RepeatFragNInfo
from DataInfo import (
    currDataset,
    datasetPath,
    matchPattern,
    commonCount,
    chrPattern,
)
from RepeatFinder import (
    findRepeatSeqs,
    integrateRepeatInfo,
    getIRComb,
    evaluateRepeat,
    generateIROutputFile,
    checkTandemRepeatExist,
    generateTROutputFile,
    generateFragmentOutputFile,
    filterRepeatInfo,
    commonRepeatFragLenTable,
)
from Evaluation.DfamEvaluation import DfamEvaluation

# importlib.reload(sys.modules['Evaluation'])


class Sequence:
    def __init__(self, cutter):
        self.fragmentLenList = []
        self.fragmentSeqList = []
        self.repeatInfoList = []
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
            self.fragmentSeqList,
            self.fragmentLenList,
            repeatFragNLenList,
            repeatFragNPositionDict,
            repeatType=2,
        )
        return self.repeatInfoList