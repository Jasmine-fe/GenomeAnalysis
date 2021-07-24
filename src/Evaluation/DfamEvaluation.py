import re
import math
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from DataInfo import chrPattern, humanChrKey
from Evaluation.RepeatEvaluation import RepeatEvaluation
from DataStructure import (
    DfamConSeqInfo,
    PositionInfo,
    DfamPositionInfo,
    refSeqSimilarityInfo,
)


class DfamEvaluation(RepeatEvaluation):
    def __init__(self, repeatInfoList):
        super().__init__(repeatInfoList)
        hitFileName = "chrX_dm6_dfam.nrph.hits"
        self.dfamPositionList = []
        self.dfamPositionLookupDic = dict()
        self.familyPositionList = []
        self.familySeqSimilarityList = []
        self.dfamBucketAmount = 10
        self.dfamEachBucketNum = 0
        self.dfamHitData = pd.read_csv(
            f"Evaluation/Source/{hitFileName}", sep="\t", header=0
        )

    def getDfamPositionList(self):
        """
        get DfamPositionInfo format data for .hit file
        """
        for idx, refRow in self.dfamHitData.iterrows():
            startIdx = refRow["ali-st"] if refRow["strand"] == "+" else refRow["ali-en"]
            endIdx = refRow["ali-en"] if refRow["strand"] == "+" else refRow["ali-st"]
            self.dfamPositionList.append(
                DfamPositionInfo(
                    refRow["family_acc"], refRow["family_name"], startIdx, endIdx
                )
            )
        self.dfamPositionList.sort(
            key=lambda tup: tup[2]
        )  # choose by key name (startIdx)
        return self.dfamPositionList

    # def dfamPositionBucketClassifier(self, bucketNum=10):
    #     self.dfamEachBucketNum = int(len(self.dfamPositionList) / bucketNum)
    #     for i in range(bucketNum):
    #         firstStartIdx, lastStartIdx = (
    #             self.dfamPositionList[self.dfamEachBucketNum * i].startIdx,
    #             self.dfamPositionList[(self.dfamEachBucketNum * (i + 1) - 1)].startIdx,
    #         )
    #         self.dfamPositionLookupDic[i] = (firstStartIdx, lastStartIdx)
    #     # handle eachBucketNum remainder
    #     self.dfamPositionLookupDic[bucketNum - 1] = (
    #         self.dfamPositionLookupDic[bucketNum - 1][0],
    #         self.dfamPositionList[-1].startIdx,
    #     )
    #     return self.dfamPositionLookupDic

    # check other match than 0, 1
    def checkRepeatMatchWithDfam(self):
        repeatMatchList = [False] * len(self.repeatPositionList)
        matchedFamilyAccList, matchedFamilyNameList = [], []
        for idx, repeatPosition in enumerate(self.repeatPositionList):
            repeatStart, repeatEnd = repeatPosition[0], repeatPosition[1]
            flag = 0
            for bucketIdx in range(self.bucketAmount):
                if (
                    self.repeatPositionLookupDic[bucketIdx][0]
                    <= repeatStart
                    <= self.repeatPositionLookupDic[bucketIdx][1]
                ):
                    for ref in self.dfamPositionList:
                        if repeatStart in range(
                            ref.startIdx, ref.endIdx
                        ) or repeatEnd in range(ref.startIdx, ref.endIdx):
                            matchedFamilyAccList.append(ref.familyAcc)
                            matchedFamilyNameList.append(ref.familyName)
                            flag = 1
                            repeatMatchList[idx] = True
                            break
                if flag == 1:
                    break
        return repeatMatchList, matchedFamilyAccList, matchedFamilyNameList

    # ----------- For specific family --------------------------------------

    def getfamilyPositionList(self, familyName):
        familyPositionIterator = filter(
            lambda x: x["familyAcc"] == familyName, self.dfamPositionList
        )
        self.familyPositionList = list(familyPositionIterator)
        return self.familyPositionList

    def getfamilySeqList(self, familyName, seq):
        self.familySeqList = []
        familyPositionList = self.getfamilyPositionList(familyName)
        for row in familyPositionList:
            self.familySeqList.append(seq[row["startIdx"] : row["endIdx"]])
        return self.familySeqList

    def refSeqSimilarity(self, consensusSeq, familySeqList):
        familySeqSimilarityList = []
        for idx, targetSeq in enumerate(familySeqList):
            alignments = pairwise2.align.localxx(consensusSeq, targetSeq)
            targetLength = len(targetSeq)
            similarityPercentage = (
                round(alignments[0].score / targetLength, 2)
                if len(alignments) > 0
                else 0
            )
            familySeqSimilarityList.append(
                refSeqSimilarityInfo(
                    hitId=idx,
                    targetSeqLength=targetLength,
                    similarityPercentage=similarityPercentage,
                )
            )
        return familySeqSimilarityList
