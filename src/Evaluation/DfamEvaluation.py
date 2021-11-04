import re
import math
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from SharedInfo import chrPattern, humanChrKey
from Evaluation.RepeatEvaluation import RepeatEvaluation
from DataStructure import (
    DfamPositionInfo,
    refSeqSimilarityInfo,
)


class DfamEvaluation(RepeatEvaluation):
    def __init__(self, repeatPositionList, hitFileName="chrX_dm6_dfam.nrph.hits"):
        super().__init__(repeatPositionList)
        self.dfamPositionList = []
        self.dfamPositionLookupDic = dict()
        self.familyPositionList = []
        self.familySeqSimilarityList = []
        self.dfamBucketAmount = 10
        self.dfamEachBucketNum = 0
        self.hitFileName = hitFileName
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

    def checkDfamMatchWithRepeat(self):
        """
        How many Dfam intersect with repeat position
        """
        dfamMatchList = [False] * len(self.dfamPositionList)
        matchedFamilyAccList, matchedFamilyNameList = [], []

        for idx, ref in enumerate(self.dfamPositionList):
            refStart, refEnd = ref.startIdx, ref.endIdx
            flag = 0
            for repeatPosition in self.repeatPositionList:
                if (
                    refStart in range(repeatPosition.startIdx, repeatPosition.endIdx)
                ) or (refEnd in range(repeatPosition.startIdx, repeatPosition.endIdx)):
                    matchedFamilyAccList.append(ref.familyAcc)
                    matchedFamilyNameList.append(ref.familyName)
                    flag = 1
                    dfamMatchList[idx] = True
                    break
        return dfamMatchList, matchedFamilyAccList, matchedFamilyNameList

    def checkRepeatMatchWithDfam(self):
        """
        How many repeats intersect with Dfam position
        """
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

    def matchRatio(self, repeatMatchList):
        totalLen = len(repeatMatchList)
        matchLen = len(list(filter(lambda x: x, repeatMatchList)))
        ratio = matchLen / totalLen
        print(f"matchCount:{matchLen}\tdfamCount:{totalLen}\tRatio:{ratio}")
        return ratio

    def familyMatchRatio(self, matchedFamilyAccList):
        matchSet = set(matchedFamilyAccList)
        dfamSet = set(list(self.dfamHitData["family_acc"]))
        matchSetLen, dfamSetLen = len(matchSet), len(dfamSet)
        ratio = round(matchSetLen / dfamSetLen, 2)
        print("Family Match Result:")
        print(f"matchCount:{matchSetLen}\tdfamCount:{dfamSetLen}\tRatio:{ratio}")
        print(f"unmatch acc: {dfamSet.difference(matchSet)}")
        return ratio

    def getUnmatchInfo(self, DRrepeatMatchList):
        unMatchDf = pd.DataFrame(columns=["index", "length", "startIdx", "endIdx"])
        for idx, value in enumerate(DRrepeatMatchList):
            if value == False:
                row = self.dfamPositionList[idx]
                LRTLength = row.endIdx - row.startIdx
                unMatchDf = unMatchDf.append(
                    {
                        "index": idx,
                        "length": LRTLength,
                        "startIdx": row.startIdx,
                        "endIdx": row.endIdx,
                    },
                    ignore_index=True,
                )
        # print statistic info
        unMatchLenSeries = pd.Series(list(unMatchDf["length"]))
        print(unMatchLenSeries.describe())
        return unMatchDf

    def generateFullLengthLTRFile(
        self,
        inputFileName="chrX_LTR_dm6_dfam.nrph.hits",
        outputFileName="chrX_FullLength_LTR_dm6_dfam.hit",
    ):
        df = pd.read_csv(f"Evaluation/Source/{inputFileName}", sep="\t", header=0)
        idx = 0
        fullLengthLTRIdxList = []
        while idx < len(df) - 1:
            if df.iloc[idx]["family_acc"] == df.iloc[idx + 1]["family_acc"]:
                fullLengthLTRIdxList.append([idx, idx + 1])
                idx += 2
            else:
                idx += 1
        flattenFullLengthLTRIdxList = [
            item for sublist in fullLengthLTRIdxList for item in sublist
        ]
        fullLengthDf = df.loc[
            flattenFullLengthLTRIdxList,
        ]
        fullLengthDf.to_csv(
            f"Evaluation/Source/{outputFileName}",
            sep="\t",
            encoding="utf-8",
            index=False,
        )

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

    def consensusSeqSimilarity(self, consensusSeq, seqList):
        seqSimilarityList = []
        for idx, targetSeq in enumerate(seqList):
            alignments = pairwise2.align.localxx(consensusSeq, targetSeq)
            targetLength = len(targetSeq)
            similarityPercentage = (
                round(alignments[0].score / targetLength, 2)
                if len(alignments) > 0
                else 0
            )
            seqSimilarityList.append(
                refSeqSimilarityInfo(
                    hitId=idx,
                    targetSeqLength=targetLength,
                    similarityPercentage=similarityPercentage,
                )
            )
        return seqSimilarityList
