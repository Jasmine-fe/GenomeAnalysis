import re
import math
import numpy as np
import pandas as pd
from DataStructure import DfamConSeqInfo, PositionInfo, DfamPositionInfo
from Evaluation.RepeatEvaluation import RepeatEvaluation


class DfamEvaluation(RepeatEvaluation):
    def __init__(self):
        self.repeatPositionList = super().repeatInfoList
        self.dfamPositionList = []

    def getDfamConSeqData(self):
        # init DfamConSeqInfo value
        dfId = ""
        seqLineCount = 1
        currLineCount = 0
        seqLength = 0
        seq = ""
        seqFlag = 0

        DfamRefDataset = open("Evaluation/Source/Dfam.embl.txt", "r")
        lines = DfamRefDataset.readlines()
        DfamRefSeqList = []
        for eachLine in lines:
            if re.search("^AC", eachLine):
                dfId = re.search("DF\d*", eachLine).group()
            elif re.search("^SQ", eachLine):
                seqFlag = 1
                currLineCount = 0
                seqLineCount = math.ceil(seqLength / 60)
                seqLength = int(re.search("\d+", eachLine).group())
            elif seqFlag == 1 and (currLineCount <= seqLineCount):
                matchList = re.findall("[a-zA-Z]+", eachLine)
                segmentSeq = "".join(matchList)
                seq += segmentSeq
                currLineCount += 1
                if currLineCount >= seqLineCount:
                    DfamRefSeqList.append(
                        DfamConSeqInfo(id=dfId, length=seqLength, seq=seq)
                    )
                    seq = ""
                    seqFlag = 0
        return DfamRefSeqList

    def generateDfamPositionData(
        self,
        readFileName,
        outputFileName,
    ):
        dfamDataset = pd.read_csv(
            f"Evaluation/Source/{readFileName}.hits",
            sep="\t",
            header=0,
            usecols=["familyAcc", "familyName", "aliSt", "aliEn", "envSt", "envEn"],
        )
        ouputFilePath = f"Evaluation/Source/{outputFileName}.txt"
        outputCols = ["familyAcc", "familyName", "startIdx", "endIdx"]
        with open(ouputFilePath, "w") as outputFile:
            outputFile.write("\t".join(outputCols) + "\n")
            for idx in range(len(dfamDataset)):
                row = dfamDataset.iloc[idx]
                li = sorted([row.aliSt, row.aliEn, row.envSt, row.envEn])
                start, end = li[0], li[-1]
                output = f"{row.familyAcc}\t{row.familyName}\t{start}\t{end}\n"
                outputFile.write(output)

    def getDfamPositionList(self, readFileName):
        dfamDataset = pd.read_csv(
            f"Evaluation/Source/{readFileName}.txt",
            sep="\t",
            header=0,
            usecols=[
                "familyAcc",
                "familyName",
                "startIdx",
                "endIdx",
            ],
        )
        for idx in range(len(dfamDataset)):
            row = dfamDataset.iloc[idx]
            self.dfamPositionList.append(
                DfamPositionInfo(
                    row.familyAcc, row.familyName, row.startIdx, row.endIdx
                )
            )
        self.dfamPositionList.sort(
            key=lambda tup: tup[2]
        )  # choose by key name (startIdx)
        return self.dfamPositionList

    def checkDfamMatch(
        self,
        readFileName,
        bucketNum=10,
    ):
        dfamDataset = pd.read_csv(
            f"Evaluation/Source/{readFileName}.txt",
            sep="\t",
            header=0,
            usecols=["familyAcc", "familyName", "startIdx", "endIdx"],
        )
        dfamDatasetMatchList = [False] * len(dfamDataset)
        eachBucketNum = int(len(self.repeatPositionList) / bucketNum)
        for idx in range(len(dfamDataset)):
            row = dfamDataset.iloc[idx]
            start, end = row.start, row.end
            flag = 0
            for bucketIdx in range(bucketNum):
                if (
                    start >= self.repeatPositionLookupDic[bucketIdx][0]
                    and start <= self.repeatPositionLookupDic[bucketIdx][1]
                ):
                    for outputRow in self.repeatPositionList[
                        eachBucketNum * bucketIdx : eachBucketNum * (bucketIdx + 1) - 1
                    ]:
                        if outputRow.startIdx in range(
                            start, end
                        ) or outputRow.endIdx in range(start, end):
                            flag = 1
                            dfamDatasetMatchList[idx] = True
                            break
                if flag:
                    break
        return dfamDatasetMatchList
