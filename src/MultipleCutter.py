import re
import sys
import pandas as pd
from SharedInfo import cutterA, cutterB
from Util.SeqUtil import parseSeqByCutter
from DataStructure import RepeatSeqOutputInfo


class MultipleCutter:
    def __init__(self, sequence, seqStateList):
        self.seqStates = {"unMatch": 0, "union": 1, "intersection": 2}
        self.sequence = sequence
        self.seqLength = len(sequence)
        self.seqStateList = seqStateList
        self.seqStateSum = [0] * len(sequence)
        self.matchStateIdxList = []
        self.matchStateRepeatInfoList = []
        self.repeatFrgmentDf = None

    def getSeqStateSum(self):
        if len(self.seqStateList) == 2:
            for idx, value in enumerate(self.seqStateSum):
                self.seqStateSum[idx] = (
                    self.seqStateList[0][idx] + self.seqStateList[1][idx]
                )
        else:
            print("Number of seqState Error")
        return self.seqStateSum

    def getSeqStateInfo(self):
        unMatchState = filter(
            lambda x: x == self.seqStates["unMatch"], self.seqStateSum
        )
        unionState = filter(
            lambda x: (x == self.seqStates["union"])
            or (x == self.seqStates["intersection"]),
            self.seqStateSum,
        )
        intersectionState = filter(
            lambda x: x == self.seqStates["intersection"], self.seqStateSum
        )
        print(
            f"chr: {self.seqLength}\nunMatch: {len(list(unMatchState))}, union:{len(list(unionState))}, intersection:{len(list(intersectionState))}"
        )
        return unMatchState, unionState, intersectionState

    def getSpecificStateIdxList(self, stateName="union"):
        self.matchStateIdxList = []
        if stateName == "intersection":
            self.matchStateIdxList = [
                idx
                for idx, state in enumerate(self.seqStateSum)
                if state == self.seqStates[stateName]
            ]
        elif stateName == "union":
            self.matchStateIdxList = [
                idx
                for idx, state in enumerate(self.seqStateSum)
                if (
                    (state == self.seqStates[stateName])
                    or (state == self.seqStates["intersection"])
                )
            ]
        return self.matchStateIdxList

    def getSpecificStatePositionList(self):
        """
        idx+baseIdx < leng(matchStateRepeatInfoList)
        idx : 400
        baseIdx = 1
        baseIdx == (matchStateRepeatInfoList[idx+baseIdx] - matchStateRepeatInfoList[idx])
        """

        self.matchStateRepeatInfoList = []
        idx = 0
        while self.matchStateIdxList[idx] < self.seqLength:
            baseCount = 1
            while (idx + baseCount < len(self.matchStateIdxList)) and (
                (self.matchStateIdxList[idx + baseCount] - self.matchStateIdxList[idx])
                == baseCount
            ):
                baseCount += 1
            else:
                startIdx = self.matchStateIdxList[idx]
                endIdx = self.matchStateIdxList[idx + baseCount - 1]
                self.matchStateRepeatInfoList.append(
                    RepeatSeqOutputInfo(
                        (
                            self.matchStateIdxList[idx + baseCount - 1]
                            - self.matchStateIdxList[idx]
                        ),
                        startIdx,
                        endIdx,
                        self.sequence[startIdx:endIdx],
                    )
                )
                if idx + baseCount >= len(self.matchStateIdxList):
                    break
            idx = idx + baseCount
        return self.matchStateRepeatInfoList

    def generateSeqStateSumFile(self, filePath):
        with open(filePath, "w") as f:
            f.write("".join(str(state) for state in self.seqStateSum))

    def cutRepeatSeqToFragment(self):
        self.matchStateRepeatInfoList.sort(key=lambda x: x.length, reverse=True)
        totalRepeat = pd.DataFrame(columns=["length", "startIdx", "endIdx", "seq"])
        for repeatInfo in self.matchStateRepeatInfoList[:1]:
            if repeatInfo.seq[: len(cutterA)] == cutterA:
                fragmentsLenList, fragmentsSeqList = parseSeqByCutter(
                    [repeatInfo.seq], cutter=cutterA
                )
                df = self.cauculateSeqPosition(
                    repeatInfo, len(cutterA), fragmentsLenList[0], fragmentsSeqList[0]
                )
            elif repeatInfo.seq[: len(cutterB)] == cutterB:
                fragmentsLenList, fragmentsSeqList = parseSeqByCutter(
                    [repeatInfo.seq], cutter=cutterB
                )
                df = self.cauculateSeqPosition(
                    repeatInfo, len(cutterB), fragmentsLenList[0], fragmentsSeqList[0]
                )
            else:
                print("Error")

            if (len(fragmentsLenList[0]) > 0) and (len(fragmentsSeqList[0]) > 0):
                totalRepeat = totalRepeat.append(df, ignore_index=True)
        self.repeatFrgmentDf = totalRepeat.loc[df["length"] != 0]
        return self.repeatFrgmentDf

    def cauculateSeqPosition(
        self, repeatInfo, cutterLen, fragmentsLenList, fragmentsSeqList
    ):
        df = pd.DataFrame(columns=["length", "startIdx", "endIdx", "seq"])
        startIdx = repeatInfo.startIdx
        for idx, fragmentLength in enumerate(fragmentsLenList):
            start = startIdx + sum(fragmentsLenList[:idx]) + cutterLen * idx
            end = start + fragmentLength
            df = df.append(
                {
                    "length": fragmentLength,
                    "startIdx": start,
                    "endIdx": end,
                    "seq": fragmentsSeqList[idx],
                },
                ignore_index=True,
            )
        return df

    def fragmentGroupbyLen(self):
        matchDfGroupByLen = self.repeatFrgmentDf.groupby(by=["length"], sort=True)
        temDf = self.repeatFrgmentDf.groupby(by=["length"]).agg({"length": "sum"})
        original_stdout = sys.stdout
        with open(f"../outputFile/SeqState/groupByLenData.txt", "w") as f:
            sys.stdout = f
            for key, row in temDf.iterrows():
                print(f"{key}:")
                for i in matchDfGroupByLen.get_group(key).index:
                    print(
                        f"({ self.repeatFrgmentDf.iloc[i]['startIdx']}, {self.repeatFrgmentDf.iloc[i]['endIdx']})\n{ self.repeatFrgmentDf.iloc[i]['seq']}"
                    )
                print("\n")
            sys.stdout = original_stdout