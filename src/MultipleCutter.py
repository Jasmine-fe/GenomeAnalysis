import re
from DataStructure import PositionInfo


class MultipleCutter:
    def __init__(self, chrLength, seqStateList):
        self.seqStates = {"unMatch": 0, "union": 1, "intersection": 2}
        self.chrLength = chrLength
        self.seqStateList = seqStateList
        self.seqStateSum = [0] * chrLength
        self.matchStateIdxList = []
        self.matchStatePositionList = []

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
            f"chr: {self.chrLength}\nunMatch: {len(list(unMatchState))}, union:{len(list(unionState))}, intersection:{len(list(intersectionState))}"
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
        idx+baseIdx < leng(matchStatePositionList)
        idx : 400
        baseIdx = 1
        baseIdx == (matchStatePositionList[idx+baseIdx] - matchStatePositionList[idx])
        """

        self.matchStatePositionList = []
        idx = 0
        while self.matchStateIdxList[idx] < self.chrLength:
            baseCount = 1
            while (idx + baseCount < len(self.matchStateIdxList)) and (
                (self.matchStateIdxList[idx + baseCount] - self.matchStateIdxList[idx])
                == baseCount
            ):
                baseCount += 1
            else:
                self.matchStatePositionList.append(
                    PositionInfo(
                        self.matchStateIdxList[idx],
                        self.matchStateIdxList[idx + baseCount - 1],
                    )
                )
                if idx + baseCount >= len(self.matchStateIdxList):
                    break
            idx = idx + baseCount
        return self.matchStatePositionList

    def generateSeqStateSumFile(self, filePath):
        with open(filePath, "w") as f:
            f.write("".join(str(state) for state in self.seqStateSum))
