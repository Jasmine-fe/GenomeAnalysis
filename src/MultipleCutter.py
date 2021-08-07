from DataStructure import PositionInfo


class MultipleCutter:
    def __init__(self, chrLength, repeatPositionList):
        self.seqStates = {"unMatch": 0, "union": 1, "intersection": 2}
        self.chrLength = chrLength
        self.repeatPositionList = repeatPositionList
        self.seqStateList = [0] * self.chrLength
        self.matchStateIdxList = []
        self.matchStatePositionList = []

    def seqStateGenerator(self):
        """
        state - 0 unmatch, 1 union, 2 intersection
        """
        for position in self.repeatPositionList:
            for base in range(position.startIdx, position.endIdx):
                self.seqStateList[base] += 1
        return self.seqStateList

    def getSeqStateInfo(self):
        unMatchState = filter(
            lambda x: x == self.seqStates["unMatch"], self.seqStateList
        )
        unionState = filter(lambda x: x == self.seqStates["union"], self.seqStateList)
        intersectionState = filter(
            lambda x: x == self.seqStates["intersection"], self.seqStateList
        )
        print(
            f"chr: {self.chrLength}\nunMatch: {len(list(unMatchState))}, union:{len(list(unionState))}, intersection:{len(list(intersectionState))}"
        )
        return unMatchState, unionState, intersectionState

    def getSpecificStateIdxList(self, stateName="union"):
        self.matchStateIdxList = [
            idx
            for idx, state in enumerate(self.seqStateList)
            if state == self.seqStates[stateName]
        ]
        return self.matchStateIdxList

    def getSpecificStatePositionList(self):
        self.matchStatePositionList = []
        idx = 0
        while idx < len(self.matchStateIdxList) - 1:
            baseCount = 1
            while (idx + baseCount < len(self.matchStateIdxList)) and (
                self.matchStateIdxList[idx + baseCount] - self.matchStateIdxList[idx]
            ) == baseCount:
                baseCount += 1
            else:
                if idx + baseCount >= len(self.matchStateIdxList):
                    break
                # elif (
                #     self.matchStateIdxList[idx]
                #     - self.matchStateIdxList[idx + baseCount]
                #     != baseCount
                # ):
                else:
                    self.matchStatePositionList.append(
                        PositionInfo(
                            self.matchStateIdxList[idx],
                            self.matchStateIdxList[idx + baseCount - 1],
                        )
                    )
                    idx = idx + baseCount
        return self.matchStatePositionList
