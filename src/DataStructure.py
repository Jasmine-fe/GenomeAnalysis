#!/usr/bin/env python
# coding: utf-8

# In[1]:
from collections import namedtuple

RepeatFragNInfo = namedtuple(
    "RepeatFragNInfo", ["fragmentLenList", "count", "position"]
)

# *PositionInfo, for position in RepeatFragNInfo
IRSPositionInfo = namedtuple(
    "IRSPositionInfo", ["chrIdx", "fragmentIdx", "baseIdx", "seq"]
)
TRSPositionInfo = namedtuple(
    "IRSPositionInfo", ["chrIdx", "fragmentIdx", "baseIdx", "seqList"]
)

RepeatEvaInfo = namedtuple("RepeatEvaInfo", ["score", "length", "mismatchRatio"])

DfamConSeqInfo = namedtuple("DfamConSeqInfo", ["id", "length", "seq"])

PositionInfo = namedtuple("PositionInfo", ["length", "startIdx", "endIdx"])

DfamPositionInfo = namedtuple(
    "PositionInfo", ["familyAcc", "familyName", "startIdx", "endIdx"]
)

refSeqSimilarityInfo = namedtuple(
    "refSeqSimilarityInfo", ["hitId", "targetSeqLength", "similarityPercentage"]
)