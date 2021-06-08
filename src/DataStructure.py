#!/usr/bin/env python
# coding: utf-8

# In[1]:
from collections import namedtuple

# In[3]:
RepeatFragNInfo = namedtuple(
    "RepeatFragNInfo", ["fragmentLenList", "count", "position"]
)
IRSPositionInfo = namedtuple(
    "IRSPositionInfo", ["chrIdx", "fragmentIdx", "baseIdx", "seq"]
)  # for position in RepeatFragNInfo
TRSPositionInfo = namedtuple(
    "IRSPositionInfo", ["chrIdx", "fragmentIdx", "baseIdx", "seqList"]
)

RepeatEvaInfo = namedtuple("RepeatEvaInfo", ["score", "length", "mismatchRatio"])
DfamRefSeqInfo = namedtuple("DfamRefSeqInfo", ["id", "length", "seq"])
