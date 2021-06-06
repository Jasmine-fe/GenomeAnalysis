#!/usr/bin/env python
# coding: utf-8

# In[1]:


from collections import namedtuple


# In[3]:


SeqRepeatInfo = namedtuple('SeqRepeatInfo', ['fragmentLenList', 'count', 'position'])
PositionInfo = namedtuple('PositionInfo', ['chrIdx', 'fragmentIdx', 'baseIdx', 'seq'])
RepeatEvaInfo = namedtuple('RepeatEvaInfo', ['score', 'length', 'mismatchRatio'])
TandemRepeatInfo = namedtuple('TandemRepeatInfo', ['fragmentLenSubset', 'position'])
DfamRefSeqInfo = namedtuple('DfamRefSeqInfo', ['id', 'length', 'seq'])

