#!/usr/bin/env python
# coding: utf-8
# In[1]:
import re
import time
import timeit
import random
import statistics
import numpy as np 
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from collections import Counter, namedtuple
from IPython.display import clear_output # clear_output(wait=True)

import import_ipynb
from DataInfo import currDataset, datasetPath, matchPattern, cutter, cutterLen, fragmentN, commonCount
# In[1]:
def parseFasta(currDataset, datasetPath, matchPattern, matchMode = False):
    print(f'...start parsing {currDataset} fasta file ...')
    start = time.time()
    fasta_sequences = SeqIO.parse(open(datasetPath), "fasta")
    seqs = matchMode and [ i.seq for i in fasta_sequences if re.search(matchPattern, i.id) ] or [ i.seq for i in fasta_sequences ]
    end = time.time()
    print(f'...cost{ end-start } sec to parse fasta file ...')
    return seqs


# In[ ]:
def parseSeq(seq, cutter):
    targetNum = seq.count(cutter)
    parseResult = seq.split(cutter)
    fragmentLength = [len(read) for read in parseResult]   
    return targetNum, parseResult, fragmentLength


# In[ ]:
def parseSeqByCutter(parseFastaSeqs): 
    start = time.time()
    print(f'...start parse seq by cutter: { cutter }')
    fragmentsLenList = []
    fragmentsSeqList = []
    for seq in parseFastaSeqs:
        targetNum, fragmentSeq, fragmentLength = parseSeq(seq, cutter)
        fragmentsLenList.append(fragmentLength)
        fragmentsSeqList.append(fragmentSeq)
    end = time.time()
    print(f'...cost { end-start } sec to cut sequence')
    return fragmentsLenList, fragmentsSeqList


# In[ ]:
def readSeqId(fileName):
    fasta_sequences = SeqIO.parse(open(fileName), "fasta")
    seqIds = [i.id for i in fasta_sequences]
    return seqIds

# In[3]:
def seqInfo(dataName, seqs):
    seqNum = len(seqs)
    sumLen = sum([len(i) for i in seqs])
    print(f'{dataName} dataset\n number of sequence:{seqNum}\n total length:{sumLen}\n')

# In[4]:
def randomSeqGenerator(length):
    seq = ''.join(random.choice('ACTG') for _ in range(length))
    return Seq(seq)
