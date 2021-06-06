#!/usr/bin/env python
# coding: utf-8

# In[1]:


import csv
import sys
import re
import time
import timeit
import random
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import statistics
get_ipython().run_line_magic('matplotlib', 'inline')

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaIterator
from collections import Counter
from itertools import islice, groupby
from collections import Counter, namedtuple
from prettytable import PrettyTable
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


# In[1]:


def exportCsvFile(fileName, rowList):
    with open(fileName, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['length','uniqueCount'])
        for row in rowList:
            writer.writerows(rowList)


# In[3]:


def seqInfo(dataName, seqs):
    seqNum = len(seqs)
    sumLen = sum([len(i) for i in seqs])
    print(f'{dataName} dataset\n number of sequence:{seqNum}\n total length:{sumLen}\n')


# In[4]:


def nGramLength(lengthList, count, N):
    nGramLegList = []
    for i in range(count-N+1):
        nGramLegList.append(sum(lengthList[i:i+N]))
    return nGramLegList


# In[4]:


def randomSeqGenerator(length):
    seq = ''.join(random.choice('ACTG') for _ in range(length))
    return Seq(a)


# In[ ]:


def unifragmentList(fragmentsLenList, fragmentsSeqList):
    flattenLenList = [item for subList in fragmentsLenList for item in subList]
    lenCounter = Counter(flattenLenList)
    flattenSeqList = [item for subList in fragmentsSeqList for item in subList]
    flattenSeqList.sort(key=len)
    uniqueLen = namedtuple('fragments',['length','uniqueCount'])
    uniqueLenList = []
    classifiedLenSeqs = ({k: list(v) for k, v in groupby(sorted(flattenSeqList, key=len), key=len)})
    for key, value in classifiedLenSeqs.items():
        tem = uniqueLen(key, len(set(value)))
        uniqueLenList.append(tem)
    sortedUniLenList = sorted(uniqueLenList, key=lambda x: x.uniqueCount, reverse=True)
    print(f"The most common length of fragments:\n {sortedUniLenList[:10]}")
    exportCsvFile("table.csv", sortedUniLenList[:100])
    return uniqueLenList


# In[ ]:


# fragment length distribution plot
def fragmentLenPlot(uniqueLenList):
    xInterval = 1
    x_axiMax = uniqueLenList[-1].length
    x = np.arange(0, x_axiMax, xInterval)
    y = [sum(map(lambda x: x.uniqueCount, uniqueLenList[j:j+xInterval])) for j in range(0, x_axiMax, xInterval) ]
    plt.figure(figsize=(20,12)) 
    plt.title("Length of unique fragments' Distribution")
    plt.xlabel("length", fontsize=20)
    plt.ylabel("occurrences", fontsize=20)
    plt.bar(x[:20000], y[:20000], color='orange')
    plt.xlim([0, 2000])
    plt.show()
    return 0

