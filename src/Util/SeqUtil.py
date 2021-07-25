import re
import time
import timeit
import random
import statistics
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter, namedtuple
from contextlib import redirect_stdout
from DataInfo import (
    currDataset,
    datasetPath,
    matchPattern,
    cutter,
    commonCount,
)


def parseFasta(currDataset, datasetPath, matchPattern, matchMode=False):
    print(f"...start parsing {currDataset} fasta file ...")
    start = time.time()
    fasta_sequences = SeqIO.parse(open(datasetPath), "fasta")
    seqs = (
        matchMode
        and [i.seq for i in fasta_sequences if re.search(matchPattern, i.id)]
        or [i.seq for i in fasta_sequences]
    )
    end = time.time()
    print(f"...cost{ end-start } sec to parse fasta file ...")
    return seqs


def parseSeq(seq, cutter):
    targetNum = seq.count(cutter)
    parseResult = seq.split(cutter)
    fragmentLength = [len(read) for read in parseResult]
    return targetNum, parseResult, fragmentLength


def parseSeqByCutter(parseFastaSeqs, cutter=cutter):
    start = time.time()
    print(f"...start parse seq by cutter: { cutter }")
    fragmentsLenList = []
    fragmentsSeqList = []
    for seq in parseFastaSeqs:
        seq = seq.upper()
        targetNum, fragmentSeq, fragmentLength = parseSeq(seq, cutter)
        fragmentsLenList.append(fragmentLength)
        fragmentsSeqList.append(fragmentSeq)
    end = time.time()
    print(f"...cost { end-start } sec to cut sequence")
    return fragmentsLenList, fragmentsSeqList


def readSeqId(fileName):
    fasta_sequences = SeqIO.parse(open(fileName), "fasta")
    seqIds = [i.id for i in fasta_sequences]
    return seqIds


def seqInfo(dataName, seqs):
    seqNum = len(seqs)
    sumLen = sum([len(i) for i in seqs])
    print(f"{dataName} dataset\n number of sequence:{seqNum}\n total length:{sumLen}\n")


def randomSeqGenerator(length):
    seq = "".join(random.choice("ACTG") for _ in range(length))
    return Seq(seq)


def seqFileGenerator(parseFastaSeq):
    """
    change parseFastaSeq [startIdx: endIdx]
    """
    with open("unMatchSeqCompare_10.txt", "w") as f:
        with redirect_stdout(f):
            print(f"Before Dfam RefSeq")
            print(f"seq10\t{parseFastaSeq[551477-600:551477]}")
            print(f"seq11\t{parseFastaSeq[558499-600:558499]}\n")

            print(f"Dfam Hit RefSeq")
            print(f"seq10\t{parseFastaSeq[551477:551961]}")
            print(f"seq11\t{parseFastaSeq[558499:558984]}\n")
