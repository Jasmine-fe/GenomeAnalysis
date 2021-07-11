import numpy as np
import pandas as pd


def getChrBaseLengthDic():
    chrSizeFilePath = "../Source/ChrBaseLength.csv"
    chrSize = pd.read_csv(chrSizeFilePath, delimiter="\t")
    chrBaseLengthDic = {}
    for idx, row in chrSize.iterrows():
        key = chr + row["Chromosome"]
        chrBaseLengthDic[key] = row["BaseLength"]
    return chrBaseLengthDic