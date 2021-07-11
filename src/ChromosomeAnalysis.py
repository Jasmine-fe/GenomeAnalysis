import numpy as np
import pandas as pd


def getChrBaseLengthDic():
    chrSizeFilePath = "../Source/ChrBaseLength.csv"
    chrSize = pd.read_csv(chrSizeFilePath, delimiter="\t")
    chrBaseLengthDic = {}
    for idx, row in chrSize.iterrows():
        key = chr + row["Chromosome"]
        value = int(row["BaseLength"].replace(",", ""))
        chrBaseLengthDic[key] = value
    return chrBaseLengthDic