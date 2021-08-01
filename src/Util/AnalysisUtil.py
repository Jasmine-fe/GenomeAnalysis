import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import Counter
from prettytable import PrettyTable


def listToSortedCounter(inputList):
    """
    In order to convert data to plot formate > ordered by length
    """
    inputCounter = Counter(inputList)
    mostCommonInputCounterList = inputCounter.most_common(
        len(inputCounter)
    )  # [(length, count)...], ordered by count
    sortedCounterList = sorted(
        mostCommonInputCounterList, key=lambda x: x[0]
    )  # [(length, count)...], ordered by length
    return sortedCounterList


def getStatisticData(inputList):
    inputSeries = pd.Series(inputList)
    print(inputSeries.describe())
    return 0


def mostCommonTable(mostCommonList, num):
    table = PrettyTable(["FragLen", "Count"])
    for i in range(num):
        fragmentLenValue = mostCommonList[i][0]
        count = mostCommonList[i][1]
        table.add_row([fragmentLenValue, count])
    print(table)
