from re import T
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from SharedInfo import colorA, colorB


def basicPlot(
    sortedCounterList,
    title="Distribution",
    xlabel="Length",
    ylabel="Count",
    xlimUpperBound=3000,
):
    df = pd.DataFrame(columns=["x", "y"], dtype=float)
    for row in sortedCounterList:
        df = df.append({"x": row[0], "y": row[1]}, ignore_index=True)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    plt.yticks(np.arange(0, 6, step=1))
    sns.set_style("whitegrid")
    sns.lineplot(data=df, x="x", y="y")
    ax.set_xlabel(xlabel, size=15)
    ax.set_ylabel(ylabel, size=15)

    xlimUpperBound and ax.set_xlim(0, xlimUpperBound)


def twoLabelBasicPlot(
    sortedCounterListA,
    sortedCounterListB,
    title="Distribution",
    xlabel="length",
    ylabel="count",
    xlimUpperBound=1000,
    labelA="",
    labelB="",
):
    xA = [i[0] for i in sortedCounterListA]
    yA = [i[1] for i in sortedCounterListA]
    xB = [i[0] for i in sortedCounterListB]
    yB = [i[1] for i in sortedCounterListB]
    fig, ax = plt.subplots()
    ax.plot(xA, yA, color="tab:blue", label=labelA)
    ax.plot(xB, yB, color="tab:orange", label=labelB)
    ax.set_title(title, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlim(0, xlimUpperBound)
    ax.legend()


def lengthScatterDistributionPlot(xList):
    fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
    sns.set_style("whitegrid")
    sns.stripplot(x=xList, linewidth=1.0, color=colorA)
    ax.set_xlabel("Length", size=15)