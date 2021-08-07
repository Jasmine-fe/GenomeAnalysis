# # tandem Repeat in N fragment
# tandemRepeatLenList = list(filter(checkTandemRepeatExist, repeatFragNLenList))
# tandemRepeatInfoList = integrateRepeatInfo(fragmentsSeqList, fragmentsLenList, tandemRepeatLenList, repeatFragNPositionDict, repeatType=1)
# generateTROutputFile(tandemRepeatInfoList, outputFileName= f"{currDatasetName}TRS", matchRatioOfSum=0.4)

# interspersed repetitive sequences
# repeatInfoList = integrateRepeatInfo(fragmentsSeqList, fragmentsLenList, repeatFragNLenList, repeatFragNPositionDict, repeatType=2)
# filterRepeatInfoList = filterRepeatInfo(repeatInfoList)
# seqPermutation = getIRComb(repeatInfoList)
# repeatEvaInfoList = generateIROutputFile(seqPermutation, outputFileName= f"{currDatasetName}_IRS", matchRatioOfSum=0.8)

# evaluation = RepeatEvaluation(repeatInfoList)
# repeatPositionList = evaluation.getRepeatPositionList()
# repeatPositionLookupDic = evaluation.positionBucketClassifier()
# outputMatchList, matchedFamilyAccList, matchedFamilyNameList = evaluation.checkOutputMatch(dfamPositionList, dfamPositionLookupDic, bucketNum= 10)

# fragmentsLenList, fragmentsSeqList = parseSeqByCutter(parseFastaSeqs)
# repeatFragNLenList, repeatFragNPositionDict = findRepeatSeqs(fragmentsLenList)
# interspersed repetitive sequences
# repeatInfoList = integrateRepeatInfo(fragmentsSeqList, fragmentsLenList, repeatFragNLenList, repeatFragNPositionDict, repeatType=2)
# filterRepeatInfoList = filterRepeatInfo(repeatInfoList)
# seqPermutation = getIRComb(repeatInfoList)
# repeatEvaInfoList = generateIROutputFile(seqPermutation, outputFileName= f"{currDatasetName}_IRS", matchRatioOfSum=0.8)

# ROO_LTR_df = pd.DataFrame(columns=['seq', 'seqLength'])
# for idx, row in enumerate(dfamPositionList):
#     seq = str(parseFastaA[0][row.startIdx: row.endIdx])
#     seqLength = row.endIdx - row.startIdx
#     ROO_LTR_df = ROO_LTR_df.append({'seq': seq, 'seqLength': seqLength}, ignore_index=True)

# Similarity plot
# def consensusSeqSimilarity(consensusSeq, seqDf):
#     print("hihi", len(seqDf))
#     seqSimilarityList = []
#     for targetSeq in seqDf:
#         alignments = pairwise2.align.localxx(consensusSeq, targetSeq)
#         targetLength = len(targetSeq)
#         similarityPercentage = round(alignments[0].score / targetLength, 2)
#         seqSimilarityList.append(similarityPercentage)
#     return seqSimilarityList

# xData = [*range(len(seqSimilarityList))]
# fig, ax = plt.subplots(figsize=(12, 5), dpi=300)
# sns.set_style("whitegrid")
# sns.stripplot(x=xData, y=seqSimilarityList, linewidth=1.0, color=colorA)
# ax.set_ylabel("Similarity", size=15)
# ax.set_ylim(0,1.2)

# ROO_LTR_df = pd.read_csv("./Evaluation/Source/ROO_LTR_Seq.csv")
# conParseFasta = parseFasta(
#     "DF0001696_ROO_LTR",
#     "./Evaluation/Source/DF0001696_ROO_LTR.fa",
#     "*",
#     matchMode=False,
# )
# consensusSeq = conParseFasta[0].upper()
# seqDf = ROO_LTR_df["seq"]
# seqSimilarityList = consensusSeqSimilarity(consensusSeq, seqDf)