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