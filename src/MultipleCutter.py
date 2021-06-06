cutterA = "GATC"
cutterB = "AAGCTT"

# split seq by cutter
fragmentLenListA = []
fragmentSeqListA = []
fragmentLenListB = []
fragmentSeqListB = []

for seq in parseFastaSeqs:
    targetNumA, fragmentSeqA, fragmentLengthA = parseSeq(seq, cutterA)
    fragmentLenListA.append(fragmentLengthA)
    fragmentSeqListA.append(fragmentSeqA)
    
    targetNumB, fragmentSeqB, fragmentLengthB = parseSeq(seq, cutterB)
    fragmentLenListB.append(fragmentLengthB)
    fragmentSeqListB.append(fragmentSeqB)


referneceSeq = parseFastaSeqs[0]
querySeq = parseFastaSeqs[0][19000000:21190000]

queryFragmentNumA, queryFragmentSeqA, queryFragmentLengthA = parseSeq(querySeq, cutterA)
queryFragmentNumB, queryFragmentSeqB, queryFragmentLengthB = parseSeq(querySeq, cutterB)
refFragmentNum, refFragmentSeq, refFragmentLength = parseSeq(referneceSeq, cutter)

N = 1
querylenprefixMatchA = queryFragmentLengthA[1:1+N]
print("querylenprefixMatchA", querylenprefixMatchA)
prefixMatchA = []
for i in range(len(refFragmentLength)-N-1):
    if refFragmentLength[i:i+N] == querylenprefixMatchA:
            prefixMatchA.append(sum(refFragmentLength[:i]))

N = 1
querylenListB = queryFragmentNumB[1:1+N]
print("querylenprefixMatchA", querylenListB)
prefixMatchB = []
for i in range(len(refFragmentLength)-N-1):
    if refFragmentLength[i:i+N] == querylenListB:
            prefixMatchB.append(sum(refFragmentLength[:i]))
prefixMatchB


N = 2
querylenList = queryFragmentLength[1:4]
prefixMatch = []
for i in range(len(refFragmentLength)-N-1):
    # if refFragmentLength[i] == querylenList[0] and refFragmentLength[i+1] == querylenList[1]:
    if refFragmentLength[i] == querylenList[0]:
        prefixMatch.append(refFragmentLength[i:i+3])
len(prefixMatch)


# Repeat Position
CutterA_positionList = []

for i in range(mostCommon):
    trv = commonRepeatInfo[i][0]
    targetDicA = positionDicA[trv]
    for j in range(len(targetDicA)):
        chrIndex = targetDicA[j][0]
        fragmentIndex = targetDicA[j][1]
        positionValue = sum(fragmentLenListA[chrIndex][:fragmentIndex])+(cutterALen*fragmentIndex)
        CutterA_positionList.append(tuple([chrIndex, positionValue]))
    
# Top 10 position
CutterA_positionList.sort(key=lambda tup: tup[0])

dA = dict()
[dA [t[0]].append(t [1]) if t [0] in list(dA.keys()) else dA.update({t [0]: [t [1]]}) for t in CutterA_positionList]
chrPositionKeys = list(dA.keys())

for i in range(len(chrPositionKeys)):
    dA[chrPositionKeys[i]].sort()
    print(f"\nChr10 : {','.join( str(v) for v in dA[chrPositionKeys[i]] )}\n")
    # print(f"{chrPositionKeys[i]}: {','.join( str(v) for v in dA[chrPositionKeys[i]] )}\n")




dB = dict()
[dB [t[0]].append(t [1]) if t [0] in list(dB.keys()) else dB.update({t [0]: [t [1]]}) for t in CutterB_positionList]
keysB = list(dB.keys())

for i in range(len(keysB)):
    dB[keysB[i]].sort()
    print(f"{keysB[i]}: {','.join( str(v) for v in dB[keysB[i]] )}\n")
    

# Analyze repeat sequence for the top 1 repeat count
topRepeatValue_ = commonRepeatInfo_[0][0]
topRepeatCount_ = commonRepeatInfo_[0][1]
topRepeatPosition_ = positionDic_[topRepeatValue_]
topOriginalRepeatSeqs_ = []
topCombinedRepeatSeqs_ = []
for i in range(len(topRepeatPosition_)):
    chrIndex, fragmentIndex = topRepeatPosition_[i][0], topRepeatPosition_[i][1]
    originalSeq = fragmentsSeqList_[chrIndex][fragmentIndex: fragmentIndex+3]
    combinedSeq = fragmentsSeqList_[chrIndex][fragmentIndex].join([fragmentsSeqList_[chrIndex][fragmentIndex+1], fragmentsSeqList_[chrIndex][fragmentIndex+2]])
    topOriginalRepeatSeqs_.append(originalSeq)
    topCombinedRepeatSeqs_.append(combinedSeq)


flattenTopRepeatSeqs_ = [seq for subseqs in topOriginalRepeatSeqs_ for seq in subseqs]
oriRepeatSeqCount_ = Counter(flattenTopRepeatSeqs_)
comRepeatSeqCount_ = Counter(topCombinedRepeatSeqs_)


k = list(dA[0])
l = list(dB[0])
kIndex, lIndex = 0, 0
tolerance = 400000
resList = []

while(kIndex < len(k) and lIndex < len(l)):
    if(abs(k[kIndex] - l[lIndex]) < tolerance):
        resList.append(tuple([k[kIndex], l[lIndex]]))
    if(kIndex+1 < len(k) and lIndex+1 < len(l)):
        if(abs(k[kIndex+1]-l[lIndex]) < abs(k[kIndex]-l[lIndex+1])):
            kIndex += 1
        else:
            lIndex += 1
    else:
        if(kIndex+1 >= len(k) and lIndex+1 < len(l)):
            lIndex += 1
        elif(lIndex+1 >= len(l) and kIndex+1 < len(k)):
            kIndex += 1
        else:
            break

resList