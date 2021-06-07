#%%
import re
import requests
from bs4 import BeautifulSoup

# %%

# Fix modify soup parser to get data output as STRBaseRepeatRef.json
def getSTRLociReatSeq(lociId):
    url = f'https://strbase.nist.gov/str_{lociId}.htm'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    matchTextList = soup.find(lambda tag:tag.name=="p" and "Repeat" in tag.text)

    print(matchTextList)
    # paragraphResult =soup.select('p')
    # repeatParagraph = paragraphResult[3]
    # print("hihi", repeatParagraph)
    reapeatEleList = []
    for ele in matchTextList:
        matchObject = re.search(r'\[\w*\]', str(ele))
        result =  matchObject and matchObject.group() or None
        if result: reapeatEleList.append(re.search(r'[a-zA-Z]+', result).group())
    return reapeatEleList

# %% 
STRLocis = ['D1S1656', 'D2S441', 'D2S1338', 'D10S1248', 'D12S391', 'D19S433', 'D22S1045']
STRLocisSeqs = [ getSTRLociReatSeq(i) for i in STRLocis[:4] ]
# STRLocisSeqs
# %%
