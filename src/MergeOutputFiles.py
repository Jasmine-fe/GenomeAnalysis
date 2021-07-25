# %%
outputFragment = ''
outputIR = ''

with open('./outputFragment.txt', 'r') as fp:
    outputFragment = fp.read()
  
with open('./outputIR.txt', 'r') as fp:
    outputIR = fp.read()
  
mergedOutput = outputFragment + "\n" + outputIR

with open ('./mergedOutput.txt', 'w') as fp:
    fp.write(mergedOutput)
# %%
