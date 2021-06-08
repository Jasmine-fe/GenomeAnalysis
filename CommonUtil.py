# In[ ]:
import csv
import matplotlib.pyplot as plt
# In[ ]:
# sortedLenCounter = [(lengthValue: count), ... ]
def fragmentLenPlot(sortedLenList):
    x, y = [], []
    for len, count in sortedLenList:
        x.append(len)
        y.append(count)
    plt.figure(figsize=(20,12)) 
    plt.title("Distribution - fragment Length", fontsize=30)
    plt.xlabel("length", fontsize=15)
    plt.ylabel("count", fontsize=15)
    plt.plot(x, y)
    return 0

# In[1]:
def exportCsvFile(fileName, rowList):
    with open(fileName, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['length','uniqueCount'])
        for row in rowList:
            writer.writerows(rowList)