#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)


# In[5]:


from DataStructure import SeqRepeatInfo


# In[1]:


# dataset
weevilDataset = 'weevil.fa'
mouseDataset = 'mm39.fa'
danioRerioDataset = 'DanioRerio.fasta'
humanDataset = 'hg38.fa'
human_chr1 = 'chr1.fa'

# match pattern
humanPattern = "^chr(?:\d*|[A-Z]*)$"
weevilPattern = "^weevil.scaffold\d+"
noPattern = "*"

# shared var
currDataset = human_chr1
datasetPath = './dataset/' + currDataset
matchPattern = noPattern
fragmentN = 5
cutter = 'GATC'
cutterLen = len(cutter)
commonCount = 5


# In[ ]:




