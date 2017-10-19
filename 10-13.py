
# coding: utf-8

# In[1]:

import sys


# In[2]:

sys.path.append("../../module")


# In[3]:

from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco


# In[4]:

f = MCS_Finder("Streptomyces")
Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[5]:

complist = [] * len(Cnlist)
for i in Cnlist:
    complist.append(Compound())
print(complist)


# In[6]:

for i, Cn in enumerate(Cnlist):
    complist[i].input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    complist[i].fit2d = False
    print(Cn, i)


# In[7]:

veclist = [0] * len(Cnlist)


# In[8]:

import time
start = time.time()
for i in range(len(Cnlist)):
    veclist[i] = kcfco.kcf_vec(complist[i])
    print(time.time() - start)
print("finish")


# In[9]:

temp1 = list(veclist[0].pandas()["str"])


# In[10]:

print(temp1)


# In[12]:

strlist = [] * len(Cnlist)
for vec in veclist:
    strlist.append(vec.pandas()["str"])


# In[13]:

print(strlist)


# In[27]:

for vecstr in strlist[0]:
    for vec in strlist:
        if vecstr not in vec.values:
            break
    else:
        print(vecstr)


# In[28]:

kyoutuu = []


# In[29]:

for vecstr in strlist[0]:
    for vec in strlist:
        if vecstr not in vec.values:
            break
    else:
        kyoutuu.append(vecstr)


# In[30]:

print(len(kyoutuu))


# In[31]:

print(kyoutuu[0])


# In[33]:

import re
temp2 = []
for i in kyoutuu:
    temp3 = i.split("-")
    length = len(temp3)
    for i in re.findall("-[0-9]", "-".join(temp3)):
        length -= 1
    temp2.append((length, temp3))


# In[34]:

print(temp2)


# In[39]:

for group in temp2:
    if group[0] > 6:
        print("-".join(group[1]), group[0])


# 共通なkcfを抜き出した  
# 新しい方のkcfはunitのおかげで小さい単位で抜きやすい。  
# いっそのこと全化合物でやってみる?

# In[ ]:



