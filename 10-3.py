
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

i = 0
for vec in veclist:
    if vec == 0:
        print("no file", i)
    else:
        print(vec, i)
    i += 1
print(len(Cnlist))


# In[10]:

temp1 = list(veclist[0].pandas()["str"])


# In[13]:

print(temp1)


# In[14]:

import re


# In[76]:

temp2 = []
for i in temp1:
    temp3 = i.split("-")
    length = len(temp3)
    for i in re.findall("-[0-9]", "-".join(temp3)):
        length -= 1
    temp2.append((length, temp3))


# In[77]:

print(re.findall("[0-9]-[0-9]", "-".join(temp3)))


# In[78]:

print(temp2)


# In[79]:

sotemp2 = sorted(temp2)
print(sotemp2)


# In[80]:

counter = 0
node_counter = dict()
for i in sotemp2:
    for k in range():
        if i[0] == k:
            if len(i[1][0]) == 2:
                counter += 1
                node_counter[str(k)] = node_counter.get(str(k), 0) + 1


# In[81]:

"1" in node_counter


# In[87]:

for i in node_counter.keys():
    if 2 == node_counter[i]:
        print(i)


# In[85]:

for i in temp2:
    if i[0] == 19:
        if len(i[1][0]) == 2:
            print("-".join(i[1]), len(i[1]))


# In[86]:

for i in temp2:
    if i[0] == 38:
        if len(i[1][0]) == 2:
            print("-".join(i[1]), len(i[1]))


# In[88]:

for i in temp2:
    if i[0] == 33:
        if len(i[1][0]) == 2:
            print("-".join(i[1]), len(i[1]))


# In[89]:

for i in temp2:
    if i[0] == 1:
        if len(i[1][0]) == 2:
            print("-".join(i[1]), len(i[1]))


# In[90]:

for i in temp2:
    if i[0] == 37:
        if len(i[1][0]) == 2:
            print("-".join(i[1]), len(i[1]))


# In[91]:

print(Cnlist)


# In[92]:

f.search(Cnlist, "Streptomyces/test.kcfs")


# In[94]:

f.kcfs2count("Streptomyces/test.kcfs", "Streptomyces/test.txt")


# In[96]:

f.split("Streptomyces/test.txt", "Streptomyces/splitedtest.kcfs")


# In[ ]:



