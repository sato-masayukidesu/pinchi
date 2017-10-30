
# coding: utf-8

# In[1]:

import os
from classes import MCS_Finder
import lxml.html
import requests
import datetime


# In[2]:

ari = []
for genus in os.listdir("test"):
    if genus == ".DS_Store" or genus == "others":
        continue
    elif os.listdir("test/" + genus) == []:
        continue
    else:
        ari.append(genus)


# In[3]:

print(len(ari))


# In[4]:

def gCfl(genus, label):
    import re
    label_list = re.split("[-()]", label)
    sep_list = re.split("[a-zA-Z][0-9]?[a-z]?", label)
    query = ""
    for i in range(label_list.count("")):
        label_list.remove("")
    for i in range(len(label_list)):
        label_list[i] += "[0-9]?[a-z]?"
    for l1, l2 in zip(sep_list, label_list):
        query += l1 + l2
    query = re.sub("\(", "[(]", query)
    query = re.sub("\)", "[)]", query)
    query += "\s"
    with open("test/" + genus + "/kcfs.kcfs")as f:
        file = f.read()
        molecule = file.split("///\n")
        Cnlist = []
        for mole in molecule:
            if re.search(query, mole) is not None:
                Cn = mole.split("\n")[0].split()[1]
                Cnlist.append(Cn)
        return Cnlist


# In[6]:

kosuu2 = dict()
for genus in ari:
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            kosuu2[(temp[1], temp[2])] = kosuu2.get((temp[1], temp[2]), 0) + int(temp[3])


# In[20]:

only = dict()
for genus in ari:
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            if kosuu2[(temp[1], temp[2])] == int(temp[3]):
                only[(temp[1], temp[2])] = int(temp[3])


# In[21]:

print(len(only))


# In[22]:

sorted(only.items(), reverse=True, key=lambda x: x[1])


# In[23]:

only2 = dict()
for genus in ari:
    if genus == "Streptomyces":
        pass
    else:
        continue
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            if kosuu2[(temp[1], temp[2])] == int(temp[3]):
                only2[(temp[1], temp[2])] = int(temp[3])


# In[24]:

print(len(only2))


# In[25]:

sorted(only2.items(), reverse=True, key=lambda x: x[1])


# In[26]:

print(len(kosuu2))


# In[29]:

total = 0
for i in kosuu2.items():
    total += i[1]
print(total)


# In[ ]:



