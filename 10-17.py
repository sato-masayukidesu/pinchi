
# coding: utf-8

# In[1]:

import sys
sys.path.append("../../module")
from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco
from rdkit import Chem


# In[2]:

f = MCS_Finder("Streptomyces")
Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[3]:

complist = []
for i in Cnlist:
    complist.append(Compound())
print(complist)


# In[4]:

pagerange = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f:
        clist = f.read().split()
        pagerange.append(clist[1])
print(pagerange)


# In[5]:

counter = 0
for z, i in enumerate(Cnlist):
    num = int(i[1:])
    for p, k in enumerate(pagerange[1:]):
        k2 = int(k[1:])
        if k2 > num:
            k3 = str(p+1)
            with open("../../../database/knapsack-kcf/KNApSAck" + k3 + ".kcf")as f:
                Clist = f.read().split("///\n")
                try:
                    for C in Clist:
                        if i == C.split()[1]:
                            molblock = kcfco.kcf_to_molblock(C)
                            # print("OK", i)
                            # print(molblock[1])
                            mol = Chem.MolFromMolBlock(molblock[1])
                            if mol is None:
                                print("None", i, z, k3)
                                if "#+" in C or "#-" in C:
                                    print("Charge in\n")
                                counter += 1
                                break
                            complist[z].input_rdkmol(mol)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# In[6]:

import time
start = time.time()
veclist = [0] * len(Cnlist)
for i in range(len(Cnlist)):
    veclist[i] = kcfco.kcf_vec(complist[i])
    print(time.time() - start)
print("finish")


# In[7]:

temp1 = list(veclist[0].pandas()["str"])
strlist = [] * len(Cnlist)
for vec in veclist:
    strlist.append(vec.pandas()["str"])


# In[8]:

kyoutuu = []
for vecstr in strlist[0]:
    for vec in strlist:
        if vecstr not in vec.values:
            break
    else:
        kyoutuu.append(vecstr)


# In[9]:

print(len(kyoutuu))


# In[10]:

import re
temp2 = []
for i in kyoutuu:
    temp3 = i.split("-")
    length = len(temp3)
    for i in re.findall("-[0-9]", "-".join(temp3)):
        length -= 1
    temp2.append((length, temp3))


# In[11]:

for group in temp2:
    if group[0] > 6:
        print("-".join(group[1]), group[0])


# In[12]:

kyoutuu2 = []
for vecstr in strlist[0]:
    no_counter = 0
    for vec in strlist:
        if vecstr not in vec.values:
            no_counter += 1
    else:
        if no_counter < len(strlist)/2:
            kyoutuu2.append(vecstr)


# In[13]:

print(len(kyoutuu2))


# In[14]:

import re
temp2 = []
for i in kyoutuu2:
    temp3 = i.split("-")
    length = len(temp3)
    for i in re.findall("-[0-9]", "-".join(temp3)):
        length -= 1
    temp2.append((length, temp3))
    
for group in temp2:
    if group[0] > 6:
        print("-".join(group[1]), group[0])


# In[16]:

set(veclist[0].pandas()["type"])


# 数に幅をもたせてkcf-sを抜き出してみた。  
# どこかで使うかも。

# In[ ]:



