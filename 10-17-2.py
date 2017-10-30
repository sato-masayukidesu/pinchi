
# coding: utf-8

# In[1]:

import sys
sys.path.append("../../module")
from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco
from rdkit.Chem import rdFMCS
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor


# In[2]:

f = MCS_Finder("Streptomyces")


# In[3]:

html = f.get_html(f.genus)


# In[4]:

print(html)


# In[5]:

get_ipython().run_cell_magic('time', '', 'Cnumber = f.get_Cnumber(html)')


# In[6]:

print(Cnumber[:-4])


# In[7]:

print(len(Cnumber))


# In[8]:

Cnlist = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f1:
        clist = f1.read().split()
        Cnlist.append(clist[1])
print(Cnlist)


# In[ ]:




# In[9]:

counter = 0
f.mol_list = []
nCnumber = []
for z, i in enumerate(Cnumber):
    num = int(i[1:])
    for p, k in enumerate(Cnlist[1:]):
        k2 = int(k[1:])
        if k2 > num:
            k3 = str(p+1)
            with open("../../../database/knapsack-kcf/KNApSAck" + k3 + ".kcf")as f2:
                Clist = f2.read().split("///\n")
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
                            rdDepictor.Compute2DCoords(mol)
                            f.mol_list.append(mol)
                            nCnumber.append(i)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# 

# In[10]:

img = Draw.MolsToGridImage(f.mol_list, legends=Cnumber[:-4], subImgSize=(400, 400))
img.save(f.genus + "/allallcomp.png")


# In[11]:

print(len(f.mol_list))


# allallcomp.pngはファイルサイズがでかすぎて????読めなかった。

# In[ ]:




# In[12]:

Cnumber = f.get_Cnlist_from_label2("C(C)(C-C(C)-C-C-C-C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C-C-C")


# In[13]:

print(Cnumber)


# In[14]:

f.mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(Cnumber):
    num = int(i[1:])
    for p, k in enumerate(Cnlist[1:]):
        k2 = int(k[1:])
        if k2 > num:
            k3 = str(p+1)
            with open("../../../database/knapsack-kcf/KNApSAck" + k3 + ".kcf")as f2:
                Clist = f2.read().split("///\n")
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
                            rdDepictor.Compute2DCoords(mol)
                            f.mol_list.append(mol)
                            nCnumber.append(i)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# In[15]:

mcs = rdFMCS.FindMCS(f.mol_list, matchValences=True, completeRingsOnly=True)
mcs_smarts = mcs.smartsString
mcs_mol = Chem.MolFromSmarts(mcs_smarts)
match_list = []
for m in f.mol_list:
    match_atoms = m.GetSubstructMatch(mcs_mol)
    match_list.append(match_atoms)
img = Draw.MolsToGridImage(f.mol_list, highlightAtomLists=match_list, legends=Cnumber, subImgSize=(400, 400))
img.save(f.genus + "/test1.png")


# また面白そうな形が釣れてしまった  
# 一個だけ二重結合の数が少なく、酸素の入っている環になっているものがあった(C00015659)

# In[16]:

dekalist = []
with open("Streptomyces/splitedcount.txt", "r")as fi:
    linelist = fi.read().split("\n")
    for line in linelist:
        if line.count("-") + line.count("(") + 1 > 6 and line.split("\t")[1] != "RING":
            dekalist.append(line)


# In[17]:

print(dekalist)


# In[18]:

print(len(dekalist))


# In[19]:

search = []
for deka in dekalist:
    if int(deka.split("\t")[3]) > 10:
        search.append(deka)


# In[20]:

print(len(search))


# In[21]:

for i in search:
    print(i)


# SKELETONの大きい奴だけ抜き出して数が少ないのも切った。  
# これに関して全部MCSを探す？

# In[22]:

d = {"0":0, "1":0, "2":0}
for i in search:
    level = len(i.split("\t")[2].split("-")[0].split("(")[0]) - 1
    d[str(level)] += 1
print(d)


# In[23]:

import re
label = "C8y(C1x)-C8y(C5x)-C8y(C8y-C8x-C8x-C8x-C8x-C8y)-C8y-C8y-C8y-C8y-C8x-C8x-C8x-C8x-C8y"
label_list = re.split("[-()]", label)
sep_list = re.split("[a-zA-Z][0-9]?[a-z]?", label)
query = ""
for i in range(label_list.count("")):
    label_list.remove("")
for i in range(len(label_list)):
    label_list[i] = label_list[i][0] + "[0-9]?[a-z]?"
for l1, l2 in zip(sep_list, label_list):
    query += l1 + l2
query = re.sub("\(", "[(]", query)
query = re.sub("\)", "[)]", query)


# In[24]:

print(query)


# In[25]:

for i in search:
    if re.search(query, i.split("\t")[2]) is not None:
        print(i)


# In[26]:

print(label_list)


# In[27]:

print(type("str"))


# 被りのやつは探さなくてもいいんじゃない？  
# なるべく大きい奴を引っ張ってきたい  
# 特徴的なtripletには何かないかな？

# In[ ]:



