
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


# In[27]:

triplet = []
with open("Streptomyces/kcfscount.txt", "r")as f:
    lines = f.read().split("\n")
    for line in lines:
        if line == "":
            continue
        if line.split("\t")[1] == "TRIPLET":
            triplet.append(line)


# In[4]:

print(triplet)


# In[5]:

print(len(triplet))


# In[7]:

d = {"1":0, "2":0, "3":0}
for triple in triplet:
    d[str(len(triple.split("\t")[2].split("-")[0]))] += 1
print(d)


# In[8]:

for triple in triplet:
    if len(triple.split("\t")[2].split("-")[0]) == 1:
        print(triple)


# In[33]:

f = MCS_Finder("Streptomyces")


# In[14]:

Cnumber = f.get_Cnlist_from_label2("C-S-S")


# In[15]:

print(Cnumber)


# In[16]:

Cnlist = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f1:
        clist = f1.read().split()
        Cnlist.append(clist[1])
print(Cnlist)


# In[17]:

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


# In[18]:

mcs = rdFMCS.FindMCS(f.mol_list, matchValences=True, completeRingsOnly=True)
mcs_smarts = mcs.smartsString
mcs_mol = Chem.MolFromSmarts(mcs_smarts)
match_list = []
for m in f.mol_list:
    match_atoms = m.GetSubstructMatch(mcs_mol)
    match_list.append(match_atoms)
img = Draw.MolsToGridImage(f.mol_list, highlightAtomLists=match_list, legends=Cnumber, subImgSize=(400, 400))
img.save(f.genus + "/testtri.png")


# 化合物は9個しか拾えなかったけどこのトリプレットは21個あった。  
# トリプレットに注目するのはアリかもしれないけどダブルカウント、トリプレカウントがあることに注意したい。

# 5/9がタンパク質のジスルフィド結合だった  
# 残りはよくわかんないの

# In[19]:

Cnumber2 = f.get_Cnlist_from_label2("N-N-N")


# In[20]:

print(Cnumber2)


# In[22]:

f.mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(Cnumber2):
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


# In[23]:

mcs = rdFMCS.FindMCS(f.mol_list, matchValences=True, completeRingsOnly=True)
mcs_smarts = mcs.smartsString
mcs_mol = Chem.MolFromSmarts(mcs_smarts)
match_list = []
for m in f.mol_list:
    match_atoms = m.GetSubstructMatch(mcs_mol)
    match_list.append(match_atoms)
img = Draw.MolsToGridImage(f.mol_list, highlightAtomLists=match_list, legends=Cnumber, subImgSize=(400, 400))
img.save(f.genus + "/testtri2.png")


# In[24]:

print(f.mol_list)


# In[25]:

img = Draw.MolsToGridImage(f.mol_list, legends=Cnumber, subImgSize=(400, 400))
img.save(f.genus + "/testtri2.png")


# mcsが存在しないため怒られた。  
# たいして面白い分子じゃなかった。  

# In[28]:

bond = []
with open("Streptomyces/kcfscount.txt", "r")as fi:
    lines = fi.read().split("\n")
    for line in lines:
        if line == "":
            continue
        if line.split("\t")[1] == "BOND":
            bond.append(line)


# In[29]:

print(bond)


# In[30]:

print(len(bond))


# In[31]:

for bonbon in bond:
    if len(bonbon.split("\t")[2].split("-")[0]) == 1:
        print(bonbon)


# In[34]:

Cnumber3 = f.get_Cnlist_from_label2("O-Z-O")
print(Cnumber3)


# In[35]:

f.mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(Cnumber3):
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


# 5/6が金属錯体だった。

# In[44]:

def get_label_test(label):
    import re
    Cnlist = []
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
    for page in ["1-9", "10-19", "20-29", "30-39", "40-49", "50-59"]:
        with open("../../../database/kcfs/KNApSAck" + page + ".kcfs")as f:
            file = f.read()
            molecule = file.split("///\n")
            for mole in molecule:
                if re.search(query, mole) is not None:
                    Cn = mole.split("\n")[0].split()[1]
                    Cnlist.append(Cn)
    return Cnlist


# In[45]:

test = get_label_test("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[46]:

print(test)


# In[47]:

print(len(test))


# In[48]:

lis = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[49]:

print(lis)


# In[50]:

sabun = []
for i in test:
    if i not in lis:
        print(i)
        sabun.append(i)


# In[51]:

print(sabun)


# In[52]:

print(len(sabun))


# In[53]:

f.mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(sabun):
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


# In[55]:

img = Draw.MolsToGridImage(f.mol_list, legends=sabun, subImgSize=(400, 400))
img.save(f.genus + "/testsabun.png")


# In[62]:

import time
genuses = []
for num in sabun:
    genuses.append(f.get_genuses(num))
    time.sleep(0.5)


# In[63]:

print(genuses)


# In[64]:

len(genuses)


# In[61]:

import os
print(os.getpid())


# In[65]:

for genus in genuses:
    print(genus)


# In[66]:

for i, genus in enumerate(genuses):
    if genus == []:
        print(sabun[i])


# 手持ちの全部のkcf-sから抜き出してみた。  
# 2個ほどページが存在しなかった。

# この単位はストレプトマイセス固有だったのか怪しいが、一つの属で20こ持っているのはこれなので固有感もある。

# In[ ]:



