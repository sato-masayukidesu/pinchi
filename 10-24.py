
# coding: utf-8

# In[1]:

import os
from classes import MCS_Finder
import lxml.html
import requests
import datetime


# In[2]:

html = requests.get("https://ja.wikipedia.org/wiki/放線菌")


# In[3]:

print(html)


# In[4]:

dom = lxml.html.fromstring(html.text)
acti = []
for i1 in range(9, 15):
    moku = dom.xpath('//*[@id="mw-content-text"]/div/ul[' + str(i1) + ']/li/ul/li')
    # print(list(moku[0]))
    for i2 in moku:
        for i3 in list(i2)[1:]:
            if i3.tag == "a":
                acti.append(list(i3)[0].text)
            elif i3.tag == "i":
                acti.append(i3.text)
            else:
                print("DAME")


# In[5]:

print(len(acti))


# 前回MCS_Finderでたくさんの属を抜いたので今回はそれの整理を行う。

# In[6]:

print(os.getcwd())


# In[7]:

print(os.listdir("test"))


# In[8]:

print(len(os.listdir("test")))


# 一個多いのは.DS_Storeのせい 

# In[9]:

ari = []
for genus in os.listdir("test"):
    if genus == ".DS_Store" or genus == "others":
        continue
    elif os.listdir("test/" + genus) == []:
        continue
    else:
        ari.append(genus)


# In[10]:

print(len(ari))


# In[11]:

print(ari)


# In[12]:

def get_all_cnumber_from_kcfs(genus):
    with open("test/" + genus + "/kcfs.kcfs")as f:
        Cnlist = []
        molecule = f.read().split("///\n")
        for mol in molecule[:-1]:
            Cn = mol.split("\n")[0].split()[1]
            Cnlist.append(Cn)
    return Cnlist


# In[13]:

get_all_cnumber_from_kcfs(ari[0])


# In[14]:

arilist = []
for genus in ari:
    arilist.append(get_all_cnumber_from_kcfs(genus))


# In[15]:

print(arilist)


# In[16]:

kosuu = dict()
for Cnlist in arilist:
    for Cn in Cnlist:
        kosuu[Cn] = kosuu.get(Cn, 0) + 1


# In[17]:

print(kosuu)


# In[18]:

print(len(kosuu))


# In[19]:

counter = 0
for i in kosuu.items():
    if i[1] > 1:
        print(i)
        counter += i[1] -1
print(len(kosuu)-counter)


# ダブってる奴は少なかった。

# 逆に共通構造や特異的な構造は見つけやすいのでは？

# In[14]:

kosuu2 = dict()
for genus in ari:
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            kosuu2[(temp[1], temp[2])] = kosuu2.get((temp[1], temp[2]), 0) + int(temp[3])


# In[15]:

print(kosuu2)


# In[16]:

sorted(kosuu2.items(), reverse=True, key=lambda x: x[1])


# 44個もある大きめのスケルトンを抜く

# In[17]:

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


# In[18]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C(C-C-C-C(C)-C-C-C(C-C)-C(C)-C)-C-C-C-C-C(C)-C-C-C-C-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[19]:

dekai = []
for genus in ari:
    Cnlist = gCfl(genus, 'C-C(C-C-C-C(C)-C-C-C(C-C)-C(C)-C)-C-C-C-C-C(C)-C-C-C-C-C')
    if Cnlist != []:
        for Cn in Cnlist:
            dekai.append(Cn)


# In[20]:

print(len(dekai), len(set(dekai)))


# 44個あったあのでかい奴は31化合物27種類に入れられていた。  
# こいつをグリッドで書いてみる

# In[21]:

Cnlist = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f1:
        clist = f1.read().split()
        Cnlist.append(clist[1])
print(Cnlist)


# In[25]:

import sys
sys.path.append("../../module")
from rdkit.Chem import rdDepictor
import kcf.converter as kcfco
from rdkit import Chem
mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(sorted(dekai)):
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
                            # rdDepictor.Compute2DCoords(mol)
                            mol_list.append(mol)
                            nCnumber.append(i)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# In[26]:

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(mol_list, legends=sorted(dekai), subImgSize=(400, 400))
img.save("test/others/test_noD.png")


# from rdkit.Chem import rdFMCS
# mcs = rdFMCS.FindMCS(mol_list)
# mcs_smarts = mcs.smartsString
# mcsMol = Chem.MolFromSmarts(mcs_smarts, mergeHs=True)
# Draw.MolToFile(Chem.Mol(mcsMol.ToBinary()), "test/others/mcstest.png", kekulize=False)

# mcsを抜こうとするとカーネルデッドする  
# 理由は不明

# 抜いたやつを眺めてるとイソプレン単位が多そうに見える。  
# カロテノイドかな？

# 共通じゃなくて単独の構造も抜きたい。

# In[27]:

for genus in ari:
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            if kosuu2[(temp[1], temp[2])] == int(temp[3]) > 10:
                print(temp)


# In[28]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C(C-C-C-C(C))-C-C-C-C-C-C-C-C(C-C-C-C)-C-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[29]:

dekai2 = []
for genus in ari:
    Cnlist = gCfl(genus, 'C-C(C-C-C-C(C))-C-C-C-C-C-C-C-C(C-C-C-C)-C-C')
    if Cnlist != []:
        for Cn in Cnlist:
            dekai2.append(Cn)


# In[30]:

print(len(dekai2))


# In[31]:

hanni = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f1:
        clist = f1.read().split()
        hanni.append(clist[1])
print(hanni)


# In[32]:

import sys
sys.path.append("../../module")
from rdkit.Chem import rdDepictor
import kcf.converter as kcfco
from rdkit import Chem
mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(sorted(dekai2)):
    num = int(i[1:])
    for p, k in enumerate(hanni[1:]):
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
                            # rdDepictor.Compute2DCoords(mol)
                            mol_list.append(mol)
                            nCnumber.append(i)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# In[37]:

print(mol_list)


# In[33]:

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(mol_list, legends=sorted(dekai2), subImgSize=(400, 400))
img.save("test/others/testonly_noD.png")


# C00016805のkcfsを見たが、スケルトンのスタート基準がよくわからない  
# また同じ構造に対して7つも表現方法がある。

# 水酸基のつきかたが微妙に違ったりする。

# 今までやっていた構造もdekaiから探してみる。

# In[34]:

for genus in ari:
    Cnlist = gCfl(genus, "C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# 数日前の結果からもわかるように固有ではないが圧倒的に大きい。

# スケルトンについてソースを読んだ方がいいかもしれない

# In[35]:

dekai3 = []
for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C(C-C-C-C-C-C(C)-C)-C-C-C-C-C-C-C(C)-C-C')
    if Cnlist != []:
        for Cn in Cnlist:
            dekai3.append(Cn)


# In[36]:

print(len(dekai3))


# In[37]:

print(dekai3 == dekai2)


# In[38]:

dekai4 = []
for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-C-C-C(C-C)-C(C-C-C)-C-C-C')
    if Cnlist != []:
        for Cn in Cnlist:
            dekai4.append(Cn)


# In[39]:

print(len(dekai4))


# In[40]:

print(dekai3 == dekai4)


# In[41]:

only = dict()
for genus in ari:
    with open("test/" + genus + "/splitedcount.txt")as f:
        units = f.read().split("\n")
        for unit in units[:-1]:
            temp = unit.split()
            if kosuu2[(temp[1], temp[2])] == int(temp[3]):
                only[(temp[1], temp[2])] = int(temp[3])


# In[42]:

sorted(only.items(), reverse=True, key=lambda x: x[1])


# In[43]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-N-C-C-N-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[44]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-N-C-C-N-C')
    if Cnlist != []:
        onlyst = Cnlist


# In[45]:

print(len(onlyst))


# In[46]:

def check(label, aite):
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
    if re.search(query, aite) is not None:
        return True
    else:
        return False


# In[48]:

import re
counter = 0
coun = 0
label = 'C-C-C-C-C-N-C-C-N-C'
with open("test/Streptomyces/kcfs.kcfs", "r")as f:
    moles = f.read().split("///\n")
    for Cn in onlyst:
        for mole in moles[:-1]:
            if Cn == mole.split()[1]:
                coun += 1
                print(Cn)
                for line in mole.split("\n"):
                    i = 0
                    if len(line.split()) < 2:
                        continue
                    elif line.split()[i] == "SKELETON":
                        # print("SKE")
                        i += 1
                    elif line.split()[i] == "RING":
                        # print("RIN")
                        i += 1
                    if check(label, line.split()[i]):
                        print(line.split()[i+1][1:-1])
                        counter += int(line.split()[i+1][1:-1])


# In[49]:

print(counter)


# In[50]:

print(coun)


# レベル0で探しているのでモノによっては2個以上出てくるCnがある。

# In[51]:

import sys
sys.path.append("../../module")
from rdkit.Chem import rdDepictor
import kcf.converter as kcfco
from rdkit import Chem
mol_list = []
counter = 0
nCnumber = []
for z, i in enumerate(sorted(onlyst)):
    num = int(i[1:])
    for p, k in enumerate(hanni[1:]):
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
                            # rdDepictor.Compute2DCoords(mol)
                            mol_list.append(mol)
                            nCnumber.append(i)
                            if "#+" in C or "#-" in C:
                                print(i, z, k3, "Charge in\n")
                            break
                except IndexError:
                    counter += 1
                    print("DAME", i, z)
            break
print(counter)


# In[52]:

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(mol_list, legends=sorted(onlyst), subImgSize=(400, 400))
img.save("test/others/testonlyst_noD.png")


# 放線菌内でのStreptomyces固有の構造が抜けた。  
# C00016546の表示がおかしいが、正しい化合物は表示されている。  
# 化合物名などを回収できるようにしたい。

# In[ ]:



