
# coding: utf-8

# In[1]:

import os
from classes import MCS_Finder
import datetime
import re
import sys
sys.path.append("../../module")
from rdkit.Chem import rdDepictor
import kcf.converter as kcfco
from rdkit import Chem

from nxrd.Compound import Compound


# In[2]:

c = Compound()


# In[3]:

hanni = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f1:
        clist = f1.read().split()
        hanni.append(clist[1])
print(hanni)


# In[4]:

i = "C00016546"
num = 16546
mol_list = []
nCnumber = []
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


# In[5]:

print(mol_list)


# In[6]:

c.input_rdkmol(mol_list[0])


# In[7]:

get_ipython().magic('matplotlib inline')
c.draw_cpd_with_labels()


# C00016564の表示がおかしいのはrdkitの二次元座標計算のせい？  
# 元から入っている二次元座標を使えば綺麗にかけた  
# 一方で元から入っているやつを使うとシクロヘキサンなどの構造が立体になっている。test/test_noD.png  
# ベンゼン環とかも模式的になって見やすくなったが、一部見にくいものができてしまった。

# C番号から化合物名を持ってきたい。

# In[8]:

import requests
import lxml.html


# In[ ]:

html = requests.get("http://kanaya.naist.jp/knapsack_jsp/information.jsp?word=C00016546")


# In[ ]:

print(html)


# In[ ]:

dom = lxml.html.fromstring(html.text)


# In[ ]:

dom.xpath('//*[@id="my_contents"]/table/tr[2]/td[1]/table/tr[1]/td')[0].text


# In[ ]:

def get_name(Cnumber):
    import requests
    import lxml.html
    html = requests.get("http://kanaya.naist.jp/knapsack_jsp/information.jsp?word=" + Cnumber)
    dom = lxml.html.fromstring(html.text)
    name = dom.xpath('//*[@id="my_contents"]/table/tr[2]/td[1]/table/tr[1]/td')[0].text
    return name


# In[46]:

print(get_name("C00016546"))


# In[10]:

ari = []
for genus in os.listdir("test"):
    if genus == ".DS_Store" or genus == "others":
        continue
    elif os.listdir("test/" + genus) == []:
        continue
    else:
        ari.append(genus)


# In[11]:

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
    query = "\s" + query
    with open("test/" + genus + "/kcfs.kcfs")as f:
        file = f.read()
        molecule = file.split("///\n")
        Cnlist = []
        for mole in molecule:
            if re.search(query, mole) is not None:
                Cn = mole.split("\n")[0].split()[1]
                Cnlist.append(Cn)
        return Cnlist


# In[12]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-N-C-C-N-C')
    if Cnlist != []:
        onlyst = Cnlist


# In[13]:

print(onlyst, len(onlyst))


# In[14]:

get_ipython().run_cell_magic('time', '', 'import time\nimport datetime\nCname = dict()\nprint(datetime.datetime.now())\nfor Cnumber in onlyst:\n    name = get_name(Cnumber)\n    Cname[Cnumber] = name\n    time.sleep(10)')


# In[57]:

print(Cname)


# In[58]:

len(Cname)


# In[61]:

sorted(Cname.items(), key=lambda x: x[0])


# 名前から抜けてそうなものが予測できる？

# In[66]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-N-C-C-N-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[72]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-C-C-N-C-C-N')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# 位置のずれたNを2個含む10員環

# In[79]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-N-C-C-N-C-C-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[80]:

for genus in ari:
    Cnlist = gCfl(genus, 'C-C-C-C-N-C-C-N-C-C')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# これらは同じカテゴリに含まれるべきではないのか  
# レベルを一段上げて絞ってみる？

# In[83]:

for genus in ari:
    Cnlist = gCfl(genus, 'C8-C8-C8-C8-C8-N5-C8-C8-N5-C8')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[84]:

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


# In[86]:

import re
counter = 0
coun = 0
label = 'C8-C8-C8-C8-C8-N5-C8-C8-N5-C8'
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



# In[87]:

for genus in ari:
    Cnlist = gCfl(genus, 'C1b-C1b')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()


# In[15]:

for genus in ari:
    Cnlist = gCfl(genus, 'C1-C1-C5-C1-C1-C1')
    if Cnlist != []:
        print(genus)
        print(Cnlist)
        print()
        onlyst = Cnlist


# In[17]:

print(len(onlyst))


# In[90]:

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
                            rdDepictor.Compute2DCoords(mol)
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


# In[91]:

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(mol_list, legends=sorted(onlyst), subImgSize=(400, 400))
img.save("test/others/testonlyst2.png")


# In[ ]:

from rdkit.Chem import rdFMCS
mcs = rdFMCS.FindMCS(mol_list)
mcs_smarts = mcs.smartsString
mcsMol = Chem.MolFromSmarts(mcs_smarts, mergeHs=True)
Draw.MolToFile(Chem.Mol(mcsMol.ToBinary()), "test/others/mcstest.png", kekulize=False)


# In[ ]:



