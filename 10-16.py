
# coding: utf-8

# In[31]:

import sys
sys.path.append("../../module")
from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco
from rdkit import Chem


# In[2]:

f = MCS_Finder("Streptomyces")


# In[3]:

html = f.get_html(f.genus)


# In[4]:

print(html)


# In[5]:

Cnumber = f.get_Cnumber(html)


# In[6]:

print(Cnumber)


# In[7]:

print(len(Cnumber))


# In[ ]:

with open("Streptomyces/kcfs.kcfs", "r")as f:
    text = f.read()
    kcf = text.split("///\n")
    molblock = kcfco.kcf_to_molblock(kcf[0])


# In[8]:

Cnlist = []
for i in range(1, 52):
    page = str(i)
    with open("../../../database/knapsack-kcf/KNApSAck" + page + ".kcf")as f:
        clist = f.read().split()
        Cnlist.append(clist[1])
print(Cnlist)


# In[9]:

len(Cnlist)


# In[10]:

int(Cnumber[0][1:])


# In[15]:

complist = []
for i in Cnumber:
    complist.append(Compound())
print(complist)


# In[16]:

print(len(complist))


# In[20]:

for z, i in enumerate(Cnumber):
    num = int(i[1:])
    for p, k in enumerate(Cnlist[1:]):
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
                            print(molblock[1])
                            complist[z].input_molblock(molblock[1])
                            break
                except IndexError:
                    print("DAME", i)
            break
    break


# In[21]:

get_ipython().magic('matplotlib inline')
complist[0].draw_cpd_with_labels()


# In[55]:

counter = 0
for z, i in enumerate(Cnumber):
    num = int(i[1:])
    for p, k in enumerate(Cnlist[1:]):
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


# In[52]:

print("sdf" in "lkjsdflkj")


# In[44]:

#print(Chem.MolFromMolBlock(molblock[1]))
print(molblock[1])
#print(C)


# In[ ]:

veclist = [0] * len(Cnumber)


# In[64]:

import time
import gc
start = time.time()
for i in range(len(Cnumber)):
    if veclist[i] == 0:
        veclist[i] = kcfco.kcf_vec(complist[i])
        print(i, time.time() - start)
        gc.collect()
        break
print("finish")


# In[38]:

print(veclist[0].pandas())


# In[28]:

print(complist[0].mol)


# In[57]:

import resource


# In[63]:

rsrc = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrc)
print(soft, hard)


# In[62]:

gc.collect()


# In[67]:

sys.getsizeof(complist[2])


# In[69]:

get_ipython().run_cell_magic('prun', '', 'start = time.time()\nfor i in range(len(Cnumber)):\n    if veclist[i] == 0:\n        veclist[i] = kcfco.kcf_vec(complist[i])\n        print(i, time.time() - start)\n        gc.collect()\n        break\nprint("finish")')


# 時間がかかりすぎるので全部やるのは難しそう？  
# 多分メモリが足りなくなる。

# ベクトル作って結果だけ全部辞書に突っ込んだらいけそう？

# 10-13の奴は全部に共通な奴だけを抜き出している。  
# 数が多めのやつを抜き出せばそれっぽい単位が見つかるのではないか？

# かかる時間の予想は4時間ほど

# 現状持っているか持っていないかでしか判断していないが2個以上持っていた場合にダブルカウントする必要があると思う。

# In[70]:

get_ipython().run_cell_magic('timeit', '', 'veclist[0] = kcfco.kcf_vec(complist[0])')


# 1loopってなんだよ

# In[71]:

sys.getsizeof(complist)


# In[72]:

sys.getsizeof(veclist)


# 参照の向こう側まではサイズを計算してくれないっぽい？

# In[ ]:



