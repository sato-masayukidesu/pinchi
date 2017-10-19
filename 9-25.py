
# coding: utf-8

# 9-22の続き  
# ベクトルを作るところから

# In[1]:

import sys


# In[2]:

sys.path.append("../../module")


# In[3]:

from nxrd.Compound import Compound


# In[4]:

from classes import MCS_Finder


# In[5]:

f = MCS_Finder("Streptomyces")
Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")
print(Cnlist)


# In[6]:

complist = []
for i in Cnlist:
    complist.append(Compound())
print(complist)
print(len(complist))


# In[7]:

for i, Cn in enumerate(Cnlist):
    complist[i].input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    complist[i].fit2d = False
    print(Cn, i)


# In[8]:

import kcf.converter as kcfco


# In[9]:

veclist = [0] * len(Cnlist)


# In[10]:

import time


# In[11]:

start = time.time()
for i in range(len(Cnlist)):
    if i == 11:
        print(i)
        veclist[i] = kcfco.kcf_vec(complist[i])
        break
print(time.time() - start)


# 割と毎回pinpass3までは一瞬でそこで時間がかかるので形が悪そう  
# pinpass3の中のall_simple_path  
# 環構造の少ないやつで試したい  
# 試すなら別のやつで

# In[ ]:



