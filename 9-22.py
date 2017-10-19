
# coding: utf-8

# In[2]:

from simnet import *


# まずは先生の作ったkcfco.similalityを試す。

# In[3]:

from classes import MCS_Finder


# In[4]:

f = MCS_Finder("Streptomyces")


# In[5]:

Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[34]:

print(Cnlist)


# In[7]:

import sys


# In[8]:

sys.path.append("../../module")


# In[9]:

from nxrd.Compound import Compound


# In[9]:

c = Compound()
for Cn in Cnlist:
    with open('KNApSAck_mol/%s.mol' % (Cn))as file:
        molblock = file.read()
        c.input_molblock(molblock)
        c.fit2d = False
        print(molblock)
        break


# In[10]:

print(c.n_atoms, c.n_bonds)


# In[10]:

get_ipython().magic('matplotlib inline')


# In[12]:

c.draw_cpd_with_labels()


# In[11]:

from rdkit import Chem
from rdkit.Chem import rdDepictor


#  Compound.input.molblockを使って化合物を入れるとself.molにmolオブジェクトが入らない

# In[35]:

c = Compound()
for Cn in Cnlist:
    c.input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    # c.fit2d = False
    print(Cn)
    break


# In[36]:

print(c.n_atoms, c.n_bonds)


# In[37]:

c.draw_cpd_with_labels()


# In[18]:

complist = [] * len(Cnlist)
for i in Cnlist:
    complist.append(Compound())


# In[19]:

print(complist)


# In[20]:

for i, Cn in enumerate(Cnlist):
    complist[i].input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    complist[i].fit2d = False
    print(Cn, i)


# In[21]:

print(len(Cnlist))


# In[22]:

import kcf.converter as kcfco


# In[23]:

kcfco.similarity(complist[0], complist[1])


# In[24]:

veclist = [0] * len(Cnlist)


# In[25]:

import time


# In[30]:

start = time.time()
for i in range(len(Cnlist)):
    veclist[i] = kcfco.kcf_vec(complist[i])
print(time.time() - start)


# バカみたいに時間がかかるのでやるときは余裕を持ってやろう!!

# In[31]:

i = 0
for vec in veclist:
    if vec == 0:
        i += 1
        print("no file", i)
    else:
        print(vec)


# In[33]:

for vec1 in veclist:
    for vec2 in veclist:
        print(kcfco.similarity(vec1, vec2))
    print("\n")


# In[39]:

veclist2 = [0] * len(Cnlist)
start = time.time()
for i in range(len(Cnlist)):
    veclist2[i] = kcfco.kcf_vec(complist[i])
    print(time.time() - start)


# In[40]:

for vec1 in veclist:
    for vec2 in veclist:
        print(kcfco.similarity(vec1, vec2))
    print("\n")


# In[42]:

similarity_dict = dict()
for vec1 in veclist:
    for vec2 in veclist:
        similarity_dict[(vec1, vec2)] = kcfco.similarity(vec1, vec2)


# In[47]:

for vec1 in veclist:
    for vec2 in veclist:
        if similarity_dict[(vec1, vec2)][0] != similarity_dict[(vec2, vec1)][0]:
            print("BAD0")
        elif similarity_dict[(vec1, vec2)][1] != similarity_dict[(vec2, vec1)][2]:
            print("BAD1")
        elif similarity_dict[(vec1, vec2)][2] != similarity_dict[(vec2, vec1)][1]:
            print("BAD2")


# In[49]:

similarity_dict2 = dict()
for vec1 in veclist:
    for vec2 in veclist:
        similarity_dict2[(vec1, vec2)] = kcfco.similarity(vec1, vec2, levels=[0])

for vec1 in veclist:
    for vec2 in veclist:
        if similarity_dict2[(vec1, vec2)][0] != similarity_dict2[(vec2, vec1)][0]:
            print("BAD0")
        if similarity_dict2[(vec1, vec2)][1] != similarity_dict2[(vec2, vec1)][2]:
            print("BAD1")
        if similarity_dict2[(vec1, vec2)][2] != similarity_dict2[(vec2, vec1)][1]:
            print("BAD2")


# In[50]:

veclist[0].pandas()


# In[ ]:



