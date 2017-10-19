
# coding: utf-8

# In[2]:

import sys


# In[3]:

sys.path.append("../../module")


# In[4]:

from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco


# In[5]:

f = MCS_Finder("Streptomyces")
Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[6]:

complist = [] * len(Cnlist)
for i in Cnlist:
    complist.append(Compound())
print(complist)


# In[7]:

for i, Cn in enumerate(Cnlist):
    complist[i].input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    complist[i].fit2d = False
    print(Cn, i)


# In[8]:

veclist = [0] * len(Cnlist)


# In[9]:

import time
start = time.time()
for i in range(len(Cnlist)):
    veclist[i] = kcfco.kcf_vec(complist[i])
    print(time.time() - start)
print("finish")


# In[10]:

i = 0
for vec in veclist:
    if vec == 0:
        print("no file", i)
    else:
        print(vec, i)
    i += 1
print(len(Cnlist))


# In[11]:

for vec1 in veclist:
    for vec2 in veclist:
        print(kcfco.similarity(vec1, vec2))
    print("\n")


# In[12]:

simlistlist = []
for vec1 in veclist:
    simlist = []
    for i in range(len(Cnlist)):
        simlist.append((kcfco.similarity(vec1, veclist[i])[0], i))
    simlistlist.append(sorted(simlist, reverse=True))


# In[13]:

for simlist in simlistlist:
    print(simlist)
    print("\n")


# In[14]:

def nxappend(G, start, end, weight):
        if not start in G.nodes():
            G.add_node(start)
        if not end in G.nodes():
            G.add_node(end)
        if not (start, end) in G.edges() and not (end, start) in G.edges():
            G.add_edge(start, end, weight=weight)
        return G


# In[15]:

import networkx as nx


# In[16]:

G = nx.Graph()
for k in range(len(Cnlist)):
    for i in range(3):
        # scorelist.append(score)
        if simlistlist[k][i][0] < 1:
            G = nxappend(G, Cnlist[k], Cnlist[simlistlist[k][i][1]], simlistlist[k][i][0])


# In[17]:

get_ipython().magic('matplotlib inline')
f.drawnx(G)
print(len(G.nodes()))
print(len(G.edges()))


# In[18]:

# f.get_simcomp(Cnlist, "SIMCOMP2/Streptomyces1")
f.make_and_draw_graph("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C", "SIMCOMP2/Streptomyces1")


# In[19]:

print(Cnlist[13])


# In[20]:

veclist[0].pandas()


# In[21]:

complist[0].draw_cpd()


# In[22]:

list(veclist[0].pandas()["str"])


# In[23]:

temp1 = list(veclist[0].pandas()["str"])


# In[24]:

temp2 = []
for i in temp1:
    temp3 = i.split("-")
    temp2.append((len(temp3), temp3))


# In[25]:

print(temp2)


# In[26]:

print(temp2.sort())


# In[27]:

print(temp2)


# In[28]:

for i in temp2:
    if i[0] == 9:
        print(i[1])


# In[43]:

for i in temp2:
    if "O5" in i[1]:
        if i[0] >= 38:
            print("-".join(i[1]), i[0])


# In[54]:

counter = 0
node_counter = dict()
for i in temp2:
    for k in range(44):
        if i[0] == k:
            if len(i[1][0]) == 2:
                counter += 1
                node_counter[str(k)] = node_counter.get(str(k), 0) + 1


# In[55]:

node_counter["38"]


# In[56]:

"1" in node_counter


# In[64]:

for i in node_counter.keys():
    if 1 == node_counter[i]:
        print(i)


# In[66]:

for i in temp2:
    if i[0] == 19 or i[0] == 27:
        if len(i[1][0]) == 2:
            print("-".join(i[1]))


# リングの部分にハイフンが余分にあるからそこを取り除きたい  
# 正規表現使ったら抜けると思う      "-[0-9]"

# In[ ]:



