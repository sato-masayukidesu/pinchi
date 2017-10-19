
# coding: utf-8

# In[1]:

import networkx as nx


# In[4]:

G = nx.Graph()


# In[5]:

G.add_path([0,1,2])


# In[6]:

G.add_edge(2,3,weight=1.1)


# In[7]:

G.edges()


# In[8]:

G.add_edge(0,2,weight=2)


# In[10]:

G.edges(2, data=True)


# In[12]:

G.edges([0,2],data=True)


# In[3]:

G=nx.path_graph(4)
shells=[[0],[1,2,3]]
pos=nx.shell_layout(G,shells)


# In[6]:

get_ipython().magic('matplotlib inline')
nx.draw_networkx(G, pos)


# In[9]:

G.neighbors(0)


# In[10]:

G.add_edge(0, 3)


# In[11]:

G.neighbors(0)


# simnet.pyを今後も使うならパラメーターを外からいじれるようにしたほうがよさそう

# Streptomyces属の代謝経路がイマイチ見つからなかった。  
# のでちょっと他のに浮気

# In[1]:

from classes import MCS_Finder


# In[2]:

f = MCS_Finder("Camellia")


# In[3]:

f.make_kcfs()


# In[3]:

Cnlist = f.get_Cnlist_from_label2("C-N-C-C-N-C-N-C-N")


# In[4]:

print(Cnlist)


# In[5]:

f.make_image("C-N-C-C-N-C-N-C-N")


# 4つは少なすぎて厳しい

# In[6]:

Cnlist = f.get_Cnlist_from_label2("C-N-C-C-N")


# In[7]:

print(Cnlist)


# これも少ない

# In[8]:

Cnlist = f.get_Cnlist_from_label2("C-C-C-N-C-N")


# In[9]:

print(Cnlist)


# 一個増えた

# In[10]:

f.make_image("C-C-C-N-C-N")


# In[11]:

f.find_MCS_grid_image(Cnlist, "kk.png")


# In[12]:

print(Cnlist)


# In[14]:

print(f.mol_list)


# In[15]:

import time
for Cnumber in Cnlist:
    print(f.get_genuses(Cnumber, other = False))
    time.sleep(2)


# In[11]:

print(genus_list)


# 構造が小さいのが原因かな？

# もっと大きいので試したい

# In[8]:

f = MCS_Finder("Tephrosia")


# In[9]:

f.make_kcfs()


# In[10]:

Cnlist = f.get_Cnlist_from_label2("C-C-C-C-C-C-C-C-C-O")


# In[12]:

print(Cnlist, len(Cnlist))


# これぐらいのやつをsimnetにかけたい。
