
# coding: utf-8

# In[1]:

from classes import MCS_Finder


# In[2]:

f = MCS_Finder("Tephrosia")


# In[3]:

Cnlist = f.get_Cnlist_from_label2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O")


# In[4]:

print(Cnlist, len(Cnlist))


# In[5]:

get_ipython().run_cell_magic('time', '', 'import time\ns = time.time()\nfor i, Cn in enumerate(Cnlist):\n    f.get_molfile([Cn])\nprint(time.time() - s)')


# In[10]:

f.make_image("C1-C1-C1-C8-C8-C8-C8-C8-C8-O")


# In[9]:

list_list = []
for Cn in Cnlist:
    list_list.append(f.get_genuses(Cn, other=False))


# In[11]:

print(list_list, len(list_list))


# In[20]:

genuses = set()
for lis in list_list:
    for i in lis:
        genuses.add(" ".join(i.split()[:2]))
        # genuses.add(i)
print(genuses, len(genuses))


# In[22]:

genuses2 = set()
for lis in list_list:
    for i in lis:
        genuses2.add(i)
print(genuses2, len(genuses2))


# 変わってないっぽい？

# 34化合物あって25種のTephrosiaに含まれている。

# In[2]:

f = MCS_Finder("Streptomyces")


# In[10]:

get_ipython().magic('matplotlib inline')
f.make_graph("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# とりあえずsimnet初期版はclassesに入った。  
# この後いろんな変数を追加する予定。  
# スコアの閾値とかmatplotlibのパラメーターとか

# In[3]:

f = MCS_Finder("Tephrosia")


# In[4]:

Cnlist = f.get_Cnlist_from_label2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O")


# In[ ]:

f.get_simcomp(Cnlist, "SIMCOMP2/Tephrosia1.txt")


# In[5]:

url = "http://rest.genome.jp/simcomp2/"
urlC = ""
for Cnumber in Cnlist:
    urlC += Cnumber + "+"
else:
    urlC = urlC[:-1]
    url += urlC + "/" + urlC + "/cutoff=0.1"


# In[6]:

print(url)


# In[7]:

import requests


# In[1]:

html = requests.get(url)


# In[ ]:

print(html)


# pythonが変なので明日やる。  
# とりあえずsimnetを書いてみる。

# In[ ]:



