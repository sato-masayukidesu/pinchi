
# coding: utf-8

# KNApSAckから化合物を持っている種族名を取り出したい。

# In[2]:

import lxml.html
import requests


# In[1]:

Cnumber = "C00017791"
genus = "Streptomyces"


# In[16]:

html = requests.get("http://kanaya.naist.jp/knapsack_jsp/information.jsp?word=" + Cnumber)


# In[17]:

dom = lxml.html.fromstring(html.text)


# In[22]:

print(dom.xpath('//*[@class="org2"]'))


# In[33]:

for element in dom.xpath('//*[@class="org2"]'):
    genus2 = element.text.split()[0]
    if genus2 == genus:
        print("OK", element.text)
    else:
        print(element.text)


# とりあえず引き抜けたし、ストレプトマイセスの判定もオッケー

# In[2]:

from classes import MCS_Finder


# In[3]:

f = MCS_Finder("Streptomyces")


# In[4]:

print(f.get_genuses(Cnumber))


# In[7]:

genuslist = f.get_genuses(Cnumber)


# In[15]:

li = []
for genus in genuslist:
    li.append(genus.replace("\xa0", " ")[:-1])
print(li)


# In[12]:

st = "test test"
li.append(st)
print(li)


# In[16]:

import imp


# In[19]:

imp.reload(classes.MCS_Finder)


# classesにも実装した。

# 今の所simcomp2で100こで2秒かかっている。  
# 1973こやるには1973*1972/2で1945378こやることになる。  
# だいたい38907秒 = 648分 = 10時間以上かかる。

# In[ ]:



