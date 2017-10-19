
# coding: utf-8

# In[1]:

from simnet import *


# In[2]:

def temp2(mode=0):
    G = nx.Graph()
    pylab.figure(figsize=(10, 10))
    scoredic = dict()
    f = MCS_Finder("Streptomyces")
    Clist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")
    for Cnumber in Clist:
        filepath = "SIMCOMP/" + Cnumber + ".txt"
        if not os.path.exists(filepath):
            continue
        with open(filepath, "r") as fi:
            scorelist = []
            score = fi.readline().split("\t")
            score[1] = float(score[1][:-1])
            if mode == 1:
                while(score[1] > 0.9):
                    scorelist.append(score)
                    if score[0] in Clist and score[1] < 1:
                        G = nxappend(G, Cnumber, score[0], score[1])
                    score = fi.readline().split("\t")
                    score[1] = float(score[1][:-1])
            elif mode == 2:
                try:
                    while(True):
                        scorelist.append(score)
                        if score[0] in Clist and score[1] < 1:
                            G = nxappend(G, Cnumber, score[0], score[1])
                        score = fi.readline().split("\t")
                        score[1] = float(score[1][:-1])
                except IndexError:
                    pass
            else:
                for i in range(3):
                    scorelist.append(score)
                    if score[0] in Clist and score[1] < 1:
                        G = nxappend(G, Cnumber, score[0], score[1])
                    score = fi.readline().split("\t")
                    score[1] = float(score[1][:-1])
            scoredic[Cnumber] = scorelist
    print(str(len(G.nodes())) + "/" + str(len(Clist)))
    print("edge:" + str(len(G.edges())))
    return G
    for g in nx.connected_component_subgraphs(G):
        nx.draw_networkx(g)
        plt.xticks([])
        plt.yticks([])
        plt.show()


# In[3]:

G = temp2()


# In[4]:

Glist = list(nx.connected_component_subgraphs(G))


# In[5]:

dic = dict()
for g in Glist:
    dic[g] = len(g.nodes())


# In[6]:

sortGlist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
print(sortGlist)


# In[7]:

sortGlist[1][0].nodes()


# In[2]:

f = MCS_Finder("Streptomyces")
Clist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[3]:

print(Clist)


# In[11]:

get_ipython().run_cell_magic('time', '', 'i = 0\nfor Cnumber in Clist:\n    filepath = "SIMCOMP2/" + Cnumber + ".txt"\n    if not os.path.exists(filepath):\n        url = "http://rest.genome.jp/simcomp/"+ Cnumber + "/knapsack"\n        i += 1\n        print(url, i)')


# In[16]:

import time
import urllib


# In[28]:

get_ipython().run_cell_magic('time', '', 'i = 0\nt = time.time()\ndic = dict()\nfor Cnumber in Clist:\n    filepath = "SIMCOMP2/" + Cnumber + ".txt"\n    if not os.path.exists(filepath):\n        amalist = Clist.copy()\n        amalist.remove(Cnumber)\n        url = "http://rest.genome.jp/simcomp2/"+ Cnumber + "/"\n        for amali in amalist:\n            url += amali + "+"\n        url = url[:-1]\n        url += "/cutoff=0.1"\n        urllib.request.urlretrieve(url, filepath)\n        i += 1\n        print(url, i)\n        print(time.time() - t)\n        print(time.strftime("%X"))\n        time.sleep(60)\n        if i == 4:\n            break')


# In[26]:

get_ipython().run_cell_magic('time', '', 'time.sleep(60)')


# SIMCOMP2を使って類似度を新しく出した。  
# これなら約2秒なので30倍くらいの短縮になる。

# In[33]:

def temp3(mode=0):
    G = nx.Graph()
    pylab.figure(figsize=(10, 10))
    scoredic = dict()
    f = MCS_Finder("Streptomyces")
    Clist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")
    for Cnumber in Clist:
        filepath = "SIMCOMP2/" + Cnumber + ".txt"
        if not os.path.exists(filepath):
            continue
        with open(filepath, "r") as fi:
            scorelist = []
            score = fi.readline().split("\t")
            score[2] = float(score[2][:-1])
            if mode == 1:
                while(score[1] > 0.9):
                    scorelist.append(score)
                    if score[0] in Clist and score[2] < 1:
                        G = nxappend(G, Cnumber, score[0], score[1])
                    score = fi.readline().split("\t")
                    score[1] = float(score[1][:-1])
            elif mode == 2:
                try:
                    while(True):
                        scorelist.append(score)
                        if score[2] < 1:
                            G = nxappend(G, Cnumber, score[1], score[2])
                        score = fi.readline().split("\t")
                        score[2] = float(score[2][:-1])
                except IndexError:
                    pass
            else:
                for i in range(3):
                    scorelist.append(score)
                    if score[0] in Clist and score[1] < 1:
                        G = nxappend(G, Cnumber, score[0], score[1])
                    score = fi.readline().split("\t")
                    score[1] = float(score[1][:-1])
            scoredic[Cnumber] = scorelist
    print(str(len(G.nodes())) + "/" + str(len(Clist)))
    print("edge:" + str(len(G.edges())))
    return G
    for g in nx.connected_component_subgraphs(G):
        nx.draw_networkx(g)
        plt.xticks([])
        plt.yticks([])
        plt.show()


# In[35]:

G1 = temp2(2)
G2 = temp3(2)


# In[37]:

print(G1.nodes(data=True))
print(G1.edges(data=True))


# In[40]:

for n1, n2 in zip(G1.nodes(data=True), G2.nodes(data=True)):
    if n1 == n2:
        print("OK")
    else:
        print("dame")
        print(n1, n2)


# In[42]:

l1 = list()
l2 = list()
for e1, e2 in zip(G1.edges(data=True), G2.edges(data=True)):
    if e1 == e2:
        print("OK")
    else:
        print("dame")
        print(e1, e2)
        l1.append(e1)
        l2.append(e2)


# In[43]:

print(l1)


# In[80]:

flag = [1] * len(l2)
k = 0
for e1 in l1:
    for i, e2 in enumerate(l2):
        if e1 == e2 and flag[i]:
            # print(e1, "OK")
            flag[i] = 0
            break
        elif e1[0] == e2[1] and e1[1] == e2[0] and e1[2] == e2[2] and flag[i]:
            flag[i] = 0
            break
        elif e1[0] == e2[1] and e1[1] == e2[0] and flag[i]:
            flag[i] = 0
            print(e1[:2])
            print(e1[2], e2[2])
            print("weight違い")
            break
        elif e1[0] == e2[0] and e1[1] == e2[1] and flag[i]:
            flag[i] = 0
            print(e1[:2])
            print(e1[2], e2[2])
            print("weight違い")
            break
    else:
        print(e1)
        k += 1
print(k)
k = 0
for i , e2 in enumerate(l2):
    if flag[i]:
        print(e2)
        k += 1
print(k)


# In[48]:

print(len(l1))


# In[49]:

print(len(G1.edges(data=True)))


# In[63]:

print(len(G2.edges(data=True)))


# In[54]:

print(len(l2))


# 数週間前に取った時と値が全然違うやつがある  
# SIMCOMP2に入ってるやつが正しそう

# In[81]:

f = MCS_Finder("Streptomyces")


# In[82]:

all_Clist = f.get_Cnumber(f.get_html(f.gene))


# In[83]:

print(len(all_Clist))


# In[84]:

temp = all_Clist[:100]
print(len(temp))


# In[98]:

get_ipython().run_cell_magic('time', '', 'i = 0\nt = time.time()\ndic = dict()\nCnumber = "C00015228"\nfilepath = "SIMCOMP2/" + Cnumber + "test.txt"\nurl = "http://rest.genome.jp/simcomp2/"+ Cnumber + "/"\nfor tem in temp:\n    url += tem + "+"\nurl = url[:-1]\nurl += "/cutoff=0.01"\nurllib.request.urlretrieve(url, filepath)\ni += 1\nprint(url, i)\nprint(time.time() - t)\nprint(time.strftime("%X"))')


# この速度なら1900こまとめてやっても平気そう？

# In[99]:

get_ipython().run_cell_magic('time', '', 'i = 0\nt = time.time()\ndic = dict()\nCnumber = "C00015228"\nfilepath = "SIMCOMP2/" + Cnumber + "test.txt"\nurl = "http://rest.genome.jp/simcomp2/"+ Cnumber + "/"\nfor tem in all_Clist:\n    url += tem + "+"\nurl = url[:-1]\nurl += "/cutoff=0.01"\n# urllib.request.urlretrieve(url, filepath)\ni += 1\nprint(url, i)\nprint(time.time() - t)\nprint(time.strftime("%X"))')


# 長すぎた。  
# 適当な長さで切って後でまとめた方が良さそう。

# simnet.pyにSIMCOMP2を使うバージョンは入れておくべき。

# In[5]:

url = "http://rest.genome.jp/simcomp2/"
urlC = ""
for Cnumber in Clist:
    urlC += Cnumber + "+"
else:
    urlC = urlC[:-1]
url += urlC + "/" + urlC + "/cutoff=0.1"
print(url)


# In[6]:

import urllib


# In[8]:

get_ipython().run_cell_magic('time', '', 'urllib.request.urlretrieve(url, "SIMCOMP2/testtesttest.txt")')


# In[ ]:



