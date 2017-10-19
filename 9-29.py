
# coding: utf-8

# In[1]:

from classes import MCS_Finder


# In[2]:

f = MCS_Finder("Tephrosia")
Cnlist = f.get_Cnlist_from_label2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O")


# In[3]:

get_ipython().run_cell_magic('time', '', 'f.get_simcomp(Cnlist, "SIMCOMP2/Tephrosia1.txt")')


# In[4]:

print(len(Cnlist))


# In[5]:

import sys
sys.path.append("../../module")


# In[6]:

import kcf.converter as kcfco


# In[7]:

from nxrd.Compound import Compound


# In[8]:

f1 = MCS_Finder("Camellia")
Cnlist = f1.get_Cnlist_from_label2("C-C-C-N-C-N")


# In[9]:

complist = []
for i in range(len(Cnlist)):
    complist.append(Compound())
print(complist)
print(len(complist))


# In[10]:

for i, Cn in enumerate(Cnlist):
    complist[i].input_molfile('KNApSAck_mol/%s.mol' % (Cn))
    complist[i].fit2d = False
    print(Cn, i)


# In[11]:

veclist = [0] * len(Cnlist)


# In[12]:

import time
start = time.time()
for i in range(len(Cnlist)):
    print(i)
    print(time.strftime("%X"))
    veclist[i] = kcfco.kcf_vec(complist[i])
print(time.time() - start)


# 約0.9秒で終わった。  

# In[44]:

for vec1 in veclist:
    for vec2 in veclist:
        print(kcfco.similarity(vec1, vec2))
    print("\n")


# In[42]:

for vec1 in veclist:
    for vec2 in veclist:
        print(kcfco.similarity(vec1, vec2, levels=[0]))
    print("\n")


# In[16]:

get_ipython().magic('matplotlib inline')
f.make_and_draw_graph("C1-C1-C1-C8-C8-C8-C8-C8-C8-O", "SIMCOMP2/Tephrosia1.txt")


# In[17]:

G = f.make_graph2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O", "SIMCOMP2/Tephrosia1.txt")


# In[18]:

import networkx as nx
for g in nx.connected_component_subgraphs(G):
    f.drawnx(g)


# In[19]:

Glist = nx.connected_component_subgraphs(G)


# In[20]:

dic = dict()
for g in Glist:
    dic[g] = len(g.nodes())
sortGlist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
print(sortGlist)


# In[21]:

Cnlist1 = sortGlist[1][0].nodes()


# In[22]:

print(Cnlist1)


# In[23]:

pos = nx.spring_layout(sortGlist[1][0])


# In[24]:

def imscatter3(pos, Cnlist1, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    artists = []
    for Cn in Cnlist1:
        filepath = "Tephrosia/" + Cn + ".png"
        image = plt.imread(filepath)
        im = OffsetImage(image, zoom=zoom)
        x, y = np.atleast_1d(pos[Cn][0], pos[Cn][1])
        for x0, y0 in zip(x, y):
            ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
            artists.append(ax.add_artist(ab))
        ax.update_datalim(np.column_stack([x, y]))
        ax.autoscale()
    return artists


# In[25]:

import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Draw
from PIL import Image
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np


# In[26]:

for Cn in Cnlist1:
    with open('KNApSAck_mol/%s.mol' % Cn)as fi:
        mol = Chem.MolFromMolBlock(fi.read())
        rdDepictor.Compute2DCoords(mol)
        filename = "Tephrosia/%s.png" % Cn
        im = Draw.MolToImage(mol, size=(500, 500))
        trans = Image.new('RGBA', im.size, (0, 0, 0, 0))
        width = im.size[0]
        height = im.size[1]
        for x in range(width):
            for y in range(height):
                pixel = im.getpixel( (x, y) )
        
                # 白なら処理しない
                if pixel[0] == 255 and pixel[1] == 255 and pixel[2] == 255:
                    continue
        
                # 白以外なら、用意した画像にピクセルを書き込み
                trans.putpixel( (x, y), pixel )
        # 透過画像を保存
        trans.save(filename)


# In[27]:

fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('w')
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
ax.spines["right"].set_color("none")  # 右消し
ax.spines["left"].set_color("none")   # 左消し
ax.spines["top"].set_color("none")    # 上消し
ax.spines["bottom"].set_color("none") # 下消し
plt.xticks([])
plt.yticks([])
imscatter3(pos, Cnlist1, zoom=0.5)
nx.draw_networkx_nodes(sortGlist[1][0], pos, with_labels=False, node_size=200, node_shape=",", node_color="w", alpha=0)
nx.draw_networkx_edges(sortGlist[1][0], pos, with_labels=False, alpha=0.3)
plt.savefig("test.png")
plt.show()


# In[28]:

from classes import MCS_Finder


# In[29]:

f = MCS_Finder("Tephrosia")
Cnlist = f.get_Cnlist_from_label2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O")
G = f.make_graph2("C1-C1-C1-C8-C8-C8-C8-C8-C8-O", "SIMCOMP2/Tephrosia1.txt")


# In[30]:

f.make_skeleton(Cnlist)


# In[31]:

import networkx as nx


# In[32]:

pos = nx.spring_layout(G)


# In[33]:

f.imscatter3(pos, Cnlist)


# In[34]:

f.plot_with_figure(G, pos, Cnlist)


# In[35]:

Glist = nx.connected_component_subgraphs(G)
dic = dict()
for g in Glist:
    dic[g] = len(g.nodes())
sortGlist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
print(sortGlist)
Cnlist1 = sortGlist[1][0].nodes()


# In[36]:

g = sortGlist[1][0]


# In[37]:

pos = nx.spring_layout(g)


# In[38]:

f.plot_with_figure(g, pos, Cnlist1)


# In[ ]:




# In[39]:

f = MCS_Finder("Streptomyces")
Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")
f.get_simcomp(Cnlist, "SIMCOMP2/Streptomyces1.txt")
G = f.make_graph2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C", "SIMCOMP2/Streptomyces1.txt")


# In[40]:

Glist = nx.connected_component_subgraphs(G)
dic = dict()
for g in Glist:
    dic[g] = len(g.nodes())
sortGlist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
print(sortGlist)
Cnlist1 = sortGlist[1][0].nodes()


# In[41]:

g = sortGlist[1][0]
pos = nx.spring_layout(g)
f.plot_with_figure(g, pos, Cnlist1)


# In[ ]:



