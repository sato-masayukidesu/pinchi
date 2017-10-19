
# coding: utf-8

# In[1]:

from simnet import *


# In[2]:

temp1()


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

print(Glist)


# In[6]:

dic = dict()
for g in Glist:
    dic[g] = len(g.nodes())


# In[7]:

print(dic)


# In[8]:

sortGlist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
print(sortGlist)


# In[9]:

for sg in sortGlist:
    nx.draw_networkx(sg[0])
    plt.xticks([])
    plt.yticks([])
    plt.show()


# In[10]:

sortGlist[1][0].nodes()


# In[11]:

sortGlist[1][0].edges()


# In[12]:

Cnlist = sortGlist[1][0].nodes()


# In[13]:

print(Cnlist)


# In[14]:

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
mollist = []
for Cn in Cnlist:
    with open('KNApSAck_mol/%s.mol' % (Cn))as f:
        mol = Chem.MolFromMolBlock(f.read())
        rdDepictor.Compute2DCoords(mol)
        mollist.append(mol)


# In[15]:

print(mollist)


# In[16]:

img = Draw.MolsToGridImage(mollist, legends=Cnlist, subImgSize=(400, 400))
img.save('Gridtest.png')


# エッジ情報が落ちてしまった

# In[17]:

pos = nx.spring_layout(sortGlist[1][0])


# In[18]:

print(pos)


# In[19]:

nx.draw_networkx(sortGlist[1][0], pos)
plt.show()


# In[20]:

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt


# In[21]:

for Cn in Cnlist:
    im = Image.open("Streptomyces/" + Cn + ".png")
    im_list = np.asarray(im)
    plt.imshow(im_list)
plt.show()


# In[22]:

im = Image.open("Streptomyces/C00015229.png")
im_list = np.asarray(im)
plt.imshow(im_list)
plt.show()


# In[23]:

plt.figure(figsize=(6, 8))
plt.axes([0.1, 0.1, 0.8, 0.8])
plt.imshow(im_list)
plt.axes([0.5, 0.6, 0.1, 0.1])
plt.imshow(im_list, alpha=0)
plt.show()


# imshowでもaxesが使えることがわかった。

# In[24]:

plt.axes([0.1, 0.1, 0.8, 0.8])
plt.xticks([])
plt.yticks([])
nx.draw_networkx(sortGlist[1][0], pos, with_labels=False)
for Cn in Cnlist:
    im=mpimg.imread("Streptomyces/" + Cn + ".png")
    lis = list(pos[Cn])
    for i in range(len(lis)):
        lis[i] = lis[i] * 0.64
        lis[i] = lis[i] + 0.1
    lis.extend([0.2, 0.2])
    plt.axes(lis)
    plt.xticks([])
    plt.yticks([])
    plt.imshow(im)
    if Cn == "C00017726":
        print(lis)

plt.show()


# In[25]:

plt.imshow(im_list, interpolation="nearest")
plt.show()


# In[26]:

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
img=mpimg.imread("Streptomyces/C00015229.png") #image to array
plt.imshow(img) #array to 2Dfigure

plt.show()


# とりあえずできたけどpngの解像度が悪すぎて見えない
# rdkitのせい

# In[27]:

with open('KNApSAck_mol/C00015228.mol'as fi:
    mol = Chem.MolFromMolBlock(fi.read())
    rdDepictor.Compute2DCoords(mol)
    filename = "Streptomyces/C00015228.png"
    Draw.MolToFile(mol, filename, size=(400, 400), transparent=True)


# In[28]:

from PIL import Image
with open('KNApSAck_mol/C00015228.mol')as fi:
    mol = Chem.MolFromMolBlock(fi.read())
    rdDepictor.Compute2DCoords(mol)
    filename = "Streptomyces/C00015228.png"
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


# 透過画像を作ることに成功した  
# これで拡大しても白い部分が邪魔にならない

# In[48]:

for Cn in Cnlist:
    with open('KNApSAck_mol/%s.mol' % Cn)as fi:
        mol = Chem.MolFromMolBlock(fi.read())
        rdDepictor.Compute2DCoords(mol)
        filename = "Streptomyces/%s.png" % Cn
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


# In[47]:

fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('w')
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
ax.spines["right"].set_color("none")  # 右消し
ax.spines["left"].set_color("none")   # 左消し
ax.spines["top"].set_color("none")    # 上消し
ax.spines["bottom"].set_color("none") # 下消し
plt.xticks([])
plt.yticks([])
nx.draw_networkx(sortGlist[1][0], pos, with_labels=False, node_size=700)
for Cn in Cnlist:
    im=mpimg.imread("Streptomyces/" + Cn + ".png")
    lis = list(pos[Cn])
    for i in range(len(lis)):
        lis[i] = lis[i] * 0.64
        lis[i] = lis[i] + 0.1 - (0.25 * 0.64)
    lis.extend([0.5, 0.5])
    ax = plt.axes(lis)
    ax.patch.set_alpha(0.1)
    ax.spines["right"].set_color("none")  # 右消し
    ax.spines["left"].set_color("none")   # 左消し
    ax.spines["top"].set_color("none")    # 上消し
    ax.spines["bottom"].set_color("none") # 下消し
    plt.xticks([])
    plt.yticks([])
    plt.imshow(im)

plt.show()


# 一応できたけど透過の具合が変な感じだから調整しないと醜い  
# あとファイルごと保存しないと結局拡大できない

# In[ ]:




# ## 9月6日分

# In[31]:

import matplotlib.path as pth


# In[32]:

pth.Path(vertices=im)


# In[ ]:

nx.draw_networkx(sortGlist[1][0], pos, with_labels=False, node_shape=im)


# matplotlib.pathはmatplotlib.path.Pathという軽い図形を描くためのクラスを描くためのもの  
# PATHを指定して表示できるわけではない

# In[33]:

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data


# In[34]:

def main():
    x = np.linspace(0, 10, 20)
    y = np.cos(x)
    image_path = 'Streptomyces/C00015228.png'
    fig, ax = plt.subplots()
    imscatter(x, y, image_path, zoom=0.5, ax=ax)
    ax.plot(x, y)
    plt.show()


# In[35]:

def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists


# In[36]:

main()


# これならmatplotlibのマーカーを画像に変更できる

# In[37]:

def imscatter2(pos, Cn, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    filepath = "Streptomyces/" + Cn + ".png"
    image = plt.imread(filepath)
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(pos[Cn][0], pos[Cn][1])
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists


# In[38]:

imscatter2(pos, "C00015228", zoom = 0.4)


# In[39]:

nx.draw_networkx(sortGlist[1][0], pos, with_labels=False, node_size=20)


# In[40]:

plt.show()


# In[41]:

def imscatter3(pos, Cnlist, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    artists = []
    for Cn in Cnlist:
        filepath = "Streptomyces/" + Cn + ".png"
        image = plt.imread(filepath)
        im = OffsetImage(image, zoom=zoom)
        x, y = np.atleast_1d(pos[Cn][0], pos[Cn][1])
        for x0, y0 in zip(x, y):
            ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
            artists.append(ax.add_artist(ab))
        ax.update_datalim(np.column_stack([x, y]))
        ax.autoscale()
    return artists


# In[46]:

fig = plt.figure(figsize=(8, 6))
fig.patch.set_facecolor('w')
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
ax.spines["right"].set_color("none")  # 右消し
ax.spines["left"].set_color("none")   # 左消し
ax.spines["top"].set_color("none")    # 上消し
ax.spines["bottom"].set_color("none") # 下消し
plt.xticks([])
plt.yticks([])
imscatter3(pos, Cnlist, zoom=0.5)
nx.draw_networkx_nodes(sortGlist[1][0], pos, with_labels=False, node_size=200, node_shape=",", node_color="w", alpha=0)
nx.draw_networkx_edges(sortGlist[1][0], pos, with_labels=False, alpha=0.3)
plt.show()


# できた。  
# 透過だけなんとかする。

# In[43]:

def imscatter4(pos, Cnlist, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    artists = []
    for Cn in Cnlist:
        filepath = "Streptomyces/" + Cn + "back.png"
        image = plt.imread(filepath)
        im = OffsetImage(image, zoom=zoom)
        x, y = np.atleast_1d(pos[Cn][0], pos[Cn][1])
        for x0, y0 in zip(x, y):
            ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
            artists.append(ax.add_artist(ab))
        ax.update_datalim(np.column_stack([x, y]))
        ax.autoscale()
        nx.draw_networkx_nodes(G=sortGlist[1][0], pos=pos, nodelist=[Cn], node_size=0)
    return artists


# In[44]:

imscatter4(pos, Cnlist, zoom=0.4)
plt.show()


# 背景が白い画像も用意した。これの余白？を切り取って貼り付けたい  
# 果たしてそれがいいことかは知らない

# パラメータの調整をして少し見やすくした。

# In[45]:

print(pos)


# In[ ]:



