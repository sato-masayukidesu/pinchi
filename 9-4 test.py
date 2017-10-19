
# coding: utf-8

# In[2]:

from classes import MCS_Finder
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit import Chem


# In[3]:

f = MCS_Finder("Streptomyces")


# In[4]:

Cnlist = f.get_Cnlist_from_label2("C(C)-C(C)-C(C-C-C-C-C-C)-C-C-C-C-C-C-C-C-C")


# In[5]:

print(Cnlist)


# In[6]:

for Cn in Cnlist:
    with open('KNApSAck_mol/%s.mol' % (Cn))as fi:
        mol = Chem.MolFromMolBlock(fi.read())
        rdDepictor.Compute2DCoords(mol)
    filename = f.gene + "/" + Cn + "back.png"
    Draw.MolToFile(mol, filename)
    break


# Cnlist内のCnの図を全て保存するためのコード

# とりあえず必要なのだけ抜く

# In[9]:

for Cn in ["C00017726", "C00026595", "C00026596", "C00015229", "C00015228"]:
    with open('KNApSAck_mol/%s.mol' % (Cn))as fi:
        mol = Chem.MolFromMolBlock(fi.read())
        rdDepictor.Compute2DCoords(mol)
    filename = f.gene + "/" + Cn + "back.png"
    Draw.MolToFile(mol, filename, size=(500, 500))


# In[22]:

get_ipython().magic('matplotlib inline')


# In[35]:

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
im = Image.open("Streptomyces/9-4test.png")
im_list = np.asarray(im)
plt.imshow(im_list)
plt.show()


# 左上の部分だけ画像を拾ってきた  
# 細かい部分は直接開いて

# In[34]:




# In[ ]:



