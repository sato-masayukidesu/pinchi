
# coding: utf-8

# In[1]:

import sys
sys.path.append("../../module")
from nxrd.Compound import Compound
from classes import MCS_Finder
import kcf.converter as kcfco


# In[2]:

f = MCS_Finder("Streptomyces")


# In[3]:

f.get_html(f.genus)


# In[ ]:



