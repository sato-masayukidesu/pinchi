
# coding: utf-8

# In[62]:

import resource


# In[63]:

rsrc = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrc)


# In[64]:

print(soft, hard)


# In[65]:

soft = 256
resource.setrlimit(rsrc,(soft,hard))


# In[66]:

k = 0
d = dict()
for i in range(100):
    k += 1
    d[str(i)] = i
print(k)


# In[67]:

print(soft, hard)


# In[68]:

print(d, len(d))


# In[69]:

print(getattr(resource.getrusage(resource.RUSAGE_SELF),'ru_utime'))


# In[70]:

print(getattr(resource.getrusage(resource.RUSAGE_SELF),'ru_ixrss'))


# In[71]:

print(getattr(resource.getrusage(resource.RUSAGE_SELF),'ru_idrss'))


# In[72]:

resource.getpagesize()


# In[73]:

import gc


# In[74]:

gc.isenabled()


# In[77]:

gc.collect()


# In[ ]:



