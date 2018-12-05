#!/usr/bin/env python
# coding: utf-8

# # First approach to VAR estimation

# ## Libraries

# In[9]:


import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.api import VAR, DynamicVAR
from statsmodels.tsa.base.datetools import dates_from_str
import matplotlib


# ## Importing data

# In[2]:


# Importing data
mdata = sm.datasets.macrodata.load_pandas().data
mdata.head()


# In[3]:


# Basic data manipulation
dates = mdata[['year', 'quarter']].astype(int).astype(str)
quarterly = dates["year"] + "Q" + dates["quarter"]
quarterly = dates_from_str(quarterly)
mdata = mdata[['realgdp','realcons','realinv']]
mdata.index = pd.DatetimeIndex(quarterly)
data = np.log(mdata).diff().dropna()
data.head()


# ## VAR

# In[4]:


# VAR
model = VAR(data)
results = model.fit(2)
results.summary()


# ## Plots

# In[11]:


# Results

results.plot();


# In[12]:


# Autocrrelation

results.plot_acorr();


# In[ ]:




