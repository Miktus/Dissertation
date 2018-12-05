#!/usr/bin/env python
# coding: utf-8

# # First approach to VAR estimation

# ## Libraries

# In[1]:


import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.api import VAR, DynamicVAR
from statsmodels.tsa.base.datetools import dates_from_str
from statsmodels.tsa.vector_ar.vecm import coint_johansen
import matplotlib
from random import random


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

# Pandas object to a time-series model expects that the index is dates

mdata.index = pd.DatetimeIndex(quarterly)
data = np.log(mdata).diff().dropna()
data.head()


# In[4]:


# Plot the data 

data.plot(subplots=True,figsize=(18, 16));


# In[5]:


# Checking stationarity

coint_johansen(data,-1,1).eig


# ## VAR

# In[6]:


# VAR

model = VAR(data)
results = model.fit(2)
results.summary()


# ## Plots

# In[7]:


# Results

results.plot();


# In[8]:


# Autocorrelation

results.plot_acorr();


# ## Lag order selection

# In[9]:


# Adding maximum lag

model.select_order(5)

# Estimate model with maximum lag according to the AIC criterion

results = model.fit(maxlags=5, ic='aic')
results.summary()


# ## Forecasting

# In[10]:


lag_order = results.k_ar

# Specify the initial value of forecast

results.forecast(data.values[-lag_order:], 5);

# Plot forecasts with confidence intervals

results.plot_forecast(10);


# ## Impulse Response Analysis

# In[11]:


# IRF

irf = results.irf(10)
irf.plot(orth=False);


# In[12]:


# Cumulated IRF

irf.plot_cum_effects(orth=False);


# ## Forecast Error Variance Decomposition (FEVD)

# In[13]:


fevd = results.fevd(5)

fevd.summary()


# In[14]:


#Plot FEVD

results.fevd(20).plot();


# # Simulated part

# In[15]:


data_sim = data.copy()

for i in range(0,len(data_sim.columns)):
    data_sim.iloc[:,i] = random()


# In[16]:


model_sim = VAR(data)
results_sim = model_sim.fit(2)
results_sim.summary()


# In[17]:


results_sim.plot();

