#!/usr/bin/env python
# coding: utf-8

# # Estimation of simulated VAR by ML

# ## Libraries

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import statsmodels.api as sm
import random as random
from statsmodels.tsa.arima_process import arma_generate_sample
import tsm as tsm
import pyflux as pf
from scipy.optimize import minimize

np.random.seed(1)


# ## Simulation of stationary time-series

# In[2]:


# Time-series dimension of simulation

nsample = 100

todays_date = datetime.datetime.now().date()
index = pd.date_range(end = todays_date, periods=nsample, freq='Q')
print(index)


# In[3]:


# Create null data frame

columns = ['A','B','C','D']
df = pd.DataFrame(index=index, columns=columns)
df = df.fillna(0) # With 0s rather than NaNs
print(df.head())


# In[4]:


# Simulation of ARMA processes

for i in np.arange(0,len(df.columns)):
    
    arparams = np.array([1, random.uniform(-0.99,0.99), random.uniform(-0.99,0.99)])
    maparams = np.array([1, random.uniform(-0.99,0.99), random.uniform(-0.99,0.99)])
    
    df.iloc[:,i] = arma_generate_sample(arparams, maparams, nsample)  

print(df.head())


# In[5]:


# Plots

df.plot(subplots=True, figsize=(30,20))


# ## VAR

# In[6]:


class VAR:
    """ 
    **** VECTOR AUTOREGRESSION (VAR) MODELS ****
    ----------
    Parameters
    data : pd.DataFrame or np.array
        Field to specify the time series data that will be used.
    lags : int
        Field to specify how many lag terms the model will have. 
    integ : int (default : 0)
        Specifies how many time to difference the dependent variables.
    target : str (pd.DataFrame) or int (np.array) (default : None)
        By default, all columns will be selected as the dependent variables.
    """

    def __init__(self,data,lags,target=None,integ=0):
        
        # Latent Variables
        self.lags = lags
        self.integ = integ
        self.target = target
        self.model_name = "VAR(" + str(self.lags) + ")"
        
        # Format the dependant variables
        self.data = data
        self.values = data.values
        self.index = data.index     
        self.data_name = data.columns.values
        self.T = data.shape[0]
        
        # Format the independent variables
        
        ## TO DO
        
        # Difference data
        X = np.transpose(self.data)
        for order in range(self.integ):
            X = np.asarray([np.diff(i) for i in X])
            self.data_name = np.asarray(["Differenced " + str(i) for i in self.data_name])
        self.data = X.T
        self.ylen = self.data_name.shape[0]
        
        """
        Y : np.array
            Contains the length-adjusted time series (accounting for lags)
        """     

        self.Y = self.data.drop(self.data.index[np.arange(0, self.lags)]).T
        
    def _design(self):
        """ Creates a design matrix
        Z : np.array
        """ 
        
        Z = np.ones(((self.ylen*self.lags+1), (self.T-self.lags)))

        row_count=1
        for lag in range(1, self.lags+1):
            for reg in range(self.ylen):
                Z[row_count, :] = self.data.ix[:,reg][(self.lags-lag):-lag]
                row_count += 1
                
        return(Z)

    def OLS(self):
        """ Creates OLS coefficient matrix
        ----------
        Parameters:
        NULL
        ----------
        Returns
        The coefficient matrix B
        """         
        
        Z = self._design()
        return np.dot(np.dot(self.Y,np.transpose(Z)),np.linalg.inv(np.dot(Z,np.transpose(Z))))
    
    def _neg_loglike(self, par):
        """ Calculate the MLE value, given the mean vector and variance matrix        
        """
        
        Z = self._design()[1:]
        
        coef = np.reshape(par[0:self.ylen**2], (self.ylen, self.ylen))
        coef_mean = par[self.ylen**2:self.ylen**2+self.ylen]
        coef_var = np.diag(par[self.ylen**2+self.ylen:])    
    
        Y_0 = (self.Y.T - coef_mean).T
        Z_0 = (Z.T - coef_mean).T 
        
        logLik = -self.Y.shape[1]*self.ylen*np.log(2*np.pi)*.5 - .5*self.Y.shape[1]*np.log(np.abs(np.linalg.det(coef_var)))- .5*np.trace(np.dot(np.dot((Y_0 - np.dot(coef,Z_0)).T,np.linalg.inv(coef_var)),Y_0 - np.dot(coef,Z_0)))
        
        return -logLik

    def MLE(self):
        """ Creates MLE coefficient matrix
        ----------
        Parameters:
        NULL
        ----------
        Returns
        The coefficient matrix MLE
        ----------
        It is based on the assumption of normality of errors
        """     
        
        cons = []
        
        for i in np.arange(self.ylen):
            cons.append(dict({'type':'ineq', 'fun': eval("lambda x: x[-" + str(i) + "]")}))
        
        # Make a list of initial parameter guesses 
        
        initParams = np.repeat(1, (self.ylen)**2 + self.ylen + self.ylen)

        # Run the minimizer
        results = minimize(self._neg_loglike, initParams, constraints = cons, method='COBYLA')

        # Print the results
        return(results.x)


# In[7]:


# Estimate VAR(1) by OLS

# In general

OLS_results = VAR(data = df, lags = 1, target = None, integ = 0).OLS()

# For more clarity

par_names_OLS = ['A Constant','A AR(1)','B to A AR(1)','C to A AR(1)','D to A AR(1)',
                 'B Constant','B AR(1)','A to B AR(1)','C to B AR(1)','D to B AR(1)',
                 'C Constant','C AR(1)','A to C AR(1)','B to C AR(1)', 'D to C AR(1)',
                 'D Constant','D AR(1)','A to D AR(1)','B to D AR(1)', 'C to D AR(1)'] 

results_OLS = pd.DataFrame(index = par_names_OLS, data = OLS_results.reshape(-1)) 
results_OLS = pd.DataFrame(index = par_names_OLS, data = OLS_results.reshape(-1))
results_OLS.columns = ['Parameters']
results_OLS


# In[8]:


# Check with the values for VAR(1) from the pyflux package

pf.VAR(df,1).fit().summary()


# In[9]:


# Estimate VAR(1) by MLE

MLE_results = VAR(data = df, lags = 1, target = None, integ = 0).MLE()
MLE_results

# For more clarity

par_names_MLE = ['A AR(1)','B to A AR(1)','C to A AR(1)','D to A AR(1)',
                 'B AR(1)','A to B AR(1)','C to B AR(1)','D to B AR(1)',
                 'C AR(1)','A to C AR(1)','B to C AR(1)', 'D to C AR(1)',
                 'D AR(1)','A to D AR(1)','B to D AR(1)', 'C to D AR(1)',
                 'A Constant','B Constant', 'C Constant', 'D Constant',
                 'Var(A)','Var(B)','Var(C)','Var(D)'] 

results_MLE = dict(zip(par_names_MLE, MLE_results.flatten()))
results_MLE


# In[10]:


# Check with the values for VAR(1) from the pyflux package

pf.VAR(df,1).fit(method='MLE').summary()

