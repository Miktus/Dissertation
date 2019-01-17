#!/usr/bin/env python
# coding: utf-8

# # Estimation of simulated VAR by ML using PyTorch

# ## Libraries

# In[1]:


import matplotlib.pyplot as plt
import torch
from torch.autograd import Variable
import math


# ## Simulation of stationary time-series

# In[2]:


# Time-series dimension of simulation

nsample = 1000
nvar = 4


# In[3]:


# Create null data frames for storing data

df = torch.empty((nsample, nvar))
rho = torch.empty((nvar))


# In[4]:


# Simulation of processes

# Start values

for i in range(0, df.shape[1]):
    df[0,i] = torch.rand(1)
    
# Rho parameters, smaller than one in absolute value
for i in range(0, df.shape[1]):
    rho[i] = torch.rand(1)

# Create the AR(1) processes
for i in range(1,df.shape[0]):
    for j in range(0, df.shape[1]):
        df[i,j] = rho[j]*df[i-1, j] + torch.randn(1)

print(df)


# In[5]:


# Plots

f, axarr = plt.subplots(df.shape[1], sharex=True, sharey=True)
f.suptitle('Basic plots')

for i in range(0, df.shape[1]):
    axarr[i].plot(df[:,i].numpy())

# Bring subplots close to each other.
f.subplots_adjust(hspace=0)

# Hide x labels and tick labels for all but bottom plot.
for ax in axarr:
    ax.label_outer()


# ## VAR

# In[6]:


class VAR:
    """ 
    **** VECTOR AUTOREGRESSION (VAR) MODELS ****
    ----------
    Parameters
    data : np.array
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
        
        # Format the independent variables
        
        def diff(x, n):  
            """ Calculate the n-th order discrete difference
            """
            new_torch = torch.zeros(list(x.shape)[0]-n)
            if n == 0:
                new_torch = x
            else: 
                for i in range(list(x.shape)[0]-n):
                    new_torch[i] = x[i] - x[i+n]
            return new_torch

        # Difference data
        
        self.data = torch.t(torch.stack([diff(i, self.integ) for i in torch.t(self.data)]))
        self.T = self.data.shape[0]
        self.ylen = self.data.shape[1]
                                    
        """
        Y : torch.array
            Contains the length-adjusted time series (accounting for lags)
        """     

        self.Y = torch.t(self.data[self.lags:,])
        
    def _design(self):
        """ Creates a design matrix
        Z : np.array
        """ 
        
        Z = torch.ones(((self.ylen*self.lags+1), (self.T-self.lags)))

        row_count=1
        for lag in range(1, self.lags+1):
            for reg in range(self.ylen):
                Z[row_count, :] = self.data[:,reg][self.lags-lag:-lag]
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
        return torch.mm(torch.mm(self.Y,torch.t(Z)),torch.inverse(torch.mm(Z,torch.t(Z))))
    
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
        
        par = Variable(torch.ones(self.lags*(self.ylen**2) + self.ylen + self.ylen), requires_grad=True)
        
        coef = torch.reshape(par[0:self.lags*self.ylen**2], (self.ylen, self.lags*self.ylen))
        coef_mean = par[self.lags*self.ylen**2:self.lags*self.ylen**2+self.ylen]
        coef_var = torch.diag(par[self.lags*self.ylen**2+self.ylen:])    
        Z = self._design()[1:]

        Y_0 = torch.t(torch.t(self.Y) - coef_mean)
        Z_0 = torch.t(torch.t(Z) - coef_mean.repeat(self.lags))
        
        learning_rate = 1e-5
        
        optimizer = torch.optim.LBFGS(params = [par], lr=learning_rate, max_iter=25)
        
        def closure():
            # Before the backward pass, use the optimizer object to zero all of the
            # gradients for the Tensors it will update (which are the learnable weights
            # of the model)
            optimizer.zero_grad()
            
            # First way (without a constant term in the likelihood function):
            loss = -(- .5*self.Y.shape[1]*torch.log(torch.abs(torch.det(coef_var))) - .5*torch.trace(torch.mm(torch.mm(torch.t(Y_0 - torch.mm(coef,Z_0)),torch.inverse(coef_var)),Y_0 - torch.mm(coef,Z_0))))

#              TO ADD: 
#             Second way:
#             dist = torch.distributions.MultivariateNormal(torch.mm(coef,Z_0), coef_var)
#             loss = -torch.mean(dist.log_prob(Y_0))
            
            # Backward pass: compute gradient of the loss with respect to model parameters
            loss.backward(retain_graph = True)
            
            print(loss)
            
            return loss
        
        # Calling the step function on an Optimizer makes an update to its parameters
        
        for i in range(50):
            optimizer.step(closure)
            
        return(par)


# In[7]:


# Estimate VAR(p) by MLE

MLE_results = VAR(data = df, lags = 3, target = None, integ = 1).MLE()
MLE_results


# In[8]:


# For more clarity, let's consider the specific case of VAR(1)

MLE_results = VAR(data = df, lags = 1, target = None, integ = 1).MLE();


# In[9]:


# Better visualisation of the results from previous cell

par_names_MLE = ['A AR(1)','B to A AR(1)','C to A AR(1)','D to A AR(1)',
                 'B AR(1)','A to B AR(1)','C to B AR(1)','D to B AR(1)',
                 'C AR(1)','A to C AR(1)','B to C AR(1)', 'D to C AR(1)',
                 'D AR(1)','A to D AR(1)','B to D AR(1)', 'C to D AR(1)',
                 'A Constant','B Constant', 'C Constant', 'D Constant',
                 'Var(A)','Var(B)','Var(C)','Var(D)'] 

results_MLE = dict(zip(par_names_MLE, MLE_results.flatten()))
results_MLE


# In[10]:


# Estimate VAR(p) by OLS

OLS_results = VAR(data = df, lags = 3, target = None, integ = 3).OLS()
print(OLS_results)


# In[11]:


# For more clarity, let's consider the specific case of VAR(1)

OLS_results = VAR(data = df, lags = 1, target = None, integ = 1).OLS()

par_names_OLS = ['A Constant','A AR(1)','B to A AR(1)','C to A AR(1)','D to A AR(1)',
                 'B Constant','B AR(1)','A to B AR(1)','C to B AR(1)','D to B AR(1)',
                 'C Constant','C AR(1)','A to C AR(1)','B to C AR(1)', 'D to C AR(1)',
                 'D Constant','D AR(1)','A to D AR(1)','B to D AR(1)', 'C to D AR(1)'] 

results_OLS =  dict(zip(par_names_OLS, OLS_results.flatten()))
results_OLS

