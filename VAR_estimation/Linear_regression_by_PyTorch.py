#!/usr/bin/env python
# coding: utf-8

# In[1]:


import torch 
from torch.autograd import Variable 


# In[2]:


x_data = Variable(torch.Tensor([[1.0], [2.0], [3.0]])) 
y_data = Variable(torch.Tensor([[2.0], [4.0], [6.0]])) 


# In[3]:


class LinearRegressionModel(torch.nn.Module): 

	def __init__(self): 
		super(LinearRegressionModel, self).__init__() 
		self.linear = torch.nn.Linear(1, 1) # One in and one out 

	def forward(self, x): 
		y_pred = self.linear(x) 
		return y_pred 


# In[4]:


# our model 
our_model = LinearRegressionModel() 


# In[6]:


criterion = torch.nn.MSELoss(reduction='sum') 
optimizer = torch.optim.SGD(our_model.parameters(), lr = 0.01) 


# In[8]:


for epoch in range(500): 

	# Forward pass: Compute predicted y by passing 
	# x to the model 
	pred_y = our_model(x_data) 

	# Compute and print loss 
	loss = criterion(pred_y, y_data) 

	# Zero gradients, perform a backward pass, 
	# and update the weights. 
	optimizer.zero_grad() 
	loss.backward() 
	optimizer.step() 
	print('epoch {}, loss {}'.format(epoch, loss.data)) 


# In[9]:


new_var = Variable(torch.Tensor([[4.0]])) 
pred_y = our_model(new_var) 
print("predict (after training)", 4, our_model(new_var).data[0][0]) 

