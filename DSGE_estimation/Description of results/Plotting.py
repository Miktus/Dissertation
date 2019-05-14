import pandas as pd
import matplotlib.pyplot as plt
from torch.distributions import uniform
import torch
# Import data

data = pd.read_csv("/Users/miktus/Documents/PSE/Dissertation/DSGE_estimation/Results/MC_All_F_100_100.csv")
# print(data)

# Optimal value

del data['Likelihood']
data.columns = ['Beta', 'Varsigma_p', 'Rho_G', 'Alpha', 'Tau']
histo = data.hist(figsize = (16,18), bins=25)

plt.show()
