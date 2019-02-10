
# # Estimation of simulated VAR by ML using PyTorch

# ## Libraries

import torch
from torch import mm, t, inverse, diag, exp, reshape, empty, log, det, trace
from torch.autograd import Variable
# import matplotlib.pyplot as plt
from plotly import tools
from plotly.offline import init_notebook_mode, iplot
import plotly.graph_objs as go
import plotly.io as pio

# Set seed

torch.manual_seed(123)

# Templates for graphs

pio.templates.default = 'plotly_dark'
init_notebook_mode(connected=True)

# Simulation of stationary time-series

# Time-series dimension of simulation

nsample = 1000
nvar = 4

# Create null data frames for storing data

df = empty((nsample, nvar))
rho = empty((nvar))

# Simulation of processes

# Start values

for i in range(0, df.shape[1]):
    df[0, i] = torch.rand(1)

# Rho parameters, smaller than one in absolute value
for i in range(0, df.shape[1]):
    rho[i] = torch.rand(1)

# Create the AR(1) processes
for i in range(1, df.shape[0]):
    for j in range(0, df.shape[1]):
        df[i, j] = rho[j] * df[i - 1, j] + torch.randn(1)

print(df)

# Plots

for i in range(0, df.shape[1]):
    to_plot = go.Scatter(
        y=df[:, i],
        x=list(range(0, df.shape[0])))
    layout = go.Layout(
        title=('Chart Number: ' + str(i)))
    figure = go.Figure(data=[to_plot], layout=layout)
    iplot(figure)

# Subplots version

fig = tools.make_subplots(rows=df.shape[1], cols=1)

for i in range(0, df.shape[1]):
    fig.append_trace(to_plot, i + 1, 1)


fig['layout'].update(height=600, width=600, title='Time series')
iplot(fig)


pio.write_image(fig, 'Graphs/Simulated_time_series.png')


# ## VAR


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

    def __init__(self, data, lags, target=None, integ=0):

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
            new_torch = torch.zeros(list(x.shape)[0] - n)
            if n == 0:
                new_torch = x
            else:
                for i in range(list(x.shape)[0] - n):
                    new_torch[i] = x[i] - x[i + n]
            return new_torch

        # Difference data

        self.data = t(torch.stack([diff(i, self.integ) for i in t(self.data)]))
        self.T = self.data.shape[0]
        self.ylen = self.data.shape[1]

        """
        Y : torch.array
            Contains the length-adjusted time series (accounting for lags)
        """

        self.Y = t(self.data[self.lags:, ])

    def _design(self):
        """ Creates a design matrix
        Z : np.array
        """

        Z = torch.ones(((self.ylen * self.lags + 1), (self.T - self.lags)))

        row_count = 1
        for lag in range(1, self.lags + 1):
            for reg in range(self.ylen):
                Z[row_count, :] = self.data[:, reg][self.lags - lag:-lag]
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
        return mm(mm(self.Y, t(Z)), inverse(mm(Z, t(Z))))

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

        par = Variable(torch.rand(self.lags * (self.ylen**2)
                                  + self.ylen + self.ylen), requires_grad=True)

        coef = reshape(par[0:self.lags * self.ylen**2],
                       (self.ylen, self.lags * self.ylen))
        coef_mean = par[self.lags * self.ylen **
                        2:self.lags * self.ylen**2 + self.ylen]
        coef_var = diag(exp(par[self.lags * self.ylen**2 + self.ylen:]))
        Z = self._design()[1:]

        Y_0 = t(t(self.Y) - coef_mean)
        Z_0 = t(t(Z) - coef_mean.repeat(self.lags))

        learning_rate = 1e-5
        n_iter = 100000

        optimizer = torch.optim.SGD(params=[par], lr=learning_rate)

        def closure():
            # Before the backward pass, use the optimizer object to zero all of the
            # gradients for the Tensors it will update (which are the learnable weights
            # of the model)
            optimizer.zero_grad()

            # Without a constant term in the likelihood function:

            loss = -(-.5 * self.Y.shape[1] * log(torch.abs(det(coef_var))) - .5 * trace(
                mm(mm(t(Y_0 - mm(coef, Z_0)), inverse(coef_var)), Y_0 - mm(coef, Z_0))))

            # Backward pass: compute gradient of the loss with respect to model parameters

            loss.backward(retain_graph=True)

            return loss, par

        # Calling the step function on an Optimizer makes an update to its parameters

        loss_vector = torch.zeros(n_iter)
        par_vector = torch.zeros(
            (n_iter, self.lags * (self.ylen**2) + self.ylen + self.ylen))

        for i in range(n_iter):
            optimizer.step(closure)
            loss_vector[i], par_vector[i, :] = closure()

        return(par, loss_vector, par_vector)

# VAR(1) by MLE


MLE_results = VAR(data=df, lags=1, target=None, integ=0).MLE()


# Better visualisation of the results from previous cell

par_names_MLE = ['A AR(1)', 'B to A AR(1)', 'C to A AR(1)', 'D to A AR(1)',
                 'B AR(1)', 'A to B AR(1)', 'C to B AR(1)', 'D to B AR(1)',
                 'C AR(1)', 'A to C AR(1)', 'B to C AR(1)', 'D to C AR(1)',
                 'D AR(1)', 'A to D AR(1)', 'B to D AR(1)', 'C to D AR(1)',
                 'A mean', 'B mean', 'C mean', 'D mean',
                 'logVar(A)', 'logVar(B)', 'logVar(C)', 'logVar(D)']

results_MLE = dict(zip(par_names_MLE, MLE_results[0].flatten()))
results_MLE

# VAR(1) by OLS

OLS_results = VAR(data=df, lags=1, target=None, integ=0).OLS()

par_names_OLS = ['A Constant', 'A AR(1)', 'B to A AR(1)', 'C to A AR(1)', 'D to A AR(1)',
                 'B Constant', 'B AR(1)', 'A to B AR(1)', 'C to B AR(1)', 'D to B AR(1)',
                 'C Constant', 'C AR(1)', 'A to C AR(1)', 'B to C AR(1)', 'D to C AR(1)',
                 'D Constant', 'D AR(1)', 'A to D AR(1)', 'B to D AR(1)', 'C to D AR(1)']

results_OLS = dict(zip(par_names_OLS, OLS_results.flatten()))
results_OLS

# Plotting loss function

loss_to_plot = go.Scatter(
    y=list(MLE_results[1].detach().numpy()),
    x=list(range(0, list(MLE_results[1].shape)[0])))
layout = go.Layout(
    title=('Loss function per iteration'))
figure = go.Figure(data=[loss_to_plot], layout=layout)
iplot(figure)

pio.write_image(figure, 'Graphs/Loss function per iteration.png')

# Plotting parameters


def plot_param(par_num):
    par_to_plot = go.Scatter(
        y=list(MLE_results[2][:, par_num].detach().numpy()),
        x=list(range(0, list(MLE_results[2].shape)[0])))
    layout = go.Layout(
        title=('Parameter ' + str(par_num) + 'per iteration'))
    par_per_iter = go.Figure(data=[par_to_plot], layout=layout)
    return iplot(par_per_iter)


plot_param(14)
