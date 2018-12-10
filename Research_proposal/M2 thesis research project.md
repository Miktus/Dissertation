M2 thesis research project
Author : Michał Miktus

Supervisor : Pablo Winant

Dissertation title: Machine Learning approach to the DSGE modeling

The thesis will consist of the implementation of the Machine Learning algorithms to the estimation
of the Dynamic Stochastic General Equilibrium model in the framework of the small open
economy with frictions on the financial and labor markets.

More specifically, the thesis will start with the comprehensive summary of the current trends in the
DSGE modeling, with particular emphasis on the importance of the inclusion of financial and
unemployment frictions. The analysis of the growing literature of business cycle fluctuations
models will allow to choose one specific model, potentially the mix of the most renowned models,
possibly extended by the author by adding another aspects or by contrary simplified in order to
grasp the main ideas and still remain relatively medium in size and complexity.

The next step will incorporate solving of the selected model and then log-linearizing its solution
around the steady steady. Subsequently, parameter estimation based on the representation of the
equilibrium from the previous stage will be performed by implementing several techniques, not to
mention the maximum likelihood or the Bayesian interference. The author will heavily use the
Python Machine Learning libraries as TensorFlow or PyTorch in order to perform flexible numerical
computations and compare their efficiency, accuracy and time complexity with the classical
DSGE estimation procedures. By direct usage of the Machine Learning tools, the author hopes to
overcome several difficulties associated with the estimation of DSGE models using likelihood-
based methods, not to the mention the stochastic singularity of the covariance matrix, preventing
classical numerical routines based on the Hessian from working properly. Moreover, the author
hopes to prove that the machine learning algorithms can be considered as a valuable alternative
to the currently heavily used Bayesian estimation methods, usually Markov Chain Monte Carlo
(MCMC) techniques.

Finally, as the above-mentioned libraries allow for the effortless calculation of the derivatives by
implementing the automatic differentiation techniques, the author will perform several stress tests
by varying different parameters values, hoping that it will help to understand the fundamental
impact of the choice of the specific parameter to the final estimation.

References:

„Introducing Financial Frictions and Unemployment into a Small Open Economy Model”,
Lawrence J. Christiano, Mathias Trabandt, Karl Walentin, Journal of Economic Dynamics &
Control, 2011
https://pytorch.org/get-started/locally/
https://www.tensorflow.org
Author’s signature: Supervisor’s signature: