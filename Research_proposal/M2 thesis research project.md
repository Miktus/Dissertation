M2 thesis research project
Author : Michał Miktus

Supervisor : Pablo Winant

Dissertation title: Machine Learning approach to the DSGE modeling

The dissertation will consist of the adaptation of the Machine Learning software stack to the estimation of the Dynamic Stochastic General Equilibrium model in the framework of the small open economy with frictions on the financial and labor markets. 

More specifically, the dissertation will start with the comprehensive summary of the current trends in the DSGE modeling, with particular emphasis on the importance of the inclusion of financial and unemployment frictions. The analysis of the growing literature of business cycle fluctuations models will allow to choose one specific model, potentially the mix of the most renowned models, possibly extended by the author by adding another aspects or by contrary simplified in order to grasp the main ideas and still remain relatively medium in size and complexity.

The next step will incorporate solving of the selected model and then log-linearizing its solution around the steady steady. Subsequently, parameter estimation based on the representation of the equilibrium from the previous stage will be performed by implementing the maximum likelihood procedure. The author will heavily use the Python Deep Learning libraries as TensorFlow or PyTorch in order to perform flexible numerical computations and compare their efficiency and time complexity with the classical DSGE estimation implementations. 

Moreover, as the above-mentioned libraries allow for the effortless calculation of the derivatives by the automatic differentiation techniques, the author will perform the sensitivity analysis by varying different parameters values, hoping that it will help to understand the fundamental impact of the choice of the specific parameter to the final estimation and its robustness.

Finally, the author hopes to come up with new efficient estimation procedure, which can benefit from the main strengths of the Deep Learning libraries, meaning the undemanding gradient calculations and vectorization.


References:

„Introducing Financial Frictions and Unemployment into a Small Open Economy Model”,
Lawrence J. Christiano, Mathias Trabandt, Karl Walentin, Journal of Economic Dynamics &
Control, 2011
https://pytorch.org/get-started/locally/
https://www.tensorflow.org
Author’s signature: Supervisor’s signature: