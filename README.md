# delta-formulation
Numerics based on delta-formulation with fully implicit time-integration.

The fully implicit time integration is reached by an iterative Newton linearisation, but the behaviour of the Newton iteration depends on the smoothness of the data.
Therefor the data is adjusted in that way that the second derivatives of all data is smooth and is negiglible. 
This step in the procedure is called *regularization*.
And the numerical scheme is chosen to be central in space, no dissipation is added to the model by the numerical method, just higher order (third-order) dispersion.
After the regularization and the central discretization the numerical method is accurate (second order), reliable (numerical are small compared to erros in the mathematical model), efficient (Newton quadratic convergence), robust (fully implicit). 
