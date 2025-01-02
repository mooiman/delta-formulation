# delta-formulation
Numerics based on delta-formulation with fully implicit time-integration.

The fully implicit time integration is reached by an iterative Newton linearisation, but the behaviour of the Newton iteration depends on the smoothness of the data.
Therefor the data is adjusted in that way that the second derivatives of all data is smooth and is negiglible. 
This step in the procedure is called *regularization*.
And also the numerical scheme should be central in space, no dissipation is added to the model by the numerical method, just dispersion.
After the regularization and the central discretization is accurate (second order), fast (Newton quadratic convergence), efficient, robust (fully implicit). 
