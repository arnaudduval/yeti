Coupling modelling
******************

Details about coupling modelling can be found in :cite:`hirschler_embedded_2019`.

.. bibliography:: ../YETI.bib


Choosing degrees
================

Let p denote the smaller degree of both subdomain displacement fields.
User should adopt the following strategy:

- for the displacement constraint, a B-Spline function λ h with degree p − 1 is defined since it is mainly related to traction forces,
- for the rotation constraint, a B-Spline function μ h with degree p − 2 is defined because it transfers a bending moment,
- same mesh refinement is chosen for both Lagrange multipliers. User should discretize these fields using as many elements as the coarsest of the subdomains over the interface.
 
Those recommandations are only based on numerical experiments. With such a choice, we never encountered instabilities in our computations. 
