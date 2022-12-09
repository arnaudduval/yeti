======
Theory
======

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Sensibilities of embedded solid element
=======================================

Problem formulation
-------------------

We consider a domain :math:`\Omega \in \mathbb{R}^3` with a boudary :math:`\Gamma`.
This domain is subjected to body force :math:`g`, boundary traction :math:`t` and prescribed
displacement :math:`g`. Its material bheviour is isotropic and linear elastic. The problem can be
expressed as:

Find the displacements :math:`u : \Omega \Rightarrow \mathbb{R}^3` such that:

.. math::

    \begin{align}
        {\rm div}~\sigma + b = 0 & {\rm \quad in~} \Omega \\
        u = g & {\rm \quad on~} \Gamma_D \\
        \sigma \cdot n = t & {\rm \quad on~} \Gamma_N
    \end{align}

With :math:`\Gamma_D \cup \Gamma_N = \Gamma` and :math:`n` is the outside normal on the bounday :math:`\Gamma`.

The stress tensor :math:`\sigma` depends on the strain tensor :math:`\varepsilon` through the relation:

.. math::

    \sigma = D \varepsilon

where :math:`D` is the constitutive matrix of the linear elatic behaviour.

Galerkin weak form gives the following equation after integration:

.. math::

    \int_{\Omega} \delta \varepsilon^T \sigma \, {\rm d} \Omega - \int_{\Omega} \delta u^T \, b \, {\rm d} \Omega - \int_{\Gamma} \delta u \, t \,  {\rm d} \Gamma = 0

We set up a discrete formulation from the Galerkin weak form, given by:

.. math::
    
    K \, u = F

where :math:`K` is the the stiffness matrix, :math:`u` the discrete displacements defined at control points and :math:`F` the load vector.
Stiffness matrix and force vector can be assembled from elementary values :math:`K^e` and :math:`F^e`.

Elementary stiffness matrix is given by:

.. math::

    K^e_{kl} = \int_{\overline{\Omega^e}} B^T_k \, D \, B_l \left| J \right| \, {\rm d}\overline{\Omega}

where :math:`k` and :math:`l` are control point indices of the element, :math:`B_i` is the strain-displacement matrix, :math:`\left| J \right|` is the Jacobian determinant and :math:`\overline{\Omega^e}` is the parametric domain of the element.

Elementary force vector is given by:

.. math::

    F^e_k = \int_{\overline{\Omega^e}} N_k \, b \, \left| J \right| \, {\rm d}\overline{\Omega} + \int_{\overline{\Gamma^e_N}} N_k \, t  \left| J \right| \, {\rm d}\overline{\Gamma}

where :math:`\overline{\Gamma}` is the parametric domain of the boundary.

Sensitivities
-------------

In order to solve an optimization problem aimed at minimizing an objective function, we want to use a gradient-based algorithm.
To do this, we seek to express the gradients analytically.

We define an objective function :math:`f`;, function of the design variables :math:`x` and the state variable :math:`u`:

.. math::

    f := g \left( x, u\left( x \right) \right)

where u is the discrete solution of the linear system :math:`Ku=F`.