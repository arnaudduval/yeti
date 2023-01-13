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

We consider a domain :math:`\Omega \in \mathbb{R}^3` with a boundary :math:`\Gamma`.
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

We define an objective function :math:`f`, function of the design variables :math:`x` and the state variable :math:`u`:

.. math::

    f := f \left( x, u\left( x \right) \right)

where u is the discrete solution of the linear system :math:`Ku=F`.

The derivative of the objectif function with respect to the design variables can be written as:

.. math::

    \frac{\mathrm{d} f}{\mathrm{d} x_i} = \frac{\partial f}{\partial x_i} + \frac{\partial f}{\partial u} \cdot \frac{\mathrm{d} u}{\mathrm{d} x_i}

We define :math:`u^*` as the solution of the adjoint problem :math:`K u^* = \frac{\partial g}{\partial u}`. The derivative of the objective function becomes:

.. math::

    \frac{\mathrm{d} f}{\mathrm{d} x_i} = \frac{\partial f}{\partial x_i} + u^* \cdot \left( \frac{\partial F}{\partial x_i} - \frac{\partial K}{\partial x_i} u \right)

From optimization model to analysis model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The derivatives of stiffness matrix :math:`K` and load vector :math:`F` with respect to the design variables are computed on the analysis model and propagated on the optimization model using the refinement operator :math:`R`, which is kept constant during the optimization process.

:math:`\tilde{Q}` and :math:`\tilde{P}` are the control points coordinates on respectively the analysis and optimization model. Derivatives can be expressed as:

.. math::

    \frac{\partial \bullet}{\partial \tilde{P}} = R^T \frac{\partial \bullet}{\partial \tilde{Q}}

Sensitivities of the objective function can be expressed as:

.. math::

    \frac{\mathrm{d} f}{\mathrm{d} x_i} = \frac{\partial f}{\partial x_i} + \frac{\partial \tilde{P}}{\partial x_i} : R^T \left(  u^* \cdot \frac{\partial F}{\partial \tilde{Q}} - u^* \cdot\frac{\partial K}{\partial \tilde{Q}} u \right)

where :math:`\frac{\partial \tilde{P}}{\partial x_i}` is an operator linking the control points of the optimization model with the design variables.

Total virtual work
~~~~~~~~~~~~~~~~~~

Let us introduce a function :math:`W` that takes as input arguments the state variable :math:`u`, the adjoint variable :math:`u^*` and the control points :math:`\tilde{Q}` of the analysis model:

.. math::

    \begin{align}
        W \left( u, u^*, \tilde{Q} \right) & = W_{\mathrm{ext}} \left( u^*, \tilde{Q} \right) + W_{\mathrm{int}} \left( u, u^*, \tilde{Q} \right) \\
        & = u^* \cdot F \left( \tilde{Q} \right) - u^* \cdot K \left( \tilde{Q} \right) u
    \end{align}

Sensitivities of the objective function becomes:

.. math::

    \frac{\mathrm{d} f}{\mathrm{d} x_i} = \frac{\partial f}{\partial x_i} + \frac{\partial \tilde{P}}{\partial x_i} : R^T \frac{\partial W}{\partial \tilde{Q}}

with:

.. math::

    \frac{\partial W}{\partial \tilde{Q}} = u^* \cdot \frac{\partial F}{\partial \tilde{Q}} - u^* \cdot\frac{\partial K}{\partial \tilde{Q}} u

Depending on the formulation of the response function :math:`f`, the partial derivative :math:`\frac{\partial f}{\partial x_i}` can be changed using the chain rule of differenciation. The sensitivities reads as:

.. math::

    \frac{\mathrm{d} f}{\mathrm{d} x_i} = \frac{\partial \tilde{P}}{\partial x_i} : R^T \left( \frac{\partial f}{\partial \tilde{Q}} + \frac{\partial W}{\partial \tilde{Q}} \right)

Differenciation of the IGA operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to to compute the derivatives of the internal and external works:

.. math::

    \begin{align}
        \frac{\partial W_{\mathrm{ext}}}{\partial \tilde{Q}} & = u^* \cdot \frac{\partial F}{\partial \tilde{Q}} \\
        \frac{\partial W_{\mathrm{int}}}{\partial \tilde{Q}} & = u^* \cdot \frac{\partial K}{\partial \tilde{Q}} u
    \end{align}


Embedded formulation
--------------------

.. list-table:: Notations
    :widths: 25 25 25
    :header-rows: 1

    * -
      - Embedded solid
      - Hull
    * - basis functions
      - :math:`R`
      - :math:`N`
    * - parameter
      - :math:`\theta`
      - :math:`\xi`
    * - physical points
      - :math:`\xi`
      - :math:`X`
    * - control points
      - :math:`P`
      - :math:`Q`

**In the embedded solid**, the coordinate :math:`\xi` is expressed as:

.. math::

    \xi = \sum_{a} R_a \left( \theta \right) P_a

And its derivative with respect to the coordinate :math:`\theta`:

.. math::

    \frac{\partial \xi}{\partial \theta} = \sum_{a} \frac{\partial R_a}{\partial \theta} P_a

and for specific directions :math:`i,j`:

.. math::
    \left( \frac{\partial \xi}{\partial \theta} \right)_{ij} = \frac{\partial \xi_i}{\partial \theta_j}

This quantity is stored in variable :code:`dxidtheta(i,j)`

Then, we derive this quantity with respect to the coordinates of a particular control point of the embedded entity :math:`P_a`:

.. math::

    \left( \frac{\partial}{\partial P_a} \left( \frac{\partial \xi}{\partial \theta} \right) \right)_{ijk} = \frac{\partial}{\partial P_{a_k}} \left( \frac{\partial \xi_i}{\partial \theta_j} \right) = \frac{\partial R_a}{\partial \theta_j} \delta_{ik}

This quantity is stored in variable :code:`DdxidthetaDP(i,j,k)`

The derivative of :math:`\frac{\partial \xi}{\partial \theta}` with respect tyo the coordinates of the hull's control points equals zero since this quantity does not depend on the control points :math:`Q`.

**In the hull**, the coordinate in the physical space is linked to the parametric coordinate :math:`\xi` by the relation:

.. math::

    X = \sum_{a} N_a \left( \xi \right) Q_a

And its derivative with respect to :math:`\xi`:

.. math::

    \frac{\partial X}{\partial \xi} = \sum_a \frac{\partial N_a}{\partial \xi} Q_a

And for specific directions :math:`i,j`:

.. math::

    \left( \frac{\partial X}{\partial \xi} \right)_{ij} = \frac{\partial X_i}{\partial \xi_j}

This quantity is stored in variable :code:`dxdxi(i,j)`

Its derivative with respect to a specific hull control point :math:`Q_a` reads:

.. math::

    \left( \frac{\partial}{\partial Q_a} \left( \frac{\partial X}{\partial \xi} \right) \right)_{ijk} = \frac{\partial}{\partial Q_{a_k}} \left( \frac{\partial X_i}{\partial \xi_j} \right) = \frac{\partial N_a}{\partial \xi_j} \delta_{ik}

This quantity is stored in variable :code:`DdxdxiDQ(i,j,k)`

To express the derivative with respect to embedded entity control point :math:`P_a`, we need to express the NURBS composition:

.. math::
    \frac{\partial X}{\partial \xi} = \sum_b \frac{\partial N_b \left( \sum_a R_a P_a \right)}{\partial \xi} Q_b

The derivative with respect to :math:`P_a` reads:

.. math::

    \left( \frac{\partial}{\partial P_a} \left( \frac{\partial X}{\partial \xi} \right)\right)_{ijk} = \frac{\partial}{\partial P_{a_k}} \left( \frac{\partial X_i}{\partial \xi_j} \right) = R_a \frac{\partial^2 X_i}{\partial \xi_j \partial \xi_k}

This quantity is stored in variable :code:`DdxdxiDP(i,j,k)`

**As a summary**

.. math::

    \begin{array}{c|c}
        \left( \frac{\partial}{\partial Q_a} \left( \frac{\partial \xi}{\partial \theta} \right) \right)_{ijk} = 0 & \left( \frac{\partial}{\partial Q_a} \left( \frac{\partial X}{\partial \xi} \right) \right)_{ijk} = \frac{\partial N_a}{\partial \xi_j} \delta_{ik} \\
        \left( \frac{\partial}{\partial P_a} \left( \frac{\partial \xi}{\partial \theta} \right) \right)_{ijk} = \frac{\partial R_a}{\partial \theta_j} \delta_{ik} & \left( \frac{\partial}{\partial P_a} \left( \frac{\partial X}{\partial \xi} \right) \right)_{ijk} = R_a \frac{\partial^2 X_i}{\partial \xi_j \partial \xi_k}
    \end{array}


Derivative of inverse mappings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this part, we express the derivative of inverse mapping :math:`\frac{\partial \xi}{\partial X}` and :math:`\frac{\partial \theta}{\partial \xi}`
with respect to a quantity named :math:`\Lambda` which can be eitehr the coordinates of control points of the hull or the embedded entity.

We start with:

.. math::

    \frac{\partial \xi}{\partial X} \cdot \frac{\partial X}{\partial \xi} = I

Which can be derived as:

.. math::

    \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial X} \right) \cdot \frac{\partial X}{\partial \xi} + \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial X}{\partial \xi} \right) = 0

Multiplying by :math:`\frac{\partial \xi}{\partial X}` gives:

.. math::

    \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial X} \right) = - \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial X}{\partial \xi} \right) \cdot \frac{\partial \xi}{\partial X}

The same reads for the derivative of :math:`\frac{\partial \theta}{\partial \xi}`:

.. math::

    \frac{\partial}{\partial \Lambda} \left( \frac{\partial \theta}{\partial \xi} \right) = - \frac{\partial \theta}{\partial \xi} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial \theta} \right) \cdot \frac{\partial \theta}{\partial \xi}

Derivative of the Jacobian determinant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several transformation to take into account:
 - Reference element space :math:`\overline{\xi}` to embedded entity parametric space :math:`\theta`
 - Embedded entity parametric space :math:`\theta` to hull parametric space :math:`\xi`
 - Hull parametric space :math:`\xi` to physical space :math:`X`

.. math::
    J = \frac{\partial X}{\partial \overline{\xi}} = \frac{\partial X}{\partial \xi} \cdot \frac{\partial \xi}{\partial \theta} \cdot \frac{\partial \theta}{\partial \overline{\xi}}

:math:`\frac{\partial \theta}{\partial \overline{\xi}}` does not depend on control points corrdinates. Thus, derivative of Jacobian determinant with respect to control points coordinates reads:

.. math::
    \frac{\partial \left| J \right|}{\partial \Lambda} = \frac{\partial}{\partial \Lambda} \left( \left| \frac{\partial X}{\partial \xi} \right| \right) \cdot \left| \frac{\partial \xi}{\partial \theta} \right| \cdot \left| \frac{\partial \theta}{\partial \overline{\xi}} \right| + \left| \frac{\partial X}{\partial \xi} \right| \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial \theta} \right) \cdot \left| \frac{\partial \theta}{\partial \overline{\xi}} \right|

Jacobi's formula give the expression of the derivative of a matrix determinant:

.. math::

    \frac{\mathrm{d}}{\mathrm{d} t} \mathrm{det} \left( A \left( t \right) \right) = \mathrm{det} \left( A \left( t \right) \right) \cdot \mathrm{tr} \left( A \left(t \right)^{-1} \cdot \frac{\mathrm{d} A \left( t \right)}{\mathrm{d} t}\right)

Applying Jocobi's formula top our case gives:

.. math::

    \begin{eqnarray}
        \frac{\partial}{\partial \Lambda} \left( \left| \frac{\partial X}{\partial \xi} \right|\right) & = & \left| \frac{\partial X}{\partial \xi} \right| \cdot \mathrm{tr} \left( \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial X}{\partial \xi} \right) \right) \\
        \frac{\partial}{\partial \Lambda} \left( \left| \frac{\partial \xi}{\partial \theta} \right|\right) & = & \left| \frac{\partial \xi}{\partial \theta} \right| \cdot \mathrm{tr} \left( \frac{\partial \theta}{\partial \xi} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial \theta} \right) \right)
    \end{eqnarray}

And the derivative of :math:`\left| J \right|` reads:

.. math::
    \begin{eqnarray}
        \frac{\partial \left| J \right|}{\partial \Lambda} & = & \left( \left| \frac{\partial X}{\partial \xi} \right| \cdot \mathrm{tr} \left( \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial X}{\partial \xi} \right) \right) \right) \cdot \left| \frac{\partial \xi}{\partial \theta} \right| \cdot \left| \frac{\partial \theta}{\partial \overline{\xi}} \right| \\
        && + \left| \frac{\partial X}{\partial \xi} \right| \cdot \left( \left| \frac{\partial \xi}{\partial \theta} \right| \cdot \mathrm{tr} \left( \frac{\partial \theta}{\partial \xi} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial \theta} \right) \right) \right) \cdot \left| \frac{\partial \theta}{\partial \overline{\xi}} \right| \\
        & = & \left| J \right| \cdot \left[ \mathrm{tr} \left( \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial X}{\partial \xi} \right) \right) + \mathrm{tr} \left( \frac{\partial \theta}{\partial \xi} \cdot \frac{\partial}{\partial \Lambda} \left( \frac{\partial \xi}{\partial \theta} \right) \right) \right]
    \end{eqnarray}

Applying this to the cases of control points :math:`P` and :math:`Q` gives:

.. math::

    \begin{eqnarray}
        \frac{\partial \left| J \right|}{\partial P} & = & \left| J \right| \cdot \left[ \mathrm{tr} \left( \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial P} \left( \frac{\partial X}{\partial \xi} \right) \right) + \mathrm{tr} \left( \frac{\partial \theta}{\partial \xi} \cdot \frac{\partial}{\partial P} \left( \frac{\partial \xi}{\partial \theta} \right) \right) \right] \\
        \frac{\partial \left| J \right|}{\partial Q} & = & \left| J \right| \cdot \mathrm{tr} \left( \frac{\partial \xi}{\partial X} \cdot \frac{\partial}{\partial Q} \left( \frac{\partial X}{\partial \xi} \right) \right)
    \end{eqnarray}