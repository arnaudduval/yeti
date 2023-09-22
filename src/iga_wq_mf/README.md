# IGA-WQ-MF analysis applied to nonlinear elastoplasticity and heat transfer problems
### Developped by J. Cornejo-Fuentes

In this repository, we develop isogeometric analysis (IGA) with weighted─quadrature (WQ) and matrix─free (MF) approaches. This code is mainly based on algorithms/theory described at:

**Papers**:

About tensor product and sum-factorization
- *Tensor decompositions and applications* by T. Kolda and B. Bader
- *Efficient matrix computation for tensor-product isogeometric analysis* by P. Antolin et al.

About iterative solvers
- *Templates for the solution of linear systems: building blocks for iterative solvers* by R. Barrett et al.
- *A note on relaxed and flexible GMRES* by L. Giraud et al.

About weighted quadrature
- *Fast formation of isogeometric Galerkin matrices by weighted quadrature* by F. Calabro et al.
- *Fast formation and assembly of finite element matrices with application to isogeometric linear elasticity* by R. Hiemstra et al.
- *Solving Poisson's equation using dataflow computing* by R. van Nieuwpoort

About matrix free
- *Matrix-free weighted quadrature for a computationally efficient isogeometric k-method* by G. Sangalli and M. Tani

About fast diagonalization
- *Preconditioners for Isogeometric Analysis thesis* by M. Montardini
- *Isogeometric preconditioners based on fast solvers for the Sylvester equation* by G. Sangalli and M. Tani

About multipatch
- *A domain decomposition method ofr Isogeometric multi-patch problems with exact local solvers* by M. Bosy et al.
- *IETI-Isogeometric tearing and interconnecting* by S. Kleiss et al.

Others
- *B free* by J. Planas, I. Romero and J.M. Sancho
- *Solution algorithms for nonlinear transient heat conduction analysis employing element-by-element iterative strategies* by J. Winget and T.J.R. Hughes

**Books**:

To understand B-splines and NURBS
- *The NURBS book* by L. Piegl
- *Isogeometric analysis: toward integration of CAD and FEA* by J. Cotrell, T.J.R. Hughes and Y. Bazilevs

To understand the Galerkin method (heat conduction and elasticity)
- *The finite element method* by T.J.R. Hughes
- *Finite element simulation of heat transfer* by J.M. Bergheau and R. Fortunier

To undestand plasticity theory
- *Computational inelasticity* by J.C. Simo and T.J.R. Hughes
- *Introduction to nonlinear finite element analysis* by N.H. Kim

## How to use it ?
This code is part of YETI project. More information at https://aduval.fr/yeti/.  