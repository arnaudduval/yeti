future module
===============

``yeti_iga.future`` is a C++17 reimplementation of YETI's core B-spline/NURBS
machinery, exposed to Python through pybind11 and built independently of the legacy
Fortran/f2py layer documented under *API reference*. It is the direction of active
development for the library, currently centered on multipatch assembly (shared control
points, propagated refinement across patch interfaces, and consistent DOF numbering).

This section documents it in two complementary parts:

.. toctree::
   :maxdepth: 2

   design
   api
