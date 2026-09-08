Changelog
=========


Version 0.1.2 - xxxx-xx-xx
--------------------------
 - Add OpenMP parallel build of stiffness matrix
 - Add GitHub workflows (build, test, publish packages on PyPI, build doc, publish on ReadTheDocs)
 - Add badges
 - Add compatibility with Python 3.13
 - Add :mod:`yeti_iga.pymfiga` subpackage: multiphysics IGA library covering explicit
   dynamics, space-time methods, fast diagonalization, and multi-model mortar assembly
 - Add optional install extras ``[viz]`` (matplotlib) and ``[fem]`` (meshpy) for
   postprocessing and FEM mesh generation respectively
 - Add 9 regression benchmarks for ``pymfiga`` (critical time step and space-time L2
   error, 4 geometries each)

Version 0.1.1 - 2025-05-07
--------------------------
- Remove extra terminal outputs
- Fix API documentation

Version 0.1.0 - 2025-04-25
--------------------------
- Initial release of YETI