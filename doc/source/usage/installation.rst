************
Installation
************

Prerequisites
*************

yeti is currently only available for Linux environment.
For usage in a Microsoft Windows environment, *Windows Subsystem for Linux (WSL) can be used*

For a Ubuntu 20.04 distribution, the following packages are required :
 - python3
 - python3-venv
 - libpython3-dev
 - git
 - gfortran
 - gcc
 - cmake (minimum version 3.21 is required, a quick guide to build cmake 2.24 is available :ref:`here <buildcmake>`.



Get sources
***********
With git :

``git clone https://lamcosplm.insa-lyon.fr/plugins/git/yeti/yeti.git``

Create virtual python environment for yeti
******************************************

Virtual environnement will be created at the root of user home directory

``python3 -m venv ~/yeti-venv``

To activate ths environment :

``source ~/yeti-venv/bin/activate``

Install required python modules with pip :

``pip install --upgrade pip``

``pip install numpy matplotlib scipy nlopt``

Compile yeti
************

``cd yeti``

Create a build directory:

``mkdir build``

``cd build``

If not already loaded, load virtual python environment

``source ~/yeti-venv/bin/avtivate``

Configure with CMake :

``cmake ..``

Build yeti :

``make -j4``

yeti library will be located in ``build/lib/python``
