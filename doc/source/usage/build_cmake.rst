..  _buildcmake:

Quick guide to build cmake (version 3.24) on Linux
==================================================

YETI requires cmake 3.21 as minimal version to be built.
Here is a quick guide to build and install cmake 3.24 on a Linux computer.
User needs to have **sudo** rights elevation.

Prerequisites
-------------

Some packages are necessary for compilation:

``sudo apt install gcc g++ gfortran libssl-dev make``

Get sources and uncompress
--------------------------
``wget https://github.com/Kitware/CMake/releases/download/v3.24.2/cmake-3.24.2.tar.gz``

``tar xzf cmake-3.24.2.tar.gz``

``cd cmake-3.24.2``

Configure, compile and install
------------------------------

``./configure``

``make``

``sudo make install``