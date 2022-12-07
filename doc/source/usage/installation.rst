************
Installation
************

.. _yeti-installation:

Prerequisites
*************

yeti is currently only available for Linux environment.
For usage in a Microsoft Windows environment, *Windows Subsystem for Linux (WSL)* can be used.
A quick tutorial for the installation of WSL can be found :ref:`here <install-wsl>`.

For a Ubuntu 20.04 distribution, the following packages are required :
 - python3
 - python3-venv
 - libpython3-dev
 - git
 - gfortran
 - gcc
 - libblas3
 - libblas-dev
 - liblapack3
 - liblapack-dev
 - cmake (minimum version 3.21 is required, a quick guide to build cmake 2.24 is available :ref:`here <buildcmake>`)


Get sources
***********
With git :

..  code-block:: bash

    git clone https://lamcosplm.insa-lyon.fr/plugins/git/yeti/yeti.git

Create virtual python environment for yeti
******************************************

Virtual environnement will be created at the root of user home directory

.. code-block:: bash

    python3 -m venv ~/yeti-venv

To activate this environment:

.. code-block:: bash

    source ~/yeti-venv/bin/activate

Install required python modules with pip:

..  code-block:: bash

    pip install --upgrade pip
    pip install numpy matplotlib scipy nlopt


For use of umfpack with scipy sparse solver:

..  code-block:: bash

    pip install wheel
    sudo apt install libsuitesparse-dev swig
    pip install scikit-umfpack


Compile yeti
************

..  code-block:: bash

    cd yeti

Create a build directory:

..  code-block:: bash

    mkdir build
    cd build

If not already loaded, load virtual python environment:

..  code-block:: bash

    source ~/yeti-venv/bin/activate

Configure with CMake:

..  code-block:: bash

    cmake ..

Build yeti:

..  code-block:: bash

    make -j4

yeti library will be located in :file:`build/lib/python`. You must add it to your ``PYTHONPATH`` environment variable:

..  code-block:: bash

    export PYTHONPATH=$PYTHONPATH:~/yeti/build/lib/python

