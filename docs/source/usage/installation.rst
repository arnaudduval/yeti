Installation
============

YETI is distributed via the pypi package repository. As it is based on Fortran code compiled with f2py and using the Fortran module feature, it has so far only been built for Linux platforms.
Contributors with the knowledge to build YETI on MacOS or Windows platforms are welcome.

Linux
-----

You can simply install Yeti by invoking pip from the command line. YETI package is available for x86_64 platforms.
::

    pip install yeti-iga

Windows
-------

YETI is currently not available for Windows platforms. However, it can be installed using Windows Subsystem for Linux (WSL).


System Requirements
~~~~~~~~~~~~~~~~~~~

- Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11
- Administrator privileges on your Windows system
- Internet connection

Step-by-Step Installation
~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Open PowerShell as Administrator**

   Right-click on the Start menu and select **"Windows Terminal (Admin)"** or **"PowerShell (Admin)"**.

2. **Install WSL**

   Run the following command to install WSL and set up the default Linux distribution:

   .. code-block:: powershell

      wsl --install

   This command performs the following actions:
   - Enables the required optional components
   - Downloads and installs the latest Linux kernel
   - Sets WSL 2 as the default
   - Installs Ubuntu by default

3. **Restart Your Computer**

   After installation, you may be prompted to restart your PC. Make sure to save your work and reboot.

4. **Set Up Your Linux User Account**

   On the first launch of the Linux distribution, you will be asked to create a user account and password.

Verify Installation
~~~~~~~~~~~~~~~~~~~

To check that WSL is installed correctly and running:

.. code-block:: powershell

   wsl --status

This will display information about the installed distributions and the WSL version in use.


References
~~~~~~~~~~

- Official WSL documentation: https://docs.microsoft.com/windows/wsl/
- WSL GitHub repository: https://github.com/microsoft/WSL

