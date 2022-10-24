..  _install-wsl:

Windows Subsystem for Linux installation
****************************************

This procedure has been tested on Windows 10, build 19044.2006

Start by opening a Windows Powershell with adminitrstor rights : 

.. image:: ../_static/install_WSL/01_launch_powershell_admin.png
   :width: 400
   :align: center

In this powershell, run the following command to install Windows Subsystem for Linux with the Ubuntu distribution : 

``wsl --install -d Ubuntu``

Force the use of WSL version 1 : 

``wsl --set-default-version 1``

Restart Windows.

For the start menu, launch Ubuntu application : 

.. image:: ../_static/install_WSL/02_launch_Ubuntu.png
   :width: 400
   :align: center

You will be asked to set a news user name and a password.

Your Windows Subsystem is now ready. You can continue by installing YETI.