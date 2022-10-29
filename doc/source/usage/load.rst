Loads
=====

Constant pressure
-----------------

In the :file:`*.inp` file, constant pressure is identified with keys ``Uxy``, where ``x`` is the index of the face (see :ref:`faces-numbering`) and ``y`` specify the load direction:

 - ``0`` for normal load
 - ``1``, ``2`` or ``3`` for load in respectively direction 1, 2 or 3

For example, a normal load of 60 per unit of surface on face 4 of elements contained in the set ``I1.EltToLoad`` is defined with : 

..  code-block::

    *Dload
    I1.EltToLoad, U40, 60.


Distributed pressure
--------------------

For distributed pressure, a pressure value should be defined for each control point of the mesh.

In the :file:`*.inp` file, distributed pressure is identified with key ``Ux4``, where ``x`` is the index of the face (see :ref:`faces-numbering`).

For example, a distributed pressure defined with a nodal distribution ``presField`` on face 6 of the elements contained in the set ``I1.EltToLoad`` is defined with : 

..  code-block::

    *Dload
    I1.EltToLoad, U64, presField

Define a nodal distribution
---------------------------

TO BE DONE
