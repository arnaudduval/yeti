Elements library
================

Element U1
----------

U1 is a solid element

Element U3
----------

U3 is a Kirchhoff-Love shell element

Element U0
----------

U0 is a mapping element

Element U30
------------

U30 is an immersed Kirchhoff-Love shell element

 - Property #2 is the ID of the corresponding mapping patch
 - Property #3 is the shell thickness

Element U00
-----------

U00 is an interface element for couling 2 domains

 - Property #2 is the ID of the domain
 - Property #3 is the ID of the Lagrange multiplier
 - Property #4 indicate if the corresponding domain is master

Element U4
----------

U4 is a Lagrange multiplier
 - Property #2 indicates if it concerns displacement DOF (0) or totation DOF (1)
