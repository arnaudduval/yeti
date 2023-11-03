Preprocessing module
====================

Temporary documentation for the content of IGAparametrization class.

``geometricSettings`` corresponds to the content of :file:`.NB` file

``mechanicalSettings`` corresponds to the content of :file:`.inp` file

``mechanicalSettings`` is a list with the following components :

 - [0] : parameters ???
 - [1] : boundary conditions
 - [2] : load
 - [3] : nodes
 - [4] : IEN connectivity
 - [5] : material properties
 - [6] : properties ???
 - [7] : tables, used for nodal distributions
 - [8] : shapeparametrization ????

``geometricSettings`` is a list with the following components :

 - [0] : dimension
 - [1] : sizes of knot vectors
 - [2] : knot vectors
 - [3] : degrees
 - [4] : knot spans
 - [5] : weights
 - [6] : number of elements per patch
 - [7] : number of nodes ???
 - [8] : number of elements (per patch ???)

Members of `IGAparametrization` class
-------------------------------------

Global model informations
~~~~~~~~~~~~~~~~~~~~~~~~~~

 - _ELT_TYPE : an array containing the element types for each patch (size : number of patchs)
 - _NBPINT : an array containing the number of integration points per element for each patch (size : number of patchs)
 - _TENSOR : an array containing the tensor type for each patch patch (size : number of patchs)
 - _mcrd : a scalar the number of coordinates of the model

Control points informations
~~~~~~~~~~~~~~~~~~~

 - _COORDS : a 2D array containing the coordinates of the control points (shape : 3, number of control points)
 - _nb_cp : a scalar containing the total number of control points

Elements informations
~~~~~~~~~~~~~~~~~~~~~

 - _IEN : a list containing the connectivity table for each patch (length : number of patchs). Each element of the list _IEN is a 2D array of shape (number of element of the patch, number of control points per element)

Material informations
~~~~~~~~~~~~~~~~~~~~~

 - _MATERIAL_PROPERTIES : a 2D array containing the material properties for each patch (shape : 3, number of patchs). The first dimension is 3 for elastic materials (Young modulus, Poisson ratio, density) but can be set higher for specific material behaviours (gradient elasticity, plasticity, ...)
 - _N_MATERIAL_PROPERTIES : an array containing the numebr of material properties for each patch (size : number of patchs)

Properties informations
~~~~~~~~~~~~~~~~~~~~~~~

 - _PROPS : a list containing arrays with the properties with real values for each patch (length : number of patchs). Each element of the _PROPS list can have a different size
 - _JPROPS : an array containing the number of properties for each patch (size : number of patchs)

Geometry informations
~~~~~~~~~~~~~~~~~~~~~

 - _dim : an array containing the dimension of each patch (size : number of patchs)
 - _Nkv : a 2D array containing the number of knots for each patch (shape : 3, number of patchs)
 - _Ukv : a list containing the knot vectors for each patch (length : number of patchs). Each element of the list _Ukv is a list of size 3 containing an array representing with the knot values of the patch in each parametric direction
 - _Jpqr : a 2D array containing the degrees of each patch in each direction (shape : 3, number of patchs)
 - _Nijk : a 2D array containing the knot span for each element in each direction (shape : 3, total number of elements)
 - _weight : a list containing the weights for each element (length : total number of elements). Each element of the list is an array whose size is number of control points of the element
 - _elementsByPatch : an array containing the number of elements for each patch (size : number of patchs)
 - _nb_patch : a scalar containg the number of patchs
 - _nnode : an array containing the number of control points per element for each patch (size : number of patchs)
 - _nb_elem : a scalar containing the total number of elements

Other informations
~~~~~~~~~~~~~~~~~~

 - _indCPbyPatch : an array containing the indices of control points for each patch (size : number of patchs). Each element of the array is an array containing the indices of the control points for a given patch (size : number of control points of the corresponding patch)




Class members (TEST)
====================

.. autoclass:: preprocessing.igaparametrization.IGAparametrization
    :members:
    :private-members:
