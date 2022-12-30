Preprocessing module
====================

Temporary documentation for the content of IGAparametrization class.

``geometricSettings`` corresponds to the content of :file:`.NB` file

``mechanicalSettings`` corresponds to the content of :file:`.inp` file

``mechanicalSettings`` is a list with the following components:

 - [0] : parameters ???
 - [1] : boundary conditions
 - [2] : load
 - [3] : nodes
 - [4] : IEN connectivity
 - [5] : material properties
 - [6] : properties ???
 - [7] : tables ???
 - [8] : shapeparametrization ????

 ``geometricSettings`` is a list with the following components:

 - [0] : dimensions of patchs
 - [1] : number of compoennts of knot vectors
 - [2] : components of knot vectors
 - [3] : degrees
 - [4] : knot spans
 - [5] : weights
 - [6] : number of elements per patch
 - [7] : number of patchs
 - [8] : number of nodes (total or per patch ???)
 - [9] : number of elements (total or per patch ???)





.. autoclass:: preprocessing.igaparametrization.IGAparametrization
    :members:
