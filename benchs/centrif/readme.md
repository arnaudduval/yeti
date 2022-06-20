Benchs for validation of centrifugal body force. Geometry is the same for all cases and reference solution was computed with Abaqus
Validation is made by comparing displacement at a corner of the geomeetry in a quasi raidal direction. Tolerance is set to 2%

 - centrif\_U1 : a single U1 patch
 - centrif\_U10 : a single U10 patch
 - centrif\_U1\_cpl\_U1\_1load : 2 U1 patchs with weak mortar coupling, centrifugal force is applied through a single instruction for all elements
 - centrif\_U1\_cpl\_U1\_2loads : 2 U1 patchs with weak mortar coupling, centrifugal force is applied through 2 instructions (one for each patch)
 - centrif\_U1\_cpl\_U10\_1load : a U1 patch and a U10 patch  with weak mortar coupling, centrifugal force is applied through a single instruction for all elements
 - centrif\_U1\_cpl\_U10\_2loads : a U1 patch and a U10 patch with weak mortar coupling, centrifugal force is applied through 2 instructions (one for each patch) 
 - centrif\_U1\_U1\_1load : 2 U1 patchs with strong coupling, centrifugal force is applied through a single instruction for all elements
 - centrif\_U1\_U1\_2loads : 2 U1 patchs with strong coupling, centrifugal force is applied through 2 instructions (one for each patch)
 - centrif\_U1\_C0 : a single U1 patch with a C0 surface, equivelent to the centrif\_U1\_U1\_1load case
