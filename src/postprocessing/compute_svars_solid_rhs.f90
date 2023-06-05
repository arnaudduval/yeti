!! Copyright 2021 Arnaud Duval

!! This file is part of Yeti.
!!
!! Yeti is free software: you can redistribute it and/or modify it under the terms 
!! of the GNU Lesser General Public License as published by the Free Software 
!! Foundation, either version 3 of the License, or (at your option) any later version.
!!
!! Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
!! PURPOSE. See the GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License along 
!! with Yeti. If not, see <https://www.gnu.org/licenses/>

!! Compute RHS for least square projection of values at Gauss point
!! Input is a displacement solution 
!! stresses and strains are computed at integrations points
!! TODO make the same routine with an input of values already existing at integration points
!! TODO RHS could be build with Gram matrix
subroutine compute_svars_solid_rhs(rhs, sol, activeElement, coords3D, ien, &            
            &   nb_elem_patch, nkv, ukv, nijk, weight, jpqr, elt_type, &
            &   tensor, props, jprops, material_properties, nbInt, nnode, &
            &   nb_patch, nb_elem, nb_cp, mcrd)
    
    use parameters
    use nurbspatch
    use embeddedMapping
    
    implicit none
    
    !! Input arguments
    !! ---------------
    integer, intent(in) :: nb_cp
    integer, intent(in) :: mcrd
    integer, intent(in) :: nb_patch
    integer, intent(in) :: nb_elem
    double precision, intent(in) :: sol     !! analysis solution
    dimension sol(mcrd, nb_cp)
    integer, intent(in) :: activeElement
    dimension activeElement(nb_elem)
    double precision, intent(in) :: coords3D
    dimension coords3D(3, nb_cp)
    
    integer, intent(in) :: ien
    dimension ien(:)
    integer, intent(in) :: nb_elem_patch
    dimension nb_elem_patch(nb_patch)
    integer, intent(in) :: nkv
    dimension nkv(3, nb_patch)
    double precision, intent(in) :: ukv
    dimension ukv(:)
    integer, intent(in) :: nijk
    dimension nijk(3, nb_elem)
    double precision, intent(in) :: weight
    dimension weight(:)
    integer, intent(in) :: jpqr
    dimension jpqr(3, nb_patch)
    character(len=*), intent(in) :: elt_type
    character(len=*), intent(in) :: tensor
    double precision, intent(in) :: props
    dimension props(:)
    integer, intent(in) :: jprops
    dimension jprops(nb_patch)
    double precision, intent(in) :: material_properties
    dimension material_properties(2, nb_patch)
    
    integer, intent(in) :: nbInt           !! Number of integration points per element for each patch
    dimension nbInt(nb_patch)
    integer, intent(in) :: nnode
    dimension nnode(nb_patch)
    
    !! Output variables
    !! ----------------
    double precision, intent(out) :: rhs    !! output RHS vector
    dimension rhs(2*2*mcrd, nb_cp)
    
    !! Local variables
    !! ---------------
    integer :: numPatch                     !! loop variable for patches
    integer :: numElem                      !! lopp variable for elements
    integer :: jelem                        !! element number counter
    double precision :: coords_elem         !! CP coordinates of current element
    dimension coords_elem(mcrd, maxval(nnode))
    integer :: sctr
    dimension sctr(maxval(nnode))
    double precision elRHS
    dimension elRHS(2*2*mcrd, maxval(nnode))
    integer :: i, icp
    
    
    jelem = 0
    !! Loop on patches
    do numPatch = 1, nb_patch
        call extractNurbsPatchGeoInfos(numPatch, Nkv, Jpqr, Nijk, Ukv, &
            &   weight, nb_elem_patch)
        call extractNurbsPatchMechInfos(NumPatch, IEN, PROPS, JPROPS,  &
            &   nnode, nb_elem_patch, ELT_TYPE,TENSOR)
        if (elt_type_patch .eq. 'U30' .or. elt_type_patch .eq. 'U10') then
            i = int(props_patch(2))
            call extractMappingInfos(i, nb_elem_patch, Nkv, Jpqr, Nijk, Ukv, &
                &   weight, IEN, PROPS, JPROPS, nnode, ELT_TYPE, TENSOR)  
        endif
        
        if (elt_type_patch .eq. 'U1' .or. & 
            &   elt_type_patch .eq. 'U2' .or. &
            &   elt_type_patch .eq. 'U3' .or. &
            &   elt_type_patch .eq. 'U10' .or. &
            &   elt_type_patch .eq. 'U30') then
            !! Loop on elements
            do numElem = 1, nb_elem_patch(numPatch)
                jelem = jelem + 1
                if (activeElement(jelem) .eq. 1) then
                    do i = 1, nnode_patch
                        coords_elem(:,i) = coords3D(:mcrd, ien_patch(i, numelem))
                    enddo 
                    
                    call extractNurbsElementInfos(numElem)
                    !! Compute elementary RHS
                    if (elt_type_patch .eq. 'U1' .or. & 
                        &   elt_type_patch .eq. 'U2' .or. &
                        &   elt_type_patch .eq. 'U3') then
                        call uSvarsRHS(elRhs(:,:nnode_patch), &
                            &   sol(:, ien_patch(:, numElem)), &
                            &   mcrd, nnode_patch, nbInt(numPatch), &
                            &   coords_elem(:,:nnode_patch),  &
                            &   elt_type_patch, tensor_patch, &
                            &   material_properties(:, numPatch))
                    endif
                    if (elt_type_patch .eq. 'U10') then
                        call uSvarsRHSembdedSol(elRhs(:,:nnode_patch), &
                            &   sol(:, ien_patch(:, numElem)), &
                            &   mcrd, nnode_patch, nnode_map, &
                            &   nb_cp, nbInt(numPatch), &
                            &   coords_elem(:,:nnode_patch), coords3D, &
                            &   tensor_patch, material_properties(:, numPatch))
                    endif
                    
                    !! Assemble to global RHS vector
                    sctr(:nnode_patch) = ien_patch(:, numElem)
                    do icp = 1, nnode_patch
                        rhs(:,sctr(icp)) = rhs(:,sctr(icp)) + elRhs(:,icp)
                    enddo
                    
                endif
            enddo
        endif
    enddo

end subroutine compute_svars_solid_rhs


!! Compute elementary contribution to RHS vector for least square projection 
!! of variables defined at integration points
subroutine uSvarsRHS(rhs, sol, mcrd, nnode, nbInt, coords, elt_type, &
                &   tensor, material_properties)
    use parameters
    implicit none
    
    !! Input arguments
    !! ---------------
    integer, intent(in) :: mcrd
    integer, intent(in) :: nnode
    integer, intent(in) :: nbInt
    double precision, intent(in) :: sol
    dimension sol(mcrd, nnode)
    double precision, intent(in) :: coords
    dimension coords(mcrd, nnode)
    character(len=*), intent(in) :: elt_type
    character(len=*), intent(in) :: tensor
    double precision, intent(in) :: material_properties
    dimension material_properties(2)
    
    !! Output variables
    !! ----------------
    double precision, intent(out) :: rhs
    dimension rhs(2*2*mcrd, nnode)
    
    !! Local variables
    !! ---------------
    double precision :: GaussPdsCoords      !! weight and coordinates of Gauss points
    dimension GaussPdsCoords(mcrd+1, nbInt)
    double precision :: R
    dimension R(nnode)
    double precision :: dRdx
    dimension dRdx(mcrd, nnode)
    double precision :: DetJac
    double precision :: dvol
    double precision :: stran
    dimension stran(2*mcrd, nbInt)
    double precision :: stress
    dimension stress(2*mcrd, nbInt)
    double precision :: ddsdde
    dimension ddsdde(2*mcrd, 2*mcrd)
    integer :: nbPtInt
    integer ntens
    integer :: icp, iPtInt, i, j, icomp
    
    rhs = zero
    stran(:,:) = zero
    ntens = 2*mcrd
    
    !! Number of integration points in each direction
    if ((mcrd .eq. 2) .or. (elt_type .eq. 'U2') .or. (elt_type .eq. 'U3')) then
        nbPtInt = int(nbInt**(1.0/2.0))
        if (nbPtInt**2. < nbInt) nbPtInt = nbPtInt + 1
    elseif (mcrd .eq. 3) then
        nbPtInt = int(nbInt**(1.0/3.0))
        if (nbPtInt**3. < nbint) nbPtInt = nbPtInt +1
    endif
    
    if (elt_type .eq. 'U1') then
        call Gauss(nbPtInt, mcrd, GaussPdsCoords, 0)
        
        !! Loop on integration points
        do iPtInt = 1, nbint
            dvol = zero
            call shap(dRdx, R, DetJac, coords, GaussPdsCoords(2:,iPtInt),mcrd)
            !! Loop on control points to compute strain
            do icp = 1, nnode
                stran(1,iPtInt) = stran(1,iPtInt) + dRdx(1,icp)*sol(1,icp)
                stran(2,iPtInt) = stran(2,iPtInt) + dRdx(2,icp)*sol(2,icp)
                stran(4,iPtInt) = stran(4,iPtInt) &
                        & + dRdx(2,icp)*sol(1,icp)  + dRdx(1,icp)*sol(2,icp)
                if (mcrd .eq. 3) then
                    stran(3,iPtInt) = stran(3,iPtInt) + dRdx(3,icp)*sol(3,icp)
                    stran(5,iPtInt) = stran(5,iPtInt) &
                        & + dRdx(3,icp)*sol(1,icp)  + dRdx(1,icp)*sol(3,icp)
                    stran(6,iPtInt) = stran(6,iPtInt) &
                        & + dRdx(3,icp)*sol(2,icp)  + dRdx(2,icp)*sol(3,icp)
                endif
            enddo
            !! Compute stress
            !! NOTE : in elasticity, call to material_lib is not necessary at 
            !! each integration point but is written here to prepare for non linear cases
            call material_lib(material_properties ,tensor, mcrd, ddsdde)
            call MulVect(ddsdde, stran(:,iPtInt), stress(:,iPtInt), ntens, ntens)
            
            !! Compute contribution to rhs
            do icp = 1, nnode
                dvol = R(icp)*GaussPdsCoords(1, iPtInt)*DetJac
                !! storage 2D : e11 e22 e33 2*e12 s11 s22 s33 s12
                !! storage 3D : e11 e22 e33 2*e12 2*e13 2*e23 s11 s22 s33 s12 s13 s23
                do icomp = 1, ntens
                    rhs(icomp, icp) = rhs(icomp, icp) + stran(icomp, iPtInt)*dvol
                    rhs(ntens+icomp, icp) = rhs(ntens+icomp, icp) + stress(icomp, iPtInt)*dvol
                enddo
            enddo
        enddo
           
    else
        !! NOT IMPLEMENTED YET FOR OTHER ELEMENTS
        write(*,*) 'This function is only implemented for elements U1'
        call exit(666) !! TODO : exit codes should be documented
    endif


end subroutine uSvarsRHS


!! Compute elementary contribution to RHS vector for least square projection 
!! of variables defined at integration points - Embedded solid
subroutine uSvarsRHSembdedSol(rhs, sol, mcrd, nnode, nnode_hull, nb_cp_all, &
                &   nbInt, coords, coords_all, &
                &   tensor, material_properties)
    use parameters
    use embeddedMapping
    
    implicit none
    
    !! Input arguments
    !! ---------------
    double precision, intent(in) :: sol  ! Computed solution
    dimension sol(mcrd, nnode)
    integer, intent(in) :: mcrd  ! Dimension
    integer, intent(in) :: nnode  ! Number of nodes (patch)
    integer, intent(in) :: nnode_hull  ! Number of nodes (hull)
    integer, intent(in) :: nb_cp_all  ! Total number of control points (all patches)
    integer, intent(in) :: nbInt  ! Integration order
    double precision, intent(in) :: coords  ! CP coordinates (patch)
    dimension coords(mcrd, nnode)
    double precision, intent(in) :: coords_all  ! CP coordinates (all patches)
    dimension coords_all(mcrd, nb_cp_all)
    character(len=*), intent(in) :: tensor  ! Type of stresses
    double precision, intent(in) :: material_properties  ! Material matrix
    dimension material_properties(2)
    
    !! Output variables
    !! ----------------
    double precision, intent(out) :: rhs
    dimension rhs(2*2*mcrd, nnode)
    
    !! Local variables
    !! ---------------
    double precision :: GaussPdsCoords  ! Weight and coordinates of Gauss points
    dimension GaussPdsCoords(mcrd+1, nbInt)
    double precision :: DetJac
    double precision :: dvol
    double precision :: stran
    dimension stran(2*mcrd, nbInt)
    double precision :: stress
    dimension stress(2*mcrd, nbInt)
    double precision :: ddsdde
    dimension ddsdde(2*mcrd, 2*mcrd)
    integer :: nbPtInt
    integer ntens
    integer :: icp, iPtInt, i, j, icomp, isave
    double precision :: coef
    
    !! Embedded Volume
    ! - NURBS basis functions
    double precision :: R, dRdTheta, dRdxi
    dimension R(nnode), dRdTheta(nnode, 3), dRdxi(3, nnode)
    ! - Parametric space
    double precision :: Theta
    dimension Theta(3)
    ! - Physical space
    double precision :: xi
    dimension xi(3)
    ! - From parametric to physical space
    double precision :: dxidTheta, dThetadxi, detdxidTheta
    dimension dxidTheta(3, 3), dThetadxi(3, 3)
    ! - From parent element to pametric space
    double precision :: dtildexidTheta, dThetadtildexi, detdThetadtildexi
    dimension dtildexidTheta(3, 3), dThetadtildexi(3, 3)
      
    !! Mapping
    ! - NURBS basis functions
    double precision :: N, dNdxi
    dimension N(nnode_hull), dNdxi(nnode_hull, 3)
    ! - From parametric to physical space
    double precision :: dXdxi, dxidX, detdXdxi, dThetadX
    dimension dXdxi(3, 3), dxidX(3, 3), dThetadX(3, 3)
    ! - From parent element to pametric space
    double precision :: dxidtildexi
    dimension dxidtildexi(3, 3)
    ! - Element infos
    double precision :: coords_hull
    dimension coords_hull(mcrd, nnode_hull)
    integer :: sctr_hull
    dimension sctr_hull(nnode_hull)

    !! Composition Mapping+Volume
    double precision :: dRdX, dRdXT
    dimension dRdX(3, nnode), dRdXT(nnode, 3)
    
    
    !! Initialization 
    !! --------------
    rhs = zero
    stran(:,:) = zero
    ntens = 2*mcrd
    
    !! Number of integration points in each direction
    nbPtInt = int(nbInt**(1.0/30.))
    if (nbPtInt**3. < nbint) nbPtInt = nbPtInt +1
    
    !! Defining Gauss points coordinates and Gauss weights
    call Gauss(nbPtInt, mcrd, GaussPdsCoords, 0)
    
    isave = 0
    
    !! Loop on integration points
    do iPtInt = 1, nbint
        
        !! 1. Embedded solid 
        !! ..................
          
        !! - Compute parametric coordinates from parent element
        Theta(:) = zero
        do j = 1, 3
            coef = GaussPdsCoords(1+j, iPtInt)
            Theta(j) = ((Ukv_elem(2, j) - Ukv_elem(1, j))*coef &
                &   + (Ukv_elem(2, j) + Ukv_elem(1, j)))*0.5d0
        enddo
        !! - Compute NURBS basis functions and derivatives of the embedded solid          
        call evalnurbs(Theta, R, dRdTheta)
        !! - Compute embedded solid physical position
        !!   NB: physical space (embedded) = parametric space (hull)
        xi(:) = zero
        do icp = 1, nnode
            xi(:) =  xi(:) + R(icp)*coords(:, icp)
        enddo
        !! - Gradient of mapping: parent element >> parameter space (embedded)
        dThetadtildexi(:, :) = zero
        do j = 1, 3
            dThetadtildexi(j, j) = 0.5d0*(Ukv_elem(2, j) - Ukv_elem(1, j))
        enddo
        !! Mapping determinant
        call MatrixDet(dThetadtildexi, detdThetadtildexi, 3)
        !! Gradient of mapping: parameter space (embedded) >> physical space (embedded)
        !!   NB: physical space (embedded) = parametric space (hull)
        dxidTheta(:, :) = zero
        do icp = 1, nnode
            dxidTheta(:, 1) = dxidTheta(:, 1) + dRdTheta(icp, 1)*coords(:, icp)
            dxidTheta(:, 2) = dxidTheta(:, 2) + dRdTheta(icp, 2)*coords(:, icp)
            dxidTheta(:, 3) = dxidTheta(:, 3) + dRdTheta(icp, 3)*coords(:, icp)
        enddo
        !! Invert mapping
        call MatrixInv(dThetadxi, dxidTheta, detdxidTheta, 3)
          
        !! 2. Hull object
        !! ..............
          
        !! - Get active element number
        call updateMapElementNumber(xi(:))
        !! - Evaluate NURBS basis functions and derivatives of the hull object
        call evalnurbs_mapping(xi(:), N(:), dNdxi(:, :))
        !! - Extract coordinates of the CPs of the hull object
        if (isave /= current_map_elem) then
            sctr_hull(:) = IEN_map(:, current_map_elem)
            do icp = 1, nnode_hull
                coords_hull(:, icp) = coords_all(:, sctr_hull(icp))
            enddo
            isave = current_map_elem
        endif
        !! - Gradient of mapping: parameter space (hull) >> physical space (hull)              
        dXdxi(:, :) = zero
        do icp = 1, nnode_hull
            dXdxi(:, 1) = dXdxi(:, 1) + dNdxi(icp, 1)*coords_hull(:, icp)
            dXdxi(:, 2) = dXdxi(:, 2) + dNdxi(icp, 2)*coords_hull(:, icp)
            dXdxi(:, 3) = dXdxi(:, 3) + dNdxi(icp, 3)*coords_hull(:, icp)
        enddo
        !! Invert mapping
        call MatrixInv(dxidX, dXdxi, detdXdxi, 3)
          
        !! 3. Composition: hull object x embedded solid
        !! ............................................
        
        !! - Intermediate mapping determinant product
        call MulMat(dThetadxi, dxidX, dThetadX, 3, 3, 3)
        !! - Basis functions composition
        call Mulmat(dRdTheta, dThetadX, dRdXT, nnode, 3, 3)
        !! - Transpose basis functions array for further usage
        dRdX(:, :) = zero
        do icp = 1, nnode
            dRdX(1, icp) = dRdXT(icp, 1)
            dRdX(2, icp) = dRdXT(icp, 2)
            dRdX(3, icp) = dRdXT(icp, 3)
        enddo
        !! - Compute product of all mappings gradients
        DetJac = detdXdxi*detdxidTheta*detdThetadtildexi
        
        !! 4. Compute strain
        !! .................
        do icp = 1, nnode
            stran(1, iPtInt) = stran(1, iPtInt) + dRdX(1, icp)*sol(1, icp)
            stran(2, iPtInt) = stran(2, iPtInt) + dRdX(2, icp)*sol(2, icp)
            stran(3, iPtInt) = stran(3, iPtInt) + dRdX(3, icp )*sol(3, icp)
            stran(4, iPtInt) = stran(4, iPtInt) &
                &   + dRdX(2, icp)*sol(1, icp)  + dRdX(1, icp)*sol(2, icp)
            stran(5, iPtInt) = stran(5, iPtInt) &
                &   + dRdX(3, icp)*sol(1, icp)  + dRdX(1, icp)*sol(3, icp)
            stran(6, iPtInt) = stran(6, iPtInt) &
                &   + dRdX(3, icp)*sol(2, icp)  + dRdX(2, icp)*sol(3, icp)
        enddo
        
        !! 5. Compute stress
        !! .................
        !! NOTE : in elasticity, call to material_lib is not necessary at 
        !! each integration point but is written here to prepare for non linear cases
        call material_lib(material_properties, tensor, mcrd, ddsdde)
        call MulVect(ddsdde, stran(:,iPtInt), stress(:,iPtInt), ntens, ntens)
        
        !! 6. Compute contribution to RHS
        !! ..............................
        
        !! Solid gradient
        dvol = GaussPdsCoords(1, iPtInt)*DetJac
        !! Compute rhs contribution
        do icp = 1, nnode
            !! storage 2D : e11 e22 e33 2*e12 s11 s22 s33 s12
            !! storage 3D : e11 e22 e33 2*e12 2*e13 2*e23 s11 s22 s33 s12 s13 s23
            do icomp = 1, ntens
                rhs(icomp, icp) = rhs(icomp, icp) + stran(icomp, iPtInt)*dvol*R(icp)
                rhs(ntens+icomp, icp) = rhs(ntens+icomp, icp) + stress(icomp, iPtInt)*dvol*R(icp)
            enddo
        enddo
    enddo

end subroutine uSvarsRHSembdedSol
