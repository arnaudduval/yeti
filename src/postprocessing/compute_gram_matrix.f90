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

!! Compute Gram matrix
!! Solid elements only
!! Author : Arnaud Duval


subroutine build_CGrammatrix(Mdata, Mrow, Mcol, activeElement, nb_data, &
        & coords3D, ien, nb_elem_patch, nkv, ukv, nijk, weight, jpqr, &
        & elt_type, props, jprops, tensor, ind_dof_free, nb_dof_free, &
        & mcrd, nbint, nb_patch, nb_elem, nnode, nb_cp, nb_dof_tot)
    use parameters
    use nurbspatch
    use embeddedMapping
    
    implicit none
    
    !! Input arguments
    !! ---------------
    
    !! NURBS
    integer, intent(in) :: nb_patch
    character(len=*), intent(in) :: elt_type
    integer, intent(in) :: ien
    dimension ien(:)
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: coords3D
    dimension coords3d(3, nb_cp)
    integer, intent(in) :: nb_elem_patch
    dimension nb_elem_patch(nb_patch)
    integer, intent(in) :: nnode
    dimension nnode(nb_patch)
    integer, intent(in) :: nb_elem
    integer, intent(in) :: nkv
    dimension nkv(3, nb_patch)
    integer, intent(in) :: jpqr
    dimension jpqr(3, nb_patch)    
    integer, intent(in) :: nijk
    dimension nijk(3, nb_elem)
    double precision, intent(in) :: weight
    dimension weight(:)
    double precision, intent(in) :: ukv
    dimension ukv(:)
    double precision, intent(in) :: props
    dimension props(:)
    integer, intent(in) :: jprops
    dimension jprops(:)
    character(len=*), intent(in) :: tensor
    integer, intent(in) :: mcrd
    integer, intent(in) :: nbInt
    dimension nbInt(nb_patch)

    !! Degrees of freedom
    integer, intent(in) :: nb_dof_free
    integer, intent(in) :: nb_dof_tot
    integer, intent(in) :: ind_dof_free
    dimension ind_dof_free(nb_dof_free)
    
    !! Storage infos
    integer, intent(in) :: nb_data
    integer, intent(in) :: activeElement
    dimension activeElement(nb_elem)
    
    !! Output variables
    !! ----------------
    double precision, intent(out) :: Mdata
    integer, intent(out) :: Mrow, Mcol
    dimension Mdata(nb_data), Mrow(nb_data), Mcol(nb_data)

    
    !! Local variables
    !! ---------------
    integer :: numPatch     !! counter for loop on patches
    integer :: numElem      !! counter for loop on elements
    integer :: count        !! counter to build coord sparse matrix
    integer :: jelem        !! index number of current element
    integer :: nnodeSum
    double precision :: coords_elem     !! CP coordinates of current element
    dimension coords_elem(mcrd, maxval(nnode))
    double precision :: cGramMatrix     !! Gram matrix for current element, already initialized to zero
    dimension cGramMatrix(maxval(nnode)*(maxval(nnode)+1)/2)
    integer :: sctr
    dimension sctr(maxval(nnode))
    integer :: i, cpi, cpj


    !! Initialization
    Mdata = zero
    Mrow = 0
    Mcol = 0
    
    !! Start assemble
    count = 1
    jelem = 0
    !! Loop on patches
    do numPatch = 1, nb_patch
        call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv, &
            &        weight,nb_elem_patch)
        call extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,  &
            &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

        if ((ELT_TYPE_patch .eq. 'U30') .or. (ELT_TYPE_patch .eq. 'U10')) then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
                &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)     
        endif
     
        nnodeSum = nnode_patch*(nnode_patch+1)/2
    
        !! Loop on elements
        do numElem = 1, nb_elem_patch(numPatch)
            jelem = jelem + 1
            
            if (activeElement(jelem) .eq. 1) then
                do i = 1, nnode_patch
                    coords_elem(:,i) = coords3D(:mcrd, ien_patch(i, numElem))
                enddo
                cGramMatrix(:) = zero
                
                call extractNurbsElementinfos(numElem)
                
                !! build elementary matrix
                if (ELT_TYPE_patch .eq. 'U1') then
                    call UGRAMMATRIX_byCP(MCRD,nnode_patch,NBINT(numPatch),    &
            &           COORDS_elem(:,:nnode_patch),ELT_TYPE_patch,        &
            &           PROPS_patch,JPROPS_patch,            &
            &           cGramMatrix(:nnodeSum))      
                else if (ELT_TYPE_patch .eq. 'U10') then
                    call UGRAMMATRIX_U10_byCP(MCRD,nnode_patch,nnode_map,nb_cp,NBINT(numPatch),    &
            &           COORDS_elem(:,:nnode_patch),coords3D,ELT_TYPE_patch,        &
            &           PROPS_patch,JPROPS_patch,            &
            &           cGramMatrix(:nnodeSum))      
                endif
                
                !! Assemble to global matrix
                sctr(:nnode_patch) = ien_patch(:, numElem)
                i = 0
                do cpj = 1, nnode_patch
                    !! cas cpi < cpj
                    do cpi = 1, cpj-1
                        i = i+1
                        Mdata(count) = cGramMatrix(i)
                        Mrow(count) = sctr(cpi) - 1
                        Mcol(count) = sctr(cpj) - 1
                        count = count+1
                    enddo
                    !! cas cpi == cpj
                    i = i+1
                    cGramMatrix(i) = cGramMatrix(i)*0.5D0
                    Mdata(count) = cGramMatrix(i)
                    Mrow(count) = sctr(cpj) - 1
                    Mcol(count) = sctr(cpj) -1
                    count = count + 1
                enddo
            endif
        
        enddo
        call FinalizenurbsPatch()
    
    enddo

end subroutine build_CGrammatrix


!! Compute elementary Gram matrix
!! do not build the matrix, just return values for each couple of control point
!! (triangular sup part only)
subroutine UGRAMMATRIX_byCP(mcrd, nnode, nbInt, coords, elt_type, &
         & PROPS,JPROPS,cGramMatrix)
    
    use parameters
    implicit none
    
    !! Input arguments
    !! ---------------
    integer, intent(in) :: mcrd
    integer, intent(in) :: nnode
    integer, intent(in) :: nbInt
    double precision, intent(in) :: coords
    dimension coords(mcrd, nnode)
    character(len=*), intent(in) :: elt_type
    integer, intent(in) :: jprops
    double precision, intent(in) :: props
    dimension props(jprops)
    
    double precision, intent(out) :: cGramMatrix
    dimension cGramMatrix(nnode*(nnode+1)/2)
    
    
    !! Local variables
    !! ---------------
    double precision :: GaussPdsCoords      !! weight and coordinates of Gauss points
    dimension GaussPdsCoords(mcrd+1, nbInt)
    integer :: nbPtInt
    double precision :: h              !! shells only, not necessary for solids
    double precision :: h3_12          !! shells only, not necessary for solids
    double precision :: detJac
    double precision :: R
    dimension R(nnode)    
    double precision :: dRdx
    dimension dRdx(mcrd, nnode)
    double precision dvol
    
    double precision :: temp1
    integer n, k, count, nodi, nodj
    
    
    
    !! Initialize matrix to zero
    cGramMatrix(:) = zero
    h = zero        !! used for shells
    h3_12 = zero    !! used for shell
    
    !! Number of integration points in each direction
    if ((mcrd .eq. 2) .or. (elt_type .eq. 'U2') .or. (elt_type .eq. 'U3')) then
        nbPtInt = int(nbInt**(1.0/2.0))
        if (nbPtInt**2. < nbInt) nbPtInt = nbPtInt + 1
    elseif (mcrd .eq. 3) then
        nbPtInt = int(nbInt**(1.0/3.0))
        if (nbPtInt**3. < nbint) nbPtInt = nbPtInt +1
    endif
    
    !! Define gauss points coordinates and Gauss weights
    if (elt_type .eq. 'U1' .or. elt_type .eq. 'U10') then
        call Gauss(nbPtInt, mcrd, GaussPdsCoords, 0)
    elseif (elt_type .eq. 'U2' .or. elt_type .eq. 'U3') then
        call Gauss(nbPtInt, 2, GaussPdsCoords(:mcrd,:), 0)
        h = props(3)
        h3_12 = h*h*h/12.D0
    endif
    
    !! Loop on integration points
    do n = 1, nbInt
        !! Compute Nurbs basis functions at integration points
        dvol = zero
        
        if (elt_type .eq. 'U1') then
            !! Solid element
            call shap(dRdx,R,DetJac,coords,GaussPdsCoords(2:,n),mcrd)
            dvol = GaussPdsCoords(1,n)*DetJac
            
            !! Assembling cGramMatrix
            count = 1
            do nodj = 1, nnode
                temp1 = dvol*R(nodj)
                do nodi = 1, nodj
                    cGramMatrix(count) = cGramMatrix(count) + temp1*R(nodi)
                    count = count+1
                enddo
            enddo
            
        elseif ((elt_type .eq. 'U2') .or. (elt_type .eq. 'U3')) then
            !! NOT IMPLEMENTED YET
            !! can be inspired from UMASSMAT routine
            write(*,*) 'This function has not been implemented yet for elements U2 and U3'
            call exit(666) !! TODO : Exit codes should be documented
            
        endif

    enddo
    
end subroutine UGRAMMATRIX_byCP

subroutine UGRAMMATRIX_U10_byCP(mcrd, nnode, nnodemap, nb_cp, nbInt, coords, coordsall, elt_type, &
            & PROPS,JPROPS,cGramMatrix)
    
    use parameters
    use embeddedMapping
    implicit none

    !! Input arguments
    !! ---------------
    integer, intent(in) :: mcrd
    integer, intent(in) :: nnode
    integer, intent(in) :: nnodemap
    integer, intent(in) :: nbInt
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: coords, coordsall
    dimension coords(mcrd, nnode), coordsall(3, nb_cp)
    character(len=*), intent(in) :: elt_type
    integer, intent(in) :: jprops
    double precision, intent(in) :: props
    dimension props(jprops)
    
    double precision, intent(out) :: cGramMatrix
    dimension cGramMatrix(nnode*(nnode+1)/2)


    !! Local variables
    !! ---------------
    double precision :: GaussPdsCoords      !! weight and coordinates of Gauss points
    dimension GaussPdsCoords(mcrd+1, nbInt)
    integer :: nbPtInt
    double precision :: detJac
    double precision :: R, dRdtheta
    dimension R(nnode), dRdtheta(nnode, 3)    
    double precision :: dRdx
    dimension dRdx(mcrd, nnode)
    double precision :: N, dNdxi
    dimension N(nnodemap), dNdxi(nnodemap, 3)
    double precision dvol
    double precision :: theta, xi
    dimension theta(3), xi(3)
    double precision :: dThetadtildexi
    dimension dThetadtildexi(3,3)
    double precision :: dxdxi, dxidx, detdxdxi
    dimension dxdxi(3,3), dxidx(3,3)
    double precision :: dthetadxi, dxidTheta, detdxidTheta
    dimension dthetadxi(3,3), dxidtheta(3,3)
    double precision detdThetadtildexi
    double precision COORDSmap
    dimension COORDSmap(MCRD, NNODEmap)
    integer :: sctr_map
    dimension sctr_map(nnodemap)


    
    double precision :: temp1, coef
    integer :: nintp, j, k, count, nodi, nodj, numcp
    integer :: isave


    !! Initialize matrix to zero
    cGramMatrix(:) = zero
    
    !! Number of integration points in each direction
    nbPtInt = int(nbInt**(1.0/30.))
    if (nbPtInt**3. < nbint) nbPtInt = nbPtInt +1

    !! Define gauss points coordinates and Gauss weights
    call Gauss(nbPtInt, mcrd, GaussPdsCoords, 0)
    isave = 0

    !! Loop on integration points
    do nintp = 1, nbInt
        !! Compute Nurbs basis functions at integration points
        dvol = zero
    
        !! Embedded solid element
        theta(:) = zero
        do j =1, 3
            coef = GaussPdsCoords(j+1, nintp)
            theta(j) = ((Ukv_elem(2, j) - Ukv_elem(1, j))*coef          &
            &            + (Ukv_elem(2, j) + Ukv_elem(1, j)))*0.5d0
        enddo

        !! Compute NURBS basis functions and derivatives of the embedded solid
        call evalnurbs(theta, R, dRdTheta)

        !! Compute embedded solid physical position
        !! (physical space (embedded) = parametric space (hull))
        xi(:) = zero
        do numcp = 1, nnode
            xi(:) = xi(:) + R(numcp)*coords(:, numcp)
        enddo

        !! Gradient of mapping: parent element >> parameter space (embedded)
        dThetadtildexi(:,:) = zero
        do j = 1, dim_patch
            dThetadtildexi(j, j) = 0.5d0*(Ukv_elem(2, j) - Ukv_elem(1, j))
        enddo

        call MatrixDet(dThetadtildexi, detdThetadtildexi, 3)

        !! Gradient of mapping: parameter space (embedded) >> physical space (embedded)
        dxidTheta(:, :) = zero
        do numCP = 1, NNODE
            dxidTheta(:, 1) = dxidTheta(:, 1) + dRdTheta(numCP, 1)*COORDS(:, numCP)
            dxidTheta(:, 2) = dxidTheta(:, 2) + dRdTheta(numCP, 2)*COORDS(:, numCP)
            dxidTheta(:, 3) = dxidTheta(:, 3) + dRdTheta(numCP, 3)*COORDS(:, numCP)
        enddo
        
        call MatrixInv(dThetadxi, dxidTheta, detdxidTheta, 3)

        !! Hull
        !! - Get active element number
        call updateMapElementNumber(xi(:))
      
        !! - Evaluate NURBS basis functions and derivatives of the hull object
        call evalnurbs_mapping(xi(:), N(:), dNdxi(:, :))

        
      
        !! - Extract coordinates of the CPs of the hull object
        if (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:, current_map_elem)
            do numCP = 1, NNODEmap
                COORDSmap(:, numCP) = COORDSall(:, sctr_map(numCP))
            enddo
            isave = current_map_elem
        endif

        !! - Gradient of mapping: parameter space (hull) >> physical space (hull)              
        dXdxi(:, :) = zero
        do numCP = 1, NNODEmap
            dXdxi(:, 1) = dXdxi(:, 1) + dNdxi(numCP, 1)*COORDSmap(:, numCP)
            dXdxi(:, 2) = dXdxi(:, 2) + dNdxi(numCP, 2)*COORDSmap(:, numCP)
            dXdxi(:, 3) = dXdxi(:, 3) + dNdxi(numCP, 3)*COORDSmap(:, numCP)
        enddo
      
        call MatrixInv(dxidX, dXdxi, detdXdxi, 3)

        !! Full Jacobian
        DetJac = detdXdxi * detdxidTheta * detdThetadtildexi
        dvol = GaussPdsCoords(1,nintp)*DetJac

        !! Assembling cGramMatrix
        count = 1
        count = 1
        do nodj = 1, nnode
            temp1 = dvol*R(nodj)
            do nodi = 1, nodj
                cGramMatrix(count) = cGramMatrix(count) + temp1*R(nodi)
                count = count+1
            enddo
        enddo
    enddo
end subroutine UGRAMMATRIX_U10_byCP
