!! Copyright 2021 Arnaud Duval
!! Copyright 2021 Marie Guerder

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

!! Write VTU file with displacement solution at nodes and least square
!! projection of values from integration points

subroutine generate_proj_vtu(filename, output_path, nb_refinement, sol, rhs, coords3D, ien, &
        &   nb_elem_patch, nkv, ukv, nijk, weight, jpqr, elt_type, tensor, props, &
        &   jprops, nnode, nb_patch, nb_elem, nb_cp, mcrd)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------
    integer, intent(in) :: nb_patch
    integer, intent(in) :: nb_elem
    integer, intent(in) :: nb_cp
    integer, intent(in) :: mcrd

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: output_path
    integer, intent(in) :: nb_refinement    !! refinement level on parametric directions
    dimension nb_refinement(3)
    double precision, intent(in) :: sol
    dimension sol(mcrd, nb_cp)
    double precision, intent(in) :: rhs
    dimension rhs(2*2*mcrd, nb_cp)
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
    integer, intent(in) :: nnode
    dimension nnode(nb_patch)


    !! Local variables
    !! ---------------
    double precision, dimension(:,:), allocatable :: vars       !! projected variables

    integer :: nb_xi
    integer :: nb_eta
    integer :: nb_zeta
    integer :: jelem
    integer :: ipatch
    integer :: ielem
    integer :: i
    integer :: ntens        !! Size of tensors
    integer :: nb_elemVTU   !! number of elements of projected FE mesh for 1 IGA element
    integer :: nb_vertice   !! number of vertices of projected FE mesh for 1 IGA element
    integer :: nb_node      !! total number of nodes of projected FE mesh
    integer :: iter
    double precision :: coords_elem
    dimension coords_elem(mcrd, maxval(nnode))
    double precision :: sol_elem
    dimension sol_elem(mcrd, maxval(nnode))
    double precision :: rhs_elem
    dimension rhs_elem(2*2*mcrd, maxval(nnode))
    integer :: nvars        !! number of values to compute at projected FE nodes
    character(len=8) :: fmt !! Char format
    character(5) :: char_ipatch  !! Char patch id.
    fmt = '(I5.5)'


    nb_xi = 2**(max(nb_refinement(1)-1, 0)) + 1
    nb_eta = 2**(max(nb_refinement(2)-1, 0)) + 1
    nb_zeta = 2**(max(nb_refinement(3)-1, 0)) + 1

    !! Write on screen
    ! write(*,*) 'Post processing ...'
    ! write(*,*) ' Map data on FEM mesh ...'
    ! write(*,*) ' Writing .vtu file ...'

    !! Compute output data
    jelem = 0

    do ipatch = 1, nb_patch



        !! Extract model info
        call extractNurbsPatchGeoInfos(ipatch, nkv, jpqr, nijk, ukv,    &
            &   weight, nb_elem_patch)
        call extractNurbsPatchMechInfos(ipatch, ien, props, jprops,   &
            &   nnode, nb_elem_patch, elt_type, tensor)

        if ((elt_type_patch .eq. 'U30') .or. (elt_type_patch .eq. 'U10')) then
            i = int(props_patch(2))
            call extractMappingInfos(i, nb_elem_patch, nkv, jpqr, nijk, ukv,    &
                &   weight, ien, props, jprops, nnode, elt_type, tensor)
        endif

        if ((elt_type_patch .eq. 'U1') .or. (elt_type_patch .eq. 'U2') &
            &   .or. (elt_type_patch .eq. 'U3') .or. (elt_type_patch .eq. 'U30') &
            &   .or. (elt_type_patch .eq. 'U10')) then

            !! Open file to write
            write (char_ipatch, fmt) ipatch  ! Converting int to str using an 'internal file'
            open(91, file=output_path // '/'// filename // '_' // char_ipatch //'.vtu', &
                &   form='formatted')
            !! Header
            write(91,*)'<VTKFile type="UnstructuredGrid"  version="0.1"   >'
            write(91,*)'<UnstructuredGrid>'

            ntens = 2*mcrd

            !! Compute size of projected FE mesh
            nb_elemVTU = 1
            nb_vertice = 1
            do i = 1, dim_patch
                nb_elemVTU = nb_elemVTU * 2**(max(nb_refinement(i) - 1, 0))
                nb_vertice = nb_vertice * (2**(max(nb_refinement(i) - 1, 0)) + 1)
            enddo
            nb_node = nb_elem_patch(ipatch)*nb_vertice
            !! Variables : coords ; disp ; strain ; stress
            nvars = 3+3+ntens+ntens
            allocate(vars(nvars*nb_vertice, nb_elem_patch(ipatch)))
            vars = 0.0D0

            !! Loop on elements
            do ielem = 1, nb_elem_patch(ipatch)
                jelem = jelem + 1

                !! Get element solution
                do i = 1, nnode_patch
                    coords_elem(:, i) = coords3D(:mcrd, ien_patch(i, ielem))
                    sol_elem(:,i) = sol(:mcrd, ien_patch(i, ielem))
                    rhs_elem(:,i) = rhs(:,ien_patch(i, ielem))
                enddo
                call extractNurbsElementInfos(ielem)

                !! Classical solid case
                if (elt_type_patch .eq. 'U1') then
                    call project_data_solid(coords_elem, sol_elem, rhs_elem, &
                        &   vars(:, ielem), nvars, nb_vertice, nb_refinement, &
                        &   mcrd, nnode(ipatch))
                !! Embedded solid case
                elseif (elt_type_patch .eq. 'U10') then
                    call project_data_embded_solid(coords_elem, coords3D, &
                        &   sol_elem, rhs_elem, vars(:, ielem), nvars, nb_vertice, &
                        &   nb_refinement, mcrd, nnode(ipatch), nnode_map, nb_cp)
                endif

            enddo
            call writepatchprojection(91, vars, nb_node, nb_elem_patch(ipatch), dim_patch, &
                    &   mcrd, nvars, nb_elemVTU, &
                &   2**max(nb_refinement(1) - 1, 0) + 1, &
                &   2**max(nb_refinement(2) - 1, 0) + 1, &
                &   2**max(nb_refinement(3) - 1, 0) + 1, &
                &   nb_vertice)

            deallocate(vars)

            !! Write end of file
            write(91,*) '</UnstructuredGrid>'
            write(91,*) '</VTKFile>'
            close(91)

            !! Write on screen
            ! write(*,*) ' --> File ' // output_path // '/' // filename // '_' // char_ipatch // '.vtu has been created.'

        endif
        call finalizeNurbsPatch()
        call deallocateMappingData()

    enddo

end subroutine generate_proj_vtu


subroutine project_data_solid(coords, sol, rhs, vars, nvars, nb_vertice, nb_ref,    &
    &   mcrd, nnode)
    !! Interpolate data at control points for an IGA element to a finite element submesh

    use parameters

    implicit none

    !! Input variables
    !! ---------------
    integer, intent(in) :: nb_ref   !! number of refinements in each direction
    dimension nb_ref(3)
    integer, intent(in) :: nnode
    integer, intent(in) :: mcrd
    integer, intent(in) :: nvars
    integer, intent(in) :: nb_vertice           !! number of projected FE nodes
    double precision, intent(in) :: coords      !! CP coordinates
    dimension coords(mcrd, nnode)
    double precision, intent(in) :: sol         !! displacement solution
    dimension sol(mcrd, nnode)
    double precision, intent(in) :: rhs         !! values at CP
    dimension rhs(2*2*mcrd, nnode)

    !! Output variables
    !! ----------------
    double precision, intent(inout) :: vars
    dimension vars(nvars*nb_vertice)



    !! Local variables
    !! ---------------
    integer :: ndofel
    integer :: nb_xi
    integer :: nb_eta
    integer :: nb_zeta
    integer :: ntens
    integer :: ixi, ieta, izeta, n, k1
    double precision :: vertice
    dimension vertice(mcrd, nb_vertice)
    double precision :: DetJac
    double precision :: R
    dimension R(nnode)
    double precision :: dRdx
    dimension dRdx(mcrd, nnode)
    integer :: offset

    vars = zero
    R = zero
    dRdx = zero

    ndofel = nnode*mcrd
    ntens = 2*mcrd

    nb_xi = 2**max(nb_ref(1) - 1, 0) + 1
    nb_eta = 2**max(nb_ref(2) - 1, 0) + 1
    nb_zeta = 2**max(nb_ref(3) - 1, 0) + 1
    if (mcrd .eq. 2) nb_zeta = 1

    do izeta = 1, nb_zeta
        do ieta = 1, nb_eta
            do ixi = 1, nb_xi
                n = (izeta-1)*nb_eta*nb_xi + (ieta-1)*nb_xi + ixi
                vertice(1,n) = two/dble(nb_xi-1)*dble(ixi-1) - one
                vertice(2,n) = two/dble(nb_eta-1)*dble(ieta-1) - one
                if (mcrd .eq. 3) then
                    vertice(3,n) = two/dble(nb_zeta-1)*dble(izeta-1) - one
                endif
                !!
                !! >>> Modify here
                !!
                !! NOTE : it could improve performance to use a routine that compute
                !!        several evaluations at the same time
                call shap(dRdx, R, DetJac, coords, vertice(1:,n),mcrd)
                !! Interpolate values
                offset = (n-1)*nvars
                do k1 = 1, nnode
                    !! coordinates
                    vars(offset+1:offset+mcrd) = &
                        &   vars(offset+1:offset+mcrd) + R(k1)*coords(:, k1)
                    !! diplacements
                    vars(offset+3+1:offset+3+mcrd) = &
                        &   vars(offset+3+1:offset+3+mcrd) + R(k1)*sol(:, k1)
                    !! strain and stress
                    vars(offset+6+1:offset+6+2*ntens) = &
                        &   vars(offset+6+1:offset+6+2*ntens) + R(k1)*rhs(:, k1)
                enddo
            enddo
        enddo
    enddo

end subroutine project_data_solid


subroutine project_data_embded_solid(coords, coords_all, sol, rhs, vars, nvars, &
    &   nb_vertice, nb_ref, mcrd, nnode, nnode_hull, nb_cp_all)
    !! Interpolate data at control points for an IGA embedded solid element
    !! to a finite element submesh

    use parameters
    use embeddedMapping

    implicit none

    !! Input variables
    !! ---------------
    integer, intent(in) :: nb_ref       !! Number of refinements in each direction
    dimension nb_ref(3)
    integer, intent(in) :: nnode        !! Number of nodes (patch)
    integer, intent(in) :: nnode_hull   !! Number of nodes (hull)
    integer, intent(in) :: nb_cp_all    !! Number of nodes (all patches)
    integer, intent(in) :: mcrd
    integer, intent(in) :: nvars
    integer, intent(in) :: nb_vertice             !! Number of projected FE nodes
    double precision, intent(in) :: coords        !! CP coordinates (patch)
    dimension coords(mcrd, nnode)
    double precision, intent(in) :: coords_all    !! CP coordinates (all patches)
    dimension coords_all(mcrd, nnode_hull)
    double precision, intent(in) :: sol           !! Displacement solution
    dimension sol(mcrd, nnode)
    double precision, intent(in) :: rhs           !! Values at CP
    dimension rhs(2*2*mcrd, nnode)

    !! Output variables
    !! ----------------
    double precision, intent(inout) :: vars
    dimension vars(nvars*nb_vertice)

    !! Local variables
    !! ---------------
    integer :: ndofel
    integer :: nb_xi
    integer :: nb_eta
    integer :: nb_zeta
    integer :: ntens
    integer :: ixi, ieta, izeta, n_p, k1, icp, idim
    double precision :: vertice
    dimension vertice(mcrd, nb_vertice)
    integer :: offset
    double precision :: coef
    !! Embedded solid
    double precision :: R                       !! Basis funcs.
    dimension R(nnode)
    double precision :: Theta                   !! Parametric coords.
    dimension Theta(3)
    double precision :: xi                      !! Physical coords.
    dimension xi(3)
    !! Hull object
    double precision :: N                       !! Basis funcs.
    dimension N(nnode_hull)
    double precision :: coords_hull             !! Coordinates
    dimension coords_hull(mcrd, nnode_hull)
    double precision :: sctr_hull               !! Connectivity
    dimension sctr_hull(nnode_hull)
    integer isave                               !! For loops

    !! Initialization
    !! --------------
    vars = zero
    R = zero
    N = zero
    isave = 0

    ndofel = nnode*mcrd
    ntens = 2*mcrd

    nb_xi = 2**max(nb_ref(1) - 1, 0) + 1
    nb_eta = 2**max(nb_ref(2) - 1, 0) + 1
    nb_zeta = 2**max(nb_ref(3) - 1, 0) + 1
    if (mcrd .eq. 2) nb_zeta = 1

    !! Computation
    !! -----------

    do izeta = 1, nb_zeta
        do ieta = 1, nb_eta
            do ixi = 1, nb_xi
                n_p = (izeta-1)*nb_eta*nb_xi + (ieta-1)*nb_xi + ixi
                vertice(1,n_p) = two/dble(nb_xi-1)*dble(ixi-1) - one
                vertice(2,n_p) = two/dble(nb_eta-1)*dble(ieta-1) - one
                if (mcrd .eq. 3) then
                    vertice(3,n_p) = two/dble(nb_zeta-1)*dble(izeta-1) - one
                endif

                !! 1. Embedded solid
                !! ..................

                !! - Compute parametric coordinates from parent element
                Theta(:) = zero
                do idim = 1, 3
                    coef = vertice(idim, n_p)
                    Theta(idim) = ((Ukv_elem(2, idim) - Ukv_elem(1, idim))*coef &
                        &   + (Ukv_elem(2, idim) + Ukv_elem(1, idim)))*0.5d0
                enddo
                !! - Compute NURBS basis functions and derivatives of the embedded solid
                call evalnurbs_noder(Theta, R)
                !! - Compute embedded solid physical position
                !!   NB: physical space (embedded) = parametric space (hull)
                xi(:) = zero
                do icp = 1, nnode
                    xi(:) =  xi(:) + R(icp)*coords(:, icp)
                enddo

                !! 2. Hull object
                !! ..............

                !! - Get active element number
                call updateMapElementNumber(xi(:))
                !! - Evaluate NURBS basis functions and derivatives of the hull object
                call evalnurbs_mapping_noder(xi(:), N(:))
                !! - Extract coordinates of the CPs of the hull object
                if (isave /= current_map_elem) then
                    sctr_hull(:) = IEN_map(:, current_map_elem)
                    do icp = 1, nnode_hull
                        coords_hull(:, icp) = coords_all(:, sctr_hull(icp))
                    enddo
                    isave = current_map_elem
                endif

                !! 3. Interpolate values
                !! .....................

                !! - Init. offset var.
                offset = (n_p-1)*nvars
                !! - Coordinates >>> need hull funcs. & coords.
                do k1 = 1, nnode_hull
                    vars(offset+1:offset+mcrd) = &
                        &   vars(offset+1:offset+mcrd) + N(k1)*coords_hull(:, k1)
                enddo
                !! - Displacement + stress & strain >>> need embedded solid funcs.
                do k1 = 1, nnode
                    !! Displacements
                    vars(offset+3+1:offset+3+mcrd) = &
                        &   vars(offset+3+1:offset+3+mcrd) + R(k1)*sol(:, k1)
                    !! Strain and stress
                    vars(offset+6+1:offset+6+2*ntens) = &
                        &   vars(offset+6+1:offset+6+2*ntens) + R(k1)*rhs(:, k1)
                enddo
            enddo
        enddo
    enddo

end subroutine project_data_embded_solid


subroutine writepatchprojection(file, vars, nb_node, nb_elem, dime, mcrd, nvars, &
        &   nb_elemVTU, nb_xi, nb_eta, nb_zeta, nb_vertice)

    use parameters

    implicit none

    !! Input parameters
    !! ----------------
    integer, intent(in) :: file         !! file unit to write
    integer, intent(in) :: nb_node
    integer, intent(in) :: nb_elem
    integer, intent(in) :: nb_elemVTU
    integer, intent(in) :: nb_vertice
    integer, intent(in) :: dime
    integer, intent(in) :: mcrd
    integer, intent(in) :: nvars
    double precision, intent(in) :: vars
    dimension vars(nvars*nb_vertice, nb_elem)
    integer, intent(in) :: nb_xi, nb_eta, nb_zeta



    !! Local variables
    !! ---------------
    integer :: ielem
    integer :: offset
    integer :: i
    integer :: comp, compt
    integer :: i_xi, i_eta, i_zeta
    double precision :: stress
    dimension stress(2*mcrd)
    double precision :: strain
    dimension strain(2*mcrd)
    integer :: ntens


    ntens = 2*mcrd


    !! Start piece
    write(file,*) '<Piece  NumberOfPoints="  ', nb_node, &
        &   '"  NumberOfCells=" ', nb_elem*nb_elemVTU, '">'

    !! Write data points
    write(file,*) '<Points>'
    write(file,*) '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" >'
    do ielem = 1, nb_elem
        offset = 1
        do i = 1, nb_vertice
            write(file,*) vars(offset:offset+2, ielem)
            offset = offset + nvars
        enddo
    enddo

    write(file,*) '</DataArray>'
    write(file,*) '</Points>'

    !! Write cell connectivity
    write(file,*) '<Cells>'
    write(file,*) '<DataArray  type="Int32"  Name="connectivity"  format="ascii">'

    comp = 0
    if (dime .eq. 2) then
        do i = 1, nb_elem
            do i_eta = 1, nb_eta-1
                do i_xi = 1, nb_xi-1
                    comp = (i-1)*nb_vertice
                    comp = comp + (i_eta-1)*nb_xi + i_xi - 1
                    write(file,*) comp, comp+1, comp+nb_xi+1, comp+nb_xi
                enddo
            enddo
        enddo
    else if (dime .eq. 3) then
        do i = 1, nb_elem
            do i_zeta = 1, nb_zeta-1
                do i_eta = 1, nb_eta-1
                    do i_xi = 1, nb_xi-1
                        comp = (i-1)*nb_vertice
                        comp = comp + (i_zeta-1)*nb_xi*nb_eta   &
                            &   + (i_eta-1)*nb_xi + i_xi-1
                        compt = comp + nb_xi*nb_eta
                        write(file,*) comp, comp+1, comp+nb_xi+1, comp+nb_xi, &
                            &   compt, compt+1, compt+nb_xi+1, compt+nb_xi
                    enddo
                enddo
            enddo
        enddo
    endif
    write(file,*)'</DataArray>'

    !! Write cell offsets
    write(file,*)'<DataArray  type="Int32"  Name="offsets"  format="ascii">'
    offset = 0
    do i = 1, nb_elem*nb_elemVTU
        offset = offset + 2**dime
        write(file,*) offset
    enddo
    write(file,*)'</DataArray>'

    !! Write cells types
    write(file,*)'<DataArray  type="UInt8"  Name="types"  format="ascii">'
    do i = 1, nb_elem*nb_elemVTU
        if (dime .eq. 2) then
            write(file,*) '9'
        elseif (dime .eq. 3) then
            write(file,*)'12'
        endif
    enddo
    write(file,*) '</DataArray>'
    write(file,*) '</Cells>'

    !! Write variables values
    write(file,*) '<PointData>'

    !! Displacement field
    write(file,*) '<DataArray  type="Float64"  Name="disp"  NumberOfComponents="3"  format="ascii">'
    do ielem = 1, nb_elem
        offset = 3
        do i = 1, nb_vertice
            write(file,*) vars(offset+1:offset+3, ielem)
            offset = offset + nvars
        enddo
    enddo
    write(file,*) '</DataArray>'

    !! Strain
    if (mcrd .eq. 2) then
        write(file,*) '<DataArray type="Float64"  Name="strain"  NumberOfComponents="3"  format="ascii">'
        do ielem = 1, nb_elem
            offset = 6
            do i = 1, nb_vertice
                strain(:ntens) = vars(offset+1:offset+ntens,ielem)
                write(file,*) strain(1), strain(2), strain(4)
                offset = offset + nvars
            enddo
        enddo
    elseif (mcrd .eq. 3) then
    write(file,*) '<DataArray type="Float64"  Name="strain"  NumberOfComponents="6"  format="ascii">'
        do ielem = 1, nb_elem
            offset = 6
            do i = 1, nb_vertice
                strain(:ntens) = vars(offset+1:offset+ntens,ielem)
                write(file,*) strain(1), strain(2), strain(3), &
                    &   strain(4), strain(5), strain(6)
                offset = offset + nvars
            enddo
        enddo
    endif
    write(file,*) '</DataArray>'

    !! stress
    if (mcrd .eq. 2) then
        write(file,*) '<DataArray type="Float64"  Name="stress"  NumberOfComponents="3"  format="ascii">'
        do ielem = 1, nb_elem
            offset = 6 + ntens
            do i = 1, nb_vertice
                stress(:ntens) = vars(offset+1:offset+ntens,ielem)
                write(file,*) stress(1), stress(2), stress(4)
                offset = offset + nvars
            enddo
        enddo
    elseif (mcrd .eq. 3) then
    write(file,*) '<DataArray type="Float64"  Name="stress"  NumberOfComponents="6"  format="ascii">'
        do ielem = 1, nb_elem
            offset = 6 + ntens
            do i = 1, nb_vertice
                stress(:ntens) = vars(offset+1:offset+ntens,ielem)
                write(file,*) stress(1), stress(2), stress(3), &
                    &   stress(4), stress(5), stress(6)
                offset = offset + nvars
            enddo
        enddo
    endif
    write(file,*) '</DataArray>'

    !! End of data
    write(file,*) '</PointData>'

    !! end of piece
    write(file,*) '</Piece>'

end subroutine writepatchprojection

