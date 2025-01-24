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

!! Write VTU file with results at coupling interface
!! ID of interface is given by input variable lgrge_patch_number which corresponds to the Lagrange U4 patch to process
!! Results are compted on refined mesh (input variable nb_refinement) of slave U00 patch

subroutine generate_coupling_vtu(filename, output_path,     &
            &   lgrge_patch_number, nb_refinement, sol, coords3D, ien, nb_elem_patch, &
            &   nkv, ukv, nijk, weight, jpqr, elt_type, tensor, props, jprops, nnode,    &
            &   mcrd, nb_patch, nb_elem, nb_cp)

    use parameters
    use nurbspatch
    use embeddedMapping


    implicit none

    !! Input arguments
    !! ---------------
    integer, intent(in) :: nb_patch     !! number of patches
    integer, intent(in) :: nb_elem      !! number of element
    integer, intent(in) :: mcrd         !! number of coordinates
    integer, intent(in) :: nb_cp        !! number of control points

    character(len=*) :: filename            !! name of VTU file to write
    character(len=*) :: output_path         !! path to output directory
    integer, intent(in) :: lgrge_patch_number   !! Index of Lagrange interface to process
    integer, intent(in) :: nb_refinement    !! refinement level on parametric directions
    dimension nb_refinement(2)
    double precision, intent(in) :: sol     !! displacement field solution
    dimension sol(mcrd, nb_cp)
    double precision, intent(in) :: coords3D    !! CP coordinates
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
    dimension jprops(:)
    integer, intent(in) :: nnode
    dimension nnode(nb_patch)

    !! Local variables
    !! ---------------
    integer :: i_patch      !! loop variable for patch index
    integer :: i_domain     !! domain patch index
    integer :: i_lgrge      !! Lagrange patch index
    integer :: i_elem       !! loop variable for element index
    integer :: i_vertice    !! loop variable for vertex index
    integer :: i_side       !! variable for loop on sides of interface
    integer :: is_master    !! 1 if domain is master, 0 otherwise
    integer :: dim_interface    !! dimension of interface$
    integer :: domains      !! indices of coupled domains patches (1 for slave, 2 for master)
    dimension domains(2)
    integer :: interfaces   !! indices of coupled interface patches (1 for slave, 2 for master)
    dimension interfaces(2)
    double precision, dimension(:,:), allocatable :: vertice    !! discretization of isoparametric element
    integer :: nb_elem_vtu      !! number of elements in vtu output
    integer :: nb_xi            !! number of points along 1st parametric directions
    integer :: nb_eta           !! number of points along 2nd parametric direction
    integer :: nb_node          !! number of vtu nodes per patch
    integer :: nb_vertice       !! number of vtu points per NURBS element
    integer :: nvar             !! number of values to evaluate per VTU points
    double precision :: coords_elem !! CP coordinates of current element
    dimension coords_elem(mcrd, maxval(nnode))
    double precision :: u_elem      !! displacement solution for current element
    dimension u_elem(mcrd, maxval(nnode))
    double precision :: R, dRdxi !! shape function and derivative
    dimension R(maxval(nnode))
    dimension dRdxi(maxval(nnode), mcrd)
    double precision :: detJac
    double precision, dimension(:,:,:), allocatable :: xibar    !! xibar parametric coordinates for each evaluation point
    double precision, dimension(:,:,:,:), allocatable :: xi       !! xi parametric coordinates for each evaluation point
    double precision, dimension(:,:,:,:), allocatable :: x        !! x physical coordinates for each evaluation point
    double precision, dimension(:,:,:,:), allocatable :: u        !! displacement field value for each evaluation point
    integer :: i, j, j_elem, i_eta, i_xi
    integer :: sctr
    dimension sctr(maxval(nnode))
    integer :: file         !! File unit
    integer :: comp         !! intermediate variable for connectivity computation
    integer :: offset       !! Offset counter when writing VTU file


    interfaces(:) = 0
    domains(:) = 0
    ! write(*,*) 'Post processing interface patch #', lgrge_patch_number
    !! Loop on patches to find corresponding interface patches
    do i_patch = 1, nb_patch
        call ExtractNurbsPatchMechInfos(i_patch, ien, props, jprops,   &
                &   nnode, nb_elem_patch, elt_type, tensor)
        if (elt_type_patch .eq. 'U00') then
            call extractNurbsPatchGeoInfos(i_patch, Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
            i_domain = int(props_patch(2))
            i_lgrge = int(props_patch(3))
            is_master = int(props_patch(4))
            if (i_lgrge .eq. lgrge_patch_number) then
                !! Current U00 patch corresponds to requested Lagrange patch
                interfaces(is_master+1) = i_patch
                domains(is_master+1) = i_domain
                ! if (is_master .eq. 0) then
                !     write(*,*) '    Slave patch found'
                ! elseif (is_master .eq. 1) then
                !     write(*,*) '    Master patch found'
                ! endif
            endif
        endif
    enddo

    !! Discretize interface (from slave side)
    call ExtractNurbsPatchMechInfos(interfaces(1), ien, props, jprops,   &
                &   nnode, nb_elem_patch, elt_type, tensor)
    dim_interface = dim_patch
    nb_xi = 2**(max(nb_refinement(1)-1, 0)) + 1
    nb_eta = 2**(max(nb_refinement(1)-1, 0)) + 1

    nb_elem_vtu = (nb_xi-1)*(nb_eta-1)
    nb_vertice = nb_xi*nb_eta
    allocate(vertice(3, nb_vertice))
    vertice(:,:) = zero
    call DiscretizeIsoparam(vertice, nb_xi, nb_eta, 1, dim_interface)
    nb_node = nb_elem_patch(interfaces(1))*nb_vertice
    allocate(xibar(dim_interface, nb_vertice, nb_elem_patch(interfaces(1))))
    allocate(xi(2, 3, nb_vertice, nb_elem_patch(interfaces(1))))
    allocate(x(2, 3, nb_vertice, nb_elem_patch(interfaces(1))))
    allocate(u(2, 3, nb_vertice, nb_elem_patch(interfaces(1))))
    xibar(:,:,:) = zero
    xi(:,:,:,:) = zero
    x(:,:,:,:) = zero
    u(:,:,:,:) = zero

    !! Compute xibar values
    do i_elem = 1, nb_elem_patch(interfaces(1))
        do i = 1, nnode_patch
            coords_elem(:,i) = coords3D(:mcrd, ien_patch(i, i_elem))
        enddo
        call ExtractNurbsElementInfos(i_elem)
        do i=1, nb_vertice
            !! Get Xibar
            do j = 1, dim_interface
                xibar(j, i, i_elem) = ((ukv_elem(2,j) - ukv_elem(1,j))*vertice(j,i)    &
                                + (ukv_elem(2,j) + ukv_elem(1,j))) *0.5D0
            enddo
        enddo
    enddo

    !! Get data from each side of interface
    do i_side = 1, 2
       if (interfaces(i_side)>0) then
        call ExtractNurbsPatchMechInfos(interfaces(i_side), ien, props, jprops,   &
                &   nnode, nb_elem_patch, elt_type, tensor)
        if (elt_type_patch .eq. 'U00') then
            call extractNurbsPatchGeoInfos(interfaces(i_side), Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
            i_domain = int(props_patch(2))
            !!i_lgrge = int(props_patch(3))
            is_master = int(props_patch(4))
            !! On balaie sur la discretisation de l'esclave
            do i_elem = 1, nb_elem_patch(interfaces(1))
                do i_vertice=1, nb_vertice
                    !! On recherche l'element auquel correspond xibar (= i_elem si esclave)
                    if (i_side .eq. 1) then
                        j_elem = i_elem
                    else
                        call updateElementNumber(xibar(:, i_vertice, i_elem))
                        j_elem = current_elem
                    endif
                    do i = 1, nnode_patch
                        coords_elem(:,i) = coords3D(:mcrd, ien_patch(i, j_elem))
                    enddo
                    call ExtractNurbsElementInfos(j_elem)

                    !! Evaluate basis function
                    call evalnurbs(xibar(:3, i_vertice, i_elem), R(:nnode_patch), dRdxi(:nnode_patch,:))

                    !! Get Xi
                    do j = 1, nnode_patch
                        xi(i_side,:mcrd, i_vertice, i_elem) = xi(i_side,:mcrd, i_vertice, i_elem)   &
                                                        &   + R(j)*coords_elem(:,j)
                    enddo
                enddo
            enddo

            !! Domain to couple
            call extractNurbsPatchGeoInfos(domains(i_side), nkv, jpqr, nijk, ukv,  &
                    &   weight, nb_elem_patch)
            call extractNurbsPatchMechInfos(domains(i_side), ien, props, jprops,   &
                    &   nnode, nb_elem_patch, elt_type, tensor)
            if (elt_type_patch .eq. 'U30') then
                !! WARNING NOT TESTED !!!
                i = int(props_patch(2))
                call extractMappingInfos(i, nb_elem_patch, nkv, jpqr, nijk, ukv,    &
                        weight, ien, props, jprops, nnode, elt_type, tensor)
            endif

            do i_elem = 1, nb_elem_patch(interfaces(1))
                do i_vertice = 1, nb_vertice
                    call updateElementNumber(xi(i_side,:, i_vertice, i_elem))
                    sctr(:nnode_patch) = ien_patch(:, current_elem)
                    do i = 1, nnode_patch
                        coords_elem(:,i) = coords3D(:mcrd, sctr(i))
                        u_elem(:,i) = sol(:mcrd, sctr(i))
                    enddo
                    !! Evaluate basis function
                    call evalnurbs(xi(i_side,:, i_vertice, i_elem), R(:nnode_patch), dRdxi(:nnode_patch,:))

                    !! Get x and u
                    do j=1, nnode_patch
                        !! TODO : pas vraiment necessaire de calculer x sur les deux faces vu qu'elles sont geometriquement identiques
                        x(i_side,:mcrd, i_vertice, i_elem) = x(i_side,:mcrd, i_vertice, i_elem) +     &
                                                &       R(j)*coords_elem(:,j)
                        u(i_side,:mcrd, i_vertice, i_elem) = u(i_side,:mcrd, i_vertice, i_elem) +     &
                                                &       R(j)*u_elem(:,j)
                    enddo
                enddo
            enddo
         endif
      endif
    enddo


    !! Write data to files
    !! for parameter space
    open(91, file=output_path // '/'//filename//'_param.vtu', form='formatted')
    !! for physical space
    open(92, file=output_path // '/'//filename//'_phys.vtu', form='formatted')

    do file = 91,92
        !! Header
        write(file,*) '<VTKFile type="UnstructuredGrid"  version="0.1"   >'
        write(file,*) '<UnstructuredGrid>'

        !! Start piece
        write(file,*) '<Piece NumberOfPoints="  ', nb_node,     &
                &   '"  NumberOfCells="  ',  &
                &   nb_elem_patch(interfaces(1))*nb_elem_vtu, '">'

        !! Write data points
        write(file, *) '<Points>'
        if (file .eq. 91) then
            write(file, *) '<DataArray  type="Float64"'//   &
                    &   '  NumberOfComponents="3"  format="ascii" >'
        elseif (file .eq. 92) then
            write(file, *) '<DataArray  type="Float64"'//   &
                    &   '  NumberOfComponents="3"  format="ascii" >'
        endif
        do i_elem = 1, nb_elem_patch(interfaces(1))
            do i = 1, nb_vertice
                if (file .eq. 91) then
                    write(file,*) xibar(:2, i, i_elem), '0.0'
                elseif (file .eq. 92) then
                    write(file,*) x(1, :3, i, i_elem)
                endif
            enddo
        enddo
        write(file, *) '</DataArray>'
        write(file, *) '</Points>'

        !! Write cells connectivity
        write(file, *) '<Cells>'
        write(file, *) '<DataArray  type="Int32"  Name="connectivity"'//    &
                    &   '  format="acsii">'
        comp = 0
        do i_elem = 1, nb_elem_patch(interfaces(1))
            do i_eta = 1, nb_eta - 1
                do i_xi = 1, nb_xi - 1
                    comp = (i_elem - 1)*nb_vertice
                    comp = comp + (i_eta-1)*nb_xi + i_xi - 1
                    write(file, *) comp, comp+1, comp+nb_xi+1, comp+nb_xi
                enddo
            enddo
        enddo
        write(file, *) '</DataArray>'

        !! Write cells offsets
        write(file, *) '<DataArray  type="Int32"  Name="offsets"'// &
                    &   '  format="ascii"> '
        offset = 0
        do i = 1, nb_elem_patch(interfaces(1))*nb_elem_vtu
            offset = offset + 2**dim_interface
            write(file, *) offset
        enddo
        write(file, *) '</DataArray>'

        !! Write cells type
        write(file, *) '<DataArray  type="UInt8"  Name="types"'//   &
                &       '  format="ascii">'
        do i = 1, nb_elem_patch(interfaces(1))*nb_elem_vtu
            write(file, *) '9'
        enddo
        write(file, *) '</DataArray>'
        write(file, *) '</Cells>'

        !! Write interest variables
        write(file, *) '<PointData>'
        !! displacement on master and slave
        do i_side = 1, 2
            if (i_side .eq. 1) then
                write(file,*)'<DataArray  type="Float64"'//     &
                    &   ' Name="disp_slave" NumberOfComponents="3" format="ascii">'
            elseif (i_side .eq. 2) then
                write(file,*)'<DataArray  type="Float64"'//     &
                &   ' Name="disp_master" NumberOfComponents="3" format="ascii">'
            endif
            do i_elem = 1, nb_elem_patch(interfaces(1))
                do i = 1, nb_vertice
                    write(file, *) u(i_side, :, i, i_elem)
                enddo
            enddo
            write(file,*) '</DataArray>'
        enddo
        !! Condition (disp master - disp slave)
        write(file,*)'<DataArray  type="Float64"'//     &
        &   ' Name="condition" NumberOfComponents="3" format="ascii">'
        do i_elem = 1, nb_elem_patch(interfaces(1))
            do i = 1, nb_vertice
                write(file, *) u(2, :, i, i_elem) - u(1, :, i, i_elem)
            enddo
        enddo
        write(file,*) '</DataArray>'



        write(file, *) '</PointData>'
        !! End of file
        write(file,*) '</Piece>'
        write(file,*) '</UnstructuredGrid>'
        write(file,*) '</VTKFile>'
    enddo

    close(91)
    close(92)
    deallocate(vertice)
    deallocate(xibar, xi, x, u)
end subroutine generate_coupling_vtu

!! Discretize isoparametric element (coordinates in [-1 1]^2 or [-1 1]^3
subroutine DiscretizeIsoparam(vertice, nb_xi, nb_eta, nb_zeta, mcrd)

    use parameters

    implicit none

    !! Input parameters
    !! ----------------
    integer, intent(in) :: nb_xi, nb_eta, nb_zeta
    integer, intent(in) :: mcrd

    !! Returns
    !! -------
    double precision, intent(out) :: vertice
    dimension vertice(3, nb_xi*nb_eta*nb_zeta)

    !! Local variables
    !! ---------------
    integer :: n
    integer :: i_xi, i_eta, i_zeta

    do i_zeta = 1, nb_zeta
        do i_eta = 1, nb_eta
            do i_xi = 1, nb_xi
                n = (i_zeta-1)*nb_eta*nb_xi + (i_eta-1)*nb_xi + i_xi
                vertice(1,n) = two/dble(nb_xi -1)*dble(i_xi -1) - one
                vertice(2,n) = two/dble(nb_eta-1)*dble(i_eta-1) - one
                if (mcrd .eq. 3) then
                    vertice(3,n) = two/dble(nb_zeta-1)*dble(i_zeta-1) - one
                endif
            enddo
        enddo
    enddo

end subroutine DiscretizeIsoparam
