!! Copyright 2023 Arnaud Duval

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
!! with Yeti. if not, see <https://www.gnu.org/licenses/>

!! Generate VTU file using Bezier cells

subroutine generate_VTU_bezier(filename, i_patch,       &
        &   sol, coords3d,              &
        &   ien, nb_elem_patch,                     &
        &   Nkv, Ukv, Nijk, weight, weight_by_cp, Jpqr,               &
        &   elt_type, tensor, props, jprops,                           &
        &   ind_cp_patch, nb_cp_patch,                    &
        &   nnode, nb_patch, nb_elem, nb_cp, mcrd)

    use parameters
    use nurbspatch

    implicit none

    !! Input arguments
    !! ---------------

    !! Output infos
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i_patch

    !! NURBS geometry
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: coords3d
    dimension coords3d(3, nb_cp)

    double precision, intent(in) :: Ukv, weight, weight_by_cp
    integer, intent(in) :: Nkv, Jpqr, Nijk, nb_cp_patch
    dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem),    &
        &       Ukv(:), weight(:), weight_by_cp(nb_cp_patch)

    !! Patches and elements
    character(len=*), intent(in) :: tensor, elt_type
    double precision, intent(in) :: props
    integer, intent(in) :: mcrd, nnode, nb_patch, nb_elem, ien, nb_elem_patch,       &
        &   jprops, ind_cp_patch
    dimension ien(:), nb_elem_patch(nb_patch), props(:), jprops(nb_patch),      &
        &   nnode(nb_patch), ind_cp_patch(nb_cp_patch)

    !! Analysis solution
    double precision, intent(in) :: sol
    dimension sol(mcrd, nb_cp)

    !! Local variables
    !! ---------------
    double precision :: coords_elem, sol_elem
    dimension :: coords_elem(mcrd, maxval(nnode))
    dimension :: sol_elem(mcrd, maxval(nnode))

    integer :: i, i_elem, i_cp, offset
    integer, allocatable :: conn_vtk(:)

    !! Printing
    write(*,*) 'Postprocessing VTU Bezier ...'
    write(*,*) 'bla bla ...'

    !! File
    open(90, file='results/'//filename//'.vtu', form='formatted')

    !! Write header
    write(90, *) '<VTKFile type="UnstructuredGrid" version="2.2"  >'
    write(90, *) '<UnstructuredGrid>'


    !! Extract patch data
    call extractNurbsPatchGeoInfos(i_patch, Nkv, Jpqr, Nijk, Ukv,       &
        &       weight, nb_elem_patch)
    call extractNurbsPatchMechInfos(i_patch, ien, props, jprops,        &
        &       nnode, nb_elem_patch, elt_type, tensor)

    if (elt_type_patch == "U1") then
        !! Only take into account solid element
        !! TODO only in 3D

        allocate(conn_vtk(nnode_patch))

        !! Start piece
        !! TODO la valeur nnode_patch est fausse ( c'est le nombre de pts par element)
        write(90, *) '<Piece NumberOfPoints="  ', nb_cp_patch,         &
            &   '"  NumberOfCells=" ', nb_elem_patch, '">'

        !! Write degrees
        write(90, *) '<CellData HigherOrderDegrees="HigherOrderDegrees">'
        write(90, *) '<DataArray  type="Int32" Name="HigherOrderDegrees" NumberOfComponents="3" format="ascii">'
        do i_elem = 1, nb_elem_patch(i_patch)
            write(90, *) (Jpqr_patch(i), i=1, 3)
        enddo
        write(90, *) '</DataArray>'
        write(90, *) '</CellData>'


        !! Write control points
        write(90, *) '<Points>'
        write(90, *) '<DataArray  type="Float64" NumberOfComponents="3"  format="ascii" >'

        do i_cp = 1, nb_cp_patch
            write(90, *) coords3d(:, ind_cp_patch(i_cp))
        enddo

        write(90,*) '</DataArray>'
        write(90,*) '</Points>'


        !! Write cells
        write(90, *) '<Cells>'

        write(90, *) '<DataArray  type="Int32"  Name="connectivity"  format="ascii">'
        do i_elem = 1, nb_elem_patch(i_patch)
            !! WARNING il faudra repasser à une numérotation locale au patch !!!!
            call ComputeBezierVTUConnectivity(conn_vtk, IEN_patch(:, i_elem), Jpqr_patch(1), Jpqr_patch(2), Jpqr_patch(3))

            write(90, *) (conn_vtk(i_cp) - 1, i_cp = 1, nnode_patch)
        enddo

        write(90, *) '</DataArray>'

        write(90, *) '<DataArray  type="Int32"  Name="offsets"  format="ascii">'
        offset = 0
        do i_elem = 1, nb_elem_patch(i_patch)
            offset = offset + nnode_patch
            write(90, *) offset
        enddo
        write(90, *) '</DataArray>'

        write(90, *) '<DataArray  type="UInt8"  Name="types"  format="ascii">'
        do i_elem = 1, nb_elem_patch(i_patch)
            write(90, *) '79'     !! 79 == type for Bezier hexahedron
        enddo
        write(90, *) '</DataArray>'

        write(90, *) '</Cells>'

        ! do i_elem = 1, nb_elem_patch(i_patch)
        !     !! extract element solution
        !     do i = 1, nnode_patch
        !         coords_elem(:, i) = coords3d(:mcrd, ien_patch(i, i_elem))
        !         sol_elem(:, i) = sol(:mcrd, ien_patch(i, i_elem))
        !     enddo
        !     call extractNurbsElementInfos(i_elem)
        ! enddo

        !! Write data at control points
        write(90, *) '<PointData RationalWeights="RationalWeights">'

        !! Write weights at control points
        write(90, *) '<DataArray  type="Float64" Name="RationalWeights" NumberOfComponents="1" format="ascii">'
        do i_cp = 1, nb_cp_patch
            write(90, *) weight_by_cp(i_cp)
        enddo
        write(90, *) '</DataArray>'

        !! Write solution at control points
        write(90,*) '<DataArray type="Float64" Name="disp" NumberOfComponents="3" format="ascii">'

        do i_cp = 1, nb_cp_patch
            write(90, *) sol(:, ind_cp_patch(i_cp))
        enddo

        write(90,*) '</DataArray>'
        write(90,*) '</PointData>'

        !! Finalize piece
        write(90, *) '</Piece>'

        deallocate(conn_vtk)
    endif

    !! Finalize file
    write(90, *) '</UnstructuredGrid>'
    write(90, *) '</VTKFile>'


end subroutine generate_vtu_bezier

subroutine ComputeBezierVTUConnectivity(conn_vtk, conn_yeti, p, q, r)
    !! Compute element connectivity for Bezier cell in VTK format from Yeti connectivity

    implicit none

    !! Inputs
    integer, intent(in) :: conn_yeti, p, q, r
    dimension conn_yeti((p+1)*(q+1)*(r+1))
    !! Output
    integer, intent(out) :: conn_vtk
    dimension conn_vtk((p+1)*(q+1)*(r+1))

    !! Local variables
    integer :: i, j, k, counter, n_cp
    integer :: ConnConvert

    n_cp = (p+1)*(q+1)*(r+1)
    !! TODO don't forget to reverse yeti nodes numbering

    !! Vertices
    !! --------
    !! i = 0, j = 0, k = 0
    conn_vtk(1) = conn_yeti(n_cp + 1 - ConnConvert(0, 0, 0, p, q, r))
    !! i = p, j = 0, k = 0
    conn_vtk(2) = conn_yeti(n_cp + 1 - ConnConvert(p, 0, 0, p, q, r))
    !! i = p, j = q, k = 0
    conn_vtk(3) = conn_yeti(n_cp + 1 - ConnConvert(p, q, 0, p, q, r))
    !! i = 0, j = q, k = 0
    conn_vtk(4) = conn_yeti(n_cp + 1 - ConnConvert(0, q, 0, p, q, r))
    !! i = 0, j = 0, k = r
    conn_vtk(5) = conn_yeti(n_cp + 1 - ConnConvert(0, 0, r, p, q, r))
    !! i = p, j = 0, k = r
    conn_vtk(6) = conn_yeti(n_cp + 1 - ConnConvert(p, 0, r, p, q, r))
    !! i = p, j = q, k = r
    conn_vtk(7) = conn_yeti(n_cp + 1 - ConnConvert(p, q, r, p, q, r))
    !! i = 0, j = q, k = r
    conn_vtk(8) = conn_yeti(n_cp + 1 - ConnConvert(0, q, r, p, q, r))

    !! Edges
    !! -----
    counter = 9
    !! Edge 1: j=0, k=0, i croissant
    do i = 1, p-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, 0, 0, p, q, r))
        counter = counter+1
    enddo
    !! Edge 2: i=p, k=0, j croissant
    do j = 1, q-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(p, j, 0, p, q, r))
        counter = counter+1
    enddo
    !! Edge 3: j=q, k=0, i croissant
    do i = 1, p-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, q, 0, p, q, r))
        counter = counter+1
    enddo
    !! Edge 4: i=0, k=0, j croissant
    do j = 1, q-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(0, j, 0, p, q, r))
        counter = counter+1
    enddo
    !! Edge 5: j=0, k=r, i croissant
    do i = 1, p-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, 0, r, p, q, r))
        counter = counter+1
    enddo
    !! Edge 6: i=p, k=r, j croissant
    do j = 1, q-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(p, j, r, p, q, r))
        counter = counter+1
    enddo
    !! Edge 7: j=q, k=r, i croissant
    do i = 1, p-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, q, r, p, q, r))
        counter = counter+1
    enddo
    !! Edge 8: i=0, k=r, j croissant
    do j = 1, q-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(0, j, r, p, q, r))
        counter = counter+1
    enddo
    !! Edge 9: i=0, j=0, k croissant
    do k = 1, r-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(0, 0, k, p, q, r))
        counter = counter+1
    enddo
    !! Edge 10: i=p, j=0, k croissant
    do k = 1, r-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(p, 0, k, p, q, r))
        counter = counter+1
    enddo
    !! Edge 11: i=p, j=q, k croissant
    do k = 1, r-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(p, q, k, p, q, r))
        counter = counter+1
    enddo
    !! Edge 12: i=0, j=q, k croissant
    do k = 1, r-1
        conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(0, q, k, p, q, r))
        counter = counter+1
    enddo


    !! Faces
    !! -----
    !! Face 1: i=0, j croisssant, puis k croissant
    do k = 1, r-1
        do j = 1, q-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(0, j, k, p, q, r))
            counter = counter+1
        enddo
    enddo
    !! Face 2: i=p, j croissant, puis k croissant
    do k = 1, r-1
        do j = 1, q-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(p, j, k, p, q, r))
            counter = counter+1
        enddo
    enddo
    !! Face 3: j=0, i croissant, puis k croissant
    do k = 1, r-1
        do i = 1, p-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, 0, k, p, q, r))
            counter = counter+1
        enddo
    enddo
    !! Face 4: j=q, i croissant, puis k croissant
    do k = 1, r-1
        do i = 1, p-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, q, k, p, q, r))
            counter = counter+1
        enddo
    enddo
    !! Face 5: k=0, i croissant, puis j croissant
    do j = 1, q-1
        do i = 1, p-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, j, 0, p, q, r))
            counter = counter+1
        enddo
    enddo
    !! Face 6 k=r, i croissant, puis j croissant
    do j = 1, q-1
        do i = 1, p-1
            conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, j, r, p, q, r))
            counter = counter+1
        enddo
    enddo

    !! Volume
    !! ------
    !! i croissant, puis j croissant, puis k croissant
    do k = 1, r-1
        do j = 1, q-1
            do i = 1, p-1
                conn_vtk(counter) = conn_yeti(n_cp + 1 - ConnConvert(i, j, k, p, q, r))
                counter = counter+1
            enddo
        enddo
    enddo

end subroutine ComputeBezierVTUConnectivity

function ConnConvert(i, j, k, p, q, r) result(idx)
    !! Compute index of a control point in Yeti convention form i, j, k indices of CP
    integer, intent(in) :: i, j, k, p, q, r
    integer :: idx

    idx = i + (p+1)*j + (p+1)*(q+1)*k +1
end function
