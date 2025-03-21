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
!! The inputs for the solution and geometry are given for a given patch that
!! has already undergone Bezier extraction

subroutine generate_VTU_bezier(filename, output_path, i_patch,       &
        &   sol, coords, weights,    &
        &   ien, Jpqr,                            &
        &   nb_cp, mcrd, nnode, nb_elem, &
        &   dim)


    use parameters
    use nurbspatch

    implicit none

    !! Input arguments
    !! ---------------

    !! Output infos
    character(len=*), intent(in) :: filename    !! name of file to write (without extension)
    character(len=*), intent(in) :: output_path !! path to output directory
    integer, intent(in) :: i_patch              !! index of patch (starts at 1)

    !! NURBS geometry (for given patch, after Bezier extraction)
    integer, intent(in) :: nb_cp, mcrd, nnode,        &
        &                  nb_elem, dim
    double precision, intent(in) :: coords,        &
        &                           weights
    dimension coords(3, nb_cp),        &
        &     weights(nb_cp)
    integer, intent(in) :: ien, Jpqr
    dimension ien(nb_elem, nnode)
    dimension Jpqr(nb_elem, dim)

    !! Bezier extracted analysis solution for given patch
    double precision, intent(in) :: sol
    dimension sol(mcrd, nb_cp)

    !! Local variables
    !! ---------------
    integer :: i, i_elem, i_cp, offset
    integer, allocatable :: conn_vtk(:)

    !! File
    open(90, file=output_path//'/'//filename//'.vtu', form='formatted')

    !! Write header
    write(90, *) '<VTKFile type="UnstructuredGrid" version="2.2"  >'
    write(90, *) '<UnstructuredGrid>'


    allocate(conn_vtk(nnode))

    !! Start piece
    write(90, *) '<Piece NumberOfPoints="  ', nb_cp,         &
        &   '"  NumberOfCells=" ', nb_elem, '">'

    !! Write degrees
    write(90, *) '<CellData HigherOrderDegrees="HigherOrderDegrees">'
    write(90, *) '<DataArray  type="Int32" Name="HigherOrderDegrees" NumberOfComponents="3" format="ascii">'
    do i_elem = 1, nb_elem
        write(90, *) (Jpqr(i_elem, i), i=1, 3)
    enddo
    write(90, *) '</DataArray>'
    write(90, *) '</CellData>'

    !! Write control points coordinates
    write(90, *) '<Points>'
    write(90, *) '<DataArray  type="Float64" NumberOfComponents="3"  format="ascii" >'

    do i_cp = 1, nb_cp
        write(90, *) coords(:, i_cp)
    enddo

    write(90,*) '</DataArray>'
    write(90,*) '</Points>'

    !! Write cells
    write(90, *) '<Cells>'

    write(90, *) '<DataArray  type="Int32"  Name="connectivity"  format="ascii">'
    do i_elem = 1, nb_elem
        call ComputeBezierVTUConnectivity(conn_vtk, ien(i_elem, :),        &
            &                             Jpqr(i_elem, 1),                 &
            &                             Jpqr(i_elem, 2),                 &
            &                             Jpqr(i_elem, 3))

        write(90, *) (conn_vtk(i), i=1, nnode)
    enddo

    write(90, *) '</DataArray>'

    write(90, *) '<DataArray  type="Int32"  Name="offsets"  format="ascii">'
    offset = 0
    do i_elem = 1, nb_elem
        offset = offset + nnode
        write(90, *) offset
    enddo
    write(90, *) '</DataArray>'

    write(90, *) '<DataArray  type="UInt8"  Name="types"  format="ascii">'
    do i_elem = 1, nb_elem
        write(90, *) '79'     !! 79 == type for Bezier hexahedron
    enddo
    write(90, *) '</DataArray>'

    write(90, *) '</Cells>'

    !! Write data at control points
    write(90, *) '<PointData RationalWeights="RationalWeights">'

    !! Write weights at control points
    write(90, *) '<DataArray  type="Float64" Name="RationalWeights" NumberOfComponents="1" format="ascii">'
    do i_cp = 1, nb_cp
        write(90, *) weights(i_cp)
    enddo
    write(90, *) '</DataArray>'

    !! Write solution at control points
    write(90,*) '<DataArray type="Float64" Name="disp" NumberOfComponents="3" format="ascii">'
    do i_cp = 1, nb_cp
        write(90, *) sol(:, i_cp)
    enddo

    write(90,*) '</DataArray>'
    write(90,*) '</PointData>'


    !! Finalize piece
    write(90, *) '</Piece>'

    !! Finalize file
    write(90, *) '</UnstructuredGrid>'
    write(90, *) '</VTKFile>'

    close(90)
    deallocate(conn_vtk)

end subroutine generate_VTU_bezier

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
    !! Compute index of a control point in Yeti convention from i, j, k indices of CP
    integer, intent(in) :: i, j, k, p, q, r
    integer :: idx

    idx = i + (p+1)*j + (p+1)*(q+1)*k +1
end function
