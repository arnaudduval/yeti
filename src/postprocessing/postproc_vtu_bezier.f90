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

        !! Write weights at control points
        write(90, *) '<PointData RationalWeights="RationalWeights">'
        write(90, *) '<DataArray  type="Float64" Name="RationalWeights" NumberOfComponents="1" format="ascii">'
        do i_cp = 1, nb_cp_patch
            write(90, *) weight_by_cp(i_cp)
        enddo
        write(90, *) '</DataArray>'
        write(90, *) '</PointData>'

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
            ! hard coding for degree 2 in all directions
            !! Attention, il faut repasser à une numérotation locale au patch !!!!
            write(90,*)  IEN_patch(27 ,i_elem)-1, IEN_patch(25 ,i_elem)-1, IEN_patch(19 ,i_elem)-1,       &
                &       IEN_patch(21 ,i_elem)-1, IEN_patch(9 ,i_elem)-1, IEN_patch(7 ,i_elem)-1,       &
                &       IEN_patch(1 ,i_elem)-1, IEN_patch(3 ,i_elem)-1, IEN_patch(26 ,i_elem)-1,       &
                &       IEN_patch(22 ,i_elem)-1, IEN_patch(20 ,i_elem)-1, IEN_patch(24 ,i_elem)-1,       &
                &       IEN_patch(8 ,i_elem)-1, IEN_patch(4 ,i_elem)-1, IEN_patch(2 ,i_elem)-1,       &
                &       IEN_patch(6 ,i_elem)-1, IEN_patch(18 ,i_elem)-1, IEN_patch(16 ,i_elem)-1,       &
                &       IEN_patch(10 ,i_elem)-1, IEN_patch(12 ,i_elem)-1, IEN_patch(15 ,i_elem)-1,       &
                &       IEN_patch(13 ,i_elem)-1, IEN_patch(17 ,i_elem)-1, IEN_patch(11 ,i_elem)-1,       &
                &       IEN_patch(23 ,i_elem)-1, IEN_patch(5 ,i_elem)-1, IEN_patch(14 ,i_elem)-1
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

        do i_elem = 1, nb_elem_patch(i_patch)
            !! extract element solution
            do i = 1, nnode_patch
                coords_elem(:, i) = coords3d(:mcrd, ien_patch(i, i_elem))
                sol_elem(:, i) = sol(:mcrd, ien_patch(i, i_elem))
            enddo
            call extractNurbsElementInfos(i_elem)
        enddo

        !! Finalize piece
        write(90, *) '</Piece>'

    endif

    !! Finalize file
    write(90, *) '</UnstructuredGrid>'
    write(90, *) '</VTKFile>'


end subroutine generate_vtu_bezier
