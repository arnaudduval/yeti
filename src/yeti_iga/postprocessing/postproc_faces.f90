!! Copyright 2022 Arnaud Duval

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


!! Generate a VTU face with info on faces indices for 3D geometries
!! WARNING : works only for 3D geometries

subroutine generate_faces_vtu(filename, output_path, nb_refinement, COORDS3D, &
    &                         IEN, nb_elem_patch, Nkv, Ukv, Nijk, weight, &
    &                         Jpqr, ELT_TYPE, &
    &                         TENSOR, PROPS, JPROPS, NNODE, &
    &                         nb_patch, nb_elem, nb_cp, MCRD)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input variables
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: output_path
    integer, intent(in) :: nb_refinement
    dimension nb_refinement(3)
    character(len=*), intent(in) :: ELT_TYPE, TENSOR
    integer, intent(in) :: nb_patch, nb_elem, nb_cp, MCRD

    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem)
    double precision, intent(in) :: Ukv, weight
    dimension Ukv(:), weight(:)
    integer, intent(in) :: IEN, nb_elem_patch, JPROPS, NNODE
    dimension IEN(:), nb_elem_patch(nb_patch), JPROPS(nb_patch), NNODE(nb_patch)
    double precision, intent(in) :: COORDS3D, PROPS
    dimension COORDS3D(3, nb_cp), PROPS(:)

    !! Local variables
    integer :: i, i_patch, i_face, i_elem, k1, i_eta, i_xi, n
    integer :: nb_eta, nb_xi
    integer :: nb_elem_vtu, nb_vertices
    double precision, allocatable :: vtu_coords(:,:)
    double precision :: vertice
    dimension vertice(3)
    double precision :: R, dRdx, DetJac
    dimension R(maxval(nnode)), dRdx(maxval(nnode), 3)
    double precision :: pt_u, pt_v
    double precision :: coords_elem
    dimension coords_elem(mcrd, maxval(nnode))
    integer :: comp, offset


    logical :: IsElemOnFace
    integer :: NbElemOnFace

    open(90, file=output_path // '/'// filename //'.vtu', form='formatted')

    !! write header
    write(90,*) '<VTKFile type="UnstructuredGrid" version="0.1">'
    write(90,*) '<UnstructuredGrid>'

    do i_patch = 1, nb_patch
        call extractNurbsPatchGeoInfos(i_patch, Nkv, Jpqr, Nijk, Ukv, weight, nb_elem_patch)
        call extractNurbsPatchMechInfos(i_patch, IEN, PROPS, JPROPS, NNODE, nb_elem_patch, ELT_TYPE, TENSOR)

        if (ELT_TYPE_patch == 'U10') then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i, nb_elem_patch, Nkv, Jpqr, Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
        endif

        if (ELT_TYPE_patch == 'U1' .or. ELT_TYPE_patch == 'U10') then
            do i_face = 1, 6
                select case(i_face)
                    case(1,2)
                        !! direction 1 : v, direction 2 : w
                        nb_xi = 2**max(nb_refinement(2), 0) + 1
                        nb_eta = 2**max(nb_refinement(3), 0) + 1
                    case(3,4)
                        !! direction 1 : u, direction 2 : w
                        nb_xi = 2**max(nb_refinement(1), 0) + 1
                        nb_eta = 2**max(nb_refinement(3), 0) + 1
                    case(5,6)
                        !! direction 1 : u, direction 2 : v
                        nb_xi = 2**max(nb_refinement(1), 0) + 1
                        nb_eta = 2**max(nb_refinement(2), 0) + 1
                end select
                nb_vertices = NbElemOnFace(i_face, Jpqr_patch, Nkv_patch) * nb_xi*nb_eta
                nb_elem_vtu = NbElemOnFace(i_face, Jpqr_patch, Nkv_patch) * (nb_xi-1)*(nb_eta-1)

                if (allocated(vtu_coords)) deallocate(vtu_coords)
                allocate(vtu_coords(3, nb_vertices))
                !! Loop on elements
                n = 0
                do i_elem = 1, nb_elem_patch(i_patch)
                    call extractNurbsElementInfos(i_elem)
                    if (IsElemOnFace(i_face, Nijk_patch(:, i_elem), Jpqr_patch, Nkv_patch)) then
                        do i = 1, nnode_patch
                            coords_elem(:,i) = coords3D(:mcrd, ien_patch(i, i_elem))
                        enddo
                        do i_eta = 1, nb_eta
                            do i_xi = 1, nb_xi
                                n = n+1
                                pt_u = two/dble(nb_xi-1) * dble(i_xi - 1) - one
                                pt_v = two/dble(nb_eta-1) * dble(i_eta - 1) - one
                                select case(i_face)
                                    case(1)
                                        vertice(1) = -one
                                        vertice(2) = pt_u
                                        vertice(3) = pt_v
                                    case(2)
                                        vertice(1) = one
                                        vertice(2) = pt_u
                                        vertice(3) = pt_v
                                    case(3)
                                        vertice(1) = pt_u
                                        vertice(2) = -one
                                        vertice(3) = pt_v
                                    case(4)
                                        vertice(1) = pt_u
                                        vertice(2) = one
                                        vertice(3) = pt_v
                                    case(5)
                                        vertice(1) = pt_u
                                        vertice(2) = pt_v
                                        vertice(3) = -one
                                    case(6)
                                        vertice(1) = pt_u
                                        vertice(2) = pt_v
                                        vertice(3) = one
                                end select
                                !! TODO : derivatives are not needed
                                call shap(dRdx, R, DetJac, coords_elem, vertice(1:),MCRD)
                                vtu_coords(:, n) = zero
                                do k1 = 1, nnode_patch
                                    vtu_coords(:mcrd, n) = vtu_coords(:mcrd, n) + R(k1)*coords_elem(:,k1)
                                enddo
                            enddo
                        enddo
                    endif
                enddo
                write(90,*) '<Piece NumberOfPoints="  ', nb_vertices,'"  NumberOfCells="  ', nb_elem_vtu, '">'

                !! Write points coordinates
                write(90,*) '<Points>'
                write(90,*) '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" >'
                do i = 1, nb_vertices
                    write(90,*) vtu_coords(:, i)
                enddo
                write(90,*) '</DataArray>'
                write(90,*) '</Points>'

                !! Write cells connectivity
                write(90,*) '<Cells>'
                write(90,*) '<DataArray  type="Int32"  Name="connectivity"  format="ascii">'
                comp = 0
                do i_elem = 1, nb_elem_vtu
                    do i_eta = 1, nb_eta-1
                        do i_xi = 1, nb_xi-1
                            comp = (i_elem-1)*(nb_xi)*(nb_eta)
                            comp = comp + (i_eta-1)*nb_xi + i_xi - 1
                            write(90,*) comp, comp+1, comp+nb_xi+1, comp+nb_xi
                        enddo
                    enddo
                enddo
                write(90,*) '</DataArray>'

                !! Write cells offset
                write(90,*)'<DataArray  type="Int32"  Name="offsets"  format="ascii"> '
                offset = 0
                do i = 1, nb_elem_vtu
                    offset = offset + 4
                    write(90,*) offset
                enddo
                write(90,*)'</DataArray>'

                !! Write cells types
                write(90,*)'<DataArray  type="UInt8"  Name="types"  format="ascii"> '
                offset = 0
                do i = 1, nb_elem_vtu
                    write(90,*) '9'
                enddo
                write(90,*)'</DataArray>'

                write(90,*) '</Cells>'

                !! Write faces values
                write(90,*) '<PointData>'
                write(90,*) '<DataArray  type="Int32"  Name="face" format="ascii"> '
                do i = 1, nb_vertices
                    write(90,*) i_face
                enddo
                write(90,*) '</DataArray>'
                write(90,*) '</PointData>'

                write(90,*) '</Piece>'
                if (allocated(vtu_coords)) deallocate(vtu_coords)
            enddo
        endif
    enddo

    write(90,*) '</UnstructuredGrid>'
    write(90,*) '</VTKFile>'
    close(90)
end subroutine generate_faces_vtu