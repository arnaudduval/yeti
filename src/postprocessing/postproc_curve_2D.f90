!! Copyright 2023 Arnaud Duval
!! Copyright 2023 Ange Gagnaire

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


!! Compute displacement and displacement gradient along a side curve of a given patch
subroutine postproc_curve_2d(sol, n_sample, i_patch, i_face,                  &
    &           x_sample, u_sample, dudx_sample, norm_sample, tan_sample, dudxi_sample, &
    &           ien, props, jprops, nnode, nb_elem_patch, elt_type, tensor,         &
    &           coords,                                                             &
    &           nkv, jpqr, nijk, ukv, weight,                                       &
    &           nb_patch, nb_elem, mcrd, nb_cp)





    use parameters
    use nurbspatch

    implicit none

    ! Input variables
    integer, intent(in) :: nb_patch
    integer, intent(in) :: nb_elem
    integer, intent(in) :: mcrd, nb_cp
    double precision, intent(in) :: sol
    dimension sol(mcrd, nb_cp)
    integer, intent(in) :: n_sample, i_patch, i_face
    double precision :: coords(3, nb_cp)

    integer, intent(in) :: ien
    dimension ien(:)
    double precision, intent(in) :: props
    dimension props(:)
    integer, intent(in) :: jprops
    dimension jprops(:)
    integer, intent(in) :: nnode
    dimension nnode(nb_patch)
    integer, intent(in) :: nb_elem_patch
    dimension nb_elem_patch(nb_patch)
    character(len=*), intent(in) :: elt_type
    character(len=*), intent(in) :: tensor
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

    double precision, intent(out) :: x_sample, u_sample, dudx_sample, norm_sample, tan_sample, dudxi_sample
    dimension x_sample(mcrd, n_sample), u_sample(mcrd, n_sample)
    dimension dudx_sample(mcrd, mcrd, n_sample), norm_sample(mcrd, n_sample), tan_sample(mcrd, n_sample)
    dimension dudxi_sample(mcrd, mcrd, n_sample)

    !! Local variables
    integer :: i_sample, icp, j, i
    double precision :: xi, x, sol_pt
    dimension xi(3), x(3), sol_pt(3)
    double precision, dimension(:), allocatable :: R
    double precision, dimension(:,:), allocatable :: dRdxi
    double precision :: coords_elem, sol_elem
    dimension coords_elem(3, maxval(nnode)), sol_elem(3, maxval(nnode))
    double precision :: dxdxi, dxidx, dudxi, dudx
    dimension dxdxi(3,3), dxidx(3,3), dudxi(3,3), dudx(3,3)
    double precision :: det_dxdxi
    double precision :: tan_vect, norm_vect
    dimension tan_vect(mcrd), norm_vect(mcrd)


    write(*, *) 'Postprocessing curve ', i_face, ' of patch ', i_patch

    !! Extract patch infos
    call extractNurbsPatchGeoInfos(i_patch, nkv, jpqr, nijk, ukv,    &
            &   weight, nb_elem_patch)
    call ExtractNurbsPatchMechInfos(i_patch, ien, props, jprops,   &
    &   nnode, nb_elem_patch, elt_type, tensor)

    allocate(R(nnode_patch), dRdxi(nnode_patch, 3))

    do i_sample = 1, n_sample
        if (i_face .eq. 1) then
            xi(1) = 0.
            xi(2) = (i_sample - 1.)/(n_sample - 1.)
            xi(3) = 0.
        else if (i_face .eq. 2) then
            xi(1) = 1.
            xi(2) = (i_sample - 1.)/(n_sample - 1.)
            xi(3) = 0.
        else if (i_face .eq. 3) then
            xi(1) = (i_sample - 1.)/(n_sample - 1.)
            xi(2) = 0.
            xi(3) = 0.
        else if (i_face .eq. 4) then
            xi(1) = (i_sample - 1.)/(n_sample - 1.)
            xi(2) = 1.
            xi(3) = 0.
        else
            write(*,*) 'Wrong face index'
            ! TODO : Raise an error
        endif

        call updateElementNumber(xi)
        do icp = 1, nnode_patch
            coords_elem(:, icp) = coords(:, IEN_patch(icp, current_elem))
            sol_elem(:, icp) = sol(:, IEN_patch(icp, current_elem))
        enddo

        ! write(*, *) 'xi', xi, 'current_elem', current_elem

        call evalnurbs(xi, R, dRdxi)

        x(:) = 0.0
        sol_pt(:) = 0.0
        dxdxi(:,:) = 0.0
        dudxi(:,:) = 0.0

        do icp = 1, nnode_patch
            x(:) = x(:) + R(icp)*coords_elem(:, icp)
            sol_pt(:) = sol_pt(:) + R(icp)*coords_elem(:, icp)
            do j = 1, 3
                dxdxi(:, j) = dxdxi(:, j) + dRdxi(icp, j)*coords_elem(:, icp)
                dudxi(:, j) = dudxi(:, j) + dRdxi(icp, j)*sol_elem(:, icp)
            enddo
        enddo

        call MatrixInv(dxidx(:mcrd, :mcrd), dxdxi(:mcrd, :mcrd), det_dxdxi, mcrd)

        call MulMat(dudxi(:mcrd,:mcrd), dxidx(:mcrd,:mcrd), dudx(:mcrd,:mcrd), mcrd, mcrd, mcrd)



        !! Compute normal vector
        if ((i_face .eq. 1).or.(i_face .eq. 2)) then
            ! parameter is xi(2)
            tan_vect(:mcrd) = dxdxi(:mcrd, 2) / norm2(dxdxi(:, 2))
        else if ((i_face .eq. 3).or.(i_face .eq. 4)) then
            ! parameter is xi(1)
            tan_vect(:mcrd) = dxdxi(:mcrd, 1) / norm2(dxdxi(:, 1))
        endif
        norm_vect(1) = tan_vect(2)
        norm_vect(2) = -tan_vect(0)

        ! TODO : rewrite code without copy (compute directly in sample values)
        x_sample(:, i_sample) = x(:mcrd)
        u_sample(:, i_sample) = sol_pt(:mcrd)
        dudx_sample(:, :, i_sample) = dudx(:mcrd, :mcrd)
        norm_sample(:, i_sample) = norm_vect(:)
        tan_sample(:, i_sample) = tan_vect(:)
        dudxi_sample(:, :, i_sample) = dudxi(:mcrd, :mcrd)

    enddo








    deallocate(R, dRdxi)

end subroutine postproc_curve_2d