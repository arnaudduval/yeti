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


!! Compute points on the surfaces involved in U5 coupling


subroutine eval_coupling_gpts(filename, npts_u, npts_v, COORDS3D, IEN, nb_elem_patch, Nkv, Ukv, Nijk, weight, Jpqr, ELT_TYPE, &
    &   PROPS, JPROPS, TENSOR, MCRD, nb_patch, nb_elem, nnode, nb_cp )
    use iso_fortran_env

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------

    character(len=*), intent(in) :: filename
    integer, intent(in) :: npts_u, npts_v

    !! NURBS geometry
    integer, intent(in) :: nb_patch, nb_elem, MCRD, nb_cp
    integer, intent(in) :: nnode, nb_elem_patch
    dimension nnode(nb_patch), nb_elem_patch(nb_patch)
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(MCRD, nb_cp)
    
    double precision, intent(in) :: Ukv, weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem),    &
        &     Ukv(:), weight(:)

    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: PROPS
    integer, intent(in) :: JPROPS, IEN
    dimension PROPS(:), JPROPS(nb_patch), IEN(:)



    !! Local Variables
    !! ---------------
    integer :: iPatch, slaveFace, slavePatch, i_embded_patch, icp, icp_map
    integer :: idim
    integer :: ipt_u, ipt_v
    logical :: isSlaveEmbded
    double precision, dimension(MCRD) :: xi_slave, xi_hull, x_phys_slave
    double precision, allocatable :: R_slave(:), N_slave(:)
    double precision, allocatable :: COORDS_elem(:,:), COORDS_elem_map(:,:)
    double precision u, v
    character(len=8) :: fmt     !! Char format
    character(5) :: char_ipatch, char_iface, char_icpl

    fmt = '(I5.5)'

    allocate(COORDS_elem(MCRD, MAXVAL(NNODE)))
    allocate(COORDS_elem_map(MCRD, MAXVAL(NNODE)))

    do iPatch = 1, nb_patch

        call extractNurbsPatchMechInfos(iPatch, IEN, PROPS, JPROPS, &
                &   NNODE, nb_elem_patch, ELT_TYPE, TENSOR)
        if (ELT_TYPE_patch .eq. 'U5') then
            slavePatch = int(PROPS_patch(4))
            slaveFace = int(PROPS_patch(5))

            !! Open file to write
            write(char_icpl, fmt) iPatch
            write(char_ipatch, fmt) slavePatch
            write(char_iface, fmt) slaveFace

            open(95, file='results/'// filename // '_' // char_icpl // '_' // char_ipatch // '_' // char_iface // '.txt')
            
            isSlaveEmbded = .false.

            call extractNurbsPatchGeoInfos(slavePatch, Nkv, Jpqr, Nijk, Ukv,    &
                    &   weight, nb_elem_patch)
            call extractNurbsPatchMechInfos(slavePatch, IEN, PROPS, JPROPS,     &
                    &   NNODE, nb_elem_patch, ELT_TYPE, TENSOR)

            !! Check if slave patch is embedded
            if (ELT_TYPE_patch .eq. 'U10') then
                isSlaveEmbded = .true.
                i_embded_patch = int(PROPS_patch(2))
                call extractMappingInfos(i_embded_patch, nb_elem_patch, Nkv, Jpqr,  &
                        &   Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
            endif

            !! Loop to compute point
            do ipt_u = 1, npts_u
                u = (one*(ipt_u-1))/(npts_u-1)
                do ipt_v = 1, npts_v
                    v = (one*(ipt_v-1))/(npts_v-1)
                    if(isSlaveEmbded) then
                        call point_on_solid_face(u,v, slaveFace, xi_slave)

                        call updateElementNumber(xi_slave(:))
                        if (allocated(R_slave)) deallocate(R_slave)
                        allocate(R_slave(nnode_patch))
                        call evalnurbs_noder(xi_slave(:), R_slave(:))
                        !! compute xi
                        xi_hull(:) = zero
                        do icp = 1, nnode_patch
                            COORDS_elem(:,icp) = coords3D(:MCRD, IEN_patch(icp, current_elem))
                            do idim = 1, MCRD
                                xi_hull(idim) = xi_hull(idim) + R_slave(icp)*COORDS_elem(idim, icp)
                            enddo
                        enddo

                        !! compute X
                        !! Get active element number
                        call updateMapElementNumber(xi_hull(:))
                        !! Evaluate functions
                        if (allocated(N_slave)) deallocate(N_slave)
                        allocate(N_slave(nnode_map))
                        call evalnurbs_mapping_noder(xi_hull(:), N_slave(:))
                        x_phys_slave(:) = zero
                        do icp_map = 1, nnode_map
                            COORDS_elem_map(:, icp_map) = COORDS3D(:MCRD, IEN_map(icp_map, current_map_elem))
                            do idim = 1, MCRD
                                x_phys_slave(idim) = x_phys_slave(idim) + N_slave(icp_map)*COORDS_elem_map(idim, icp_map)
                            enddo
                        enddo
                    else
                        call point_on_solid_face(u,v, slaveFace, xi_slave)
                        call updateElementNumber(xi_slave(:))
                        if (allocated(R_slave)) deallocate(R_slave)
                        allocate(R_slave(nnode_patch))
                        call evalnurbs_noder(xi_slave(:), R_slave(:))
                        
                        !! compute X phys
                        x_phys_slave(:) = zero
                        do icp = 1, nnode_patch
                            COORDS_elem(:,icp) = coords3D(:MCRD, IEN_patch(icp, current_elem))
                            do idim = 1, MCRD
                                x_phys_slave(idim) = x_phys_slave(idim) + R_slave(icp)*COORDS_elem(idim, icp)
                            enddo
                        enddo
                    endif
                    write(95, *) u, v, x_phys_slave(:MCRD)
                enddo
            enddo

            close(95)
        endif


    enddo



end subroutine eval_coupling_gpts