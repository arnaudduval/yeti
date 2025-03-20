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


!! Build sparse coupling matrix with given integration points position
!! Returns
!!      Cdata : vector containing values of elementary matrices
!!      Crow : vector containing row indices of non zero values
!!      Ccol : vector containing column indices of non zero values

!! Coupling with integration points position computed in master patch space (same parametric space as coupling patch)
!! and then projected on the slave patch

!! Integration is made in the master space

!! WARNING only solid 3D case with coupling between two domains is handled (soon available in 2D !!!)

subroutine cplg_matrixU5(output_path, nb_data, &
    &   COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,   &
    &   PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,ind_dof_free,   &
    &   nb_dof_free,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp,    &
    &   nb_dof_tot, order, Cdata,Crow,Ccol)

    use ISO_FORTRAN_ENV

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------

    !! Output path for debug files (projection)
    character(len=*) :: output_path

    !! NURBS geometry
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3,nb_cp)

    double precision, intent(in) :: Ukv, weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),   &
        &     Ukv(:),weight(:)


    !! Patches and Elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
    integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,IEN,   &
        &   nb_elem_patch, JPROPS
    dimension MATERIAL_PROPERTIES(2,nb_patch),  &
        &   PROPS(:),   &
        &   NNODE(nb_patch),    &
        &   IEN(:), &
        &   nb_elem_patch(nb_patch),    &
        &   JPROPS(nb_patch),   &
        &   NBINT(nb_patch)


    !! Degree Of Freedom
    integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
    dimension ind_dof_free(nb_dof_tot)


    !! Storage INFOS
    integer(kind=8), intent(in) :: nb_data

    !! Integration order
    integer, intent(in) :: order

    !! Output variables
    !! ----------------
    integer,          intent(out) :: Crow,Ccol
    double precision, intent(out) :: Cdata
    dimension Cdata(nb_data),Crow(nb_data),Ccol(nb_data)




    !! Local variables
    !! ---------------
    integer(kind=8) :: count
    integer :: iPatch, igps, inode, idim, ielface, idxGP, idom
    integer :: i, j, k, l, n
    integer :: idof, jdof, icp, jcp
    integer :: ielem
    integer :: masterPatch, masterFace, slavePatch, slaveFace, domPatch
    integer, dimension(:,:), allocatable :: saveIEN
    integer, dimension(:), allocatable :: saveEL, saveEM, saveES
    double precision :: factor
    double precision, dimension(:,:), allocatable :: dRdxi
    double precision, dimension(mcrd, mcrd) :: dxdxi
    double precision :: DetJac

    !!  Extract infos
    integer, allocatable          :: sctr(:), sctr_l(:)
    double precision, allocatable :: COORDS_elem(:,:)

    !! Coupling matrix assembly
    double precision, allocatable :: CMAT(:,:)

    !! Integration points
    integer :: nbPtInt, nb_gps
    double precision, dimension(:,:), allocatable :: GaussPdsCoords
    double precision, dimension(:,:), allocatable :: xi_master, xi_slave, xi_interface
    double precision, dimension(:), allocatable :: weightGP
    double precision, dimension(:,:), allocatable :: x_phys, x_phys_slave
    double precision, dimension(:,:), allocatable :: R_master, R_slave, R_lgrge
    logical :: IsElemOnFace

    !! ------------------------ NEW
    !! Manage embedded entities
    integer :: icp_map, i_embded_patch
    double precision, dimension(:,:), allocatable :: COORDS_elem_map
    double precision, dimension(:,:), allocatable :: xi_master_4embded, xi_slave_4embded, &
        &   xi_slave_2checkproj, x_phys_slave_hull
    double precision, dimension(:,:), allocatable :: R_master_4embded, N_slave
    logical :: IsMasterEmbded, IsSlaveEmbded
    !! ------------------------

    !! Manage infos from projection algorithm
    integer :: info

    !!! CHECK PROJECTION --->
    character(len=8) :: fmt
    character(5) :: char_iPatch
    fmt = '(I5.5)'
    !!! <--- CHECK PROJECTION

    !! Allocations
    allocate(sctr(MAXVAL(NNODE)), sctr_l(MAXVAL(NNODE)))
    allocate(COORDS_elem(MCRD,MAXVAL(NNODE)))
    allocate(CMAT(MCRD, maxval(NNODE)**2))

    CMAT(:,:) = zero
    !! Start assembly
    count=1
    do iPatch = 1, nb_patch

        !!! CHECK PROJECTION --->
        write (char_iPatch, fmt) iPatch  ! Converting integer to string using an 'internal file'
        !!! <--- CHECK PROJECTION

        call extractNurbsPatchMechInfos(iPatch, IEN, PROPS, JPROPS, &
                &   NNODE, nb_elem_patch, ELT_TYPE, TENSOR)

        if(ELT_TYPE_patch .eq. 'U5') then
            call extractNurbsPatchGeoInfos(iPatch, Nkv, Jpqr, Nijk, Ukv,    &
                    &       weight, nb_elem_patch)

            masterPatch = int(PROPS_patch(2))
            masterFace = int(PROPS_patch(3))
            slavePatch = int(PROPS_patch(4))
            slaveFace = int(PROPS_patch(5))

            IsMasterEmbded = .false.
            IsSlaveEmbded = .false.

            !! Print info.
            ! write(*,'(A)') '--------------------'
            ! write(*,'(A, I2, A)') 'Patch ', iPatch, ' of type U5'
            ! write(*,'(A, I1, A, I2)') '  > coupling of face no.', masterFace, ' of patch ', &
            !     &   masterPatch
            ! write(*,'(A, I1, A, I2)') '  > on face no.', slaveFace, ' of patch ', slavePatch


            !! Compute integration points position
            !! -----------------------------------

            !! Integration order
            nbPtInt = max(maxval(Jpqr(:,masterPatch)),maxval(Jpqr(:,slavePatch))) + &
                &   maxval(Jpqr(:,iPatch))
            if (nbPtInt .lt. order) then
                nbPtInt = order
            endif
            !! Integration points on master patch
            call extractNurbsPatchGeoInfos(masterPatch, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(masterPatch,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            !! Allocate array for shape function derivatives
            if (allocated(dRdxi)) deallocate(dRdxi)
            allocate(dRdxi(nnode_patch, 3))

            !! Check if master patch is embedded
            if (ELT_TYPE_patch .eq. 'U10') then
                IsMasterEmbded = .true.
                i_embded_patch = int(PROPS_patch(2))
                call extractMappingInfos(i_embded_patch, nb_elem_patch, Nkv, Jpqr,   &
                        &   Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
            endif

            select case(masterFace)
                case(1,2)
                    !! v and w direction
                    !! WARNING : THIS MAY NOT WORK PROPERLY WITH REPEATED KNOT INSIDE KNOT VECTOR
                    !!           (CONTINUITY DROP)
                    if (dim_patch .eq. 3) then
                        nb_gps = (nbPtInt**(dim_patch-1)) * (Nkv_patch(2)-2*Jpqr_patch(2)-1) * &
                            &   (Nkv_patch(3)-2*Jpqr_patch(3)-1)
                    elseif (dim_patch .eq. 2) then
                        nb_gps = (nbPtInt**(dim_patch-1)) * (Nkv_patch(2)-2*Jpqr_patch(2)-1)
                    endif
                case(3,4)
                    !! u and w direction
                    !! WARNING : THIS MAY NOT WORK PROPERLY WITH REPEATED KNOT INSIDE KNOT VECTOR
                    !!           (CONTINUITY DROP)
                    if (dim_patch .eq. 3) then
                        nb_gps = (nbPtInt**(dim_patch-1)) * (Nkv_patch(1)-2*Jpqr_patch(1)-1) * &
                            &   (Nkv_patch(3)-2*Jpqr_patch(3)-1)
                    elseif (dim_patch .eq. 2) then
                        nb_gps = (nbPtInt**(dim_patch-1)) * (Nkv_patch(1)-2*Jpqr_patch(1)-1)
                    endif
                case(5,6)
                    !! u and v direction
                    !! WARNING : THIS MAY NOT WORK PROPERLY WITH REPEATED KNOT INSIDE KNOT VECTOR
                    !!           (CONTINUITY DROP)
                    if (dim_patch .lt. 3) then
                        write(*,*) "Can not define face ", masterFace, " in dimension ", dim_patch
                        call exit(123)
                    endif
                    nb_gps = (nbPtInt**(dim_patch-1)) * (Nkv_patch(1)-2*Jpqr_patch(1)-1) * &
                        &   (Nkv_patch(2)-2*Jpqr_patch(2)-1)
            end select

            if (allocated(GaussPdsCoords)) deallocate(GaussPdsCoords)
            if (allocated(weightGP)) deallocate(weightGP, xi_master, xi_slave, &
                &   xi_interface, saveEM, saveES)

            ! TODO : maybe there are too much value allocated (use nbPtInt**(dim_patch-1) instead ?) => to verify
            allocate(GaussPdsCoords(dim_patch+1, nbPtInt**(dim_patch)))
            !! TODO : use dimension dim_patch_interface instead of 3 to store Gauss points informations
            allocate(weightGP(nb_gps), xi_master(dim_patch, nb_gps), xi_slave(dim_patch, nb_gps), &
                &   xi_interface(3, nb_gps))
            !! TODO : allocated xi_interface(dim_patch - 1, nb_gps)
            allocate(saveEM(nb_gps), saveES(nb_gps))

            call Gauss(nbPtInt, dim_patch, GaussPdsCoords, masterFace)

            ielface = 1
            do ielem = 1, nb_elem_patch(masterPatch)
                call extractNurbsElementInfos(ielem)
                if(IsElemOnFace(masterFace, Nijk_patch(:,ielem), Jpqr_patch, Nkv_patch, dim_patch)) then
                    do igps = 1, nbPtInt**(dim_patch-1)
                        idxGP = (ielface-1)*(nbPtInt**(dim_patch-1))+igps
                        weightGP(idxGP) = GaussPdsCoords(1, igps)
                        do idim = 1, dim_patch
                            xi_master(idim,idxGP) =    &
                                & ((Ukv_elem(2, idim) - Ukv_elem(1,idim)) * &
                                &   GaussPdsCoords(idim+1,igps) + &
                                &   (Ukv_elem(2, idim) + Ukv_elem(1,idim))) * 0.5D0

                            !! 2D Jacobian for surface integration
                            !!    Add contribution only if we are not on the parametric direction
                            !!    corresponding to the interface
                            !!    (e.g., surfaces 3, 4, 5, 6 if masterFace == 1 or 2)
                            if ((2*idim .ne. masterFace) .and. (2*idim-1 .ne. masterFace)) then
                                weightGP(idxGP) = weightGP(idxGP) * &
                                    &   (Ukv_elem(2, idim) - Ukv_elem(1,idim)) * 0.5D0
                            endif
                        enddo
                        ! write(*,*) idxGP, xi_master(:, idxGP)
                    enddo
                    ielface = ielface + 1
                endif
            enddo


            !! Compute coordinates in physical space
            !! -------------------------------------

            !! Data allocation
            !! - Physical coordinates
            if (allocated(x_phys)) deallocate(x_phys)
            allocate(x_phys(MCRD, nb_gps))
            x_phys(:,:) = zero
            !! - Basis functions - master
            if (allocated(R_master)) deallocate(R_master)
            allocate(R_master(nnode_patch, nb_gps))

            if (IsMasterEmbded) then
                !! - Intermediate parametric coordinates
                if (allocated(xi_master_4embded)) deallocate(xi_master_4embded)
                allocate(xi_master_4embded(MCRD, nb_gps))
                xi_master_4embded(:,:) = zero
                !! - Intermediate basis functions
                if (allocated(R_master_4embded)) deallocate(R_master_4embded)
                allocate(R_master_4embded(nnode_map, nb_gps))
                !! - Coordsmap
                if (allocated(COORDS_elem_map)) deallocate(COORDS_elem_map)
                allocate(COORDS_elem_map(MCRD,nnode_map))
            endif

            !!! CHECK PROJECTION --->
            !! Results
            ! open(11, file=output_path // '/verif_proj_patch'// trim(char_iPatch) //'.txt', form='formatted')
            ! write(11, *) '# Physical coordinates - master side'
            !! Warnings
            ! open(12, file=output_path // '/warnings_proj_patch'// trim(char_iPatch) //'.txt', form='formatted')
            ! write(12, *) 'Patch ' // trim(char_iPatch)
            !!! <--- CHECK PROJECTION

            !! Computation
            do igps = 1, nb_gps
                call updateElementNumber(xi_master(:,igps))
                saveEM(igps) = current_elem
                !! Differenciate cases if master is embedded or not
                if (IsMasterEmbded) then
                    call evalnurbs_noder(xi_master(:, igps), R_master(:,igps))
                    do icp = 1, nnode_patch
                        COORDS_elem(:,icp) = COORDS3D(:MCRD,IEN_patch(icp,current_elem))
                        do idim = 1, MCRD
                            xi_master_4embded(idim, igps) = xi_master_4embded(idim, igps) +  &
                                &   R_master(icp,igps)*COORDS_elem(idim,icp)
                        enddo
                    enddo
                else
                    call evalnurbs_noder(xi_master(:, igps), R_master(:,igps))
                    do icp = 1, nnode_patch
                        COORDS_elem(:,icp) = COORDS3D(:MCRD,IEN_patch(icp,current_elem))
                        do idim = 1, MCRD
                            x_phys(idim, igps) = x_phys(idim, igps) + &
                                &   R_master(icp,igps)*COORDS_elem(idim,icp)
                        enddo
                    enddo

                    !! TEMPORARY FIX
                    !! DetJac is needed
                    !! derivatives of function must be computed
                    call evalnurbs(xi_master(:, igps), R_master(:, igps), dRdxi(:,:))
                    dxdxi = 0.0
                    do icp = 1, nnode_patch
                        COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
                        do idim = 1, MCRD
                            dxdxi(idim, :) = dxdxi(idim, :) + &
                                &  dRdxi(icp, :)*COORDS_elem(idim, icp)
                        enddo
                    enddo
                    !! Compute surface Jacobian depending on considered surface
                    select case (masterFace + mcrd*10)
                        case(21,22)
                            DetJac = sqrt(dxdxi(1,2)**2.0 + dxdxi(2,2)**2.0)
                        case(23,24)
                            DetJac = sqrt(dxdxi(1,1)**2.0 + dxdxi(2,1)**2.0)
                        case(31,32)               !Face 1 et 2
                            call SurfElem(dxdxi(:,2),dxdxi(:,3),DetJac)
                        case(33,34)               !Face 3 et 4
                            call SurfElem(dxdxi(:,3),dxdxi(:,1),DetJac)
                        case(35,36)               !Face 5 et 6
                            call SurfElem(dxdxi(:,1),dxdxi(:,2),DetJac)
                    end select
                    weightGP(igps) = weightGP(igps)*DetJac

                    !! END OF TEMPORARY FIX

                endif
                !! Manage embedded entities
                if (IsMasterEmbded) then
                    !! Get active element number
                    call updateMapElementNumber(xi_master_4embded(:, igps))
                    !! Evaluate functions and derivatives
                    call evalnurbs_mapping_noder(xi_master_4embded(:, igps), &
                        &   R_master_4embded(:, igps))
                    !! Extract coordinates of the CPs of the mapping & compute phys. coords.
                    do icp_map = 1, nnode_map
                        COORDS_elem_map(:, icp_map) = COORDS3D(:MCRD, &
                            &   IEN_map(icp_map, current_map_elem))
                        do idim = 1, MCRD
                            x_phys(idim, igps) = x_phys(idim, igps) +   &
                                &   R_master_4embded(icp_map, igps)*COORDS_elem_map(idim, icp_map)
                        enddo
                    enddo
                endif
                !!! CHECK PROJECTION --->
                ! write(11, *) x_phys(:, igps)
                !!! <--- CHECK PROJECTION
            enddo


            !! Projection on slave patch
            !! -------------------------
            call extractNurbsPatchGeoInfos(slavePatch, Nkv,Jpqr,Nijk,Ukv,    &
                &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(slavePatch,IEN,PROPS,JPROPS, &
                &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            !! Data allocation
            !! - Basis functions - slave
            if (allocated(R_slave)) deallocate(R_slave)
            allocate(R_slave(nnode_patch, nb_gps))

            !!! CHECK PROJECTION --->
            !! - Physical coordinates - slave
            if (allocated(x_phys_slave)) deallocate(x_phys_slave)
            allocate(x_phys_slave(MCRD, nb_gps))
            x_phys_slave(:,:) = zero
            !! <--- CHECK PROJECTION

            ! Check if slave patch is embedded
            if (ELT_TYPE_patch .eq. 'U10') then
                IsSlaveEmbded = .true.
                i_embded_patch = int(PROPS_patch(2))
                call extractMappingInfos(i_embded_patch, nb_elem_patch, Nkv, Jpqr,   &
                        &   Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
                !! Allocate vars.
                !! - Xi slave
                if (allocated(xi_slave_4embded)) deallocate(xi_slave_4embded)
                allocate(xi_slave_4embded(3, nb_gps))
                xi_slave_4embded(:,:) = zero
                !! - Coordsmap
                if (allocated(COORDS_elem_map)) deallocate(COORDS_elem_map)
                allocate(COORDS_elem_map(MCRD,nnode_map))
                !!! CHECK PROJECTION --->
                if (allocated(xi_slave_2checkproj)) deallocate(xi_slave_2checkproj)
                allocate(xi_slave_2checkproj(3, nb_gps))
                xi_slave_2checkproj(:,:) = zero
                if (allocated(N_slave)) deallocate(N_slave)
                allocate(N_slave(nnode_map, nb_gps))
                if (allocated(x_phys_slave_hull)) deallocate(x_phys_slave_hull)
                allocate(x_phys_slave_hull(3, nb_gps))
                x_phys_slave_hull(:,:) = zero
                !!! <--- CHECK PROJECTION
            endif

            !! Embedded slave case
            if (IsSlaveEmbded) then
                !! Project points on hull
                do igps = 1, nb_gps
                    call point_inversion_surface(x_phys(:, igps), slaveFace, COORDS3D, &
                        &   nb_cp, .true., xi_slave_4embded(:,igps), info)
                    !! Projection info. message
                    if (info .ne. 0) then
                        ! write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                        !     &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
                    endif
                enddo
                !!! CHECK PROJECTION --->
                ! write(11, *) '# Physical coordinates - slave side (hull)'
                N_slave(:, :) = zero
                do igps = 1, nb_gps
                    !! Get active element number
                    call updateMapElementNumber(xi_slave_4embded(:, igps))
                    !! Evaluate functions and derivatives
                    call evalnurbs_mapping_noder(xi_slave_4embded(:, igps), N_slave(:, igps))
                    do icp_map = 1, nnode_map
                        !! Extract coordinates of the CPs of the mapping
                        COORDS_elem_map(:, icp_map) = &
                            &   COORDS3D(:MCRD, IEN_map(icp_map, current_map_elem))
                        do idim = 1, MCRD
                            !! Compute phys. coords.
                            x_phys_slave_hull(idim, igps) = x_phys_slave_hull(idim, igps) +   &
                                &   N_slave(icp_map, igps)*COORDS_elem_map(idim, icp_map)
                        enddo
                    enddo
                    ! write(11, *) x_phys_slave_hull(:, igps)
                enddo
                ! write(11, *) '# Physical coordinates - slave side (embedded)'
                !!! <--- CHECK PROJECTION
                do igps = 1, nb_gps  ! Project points on embedded entity
                    call point_inversion_surface(xi_slave_4embded(:, igps), slaveFace, COORDS3D, &
                        &   nb_cp, .false., xi_slave(:,igps), info)
                    !! Projection info. message
                    if (info .ne. 0) then
                        ! write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                        !     &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
                    endif
                    !! New search of current element because it may have changed at last
                    !!     projection iteration (rare)
                    call updateElementNumber(xi_slave(:,igps))
                    saveES(igps) = current_elem
                    call evalnurbs_noder(xi_slave(:, igps), R_slave(:,igps))
                    !!! CHECK PROJECTION --->
                    N_slave(:, :) = zero
                    !! Compute Xi
                    do icp = 1, nnode_patch
                        COORDS_elem(:,icp) = COORDS3D(:MCRD,IEN_patch(icp,current_elem))
                        do idim = 1, MCRD
                            xi_slave_2checkproj(idim, igps) = xi_slave_2checkproj(idim, igps) +  &
                                &   R_slave(icp,igps)*COORDS_elem(idim,icp)
                        enddo
                    enddo
                    !! Compute X
                    !! Get active element number
                    call updateMapElementNumber(xi_slave_2checkproj(:, igps))
                    !! Evaluate functions and derivatives
                    call evalnurbs_mapping_noder(xi_slave_2checkproj(:, igps), N_slave(:, igps))
                    do icp_map = 1, nnode_map
                        !! Extract coordinates of the CPs of the mapping
                        COORDS_elem_map(:, icp_map) = &
                            &   COORDS3D(:MCRD, IEN_map(icp_map, current_map_elem))
                        do idim = 1, MCRD
                            !! Compute phys. coords.
                            x_phys_slave(idim, igps) = x_phys_slave(idim, igps) +   &
                                &   N_slave(icp_map, igps)*COORDS_elem_map(idim, icp_map)
                        enddo
                    enddo
                    ! write(11, *) x_phys_slave(:, igps)
                    !!! <--- CHECK PROJECTION
                enddo
            !! Classical case: project points on surface
            else
                !!! CHECK PROJECTION --->
                ! write(11, *) '# Physical coordinates - slave side'
                !!! <--- CHECK PROJECTION
                do igps = 1, nb_gps
                    if (dim_patch .eq. 3) then
                        call point_inversion_surface(x_phys(:, igps), slaveFace, COORDS3D, nb_cp,   &
                            &   .false., xi_slave(:,igps), info)
                    elseif (dim_patch .eq. 2) then
                        call point_inversion_plane_curve(x_phys(:, igps), slaveFace, COORDS3D,      &
                            &   nb_cp, .false., xi_slave(:, igps), info)
                    else
                        write(*,*) "Dimension ", dim_patch, " not available"
                        call exit(123)
                    endif
                    !! Projection info. message
                    if (info .ne. 0) then
                        ! write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                        !     &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
                    endif
                    !! New search of current element because it may have changed at last
                    !!     projection iteration (rare)
                    call updateElementNumber(xi_slave(:,igps))
                    saveES(igps) = current_elem
                    call evalnurbs_noder(xi_slave(:, igps), R_slave(:,igps))
                    !!! CHECK PROJECTION --->
                    do icp = 1, nnode_patch
                        COORDS_elem(:,icp) = COORDS3D(:MCRD,IEN_patch(icp,current_elem))
                        do idim = 1, MCRD
                            x_phys_slave(idim, igps) = x_phys_slave(idim, igps) + &
                                &   R_slave(icp,igps)*COORDS_elem(idim,icp)
                        enddo
                    enddo
                    ! write(11, *) x_phys_slave(:, igps)
                    !!! <--- CHECK PROJECTION
                enddo
            endif

            !!! CHECK PROJECTION --->
            ! close(11)
            ! close(12)
            !!! <--- CHECK PROJECTION


            !! Define integration points on Lagrange patch
            !! -------------------------------------------
            if (dim_patch .eq. 3) then
                select case(masterFace)
                    case(1,2)
                        xi_interface(1,:) = xi_master(2,:)
                        xi_interface(2,:) = xi_master(3,:)
                    case(3,4)
                        xi_interface(1,:) = xi_master(1,:)
                        xi_interface(2,:) = xi_master(3,:)
                    case(5,6)
                        xi_interface(1,:) = xi_master(1,:)
                        xi_interface(2,:) = xi_master(2,:)
                end select
                xi_interface(3,:) = zero
            elseif (dim_patch .eq. 2) then
                select case(masterFace)
                    case(1,2)
                        xi_interface(1,:) = xi_master(2,:)
                    case(3,4)
                        xi_interface(1,:) = xi_master(1,:)
                end select
                xi_interface(2,:) = zero
                xi_interface(3,:) = zero
            endif

            call extractNurbsPatchMechInfos(iPatch, IEN, PROPS, JPROPS, &
                &   NNODE, nb_elem_patch, ELT_TYPE, TENSOR)
            call extractNurbsPatchGeoInfos(iPatch, Nkv, Jpqr, Nijk, Ukv,    &
                &       weight, nb_elem_patch)
            if (allocated(R_lgrge)) deallocate(R_lgrge)
            if (allocated(saveEL)) deallocate(saveEL)
            if (allocated(saveIEN)) deallocate(saveIEN)
            allocate(R_lgrge(nnode_patch, nb_gps))
            allocate(saveEL(nb_gps))
            allocate(saveIEN(nnode_patch, nb_elem_patch(iPatch)))

            saveIEN(:,:) = IEN_patch(:,:)

            do igps = 1, nb_gps
                call updateElementNumber(xi_interface(:,igps))
                call evalLgrge(xi_interface(:,igps), R_lgrge(:,igps))
                saveEL(igps) = current_elem
            enddo

            !! Fill coupling matrix
            !! --------------------

            do idom = 1, 2      !! 1 : master, 2 : slave
                if (idom .eq. 1) then
                    domPatch = masterPatch
                    factor = 1.D0
                else
                    domPatch = slavePatch
                    factor = -1.D0
                endif
                call extractNurbsPatchGeoInfos(domPatch, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
                call extractNurbsPatchMechInfos(domPatch,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
                do igps = 1, nb_gps
                    sctr_l(:NNODE(iPatch)) = saveIEN(:,saveEL(igps))
                    n = nnode_patch*NNODE(iPatch)
                    if(idom .eq. 1) then
                        sctr(:nnode_patch) = IEN_patch(:,saveEM(igps))
                        call cplingdispU5(R_lgrge(:, igps), R_master(:, igps), weightGP(igps), &
                        & NNODE(iPatch), nnode_patch, MCRD, CMAT(:,:n))
                    else
                        sctr(:nnode_patch) = IEN_patch(:,saveES(igps))
                        call cplingdispU5(R_lgrge(:, igps), R_slave(:, igps), weightGP(igps), &
                        & NNODE(iPatch), nnode_patch, MCRD, CMAT(:,:n))
                    endif
                    i = 0
                    do jcp = 1, NNODE(iPatch)
                        jdof = (sctr_l(jcp)-1)*MCRD
                        do icp = 1, nnode_patch
                            idof = (sctr(icp)-1)*MCRD
                            i = i+1
                            do k = 1, MCRD
                                Cdata(count) = factor*CMAT(k,i)
                                Crow(count) = idof + k - 1
                                Ccol(count) = jdof + k - 1
                                count = count+1
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        endif
    enddo

    !! Deallocations
    deallocate(sctr, sctr_l)
    deallocate(COORDS_elem)
    deallocate(CMAT)

end subroutine cplg_matrixU5

!! Compute coupling terms for coupling element U5
subroutine cplingdispU5(Rl, Rd, detJ, NNODE_l, NNODE_d, MCRD, CMAT)

    use parameters

    implicit none

    !! Input arguments
    !! Interpolation functions of Lagrange multiplier at given point
    double precision, intent(in) :: Rl
    dimension Rl(NNODE_l)
    !! Interpolation function of domain at given point
    double precision, intent(in) :: Rd
    dimension Rd(NNODE_d)
    !! Jacobian and Gauss weight
    double precision, intent(in) :: detJ
    !! Number of nodes per element for Lagrange patch
    integer, intent(in) :: NNODE_l
    !! Number of nodes per element for coupled domain
    integer, intent(in) :: NNODE_d
    !! Number of DOF per node
    integer, intent(in) :: MCRD

    !! Output variables
    double precision, intent(out) :: CMAT
    dimension CMAT(MCRD, NNODE_d*NNODE_l)
    integer :: count, i, j, k

    CMAT(:,:) = zero

    ! Assembling
    !! TODO Use 3 dimension for CMAT for sake of simplicity : CMAT(MCRD, NNODE_d, NNODE_l)
    !! TODO 2 : First direction with size MCRD is not necessary : it contains the same values
    count = 1
    do j = 1, NNODE_l
        do i = 1, NNODE_d
            do k = 1, MCRD
                CMAT(k,count) = Rd(i)*Rl(j)*detJ
            enddo
            count = count+1
        enddo
    enddo
end subroutine cplingdispU5
