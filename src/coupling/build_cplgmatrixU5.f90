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

subroutine cplg_matrixU5(nb_data, &
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
    integer :: iPatch, igps, inode, idim, ielface, idxGP, idom, iknotint
    integer :: i, j, k, l, n
    integer :: idof, jdof, icp, jcp
    integer :: ielem, nb_elem_slaveonface, nb_knot
    integer :: masterPatch, masterFace, slavePatch, slaveFace, domPatch, derivCoupling
    integer, dimension(:,:), allocatable :: saveIEN
    integer, dimension(:), allocatable :: saveEL, saveEM, saveES, saveESK
    double precision :: factor
    double precision, dimension(:,:,:), allocatable :: dRdxi_master, dRdxi_slave
    double precision, dimension(:, :, :), allocatable :: dxdxi_master, dxdxi_slave
    double precision :: DetJac
    double precision :: u0

    !!  Extract infos
    integer, allocatable          :: sctr(:), sctr_l(:)
    double precision, allocatable :: COORDS_elem(:,:)

    !! Coupling matrix assembly
    double precision, allocatable :: CMAT(:,:)

    !! Integration points
    integer :: nbPtInt, nb_gps
    double precision, dimension(:,:), allocatable :: GaussPdsCoords
    double precision, dimension(:,:), allocatable :: xi_master, xi_slave, xi_interface, xi_slave_knot
    double precision, dimension(:), allocatable :: weightGP
    double precision, dimension(:), allocatable :: dist_elem_gps
    double precision, dimension(:,:), allocatable :: x_phys, x_phys_slave, x_phys_slave_knot
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

    ! write(*,*)'IEN', IEN
    write(*,*)'MCRD', MCRD
    write(*,*)'nb_patch', nb_patch
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
            derivCoupling = int(PROPS_patch(6))

            IsMasterEmbded = .false.
            IsSlaveEmbded = .false.

            !! Print info.
            write(*,'(A)') '--------------------'
            write(*,'(A, I2, A)') 'Patch ', iPatch, ' of type U5'
            write(*,'(A, I1, A, I2)') '  > coupling of face no.', masterFace, ' of patch ', &
                &   masterPatch
            write(*,'(A, I1, A, I2)') '  > on face no.', slaveFace, ' of patch ', slavePatch


            !! Compute slave patch knots position (ONLY FOR SLAVE PATCH DIM = 2 )
            !! ----------------------------------
            

            call extractNurbsPatchGeoInfos(slavePatch, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(slavePatch,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            if (dim_patch .eq. 2) then    
                write(*,*) 'Ukv2 slavepatch : ', Ukv2_patch

                ! where(Nkv_patch(1,:) .eq. maxval(Nkv_patch(1,:))) 
                ! endwhere 
                    
                nb_elem_slaveonface = 0

                do ielem = 1, nb_elem_patch(slavePatch)
                    call extractNurbsElementInfos(ielem)
                    if(IsElemOnFace(slaveFace, Nijk_patch(:,ielem), Jpqr_patch, Nkv_patch, dim_patch)) then
                        nb_elem_slaveonface = 1 + nb_elem_slaveonface
                        write(*,'(A)') '--------------------'
                        write(*,'(A, I2, A)') 'Element on face : ', ielem
                        write(*,*) 'Nijk : ', Nijk_patch(:,ielem)
                        write(*,*) 'Nj : ', Nijk_patch(2,ielem)
                        write(*,*) 'xi knot : ', Ukv2_patch(Nijk_patch(2,ielem))
                        write(*,*) 'xi knot 2 : ', Ukv2_patch(Nijk_patch(2,ielem)+1)
                        write(*,*) 'max : ', maxval(Nijk_patch(1,:))
                        write(*,*) 'nb_elem_slaveonface : ', nb_elem_slaveonface
                    endif
                enddo

                if (allocated(xi_slave_knot)) deallocate(xi_slave_knot)
                allocate(xi_slave_knot(dim_patch, nb_elem_slaveonface*2))

                nb_knot = 0
                
                select case(slaveFace)
                    case(1,2)
                        do ielem = 1, nb_elem_patch(slavePatch)
                            call extractNurbsElementInfos(ielem)
                            if(IsElemOnFace(slaveFace, Nijk_patch(:,ielem), Jpqr_patch, Nkv_patch, dim_patch)) then
                                nb_knot= 1 + nb_knot
                                xi_slave_knot(2,nb_knot) = Ukv2_patch(Nijk_patch(2,ielem))
                                nb_knot = 1 + nb_knot
                                xi_slave_knot(2,nb_knot) = Ukv2_patch(Nijk_patch(2,ielem)+1)
                            endif
                        enddo
                    case(3,4)
                        do ielem = 1, nb_elem_patch(slavePatch)
                            call extractNurbsElementInfos(ielem)
                            if(IsElemOnFace(slaveFace, Nijk_patch(:,ielem), Jpqr_patch, Nkv_patch, dim_patch)) then
                                nb_knot= 1 + nb_knot
                                xi_slave_knot(1,nb_knot) = Ukv1_patch(Nijk_patch(1,ielem))
                                nb_knot = 1 + nb_knot
                                xi_slave_knot(1,nb_knot) = Ukv1_patch(Nijk_patch(1,ielem)+1)
                            endif
                        enddo
                end select

                select case(slaveFace)
                    case(1)
                        xi_slave_knot(1, :) = 0
                    case(2)
                        xi_slave_knot(1, :) = 1
                    case(3)
                        xi_slave_knot(2, :) = 0
                    case(4)
                        xi_slave_knot(2, :) = 1
                end select

                write(*,*) 'xi_slave_knot1 : ', xi_slave_knot(1,:)
                write(*,*) 'xi_slave_knot2 : ', xi_slave_knot(2,:)

                !! Data allocation
                !! - Physical coordinates
                if (allocated(x_phys_slave_knot)) deallocate(x_phys_slave_knot)
                allocate(x_phys_slave_knot(MCRD, nb_knot))
                x_phys_slave_knot(:,:) = zero

                
                !! - Basis functions - slave
                if (allocated(R_slave)) deallocate(R_slave)
                allocate(R_slave(nnode_patch, nb_knot))

                if (allocated(saveESK)) deallocate(saveESK)
                allocate(saveESK(nb_knot))

                if (allocated(dist_elem_gps)) deallocate(dist_elem_gps)
                allocate(dist_elem_gps(nb_elem_slaveonface))

                write(*,*) '----------------------------------------- '

                !! Computation
                do iknotint = 1, nb_knot
                    call updateElementNumber(xi_slave_knot(:,iknotint))
                    saveESK(iknotint) = current_elem
                    call evalnurbs_noder(xi_slave_knot(:, iknotint), R_slave(:,iknotint))
                    do icp = 1, nnode_patch
                        COORDS_elem(:,icp) = COORDS3D(:MCRD,IEN_patch(icp,current_elem))
                        do idim = 1, MCRD
                            x_phys_slave_knot(idim, iknotint) = x_phys_slave_knot(idim, iknotint) + &
                                &   R_slave(icp,iknotint)*COORDS_elem(idim,icp)
                        enddo
                    enddo
                enddo
                write(*,*) 'x_phys_slave_knot : ', x_phys_slave_knot
            endif   

            !! Compute integration points position
            !! -----------------------------------

            !! Integration order
            nbPtInt = max(maxval(Jpqr(:,masterPatch)),maxval(Jpqr(:,slavePatch))) + &
                &   maxval(Jpqr(:,iPatch))
            if (nbPtInt .lt. order) then
                nbPtInt = order
            endif
            write(*,*) "Nb de points d'integration",  nbPtInt
            !! Integration points on master patch
            call extractNurbsPatchGeoInfos(masterPatch, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(masterPatch,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            !! Check if master patch is embedded
            if (ELT_TYPE_patch .eq. 'U10') then
                IsMasterEmbded = .true.
                i_embded_patch = int(PROPS_patch(2))
                call extractMappingInfos(i_embded_patch, nb_elem_patch, Nkv, Jpqr,   &
                        &   Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
            endif
            write(*,*) 'dim_patch', dim_patch
            write(*,*) 'Nkv_patch', Nkv_patch
            write(*,*) 'Jpqr_patch', Jpqr_patch
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
            write(*,*) 'nb_gps' , nb_gps
            write(*,*) (nbPtInt**(dim_patch-1)), (Nkv_patch(2)-2*Jpqr_patch(2)-1)
            write(*,*) 'nkv' , Nkv_patch(1), Nkv_patch(2)
            ! write(*,*) 'ATTENTION VALEUR DE NB_GPS FORCEE'
            ! if (derivCoupling .eq. 0) then
            !     nb_gps = 30
            ! endif

            ! if (derivCoupling .eq. 1) then
            !     nb_gps = 24
            ! endif

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
            write(*,*) "MCRD = ", MCRD
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

            ! Correctif nb_gps temporaire 
            write(*,*) 'idxGP', idxGP
            nb_gps = idxGP

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
                        !! Allocate array for shape function derivatives
            if (allocated(dRdxi_master)) deallocate(dRdxi_master)
            allocate(dRdxi_master(nnode_patch, 3, nb_gps))
            if (allocated(dRdxi_slave)) deallocate(dRdxi_slave)
            allocate(dRdxi_slave(nnode_patch, 3, nb_gps))
            if (allocated(dxdxi_master)) deallocate(dxdxi_master)
            allocate(dxdxi_master(3, 3, nb_gps))
            if (allocated(dxdxi_slave)) deallocate(dxdxi_slave)
            allocate(dxdxi_slave(3, 3, nb_gps))


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
            open(11, file='results/verif_proj_patch'// trim(char_iPatch) //'.txt', form='formatted')
            write(11, *) '# Physical coordinates - master side'
            !! Warnings
            open(12, file='results/warnings_proj_patch'// trim(char_iPatch) //'.txt', form='formatted')
            write(12, *) 'Patch ' // trim(char_iPatch)
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
                        ! write(*,*) 'COORDS_elem', COORDS_elem(:,icp)
                        do idim = 1, MCRD
                            x_phys(idim, igps) = x_phys(idim, igps) + &
                                &   R_master(icp,igps)*COORDS_elem(idim,icp)
                        enddo
                    enddo
                    ! write(*,*) 'COORDS_elem', COORDS_elem(:,icp)
                    
                    !! TEMPORARY FIX
                    !! DetJac is needed
                    !! derivatives of function must be computed
                    call evalnurbs(xi_master(:, igps), R_master(:, igps), dRdxi_master(:,:,igps))
                    dxdxi_master(:,:,igps) = 0.0
                    do icp = 1, nnode_patch
                        COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
                        do idim = 1, MCRD
                            dxdxi_master(idim, :, igps) = dxdxi_master(idim, :, igps) + &
                                &  dRdxi_master(icp, :, igps)*COORDS_elem(idim, icp)
                        enddo
                    enddo
                    !! Compute surface Jacobian depending on considered surface
                    select case (masterFace + mcrd*10)
                        case(21,22)
                            DetJac = sqrt(dxdxi_master(1,2,igps)**2.0 + dxdxi_master(2,2,igps)**2.0)
                        case(23,24)
                            DetJac = sqrt(dxdxi_master(1,1,igps)**2.0 + dxdxi_master(2,1,igps)**2.0)
                        case(31,32)               !Face 1 et 2
                            call SurfElem(dxdxi_master(:,2,igps),dxdxi_master(:,3, igps),DetJac)
                        case(33,34)               !Face 3 et 4
                            call SurfElem(dxdxi_master(:,3,igps),dxdxi_master(:,1,igps),DetJac)
                        case(35,36)               !Face 5 et 6
                            call SurfElem(dxdxi_master(:,1,igps),dxdxi_master(:,2,igps),DetJac)
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
                write(11, *) x_phys(:, igps)
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
                        write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                            &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
                    endif
                enddo
                !!! CHECK PROJECTION --->
                write(11, *) '# Physical coordinates - slave side (hull)'
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
                    write(11, *) x_phys_slave_hull(:, igps)
                enddo
                write(11, *) '# Physical coordinates - slave side (embedded)'
                !!! <--- CHECK PROJECTION
                do igps = 1, nb_gps  ! Project points on embedded entity
                    call point_inversion_surface(xi_slave_4embded(:, igps), slaveFace, COORDS3D, &
                        &   nb_cp, .false., xi_slave(:,igps), info)
                    !! Projection info. message
                    if (info .ne. 0) then
                        write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                            &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
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
                    write(11, *) x_phys_slave(:, igps)
                    !!! <--- CHECK PROJECTION
                enddo

            !! Classical case: project points on surface
            else
                !!! CHECK PROJECTION --->
                write(11, *) '# Physical coordinates - slave side'
                write(*,*) 'nb_gps', nb_gps
                !!! <--- CHECK PROJECTION
                do igps = 1, nb_gps
                    write(*,*) '********************', igps
                    write(*,*) 'x_phys', x_phys(:,igps)
                    ! dist_elem_gps(:) = zero
                    ! do ielem = 1, nb_elem_slaveonface
                    !     write(*,*) 'distn', norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2-1)) + &
                    !         & norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2)) 
                    !     dist_elem_gps(ielem) = norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2-1)) + &
                    !         & norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2)) 
                    ! enddo
                    ! write(*,*) 'dist', minval(dist_elem_gps)

                    ! do ielem = 1, nb_elem_slaveonface
                    !     if (dist_elem_gps(ielem) .eq. minval(dist_elem_gps)) then
                    !         u0 = (xi_slave_knot(2,ielem*2-1) + xi_slave_knot(2,ielem*2))/2
                    !     endif
                    ! enddo
                    ! write(*,*) 'u0', u0
                    if (dim_patch .eq. 3) then
                        call point_inversion_surface(x_phys(:, igps), slaveFace, COORDS3D, nb_cp,   &
                            &   .false., xi_slave(:,igps), info)
                    elseif (dim_patch .eq. 2) then
                        dist_elem_gps(:) = zero
                        do ielem = 1, nb_elem_slaveonface
                            ! write(*,*) 'distn', norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2-1)) + &
                            !     & norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2)) 
                            dist_elem_gps(ielem) = norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2-1)) + &
                                & norm2(x_phys(:,igps) - x_phys_slave_knot(:,ielem*2)) 
                        enddo
                        write(*,*) 'dist', minval(dist_elem_gps)

                        do ielem = 1, nb_elem_slaveonface
                            if (dist_elem_gps(ielem) .eq. minval(dist_elem_gps)) then
                                u0 = (xi_slave_knot(2,ielem*2-1) + xi_slave_knot(2,ielem*2))/2
                            endif
                        enddo
                        write(*,*) 'u0', u0

                        call point_inversion_plane_curve(x_phys(:, igps), slaveFace, COORDS3D,      &
                            &   nb_cp, .false., u0, xi_slave(:, igps), info)
                    else
                        write(*,*) "Dimension ", dim_patch, " not available"
                        call exit(123)
                    endif
                    !! Projection info. message
                    if (info .ne. 0) then
                        write(12,'(A, /, A, I5.1, /, A, I1)') "======== WARNING ========", &
                            &   "Gauss point nb.", igps, "   >> projection exit with info = ", info
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

                    !! TODO : dxdxi pour l'esclave  
                    call evalnurbs(xi_slave(:, igps), R_slave(:, igps), dRdxi_slave(:,:,igps))
                    dxdxi_slave(:,:,igps) = 0.0
                    do icp = 1, nnode_patch
                        COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
                        do idim = 1, MCRD
                            dxdxi_slave(idim, :, igps) = dxdxi_slave(idim, :, igps) + &
                                &  dRdxi_slave(icp, :, igps)*COORDS_elem(idim, icp)
                        enddo
                    enddo

                    write(11, *) x_phys_slave(:, igps)
                    write(*,*) x_phys(:, igps)
                    ! write(*,*) nb_gps
                    write(*,*) xi_slave(:, igps)
                    !!! <--- CHECK PROJECTION
                enddo
            endif

            !!! CHECK PROJECTION --->
            close(11)
            close(12)
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
            ! write(*,*) 'NNODE(iPatch)', NNODE(iPatch)
            !! Fill coupling matrix
            !! --------------------

            do idom = 1, 2     !! 1 : master, 2 : slave !! 
                if (idom .eq. 1) then
                    domPatch = masterPatch
                    factor = 1.D0
                else
                    domPatch = slavePatch
                    factor = -1.D0
                endif
                if (derivCoupling .eq. 1) then
                    factor = 1.D0
                endif
                call extractNurbsPatchGeoInfos(domPatch, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
                call extractNurbsPatchMechInfos(domPatch,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
                do igps = 1, nb_gps
                ! do igps = 1 , 1
                    ! write(*,*) 'Points de Gauss', igps
                    sctr_l(:NNODE(iPatch)) = saveIEN(:,saveEL(igps))
                    ! write(*,*) 'sctr_l', sctr_l !! Allocation ?????????
                    n = nnode_patch*NNODE(iPatch)
                    if(idom .eq. 1) then
                        sctr(:nnode_patch) = IEN_patch(:,saveEM(igps))
                        ! write(*,*) 'sctrM', sctr
                        if (derivCoupling .eq. 0) then
                            call cplingdispU5(R_lgrge(:, igps), R_master(:, igps), weightGP(igps), &
                            & NNODE(iPatch), nnode_patch, MCRD, CMAT(:,:n))
                        elseif (derivCoupling .eq. 1) then 
                            call cplingddispU5(R_lgrge(:, igps), dRdxi_master(:, :, igps), & 
                            & dxdxi_master(:(dim_patch),:(dim_patch), igps), weightGP(igps), &
                            & NNODE(iPatch), nnode_patch, dim_patch, dim_patch, masterFace, CMAT(:,:n))
                        endif
                    else
                        sctr(:nnode_patch) = IEN_patch(:,saveES(igps))
                        ! write(*,*) 'sctrS', sctr
                        if (derivCoupling .eq. 0) then
                            call cplingdispU5(R_lgrge(:, igps), R_slave(:, igps), weightGP(igps), &
                            & NNODE(iPatch), nnode_patch, MCRD, CMAT(:,:n))
                        elseif (derivCoupling .eq. 1) then
                            call cplingddispU5(R_lgrge(:, igps), dRdxi_slave(:, :, igps), & !! dim_patch+1 ??????
                            & dxdxi_slave(:(dim_patch),:(dim_patch), igps), weightGP(igps), &
                            & NNODE(iPatch), nnode_patch, dim_patch, dim_patch, slaveFace, CMAT(:,:n))
                        endif
                    endif

                    ! ! Calcul du gradient sous la forme [dUx/dy dUy/dy dUx/dx dUy/dx]
                    ! i = 0
                    ! do jcp = 1, NNODE(iPatch)
                    !     jdof = (sctr_l(jcp)-1)*MCRD
                    !     do icp = 1, nnode_patch
                    !         idof = (sctr(icp)-1)*MCRD
                    !         i = i+1
                    !         Cdata(count) = factor*CMAT(jcp,i)
                    !         Crow(count) = idof 
                    !         Ccol(count) = jdof 
                    !         count = count+1
                    !         if (CMAT(jcp,i) .ne. 0) then
                    !             write(*,*) 'idof', idof 
                    !             write(*,*) 'jdof', jdof 
                    !             write(*,*) 'icp', icp
                    !             write(*,*) 'jcp', jcp
                    !             write(*,*) 'count', count
                    !             write(*,*) 'CMAT', CMAT(jcp,i)
                    !             write(*,*) '*************************************'
                    !         endif 
                    !         Cdata(count) = factor*CMAT(jcp,i)
                    !         Crow(count) = idof + 1
                    !         Ccol(count) = jdof + 1
                    !         count = count+1
                    !     enddo
                    ! enddo

                    ! couplage
                    i = 0
                    do jcp = 1, NNODE(iPatch)
                        jdof = (sctr_l(jcp)-1)*MCRD ! lagrange
                        do icp = 1, nnode_patch
                            idof = (sctr(icp)-1)*MCRD  ! déplacement
                            i = i+1
                            do k = 1, MCRD
                                Cdata(count) = factor*CMAT(k,i)
                                Crow(count) = idof + k - 1
                                Ccol(count) = jdof + k - 1
                                ! write(*,*) 'count', count
                                ! write(*,*) 'Crow', Crow(count)
                                ! write(*,*) 'Ccol', Ccol(count)
                                ! write(*,*) 'Cdata', Cdata(count)
                                count = count+1
                            enddo
                        enddo
                    enddo

                enddo
            enddo
        endif
    enddo

    ! write(*,*), 'Cdata', Cdata
    ! write(*,*), 'Crow', Crow
    ! write(*,*), 'Ccol', Ccol

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

subroutine cplingddispU5(Rl, dRdxid, dxdxi, detJ, NNODE_l, NNODE_d, MCRD, dim_patchd, i_face, CMAT)

    use parameters

    implicit none

    !! Input arguments 
    !! Interpolation functions of Lagrange multiplier at given point
    double precision, intent(in) :: Rl
    dimension Rl(NNODE_l)
    !! derivative of interpolation function of domain at given point
    double precision, intent(in) :: dRdxid
    dimension dRdxid(NNODE_d, dim_patchd)
    !! derivative of phys coord with recpect to parametric coord
    double precision, intent(in) :: dxdxi
    dimension dxdxi(MCRD, dim_patchd)

    
    !! Jacobian and Gauss weight
    double precision, intent(in) :: detJ
    !! Number of nodes per element for Lagrange patch
    integer, intent(in) :: NNODE_l
    !! Number of nodes per element for coupled domain
    integer, intent(in) :: NNODE_d
    !! Number of DOF per node
    integer, intent(in) :: MCRD
    !! dimension of patch 
    integer, intent(in) :: dim_patchd
    integer, intent(in) :: i_face

    !! Output variables
    double precision, intent(out) :: CMAT
    dimension CMAT(MCRD, NNODE_d*NNODE_l)

    !! Local variables
    integer :: count, i, j, k, l, m
    !! derivative of parametric coord with respect to phys coord
    double precision, dimension(dim_patchd, MCRD) :: dxidx
    !! pas utilisé
    double precision :: det_dxdxi 
    double precision, dimension(MCRD) :: tan_vect, norm_vect

    ! write(*,*) 'Rl', Rl, 'detJ', detJ, 'i_face', i_face
    ! write(*,*) 'NNODE_l', NNODE_l, 'NNODE_d', NNODE_d, 'MCRD', MCRD, 'dim_patchd', dim_patchd

    CMAT(:,:) = zero

    ! Assembling
    !! TODO Use 3 dimension for CMAT for sake of simplicity : CMAT(MCRD, NNODE_d, NNODE_l)
    !! TODO 2 : First direction with size MCRD is not necessary : it contains the same values
    dxidx(:,:) = zero
    call MatrixInv(dxidx(:2,:2), dxdxi(:2,:2), det_dxdxi, 2)



    !! Compute normal vector
    select case (i_face)
        case(1)
            tan_vect = dxdxi(:, 2) / norm2(dxdxi(:, 2))
        case(2)   
            tan_vect = - dxdxi(:, 2) / norm2(dxdxi(:, 2))
        case(3)
            tan_vect = dxdxi(:, 1) / norm2(dxdxi(:, 1))
        case(4)
            tan_vect = - dxdxi(:, 1) / norm2(dxdxi(:, 1))
    end select

    norm_vect(1) = tan_vect(2)
    norm_vect(2) = -tan_vect(1)
    norm_vect(3) = 0
    
    ! write(*,*)"Rl",Rl
    ! write(*,*)"dRdxid",dRdxid
    ! write(*,*)"dxidx",dxidx
    ! write(*,*)"detJ",detJ

    ! write(*,*)"norm_vect", norm_vect

    ! write(*,*) 'dRdxid'
    ! do i = 1, NNODE_d
    !     write(*,*) 'NNODE_d', i 
    !     write(*,*) (dRdxid(i, l), l=1, dim_patchd)
    ! enddo   

    ! write(*,*) 'dxidx'
    ! do l = 1, dim_patchd
    !     write(*,*) (dxidx(l,k), k=1, MCRD)
    ! enddo   

    ! write(*,*) 'norm_vect'
    ! do l = 1, MCRD
    !     write(*,*) norm_vect(l)
    ! enddo   

    ! write(*,*) 'dRdxid*norm'
    ! do i = 1, NNODE_d
    !     write(*,*) 'NNODE_d', i 
    !     do l = 1, dim_patchd
    !         write(*,*) 'dim_patchd', l
    !         write(*,*) dRdxid(i,l)*norm_vect(l)
    !     enddo
    ! enddo   

    count = 1 
    do j = 1, NNODE_l
        do i = 1, NNODE_d
            do l = 1, dim_patchd
                do m = 1, MCRD
                    do k = 1, MCRD
                        CMAT(k,count) = CMAT(k,count) + Rl(j)*dRdxid(i,l)* &
                        & dxidx(l,m)*norm_vect(m)*detJ
                    enddo
                enddo
            enddo
            count = count+1
        enddo
    enddo

    ! count = 1 
    ! do j = 1, NNODE_l
    !     do i = 1, NNODE_d
    !         do l = 1, dim_patchd
    !             do m = 1, MCRD
    !                 do k = 1, MCRD
    !                     CMAT(k,count) = CMAT(k,count) + Rl(j)*dRdxid(i,l)* &
    !                     & dxidx(l,m)*norm_vect(m)*detJ
    !                 enddo
    !             enddo
    !         enddo
    !         write(*,*) 'count', count
    !         write(*,*) 'CMAT', CMAT(:,count)
    !         write(*,*) '*************************************'
    !         count = count+1
    !     enddo
    ! enddo

    ! count = 1 
    ! do j = 1, NNODE_l
    !     do i = 1, NNODE_d
    !         do k = 1, MCRD
    !             CMAT(k,count) = CMAT(k,count) + dRdxid(i,k)
    !             if (dRdxid(i,k) .ne. 0) then
    !                 write(*,*) 'j', j 
    !                 write(*,*) 'i', i 
    !                 write(*,*) 'k', k
    !                 write(*,*) 'count', count
    !                 write(*,*) '*************************************'
    !             endif 
    !         enddo
    !         write(*,*) 'CMAT', CMAT(:,count)
    !         count = count+1
    !     enddo
    ! enddo

    ! count = 1 
    ! do j = 1, NNODE_l
    !     do i = 1, NNODE_d
    !         do l = 1, dim_patchd
    !             do k = 1, MCRD
    !                 CMAT(k,count) = CMAT(k,count) + dRdxid(i,l)* &
    !                 & norm_vect(l)
    !             enddo    
    !         enddo
    !         write(*,*) 'count', count
    !         write(*,*) 'CMAT', CMAT(:,count)
    !         write(*,*) '*************************************'
    !         count = count+1
    !     enddo
    ! enddo


    ! write(*,*)"CMAT",CMAT
endsubroutine cplingddispU5
    