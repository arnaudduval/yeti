!! Copyright 2018-2019 Thibaut Hirschler
!! Copyright 2020-2023 Arnaud Duval
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


subroutine rigid_motion(u,  &
            &   patch_idx, translation, rotation_center, rotation_angle, &
            &   activeElement,nb_data, &
            &   COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight, &
            &   Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,n_mat_props,RHO,TENSOR,  &
            &   indDLoad,JDLType,ADLMAG,load_target_nbelem, &
            &   load_additionalInfos, nb_load_additionalInfos, bc_values,   &
            &   nb_bc,bc_target, &
            &   bc_target_nbelem,ind_dof_free,nb_dof_free,MCRD,NBINT,nb_load,   &
            &   nb_patch,nb_elem,nnode,nb_cp,nb_dof_tot,nodal_dist, &
            &   nb_n_dist, nb_cp_n_dist, n_mat_props_max)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------
    double precision, intent(in) :: translation, rotation_center, rotation_angle
    integer, intent(in) :: patch_idx
    dimension translation(2), rotation_center(2)

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
    double precision, intent(in) :: MATERIAL_PROPERTIES,RHO,PROPS
    integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,IEN,   &
        &     nb_elem_patch, JPROPS
    integer, intent(in) :: n_mat_props      !! Number of material properties per patch
    dimension n_mat_props(nb_patch)
    integer, intent(in) :: n_mat_props_max  !! Maximum value of n_mat_props
    dimension MATERIAL_PROPERTIES(n_mat_props_max,nb_patch),  &
        &     RHO(nb_patch),   &
        &     PROPS(:),    &
        &     NNODE(nb_patch), &
        &     IEN(:),  &
        &     nb_elem_patch(nb_patch), &
        &     JPROPS(nb_patch),    &
        &     NBINT(nb_patch)


    !! Loads
    double precision, intent(in) :: ADLMAG,load_additionalInfos, &
        &     nodal_dist
    integer, intent(in) :: nb_load,indDLoad,JDLType,load_target_nbelem
    integer, intent(in) :: nb_n_dist, nb_cp_n_dist
    integer, intent(in) :: nb_load_additionalInfos
    dimension ADLMAG(nb_load),  &
        &     load_additionalInfos(:),  &
        &     nodal_dist(nb_n_dist, nb_cp_n_dist),  &
        &     indDLoad(:),  &
        &     JDLType(nb_load), &
        &     load_target_nbelem(nb_load),   &
        &     nb_load_additionalInfos(:)


    !! Boundary Conditions
    double precision, intent(in) :: bc_values
    integer, intent(in) :: nb_bc,bc_target,bc_target_nbelem
    dimension bc_values(2,nb_bc),   &
        &     bc_target(:),         &
        &     bc_target_nbelem(nb_bc)


    !! Degrees Of Freedom
    integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
    dimension ind_dof_free(nb_dof_tot)


    !! Storage infos
    integer, intent(in) :: nb_data,activeElement
    dimension activeElement(nb_elem)


    !! Output variables
    !! ----------------

    !! linear system to solve
    double precision, intent(out) :: u
    dimension u(nb_dof_tot)


    !! Local variables
    !! ---------------

    !! for UELMAT routine
    integer :: NDOFEL
    double precision :: COORDS_elem, RHS, AMATRX, MAT_patch, &
        &               n_dist_elem
    dimension COORDS_elem(MCRD,MAXVAL(NNODE)),RHS(MCRD*MAXVAL(NNODE)),  &
        &     AMATRX(MCRD,MCRD,MAXVAL(NNODE)*(MAXVAL(NNODE)+1)/2),      &
        &     MAT_patch(maxval(n_mat_props)), n_dist_elem(nb_n_dist,MAXVAL(NNODE))

    !! Global stiffness matrix and force vector
    integer :: num_elem,i,j,JELEM,Numpatch, sctr,num_load,  &
        &      num_cp, ddl
    dimension sctr(MAXVAL(NNODE))

    double precision :: rotation_matrix, rotation
    dimension rotation_matrix(2,2), rotation(2)

    !! Integers
    integer :: n,m,dofi,dofj,ddli,ddlj,cpi,cpj, kk,ll, nnodeSum,count


    !! Initialisation

    rotation_matrix(1,1) = cos(rotation_angle)
    rotation_matrix(1,2) = -sin(rotation_angle)
    rotation_matrix(2,1) = sin(rotation_angle)
    rotation_matrix(2,2) = cos(rotation_angle)
    
    !! Initialize K to zero and F to concentrated loads
    u     = zero

    !! Assembly

    count = 1
    JELEM = 0
    do NumPatch = 1,nb_patch
        if (Numpatch == patch_idx) then
            call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv, &
                &        weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,  &
                &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            if ((ELT_TYPE_patch .eq. 'U30') .or. (ELT_TYPE_patch .eq. 'U10')) then
                i = int(PROPS_patch(2))
                call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
                    &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)
            endif

            NDOFEL   = nnode_patch*MCRD
            nnodeSum = nnode_patch*(nnode_patch+1)/2

            !! Loop on elments
            do num_elem = 1,nb_elem_patch(NumPatch)
                JELEM = JELEM + 1
                ! write(*,*) "*****Element ", JELEM, "********"
                if (activeElement(JELEM)==1) then

                    do i = 1,nnode_patch
                        ! write(*,*) i
                        COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
                        n_dist_elem(:,i) = nodal_dist(:,IEN_patch(i,num_elem))
                        ! write(*,*) COORDS_elem(:,i)
                    enddo
                    call extractNurbsElementInfos(num_elem)

                    sctr(:nnode_patch) = IEN_patch(:,num_elem)
                    
                    !! Update Load Vector
                    do i = 1,nnode_patch
                        ddli  = (sctr(i)-1)*MCRD
                        ! write(*,*), ddli, matmul(rotation_matrix,COORDS3D(:MCRD,sctr(i)))
                        ! write(*,*), ddli, matmul(rotation_matrix,COORDS3D(:MCRD,sctr(i))) - COORDS3D(:MCRD,sctr(i))
                        ! rotation = matmul(rotation_matrix,(COORDS3D(:MCRD,sctr(i)) - rotation_center)) & 
                        ! & - COORDS3D(:MCRD,sctr(i)) + rotation_center ! Grande rotations
                        rotation(1) = rotation_angle * (rotation_center(2) - COORDS3D(2,sctr(i))) ! petites rotations
                        rotation(2) = rotation_angle * (-rotation_center(1) + COORDS3D(1,sctr(i)))
                        ! write(*,*), ddli, rotation
                        do kk = 1,MCRD
                            u(ddli+kk) = translation(kk) + rotation(kk)
                            ! write(*,*), kk, rotation(kk)
                            !! + COORDS3D(:MCRD,sctr(i)) translation uniquement 
                            !! METTRE A JOUR COORDS 3D VALEUR INITIALE PLUS U (ou decaler les patchs ?)
                        enddo
                    enddo

                endif

            enddo !! End of loop on element

            call deallocateMappingData()
            call finalizeNurbsPatch()
        endif

    enddo ! End of loop on patch

end subroutine rigid_motion
