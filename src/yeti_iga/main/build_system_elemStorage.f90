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

!! Build global sparse stiffness matrix for scipy.sparse
!! Inputs are handled by Python object IGAparametrization
!! Outputs : Kdata : vector containing element matrices
!! Krow : matrix line indices
!! Kcol : matrix column indices

subroutine sys_linmat_lindef_static(Kdata, Krow, Kcol, F,  &
        &   activeElement, nb_data, &
        &   COORDS3D, IEN, nb_elem_patch, Nkv, Ukv, Nijk, weight, &
        &   Jpqr, ELT_TYPE, PROPS, JPROPS, MATERIAL_PROPERTIES, n_mat_props, RHO, TENSOR,  &
        &   indDLoad, JDLType, ADLMAG, load_target_nbelem, &
        &   load_additionalInfos, nb_load_additionalInfos, bc_values,   &
        &   nb_bc, bc_target, &
        &   bc_target_nbelem, ind_dof_free, nb_dof_free, MCRD, NBINT, nb_load,   &
        &   nb_patch, nb_elem, nnode, nb_cp, nb_dof_tot, nodal_dist, &
        &   nb_n_dist, nb_cp_n_dist, n_mat_props_max)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------

    !! NURBS geometry
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3, nb_cp)

    double precision, intent(in) :: Ukv, weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem),   &
        &     Ukv(:), weight(:)


    !! Patches and Elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: MATERIAL_PROPERTIES, RHO, PROPS
    integer, intent(in) :: MCRD, NNODE, nb_patch, nb_elem, NBINT, IEN,   &
        &     nb_elem_patch, JPROPS
    integer, intent(in) :: n_mat_props      !! Number of material properties per patch
    dimension n_mat_props(nb_patch)
    integer, intent(in) :: n_mat_props_max  !! Maximum value of n_mat_props
    dimension MATERIAL_PROPERTIES(n_mat_props_max, nb_patch),  &
        &     RHO(nb_patch),   &
        &     PROPS(:),    &
        &     NNODE(nb_patch), &
        &     IEN(:),  &
        &     nb_elem_patch(nb_patch), &
        &     JPROPS(nb_patch),    &
        &     NBINT(nb_patch)


    !! Loads
    double precision, intent(in) :: ADLMAG, load_additionalInfos, &
        &     nodal_dist
    integer, intent(in) :: nb_load, indDLoad, JDLType, load_target_nbelem
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
    integer, intent(in) :: nb_bc, bc_target, bc_target_nbelem
    dimension bc_values(2, nb_bc),   &
        &     bc_target(:),         &
        &     bc_target_nbelem(nb_bc)


    !! Degrees Of Freedom
    integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
    dimension ind_dof_free(nb_dof_tot)


    !! Storage infos
    integer, intent(in) :: nb_data, activeElement
    dimension activeElement(nb_elem)


    !! Output variables
    !! ----------------

    !! linear system to solve
    integer,          intent(out) :: Krow, Kcol
    double precision, intent(out) :: Kdata, F
    dimension Kdata(nb_data), Krow(nb_data), Kcol(nb_data), F(nb_dof_tot)


    !! Local variables
    !! ---------------

    !! for UELMAT routine
    integer :: NDOFEL
    double precision :: COORDS_elem, RHS, AMATRX, MAT_patch, &
        &               n_dist_elem
    dimension COORDS_elem(MCRD, MAXVAL(NNODE)), RHS(MCRD * MAXVAL(NNODE)),  &
        &     AMATRX(MCRD, MCRD, MAXVAL(NNODE) * (MAXVAL(NNODE) + 1) / 2),      &
        &     MAT_patch(maxval(n_mat_props)), n_dist_elem(nb_n_dist, MAXVAL(NNODE))

    !! Global stiffness matrix and force vector
    integer :: num_elem, i, j, JELEM, Numpatch, sctr, num_load,  &
        &      num_cp, ddl
    dimension sctr(MAXVAL(NNODE))

    !! Integers
    integer :: n, m, dofi, dofj, ddli, ddlj, cpi, cpj, kk, ll, nnodeSum, count


    !! Initialisation

    !! Initialize K to zero and F to concentrated loads
    Kdata = zero
    Krow  = 0
    Kcol  = 0
    F     = zero
    kk    = 0
    do num_load = 1, nb_load
        i = JDLType(num_load)
        if (i / 10 < 1) then
            do num_cp = 1, load_target_nbelem(num_load)
                ddl = (indDLoad(kk + num_cp) - 1) * MCRD + i
                F(ddl) = ADLMAG(num_load)
            end do
        end if
        kk = kk + load_target_nbelem(num_load)
    end do

    !! Assembly

    count = 1
    JELEM = 0
    do NumPatch = 1, nb_patch

        call extractNurbsPatchGeoInfos(NumPatch, Nkv, Jpqr, Nijk, Ukv, &
            &        weight, nb_elem_patch)
        call extractNurbsPatchMechInfos(NumPatch, IEN, PROPS, JPROPS,  &
            &        NNODE, nb_elem_patch, ELT_TYPE, TENSOR)

        if ((ELT_TYPE_patch == 'U30') .or. (ELT_TYPE_patch == 'U10')) then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i, nb_elem_patch, Nkv, Jpqr, Nijk, Ukv, &
                &           weight, IEN, PROPS, JPROPS, NNODE, ELT_TYPE, TENSOR)
        end if

        NDOFEL   = nnode_patch * MCRD
        nnodeSum = nnode_patch * (nnode_patch + 1) / 2

        !! Loop on elments
        do num_elem = 1, nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            ! write(*,*) "*****Element ", JELEM, "********"
            if (activeElement(JELEM) == 1) then

                do i = 1, nnode_patch
                    COORDS_elem(:, i) = COORDS3D(:MCRD, IEN_patch(i, num_elem))
                    n_dist_elem(:, i) = nodal_dist(:, IEN_patch(i, num_elem))
                end do
                call extractNurbsElementInfos(num_elem)


                !! Compute elementary matrix and load vector
                RHS    = zero
                AMATRX = zero
                MAT_patch(:n_mat_props(NumPatch)) = MATERIAL_PROPERTIES(:n_mat_props(NumPatch), NumPatch)
                if (ELT_TYPE_patch == 'U0') then
                    !! Void element --> do nothing

                elseif (ELT_TYPE_patch == 'U1') then
                    !! Solid element
                    call UELMAT_byCP(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, MAT_patch(:2), RHO(NumPatch), nb_load,   &
                        &   indDLoad, load_target_nbelem, JDLType, ADLMAG,     &
                        &   load_additionalInfos, SIZE(load_additionalInfos),&
                        &   nb_load_additionalInfos,        &
                        &   n_dist_elem, nb_n_dist, RHS(:NDOFEL),             &
                        &   AMATRX(:,:,:nnodeSum))
                elseif (ELT_TYPE_patch == 'U99') then
                    !! Solid element with high order elastic behaviour law
                    !! Name U99 is set for temporary development
                    !! For this type of element we take into account material
                    !! behaviour with more than 2 material properties
                    !! TODO : implement material behaviour with more than 2 properties
                    !! for other element types
                    !! NOTE : for this element, RHO (density) is not used
                    call UELMAT_HO_byCP(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, MAT_patch, n_mat_props(NumPatch), RHO(NumPatch), nb_load,   &
                        &   indDLoad, load_target_nbelem, JDLType, ADLMAG, &
                        &   load_additionalInfos, SIZE(load_additionalInfos),&
                        &   n_dist_elem, nb_n_dist, RHS(:NDOFEL),             &
                        &   AMATRX(:,:,:nnodeSum))
                elseif (ELT_TYPE_patch == 'U98') then
                    !! Same as U99 but with only 1st gradient of strain taken into account
                    call UELMAT_HO_byCP_1stG(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, MAT_patch, n_mat_props(NumPatch), RHO(NumPatch), nb_load,   &
                        &   indDLoad, load_target_nbelem, JDLType, ADLMAG, &
                        &   load_additionalInfos, SIZE(load_additionalInfos),&
                        &   n_dist_elem, nb_n_dist, RHS(:NDOFEL),             &
                        &   AMATRX(:,:,:nnodeSum))
                elseif (ELT_TYPE_patch == 'U3') then
                    !! Shell element
                    call UELMAT3_byCP(NDOFEL, MCRD, nnode_patch, JELEM,        &
                        &   NBINT(NumPatch), COORDS_elem, TENSOR_patch, MAT_patch, &
                        &   PROPS_patch, JPROPS_patch, nb_load, indDLoad,          &
                        &   load_target_nbelem, JDLType, ADLMAG, RHS(:NDOFEL),     &
                        &   AMATRX(:,:,:nnodeSum))

                elseif (ELT_TYPE_patch == 'U10') then
                    !! Embedded solid element
                    call UELMAT10(NDOFEL, MCRD, nnode_patch, nnode_map, nb_cp,      &
                        &   JELEM, NBINT(NumPatch), COORDS_elem(:,:nnode_patch),  &
                        &   COORDS3D, TENSOR_patch, MAT_patch, RHO(NumPatch),      &
                        &   nb_load, indDLoad, load_target_nbelem, JDLType, ADLMAG, &
                        &   load_additionalInfos, SIZE(load_additionalInfos),    &
                        &   nb_load_additionalInfos,                            &
                        &   n_dist_elem, nb_n_dist,                              &
                        &   RHS(:NDOFEL), AMATRX(:,:,:nnodeSum))

                elseif (ELT_TYPE_patch == 'U30') then
                    !! embedded shell element
                    call UELMAT30(NDOFEL, MCRD, nnode_patch, nnode_map, nb_cp,      &
                        &   JELEM, NBINT(NumPatch), COORDS_elem, COORDS3D,         &
                        &   TENSOR_patch, MAT_patch, PROPS_patch, JPROPS_patch,    &
                        &   nb_load, indDLoad, load_target_nbelem, JDLType, ADLMAG, &
                        &   RHS(:NDOFEL), AMATRX(:,:,:nnodeSum) )

                else
                    write(*,*) 'Element'// ELT_TYPE_patch //' not available.'
                end if


                !! Assemble AMATRX to global stiffness matrix K
                sctr(:nnode_patch) = IEN_patch(:, num_elem)

                i = 0
                do cpj = 1, nnode_patch
                    dofj = (sctr(cpj) - 1) * MCRD

                    !! case cpi < cpj
                    do cpi = 1, cpj - 1
                        dofi = (sctr(cpi) - 1) * MCRD
                        i   = i + 1
                        do ll = 1, MCRD
                            do kk = 1, MCRD
                                Kdata(count) = AMATRX(kk, ll, i)
                                Krow( count) = dofi + kk - 1
                                Kcol( count) = dofj + ll - 1
                                count = count + 1
                            end do
                        end do
                    end do

                    !! case cpi == cpj
                    i = i + 1
                    do ll = 1, MCRD
                        AMATRX(ll, ll, i) = AMATRX(ll, ll, i) * 0.5d0
                        do kk = 1, ll
                            Kdata(count) = AMATRX(kk, ll, i)
                            Krow( count) = dofj + kk - 1
                            Kcol( count) = dofj + ll - 1
                            count = count + 1
                        end do
                    end do

                end do

                !! Update Load Vector
                dofi = 0
                do i = 1, nnode_patch
                    ddli  = (sctr(i) - 1) * MCRD
                    do kk = 1, MCRD
                        F(ddli + kk) = F(ddli + kk) + RHS(dofi + kk)
                    end do
                    dofi = dofi + MCRD
                end do

            end if

        end do !! End of loop on element

        call deallocateMappingData()
        call finalizeNurbsPatch()


    end do ! End of loop on patch

end subroutine sys_linmat_lindef_static
