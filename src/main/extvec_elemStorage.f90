

subroutine build_extvector(F,  &
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
    double precision, intent(out) ::  F
    dimension F(nb_dof_tot)


    !! Local variables
    !! ---------------

    !! for UELMAT routine
    integer :: NDOFEL
    double precision :: COORDS_elem, RHS,n_dist_elem
    dimension COORDS_elem(MCRD, MAXVAL(NNODE)), RHS(MCRD * MAXVAL(NNODE)), n_dist_elem(nb_n_dist, MAXVAL(NNODE))

    !! Global stiffness matrix and force vector
    integer :: num_elem, i, j, JELEM, Numpatch, sctr, num_load,  &
        &      num_cp, ddl
    dimension sctr(MAXVAL(NNODE))

    !! Integers
    integer :: n, m, dofi, dofj, ddli, ddlj, cpi, cpj, kk, ll, nnodeSum, count


    !! Initialisation

    !! Initialize K to zero and F to concentrated loads
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

        NDOFEL   = nnode_patch * MCRD
        nnodeSum = nnode_patch * (nnode_patch + 1) / 2

        !! Loop on elments
        do num_elem = 1, nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            if (activeElement(JELEM) == 1) then

                do i = 1, nnode_patch
                    COORDS_elem(:, i) = COORDS3D(:MCRD, IEN_patch(i, num_elem))
                    n_dist_elem(:, i) = nodal_dist(:, IEN_patch(i, num_elem))
                end do
                call extractNurbsElementInfos(num_elem)


                !! Compute elementary matrix and load vector
                RHS    = zero
                if (ELT_TYPE_patch == 'U0') then
                    !! Void element --> do nothing

                elseif (ELT_TYPE_patch == 'U1') then
                    !! Solid element
                    call UELVEC_byCP(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, RHO(NumPatch), nb_load,   &
                        &   indDLoad, load_target_nbelem, JDLType, ADLMAG,     &
                        &   load_additionalInfos, SIZE(load_additionalInfos),&
                        &   nb_load_additionalInfos,        &
                        &   n_dist_elem, nb_n_dist, RHS(:NDOFEL))
                else
                    write(*,*) 'Element'// ELT_TYPE_patch //' not available.'
                end if

                sctr(:nnode_patch) = IEN_patch(:, num_elem)

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

end subroutine build_extvector

