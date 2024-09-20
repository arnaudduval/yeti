

subroutine compute_intVector(FINT,  &
        &   Uk,activeElement,COORDS3D, &
        &   Nkv,Ukv,Nijk,weight,Jpqr,nb_elem_patch, &
        &   NBINT,IEN,TENSOR, ELT_TYPE,MATERIAL_PROPERTIES,n_mat_props, PROPS,JPROPS,   &
        &   MCRD, NNODE, &
        &   nb_data,n_mat_props_max,nb_patch,nb_cp,nb_elem,nb_dof_tot)

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
    double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS

    integer, intent(in) :: MCRD, NNODE, nb_patch, nb_elem, NBINT, IEN,   &
        &     nb_elem_patch, JPROPS
    integer, intent(in) :: n_mat_props      !! Number of material properties per patch
    dimension n_mat_props(nb_patch)
    integer, intent(in) :: n_mat_props_max  !! Maximum value of n_mat_props
    dimension MATERIAL_PROPERTIES(n_mat_props_max, nb_patch),  &
        &     PROPS(:),    &
        &     NNODE(nb_patch), &
        &     IEN(:),  &
        &     nb_elem_patch(nb_patch), &
        &     JPROPS(nb_patch),    &
        &     NBINT(nb_patch)



    ! Degrees Of Freedom
    integer, intent(in) :: nb_dof_tot
    double precision, intent(in) :: Uk
    dimension Uk(nb_dof_tot)


    ! !! Storage infos
    integer, intent(in) :: nb_data
    integer, intent(in) :: activeElement
    dimension activeElement(nb_elem)

    !! Output variables
    !! ----------------

    !! linear system to solve
    double precision, intent(out) :: FINT
    dimension FINT(nb_dof_tot)


    !! Local variables
    !! ---------------

    !! for UELMAT routine
    integer :: NDOFEL
    double precision :: COORDS_elem, FINT_elem, MAT_patch
    dimension COORDS_elem(MCRD, MAXVAL(NNODE)),  &
        &     FINT_elem(MCRD, MAXVAL(NNODE)),      &
        &     MAT_patch(maxval(n_mat_props))

    double precision :: Uelem
    dimension Uelem(MCRD, MAXVAL(NNODE))

    !! Global stiffness matrix and force vector
    integer :: num_elem, i, j, JELEM, Numpatch, sctr, num_load,  &
        &      num_cp, ddl
    dimension sctr(MAXVAL(NNODE))

    !! Integers
    integer :: n, m, dofi, dofj, cpi, cpj, kk, ll, nnodeSum, count

    !! Initialisation
    Fint = zero


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
                sctr(:nnode_patch) = IEN_patch(:, num_elem)
                do i = 1, nnode_patch
                    dofi = (sctr(i) - 1) * MCRD
                    COORDS_elem(:, i) = COORDS3D(:MCRD, sctr(i))
                    do kk=1,MCRD
                        Uelem(kk,i) = Uk(dofi + kk)
                    end do
                end do
                call extractNurbsElementInfos(num_elem)

                !! Compute elementary matrix 
                FINT_elem = zero
                MAT_patch(:n_mat_props(NumPatch)) = MATERIAL_PROPERTIES(:n_mat_props(NumPatch), NumPatch)
                if (ELT_TYPE_patch == 'U1') then
                    !! Solid element
                    call UELINTVEC_byCP(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, MAT_patch(:2), Uelem,            &
                        &   FINT_elem(:,:nnodeSum))

                else
                    write(*,*) 'Element'// ELT_TYPE_patch //' not available.'
                end if

            
                !! Assemble FINT to global vector
                sctr(:nnode_patch) = IEN_patch(:, num_elem)
                do i = 1, nnode_patch
                    dofi = (sctr(i) - 1) * MCRD
                    do kk = 1, MCRD
                        FINT(dofi + kk) = FINT(dofi + kk) + FINT_elem(kk,i)
                    end do
                end do

            end if

        end do !! End of loop on element

        call deallocateMappingData()
        call finalizeNurbsPatch()


    end do ! End of loop on patch

end subroutine compute_intVector

