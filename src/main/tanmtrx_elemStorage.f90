

subroutine build_tanMatrix(KTdata, KTrow, KTcol,  &
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
    integer,          intent(out) :: KTrow, KTcol
    double precision, intent(out) :: KTdata
    dimension KTdata(nb_data), KTrow(nb_data), KTcol(nb_data)


    !! Local variables
    !! ---------------

    !! for UELMAT routine
    integer :: NDOFEL
    double precision :: COORDS_elem, AMATRX, MAT_patch
    dimension COORDS_elem(MCRD, MAXVAL(NNODE)),  &
        &     AMATRX(MCRD, MCRD, MAXVAL(NNODE) * (MAXVAL(NNODE) + 1) / 2),      &
        &     MAT_patch(maxval(n_mat_props))

    double precision :: Uelem
    dimension Uelem(MCRD, MAXVAL(NNODE) * (MAXVAL(NNODE) + 1) / 2)

    !! Global stiffness matrix and force vector
    integer :: num_elem, i, j, JELEM, Numpatch, sctr, num_load,  &
        &      num_cp, ddl
    dimension sctr(MAXVAL(NNODE))

    !! Integers
    integer :: n, m, dofi, dofj, cpi, cpj, kk, ll, nnodeSum, count

    !! Initialisation

    !! Initialize K to zero and F to concentrated loads
    KTdata = zero
    KTrow  = 0
    KTcol  = 0

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
                    do j=1,MCRD
                        Uelem(j,i) = Uk((IEN_patch(i, num_elem)-1)*MCRD+j)
                    end do
                end do
                call extractNurbsElementInfos(num_elem)

                !! Compute elementary matrix 
                AMATRX = zero
                MAT_patch(:n_mat_props(NumPatch)) = MATERIAL_PROPERTIES(:n_mat_props(NumPatch), NumPatch)
                if (ELT_TYPE_patch == 'U1') then
                    !! Solid element
                    call UELTANMAT_byCP(NDOFEL, MCRD, nnode_patch, JELEM, &
                        &   NBINT(NumPatch), COORDS_elem(:,:nnode_patch),    &
                        &   TENSOR_patch, MAT_patch(:2), Uelem,            &
                        &   AMATRX(:,:,:nnodeSum))

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
                                KTdata(count) = AMATRX(kk, ll, i)
                                KTrow( count) = dofi + kk - 1
                                KTcol( count) = dofj + ll - 1
                                count = count + 1
                            end do
                        end do
                    end do

                    !! case cpi == cpj
                    i = i + 1
                    do ll = 1, MCRD
                        AMATRX(ll, ll, i) = AMATRX(ll, ll, i) * 0.5d0
                        do kk = 1, ll
                            KTdata(count) = AMATRX(kk, ll, i)
                            KTrow( count) = dofj + kk - 1
                            KTcol( count) = dofj + ll - 1
                            count = count + 1
                        end do
                    end do

                end do

            end if

        end do !! End of loop on element

        call deallocateMappingData()
        call finalizeNurbsPatch()


    end do ! End of loop on patch

end subroutine build_tanMatrix

