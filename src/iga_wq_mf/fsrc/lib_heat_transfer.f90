! ==========================
! module :: Heat transfer 
! author :: Joaquin Cornejo
! ==========================

subroutine clean_dirichlet_1dim(nr, A, ndod, dod)
    !! Set to 0 (Dirichlet condition) the values of an array using the dod indices 

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ndod
    double precision, intent(inout) :: A
    dimension :: A(nr)

    integer, intent(in) :: dod
    dimension :: dod(ndod)

    ! Clean array
    A(dod) = 0.d0 

end subroutine clean_dirichlet_1dim

subroutine create_block_L(nr, ndod, dod, indi_L, indj_L, L, indi_LT, indj_LT, LT)
    !! Creates the block matrix L of size nrxnc. 
    !! This matrix satisfies : L M L' = Mnn (n are the free control points)
    !! Returns L and LT in CSR format
    !! In this matrix, the number of non zero values is equal to the number of rows (nc = nr + ndod, nr = ndof)

    implicit none
    ! Input/output data
    ! -----------------
    integer, intent(in) :: nr, ndod, dod
    dimension :: dod(ndod)

    integer, intent(out) :: indi_L, indj_L, indi_LT, indj_LT
    dimension :: indi_L(nr+1), indj_L(nr), indi_LT(nr+ndod+1), indj_LT(nr)
    double precision, intent(out) :: L, LT
    dimension :: L(nr), LT(nr)
    
    ! Local data
    ! ----------
    integer :: i, j, nc, indi_coo, indj_coo
    dimension :: indi_coo(nr), indj_coo(nr)
    integer, allocatable, dimension(:) :: dof
    double precision :: data_coo
    dimension :: data_coo(nr, 1)

    ! Initialize 
    nc = nr + ndod
    if (any(dod.ge.nc+1)) stop 'Problem creating C'
    
    ! Get dof as complement of dod
    allocate(dof(nr))
    dof = 1; i = 1; j = 1
    do while ((j.le.nc).and.(i.le.nr))
        if (any(dod.eq.j)) then
            continue
        else
            dof(i) = j
            i = i + 1 
        end if
        j = j + 1
    end do

    ! Get COO format
    do i = 1, nr ! ndof = nr
        indi_coo(i) = i
        indj_coo(i) = dof(i)
        data_coo(i, 1) = 1.d0 
    end do

    ! Convert L in COO to CSR format
    call coo2csr(1, nr, nr, data_coo, indi_coo, indj_coo, L, indj_L, indi_L)

    ! Convert L' in COO to CSR format
    call coo2csr(1, nc, nr, data_coo, indj_coo, indi_coo, LT, indj_LT, indi_LT)

end subroutine create_block_L

module heat_transfer

    implicit none
    integer, parameter :: dimen = 3
    double precision, parameter :: span_tol = 1e-8

    contains

    subroutine interpolate_temperature_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_B_u, data_B_v, data_B_w, T_ctrlpts, T_interp)
        !! Computes strain in 3D (from parametric space to physical space)
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------  
        integer, intent(in) ::  nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
        double precision, intent(in) :: T_ctrlpts
        dimension :: T_ctrlpts(nr_u*nr_v*nr_w)

        double precision, intent(out) :: T_interp
        dimension :: T_interp(nc_u*nc_v*nc_w)

        ! Local data
        !-----------
        integer :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        ! Convert CSR to CSC
        call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
        call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
        call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

        ! Interpolation
        call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                            T_ctrlpts, T_interp)

    end subroutine interpolate_temperature_3d

    subroutine compute_heat_properties(nr, table_conductivity, table_capacity, &
                                    nnz, temperature, conductivity, capacity)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nr, nnz
        double precision, intent(in) :: table_conductivity, table_capacity, temperature
        dimension ::  table_conductivity(nr, 2), table_capacity(nr, 2), temperature(nnz)

        double precision, intent(out) :: conductivity, capacity
        dimension :: conductivity(nnz), capacity(nnz)

        call linear_interpolation(nr, table_conductivity, nnz, temperature, conductivity, span_tol)
        call linear_interpolation(nr, table_capacity, nnz, temperature, capacity, span_tol)

    end subroutine compute_heat_properties

    subroutine compute_heat_coefficients(nc_total, conductivity, capacity, invJ, detJ, Kcoefs, Ccoefs)
        !! Computes the coefficients to use to calculate source vector and tangent heat matrix
        !! By the moment, it only considers isotropic materials

        implicit none 
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nc_total
        double precision, intent(in) :: conductivity, capacity, invJ, detJ
        dimension :: conductivity(nc_total), capacity(nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: Kcoefs, Ccoefs
        dimension :: Kcoefs(dimen, dimen, nc_total), Ccoefs(nc_total)

        ! Local data
        ! ----------
        double precision :: invJJt, detJJt
        dimension :: invJJt(dimen, dimen)
        integer :: k

        do k = 1, nc_total

            ! Initialize
            invJJt = invJ(:, :, k)
            detJJt = detJ(k)

            ! Compute coefficients
            Kcoefs(:, :, k) = conductivity(k) * matmul(invJJt, transpose(invJJt)) * detJJt
            Ccoefs(k) = capacity(k) * detJJt

        end do

    end subroutine compute_heat_coefficients

end module heat_transfer
