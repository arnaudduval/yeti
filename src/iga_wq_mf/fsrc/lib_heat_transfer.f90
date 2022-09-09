! ==========================
! module :: Heat transfer 
! author :: Joaquin Cornejo
! ==========================

module heat_transfer

    implicit none
    integer, parameter :: dimen = 3! By now only consider 3D 
    double precision, parameter :: span_tol = 1e-8

    contains

    subroutine interpolate_temperature_field_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                            data_B_u, data_B_v, data_B_w, Tin, Tout)

        implicit none 
        ! Input/ output
        ! --------------------  
        integer, intent(in) ::  nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
        double precision, intent(in) :: Tin
        dimension :: Tin(nr_u*nr_v*nr_w)

        double precision, intent(out) :: Tout
        dimension :: Tout(nc_u*nc_v*nc_w)

        ! Local data
        !-----------------
        integer :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        ! Initialize
        call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
        call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
        call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

        ! Interpolation
        call sumproduct3d_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                                    nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                    nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                                    Tin, Tout)

    end subroutine interpolate_temperature_field_3d

    subroutine compute_heat_properties(nbpts, table_conductivity, table_capacity, nnz, temperature, KK, CC)

        implicit none 
        ! Input / output data
        ! ---------------------
        integer, intent(in) :: nbpts, nnz
        double precision, intent(in) :: table_conductivity, table_capacity, temperature
        dimension ::  table_conductivity(nbpts, 2), table_capacity(nbpts, 2), temperature(nnz)

        double precision, intent(out) :: KK, CC
        dimension :: KK(nnz), CC(nnz)

        call linear_interpolation(nbpts, table_conductivity, nnz, temperature, KK, span_tol)
        call linear_interpolation(nbpts, table_capacity, nnz, temperature, CC, span_tol)

    end subroutine compute_heat_properties

    subroutine compute_heat_coefficients(nc_total, KK, CC, invJ, detJ, Kcoefs, Ccoefs)
        !! Computes the coefficients to use in internal "force" vector and tangent heat matrix
        !! By now, only isotropic materials

        implicit none 
        ! Input / output 
        ! --------------------
        integer, intent(in) :: nc_total
        double precision, intent(in) :: KK, CC, invJ, detJ
        dimension :: KK(nc_total), CC(nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: Kcoefs, Ccoefs
        dimension :: Kcoefs(dimen, dimen, nc_total), Ccoefs(nc_total)

        ! Local data
        ! -------------
        double precision :: invJJt, detJJt
        dimension :: invJJt(dimen, dimen)
        integer :: k

        do k = 1, nc_total

            ! Initialize
            invJJt = invJ(:, :, k)
            detJJt = detJ(k)

            Kcoefs(:, :, k) = KK(k) * matmul(invJJt, transpose(invJJt)) * detJJt
            Ccoefs(k) = CC(k) * detJJt

        end do

    end subroutine compute_heat_coefficients

end module heat_transfer
