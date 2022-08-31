! ==========================
! module :: Heat transfer 
! author :: Joaquin Cornejo
! ==========================

module heat_transfer

    implicit none
    integer, parameter :: dimen = 3! By now only consider 3D 
    double precision, parameter :: span_tol = 1e-8

    contains

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
