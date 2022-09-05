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

! module wq_mf_pretreatment

!     implicit none

!     type :: paraprop

!         ! Inputs 
!         ! ------------
!         integer, dimension(:), pointer :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
!         double precision, dimension(:, :), pointer ::  data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w

!         ! Local 
!         ! -----------
!         integer :: nnz_u, nnz_v, nnz_w, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nr_total, nc_total

!         ! Outputs
!         ! -----------
!         integer, dimension(:), pointer :: indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w
!         double precision, dimension(:, :), pointer ::  data_BT_u, data_BT_v, data_BT_w

!     end type paraprop

!     contains

!     subroutine erase_row_direction(nr_in, nnz_in, indi_in, indj_in, data_in, table, &
!                                     nr_out, nnz_out, indi_out, indj_out, data_out)

!         implicit none
!         ! Input / output
!         ! ----------------
!         integer, intent(in) :: nr_in, nnz_in
!         integer, intent(in) :: indi_in, indj_in
!         dimension :: indi_in(nr_in+1), indj_in(nnz_in)
!         double precision, intent(in) :: data_in
!         dimension :: data_in(nnz_in, 6)
!         integer, intent(in) :: table
!         dimension :: table(2)

!         integer, intent(inout) :: nr_out, nnz_out
!         integer, dimension(*) :: indi_out, indj_out
!         double precision, dimension(*) :: data_out

!         ! Local data
!         ! --------------
!         logical :: doerase
!         integer :: nr2er, i
!         integer, allocatable, dimension(:) :: rows2er

!         ! Gat rows to be erased
!         nr2er =sum(table)
!         doerase = .true.
!         if (nr2er.eq.0) then 
!             doerase = .false.
!         else if (nr2er.eq.2) then
!             allocate(rows2er(2))
!             rows2er = (/1, nr_in/)
!         else if (nr2er.eq.1) then
!             allocate(rows2er(1))
!             if (table(1).eq.1) then
!                 rows2er = (/1/)
!             else
!                 rows2er = (/nr_in/)
!             end if
!         else
!             stop 'Something when wrong in table'
!         end if

!         ! Get final rows
!         nr_out = nr_in - nr2er

!         if (nnz_out.le.0) then
!             ! Get number of rows and nnz
!             call erase_rows_csr(nr2er, rows2er, nr_in, 6, nnz_in, indi_in, indj_in, data_in, &
!                                 nnz_out, indi_out, indj_out, data_out)

!         else
!             ! Erase rows
!             if (doerase) then 
!                 call erase_rows_csr(nr2er, rows2er, nr_in, 6, nnz_in, indi_in, indj_in, data_in, &
!                                 nnz_out, indi_out, indj_out, data_out)
!             else
!                 do i = 1, nr_in+1
!                     indi_out(i) = indi_in(i)
!                 end do

!                 do i = 1, nnz_in
!                     indj_out(i) = indj_in(i)
!                 end do

!             end if

!         end if

!     end subroutine erase_row_direction

!     subroutine paraprop_initialize(obj, nr_ui, nc_ui, nr_vi, nc_vi, nr_wi, nc_wi, nnz_ui, nnz_vi, nnz_wi, &
!                                 indi_ui, indj_ui, indi_vi, indj_vi, indi_wi, indj_wi, &
!                                 data_B_ui, data_B_vi, data_B_wi, data_W_ui, data_W_vi, data_W_wi, table)

!         implicit none
!         ! Input / output
!         ! ---------------
!         integer, parameter :: d = 3
!         integer, intent(in) :: nr_ui, nc_ui, nr_vi, nc_vi, nr_wi, nc_wi, nnz_ui, nnz_vi, nnz_wi
!         integer, intent(in) :: indi_ui, indj_ui, indi_vi, indj_vi, indi_wi, indj_wi
!         dimension ::    indi_ui(nr_ui+1), indj_ui(nnz_ui), &
!                         indi_vi(nr_vi+1), indj_vi(nnz_vi), &
!                         indi_wi(nr_wi+1), indj_wi(nnz_wi)
!         double precision, intent(in) :: data_B_ui, data_W_ui, data_B_vi, data_W_vi, data_B_wi, data_W_wi
!         dimension ::    data_B_ui(nnz_ui, 2), data_W_ui(nnz_ui, 4), &
!                         data_B_vi(nnz_vi, 2), data_W_vi(nnz_vi, 4), &
!                         data_B_wi(nnz_wi, 2), data_W_wi(nnz_wi, 4)
!         integer, intent(in) :: table
!         dimension :: table(d, 2)
!         type(paraprop), pointer :: obj

!         ! Local data
!         ! ----------------
!         integer :: nr_in, nnz_in, nr_out, nnz_out
!         dimension :: nr_in(d), nnz_in(d)
!         integer, allocatable, dimension(:) :: indi_out, indj_out
!         double precision, allocatable, dimension(:, :) :: a_in, a_out

!         ! Stock data
!         allocate(obj)
!         nr_in = (/nr_ui, nr_vi, nr_wi/)
!         nnz_in = (/nnz_ui, nnz_vi, nnz_wi/)

!         ! Erase data in direction u 
!         ! -------------------------
!         nr_out = -1; nnz_out = -1
!         call erase_row_direction(nr_in(1), nnz_in(1), indi_ui, indj_ui, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         allocate(indi_out(nr_out+1), indj_out(nnz_out), a_out(nnz_out, 6), a_in(nnz_in(1), 6))
!         a_in(:, :2) = data_B_ui
!         a_in(:, 3:) = data_W_ui
!         call erase_row_direction(nr_in(1), nnz_in(1), indi_ui, indj_ui, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         obj%nnz_u = nnz_out
!         obj%nr_u = nr_out
!         obj%nc_u = nc_ui
!         allocate(obj%indi_u(obj%nr_u+1), obj%indj_u(obj%nnz_u), obj%data_B_u(obj%nnz_u, 2), obj%data_W_u(obj%nnz_u, 4))
!         obj%indi_u = indi_out
!         obj%indj_u = indj_out
!         obj%data_B_u = a_out(:,:2)
!         obj%data_W_u = a_out(:,3:)

!         ! Erase data in direction v
!         ! -------------------------
!         nr_out = -1; nnz_out = -1
!         call erase_row_direction(nr_in(2), nnz_in(2), indi_vi, indj_vi, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         allocate(indi_out(nr_out+1), indj_out(nnz_out), a_out(nnz_out, 6), a_in(nnz_in(1), 6))
!         a_in(:, :2) = data_B_vi
!         a_in(:, 3:) = data_W_vi
!         call erase_row_direction(nr_in(2), nnz_in(2), indi_vi, indj_vi, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         obj%nnz_v = nnz_out
!         obj%nr_v = nr_out
!         obj%nc_v = nc_vi
!         allocate(obj%indi_v(obj%nr_v+1), obj%indj_v(obj%nnz_v), obj%data_B_v(obj%nnz_v, 2), obj%data_W_v(obj%nnz_v, 4))
!         obj%indi_v = indi_out
!         obj%indj_v = indj_out
!         obj%data_B_v = a_out(:,:2)
!         obj%data_W_v = a_out(:,3:)

!         ! Erase data in direction w
!         ! -------------------------
!         nr_out = -1; nnz_out = -1
!         call erase_row_direction(nr_in(3), nnz_in(3), indi_wi, indj_wi, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         allocate(indi_out(nr_out+1), indj_out(nnz_out), a_out(nnz_out, 6), a_in(nnz_in(3), 6))
!         a_in(:, :2) = data_B_wi
!         a_in(:, 3:) = data_W_wi
!         call erase_row_direction(nr_in(3), nnz_in(3), indi_wi, indj_wi, a_in, table, &
!                                 nr_out, nnz_out, indi_out, indj_out, a_out)
!         obj%nnz_w = nnz_out
!         obj%nr_w = nr_out
!         obj%nc_w = nc_wi
!         allocate(obj%indi_w(obj%nr_w+1), obj%indj_w(obj%nnz_w), obj%data_B_w(obj%nnz_w, 2), obj%data_W_w(obj%nnz_w, 4))
!         obj%indi_w = indi_out
!         obj%indj_w = indj_out
!         obj%data_B_w = a_out(:,:2)
!         obj%data_W_w = a_out(:,3:)

!     end subroutine paraprop_initialize

! end module wq_mf_pretreatment