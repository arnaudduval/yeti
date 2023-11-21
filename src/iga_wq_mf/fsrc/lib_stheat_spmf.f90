module matrixfreestheat

    implicit none
    type stthermomat
        integer :: dimen_sp
        double precision, dimension(:), pointer :: Cprop=>null(), detJ=>null(), detG=>null()
        double precision, dimension(:, :, :), pointer :: Kprop=>null(), invJ=>null()
        integer :: ncols_sp, ncols_tm, ncols_total
    end type stthermomat

contains

    subroutine setup_geometry(mat, dimen_sp, nnz_sp, nnz_tm, invJ, detJ, detG)
        !! Points to the data of the inverse and determinant of the Jacobian. 

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: dimen_sp, nnz_sp, nnz_tm
        double precision, target, intent(in) :: invJ, detJ, detG
        dimension :: invJ(dimen_sp, dimen_sp, nnz_sp), detJ(nnz_sp), detG(nnz_tm)

        mat%dimen_sp = dimen_sp
        mat%invJ => invJ
        mat%detJ => detJ
        mat%detG => detG
        mat%ncols_sp = nnz_sp
        mat%ncols_tm = nnz_tm
        mat%ncols_total = nnz_sp*nnz_tm

    end subroutine setup_geometry

    subroutine setup_conductivityprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(mat%dimen_sp, mat%dimen_sp, nnz)

        mat%Kprop => prop
        if (nnz.ne.mat%ncols_total) stop 'Size problem'

    end subroutine setup_conductivityprop

    subroutine setup_capacityprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Cprop => prop
        if (nnz.ne.mat%ncols_total) stop 'Size problem'

    end subroutine setup_capacityprop

    subroutine mf_partialt_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_t, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indj_T_u, indj_T_v, indi_T_t, indj_T_t
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_t(nc_t+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_t, indj_u, indj_v, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        double precision :: tmp
        dimension :: tmp(nc_total)
        integer :: i, j, k

        if (nr_total.ne.nr_u*nr_v*nr_t) stop 'Size problem'

        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                            array_in, tmp)

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                tmp(k) = tmp(k)*mat%Cprop(k)*mat%detJ(i)
            end do
        end do

        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 3), &
                            tmp, array_out)
        
    end subroutine mf_partialt_u_v_3d

    subroutine mf_partialt_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_w, data_W_t, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w, indi_T_t, indj_T_t
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        double precision :: tmp
        dimension :: tmp(nc_total)
        integer :: i, j, k

        if (nr_total.ne.nr_u*nr_v*nr_w*nr_t) stop 'Size problem'

        call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                            array_in, tmp)

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                tmp(k) = tmp(k)*mat%Cprop(k)*mat%detJ(i)
            end do
        end do

        call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 3), &
                            tmp, array_out)
        
    end subroutine mf_partialt_u_v_4d

    subroutine mf_gradx_u_gradx_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_t, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indj_T_u, indj_T_v, indi_T_t, indj_T_t
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_t(nc_t+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_t, indj_u, indj_v, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        double precision :: tmp_0, tmp_1, tmp_2, coefs
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, dimen, nc_total)
        integer :: i, j, k, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_t) stop 'Size problem'

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                coefs(:, :, k) = matmul(mat%invJ(:, :, i), matmul(mat%Kprop(:, :, k), &
                                    transpose(mat%invJ(:, :, i))))*mat%detJ(i)*mat%detG(j)
            end do
        end do

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                                array_in, tmp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                tmp_1 = tmp_0*coefs(i, j, :)
                call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                    tmp_1, tmp_2)
                array_out = array_out + tmp_2
            end do
        end do

    end subroutine mf_gradx_u_gradx_v_3d

    subroutine mf_gradx_u_gradx_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_w, data_W_t, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w, indi_T_t, indj_T_t
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        double precision :: tmp_0, tmp_1, tmp_2, coefs
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, dimen, nc_total)
        integer :: i, j, k, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_w*nr_t) stop 'Size problem'

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                coefs(:, :, k) = matmul(mat%invJ(:, :, i), matmul(mat%Kprop(:, :, k), &
                                    transpose(mat%invJ(:, :, i))))*mat%detJ(i)*mat%detG(j)
            end do
        end do

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                                nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                                array_in, tmp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                tmp_1 = tmp_0*coefs(i, j, :)
                call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                    nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                    tmp_1, tmp_2)
                array_out = array_out + tmp_2
            end do
        end do

    end subroutine mf_gradx_u_gradx_v_4d

end module matrixfreestheat
