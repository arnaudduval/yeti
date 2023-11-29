module matrixfreestheat

    implicit none
    type stthermomat
        integer :: dimen_sp
        double precision :: Cmean
        double precision, dimension(:), allocatable :: Kmean
        double precision, dimension(:), pointer :: Cprop=>null(), Cdersprop=>null(), detJ=>null(), detG=>null()
        double precision, dimension(:, :), pointer :: Kdersprop=>null()
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

        allocate(mat%Kmean(mat%dimen_sp))
        mat%Kmean = 1.d0; mat%Cmean = 1.d0

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

    subroutine setup_conductivityDersprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(mat%dimen_sp, nnz)

        mat%Kdersprop => prop
        if (nnz.ne.mat%ncols_total) stop 'Size problem'

    end subroutine setup_conductivityDersprop

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

    subroutine setup_capacityDersprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Cdersprop => prop
        if (nnz.ne.mat%ncols_total) stop 'Size problem'

    end subroutine setup_capacityDersprop

    subroutine compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
        
        use separatevariables
        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nc_list
        dimension :: nc_list(mat%dimen_sp+1)

        double precision, intent(out) :: univMcoefs(mat%dimen_sp+1, maxval(nc_list)), &
                                        univKcoefs(mat%dimen_sp+1, maxval(nc_list))

        ! Local data
        ! ----------
        type(sepoperator) :: oper
        integer :: i, j, k, l
        logical, allocatable, dimension(:) :: update
        double precision, allocatable, dimension(:, :) :: CC
        double precision :: tensor(mat%dimen_sp, mat%dimen_sp)

        allocate(CC(mat%dimen_sp+1, mat%ncols_total), update(mat%dimen_sp+1)); update = .true.
        call initialize_operator(oper, size(nc_list), nc_list, update)

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                tensor = matmul(mat%invJ(:, :, i), matmul(mat%Kprop(:, :, k), &
                                    transpose(mat%invJ(:, :, i))))*mat%detJ(i)*mat%detG(j)
                do l = 1, mat%dimen_sp
                    CC(l, k) = tensor(l, l)
                end do
                CC(mat%dimen_sp+1, k) = mat%Cprop(k)*mat%detJ(i)
            end do
        end do

        if (mat%dimen_sp.eq.2) then
            call separatevariables_3d(oper, CC)
        else if (mat%dimen_sp.eq.3) then
            call separatevariables_4d(oper, CC)
        end if
        univMcoefs = oper%univmasscoefs; univKcoefs = oper%univstiffcoefs
    
    end subroutine compute_separationvariables

    subroutine compute_mean(mat, nclist)
        !! Computes the average of the material properties

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nclist
        dimension :: nclist(mat%dimen_sp+1)

        ! Local data
        ! ----------
        integer, parameter :: NP = 3
        integer :: i, c_sp, c_tm, cp, gsp, gtm, gp, indlist_sp(mat%dimen_sp, NP), sample_sp(NP**mat%dimen_sp), sample_tm(NP)
        double precision, dimension(:, :), allocatable :: CC_K
        double precision, dimension(:), allocatable :: CC_M
        double precision :: tensor(mat%dimen_sp, mat%dimen_sp)
        
        do i = 1, mat%dimen_sp
            indlist_sp(i, :) = (/1, int((nclist(i) + 1)/2), nclist(i)/)
        end do
        call indices2list(mat%dimen_sp, NP, indlist_sp, nclist(:mat%dimen_sp), sample_sp)
        sample_tm = (/1, int((nclist(mat%dimen_sp+1) + 1)/2), nclist(mat%dimen_sp+1)/)
        
        allocate(CC_K(mat%dimen_sp, size(sample_sp)*size(sample_tm)))
        do c_tm = 1, size(sample_tm)
            do c_sp = 1, size(sample_sp)
                gtm = sample_tm(c_tm)
                gsp = sample_sp(c_sp)
                gp  = gsp + (gtm - 1)*mat%ncols_sp
                cp  = c_sp + (c_tm - 1)*size(sample_sp)
                tensor = matmul(mat%invJ(:, :, gsp), matmul(mat%Kprop(:, :, gp), &
                                transpose(mat%invJ(:, :, gsp))))*mat%detJ(gsp)*mat%detG(gtm)
                do i = 1, mat%dimen_sp
                    CC_K(i, cp) = tensor(i, i)
                end do
            end do
        end do

        do i = 1, mat%dimen_sp
            if (mat%dimen_sp.eq.2) then
                call trapezoidal_rule_3d(NP, NP, NP, CC_K(i, :), mat%Kmean(i))
            else if (mat%dimen_sp.eq.3) then
                call trapezoidal_rule_4d(NP, NP, NP, NP, CC_K(i, :), mat%Kmean(i))
            end if
        end do

        allocate(CC_M(size(sample_sp)*size(sample_tm)))
        do c_tm = 1, size(sample_tm)
            do c_sp = 1, size(sample_sp)
                gtm = sample_tm(c_tm)
                gsp = sample_sp(c_sp)
                gp  = gsp + (gtm - 1)*mat%ncols_sp
                cp  = c_sp + (c_tm - 1)*size(sample_sp)
                CC_M(cp) = mat%Cprop(gp)*mat%detJ(gsp)
            end do
        end do

        if (mat%dimen_sp.eq.2) then
            call trapezoidal_rule_3d(NP, NP, NP, CC_M, mat%Cmean)
        else if (mat%dimen_sp.eq.3) then
            call trapezoidal_rule_4d(NP, NP, NP, NP, CC_M, mat%Cmean)
        end if

    end subroutine compute_mean

    subroutine mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
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
                tmp(k) = tmp(k)*mat%Cdersprop(k)*mat%detJ(i)*mat%detG(j)
            end do
        end do

        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                            tmp, array_out)
        
    end subroutine mf_u_v_3d

    subroutine mf_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
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
                tmp(k) = tmp(k)*mat%Cdersprop(k)*mat%detJ(i)*mat%detG(j)
            end do
        end do

        call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                            tmp, array_out)
        
    end subroutine mf_u_v_4d

    subroutine mf_u_partialt_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
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
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 2), & 
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
                            nnz_t, indi_t, indj_t, data_W_t(:, 2), &
                            tmp, array_out)
        
    end subroutine mf_u_partialt_v_3d

    subroutine mf_u_partialt_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
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
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 2), & 
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
                            nnz_t, indi_t, indj_t, data_W_t(:, 2), &
                            tmp, array_out)
        
    end subroutine mf_u_partialt_v_4d

    subroutine mf_gradx_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
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
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, nc_total)
        integer :: i, j, k, alpha, zeta
        dimension :: alpha(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_t) stop 'Size problem'

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                coefs(:, k) = matmul(mat%invJ(:, :, i), mat%Kdersprop(:, k))*mat%detJ(i)*mat%detG(j)
            end do
        end do

        array_out = 0.d0
        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), & 
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                            array_in, tmp_0)
        do i = 1, dimen
            alpha = 1; alpha(i) = 2
            zeta  = 1 + (alpha - 1)*2
            tmp_1 = tmp_0*coefs(i, :)
            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                tmp_1, tmp_2)
            array_out = array_out + tmp_2
        end do

    end subroutine mf_gradx_u_v_3d

    subroutine mf_gradx_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
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
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, nc_total)
        integer :: i, j, k, alpha, zeta
        dimension :: alpha(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_w*nr_t) stop 'Size problem'

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                coefs(:, k) = matmul(mat%invJ(:, :, i), mat%Kdersprop(:, k))*mat%detJ(i)*mat%detG(j)
            end do
        end do

        array_out = 0.d0
        call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), & 
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                            array_in, tmp_0)
        do i = 1, dimen
            alpha = 1; alpha(i) = 2
            zeta  = 1 + (alpha - 1)*2
            tmp_1 = tmp_0*coefs(i, :)
            call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                tmp_1, tmp_2)
            array_out = array_out + tmp_2
        end do

    end subroutine mf_gradx_u_v_4d

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
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, & 
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

module stheatsolver2

    use matrixfreestheat
    use datastructure
    type stcgsolver
        integer :: matrixfreetype = 1, dimen = 3
        double precision :: Cmean = 1.d0
        type(structure) :: temp_struct
    end type stcgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_t, indj_T_u, indj_T_v, indj_T_t
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
        double precision :: tmp(nr_total)

        call mf_u_partialt_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                    nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                    indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                                    indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                                    data_W_t, array_in, array_out)

        call mf_gradx_u_gradx_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_t, array_in, tmp)

        array_out = array_out + tmp

        if (solv%matrixfreetype.eq.1) then
            return
        else if (solv%matrixfreetype.eq.2) then
            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_t, array_in, tmp)
            array_out = array_out + tmp
            call mf_gradx_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_W_u, data_W_v, &
                            data_W_t, array_in, tmp)
            array_out = array_out + tmp
        else
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, &
                data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, data_W_t, table, Kmean)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_u, indi_v, indi_t, indj_u, indj_v, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, data_W_t
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_t(nnz_t, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2)
        double precision, intent(in) :: Kmean
        dimension :: Kmean(solv%dimen)

        call init_3datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                            indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, data_W_t)
        call update_datastructure(solv%temp_struct, solv%dimen, table)
        solv%temp_struct%isspacetime = .true.
        call space_eigendecomposition(solv%temp_struct, solv%dimen, Kmean)
        call time_schurdecomposition(solv%temp_struct)
    
    end subroutine initializefastdiag

    subroutine solve_schurtriangular__(solv, nr, coefs, b, x)
        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr
        double precision, intent(in) :: coefs(2)
        double precision, intent(in) :: b(nr)
        double precision, intent(out) :: x(nr)
        
        ! Local data
        ! ----------
        double complex :: uptrg(nr, nr), Vb(nr), sol(nr), Usol(nr)

        Vb = matmul(transpose(conjg(solv%temp_struct%Lschur_tm)), dcmplx(b))
        uptrg = dcmplx(coefs(1))*solv%temp_struct%Luptr_tm + dcmplx(coefs(2))*solv%temp_struct%Ruptr_tm
        call solve_complex_uppertriangular_system(nr, uptrg, Vb, sol)
        Usol = matmul(solv%temp_struct%Rschur_tm, sol)
        x = realpart(Usol)

    end subroutine solve_schurtriangular__

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_t, i, j, l
        double precision, allocatable, dimension(:, :) :: identity
        double precision, allocatable, dimension(:) :: tmp, tmp1, tmp2, tmp3, btmp, sol
        integer, pointer, dimension(:, :, :) :: dof
        integer, allocatable, dimension(:) :: indices 
        double precision :: eigval
        integer :: nrows(solv%dimen)

        array_out = 0.d0

        nrows = solv%temp_struct%nrows
        nr_u = nrows(1); nr_v = nrows(2); nr_t = nrows(3)
        allocate(dof(nr_u, nr_v, nr_t))
        dof  = reshape(solv%temp_struct%dof, shape=(/nr_u, nr_v, nr_t/))

        ! Compute (Ut x Uv x Uu)'.array_in
        allocate(tmp(nr_u*nr_v*nr_t), identity(nr_t, nr_t))
        call create_identity(nr_t, identity)

        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_t, nr_t, &
                        transpose(solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        identity, array_in(solv%temp_struct%dof), tmp)

        allocate(tmp1(nr_total))
        tmp1 = 0.d0; tmp1(solv%temp_struct%dof) = tmp
        deallocate(tmp)

        allocate(btmp(nr_t), sol(nr_t), indices(nr_t), tmp2(nr_total))
        tmp2 = 0.d0
        do j = 1, nr_v
            do i = 1, nr_u
                l = i + (j - 1)*nr_u
                eigval = solv%temp_struct%diageigval_sp(l)
                indices = dof(i, j, :)
                btmp = tmp1(indices)
                call solve_schurtriangular__(solv, nr_t, (/solv%Cmean, eigval/), btmp, sol)
                tmp2(indices) = sol
            end do
        end do
        deallocate(tmp1)

        ! Compute (Ut x Uv x Uu).array_tmp
        allocate(tmp3(nr_u*nr_v*nr_t))
        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_t, nr_t, &
                        solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        identity, tmp2(solv%temp_struct%dof), tmp3)
        array_out(solv%temp_struct%dof) = tmp3
        deallocate(tmp2, tmp3)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(nnz)

        call set2zero(solv%temp_struct, nnz, array)

    end subroutine clear_dirichlet

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), &
                        s(nr_total), Ap(nr_total), As(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call clear_dirichlet(solv, nr_total, r) 
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, p, Ap)

            call clear_dirichlet(solv, nr_total, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, s, As)
            call clear_dirichlet(solv, nr_total, As)
            omega = dot_product(As, s)/dot_product(As, As)
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_t(nnz_t, 4)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
                        ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_W_u, data_W_v, data_W_t, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
    
            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module stheatsolver2

module stheatsolver3

    use matrixfreestheat
    use datastructure
    type stcgsolver
        integer :: matrixfreetype=1, dimen = 4
        double precision :: Cmean = 1.d0
        type(structure) :: temp_struct
    end type stcgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t, indj_T_u, indj_T_v, indj_T_w, indj_T_t
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
        double precision :: tmp(nr_total)

        call mf_u_partialt_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                    nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                    indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                    indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                    data_W_w, data_W_t, array_in, array_out)

        call mf_gradx_u_gradx_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_w, data_W_t, array_in, tmp)

        array_out = array_out + tmp

        if (solv%matrixfreetype.eq.1) then
            return
        else if (solv%matrixfreetype.eq.2) then
            call mf_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_w, data_W_t, array_in, tmp)

            array_out = array_out + tmp
            call mf_gradx_u_v_4d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                indi_T_w, indj_T_w, indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, &
                                data_W_w, data_W_t, array_in, tmp)

            array_out = array_out + tmp
        else
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, &
                data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t, table)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), data_B_t(nnz_t, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2)
        
        ! Local data
        ! ----------
        double precision :: dummymean(solv%dimen)
        call init_4datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t)
        call update_datastructure(solv%temp_struct, solv%dimen, table)
        solv%temp_struct%isspacetime = .true.; dummymean = 1.d0
        call space_eigendecomposition(solv%temp_struct, solv%dimen, dummymean)
        call time_schurdecomposition(solv%temp_struct)

    end subroutine initializefastdiag

    subroutine solve_schurtriangular__(solv, nr, coefs, b, x)
        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr
        double precision, intent(in) :: coefs(2)
        double precision, intent(in) :: b(nr)
        double precision, intent(out) :: x(nr)
        
        ! Local data
        ! ----------
        double complex :: uptrg(nr, nr), Vb(nr), sol(nr), Usol(nr)

        Vb = matmul(transpose(conjg(solv%temp_struct%Lschur_tm)), dcmplx(b))
        uptrg = dcmplx(coefs(1))*solv%temp_struct%Luptr_tm + dcmplx(coefs(2))*solv%temp_struct%Ruptr_tm
        call solve_complex_uppertriangular_system(nr, uptrg, Vb, sol)
        Usol = matmul(solv%temp_struct%Rschur_tm, sol)
        x = realpart(Usol)

    end subroutine solve_schurtriangular__

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w, nr_t, i, j, k, l
        double precision, allocatable, dimension(:, :) :: identity
        double precision, allocatable, dimension(:) :: tmp, tmp1, tmp2, tmp3, btmp, sol
        integer, allocatable, dimension(:, :, :, :) :: dof
        integer, allocatable, dimension(:) :: indices 
        double precision :: eigval
        integer :: nrows(solv%dimen)

        array_out = 0.d0

        nrows = solv%temp_struct%nrows
        nr_u = nrows(1); nr_v = nrows(2); nr_w = nrows(3); nr_t = nrows(4)
        allocate(dof(nr_u, nr_v, nr_w, nr_t))
        dof  = reshape(solv%temp_struct%dof, shape=(/nr_u, nr_v, nr_w, nr_t/))

        ! Compute (Ut x Uw x Uv x Uu)'.array_in
        allocate(tmp(nr_u*nr_v*nr_w*nr_t), identity(nr_t, nr_t))
        call create_identity(nr_t, identity)

        call sumfacto4d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        transpose(solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                        identity, array_in(solv%temp_struct%dof), tmp)

        allocate(tmp1(nr_total))
        tmp1 = 0.d0; tmp1(solv%temp_struct%dof) = tmp
        deallocate(tmp)

        allocate(btmp(nr_t), sol(nr_t), indices(nr_t), tmp2(nr_total))
        tmp2 = 0.d0
        do k = 1, nr_w
            do j = 1, nr_v
                do i = 1, nr_u
                    l = i + (j - 1)*nr_u + (k - 1)*nr_u*nr_v
                    eigval = solv%temp_struct%diageigval_sp(l)
                    indices = dof(i, j, k, :)
                    btmp = tmp1(indices)
                    call solve_schurtriangular__(solv, nr_t, (/solv%Cmean, eigval/), btmp, sol)
                    tmp2(indices) = sol
                end do
            end do
        end do  
        deallocate(tmp1)

        ! Compute (Ut x Uw x Uv x Uu).array_tmp
        allocate(tmp3(nr_u*nr_v*nr_w*nr_t))
        call sumfacto4d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), &
                        identity, tmp2(solv%temp_struct%dof), tmp3)
        array_out(solv%temp_struct%dof) = tmp3
        deallocate(tmp2, tmp3)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(nnz)

        call set2zero(solv%temp_struct, nnz, array)

    end subroutine clear_dirichlet

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_w, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, Ap, As
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), &
                        s(nr_total), Ap(nr_total), As(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call clear_dirichlet(solv, nr_total, r) 
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return

        do iter = 1, nbIterPCG
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, p, Ap)

            call clear_dirichlet(solv, nr_total, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, s, As)
            call clear_dirichlet(solv, nr_total, As)
            omega = dot_product(As, s)/dot_product(As, As)
            x = x + alpha*p + omega*s
            r = s - omega*As

            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Ap)
            rsold = rsnew
        end do

    end subroutine BiCGSTAB

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indi_T_t
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), indi_T_t(nc_t+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w, indj_T_t
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w), indj_T_t(nnz_t)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w, data_BT_t
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2), data_BT_t(nnz_t, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
        integer, intent(in) :: indj_u, indj_v, indj_w, indj_t
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w, data_W_t
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4), data_W_t(nnz_t, 4)

        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
                        ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
        integer :: iter

        x = 0.d0; r = b
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return
        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            resPCG(iter+1) = norm2(r)/normb
            if (resPCG(iter+1).le.threshold) exit
    
            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module stheatsolver3