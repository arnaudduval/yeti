module matrixfreeheat

    implicit none
    type thermomat
        integer :: dimen
        logical :: isLumped = .false.
        double precision :: Cmean, scalars(2) = (/1.d0, 1.d0/)
        double precision, dimension(:), allocatable :: Kmean
        double precision, dimension(:), pointer :: Hprop=>null(), Cprop=>null(), detJ=>null()
        double precision, dimension(:, :, :), pointer :: Kprop=>null(), invJ=>null()
        integer :: ncols_sp
    end type thermomat

contains

    subroutine setup_geometry(mat, dimen, nnz, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 

        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: dimen, nnz
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(dimen, dimen, nnz), detJ(nnz)

        mat%dimen = dimen
        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nnz
        allocate(mat%Kmean(mat%dimen))
        mat%Kmean = 1.d0; mat%Cmean = 1.d0

    end subroutine setup_geometry

    subroutine setup_conductivityprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(mat%dimen, mat%dimen, nnz)

        mat%Kprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_conductivityprop

    subroutine setup_capacityprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Cprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_capacityprop

    subroutine setup_thmchcoupledprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Hprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_thmchcoupledprop

    subroutine compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
        
        use separatevariables
        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nc_list
        dimension :: nc_list(mat%dimen)

        double precision, intent(out) :: univMcoefs(mat%dimen, maxval(nc_list)), univKcoefs(mat%dimen, maxval(nc_list))

        ! Local data
        ! ----------
        type(sepoperator) :: oper
        integer :: i, gp
        integer, allocatable, dimension(:) :: nc_list_t
        logical, allocatable, dimension(:) :: update
        double precision, allocatable, dimension(:, :) :: CC
        double precision :: tensor(mat%dimen, mat%dimen)

        if (.not.associated(mat%Cprop)) then
            allocate(CC(mat%dimen, mat%ncols_sp), update(mat%dimen), nc_list_t(mat%dimen))
            update = .true.; nc_list_t = nc_list
            call initialize_operator(oper, mat%dimen, nc_list_t, update)
            
            do gp = 1, mat%ncols_sp
                tensor = matmul(mat%invJ(:, :, gp), matmul(mat%Kprop(:, :, gp), &
                        transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
                do i = 1, mat%dimen
                    CC(i, gp) = tensor(i, i)
                end do
            end do

            if (mat%dimen.eq.2) then
                call separatevariables_2d(oper, CC)
            else if (mat%dimen.eq.3) then
                call separatevariables_3d(oper, CC)
            end if

            univMcoefs = oper%univmasscoefs; univKcoefs = oper%univstiffcoefs

        else
            allocate(CC(mat%dimen+1, mat%ncols_sp), update(mat%dimen+1), nc_list_t(mat%dimen+1))
            update = .true.; update(mat%dimen+1) = .false.
            nc_list_t(:mat%dimen) = nc_list; nc_list_t(mat%dimen+1) = 1
            call initialize_operator(oper, mat%dimen+1, nc_list_t, update)

            do gp = 1, mat%ncols_sp
                tensor = matmul(mat%invJ(:, :, gp), matmul(mat%Kprop(:, :, gp), &
                        transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
                do i = 1, mat%dimen
                    CC(i, gp) = tensor(i, i)
                end do
                CC(mat%dimen+1, gp) = mat%Cprop(gp)*mat%detJ(gp)
            end do

            if (mat%dimen.eq.2) then
                call separatevariables_3d(oper, CC)
            else if (mat%dimen.eq.3) then
                call separatevariables_4d(oper, CC)
            end if
            univMcoefs = oper%univmasscoefs(:mat%dimen, :); univKcoefs = oper%univstiffcoefs(:mat%dimen, :)
        end if
        
    end subroutine compute_separationvariables

    subroutine compute_mean(mat, nclist)
        !! Computes the average of the material properties

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nclist
        dimension :: nclist(mat%dimen)

        ! Local data
        ! ----------
        integer, parameter :: NP = 3
        integer :: i, c, gp, sample(NP**mat%dimen), indlist(mat%dimen, NP)
        double precision, dimension(:, :), allocatable :: CC_K
        double precision, dimension(:), allocatable :: CC_M
        double precision :: tensor(mat%dimen, mat%dimen)
        
        ! Select sample
        do i = 1, mat%dimen 
            indlist(i, :) =  (/1, int((nclist(i) + 1)/2), nclist(i)/)
        end do

        call indices2list(mat%dimen, NP, indlist, nclist, sample)
        
        allocate(CC_K(mat%dimen, size(sample)))
        do c = 1, size(sample)
            gp = sample(c)
            tensor = matmul(mat%invJ(:, :, gp), matmul(mat%Kprop(:, :, gp), &
                            transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
            do i = 1, mat%dimen
                CC_K(i, c) = tensor(i, i)
            end do
        end do
        do i = 1, mat%dimen
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(NP, NP, CC_K(i, :), mat%Kmean(i))
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(NP, NP, NP, CC_K(i, :), mat%Kmean(i))
            end if
        end do

        if (associated(mat%Cprop)) then
            allocate(CC_M(size(sample)))
            do c = 1, size(sample)
                gp = sample(c)
                CC_M(c) = mat%Cprop(gp)*mat%detJ(gp)
            end do
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(NP, NP, CC_M, mat%Cmean)
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(NP, NP, NP, CC_M, mat%Cmean)
            end if
        end if   

    end subroutine compute_mean

    subroutine mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_T_u, indi_T_v, indj_T_u, indj_T_v
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        double precision :: tmp_in, array_tmp
        dimension :: tmp_in(nr_total), array_tmp(nc_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        tmp_in = array_in; if (mat%isLumped) tmp_in = 1.d0

        call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            tmp_in, array_tmp)

        array_tmp = array_tmp*mat%Cprop*mat%detJ

        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            array_tmp, array_out)
        if (mat%isLumped) array_out = array_out*array_in
        
    end subroutine mf_u_v_2d

    subroutine mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        double precision :: tmp_in, tmp
        dimension :: tmp_in(nr_total), tmp(nc_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        tmp_in = array_in; if (mat%isLumped) tmp_in = 1.d0

        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                            tmp_in, tmp)

        tmp = tmp*mat%Cprop*mat%detJ

        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            tmp, array_out)
        if (mat%isLumped) array_out = array_out*array_in
        
    end subroutine mf_u_v_3d

    subroutine mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        double precision :: tmp_0, tmp_1, tmp_2, coefs
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, dimen, nc_total)
        integer :: i, j, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        do i = 1, nc_total
            coefs(:, :, i) = matmul(mat%invJ(:, :, i), matmul(mat%Kprop(:, :, i), &
                            transpose(mat%invJ(:, :, i))))*mat%detJ(i)
        end do

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                array_in, tmp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                tmp_1 = tmp_0*coefs(i, j, :)
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    tmp_1, tmp_2)
                array_out = array_out + tmp_2
            end do
        end do
        
    end subroutine mf_gradu_gradv_2d

    subroutine mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        double precision :: tmp_0, tmp_1, tmp_2, coefs
        dimension :: tmp_0(nc_total), tmp_1(nc_total), tmp_2(nr_total), coefs(dimen, dimen, nc_total)
        integer :: i, j, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'

        do i = 1, nc_total
            coefs(:, :, i) = matmul(mat%invJ(:, :, i), matmul(mat%Kprop(:, :, i), &
                                transpose(mat%invJ(:, :, i))))*mat%detJ(i)
        end do

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                                array_in, tmp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                tmp_1 = tmp_0*coefs(i, j, :)
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                    tmp_1, tmp_2)
                array_out = array_out + tmp_2
            end do
        end do
        
    end subroutine mf_gradu_gradv_3d

    subroutine mf_uv_gradugradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(nr_total)

        call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_tmp)

        array_out = array_out + mat%scalars(2)*array_tmp
        
    end subroutine mf_uv_gradugradv_2d

    subroutine mf_uv_gradugradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(nr_total)

        call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_tmp)

        array_out = array_out + mat%scalars(2)*array_tmp
        
    end subroutine mf_uv_gradugradv_3d

    subroutine mf_gradu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        integer :: i, k, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        array_out = 0.d0  
        do i = 1, dimen
            beta = 1
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in(i, :), t1) 
            t1 = t1*mat%Hprop*mat%detJ
            do k = 1, dimen
                alpha = 1; alpha(k) = 2; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(k, i, :)
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), t2, t3)
                array_out = array_out + t3
            end do
        end do
            
    end subroutine mf_gradu_tv_2d

    subroutine mf_gradu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        integer :: i, k, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        array_out = 0.d0  
        do i = 1, dimen
            beta = 1
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in(i, :), t1)  
            t1 = t1*mat%Hprop*mat%detJ
            do k = 1, dimen
                alpha = 1; alpha(k) = 2; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(k, i, :)
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                            nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t2, t3)
                array_out = array_out + t3
            end do
        end do
                    
    end subroutine mf_gradu_tv_3d

end module matrixfreeheat
