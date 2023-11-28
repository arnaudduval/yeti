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


module heatsolver2

    use matrixfreeheat
    use datastructure
    type cgsolver
        logical :: withdiag = .true.
        integer :: matrixfreetype = 1, dimen = 2
        type(structure) :: temp_struct
    end type cgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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

        if (solv%matrixfreetype.eq.1) then

            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        else 
            call mf_uv_gradugradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, &
                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                data_B_u, data_B_v, data_W_u, data_W_v, table_dirichlet, Kmean)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_B_u, data_B_v, data_W_u, data_W_v
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
        logical, intent(in) :: table_dirichlet
        dimension :: table_dirichlet(solv%dimen, 2)
        double precision, intent(in) :: Kmean
        dimension :: Kmean(solv%dimen)

        call init_2datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v)
        call update_datastructure(solv%temp_struct, solv%dimen, table_dirichlet)
        call space_eigendecomposition(solv%temp_struct, solv%dimen, Kmean)
    
    end subroutine initializefastdiag

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v
        double precision, allocatable, dimension(:) :: tmp, tmp2

        ! Compute (Uw x Uv x Uu)'.array_in
        nr_u = solv%temp_struct%nrows(1)
        nr_v = solv%temp_struct%nrows(2)
        allocate(tmp(nr_u*nr_v))
        call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, transpose(solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                transpose(solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), array_in(solv%temp_struct%dof), tmp)

        if (solv%withdiag) then
            tmp = tmp/solv%temp_struct%diageigval_sp
        end if
    
        ! Compute (Uv x Uu).array_tmp
        allocate(tmp2(nr_u*nr_v))
        call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), tmp, tmp2)
        array_out(solv%temp_struct%dof) = tmp2            
        deallocate(tmp, tmp2)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(nnz)

        call set2zero(solv%temp_struct, nnz, array)

    end subroutine clear_dirichlet

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, p, Ap)
            call clear_dirichlet(solv, nr_total, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                    nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                    data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                    data_W_u, data_W_v, s, As)
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

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)

            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, stilde, Astilde)
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

    subroutine LOBPCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ishigher, nbIterPCG, threshold, eigenvec, eigenval)
        !! Using LOBPCG algorithm to compute the stability of the transient heat problem
        
        implicit none
        ! Input / output data
        ! -------------------
        integer, parameter :: d = 3
        type(cgsolver) :: solv
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

        logical, intent(in) :: ishigher
        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(nr_total)

        ! Local data
        ! ----------
        integer :: k, ii
        double precision, dimension(d, nr_total) :: RM1, RM2, RM3
        double precision, dimension(d, d) :: AA1, BB1
        double precision, dimension(d) :: delta
        double precision, dimension(nr_total) :: u, v, g, gtil, p, tmp
        double precision :: q, norm
        double precision, allocatable, dimension(:) :: ll
        double precision, allocatable, dimension(:, :) :: qq
    
        call random_number(eigenvec)
        call clear_dirichlet(solv, nr_total, eigenvec)
        norm = norm2(eigenvec)
        eigenvec = eigenvec/norm

        call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, u)
        call clear_dirichlet(solv, nr_total, u)

        q = sqrt(dot_product(eigenvec, u))
        eigenvec = eigenvec/q; u = u/q
        call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, v)
        call clear_dirichlet(solv, nr_total, v)
        eigenval = dot_product(eigenvec, v)
        p = 0.d0
        norm = 1.d0

        do k = 1, nbIterPCG
            if (norm.le.threshold) return
    
            g = v - eigenval*u
            norm = norm2(g)
            call applyfastdiag(solv, nr_total, g, gtil)
            call clear_dirichlet(solv, nr_total, gtil)
            g = gtil

            RM1(1, :) = eigenvec; RM1(2, :) = -g; RM1(3, :) = p
            RM2(1, :) = v; RM3(1, :) = u;
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v,  -g, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(2, :) = tmp
            
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(3, :) = tmp
                        
            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, -g, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(2, :) = tmp

            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(3, :) = tmp
            
            call rayleigh_submatrix(d, nr_total, RM1, RM2, AA1); AA1 = 0.5d0*(AA1 + transpose(AA1))
            call rayleigh_submatrix(d, nr_total, RM1, RM3, BB1); BB1 = 0.5d0*(BB1 + transpose(BB1))

            if (k.eq.1) then
                allocate(ll(d-1), qq(d-1, d-1))
                call compute_eigdecomp_pdr(size(ll), AA1(:d-1, :d-1), BB1(:d-1, :d-1), ll, qq)
            else
                allocate(ll(d), qq(d, d))
                call compute_eigdecomp_pdr(size(ll), AA1, BB1, ll, qq)
            end if
            
            if (ishigher) then 
                eigenval = maxval(ll); ii = maxloc(ll, dim=1)
            else
                eigenval = minval(ll); ii = minloc(ll, dim=1)
            end if

            delta = 0.d0
            if (k.eq.1) then
                delta(:2) = qq(:, ii)
            else
                delta = qq(:, ii)
            end if
    
            p = -g*delta(2) + p*delta(3)
            eigenvec = eigenvec*delta(1) + p
            call mf_u_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v,&
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, eigenvec, u)
            call clear_dirichlet(solv, nr_total, u)
            
            q = sqrt(dot_product(eigenvec, u))
            eigenvec = eigenvec/q; u = u/q
            call mf_gradu_gradv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, eigenvec, v)
            call clear_dirichlet(solv, nr_total, v)
            
            norm = norm2(g)
            deallocate(ll, qq)
        end do

    end subroutine LOBPCGSTAB

end module heatsolver2

module heatsolver3

    use matrixfreeheat
    use datastructure
    type cgsolver
        logical :: withdiag = .true.
        integer :: matrixfreetype = 1, dimen = 3
        type(structure) :: temp_struct
    end type cgsolver

contains

    subroutine matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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

        if (solv%matrixfreetype.eq.1) then

            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.3) then

            call mf_uv_gradugradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else 
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, Kmean)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2), &
                    data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2)
        double precision, intent(in) :: Kmean
        dimension :: Kmean(solv%dimen)

        call init_3datastructure(solv%temp_struct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)
        call update_datastructure(solv%temp_struct, solv%dimen, table)
        call space_eigendecomposition(solv%temp_struct, solv%dimen, Kmean)
    
    end subroutine initializefastdiag

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w
        double precision, allocatable, dimension(:) :: tmp, tmp2

        array_out = 0.d0

        ! Compute (Uw x Uv x Uu)'.array_in
        nr_u = solv%temp_struct%nrows(1)
        nr_v = solv%temp_struct%nrows(2)
        nr_w = solv%temp_struct%nrows(3)
        allocate(tmp(nr_u*nr_v*nr_w))

        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        transpose(solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        transpose(solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
        array_in(solv%temp_struct%dof), tmp)

        tmp = tmp/solv%temp_struct%diageigval_sp

        ! Compute (Uw x Uv x Uu).array_tmp
        allocate(tmp2(nr_u*nr_v*nr_w))
        call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        solv%temp_struct%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%temp_struct%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        solv%temp_struct%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), tmp, tmp2)
        array_out(solv%temp_struct%dof) = tmp2
        deallocate(tmp, tmp2)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(nnz)

        call set2zero(solv%temp_struct, nnz, array)

    end subroutine clear_dirichlet

    subroutine BiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Ap)
            call clear_dirichlet(solv, nr_total, Ap)
            alpha = rsold/dot_product(Ap, rhat)
            s = r - alpha*Ap

            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, As)
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

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
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
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
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

    subroutine LOBPCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ishigher, nbIterPCG, threshold, eigenvec, eigenval)
        !! Using LOBPCG algorithm to compute the stability of the transient heat problem
        
        implicit none
        ! Input / output data
        ! -------------------
        integer, parameter :: d = 3
        type(cgsolver) :: solv
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
        
        logical, intent(in) :: ishigher
        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(nr_total)

        ! Local data
        ! ----------
        integer :: k, ii
        double precision, dimension(d, nr_total) :: RM1, RM2, RM3
        double precision, dimension(d, d) :: AA1, BB1
        double precision, dimension(d) :: delta
        double precision, dimension(nr_total) :: u, v, g, gtil, p, tmp
        double precision :: q, norm
        double precision, allocatable, dimension(:) ::  ll
        double precision, allocatable, dimension(:, :) ::  qq

        call random_number(eigenvec)
        call clear_dirichlet(solv, nr_total, eigenvec)
        norm = norm2(eigenvec)
        eigenvec = eigenvec/norm

        call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, u)
        call clear_dirichlet(solv, nr_total, u)
        
        q = sqrt(dot_product(eigenvec, u))
        eigenvec = eigenvec/q; u = u/q
        call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, v)
        call clear_dirichlet(solv, nr_total, v)
        
        eigenval = dot_product(eigenvec, v)
        p = 0.d0
        norm = 1.d0

        do k = 1, nbIterPCG
            if (norm.le.threshold) return
    
            g = v - eigenval*u
            norm = norm2(g)
            call applyfastdiag(solv, nr_total, g, gtil)
            call clear_dirichlet(solv, nr_total, gtil)
            g = gtil

            RM1(1, :) = eigenvec; RM1(2, :) = -g; RM1(3, :) = p
            RM2(1, :) = v; RM3(1, :) = u;
            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, -g, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(2, :) = tmp

            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(3, :) = tmp

            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, -g, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(2, :) = tmp
            
            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(3, :) = tmp
            
            call rayleigh_submatrix(d, nr_total, RM1, RM2, AA1); AA1 = 0.5d0*(AA1 + transpose(AA1))
            call rayleigh_submatrix(d, nr_total, RM1, RM3, BB1); BB1 = 0.5d0*(BB1 + transpose(BB1))

            if (k.eq.1) then
                allocate(ll(d-1), qq(d-1, d-1))
                call compute_eigdecomp_pdr(size(ll), AA1(:d-1, :d-1), BB1(:d-1, :d-1), ll, qq)
            else
                allocate(ll(d), qq(d, d))
                call compute_eigdecomp_pdr(size(ll), AA1, BB1, ll, qq)
            end if
            
            if (ishigher) then 
                eigenval = maxval(ll); ii = maxloc(ll, dim=1)
            else
                eigenval = minval(ll); ii = minloc(ll, dim=1)
            end if

            delta = 0.d0
            if (k.eq.1) then
                delta(:2) = qq(:, ii)
            else
                delta = qq(:, ii)
            end if
    
            p = -g*delta(2) + p*delta(3)
            eigenvec = eigenvec*delta(1) + p
            call mf_u_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, eigenvec, u)
            call clear_dirichlet(solv, nr_total, u)
            
            q = sqrt(dot_product(eigenvec, u))
            eigenvec = eigenvec/q; u = u/q
            call mf_gradu_gradv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, eigenvec, v)
            call clear_dirichlet(solv, nr_total, v)
            
            norm = norm2(g)
            deallocate(ll, qq)
        end do

    end subroutine LOBPCGSTAB

end module heatsolver3