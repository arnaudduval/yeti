module matrixfreeheat
    use omp_lib
    use structured_data
    implicit none
    type thermomat
        integer :: dimen, ncols_sp
        logical :: isLumped = .false.
        double precision :: scalars(2) = (/1.d0, 1.d0/)
        ! Material properties
        double precision :: Cmean
        double precision, dimension(:), allocatable :: Kmean
        double precision, dimension(:), pointer :: detJ=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:), allocatable :: Cprop, Hprop
        double precision, dimension(:, :, :), allocatable :: Kprop
    end type thermomat

contains

    subroutine setup_geometry(mat, dimen, nnz, invJ, detJ)
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

        ! Local data 
        ! ----------
        integer :: i!, nb_tasks

        if ((.not.associated(mat%invJ)).or.(.not.associated(mat%detJ))) stop 'Define geometry'
        if (nnz.ne.mat%ncols_sp) stop 'Size problem'
        allocate(mat%Kprop(mat%dimen, mat%dimen, nnz))

        !!!$OMP PARALLEL
        !!nb_tasks = omp_get_num_threads()
        !!!$OMP DO SCHEDULE(STATIC, mat%ncols_sp/nb_tasks) 
        do i = 1, mat%ncols_sp
            mat%Kprop(:, :, i) = matmul(mat%invJ(:, :, i), matmul(prop(:, :, i), &
                    transpose(mat%invJ(:, :, i))))*mat%detJ(i)
        end do
        !!!$OMP END DO
        !!!$OMP END PARALLEL
    end subroutine setup_conductivityprop

    subroutine setup_capacityprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)
        integer :: i!, nb_tasks

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_sp) stop 'Size problem'
        allocate(mat%Cprop(nnz))

        !!!$OMP PARALLEL
        !!nb_tasks = omp_get_num_threads()
        !!!$OMP DO SCHEDULE(STATIC, nnz/nb_tasks)
        do i = 1, nnz
            mat%Cprop(i) = prop(i)*mat%detJ(i)
        end do
        !!!$OMP END DO
        !!!$OMP END PARALLEL
    end subroutine setup_capacityprop

    subroutine setup_thmchcoupledprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)
        integer :: i!, nb_tasks

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_sp) stop 'Size problem'
        allocate(mat%Hprop(nnz))
        !!!$OMP PARALLEL
        !nb_tasks = omp_get_num_threads()
        !!!$OMP DO SCHEDULE(STATIC, nnz/nb_tasks)
        do i = 1, nnz
            mat%Hprop(i) = prop(i)*mat%detJ(i)
        end do
        !!!$OMP END DO
        !!!$OMP END PARALLEL
    end subroutine setup_thmchcoupledprop

    subroutine compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
        use separatevariables
        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        integer, intent(in) :: nc_list
        dimension :: nc_list(mat%dimen)

        double precision, intent(out) :: univMcoefs(mat%dimen, maxval(nc_list)), &
                                        univKcoefs(mat%dimen, maxval(nc_list))

        ! Local data
        ! ----------
        type(sepoperator) :: oper
        integer :: i, gp
        integer, allocatable, dimension(:) :: nc_list_t
        logical, allocatable, dimension(:) :: update
        double precision, allocatable, dimension(:, :) :: CC

        if (.not.allocated(mat%Cprop)) then
            allocate(CC(mat%dimen, mat%ncols_sp), update(mat%dimen), nc_list_t(mat%dimen))
            update = .true.; nc_list_t = nc_list
            call initialize_operator(oper, mat%dimen, nc_list_t, update)
            
            do gp = 1, mat%ncols_sp
                do i = 1, mat%dimen
                    CC(i, gp) = mat%Kprop(i, i, gp)
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
                do i = 1, mat%dimen
                    CC(i, gp) = mat%Kprop(i, i, gp)
                end do
                CC(mat%dimen+1, gp) = mat%Cprop(gp)
            end do

            if (mat%dimen.eq.2) then
                call separatevariables_3d(oper, CC)
            else if (mat%dimen.eq.3) then
                call separatevariables_4d(oper, CC)
            end if
            univMcoefs = oper%univmasscoefs(:mat%dimen, :); univKcoefs = oper%univstiffcoefs(:mat%dimen, :)
        end if
        
    end subroutine compute_separationvariables

    subroutine compute_variablesmean(mat, nclist)
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
        
        ! Select sample
        do i = 1, mat%dimen 
            indlist(i, :) =  (/1, int((nclist(i) + 1)/2), nclist(i)/)
        end do

        call indices2list(mat%dimen, NP, indlist, nclist, sample)
        
        allocate(CC_K(mat%dimen, size(sample)))
        do c = 1, size(sample)
            gp = sample(c)
            do i = 1, mat%dimen
                CC_K(i, c) = mat%Kprop(i, i, gp)
            end do
        end do
        do i = 1, mat%dimen
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(NP, NP, CC_K(i, :), mat%Kmean(i))
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(NP, NP, NP, CC_K(i, :), mat%Kmean(i))
            end if
        end do

        if (allocated(mat%Cprop)) then
            allocate(CC_M(size(sample)))
            do c = 1, size(sample)
                gp = sample(c)
                CC_M(c) = mat%Cprop(gp)
            end do
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(NP, NP, CC_M, mat%Cmean)
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(NP, NP, NP, CC_M, mat%Cmean)
            end if
        end if   

    end subroutine compute_variablesmean

    subroutine mf_u_v(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        double precision :: tmp_in, tmp
        dimension :: tmp_in(nr_total), tmp(basisdata%nc_total)
        integer :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indjT_u, indjT_v, indjT_w
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'

        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        allocate(indi_u(nr_u+1), indj_u(nnz_u), data_W_u(nnz_u, 4), indiT_u(nc_u+1), indjT_u(nnz_u), data_BT_u(nnz_u, 2))
        indi_u = basisdata%indi(1, 1:nr_u+1); indj_u = basisdata%indj(1, 1:nnz_u)
        data_W_u = basisdata%data_bw(1, 1:nnz_u, 3:6)
        indiT_u = basisdata%indiT(1, 1:nc_u+1); indjT_u = basisdata%indjT(1, 1:nnz_u)
        data_BT_u = basisdata%data_bwT(1, 1:nnz_u, 1:2)
        allocate(indi_v(nr_v+1), indj_v(nnz_v), data_W_v(nnz_v, 4), indiT_v(nc_v+1), indjT_v(nnz_v), data_BT_v(nnz_v, 2))
        indi_v = basisdata%indi(2, 1:nr_v+1); indj_v = basisdata%indj(2, 1:nnz_v)
        data_W_v = basisdata%data_bw(2, 1:nnz_v, 3:6)
        indiT_v = basisdata%indiT(2, 1:nc_v+1); indjT_v = basisdata%indjT(2, 1:nnz_v)
        data_BT_v = basisdata%data_bwT(2, 1:nnz_v, 1:2)
        if (basisdata%dimen.eq.3) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        tmp_in = array_in; if (mat%isLumped) tmp_in = 1.d0
        if (basisdata%dimen.eq.2) then
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), &
                                tmp_in, tmp)
        else if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), &
                                nnz_w, indiT_w, indjT_w, data_BT_w(:, 1), &
                                tmp_in, tmp)
        end if

        !!!$OMP PARALLEL
        !!!$OMP WORKSHARE 
        tmp = tmp*mat%Cprop
        !!!$OMP END WORKSHARE
        !!!$OMP END PARALLEL

        if (basisdata%dimen.eq.2) then
            call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                tmp, array_out)
        else if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                tmp, array_out)
        end if
        if (mat%isLumped) array_out = array_out*array_in
        
    end subroutine mf_u_v

    subroutine mf_gradu_gradv(mat, basisdata, nr_total, array_in, array_out)
        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        double precision :: tmp_0, tmp_1, tmp_2
        dimension :: tmp_0(basisdata%nc_total), tmp_1(basisdata%nc_total), tmp_2(nr_total)
        integer :: i, j, alpha, beta, zeta
        dimension :: alpha(basisdata%dimen), beta(basisdata%dimen), zeta(basisdata%dimen)
        integer :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indjT_u, indjT_v, indjT_w
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'

        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        allocate(indi_u(nr_u+1), indj_u(nnz_u), data_W_u(nnz_u, 4), indiT_u(nc_u+1), indjT_u(nnz_u), data_BT_u(nnz_u, 2))
        indi_u = basisdata%indi(1, 1:nr_u+1); indj_u = basisdata%indj(1, 1:nnz_u)
        data_W_u = basisdata%data_bw(1, 1:nnz_u, 3:6)
        indiT_u = basisdata%indiT(1, 1:nc_u+1); indjT_u = basisdata%indjT(1, 1:nnz_u)
        data_BT_u = basisdata%data_bwT(1, 1:nnz_u, 1:2)
        allocate(indi_v(nr_v+1), indj_v(nnz_v), data_W_v(nnz_v, 4), indiT_v(nc_v+1), indjT_v(nnz_v), data_BT_v(nnz_v, 2))
        indi_v = basisdata%indi(2, 1:nr_v+1); indj_v = basisdata%indj(2, 1:nnz_v)
        data_W_v = basisdata%data_bw(2, 1:nnz_v, 3:6)
        indiT_v = basisdata%indiT(2, 1:nc_v+1); indjT_v = basisdata%indjT(2, 1:nnz_v)
        data_BT_v = basisdata%data_bwT(2, 1:nnz_v, 1:2)
        if (basisdata%dimen.eq.3) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        array_out = 0.d0
        do j = 1, basisdata%dimen
            beta = 1; beta(j) = 2
            if (basisdata%dimen.eq.2) then
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), &
                                    array_in, tmp_0)
            else if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), &
                                    nnz_w, indiT_w, indjT_w, data_BT_w(:, beta(3)), &
                                    array_in, tmp_0)
            end if
            do i = 1, basisdata%dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                !!!$OMP PARALLEL
                !!!$OMP WORKSHARE 
                tmp_1 = tmp_0*mat%Kprop(i, j, :)
                !!!$OMP END WORKSHARE
                !!!$OMP END PARALLEL
                
                if (basisdata%dimen.eq.2) then
                    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        tmp_1, tmp_2)
                else if (basisdata%dimen.eq.3) then
                    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                        tmp_1, tmp_2)
                end if
                array_out = array_out + tmp_2
            end do
        end do

    end subroutine mf_gradu_gradv

    subroutine mf_uv_gradugradv(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp1, array_tmp2
        dimension :: array_tmp1(nr_total), array_tmp2(nr_total)

        !!!$OMP PARALLEL NUM_THREADS(omp_get_num_threads())
        !!!$OMP SINGLE 
        call mf_u_v(mat, basisdata, nr_total, array_in, array_tmp1)
        !!!$OMP END SINGLE NOWAIT

        !!!$OMP SINGLE 
        call mf_gradu_gradv(mat, basisdata, nr_total, array_in, array_tmp2)
        !!!$OMP END SINGLE NOWAIT
        !!!$OMP END PARALLEL

        !!!$OMP PARALLEL
        !!!$OMP WORKSHARE
        array_out = mat%scalars(1)*array_tmp1 + mat%scalars(2)*array_tmp2
        !!!$OMP END WORKSHARE
        !!!$OMP END PARALLEL

    end subroutine mf_uv_gradugradv

    subroutine mf_gradu_tv(mat, basisdata, nr_total, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        type(thermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total

        double precision, intent(in) :: array_in
        dimension :: array_in(basisdata%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        integer :: i, k, alpha, beta, zeta
        dimension :: alpha(basisdata%dimen), beta(basisdata%dimen), zeta(basisdata%dimen)
        double precision :: t1, t2, t3
        dimension :: t1(basisdata%nc_total), t2(basisdata%nc_total), t3(nr_total)
        integer :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indjT_u, indjT_v, indjT_w
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'

        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        allocate(indi_u(nr_u+1), indj_u(nnz_u), data_W_u(nnz_u, 4), indiT_u(nc_u+1), indjT_u(nnz_u), data_BT_u(nnz_u, 2))
        indi_u = basisdata%indi(1, 1:nr_u+1); indj_u = basisdata%indj(1, 1:nnz_u)
        data_W_u = basisdata%data_bw(1, 1:nnz_u, 3:6)
        indiT_u = basisdata%indiT(1, 1:nc_u+1); indjT_u = basisdata%indjT(1, 1:nnz_u)
        data_BT_u = basisdata%data_bwT(1, 1:nnz_u, 1:2)
        allocate(indi_v(nr_v+1), indj_v(nnz_v), data_W_v(nnz_v, 4), indiT_v(nc_v+1), indjT_v(nnz_v), data_BT_v(nnz_v, 2))
        indi_v = basisdata%indi(2, 1:nr_v+1); indj_v = basisdata%indj(2, 1:nnz_v)
        data_W_v = basisdata%data_bw(2, 1:nnz_v, 3:6)
        indiT_v = basisdata%indiT(2, 1:nc_v+1); indjT_v = basisdata%indjT(2, 1:nnz_v)
        data_BT_v = basisdata%data_bwT(2, 1:nnz_v, 1:2)
        if (basisdata%dimen.eq.3) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        array_out = 0.d0
        do i = 1, basisdata%dimen
            beta = 1
            if (basisdata%dimen.eq.2) then
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), &
                                    array_in(i, :), t1)  
            else if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), &
                                    nnz_w, indiT_w, indjT_w, data_BT_w(:, beta(3)), &
                                    array_in(i, :), t1)  
            end if

            !!!$OMP PARALLEL
            !!!$OMP WORKSHARE
            t1 = t1*mat%Hprop
            !!!$OMP END WORKSHARE
            !!!$OMP END PARALLEL
            
            do k = 1, basisdata%dimen
                alpha = 1; alpha(k) = 2; zeta = beta + (alpha - 1)*2

                !!!$OMP PARALLEL
                !!!$OMP WORKSHARE
                t2 = t1*mat%invJ(k, i, :)
                !!!$OMP END WORKSHARE
                !!!$OMP END PARALLEL
                
                if (basisdata%dimen.eq.2) then
                    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        t2, t3)
                else if (basisdata%dimen.eq.3) then
                    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                        t2, t3)
                end if
                array_out = array_out + t3
            end do
        end do
                    
    end subroutine mf_gradu_tv

end module matrixfreeheat

module heatsolver
    use omp_lib
    use matrixfreeheat
    use structured_data
    type cgsolver
        logical :: withdiag = .true., applyfd = .true.
        integer :: matrixfreetype = 1
        type(basis_data), pointer :: globsyst=>null()
        type(reduced_system), pointer :: redsyst=>null()
    end type cgsolver

contains

    subroutine matrixfree_matvec(solv, mat, nr_total, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        if (solv%matrixfreetype.eq.1) then
            call mf_u_v(mat, solv%globsyst, nr_total, array_in, array_out)
        else if (solv%matrixfreetype.eq.2) then
            call mf_gradu_gradv(mat, solv%globsyst, nr_total, array_in, array_out)
        else if (solv%matrixfreetype.eq.3) then
            call mf_uv_gradugradv(mat, solv%globsyst, nr_total, array_in, array_out)
        else 
            stop 'Not coded'
        end if

    end subroutine matrixfree_matvec

    subroutine initialize_solver(solv, globsyst, redsyst)
        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(basis_data), target :: globsyst
        type(reduced_system), target :: redsyst
        solv%globsyst => globsyst
        solv%redsyst => redsyst
    end subroutine initialize_solver

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
        integer :: nr_u, nr_v, nr_w
        double precision, allocatable, dimension(:) :: tmp, tmp2

        if (.not.solv%applyfd) then
            array_out = array_in
            return
        end if
        nr_u = solv%redsyst%basisdata%nrows(1)
        nr_v = solv%redsyst%basisdata%nrows(2)
        nr_w = 1
        if (solv%globsyst%dimen.eq.3) nr_w = solv%redsyst%basisdata%nrows(3)

        array_out = 0.d0
        allocate(tmp(nr_u*nr_v*nr_w))
        if (solv%globsyst%dimen.eq.2) then
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, &
                        transpose(solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        array_in(solv%redsyst%dof), tmp)
        else if (solv%globsyst%dimen.eq.3) then
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        transpose(solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        transpose(solv%redsyst%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                        array_in(solv%redsyst%dof), tmp)
        end if

        if (solv%withdiag) then
            !!!$OMP PARALLEL 
            !!!$OMP WORKSHARE
            tmp = tmp/solv%redsyst%diageigval_sp
            !!!$OMP END WORKSHARE
            !!!$OMP END PARALLEL
        end if

        allocate(tmp2(nr_u*nr_v*nr_w))
        if (solv%globsyst%dimen.eq.2) then
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, &
                        solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), tmp, tmp2)
        else if (solv%globsyst%dimen.eq.3) then
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                        solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        solv%redsyst%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), tmp, tmp2)
        end if

        array_out(solv%redsyst%dof) = tmp2
        deallocate(tmp, tmp2)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, size_inout, array_inout)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: size_inout
        double precision, intent(inout) :: array_inout(size_inout)
        call set2zero(solv%redsyst, size_inout, array_inout)
    end subroutine clear_dirichlet

    subroutine PBiCGSTAB(solv, mat, nr_total, iterations, threshold, b, x, residual)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, iterations
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, residual
        dimension :: x(nr_total), residual(iterations+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension :: r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
                    ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
        integer :: k

        x = 0.d0; r = b; residual = 0.d0
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        rsold = dot_product(r, rhat); normb = norm2(r)
        if (normb.le.1.d-14) return
        residual(1) = 1.d0

        do k = 1, iterations
            call applyfastdiag(solv, nr_total, p, ptilde)
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_matvec(solv, mat, nr_total, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            alpha = rsold/dot_product(Aptilde, rhat)
            s = r - alpha*Aptilde
            x = x + alpha*ptilde
            if (norm2(s).le.max(threshold*normb, 1.d-14)) exit

            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_matvec(solv, mat, nr_total, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
            omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
            r = s - omega*Astilde
            x = x + omega*stilde
            
            if (norm2(r).le.max(threshold*normb, 1.d-14)) exit
            residual(k+1) = norm2(r)/normb

            rsnew = dot_product(r, rhat)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

    subroutine RQMIN(solv, mat, nr_total, ishigher, iterations, threshold, eigenvec, eigenval)
        !! Using RQMIN algorithm to compute the stability of the transient heat problem
        
        implicit none
        ! Input / output data
        ! -------------------
        integer, parameter :: sizemat = 2
        type(cgsolver) :: solv
        type(thermomat) :: mat
        integer, intent(in) :: nr_total, iterations
        logical, intent(in) :: ishigher
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(nr_total)

        ! Local data
        ! ----------
        integer :: k, j, ii
        double precision, dimension(sizemat) :: ll
        double precision, dimension(sizemat, nr_total) :: RM1, RM2, RM3
        double precision, dimension(sizemat, sizemat) :: AA1, BB1, qq
        double precision, dimension(nr_total) :: u, v, g, gtil, p, tmp, Mg, Mgtil
        double precision :: q, gnorm, gnorm0, delta

        call random_number(eigenvec)
        call clear_dirichlet(solv, nr_total, eigenvec)

        call mf_u_v(mat, solv%globsyst, nr_total, eigenvec, u)
        call clear_dirichlet(solv, nr_total, u)
        
        q = sqrt(dot_product(eigenvec, u))
        eigenvec = eigenvec/q; u = u/q
        call mf_gradu_gradv(mat, solv%globsyst, nr_total, eigenvec, v)
        call clear_dirichlet(solv, nr_total, v)
        
        eigenval = dot_product(eigenvec, v)
        g = eigenvec; gnorm0 = norm2(g); gnorm = gnorm0

        do k = 1, iterations
            if (gnorm.le.threshold*gnorm0) return
            gtil = g
            tmp = 2*(v - eigenval*u)
            call applyfastdiag(solv, nr_total, tmp, g)
            call clear_dirichlet(solv, nr_total, g)

            if (k.eq.1) then
                p = -g
            else
                call mf_u_v(mat, solv%globsyst, nr_total, g, Mg)
                call mf_u_v(mat, solv%globsyst, nr_total, gtil, Mgtil)
                p = -g + dot_product(g, Mg)/dot_product(gtil, Mgtil)*p
            end if

            RM1(1, :) = eigenvec; RM1(2, :) = p
            RM2(1, :) = v; RM3(1, :) = u;
            call mf_gradu_gradv(mat, solv%globsyst, nr_total, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(2, :) = tmp

            call mf_u_v(mat, solv%globsyst, nr_total, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(2, :) = tmp
            
            call rayleigh_submatrix(sizemat, nr_total, RM1, RM2, AA1)
            call rayleigh_submatrix(sizemat, nr_total, RM1, RM3, BB1)
            
            call compute_eigdecomp_pdr(sizemat, AA1, BB1, ll, qq)
            do j = 1, sizemat
                if ((ll(j).lt.0.d0)) ll(j) = 0.d0
            end do
            
            if (ishigher) then 
                eigenval = maxval(ll); ii = maxloc(ll, dim=1)
            else
                eigenval = minval(ll); ii = minloc(ll, dim=1)
            end if

            if (abs(qq(1, ii)).gt.1.d-8) delta = qq(2, ii)/qq(1, ii)
    
            eigenvec = eigenvec + delta*p
            call mf_u_v(mat, solv%globsyst, nr_total, eigenvec, u)
            call clear_dirichlet(solv, nr_total, u)
            
            q = sqrt(dot_product(eigenvec, u))
            eigenvec = eigenvec/q; u = u/q
            call mf_gradu_gradv(mat, solv%globsyst, nr_total, eigenvec, v)
            call clear_dirichlet(solv, nr_total, v)            
            gnorm = norm2(g)
        end do

    end subroutine RQMIN

end module heatsolver