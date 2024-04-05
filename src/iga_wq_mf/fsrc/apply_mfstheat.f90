module matrixfreestheat
    use omp_lib
    use structured_data
    implicit none
    type stthermomat
        integer :: dimen, dimen_sp, ncols_sp, ncols_tm, ncols_total
        ! Material properties
        double precision :: Cmean
        double precision, dimension(:), allocatable :: Kmean
        double precision, dimension(:), pointer ::detJ=>null(), detG=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:), allocatable ::  Cprop, Cdersprop
        double precision, dimension(:, :), allocatable :: Kdersprop
        double precision, dimension(:, :, :), allocatable :: Kprop
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

        mat%dimen = dimen_sp + 1
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

        ! Local data 
        ! ----------
        integer :: i, j, k, nb_tasks

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_total) stop 'Size problem'
        allocate(mat%Kprop(mat%dimen_sp, mat%dimen_sp, nnz))

        !$OMP PARALLEL PRIVATE(k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, mat%ncols_total/nb_tasks) 
        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                mat%Kprop(:, :, k) = matmul(mat%invJ(:, :, i), matmul(prop(:, :, k), &
                                    transpose(mat%invJ(:, :, i))))*mat%detJ(i)*mat%detG(j)
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end subroutine setup_conductivityprop

    subroutine setup_conductivityDersprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(mat%dimen_sp, nnz)

        ! Local data 
        ! ----------
        integer :: i, j, k, nb_tasks

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_total) stop 'Size problem'
        allocate(mat%Kdersprop(mat%dimen_sp, nnz))

        !$OMP PARALLEL PRIVATE(k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, mat%ncols_total/nb_tasks) 
        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                mat%Kdersprop(:, k) = matmul(mat%invJ(:, :, i), prop(:, k))*mat%detJ(i)*mat%detG(j)
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end subroutine setup_conductivityDersprop

    subroutine setup_capacityprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        ! Local data 
        ! ----------
        integer :: i, j, k, nb_tasks

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_total) stop 'Size problem'
        allocate(mat%Cprop(nnz))
        
        !$OMP PARALLEL PRIVATE(k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, mat%ncols_total/nb_tasks) 
        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                mat%Cprop(k) = prop(k)*mat%detJ(i)
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine setup_capacityprop

    subroutine setup_capacityDersprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        ! Local data 
        ! ----------
        integer :: i, j, k, nb_tasks

        if ((.not.associated(mat%detJ)).or.(.not.associated(mat%detG))) stop 'Define geometry'
        if (nnz.ne.mat%ncols_total) stop 'Size problem'
        allocate(mat%Cdersprop(nnz))

        !$OMP PARALLEL PRIVATE(k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, mat%ncols_total/nb_tasks) 
        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                mat%Cdersprop(k) = prop(k)*mat%detJ(i)*mat%detG(j)
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

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

        allocate(CC(mat%dimen_sp+1, mat%ncols_total), update(mat%dimen_sp+1)); update = .true.
        call initialize_operator(oper, size(nc_list), nc_list, update)

        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                k = i + (j-1)*mat%ncols_sp
                do l = 1, mat%dimen_sp
                    CC(l, k) = mat%Kprop(l, l, k)
                end do
                CC(mat%dimen_sp+1, k) = mat%Cprop(k)
            end do
        end do

        if (mat%dimen_sp.eq.2) then
            call separatevariables_3d(oper, CC)
        else if (mat%dimen_sp.eq.3) then
            call separatevariables_4d(oper, CC)
        end if
        univMcoefs = oper%univmasscoefs; univKcoefs = oper%univstiffcoefs
    
    end subroutine compute_separationvariables

    subroutine compute_variablesmean(mat, nclist)
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
                do i = 1, mat%dimen_sp
                    CC_K(i, cp) = mat%Kprop(i, i, gp)
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
                CC_M(cp) = mat%Cprop(gp)
            end do
        end do

        if (mat%dimen_sp.eq.2) then
            call trapezoidal_rule_3d(NP, NP, NP, CC_M, mat%Cmean)
        else if (mat%dimen_sp.eq.3) then
            call trapezoidal_rule_4d(NP, NP, NP, NP, CC_M, mat%Cmean)
        end if

    end subroutine compute_variablesmean

    subroutine mf_u_v(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        integer :: pos_tm
        double precision :: tmp
        dimension :: tmp(basisdata%nc_total)
        integer :: nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w, data_W_t
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indiT_t, indjT_u, indjT_v, indjT_w, indjT_t
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w, data_BT_t

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'
        pos_tm = mat%dimen
        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        nr_t = basisdata%nrows(pos_tm); nc_t = basisdata%ncols(pos_tm); nnz_t = basisdata%nnzs(pos_tm)
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
        allocate(indi_t(nr_t+1), indj_t(nnz_t), data_W_t(nnz_t, 4), indiT_t(nc_t+1), indjT_t(nnz_t), data_BT_t(nnz_t, 2))
        indi_t = basisdata%indi(pos_tm, 1:nr_t+1); indj_t = basisdata%indj(pos_tm, 1:nnz_t)
        data_W_t = basisdata%data_bw(pos_tm, 1:nnz_t, 3:6)
        indiT_t = basisdata%indiT(pos_tm, 1:nc_t+1); indjT_t = basisdata%indjT(pos_tm, 1:nnz_t)
        data_BT_t = basisdata%data_bwT(pos_tm, 1:nnz_t, 1:2)
        if (basisdata%dimen.eq.4) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                array_in, tmp)
        else if (basisdata%dimen.eq.4) then
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_w, indiT_w, indjT_w, data_BT_w(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                array_in, tmp)
        end if

        !$OMP PARALLEL
        !$OMP WORKSHARE
        tmp = tmp*mat%Cdersprop
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
        
        if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                                tmp, array_out)
        else if (basisdata%dimen.eq.4) then
            call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                                tmp, array_out)
        end if
        
    end subroutine mf_u_v

    subroutine mf_u_partialt_v(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! ----------
        integer :: pos_tm
        double precision :: tmp
        dimension :: tmp(basisdata%nc_total)        
        integer :: nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w, data_W_t
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indiT_t, indjT_u, indjT_v, indjT_w, indjT_t
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w, data_BT_t

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'
        pos_tm = mat%dimen
        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        nr_t = basisdata%nrows(pos_tm); nc_t = basisdata%ncols(pos_tm); nnz_t = basisdata%nnzs(pos_tm)
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
        allocate(indi_t(nr_t+1), indj_t(nnz_t), data_W_t(nnz_t, 4), indiT_t(nc_t+1), indjT_t(nnz_t), data_BT_t(nnz_t, 2))
        indi_t = basisdata%indi(pos_tm, 1:nr_t+1); indj_t = basisdata%indj(pos_tm, 1:nnz_t)
        data_W_t = basisdata%data_bw(pos_tm, 1:nnz_t, 3:6)
        indiT_t = basisdata%indiT(pos_tm, 1:nc_t+1); indjT_t = basisdata%indjT(pos_tm, 1:nnz_t)
        data_BT_t = basisdata%data_bwT(pos_tm, 1:nnz_t, 1:2)
        if (basisdata%dimen.eq.4) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 2), & 
                                array_in, tmp)
        else if (basisdata%dimen.eq.4) then
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_w, indiT_w, indjT_w, data_BT_w(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 2), & 
                                array_in, tmp)
        end if

        !$OMP PARALLEL
        !$OMP WORKSHARE
        tmp = tmp*mat%Cprop
        !$OMP END WORKSHARE
        !$OMP END PARALLEL

        if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_t, indi_t, indj_t, data_W_t(:, 2), &
                                tmp, array_out)
        else if (basisdata%dimen.eq.4) then
            call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                nnz_t, indi_t, indj_t, data_W_t(:, 2), &
                                tmp, array_out)
        end if
        
    end subroutine mf_u_partialt_v

    subroutine mf_gradx_u_v(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        integer :: pos_tm, i, alpha, zeta 
        dimension :: alpha(basisdata%dimen), zeta(basisdata%dimen)
        double precision :: tmp_0, tmp_1, tmp_2
        dimension :: tmp_0(basisdata%nc_total), tmp_1(basisdata%nc_total), tmp_2(basisdata%nr_total)

        integer :: nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w, data_W_t
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indiT_t, indjT_u, indjT_v, indjT_w, indjT_t
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w, data_BT_t

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'
        pos_tm = mat%dimen
        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        nr_t = basisdata%nrows(pos_tm); nc_t = basisdata%ncols(pos_tm); nnz_t = basisdata%nnzs(pos_tm)
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
        allocate(indi_t(nr_t+1), indj_t(nnz_t), data_W_t(nnz_t, 4), indiT_t(nc_t+1), indjT_t(nnz_t), data_BT_t(nnz_t, 2))
        indi_t = basisdata%indi(pos_tm, 1:nr_t+1); indj_t = basisdata%indj(pos_tm, 1:nnz_t)
        data_W_t = basisdata%data_bw(pos_tm, 1:nnz_t, 3:6)
        indiT_t = basisdata%indiT(pos_tm, 1:nc_t+1); indjT_t = basisdata%indjT(pos_tm, 1:nnz_t)
        data_BT_t = basisdata%data_bwT(pos_tm, 1:nnz_t, 1:2)
        if (basisdata%dimen.eq.4) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        array_out = 0.d0
        if (basisdata%dimen.eq.3) then
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                array_in, tmp_0)
        else if (basisdata%dimen.eq.4) then
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                                nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                nnz_w, indiT_w, indjT_w, data_BT_w(:, 1), & 
                                nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                array_in, tmp_0)
        end if
        
        do i = 1, mat%dimen_sp
            alpha = 1; alpha(i) = 2
            zeta  = 1 + (alpha - 1)*2

            !$OMP PARALLEL
            !$OMP WORKSHARE
            tmp_1 = tmp_0*mat%Kdersprop(i, :)
            !$OMP END WORKSHARE
            !$OMP END PARALLEL
            
            if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                                    tmp_1, tmp_2)
            else if (basisdata%dimen.eq.4) then
                call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                    nnz_t, indi_t, indj_t, data_W_t(:, 1), &
                                    tmp_1, tmp_2)
            end if
            array_out = array_out + tmp_2
        end do

    end subroutine mf_gradx_u_v

    subroutine mf_gradx_u_gradx_v(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data
        ! -------------------
        type(stthermomat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data 
        ! -----------
        integer :: pos_tm, i, j, alpha, beta, zeta
        dimension :: alpha(basisdata%dimen), beta(basisdata%dimen), zeta(basisdata%dimen)
        double precision :: tmp_0, tmp_1, tmp_2
        dimension :: tmp_0(basisdata%nc_total), tmp_1(basisdata%nc_total), tmp_2(basisdata%nr_total)
        integer :: nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w, data_W_t
        integer, dimension(:), allocatable :: indiT_u, indiT_v, indiT_w, indiT_t, indjT_u, indjT_v, indjT_w, indjT_t
        double precision, dimension(:, :), allocatable :: data_BT_u, data_BT_v, data_BT_w, data_BT_t

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'
        pos_tm = mat%dimen
        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        nr_t = basisdata%nrows(pos_tm); nc_t = basisdata%ncols(pos_tm); nnz_t = basisdata%nnzs(pos_tm)
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
        allocate(indi_t(nr_t+1), indj_t(nnz_t), data_W_t(nnz_t, 4), indiT_t(nc_t+1), indjT_t(nnz_t), data_BT_t(nnz_t, 2))
        indi_t = basisdata%indi(pos_tm, 1:nr_t+1); indj_t = basisdata%indj(pos_tm, 1:nnz_t)
        data_W_t = basisdata%data_bw(pos_tm, 1:nnz_t, 3:6)
        indiT_t = basisdata%indiT(pos_tm, 1:nc_t+1); indjT_t = basisdata%indjT(pos_tm, 1:nnz_t)
        data_BT_t = basisdata%data_bwT(pos_tm, 1:nnz_t, 1:2)
        if (basisdata%dimen.eq.4) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4), indiT_w(nc_w+1), indjT_w(nnz_w), data_BT_w(nnz_w, 2))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
            indiT_w = basisdata%indiT(3, 1:nc_w+1); indjT_w = basisdata%indjT(3, 1:nnz_w)
            data_BT_w = basisdata%data_bwT(3, 1:nnz_w, 1:2)
        end if

        array_out = 0.d0
        do j = 1, mat%dimen_sp
            beta = 1; beta(j) = 2
            if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                    nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                    array_in, tmp_0)
            else if (basisdata%dimen.eq.4) then
                call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                    nnz_w, indiT_w, indjT_w, data_BT_w(:, beta(3)), & 
                                    nnz_t, indiT_t, indjT_t, data_BT_t(:, 1), & 
                                    array_in, tmp_0)
            end if

            do i = 1, mat%dimen_sp
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2

                !$OMP PARALLEL
                !$OMP WORKSHARE
                tmp_1 = tmp_0*mat%Kprop(i, j, :)
                !$OMP END WORKSHARE
                !$OMP END PARALLEL
                
                if (basisdata%dimen.eq.3) then
                    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_t, indi_t, indj_t, data_W_t(:, 1), &  
                                        tmp_1, tmp_2)
                else if (basisdata%dimen.eq.4) then
                    call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, & 
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                        nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                        tmp_1, tmp_2)
                end if
                array_out = array_out + tmp_2
            end do
        end do
    end subroutine mf_gradx_u_gradx_v

end module matrixfreestheat

module stheatsolver
    use omp_lib
    use matrixfreestheat
    use structured_data
    type stcgsolver
        integer :: matrixfreetype = 1
        logical :: applyfd = .true.
        type(basis_data), pointer :: globsyst=>null()
        type(reduced_system), pointer :: redsyst=>null()
        double precision :: scalarleft = 1.d0
    end type stcgsolver

contains

    subroutine matrixfree_matvec(solv, mat, nr_total, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)
        
        ! Local data
        ! ----------
        double precision :: tmp1(nr_total), tmp2(nr_total)

        !$OMP PARALLEL NUM_THREADS(omp_get_num_threads())
        !$OMP SINGLE 
        call mf_u_partialt_v(mat, solv%globsyst, nr_total, array_in, tmp1)
        !$OMP END SINGLE NOWAIT

        !$OMP SINGLE 
        call mf_gradx_u_gradx_v(mat, solv%globsyst, nr_total, array_in, tmp2)
        !$OMP END SINGLE NOWAIT
        !$OMP END PARALLEL

        !$OMP PARALLEL
        !$OMP WORKSHARE
        array_out = tmp1 + tmp2
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
        
        if (solv%matrixfreetype.eq.1) return 
        if (solv%matrixfreetype.eq.2) then
            !$OMP PARALLEL NUM_THREADS(omp_get_num_threads())
            !$OMP SINGLE 
            call mf_u_v(mat, solv%globsyst, nr_total, array_in, tmp1)
            !$OMP END SINGLE NOWAIT

            !$OMP SINGLE
            call mf_gradx_u_v(mat, solv%globsyst, nr_total, array_in, tmp2)
            !$OMP END SINGLE NOWAIT
            !$OMP END PARALLEL

            !$OMP PARALLEL
            !$OMP WORKSHARE
            array_out = array_out + tmp1 + tmp2
            !$OMP END WORKSHARE
            !$OMP END PARALLEL
            
        else
            stop 'Not coded'
        end if

    end subroutine matrixfree_matvec

    subroutine initialize_solver(solv, globsyst, redsyst)
        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(basis_data), target :: globsyst
        type(reduced_system), target :: redsyst
        solv%globsyst => globsyst
        solv%redsyst => redsyst
    end subroutine initialize_solver

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
        Vb = matmul(transpose(conjg(solv%redsyst%Lschur_tm)), dcmplx(b))
        uptrg = dcmplx(coefs(1))*solv%redsyst%Luptr_tm + dcmplx(coefs(2))*solv%redsyst%Ruptr_tm
        call solve_complex_uppertriangular_system(nr, uptrg, Vb, sol)
        Usol = matmul(solv%redsyst%Rschur_tm, sol)
        x = realpart(Usol)
    end subroutine solve_schurtriangular__

    subroutine applyfastdiag(solv, nr_total, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! Applied to steady heat problems
        !! by G. Sanaglli and M. Tani
        
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
        integer :: pos_tm, i, j, k, l, nb_tasks
        integer :: nr_u, nr_v, nr_w, nr_t
        double precision, allocatable, dimension(:, :) :: identity
        double precision, allocatable, dimension(:) :: tmp1, tmp3, btmp, stmp
        double precision, allocatable, dimension(:, :, :, :) :: tmp2
        integer, allocatable, dimension(:, :, :, :) :: dof
        double precision :: eigval
        
        if (.not.solv%applyfd) then
            array_out = array_in
            return
        end if
        pos_tm = solv%globsyst%dimen
        nr_u = solv%redsyst%basisdata%nrows(1)
        nr_v = solv%redsyst%basisdata%nrows(2)
        nr_t = solv%redsyst%basisdata%nrows(pos_tm)
        nr_w = 1; if (solv%globsyst%dimen.eq.4) nr_w = solv%redsyst%basisdata%nrows(3)

        allocate(identity(nr_t, nr_t))
        call create_identity(nr_t, identity)
        allocate(dof(nr_u, nr_v, nr_w, nr_t))
        dof = reshape(solv%redsyst%dof, shape=(/nr_u, nr_v, nr_w, nr_t/))

        array_out = 0.d0
        allocate(tmp1(nr_u*nr_v*nr_w*nr_t))
        if (solv%globsyst%dimen.eq.3) then
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_t, nr_t, &
                        transpose(solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        identity, array_in(solv%redsyst%dof), tmp1)
        else if (solv%globsyst%dimen.eq.4) then
            call sumfacto4d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        transpose(solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                        transpose(solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                        transpose(solv%redsyst%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                        identity, array_in(solv%redsyst%dof), tmp1)
        end if 

        allocate(tmp2(nr_u, nr_v, nr_w, nr_t))
        tmp2 = reshape(tmp1, shape=(/nr_u, nr_v, nr_w, nr_t/))
        deallocate(tmp1)
        allocate(btmp(nr_t), stmp(nr_t))

        !$OMP PARALLEL PRIVATE(l, eigval, btmp, stmp)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nr_w*nr_v*nr_u/nb_tasks) 
        do k = 1, nr_w
            do j = 1, nr_v
                do i = 1, nr_u
                    l = i + (j - 1)*nr_u + (k - 1)*nr_u*nr_v
                    eigval = solv%redsyst%diageigval_sp(l)
                    btmp = tmp2(i, j, k, :)
                    call solve_schurtriangular__(solv, nr_t, (/solv%scalarleft, eigval/), btmp, stmp)
                    tmp2(i, j, k, :) = stmp
                end do
            end do
        end do  
        !$OMP END DO
        !$OMP END PARALLEL
        deallocate(stmp, btmp)
        allocate(tmp1(nr_u*nr_v*nr_w*nr_t))
        tmp1 = pack(tmp2, .true.)
        deallocate(tmp2)

        allocate(tmp3(nr_u*nr_v*nr_w*nr_t))
        if (solv%globsyst%dimen.eq.3) then
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_t, nr_t, &
                        solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        identity, tmp1, tmp3)
        else if (solv%globsyst%dimen.eq.4) then
            call sumfacto4d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, nr_t, nr_t, &
                        solv%redsyst%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                        solv%redsyst%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                        solv%redsyst%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), &
                        identity, tmp1, tmp3)
        end if 
        array_out(solv%redsyst%dof) = tmp3
        deallocate(tmp1, tmp3)

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, size_inout, array_inout)
        implicit none
        ! Input / output  data 
        !---------------------
        type(stcgsolver) :: solv
        integer, intent(in) :: size_inout
        double precision, intent(inout) :: array_inout(size_inout)
        call set2zero(solv%redsyst, size_inout, array_inout)
    end subroutine clear_dirichlet

    subroutine PBiCGSTAB(solv, mat, nr_total, iterations, threshold, b, x, residual)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total
        integer, intent(in) :: iterations
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, residual
        dimension :: x(nr_total), residual(iterations+1)

        ! Local data
        ! -----------
        double precision :: rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
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

    subroutine PGMRES(solv, mat, nr_total, nbRestarts, iterations, threshold, b, x, residual)

        implicit none
        ! Input / output data
        ! -------------------
        type(stcgsolver) :: solv
        type(stthermomat) :: mat
        integer, intent(in) :: nr_total
        integer, intent(in) :: nbRestarts, iterations
        double precision, intent(in) :: threshold, b
        dimension :: b(nr_total)
        
        double precision, intent(out) :: x, residual
        dimension :: x(nr_total), residual(nbRestarts*(iterations+1))

        ! Local data
        ! -----------
        double precision :: H, V, Z, beta, e1, y
        dimension :: H(iterations+1, iterations), V(iterations+1, nr_total), &
                    Z(iterations+1, nr_total), beta(nbRestarts+1), e1(iterations+1), y(iterations)
        double precision :: r(nr_total), w(nr_total), Ax(nr_total), Pv(nr_total)
        integer :: i, j, k, c

        e1 = 0.d0; beta = 0.d0; x = 0.d0; residual = 0.d0; c = 1
        do k = 1, nbRestarts
            H = 0.d0; V = 0.d0; Z = 0.d0; y = 0.d0
            call matrixfree_matvec(solv, mat, nr_total, x, Ax)
            r = b - Ax
            call clear_dirichlet(solv, nr_total, r)
            beta(k) = norm2(r)
            if ((k.eq.1).and.(beta(1).le.1.d-14)) return
            residual(c) = beta(k)/beta(1); c = c + 1
            if (beta(k).le.max(threshold*beta(1), 1.d-14)) exit
            V(1, :) = r/beta(1)
            e1(1) = beta(k)

            do j = 1, iterations
                call applyfastdiag(solv, nr_total, V(j, :), Pv)
                call clear_dirichlet(solv, nr_total, Pv)
                Z(j, :) = Pv
                call matrixfree_matvec(solv, mat, nr_total, Pv, w)
                call clear_dirichlet(solv, nr_total, w)
                do i = 1, j
                    H(i, j) = dot_product(w, V(i, :))
                    w = w - H(i, j)*V(i, :)
                end do
                H(j+1, j) = norm2(w)
                if (abs(H(j+1, j)).gt.1e-12) then
                    V(j+1, :) = w/H(j+1, j)
                end if
                call solve_linear_system(j+1, j, H(:j+1, :j), e1(:j+1), y(:j))
                beta(k+1) = norm2(matmul(H(:j+1, :j), y(:j)) - e1(:j+1))
                if (beta(k+1).le.max(threshold*beta(1), 1.d-14)) exit
                residual(c) = beta(k+1)/beta(1); c = c + 1
            end do
            x = x + matmul(y(:j), Z(:j, :))
        end do

    end subroutine PGMRES

end module stheatsolver

! subroutine mf_u_v(mat, basisdata, nr_total, array_in, array_out)

!     implicit none 
!     ! Input / output data
!     ! -------------------
!     type(stthermomat) :: mat
!     type(basis_data) :: basisdata
!     integer, intent(in) :: nr_total
!     double precision, intent(in) :: array_in
!     dimension :: array_in(nr_total)

!     double precision, intent(out) :: array_out
!     dimension :: array_out(nr_total)

!     ! Local data 
!     ! ----------
!     integer :: pos_tm, i, j, k
!     double precision :: tmp
!     dimension :: tmp(basisdata%nc_total)
!     integer :: nr_u, nr_v, nr_w, nr_t, nc_u, nc_v, nc_w, nc_t
!     double precision, dimension(:, :, :), allocatable :: BT_u, BT_v, BT_w, BT_t, W_u, W_v, W_w, W_t

!     if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
!     if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'

!     pos_tm = mat%dimen
!     nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1)
!     nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2)
!     nr_t = basisdata%nrows(pos_tm); nc_t = basisdata%ncols(pos_tm)
!     allocate(BT_u(nc_u, nr_u, 2), W_u(nr_u, nc_u, 4), BT_v(nc_v, nr_v, 2), W_v(nr_v, nc_v, 4), &
!             BT_t(nc_t, nr_t, 2), W_t(nr_t, nc_t, 4))
!     BT_u = basisdata%BTdense(1, 1:nc_u, 1:nr_u, :); W_u = basisdata%Wdense(1, 1:nr_u, 1:nc_u, :)
!     BT_v = basisdata%BTdense(2, 1:nc_v, 1:nr_v, :); W_v = basisdata%Wdense(2, 1:nr_v, 1:nc_v, :)
!     BT_t = basisdata%BTdense(pos_tm, 1:nc_t, 1:nr_t, :); W_t = basisdata%Wdense(pos_tm, 1:nr_t, 1:nc_t, :)
!     if (basisdata%dimen.eq.4) then
!         nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3)
!         allocate(BT_w(nc_w, nr_w, 2), W_w(nr_w, nc_w, 4))
!         BT_w = basisdata%BTdense(3, 1:nc_w, 1:nr_w, :); W_w = basisdata%Wdense(3, 1:nr_w, 1:nc_w, :)
!     end if

!     if (basisdata%dimen.eq.3) then
!         call sumfacto3d_dM(nc_u, nr_u, nc_v, nr_v, nc_t, nr_t, &
!                             BT_u(:, :, 1), BT_v(:, :, 1), BT_t(:, :, 1), & 
!                             array_in, tmp)
!     else if (basisdata%dimen.eq.4) then
!         call sumfacto4d_dM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, nc_t, nr_t, &
!                             BT_u(:, :, 1), BT_v(:, :, 1), BT_w(:, :, 1), BT_t(:, :, 1), & 
!                             array_in, tmp)
!     end if


!     if (basisdata%dimen.eq.3) then
!         call sumfacto3d_dM(nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
!                         W_u(:, :, 1), W_v(:, :, 1), W_t(:, :, 1), &
!                         tmp, array_out)
!     else if (basisdata%dimen.eq.4) then
!         call sumfacto4d_dM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
!                         W_u(:, :, 1), W_v(:, :, 1), W_w(:, :, 1), W_t(:, :, 1), &
!                         tmp, array_out)
!     end if
    
! end subroutine mf_u_v
