module matrixfreeplasticity
    use structured_data
    implicit none
    type :: mecamat
    
        integer :: dimen, nvoigt, ncols_sp
        logical :: isLumped = .false., isElastic = .true., withNN = .false., withBB = .false.
        double precision :: scalars(2) = (/1.d0, 1.d0/)
        ! Material properties
        double precision, dimension(:, :), allocatable :: Smean
        double precision, dimension(:), allocatable :: Mmean
        double precision, dimension(:), pointer :: detJ=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:), allocatable :: Mprop, Hprop
        double precision, dimension(:, :), allocatable :: CepArgs
        double precision, dimension(:, :), pointer :: NN=>null(), BB=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn, JJbb
    
    end type mecamat

contains

    subroutine setup_geometry(mat, dimen, nnz, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: dimen, nnz
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(dimen, dimen, nnz), detJ(nnz)

        mat%dimen = dimen
        mat%nvoigt = dimen*(dimen+1)/2
        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nnz

        allocate(mat%Smean(mat%dimen, mat%dimen), mat%Mmean(mat%dimen))
        mat%Smean = 1.d0; mat%Mmean = 1.d0

    end subroutine setup_geometry

    subroutine setup_massprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_sp) stop 'Size problem'
        allocate(mat%Mprop(nnz))
        mat%Mprop = prop*mat%detJ
    end subroutine setup_massprop

    subroutine setup_thmchcoupledprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        if (.not.associated(mat%detJ)) stop 'Define geometry'
        if (nnz.ne.mat%ncols_sp) stop 'Size problem'
        allocate(mat%Hprop(nnz))
        mat%Hprop = prop*mat%detJ
    end subroutine setup_thmchcoupledprop

    subroutine setup_mechanicalArguments(mat, nbrows, mechArgs)
        implicit none
        ! Input / output data
        ! -------------------
        double precision, parameter :: threshold = 1.d-12
        type(mecamat) :: mat
        integer, intent(in) :: nbrows
        double precision, target, intent(in) :: mechArgs
        dimension :: mechArgs(nbrows, mat%ncols_sp)

        ! Local data
        ! ----------
        double precision :: tensor(mat%dimen, mat%dimen)
        integer :: i

        if (nbrows.gt.2) mat%isElastic = .false.
        if (mat%isElastic) then
            if (nbrows.lt.2) stop 'Size problem'
            allocate(mat%CepArgs(2, mat%ncols_sp))
            mat%CepArgs = mechArgs(:2, :)
        else
            if (nbrows.lt.4+2*mat%nvoigt) stop 'Size problem'
            allocate(mat%CepArgs(4, mat%ncols_sp))
            mat%CepArgs = mechArgs(:4, :)  
            mat%NN      => mechArgs(5:4+mat%nvoigt, :)
            mat%BB      => mechArgs(5+mat%nvoigt:4+2*mat%nvoigt, :)

            if (.not.allocated(mat%JJnn)) allocate(mat%JJnn(mat%dimen, mat%dimen, mat%ncols_sp))
            if (.not.allocated(mat%JJbb)) allocate(mat%JJbb(mat%dimen, mat%dimen, mat%ncols_sp))
            mat%JJnn = 0.d0; mat%JJbb = 0.d0
            if (any(abs(mat%NN).gt.threshold)) then
                mat%withNN = .true.
                do i = 1, mat%ncols_sp
                    call array2symtensor(mat%dimen, mat%nvoigt, mat%NN(:, i), tensor)
                    mat%JJnn(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
                end do
            end if
            if (any(abs(mat%BB).gt.threshold)) then
                mat%withBB = .true.
                do i = 1, mat%ncols_sp
                    call array2symtensor(mat%dimen, mat%nvoigt, mat%BB(:, i), tensor)
                    mat%JJbb(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
                end do
            end if
        end if

        do i = 1, size(mat%CepArgs, dim=1)
            mat%CepArgs(i, :) = mat%CepArgs(i, :)*mat%detJ
        end do

    end subroutine setup_mechanicalArguments

    subroutine setup_jacobienjacobien(mat)
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat

        ! Local data
        ! ----------
        integer :: i
        
        if (.not.allocated(mat%JJjj)) allocate(mat%JJjj(mat%dimen, mat%dimen, mat%ncols_sp))
        mat%JJjj = 0.d0
        do i = 1, mat%ncols_sp
            mat%JJjj(:, :, i) = matmul(mat%invJ(:, :, i), transpose(mat%invJ(:, :, i)))
        end do
    end subroutine setup_jacobienjacobien

    subroutine compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
        
        use separatevariables
        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nc_list
        dimension :: nc_list(mat%dimen)

        double precision, intent(out) :: univMcoefs(mat%dimen, mat%dimen, maxval(nc_list)), &
                                        univKcoefs(mat%dimen, mat%dimen, maxval(nc_list))

        ! Local data
        ! ----------
        type(sepoperator) :: oper
        integer :: i, j, k, gp
        integer, allocatable, dimension(:) :: nc_list_t
        logical, allocatable, dimension(:) :: update
        double precision, allocatable, dimension(:, :) :: CC
        double precision :: DD, NN, TNN, BB, TBB, tensor
        dimension :: DD(mat%dimen, mat%dimen), NN(mat%nvoigt), TNN(mat%dimen, mat%dimen), &
                    BB(mat%nvoigt), TBB(mat%dimen, mat%dimen), tensor(mat%dimen, mat%dimen)

        if (.not.allocated(mat%Mprop)) then
            allocate(CC(mat%dimen, mat%ncols_sp), update(mat%dimen), nc_list_t(mat%dimen))
            update = .true.; nc_list_t = nc_list
            call initialize_operator(oper, mat%dimen, nc_list_t, update)
            do i = 1, mat%dimen
                do gp = 1, mat%ncols_sp
                    
                    ! Elastic
                    DD = 0.d0
                    DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
                    do k = 1, mat%dimen
                        DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
                    end do

                    ! Plastic
                    if (.not.mat%isElastic) then
                        NN = mat%NN(:, gp); BB = mat%BB(:, gp)
                        call array2symtensor(mat%dimen, size(NN), NN, TNN)
                        call array2symtensor(mat%dimen, size(BB), BB, TBB)
                        do j = 1, mat%dimen
                            do k = 1, mat%dimen
                                DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
                                            + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
                            end do
                        end do
                    end if
                    tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))
                    do j = 1, mat%dimen
                        CC(j, gp) = tensor(j, j)
                    end do
                end do

                if (mat%dimen.eq.2) then
                    call separatevariables_2d(oper, CC)
                else if (mat%dimen.eq.3) then
                    call separatevariables_3d(oper, CC)
                end if
                univMcoefs(i, :, :) = oper%univmasscoefs; univKcoefs(i, :, :) = oper%univstiffcoefs
            end do
        else
            allocate(CC(mat%dimen+1, mat%ncols_sp), update(mat%dimen+1), nc_list_t(mat%dimen+1))
            update = .true.; update(mat%dimen+1) = .false.
            nc_list_t(:mat%dimen) = nc_list; nc_list_t(mat%dimen+1) = 1
            call initialize_operator(oper, mat%dimen+1, nc_list_t, update)

            do i = 1, mat%dimen
                do gp = 1, mat%ncols_sp
                    
                    ! Elastic
                    DD = 0.d0
                    DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
                    do k = 1, mat%dimen
                        DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
                    end do
                    
                    ! Plastic
                    if (.not.mat%isElastic) then
                        NN = mat%NN(:, gp); BB = mat%BB(:, gp)
                        call array2symtensor(mat%dimen, size(NN), NN, TNN)
                        call array2symtensor(mat%dimen, size(BB), BB, TBB)
                        do j = 1, mat%dimen
                            do k = 1, mat%dimen
                                DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
                                            + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
                            end do
                        end do
                    end if
                    tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))
                    do j = 1, mat%dimen
                        CC(j, gp) = tensor(j, j)
                    end do
                    CC(mat%dimen+1, gp) = mat%Mprop(gp)
                end do

                if (mat%dimen.eq.2) then
                    call separatevariables_3d(oper, CC)
                else if (mat%dimen.eq.3) then
                    call separatevariables_4d(oper, CC)
                end if
                univMcoefs(i, :, :) = oper%univmasscoefs(:mat%dimen, :); univKcoefs(i, :, :) = oper%univstiffcoefs(:mat%dimen, :)
            end do
            
        end if
    end subroutine compute_separationvariables

    subroutine compute_variablesmean(mat, nclist)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nclist
        dimension :: nclist(mat%dimen)

        ! Local data
        ! ----------
        integer, parameter :: NP = 3
        integer :: i, j, k, c, gp, sample(NP**mat%dimen), indlist(mat%dimen, NP)
        double precision :: DD, TNN, NN, TBB, BB, tensor
        dimension :: DD(mat%dimen, mat%dimen), TNN(mat%dimen, mat%dimen), NN(mat%nvoigt), &
                    TBB(mat%dimen, mat%dimen), BB(mat%nvoigt), tensor(mat%dimen, mat%dimen)
        double precision, allocatable, dimension(:, :) :: Scoefs
        double precision, allocatable, dimension(:) :: Mcoefs

        ! Select sample
        do i = 1, mat%dimen 
            indlist(i, :) =  (/1, int((nclist(i) + 1)/2), nclist(i)/)
        end do

        call indices2list(mat%dimen, NP, indlist, nclist, sample)
        
        allocate(Scoefs(mat%dimen, size(sample)))
        do i = 1, mat%dimen
            do c = 1, size(sample)
                gp = sample(c)

                ! Elastic
                DD = 0.d0
                DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
                do k = 1, mat%dimen
                    DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
                end do
                
                ! Plastic
                if (.not.mat%isElastic) then
                    NN = mat%NN(:, gp); BB = mat%NN(:, gp)
                    call array2symtensor(mat%dimen, size(NN), NN, TNN)
                    call array2symtensor(mat%dimen, size(BB), BB, TBB)
                    do j = 1, mat%dimen
                        do k = 1, mat%dimen
                            DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
                                        + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
                        end do
                    end do
                end if
                
                tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))
                do j = 1, mat%dimen
                    Scoefs(j, c) = tensor(j, j)
                end do 
            end do
    
            do j = 1, mat%dimen
                if (mat%dimen.eq.2) then
                    call trapezoidal_rule_2d(3, 3, Scoefs(j, :), mat%Smean(i, j))
                else if (mat%dimen.eq.3) then
                    call trapezoidal_rule_3d(3, 3, 3, Scoefs(j, :), mat%Smean(i, j))
                end if
            end do   
        end do

        if (allocated(mat%Mprop)) then
            allocate(Mcoefs(size(sample)))
            do c = 1, size(sample)
                gp = sample(c)
                Mcoefs(c) = mat%Mprop(gp)
            end do
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(3, 3, Mcoefs, tensor(1, 1))
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(3, 3, 3, Mcoefs, tensor(1, 1))
            end if 
            mat%Mmean = tensor(1, 1)
        end if

    end subroutine compute_variablesmean

    subroutine mf_tu_tv(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        type(mecamat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(basisdata%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(basisdata%dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i
        double precision :: tmp_in, array_tmp, array_tmp2
        dimension :: tmp_in(nr_total), array_tmp(basisdata%nc_total), array_tmp2(nr_total)
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
        do i = 1, mat%dimen
            tmp_in = array_in(i, :); if (mat%isLumped) tmp_in = 1.d0
            if (basisdata%dimen.eq.2) then
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), &
                                    tmp_in, array_tmp)
            else if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, 1), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, 1), & 
                                    nnz_w, indiT_w, indjT_w, data_BT_w(:, 1), & 
                                    tmp_in, array_tmp)
            end if
            
            array_tmp = array_tmp*mat%Mprop

            if (basisdata%dimen.eq.2) then
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                    array_tmp, array_tmp2)
            else if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                    array_tmp, array_tmp2)
            end if
            array_out(i, :) = array_tmp2; if (mat%isLumped) array_out(i, :) = array_tmp2*array_in(i, :)
        end do

    end subroutine mf_tu_tv

    subroutine mf_gradtu_gradtv(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        type(mecamat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total

        double precision, intent(in) :: array_in
        dimension :: array_in(basisdata%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(basisdata%dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, j, k, l, m, alpha, beta, zeta, nbCepArgs = 2
        dimension :: alpha(basisdata%dimen), beta(basisdata%dimen), zeta(basisdata%dimen)
        double precision, allocatable, dimension(:) :: t4, t5, t6
        double precision, allocatable, dimension(:, :) :: kt1
        double precision :: t1, t2, t3, t7, t8, t9
        dimension :: t1(basisdata%nc_total), t2(basisdata%nc_total), t3(basisdata%nc_total), &
                        t7(basisdata%nc_total), t8(basisdata%nc_total), t9(basisdata%nr_total)
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
                
        if (.not.mat%isElastic) then 
            nbCepArgs = 4
            allocate(t4(basisdata%nc_total), t5(basisdata%nc_total), t6(basisdata%nc_total))
        end if
        allocate(kt1(nbCepArgs, basisdata%nc_total))

        array_out = 0.d0
        do j = 1, mat%dimen
            do m = 1, mat%dimen
                beta = 1; beta(m) = 2
                if (basisdata%dimen.eq.2) then
                    call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                        nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                        nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                        array_in(j, :), t1) 
                else if (basisdata%dimen.eq.3) then
                    call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                        nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                        nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                        nnz_w, indiT_w, indjT_w, data_BT_w(:, beta(3)), &
                                        array_in(j, :), t1) 
                end if

                do k = 1, nbCepArgs
                    kt1(k, :) = mat%CepArgs(k, :)*t1
                end do

                t2 = kt1(1, :)*mat%invJ(m, j, :)
                if (.not.mat%isElastic) then
                    t4 = 0.d0; t5 = 0.d0; t6 = 0.d0
                    if (mat%withNN) t4 = kt1(3, :)*mat%JJnn(m, j, :)
                    if (mat%withBB) t5 = kt1(4, :)*mat%JJbb(m, j, :)
                    if (mat%withNN) t6 = kt1(4, :)*mat%JJnn(m, j, :)
                end if

                do i = 1, mat%dimen
                    t3 = kt1(2, :)*mat%invJ(m, i, :)
                    t9 = 0.d0

                    do l = 1, mat%dimen
                        alpha = 1; alpha(l) = 2
                        zeta  = beta + (alpha - 1)*2
                        
                        t7 = t2*mat%invJ(l, i, :) + t3*mat%invJ(l, j, :)
                        if (.not.mat%isElastic) then
                            if (mat%withNN) t7 = t7 + t4*mat%JJnn(l, i, :) - t5*mat%JJnn(l, i, :)  
                            if (mat%withBB) t7 = t7 + t6*mat%JJbb(l, i, :)
                        end if

                        if (i.eq.j) t7 = t7 + kt1(2, :)*mat%JJjj(l, m, :)
                        if (basisdata%dimen.eq.2) then
                            call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                                t7, t8) 
                        else if (basisdata%dimen.eq.3) then
                            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                                t7, t8)
                        end if
                        
                        t9 = t9 + t8
                    end do

                    array_out(i, :) = array_out(i, :) + t9
                end do
            end do
        end do
            
    end subroutine mf_gradtu_gradtv

    subroutine mf_tutv_gradtugradtv(mat, basisdata, nr_total, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        type(mecamat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total

        double precision, intent(in) :: array_in
        dimension :: array_in(basisdata%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(basisdata%dimen, nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(basisdata%dimen, nr_total)

        call mf_tu_tv(mat, basisdata, nr_total, array_in, array_out)
        array_out = mat%scalars(1)*array_out
        call mf_gradtu_gradtv(mat, basisdata, nr_total, array_in, array_tmp)
        array_out = array_out + mat%scalars(2)*array_tmp

    end subroutine mf_tutv_gradtugradtv

    subroutine mf_u_gradtv(mat, basisdata, nr_total, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        type(mecamat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(basisdata%dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, l, alpha, beta, zeta
        dimension :: alpha(basisdata%dimen), beta(basisdata%dimen), zeta(basisdata%dimen)
        double precision :: t1, t2, t3
        dimension :: t1(basisdata%nc_total), t2(basisdata%nc_total), t3(basisdata%nr_total)
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
        do l = 1, basisdata%dimen
            beta = 1; beta(l) = 2
            if (basisdata%dimen.eq.2) then
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                    array_in, t1) 
            else if (basisdata%dimen.eq.3) then
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                    nnz_u, indiT_u, indjT_u, data_BT_u(:, beta(1)), &
                                    nnz_v, indiT_v, indjT_v, data_BT_v(:, beta(2)), & 
                                    nnz_w, indiT_w, indjT_w, data_BT_w(:, beta(3)), &
                                    array_in, t1) 
            end if
            
            t1 = t1*mat%Hprop
            do i = 1, basisdata%dimen
                alpha = 1; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(l, i, :)

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
                
                array_out(i, :) = array_out(i, :) + t3
            end do
        end do

    end subroutine mf_u_gradtv

    subroutine intforce(mat, basisdata, nr_total, nc_total, stress, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_total, nc_total
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(basisdata%dimen, nr_total)

        ! Local data
        ! ----------
        double precision :: Tstress, t1, t2
        dimension :: Tstress(basisdata%dimen, basisdata%dimen), &
                    t1(basisdata%dimen, basisdata%dimen, nc_total), t2(nr_total)
        integer :: i, k, alpha(basisdata%dimen), zeta(basisdata%dimen)
        integer :: nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, dimension(:), allocatable :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        double precision, dimension(:, :), allocatable :: data_W_u, data_W_v, data_W_w

        if (nr_total.ne.basisdata%nr_total) stop 'Size problem'
        if (nc_total.ne.basisdata%nc_total) stop 'Size problem'
        if (mat%dimen.ne.basisdata%dimen) stop 'Dimension problem'

        nr_u = basisdata%nrows(1); nc_u = basisdata%ncols(1); nnz_u = basisdata%nnzs(1)
        nr_v = basisdata%nrows(2); nc_v = basisdata%ncols(2); nnz_v = basisdata%nnzs(2)
        allocate(indi_u(nr_u+1), indj_u(nnz_u), data_W_u(nnz_u, 4))
        indi_u = basisdata%indi(1, 1:nr_u+1); indj_u = basisdata%indj(1, 1:nnz_u)
        data_W_u = basisdata%data_bw(1, 1:nnz_u, 3:6)
        allocate(indi_v(nr_v+1), indj_v(nnz_v), data_W_v(nnz_v, 4))
        indi_v = basisdata%indi(2, 1:nr_v+1); indj_v = basisdata%indj(2, 1:nnz_v)
        data_W_v = basisdata%data_bw(2, 1:nnz_v, 3:6)
        if (basisdata%dimen.eq.3) then
            nr_w = basisdata%nrows(3); nc_w = basisdata%ncols(3); nnz_w = basisdata%nnzs(3)
            allocate(indi_w(nr_w+1), indj_w(nnz_w), data_W_w(nnz_w, 4))
            indi_w = basisdata%indi(3, 1:nr_w+1); indj_w = basisdata%indj(3, 1:nnz_w)
            data_W_w = basisdata%data_bw(3, 1:nnz_w, 3:6)
        end if

        do i = 1, nc_total
            call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
            t1(:, :, i) = matmul(mat%invJ(:, :, i), Tstress)*mat%detJ(i)
        end do
        
        array_out = 0.d0
        do i = 1, basisdata%dimen
            do k = 1, basisdata%dimen
                alpha = 1; alpha(k) = 2
                zeta  = 1 + (alpha - 1)*2
                if (basisdata%dimen.eq.2) then
                    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        t1(k, i, :), t2)
                else if (basisdata%dimen.eq.3) then
                    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                        nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                        nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                        nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                        t1(k, i, :), t2)
                end if
                array_out(i, :) = array_out(i, :) + t2
            end do
        end do

    end subroutine intforce

end module matrixfreeplasticity

module plasticitysolver

    use matrixfreeplasticity
    use structured_data

    type cgsolver
        logical :: withdiag = .true., applyfd = .true.
        integer :: matrixfreetype = 2
        type(basis_data), pointer :: globsyst=>null()
        type(reduced_system), dimension(:), pointer :: redsyst=>null()
    end type cgsolver

contains

    subroutine matrixfree_matvec(solv, mat, nr_total, array_in, array_out)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat)  :: mat
        integer, intent(in) :: nr_total
        double precision, intent(in) :: array_in
        dimension :: array_in(solv%globsyst%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(solv%globsyst%dimen, nr_total)

        if (solv%matrixfreetype.eq.1) then
            call mf_tu_tv(mat, solv%globsyst, nr_total, array_in, array_out)
        else if (solv%matrixfreetype.eq.2) then
            call mf_gradtu_gradtv(mat, solv%globsyst, nr_total, array_in, array_out)
        else if (solv%matrixfreetype.eq.3) then
            call mf_tutv_gradtugradtv(mat,  solv%globsyst, nr_total, array_in, array_out)
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
        type(reduced_system), dimension(:), target :: redsyst
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
        dimension :: array_in(solv%globsyst%dimen, nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(solv%globsyst%dimen, nr_total)

        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w, i
        double precision, allocatable, dimension(:) :: tmp1, tmp2

        if (.not.solv%applyfd) then
            array_out = array_in
            return
        end if

        array_out = 0.d0
        do i = 1, solv%globsyst%dimen

            nr_u = solv%redsyst(i)%basisdata%nrows(1)
            nr_v = solv%redsyst(i)%basisdata%nrows(2)
            nr_w = 1
            if (solv%globsyst%dimen.eq.3) nr_w = solv%redsyst(i)%basisdata%nrows(3)

            allocate(tmp1(nr_u*nr_v*nr_w))
            if (solv%globsyst%dimen.eq.2) then
                call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, &
                            transpose(solv%redsyst(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                            transpose(solv%redsyst(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                            array_in(i, solv%redsyst(i)%dof), tmp1)
            else if (solv%globsyst%dimen.eq.3) then                
                call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            transpose(solv%redsyst(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                            transpose(solv%redsyst(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                            transpose(solv%redsyst(i)%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                            array_in(i, solv%redsyst(i)%dof), tmp1)
            end if           

            if (solv%withdiag) then
                tmp1 = tmp1/solv%redsyst(i)%diageigval_sp
            end if

            allocate(tmp2(nr_u*nr_v*nr_w))
            if (solv%globsyst%dimen.eq.2) then
                call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, &
                            solv%redsyst(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                            solv%redsyst(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), tmp1, tmp2)
            else if (solv%globsyst%dimen.eq.3) then                
                call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            solv%redsyst(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                            solv%redsyst(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                            solv%redsyst(i)%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), tmp1, tmp2)
            end if

            array_out(i, solv%redsyst(i)%dof) = tmp2
            deallocate(tmp1, tmp2)
        end do

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, size_inout, array_inout)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: size_inout
        double precision, intent(inout) :: array_inout(solv%globsyst%dimen, size_inout)

        ! Local data
        ! ----------
        integer :: i
        do i = 1, solv%globsyst%dimen
            call set2zero(solv%redsyst(i), size_inout, array_inout(i, :))
        end do
    end subroutine clear_dirichlet

    subroutine PBiCGSTAB(solv, mat, nr_total, iterations, threshold, b, x, residual)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, iterations
        double precision, intent(in) :: threshold, b
        dimension :: b(solv%globsyst%dimen, nr_total)
        
        double precision, intent(out) :: x, residual
        dimension :: x(solv%globsyst%dimen, nr_total), residual(iterations+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(solv%globsyst%dimen, nr_total), rhat(solv%globsyst%dimen, nr_total), &
                        p(solv%globsyst%dimen, nr_total), s(solv%globsyst%dimen, nr_total), &
                        ptilde(solv%globsyst%dimen, nr_total), Aptilde(solv%globsyst%dimen, nr_total), &
                        Astilde(solv%globsyst%dimen, nr_total), stilde(solv%globsyst%dimen, nr_total)
        integer :: k

        x = 0.d0; r = b; residual = 0.d0
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        call block_dot_product(solv%globsyst%dimen, nr_total, r, rhat, rsold)
        normb = norm2(r) 
        if (normb.le.1.d-14) return
        residual(1) = 1.d0
        do k = 1, iterations
            call applyfastdiag(solv, nr_total, p, ptilde) 
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_matvec(solv, mat, nr_total, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            call block_dot_product(solv%globsyst%dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_matvec(solv, mat, nr_total, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
            
            call block_dot_product(solv%globsyst%dimen, nr_total, Astilde, s, prod)
            call block_dot_product(solv%globsyst%dimen, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            if (norm2(r).le.max(threshold*normb, 1.d-14)) exit
            residual(k+1) = norm2(r)/normb
            call block_dot_product(solv%globsyst%dimen, nr_total, r, rhat, rsnew)
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
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, iterations        
        logical, intent(in) :: ishigher
        double precision, intent(in) :: threshold
        
        double precision, intent(out) :: eigenvec, eigenval
        dimension :: eigenvec(solv%globsyst%dimen, nr_total)

        ! Local data
        ! ----------
        integer :: k, j, ii
        double precision, dimension(sizemat) :: ll
        double precision, dimension(sizemat, solv%globsyst%dimen, nr_total) :: RM1, RM2, RM3
        double precision, dimension(sizemat, sizemat) :: AA1, BB1, qq
        double precision, dimension(solv%globsyst%dimen, nr_total) :: u, v, g, gtil, p, tmp, Mg, Mgtil
        double precision :: q, gnorm, gnorm0, delta, prod1, prod2

        call random_number(eigenvec)
        call clear_dirichlet(solv, nr_total, eigenvec)

        call mf_tu_tv(mat, solv%globsyst, nr_total, eigenvec, u)
        call clear_dirichlet(solv, nr_total, u)
        
        call block_dot_product(solv%globsyst%dimen, nr_total, eigenvec, u, prod1)
        q = sqrt(prod1)
        eigenvec = eigenvec/q; u = u/q
        call mf_gradtu_gradtv(mat, solv%globsyst, nr_total, eigenvec, v)
        call clear_dirichlet(solv, nr_total, v)
        
        call block_dot_product(solv%globsyst%dimen, nr_total, eigenvec, v, eigenval)
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
                call mf_tu_tv(mat, solv%globsyst, nr_total, g, Mg)
                call mf_tu_tv(mat, solv%globsyst, nr_total, gtil, Mgtil)
                call block_dot_product(solv%globsyst%dimen, nr_total, g, Mg, prod1)
                call block_dot_product(solv%globsyst%dimen, nr_total, gtil, Mgtil, prod2)
                p = -g + prod1/prod2*p
            end if

            RM1(1, :, :) = eigenvec; RM1(2, :, :) = p
            RM2(1, :, :) = v; RM3(1, :, :) = u;
            call mf_gradtu_gradtv(mat, solv%globsyst, nr_total, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM2(2, :, :) = tmp

            call mf_tu_tv(mat, solv%globsyst, nr_total, p, tmp)
            call clear_dirichlet(solv, nr_total, tmp)
            RM3(2, :, :) = tmp
            
            call rayleigh_submatrix2(sizemat, solv%globsyst%dimen, nr_total, RM1, RM2, AA1)
            call rayleigh_submatrix2(sizemat, solv%globsyst%dimen, nr_total, RM1, RM3, BB1)
            
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
            call mf_tu_tv(mat, solv%globsyst, nr_total, eigenvec, u)
            call clear_dirichlet(solv, nr_total, u)
            
            call block_dot_product(solv%globsyst%dimen, nr_total, eigenvec, u, prod1)
            q = sqrt(prod1)
            eigenvec = eigenvec/q; u = u/q
            call mf_gradtu_gradtv(mat, solv%globsyst, nr_total, eigenvec, v)
            call clear_dirichlet(solv, nr_total, v)            
            gnorm = norm2(g)
        end do

    end subroutine RQMIN

end module plasticitysolver
