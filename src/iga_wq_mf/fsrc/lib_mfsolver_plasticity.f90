subroutine block_dot_product(nm, nr, A, B, result)
    !! Computes dot product of A and B. Both are actually vectors arranged following each dimension
    !! Vector A is composed of [Au, Av, Aw] and B of [Bu, Bv, Bw]. 
    !! Dot product A.B = Au.Bu + Av.Bv + Aw.Bw 

    implicit none
    ! Input/ output data
    ! ------------------
    integer, intent(in) :: nm, nr
    double precision, intent(in) :: A, B
    dimension :: A(nm, nr), B(nm, nr)

    double precision :: result

    ! Local data
    ! ----------
    integer :: i
    double precision :: tmp

    result = 0.d0
    do i = 1, nm 
        tmp = dot_product(A(i, :), B(i, :))
        result = result + tmp
    end do

end subroutine block_dot_product

module matrixfreeplasticity

    implicit none
    type :: mecamat
    
        integer :: dimen, nvoigt
        logical :: isLumped = .false., isElastic = .true.
        double precision :: scalars(2) = (/1.d0, 1.d0/)
        double precision, dimension(:), pointer :: detJ=>null(), Mprop=>null(), Hprop=>null()
        double precision, dimension(:, :), pointer :: CepArgs=>null(), NN=>null(), BB=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn, JJbb
        double precision, dimension(:, :), allocatable :: Smean
        double precision, dimension(:), allocatable :: Mmean
        integer :: ncols_sp
    
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

        mat%Mprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_massprop

    subroutine setup_thmchcoupledprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Hprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_thmchcoupledprop

    subroutine setup_mechanicalArguments(mat, nbrows, mechArgs)
        !! Points to data of the mechanical behavior 

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
            mat%CepArgs => mechArgs(:2, :)
        else
            mat%CepArgs => mechArgs(:4, :)    
            mat%NN      => mechArgs(5:5+mat%nvoigt, :)
            mat%BB      => mechArgs(5+mat%nvoigt:5+2*mat%nvoigt, :)

            if (.not.allocated(mat%JJnn)) allocate(mat%JJnn(mat%dimen, mat%dimen, mat%ncols_sp))
            if (.not.allocated(mat%JJbb)) allocate(mat%JJbb(mat%dimen, mat%dimen, mat%ncols_sp))
            mat%JJnn = 0.d0; mat%JJbb = 0.d0
            if (all(abs(mat%NN).le.threshold)) return
            do i = 1, mat%ncols_sp
                call array2symtensor(mat%dimen, mat%nvoigt, mat%NN(:, i), tensor)
                mat%JJnn(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
                call array2symtensor(mat%dimen, mat%nvoigt, mat%BB(:, i), tensor)
                mat%JJbb(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
            end do
        end if

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

        if (.not.associated(mat%Mprop)) then
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
                    if (mat%isElastic.eqv..false.) then
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
                    tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
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
                    if (mat%isElastic.eqv..false.) then
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
                    tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
                    do j = 1, mat%dimen
                        CC(j, gp) = tensor(j, j)
                    end do
                    CC(mat%dimen+1, gp) = mat%Mprop(gp)*mat%detJ(gp)
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

    subroutine compute_mean(mat, nclist)
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
                if (mat%isElastic.eqv..false.) then
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
                
                tensor = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
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

        if (associated(mat%Mprop)) then
            allocate(Mcoefs(size(sample)))
            do c = 1, size(sample)
                gp = sample(c)
                Mcoefs(c) = mat%Mprop(gp)*mat%detJ(gp)
            end do
            if (mat%dimen.eq.2) then
                call trapezoidal_rule_2d(3, 3, Mcoefs, tensor(1, 1))
            else if (mat%dimen.eq.3) then
                call trapezoidal_rule_3d(3, 3, 3, Mcoefs, tensor(1, 1))
            end if 
            mat%Mmean = tensor(1, 1)
        end if

    end subroutine compute_mean

    subroutine mf_tu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i
        double precision :: tmp_in, array_tmp, array_tmp2
        dimension :: tmp_in(nr_total), array_tmp(nc_total), array_tmp2(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        do i = 1, dimen
            tmp_in = array_in(i, :); if (mat%isLumped) tmp_in = 1.d0

            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                tmp_in, array_tmp)

            array_tmp = array_tmp*mat%Mprop*mat%detJ

            call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                array_tmp, array_tmp2)

            array_out(i, :) = array_tmp2; if (mat%isLumped) array_out(i, :) = array_tmp2*array_in(i, :)

        end do

    end subroutine mf_tu_tv_2d

    subroutine mf_tu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i
        double precision :: tmp_in, array_tmp, array_tmp2
        dimension :: tmp_in(nr_total), array_tmp(nc_total), array_tmp2(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'

        do i = 1, dimen
            tmp_in = array_in(i, :); if (mat%isLumped) tmp_in = 1.d0

            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                                tmp_in, array_tmp)

            array_tmp = array_tmp*mat%Mprop*mat%detJ

            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                array_tmp, array_tmp2)
            array_out(i, :) = array_tmp2; if (mat%isLumped) array_out(i, :) = array_tmp2*array_in(i, :)
        end do
            
    end subroutine mf_tu_tv_3d

    subroutine mf_gradtu_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 2D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, j, k, l, m, alpha, beta, zeta, nbCepArgs = 4
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision, allocatable, dimension(:) :: t4, t5, t6
        double precision, allocatable, dimension(:, :) :: kt1
        double precision :: t1, t2, t3, t7, t8, t9
        dimension :: t1(nc_total), t2(nc_total), t3(nc_total), &
                        t7(nc_total), t8(nr_total), t9(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        array_out = 0.d0

        if (mat%isElastic) then 
            nbCepArgs = 2
        else
            allocate(t4(nc_total), t5(nc_total), t6(nc_total))
        end if
        allocate(kt1(nbCepArgs, nc_total))

        do j = 1, dimen
            do m = 1, dimen
                beta = 1; beta(m) = 2
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in(j, :), t1) 

                do k = 1, nbCepArgs
                    kt1(k, :) = mat%CepArgs(k, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(m, j, :)
                if (mat%isElastic.eqv..false.) then
                    t4 = kt1(3, :)*mat%JJnn(m, j, :)
                    t5 = kt1(4, :)*mat%JJbb(m, j, :)
                    t6 = kt1(4, :)*mat%JJnn(m, j, :)
                end if


                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(m, i, :)
                    t9 = 0.d0

                    do l = 1, dimen
                        alpha = 1; alpha(l) = 2
                        zeta  = beta + (alpha - 1)*2
                        
                        t7 = t2*mat%invJ(l, i, :) + t3*mat%invJ(l, j, :)
                        if (mat%isElastic.eqv..false.) then
                            t7 = t7 + t4*mat%JJnn(l, i, :) - t5*mat%JJnn(l, i, :) + t6*mat%JJbb(l, i, :)
                        end if
                        if (i.eq.j) t7 = t7 + kt1(2, :)*mat%JJjj(l, m, :)
                        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                t7, t8)
                        t9 = t9 + t8
                    end do
            
                    array_out(i, :) = array_out(i, :) + t9
                end do
            end do
        end do
            
    end subroutine mf_gradtu_gradtv_2d

    subroutine mf_gradtu_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, j, k, l, m, alpha, beta, zeta, nbCepArgs = 4
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision, allocatable, dimension(:) :: t4, t5, t6
        double precision, allocatable, dimension(:, :) :: kt1
        double precision :: t1, t2, t3, t7, t8, t9
        dimension :: t1(nc_total), t2(nc_total), t3(nc_total), &
                        t7(nc_total), t8(nc_total), t9(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        array_out = 0.d0
        
        if (mat%isElastic) then 
            nbCepArgs = 2
        else
            allocate(t4(nc_total), t5(nc_total), t6(nc_total))
        end if
        allocate(kt1(nbCepArgs, nc_total))

        do j = 1, dimen
            do m = 1, dimen
                beta = 1; beta(m) = 2
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in(j, :), t1) 

                do k = 1, nbCepArgs
                    kt1(k, :) = mat%CepArgs(k, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(m, j, :)
                if (mat%isElastic.eqv..false.) then
                    t4 = kt1(3, :)*mat%JJnn(m, j, :)
                    t5 = kt1(4, :)*mat%JJbb(m, j, :)
                    t6 = kt1(4, :)*mat%JJnn(m, j, :)
                end if

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(m, i, :)
                    t9 = 0.d0

                    do l = 1, dimen
                        alpha = 1; alpha(l) = 2
                        zeta  = beta + (alpha - 1)*2
                        
                        t7 = t2*mat%invJ(l, i, :) + t3*mat%invJ(l, j, :)
                        if (mat%isElastic.eqv..false.) then
                            t7 = t7 + t4*mat%JJnn(l, i, :) - t5*mat%JJnn(l, i, :) + t6*mat%JJbb(l, i, :)
                        end if
                        if (i.eq.j) t7 = t7 + kt1(2, :)*mat%JJjj(l, m, :)
                        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t8, t9)
                        t9 = t9 + t8
                    end do

                    array_out(i, :) = array_out(i, :) + t9
                end do
            end do
        end do
            
    end subroutine mf_gradtu_gradtv_3d

    subroutine mf_tutv_gradtugradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 2D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(dimen, nr_total)

        call mf_tu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradtu_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_tmp)

        array_out = array_out + mat%scalars(2)*array_tmp

    end subroutine mf_tutv_gradtugradtv_2d

    subroutine mf_tutv_gradtugradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(dimen, nr_total)

        call mf_tu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradtu_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, array_in, array_tmp)
        
        array_out = array_out + mat%scalars(2)*array_tmp

    end subroutine mf_tutv_gradtugradtv_3d

    subroutine mf_u_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)
        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        array_out = 0.d0     
        do l = 1, dimen
            beta = 1; beta(l) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in, t1) 
            t1 = t1*mat%Hprop*mat%detJ
            do i = 1, dimen
                alpha = 1; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(l, i, :)
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), t2, t3)
                array_out(i, :) = array_out(i, :) + t3
            end do
        end do
            
    end subroutine mf_u_gradtv_2d

    subroutine mf_u_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
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
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        array_out = 0.d0     
        do l = 1, dimen
            beta = 1; beta(l) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in, t1) 
            t1 = t1*mat%Hprop*mat%detJ
            do i = 1, dimen
                alpha = 1; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(l, i, :)
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                            nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t2, t3)
                array_out(i, :) = array_out(i, :) + t3
            end do
        end do
            
    end subroutine mf_u_gradtv_3d

    subroutine intforce_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, stress, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ----------
        double precision :: Tstress, t1, t2
        dimension :: Tstress(dimen, dimen), t1(dimen, dimen, nc_total), t2(nr_total)
        integer :: i, k, alpha(dimen), zeta(dimen)
        
        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        
        do i = 1, nc_total
            call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
            t1(:, :, i) = matmul(mat%invJ(:, :, i), Tstress)*mat%detJ(i)
        end do
        
        array_out = 0.d0
        do i = 1, dimen
            do k = 1, dimen
                alpha = 1; alpha(k) = 2
                zeta  = 1 + (alpha - 1)*2

                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                t1(k, i, :), t2)
                array_out(i, :) = array_out(i, :) + t2
            end do
        end do

    end subroutine intforce_2d

    subroutine intforce_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, stress, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ----------
        double precision :: Tstress, t1, t2
        dimension :: Tstress(dimen, dimen), t1(dimen, dimen, nc_total), t2(nr_total)
        integer :: i, k, alpha(dimen), zeta(dimen)
        
        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        
        do i = 1, nc_total
            call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
            t1(:, :, i) = matmul(mat%invJ(:, :, i), Tstress)*mat%detJ(i)
        end do
        
        array_out = 0.d0
        do i = 1, dimen
            do k = 1, dimen
                alpha = 1; alpha(k) = 2
                zeta  = 1 + (alpha - 1)*2

                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                t1(k, i, :), t2)
                array_out(i, :) = array_out(i, :) + t2
            end do
        end do

    end subroutine intforce_3d

end module matrixfreeplasticity

module plasticitysolver2

    use matrixfreeplasticity
    use datastructure
    
    type cgsolver
        logical :: withdiag = .true., applyfd = .true.
        integer :: matrixfreetype = 2, dimen = 2
        type(structure) :: disp_struct(2)
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
        type(mecamat) :: mat
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
        dimension :: array_in(solv%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        if (solv%matrixfreetype.eq.1) then

            call mf_tu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)


        else if (solv%matrixfreetype.eq.2) then

            call mf_gradtu_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                                data_W_u, data_W_v, array_in, array_out)

        else if (solv%matrixfreetype.eq.3) then

            call mf_tutv_gradtugradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                                data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                                data_W_u, data_W_v, array_in, array_out)

        else
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, table, mean)

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
        logical, intent(in) :: table
        dimension :: table(solv%dimen, 2, solv%dimen)
        double precision, intent(in) :: mean
        dimension :: mean(solv%dimen, solv%dimen)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call init_2datastructure(solv%disp_struct(i), nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                    data_B_u, data_B_v, data_W_u, data_W_v)
            call update_datastructure(solv%disp_struct(i), solv%dimen, table(:, :, i))
            if (.not.solv%applyfd) cycle
            call space_eigendecomposition(solv%disp_struct(i), solv%dimen, mean(i, :))
        end do
        
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
        dimension :: array_in(solv%dimen, nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)
    
        ! Local data
        ! ----------
        integer :: nr_u, nr_v, i
        double precision, allocatable, dimension(:) :: tmp, tmp2

        if (.not.solv%applyfd) then
            array_out = array_in
            return
        end if

        array_out  = 0.d0
        do i = 1, solv%dimen
            ! Compute (Uv x Uu)'.array_in
            nr_u = solv%disp_struct(i)%nrows(1)
            nr_v = solv%disp_struct(i)%nrows(2)
            allocate(tmp(nr_u*nr_v))
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, transpose(solv%disp_struct(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                    transpose(solv%disp_struct(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), array_in(i, solv%disp_struct(i)%dof), tmp)
            
            if (solv%withdiag) then
                tmp = tmp/solv%disp_struct(i)%diageigval_sp
            end if

            ! Compute (Uv x Uu).array_tmp
            allocate(tmp2(nr_u*nr_v))
            call sumfacto2d_dM(nr_u, nr_u, nr_v, nr_v, solv%disp_struct(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                    solv%disp_struct(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), tmp, tmp2)
            array_out(i, solv%disp_struct(i)%dof) = tmp2            
            deallocate(tmp, tmp2)
        end do

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(solv%dimen, nnz)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call set2zero(solv%disp_struct(i), nnz, array(i, :))
        end do

    end subroutine clear_dirichlet

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, &
                        nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
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
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), & 
                        s(solv%dimen, nr_total), ptilde(solv%dimen, nr_total), Aptilde(solv%dimen, nr_total), &
                        Astilde(solv%dimen, nr_total), stilde(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b; resPCG = 0.d0
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb = norm2(r)
        if (normb.le.1.d-14) return
        resPCG(1) = 1.d0

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde)
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            call block_dot_product(solv%dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, &
                        nnz_u, nnz_v, indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)            
            call block_dot_product(solv%dimen, nr_total, Astilde, s, prod)
            call block_dot_product(solv%dimen, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            if (norm2(r).le.max(threshold*normb, 1.d-14)) exit
            resPCG(iter+1) = norm2(r)/normb
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module plasticitysolver2

module plasticitysolver3

    use matrixfreeplasticity
    use datastructure

    type cgsolver
        logical :: withdiag = .true., applyfd = .true.
        integer :: matrixfreetype = 2, dimen = 3
        type(structure) :: disp_struct(3)
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
        type(mecamat)  :: mat
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
        dimension :: array_in(solv%dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        if (solv%matrixfreetype.eq.1) then

            call mf_tu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        else if (solv%matrixfreetype.eq.2) then

            call mf_gradtu_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, array_in, array_out)
        
        else if (solv%matrixfreetype.eq.3) then

            call mf_tutv_gradtugradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, array_in, array_out)
                            
        else
            stop 'Not coded'
        end if

    end subroutine matrixfree_spMdV

    subroutine initializefastdiag(solv, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, mean)

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
        dimension :: table(solv%dimen, 2, solv%dimen)
        double precision, intent(in) :: mean
        dimension :: mean(solv%dimen, solv%dimen)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call init_3datastructure(solv%disp_struct(i), nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)
            call update_datastructure(solv%disp_struct(i), solv%dimen, table(:, :, i))
            if (.not.solv%applyfd) cycle
            call space_eigendecomposition(solv%disp_struct(i), solv%dimen, mean(i, :))
        end do
        
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
        dimension :: array_in(solv%dimen, nr_total)
    
        double precision, intent(out) :: array_out
        dimension :: array_out(solv%dimen, nr_total)

        ! Local data
        ! ----------
        integer :: nr_u, nr_v, nr_w, i
        double precision, allocatable, dimension(:) :: tmp, tmp2

        if (.not.solv%applyfd) then
            array_out = array_in
            return
        end if

        array_out = 0.d0
        do i = 1, solv%dimen
            ! Compute (Uw x Uv x Uu)'.array_in
            nr_u = solv%disp_struct(i)%nrows(1)
            nr_v = solv%disp_struct(i)%nrows(2)
            nr_w = solv%disp_struct(i)%nrows(3)
            allocate(tmp(nr_u*nr_v*nr_w))

            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            transpose(solv%disp_struct(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u)), &
                            transpose(solv%disp_struct(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v)), &
                            transpose(solv%disp_struct(i)%eigvec_sp_dir(3, 1:nr_w, 1:nr_w)), &
                            array_in(i, solv%disp_struct(i)%dof), tmp)

            if (solv%withdiag) then
                tmp = tmp/solv%disp_struct(i)%diageigval_sp
            end if

            ! Compute (Uw x Uv x Uu).array_tmp
            allocate(tmp2(nr_u*nr_v*nr_w))
            call sumfacto3d_dM(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            solv%disp_struct(i)%eigvec_sp_dir(1, 1:nr_u, 1:nr_u), &
                            solv%disp_struct(i)%eigvec_sp_dir(2, 1:nr_v, 1:nr_v), &
                            solv%disp_struct(i)%eigvec_sp_dir(3, 1:nr_w, 1:nr_w), tmp, tmp2)
            array_out(i, solv%disp_struct(i)%dof) = tmp2
            deallocate(tmp, tmp2)
        end do

    end subroutine applyfastdiag

    subroutine clear_dirichlet(solv, nnz, array)
        implicit none
        ! Input / output  data 
        !---------------------
        type(cgsolver) :: solv
        integer, intent(in) :: nnz
        double precision, intent(inout) :: array(solv%dimen, nnz)

        ! Local data
        ! ----------
        integer :: i

        do i = 1, solv%dimen
            call set2zero(solv%disp_struct(i), nnz, array(i, :))
        end do

    end subroutine clear_dirichlet

    subroutine PBiCGSTAB(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, nbIterPCG, threshold, b, x, resPCG)

        implicit none
        ! Input / output data
        ! -------------------
        type(cgsolver) :: solv
        type(mecamat) :: mat
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
        dimension :: b(solv%dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(solv%dimen, nr_total), resPCG(nbIterPCG+1)

        ! Local data
        ! -----------
        double precision :: prod, prod2, rsold, rsnew, alpha, omega, beta, normb
        double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
        dimension ::    r(solv%dimen, nr_total), rhat(solv%dimen, nr_total), p(solv%dimen, nr_total), & 
                        s(solv%dimen, nr_total), ptilde(solv%dimen, nr_total), Aptilde(solv%dimen, nr_total), &
                        Astilde(solv%dimen, nr_total), stilde(solv%dimen, nr_total)
        integer :: iter

        x = 0.d0; r = b; resPCG = 0.d0
        call clear_dirichlet(solv, nr_total, r)
        rhat = r; p = r
        call block_dot_product(solv%dimen, nr_total, r, rhat, rsold)
        normb = norm2(r) 
        if (normb.le.1.d-14) return
        resPCG(1) = 1.d0

        do iter = 1, nbIterPCG
            call applyfastdiag(solv, nr_total, p, ptilde) 
            call clear_dirichlet(solv, nr_total, ptilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call clear_dirichlet(solv, nr_total, Aptilde)
            call block_dot_product(solv%dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde
            
            call applyfastdiag(solv, nr_total, s, stilde)
            call clear_dirichlet(solv, nr_total, stilde)
            call matrixfree_spMdV(solv, mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call clear_dirichlet(solv, nr_total, Astilde)
            
            call block_dot_product(solv%dimen, nr_total, Astilde, s, prod)
            call block_dot_product(solv%dimen, nr_total, Astilde, Astilde, prod2)

            omega = prod/prod2
            x = x + alpha*ptilde + omega*stilde
            r = s - omega*Astilde    
            
            if (norm2(r).le.max(threshold*normb, 1.d-14)) exit
            resPCG(iter+1) = norm2(r)/normb
            call block_dot_product(solv%dimen, nr_total, r, rhat, rsnew)
            beta = (alpha/omega)*(rsnew/rsold)
            p = r + beta*(p - omega*Aptilde)
            rsold = rsnew
        end do

    end subroutine PBiCGSTAB

end module plasticitysolver3
