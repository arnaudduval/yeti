module matrixfreeplasticity

    implicit none
    type :: mecamat
    
        integer :: dimen, nvoigt
        double precision :: elasticmodulus, poissonratio, elasticlimit
        double precision, dimension(:), pointer :: detJ=>null()
        double precision, dimension(:, :), pointer :: CepArgs=>null(), NN=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn
        double precision, dimension(:, :), allocatable :: mean
        double precision :: lambda, mu, bulk
        integer :: ncols_sp
    
    end type mecamat

contains

    subroutine initialize_mecamat(mat, dimen, elasticmodulus, poissonratio, elasticlimit)
        !! Creates a material using combined isotropic/kinematic hardening theory 

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: dimen
        double precision, intent(in) :: elasticmodulus, poissonratio, elasticlimit
        type(mecamat) :: mat

        mat%dimen  = dimen
        mat%nvoigt = dimen*(dimen+1)/2

        ! Compute mechanical properties
        mat%elasticmodulus = elasticmodulus
        mat%poissonratio   = poissonratio
        mat%elasticlimit   = elasticlimit
        mat%lambda = poissonratio*elasticmodulus/((1+poissonratio)*(1-2*poissonratio))
        mat%mu     = elasticmodulus/(2*(1+poissonratio))
        mat%bulk   = mat%lambda + 2.d0/3.d0*mat%mu
        allocate(mat%mean(mat%dimen, mat%dimen))
        mat%mean   = 1.d0

    end subroutine initialize_mecamat

    subroutine setup_geometry(mat, nnz, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(mat%dimen, mat%dimen, nnz), detJ(nnz)

        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nnz

    end subroutine setup_geometry

    subroutine setup_jacobiennormal(mat, MechArgs)
        !! Points to data of the mechanical behavior 

        implicit none
        ! Input / output data
        ! -------------------
        double precision, parameter :: threshold = 1.d-12
        type(mecamat) :: mat
        double precision, target, intent(in) :: MechArgs
        dimension :: MechArgs(mat%nvoigt+3, mat%ncols_sp)

        ! Local data
        ! ----------
        double precision :: TNN(mat%dimen, mat%dimen)
        integer :: i

        mat%CepArgs => MechArgs(:3, :)
        mat%NN      => MechArgs(4:, :)
        if (.not.allocated(mat%JJnn)) allocate(mat%JJnn(mat%dimen, mat%dimen, mat%ncols_sp))
        mat%JJnn = 0.d0
        if (all(abs(mat%NN).lt.threshold)) return
        do i = 1,  mat%ncols_sp
            call array2symtensor(mat%dimen, mat%nvoigt, mat%NN(:, i), TNN)
            mat%JJnn(:, :, i) = matmul(mat%invJ(:, :, i), TNN)
        end do
    end subroutine setup_jacobiennormal

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

    subroutine compute_separationvariables_ijblock(i, j, dimen, nvoigt, nc_list, mat, Mcoefs, Kcoefs)
        
        use separatevariables
        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: i, j, dimen, nvoigt, nc_list
        dimension :: nc_list(dimen)

        double precision, intent(out) :: Mcoefs(dimen, maxval(nc_list)), Kcoefs(dimen, maxval(nc_list))

        ! Local data
        ! ----------
        type(sepoperator) :: oper
        logical :: update(dimen)
        integer :: k, l, m, genpos
        double precision :: DD, coefs, NN, TNN
        dimension :: DD(dimen, dimen), coefs(dimen, dimen, mat%ncols_sp), NN(nvoigt), TNN(dimen, dimen)

        do genpos = 1, mat%ncols_sp
            NN = mat%NN(:, genpos)
            call array2symtensor(dimen, size(NN), NN, TNN)

            DD = 0.d0
            DD(i, j) = DD(i, j) + mat%CepArgs(1, genpos)
            DD(j, i) = DD(j, i) + mat%CepArgs(2, genpos)
            if (i.eq.j) then
                do k = 1, dimen
                    DD(k, k) = DD(k, k) + mat%CepArgs(2, genpos)
                end do
            end if
            do l = 1, dimen
                do m = 1, dimen
                    DD(l, m) = DD(l, m) + mat%CepArgs(3, genpos)*TNN(i, l)*TNN(j, m)
                end do
            end do
            coefs(:, :, genpos) = matmul(mat%invJ(:, :, genpos), matmul(DD, transpose(mat%invJ(:, :, genpos))))
        end do

        update = .true.
        call initialize_operator(oper, dimen, nc_list, update)
        if (dimen.eq.2) then
            call separatevariables_2d(oper, coefs)
        else if (dimen.eq.3) then
            call separatevariables_3d(oper, coefs)
        end if
        Mcoefs = oper%Mcoefs; Kcoefs = oper%Kcoefs

    end subroutine compute_separationvariables_ijblock

    subroutine compute_mean_ijblock(i, j, dimen, nvoigt, mat, samplesize, sample, mean)
        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: i, j, dimen, nvoigt, samplesize, sample
        dimension :: sample(samplesize)

        double precision, intent(out) :: mean(dimen)

        ! Local data
        ! ----------
        integer :: k, l, m, n, c, genpos
        double precision :: DD, coefs, NN, TNN
        dimension :: DD(dimen, dimen), coefs(dimen, dimen, samplesize), NN(nvoigt), TNN(dimen, dimen)

        do c = 1, samplesize
            genpos = sample(c)
            NN = mat%NN(:, genpos)
            call array2symtensor(dimen, size(NN), NN, TNN)

            DD = 0.d0
            DD(i, j) = DD(i, j) + mat%CepArgs(1, genpos)
            DD(j, i) = DD(j, i) + mat%CepArgs(2, genpos)
            if (i.eq.j) then
                do k = 1, dimen
                    DD(k, k) = DD(k, k) + mat%CepArgs(2, genpos)
                end do
            end if
            do l = 1, dimen
                do m = 1, dimen
                    DD(l, m) = DD(l, m) + mat%CepArgs(3, genpos)*TNN(i, l)*TNN(j, m)
                end do
            end do
            coefs(:, :, c) = matmul(mat%invJ(:, :, genpos), matmul(DD, transpose(mat%invJ(:, :, genpos))))
        end do

        do n = 1, dimen
            if (dimen.eq.2) then
                call trapezoidal_rule_2d(3, 3, coefs(n, n, :), mean(n))
            else if (dimen.eq.3) then
                call trapezoidal_rule_3d(3, 3, 3, coefs(n, n, :), mean(n))
            end if
        end do   

    end subroutine compute_mean_ijblock

    subroutine compute_mean_diagblocks(mat, dimen, nvoigt, nclist)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: dimen, nvoigt, nclist
        dimension :: nclist(dimen)

        ! Local data
        ! ----------
        integer :: i, j, k, c, genpos, pos, ind(3)
        integer, dimension(:), allocatable :: sample
        integer, dimension(:, :), allocatable :: indlist
        double precision :: mean(dimen)
        
        if (product(nclist).ne.mat%ncols_sp) stop 'Wrong dimensions'
        allocate(indlist(dimen, 3), sample(3**dimen))
        do i = 1, dimen 
            pos = int((nclist(i) + 1)/2); ind = (/1, pos, nclist(i)/)
            indlist(i, :) = ind
        end do
    
        ! Select a set of coefficients
        c = 1
        if (dimen.eq.2) then
            do j = 1, 3
                do i = 1, 3
                    genpos = indlist(1, i) + (indlist(2, j) - 1)*nclist(1)
                    sample(c) = genpos
                    c = c + 1
                end do
            end do
        else if (dimen.eq.3) then
            do k = 1, 3
                do j = 1, 3
                    do i = 1, 3
                        genpos = indlist(1, i) + (indlist(2, j) - 1)*nclist(1) + (indlist(3, k) - 1)*nclist(1)*nclist(2)
                        sample(c) = genpos
                        c = c + 1
                    end do
                end do
            end do
        else
            stop 'Not possible compute mean diagonal blocks'
        end if

        do i = 1, dimen
            call compute_mean_ijblock(i, i, dimen, nvoigt, mat, size(sample), sample, mean)
            mat%mean(i, :) = mean
        end do

    end subroutine compute_mean_diagblocks

    subroutine mf_stiffness_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 3D where S is stiffness matrix
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
        integer :: i, j, k, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: kt1, t1, t2, t3, t4, t5, t6, t7
        dimension ::    kt1(3, nc_total), t1(nc_total), t2(nc_total), t3(nc_total), &
                        t4(nc_total), t5(nc_total), t6(nr_total), t7(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Number of rows not equal'
        array_out = 0.d0       
        do j = 1, dimen
            do l = 1, dimen
                beta = 1; beta(l) = 2
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in(j, :), t1) 

                do i = 1, 3
                    kt1(i, :) = mat%CepArgs(i, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(l, j, :)
                t4 = kt1(3, :)*mat%JJnn(l, j, :)

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(l, i, :)
                    t7 = 0.d0

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta  = beta + (alpha - 1)*2
                        t5    = t2*mat%invJ(k, i, :) + t3*mat%invJ(k, j, :) + t4*mat%JJnn(k, i, :)
                        if (i.eq.j) t5 = t5 + kt1(2, :)*mat%JJjj(k, l, :)
                        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                t5, t6)
                        t7 = t7 + t6
                    end do

                    array_out(i, :) = array_out(i, :) + t7
                end do
            end do
        end do
            
    end subroutine mf_stiffness_2d

    subroutine mf_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes S.u in 3D where S is stiffness matrix
        !! IN CSR FORMAT

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
        integer :: i, j, k, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: kt1, t1, t2, t3, t4, t5, t6, t7
        dimension ::    kt1(3, nc_total), t1(nc_total), t2(nc_total), t3(nc_total), &
                        t4(nc_total), t5(nc_total), t6(nr_total), t7(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Number of rows not equal'
        array_out = 0.d0       
        do j = 1, dimen
            do l = 1, dimen
                beta = 1; beta(l) = 2
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in(j, :), t1) 

                do i = 1, 3
                    kt1(i, :) = mat%CepArgs(i, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(l, j, :)
                t4 = kt1(3, :)*mat%JJnn(l, j, :)

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(l, i, :)
                    t7 = 0.d0

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta  = beta + (alpha - 1)*2
                        t5    = t2*mat%invJ(k, i, :) + t3*mat%invJ(k, j, :) + t4*mat%JJnn(k, i, :)
                        if (i.eq.j) t5 = t5 + kt1(2, :)*mat%JJjj(k, l, :)
                        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t5, t6)
                        t7 = t7 + t6
                    end do

                    array_out(i, :) = array_out(i, :) + t7
                end do
            end do
        end do
            
    end subroutine mf_stiffness_3d

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
        
        if (nr_total.ne.nr_u*nr_v) stop 'Number of rows not equal'
        
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
        
        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Number of rows not equal'
        
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