subroutine eigendecomp_plasticity_2d(nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, table, ddl, U_u, U_v, D_u, D_v)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 2
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, ddl
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    integer, intent(in) :: table
    dimension :: table(dimen, 2, dimen)

    double precision, intent(out) :: U_u, U_v
    dimension :: U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl)
    double precision, intent(out) :: D_u, D_v
    dimension :: D_u(nr_u, ddl), D_v(nr_v, ddl)

    ! Local data
    ! -----------
    integer :: i, j, c
    double precision, allocatable, dimension(:) :: Mdiag_u, Mdiag_v, Kdiag_u, Kdiag_v, &
                                                    Mcoef_u, Mcoef_v, Kcoef_u, Kcoef_v
    if (.not.any((/dimen, dimen**2/).eq.ddl)) stop 'Not possible eigen decomposition'
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v)) 
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0

    if (ddl.eq.dimen) then
        do i = 1, dimen
            call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                    data_B_u, data_W_u, table(1, :, i), table(1, :, i), &
                                    D_u(:, i), U_u(:, :, i), Kdiag_u, Mdiag_u)

            call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                    data_B_v, data_W_v, table(2, :, i), table(2, :, i), &
                                    D_v(:, i), U_v(:, :, i), Kdiag_v, Mdiag_v) 
        end do
    else
        do i = 1, dimen
            do j = 1, dimen
                c = j + (i-1)*dimen
                call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                        data_B_u, data_W_u, table(1, :, i), table(1, :, j), &
                                        D_u(:, c), U_u(:, :, c), Kdiag_u, Mdiag_u)

                call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                        data_B_v, data_W_v, table(2, :, i), table(2, :, j), &
                                        D_v(:, c), U_v(:, :, c), Kdiag_v, Mdiag_v)

            end do
        end do
    end if

    deallocate(Mdiag_u, Mdiag_v, Kdiag_u, Kdiag_v)
    deallocate(Mcoef_u, Mcoef_v, Kcoef_u, Kcoef_v)
    
end subroutine eigendecomp_plasticity_2d

subroutine eigendecomp_plasticity_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, table, &
                                ddl, U_u, U_v, U_w, D_u, D_v, D_w)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, ddl
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    integer, intent(in) :: table
    dimension :: table(dimen, 2, dimen)

    double precision, intent(out) :: U_u, U_v, U_w
    dimension :: U_u(nr_u, nr_u, ddl), U_v(nr_v, nr_v, ddl), U_w(nr_w, nr_w, ddl)
    double precision, intent(out) :: D_u, D_v, D_w
    dimension :: D_u(nr_u, ddl), D_v(nr_v, ddl), D_w(nr_w, ddl)

    ! Local data
    ! -----------
    integer :: i, j, c
    double precision, allocatable, dimension(:) :: Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w, &
                                                    Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w

    if (.not.any((/dimen, dimen**2/).eq.ddl)) stop 'Not possible eigen decomposition'
    allocate(Kdiag_u(nr_u), Mdiag_u(nr_u), Kdiag_v(nr_v), Mdiag_v(nr_v), Kdiag_w(nr_w), Mdiag_w(nr_w))
    allocate(Mcoef_u(nc_u), Kcoef_u(nc_u), Mcoef_v(nc_v), Kcoef_v(nc_v), Mcoef_w(nc_w), Kcoef_w(nc_w)) 
    Mcoef_u = 1.d0; Kcoef_u = 1.d0
    Mcoef_v = 1.d0; Kcoef_v = 1.d0
    Mcoef_w = 1.d0; Kcoef_w = 1.d0

    if (ddl.eq.dimen) then
        do i = 1, dimen
            call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                    data_B_u, data_W_u, table(1, :, i), table(1, :, i), &
                                    D_u(:, i), U_u(:, :, i), Kdiag_u, Mdiag_u)

            call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                    data_B_v, data_W_v, table(2, :, i), table(2, :, i), &
                                    D_v(:, i), U_v(:, :, i), Kdiag_v, Mdiag_v)

            call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                    data_B_w, data_W_w, table(3, :, i), table(3, :, i), &
                                    D_w(:, i), U_w(:, :, i), Kdiag_w, Mdiag_w) 
        end do
    else
        do i = 1, dimen
            do j = 1, dimen
                c = j + (i-1)*dimen
                call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                                        data_B_u, data_W_u, table(1, :, i), table(1, :, j), &
                                        D_u(:, c), U_u(:, :, c), Kdiag_u, Mdiag_u)

                call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                                        data_B_v, data_W_v, table(2, :, i), table(2, :, j), &
                                        D_v(:, c), U_v(:, :, c), Kdiag_v, Mdiag_v)

                call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                                        data_B_w, data_W_w, table(3, :, i), table(3, :, j), &
                                        D_w(:, c), U_w(:, :, c), Kdiag_w, Mdiag_w) 
            end do
        end do
    end if

    deallocate(Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)
    deallocate(Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w)

end subroutine eigendecomp_plasticity_3d

module matrixfreeplasticity

    implicit none
    type :: mecamat
    
        ! Inputs 
        integer :: dimen, nvoigt
        double precision :: elasticmodulus, poissonratio, elasticlimit
        double precision, dimension(:), pointer :: detJ=>null()
        double precision, dimension(:, :), pointer :: CepArgs=>null(), NN=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn
        double precision, dimension(:, :, :), allocatable :: mean
        
        ! Outputs
        double precision :: lambda, mu, bulk

        ! Local
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
        type(mecamat), pointer :: mat

        allocate(mat)
        mat%dimen  = dimen
        mat%nvoigt = dimen*(dimen+1)/2

        ! Compute mechanical properties
        mat%elasticmodulus = elasticmodulus
        mat%poissonratio   = poissonratio
        mat%elasticlimit   = elasticlimit
        mat%lambda  = poissonratio*elasticmodulus/((1+poissonratio)*(1-2*poissonratio))
        mat%mu      = elasticmodulus/(2*(1+poissonratio))
        mat%bulk    = mat%lambda + 2.d0/3.d0*mat%mu
        allocate(mat%mean(mat%dimen, mat%dimen, mat%dimen))
        mat%mean    = 1.d0

    end subroutine initialize_mecamat

    subroutine setup_geometry(mat, nc, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        integer, intent(in) :: nc
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(mat%dimen, mat%dimen, nc), detJ(nc)

        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nc

    end subroutine setup_geometry

    subroutine setup_jacobiennormal(mat, CepArgs)
        !! Points to data of the mechanical behavior 

        implicit none
        ! Input / output data
        ! -------------------
        double precision, parameter :: threshold = 1.d-12
        type(mecamat), pointer :: mat
        double precision, target, intent(in) :: CepArgs
        dimension :: CepArgs(mat%nvoigt+3, mat%ncols_sp)

        ! Local data
        ! ----------
        double precision :: TNN(mat%dimen, mat%dimen)
        integer :: i

        mat%CepArgs => CepArgs(:3, :)
        mat%NN      => CepArgs(4:, :)
        if (.not.allocated(mat%JJnn)) allocate(mat%JJnn(mat%dimen, mat%dimen, mat%ncols_sp))
        mat%JJnn = 0.d0
        if (any(abs(mat%NN).gt.threshold)) then
            do i = 1,  mat%ncols_sp
                call array2symtensor(mat%dimen, mat%nvoigt, mat%NN(:, i), TNN)
                mat%JJnn(:, :, i) = matmul(mat%invJ(:, :, i), TNN)
            end do
        end if
    end subroutine setup_jacobiennormal

    subroutine setup_jacobienjacobien(mat)
        
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat

        ! Local data
        ! ----------
        integer :: i
        
        if (.not.allocated(mat%JJjj)) allocate(mat%JJjj(mat%dimen, mat%dimen, mat%ncols_sp))
        mat%JJjj = 0.d0
        do i = 1, mat%ncols_sp
            mat%JJjj(:, :, i) = matmul(mat%invJ(:, :, i), transpose(mat%invJ(:, :, i)))
        end do

    end subroutine setup_jacobienjacobien

    subroutine compute_mean_ijblock(i, j, dimen, nvoigt, mat, samplesize, sample, mean)
        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
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
            NN = mat%NN(genpos, :)
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
            call trapezoidal_rule_2d(3, 3, coefs(n, n, :), mean(n))
        end do   

    end subroutine compute_mean_ijblock

    subroutine compute_mean_diagblocks(mat, dimen, nvoigt, nclist)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
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
            mat%mean(i, i, :) = mean
        end do

    end subroutine compute_mean_diagblocks

    subroutine compute_mean_allblocks(mat, dimen, nvoigt, nclist)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
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
            stop 'Not possible compute mean all blocks'
        end if

        do i = 1, dimen
            do j = 1, dimen
                call compute_mean_ijblock(i, j, dimen, nvoigt, mat, size(sample), sample, mean)
                mat%mean(i, j, :) = mean
            end do
        end do

    end subroutine compute_mean_allblocks

    subroutine mf_wq_stiffness_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 3D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat), pointer :: mat
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
            
    end subroutine mf_wq_stiffness_2d

    subroutine mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes S.u in 3D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat), pointer :: mat
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
            
    end subroutine mf_wq_stiffness_3d

    subroutine wq_intforce_2d(mat, stress, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat), pointer :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

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

    end subroutine wq_intforce_2d

    subroutine wq_intforce_3d(mat, stress, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat), pointer :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

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

    end subroutine wq_intforce_3d

end module matrixfreeplasticity