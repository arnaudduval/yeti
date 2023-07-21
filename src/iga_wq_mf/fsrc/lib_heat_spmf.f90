subroutine eigendecomp_heat_2d(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, &
                                Mcoef_u, Mcoef_v, Kcoef_u, Kcoef_v, mean, isDiag, &
                                U_u, U_v, Deigen, Mdiag_u, Mdiag_v, Kdiag_u, Kdiag_v)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

    double precision, intent(in) :: Mcoef_u, Mcoef_v, Kcoef_u, Kcoef_v
    dimension ::    Mcoef_u(nc_u), Mcoef_v(nc_v), &
                    Kcoef_u(nc_u), Kcoef_v(nc_v)
    double precision, intent(in) :: mean
    dimension :: mean(2)
    logical, intent(in) :: isDiag

    double precision, intent(out) :: U_u, U_v, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), Deigen(nr_u*nr_v)
    double precision, intent(out) :: Mdiag_u, Mdiag_v, Kdiag_u, Kdiag_v
    dimension :: Mdiag_u(nr_u), Mdiag_v(nr_v), Kdiag_u(nr_u), Kdiag_v(nr_v)

    ! Local data
    ! -----------
    double precision :: D_u, D_v
    dimension :: D_u(nr_u), D_v(nr_v)
    double precision, dimension(:), allocatable :: I_u, I_v

    ! Eigen decomposition
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u, data_W_u, (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)

    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v, data_W_v, (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)  

    ! Find diagonal
    if (.not.isDiag) return
    allocate(I_u(nr_u), I_v(nr_v))
    I_u = 1.d0; I_v = 1.d0
    call find_parametric_diag_2d(nr_u, nr_v, I_u, I_v, D_u, D_v, mean, Deigen)
    deallocate(I_u, I_v)

end subroutine eigendecomp_heat_2d

subroutine eigendecomp_heat_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, mean, isDiag, &
                                U_u, U_v, U_w, Deigen, Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

    double precision, intent(in) :: Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w
    dimension ::    Mcoef_u(nc_u), Mcoef_v(nc_v), Mcoef_w(nc_w), &
                    Kcoef_u(nc_u), Kcoef_v(nc_v), Kcoef_w(nc_w)
    double precision, intent(in) :: mean
    dimension :: mean(3)
    logical, intent(in) :: isDiag

    double precision, intent(out) :: U_u, U_v, U_w, Deigen
    dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), Deigen(nr_u*nr_v*nr_w)
    double precision, intent(out) :: Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    dimension :: Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w)

    ! Local data
    ! -----------
    double precision :: D_u, D_v, D_w
    dimension :: D_u(nr_u), D_v(nr_v), D_w(nr_w)
    double precision, dimension(:), allocatable :: I_u, I_v, I_w

    ! Eigen decomposition
    call eigen_decomposition(nr_u, nc_u, Mcoef_u, Kcoef_u, nnz_u, indi_u, indj_u, &
                            data_B_u, data_W_u, (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)

    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v, data_W_v, (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)  

    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w, data_W_w, (/0, 0/), D_w, U_w, Kdiag_w, Mdiag_w)   

    ! Find diagonal
    if (.not.isDiag) return
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, mean, Deigen)
    deallocate(I_u, I_v, I_w)

end subroutine eigendecomp_heat_3d

module matrixfreeheat

    implicit none
    type thermomat
        ! Inputs 
        integer :: dimen = 3
        double precision :: scalars(2) = (/1.d0, 1.d0/)

        double precision, dimension(:), pointer :: Cprop=>null(), Ccoefs=>null(), detJ=>null(), detG=>null()
        double precision, dimension(:, :, :), pointer :: Kprop=>null(), Kcoefs=>null(), invJ=>null()
        double precision, dimension(:), allocatable :: mean

        ! Local
        integer :: ncols_sp, ncols_tm

    end type thermomat

contains

    subroutine setup_timediscret(mat, nnz, detG)
        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) :: detG
        dimension :: detG(nnz)

        mat%detG => detG
        mat%ncols_tm = nnz

    end subroutine setup_timediscret

    subroutine setup_geometry(mat, nnz, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(mat%dimen, mat%dimen, nnz), detJ(nnz)

        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nnz

    end subroutine setup_geometry

    subroutine setup_conductivitycoefs(mat, nnz, coefs)

        implicit none
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  coefs
        dimension :: coefs(mat%dimen, mat%dimen, nnz)
        mat%Kcoefs => coefs
        mat%ncols_sp = nnz

    end subroutine setup_conductivitycoefs

    subroutine setup_capacitycoefs(mat, nnz, coefs)

        implicit none
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  coefs
        dimension :: coefs(nnz)
        mat%Ccoefs => coefs
        mat%ncols_sp = nnz

    end subroutine setup_capacitycoefs

    subroutine update_conductivitycoefs(mat, info)
        !! Computes conductivity coefficients coef = J^-1 lambda detJ J^-T
        
        use omp_lib
        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
        integer, intent(out) :: info

        ! Local data
        ! ----------
        integer :: i, nnz, nb_tasks
        double precision :: invJt, Kt
        dimension :: invJt(mat%dimen, mat%dimen), Kt(mat%dimen, mat%dimen)  
        
        info = 1
        nnz  = mat%ncols_sp
    
        if (size(mat%Kprop, dim=3).eq.1) then 
    
            !$OMP PARALLEL PRIVATE(invJt, Kt)
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute K = invJ * prop * detJ * invJ'
                invJt = mat%invJ(:, :, i)
                Kt = mat%detJ(i) * matmul(invJt, mat%Kprop(:, :, 1)) 
                mat%Kcoefs(:, :, i) = matmul(Kt, transpose(invJt))
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
    
        else if (size(mat%Kprop, dim=3).eq.nnz) then 
    
            !$OMP PARALLEL PRIVATE(invJt, Kt)
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute K = invJ * prop * detJ * invJ'
                invJt = mat%invJ(:, :, i)
                Kt = mat%detJ(i) * matmul(invJt, mat%Kprop(:, :, i)) 
                mat%Kcoefs(:, :, i) = matmul(Kt, transpose(invJt))
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else 
            info = 0
            print*, "Error computing thermal coefficient (Conductivity)"
        end if
    
    end subroutine update_conductivitycoefs

    subroutine update_capacitycoefs(mat, info)
        !! Computes capacity coefficient coef = sigma * detJ
        
        use omp_lib
        implicit none 
        ! Input / output data
        ! -------------------  
        type(thermomat), pointer :: mat    
        integer, intent(out) :: info
    
        ! Local data
        ! ----------
        integer :: i, nnz, nb_tasks
    
        info = 1
        nnz  = mat%ncols_sp
    
        if (size(mat%Cprop).eq.1) then 
    
            !$OMP PARALLEL
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute C = detJ  * prop
                mat%Ccoefs(i) = mat%detJ(i) * mat%Cprop(1)
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
    
        else if (size(mat%Cprop).eq.nnz) then
    
            !$OMP PARALLEL
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute C = detJ  * prop
                mat%Ccoefs(i) = mat%detJ(i) * mat%Cprop(i)
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
    
        else
            info = 0
            print*, "Error computing thermal coefficient (Capacity) "
        end if
    
    end subroutine update_capacitycoefs

    subroutine compute_mean_3d(mat, nc_u, nc_v, nc_w)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: samplesize = 3**3
        type(thermomat), pointer :: mat
        integer, intent(in) :: nc_u, nc_v, nc_w

        ! Local data
        ! ----------
        integer :: i, j, k, l, genPos, pos
        integer :: ind_u, ind_v, ind_w, sample
        dimension :: ind_u(3), ind_v(3), ind_w(3), sample(samplesize)
        
        if (nc_u*nc_v*nc_w.ne.mat%ncols_sp) stop 'Wrong dimensions'
        pos = int((nc_u+1)/2); ind_u = (/1, pos, nc_u/)
        pos = int((nc_v+1)/2); ind_v = (/1, pos, nc_v/)
        pos = int((nc_w+1)/2); ind_w = (/1, pos, nc_w/)
    
        ! Select a set of coefficients
        l = 1
        do k = 1, 3
            do j = 1, 3
                do i = 1, 3
                    genPos = ind_u(i) + (ind_v(j) - 1)*nc_u + (ind_w(k) - 1)*nc_u*nc_v
                    sample(l) = genPos
                    l = l + 1
                end do
            end do
        end do

        if (.not.allocated(mat%mean)) allocate(mat%mean(4))
        mat%mean = 0.d0
        if (associated(mat%Kcoefs)) then
            call trapezoidal_rule_3d(3, 3, 3, mat%Kcoefs(1, 1, sample), mat%mean(1))
            call trapezoidal_rule_3d(3, 3, 3, mat%Kcoefs(2, 2, sample), mat%mean(2))
            call trapezoidal_rule_3d(3, 3, 3, mat%Kcoefs(3, 3, sample), mat%mean(3))
        end if
        if (associated(mat%Ccoefs)) then
            call trapezoidal_rule_3d(3, 3, 3, mat%Ccoefs(sample), mat%mean(4))
        end if   

    end subroutine compute_mean_3d

    !! Matrix free 3d functions

    subroutine mf_wq_capacity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes C.u where C is the capacity matrix in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
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
        double precision :: array_temp
        dimension :: array_temp(nc_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'

        call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                            array_in, array_temp)

        array_temp = array_temp * mat%Ccoefs

        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            array_temp, array_out)
        
    end subroutine mf_wq_capacity_3d

    subroutine mf_wq_conductivity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes K.u where K is conductivity matrix in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(thermomat), pointer :: mat
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
        double precision :: array_temp_0, array_temp_1, array_temp_2
        dimension :: array_temp_0(nc_total), array_temp_1(nc_total), array_temp_2(nc_total)
        integer :: i, j, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                                array_in, array_temp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                array_temp_1 = array_temp_0 * mat%Kcoefs(i, j, :)
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), & 
                                    array_temp_1, array_temp_2)
                array_out = array_out + array_temp_2
            end do
        end do
        
    end subroutine mf_wq_conductivity_3d

    subroutine mf_wq_condcap_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes (alpha*C + beta*K).u where C and K are capacity and conductivity matrices respectively in 3D case
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
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
        double precision :: array_temp
        dimension :: array_temp(nr_total)

        call mf_wq_capacity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_wq_conductivity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_temp)

        array_out = array_out + mat%scalars(2)*array_temp
        
    end subroutine mf_wq_condcap_3d

    subroutine mf_wq_spacetimeheat_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                    nnz_u, nnz_v, nnz_w, nnz_t,indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                    indi_T_t, indj_T_t, data_BT_u, data_BT_v, data_BT_w, data_BT_t, indi_u, indj_u, &
                                    indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, data_W_u, data_W_v, data_W_w, data_W_t, &
                                    array_in, array_out)
        !! Computes K.u where K is conductivity matrix in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(thermomat), pointer :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t
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
        ! -----------
        double precision :: array_temp_0, array_temp_1, array_temp_2
        dimension :: array_temp_0(nc_total), array_temp_1(nc_total), array_temp_2(nc_total)
        integer :: i, j, k, l, gen, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v*nr_w*nr_t) stop 'Size problem'

        call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                            nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 2), & 
                            array_in, array_temp_0)
        
        do j = 1, mat%ncols_tm
            do i = 1, mat%ncols_sp
                gen = i + (j-1)*mat%ncols_sp
                array_temp_0(gen) = array_temp_0(gen) * mat%Ccoefs(i)
            end do
        end do

        call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                            nnz_t, indi_t, indj_t, data_W_t(:, 2), &
                            array_temp_0, array_out)

        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto4d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                                nnz_t, indi_T_t, indj_T_t, data_BT_t(:, 1), & 
                                array_in, array_temp_0)

            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2

                do l = 1, mat%ncols_tm
                    do k = 1, mat%ncols_sp
                        gen = l + (k-1)*mat%ncols_sp
                        array_temp_1(gen) = array_temp_0(gen) * mat%Kcoefs(i, j, k) * mat%detG(l)
                    end do
                end do

                call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                    nnz_t, indi_t, indj_t, data_W_t(:, 1), & 
                                    array_temp_1, array_temp_2)

                array_out = array_out + array_temp_2
            end do
        end do
        
    end subroutine mf_wq_spacetimeheat_3d

    !! Matrix free 2d functions 
    subroutine mf_wq_capacity_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes C.u where C is the capacity matrix in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
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
        double precision :: array_temp
        dimension :: array_temp(nc_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'

        call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                            nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                            nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                            array_in, array_temp)

        array_temp = array_temp * mat%Ccoefs

        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                            nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                            nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                            array_temp, array_out)
        
    end subroutine mf_wq_capacity_2d

    subroutine mf_wq_conductivity_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes K.u where K is conductivity matrix in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(thermomat), pointer :: mat
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
        double precision :: array_temp_0, array_temp_1, array_temp_2
        dimension :: array_temp_0(nc_total), array_temp_1(nc_total), array_temp_2(nc_total)
        integer :: i, j, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'

        array_out = 0.d0
        do j = 1, dimen
            beta = 1; beta(j) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                                array_in, array_temp_0)
            do i = 1, dimen
                alpha = 1; alpha(i) = 2
                zeta = beta + (alpha - 1)*2
                array_temp_1 = array_temp_0 * mat%Kcoefs(i, j, :)
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                    nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                    nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                    array_temp_1, array_temp_2)
                array_out = array_out + array_temp_2
            end do
        end do
        
    end subroutine mf_wq_conductivity_2d

end module matrixfreeheat

subroutine compute_transientheat_cond(nnz, Kcoefs, Ccoefs, Kmean, Cmean, kappa)

    implicit none
    ! Input /  output data
    ! --------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nnz
    double precision, intent(in) :: Kcoefs, Ccoefs
    dimension :: Kcoefs(dimen, dimen, nnz), Ccoefs(nnz)
    double precision, intent(in) :: Kmean, Cmean
    dimension :: Kmean(dimen)

    double precision, intent(out) :: kappa
    
    ! Local data
    ! ----------
    integer :: i, j, k
    double precision :: eye, KK, Kmean2, rho, x
    dimension :: eye(dimen, dimen), KK(dimen, dimen), Kmean2(dimen), rho(dimen), x(dimen, dimen)
    double precision :: supKold, infKold, supCold, infCold, supKnew, infKnew, Cnew, supK, infK, supC, infC

    eye = 0.d0
    do i = 1, dimen
        eye(i, i) = 1.d0
    end do

    do i = 1, dimen
        kmean2(i) = 1.0/sqrt(Kmean(i))
    end do

    KK = Kcoefs(:, :, 1)
    do i = 1, dimen
        do j = 1, dimen
            KK(i, j) = KK(i, j) * Kmean2(i) * Kmean2(j)
        end do
    end do

    call compute_geneigs(dimen, KK, eye, rho, x)
    supKold = maxval(rho); infKold = minval(rho)
    supCold = Ccoefs(1)/Cmean; infCold = Ccoefs(1)/Cmean

    do k = 2, nnz
        
        KK = Kcoefs(:, :, k)
        do i = 1, dimen
            do j = 1, dimen
                KK(i, j) = KK(i, j) * Kmean2(i) * Kmean2(j)
            end do
        end do

        call compute_geneigs(dimen, KK, eye, rho, x)
        supKnew = maxval(rho); infKnew = minval(rho)
        supK = max(supKold, supKnew); infK = min(infKold, infKnew)
        Cnew = Ccoefs(k)/Cmean
        supC = max(supCold, Cnew); infC = min(infCold, Cnew)

        supKold = supK; infKold = infK
        supCold = supC; infCold = infC

    end do
    
    kappa = max(supK, supC)/min(infK, infC)

end subroutine compute_transientheat_cond

subroutine compute_steadyheat_cond(nnz, Kcoefs, Kmean, kappa)
    
    implicit none
    ! Input /  output data
    ! --------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nnz
    double precision, intent(in) :: Kcoefs
    dimension :: Kcoefs(dimen, dimen, nnz)
    double precision, intent(in) :: Kmean
    dimension :: Kmean(dimen)

    double precision, intent(out) :: kappa
    
    ! Local data
    ! ----------
    integer :: i, j, k
    double precision :: eye, KK, Kmean2, rho, x
    dimension :: eye(dimen, dimen), KK(dimen, dimen), Kmean2(dimen), rho(dimen), x(dimen, dimen)
    double precision :: supKold, infKold, supKnew, infKnew, supK, infK

    eye = 0.d0
    do i = 1, dimen
        eye(i, i) = 1.d0
    end do

    do i = 1, dimen
        kmean2(i) = 1.0/sqrt(Kmean(i))
    end do

    KK = Kcoefs(:, :, 1)
    do i = 1, dimen
        do j = 1, dimen
            KK(i, j) = KK(i, j) * Kmean2(i) * Kmean2(j)
        end do
    end do

    call compute_geneigs(dimen, KK, eye, rho, x)
    supKold = maxval(rho); infKold = minval(rho)

    do k = 2, nnz
        
        KK = Kcoefs(:, :, k)
        do i = 1, dimen
            do j = 1, dimen
                KK(i, j) = KK(i, j) * Kmean2(i) * Kmean2(j)
            end do
        end do

        call compute_geneigs(dimen, KK, eye, rho, x)
        supKnew = maxval(rho); infKnew = minval(rho)
        supK = max(supKold, supKnew); infK = min(infKold, infKnew)

        supKold = supK; infKold = infK

    end do
    
    kappa = supK/infK

end subroutine compute_steadyheat_cond 