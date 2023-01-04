! ==========================
! module :: Heat transfer 
! author :: Joaquin Cornejo
! ==========================

subroutine eigendecomp_heat_3d(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                Mcoef_u, Mcoef_v, Mcoef_w, Kcoef_u, Kcoef_v, Kcoef_w, mean, doDeigen, &
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
    dimension ::    Mcoef_u(nr_u), Mcoef_v(nr_v), Mcoef_w(nr_w), &
                    Kcoef_u(nr_u), Kcoef_v(nr_v), Kcoef_w(nr_w)
    double precision, intent(in) :: mean
    dimension :: mean(3)
    logical, intent(in) :: doDeigen

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
                            data_B_u(:, 1), data_W_u(:, 1), data_B_u(:, 2), &
                            data_W_u(:, 4), (/0, 0/), D_u, U_u, Kdiag_u, Mdiag_u)

    call eigen_decomposition(nr_v, nc_v, Mcoef_v, Kcoef_v, nnz_v, indi_v, indj_v, &
                            data_B_v(:, 1), data_W_v(:, 1), data_B_v(:, 2), &
                            data_W_v(:, 4), (/0, 0/), D_v, U_v, Kdiag_v, Mdiag_v)  

    call eigen_decomposition(nr_w, nc_w, Mcoef_w, Kcoef_w, nnz_w, indi_w, indj_w, &
                            data_B_w(:, 1), data_W_w(:, 1), data_B_w(:, 2), &
                            data_W_w(:, 4), (/0, 0/), D_w, U_w, Kdiag_w, Mdiag_w)   

    ! Find diagonal
    if (.not.doDeigen) return
    allocate(I_u(nr_u), I_v(nr_v), I_w(nr_w))
    I_u = 1.d0; I_v = 1.d0; I_w = 1.d0
    call find_parametric_diag_3d(nr_u, nr_v, nr_w, I_u, I_v, I_w, D_u, D_v, D_w, mean, Deigen)
    deallocate(I_u, I_v, I_w)

end subroutine eigendecomp_heat_3d

subroutine interpolate_isotropic_prop(nr, table_prop, nnz, temperature, prop)

    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: span_tol = 1d-8
    integer, intent(in) :: nr, nnz
    double precision, intent(in) :: table_prop, temperature
    dimension :: table_prop(nr, 2), temperature(nnz)

    double precision, intent(out) :: prop
    dimension :: prop(nnz)

    call linear_interpolation(nr, table_prop, nnz, temperature, prop, span_tol)

end subroutine interpolate_isotropic_prop

subroutine interpolate_temperature(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, & 
                            data_B_u, data_B_v, data_B_w, T_ctrlpts, T_interp)
    !! Computes strain in 3D (from parametric space to physical space)
    !! IN CSR FORMAT

    implicit none 
    ! Input / output data
    ! -------------------  
    integer, intent(in) ::  nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
    integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w
    dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
    double precision, intent(in) :: T_ctrlpts
    dimension :: T_ctrlpts(nr_u*nr_v*nr_w)

    double precision, intent(out) :: T_interp
    dimension :: T_interp(nc_u*nc_v*nc_w)

    ! Local data
    !-----------
    integer :: indi_T_u, indi_T_v, indi_T_w
    dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
    integer :: indj_T_u, indj_T_v, indj_T_w
    dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
    double precision :: data_BT_u, data_BT_v, data_BT_w
    dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

    call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
    call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
    call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

    call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), &
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), &
                        T_ctrlpts, T_interp)

end subroutine interpolate_temperature

module heat_spmf

    implicit none

    type thermomat
        ! Inputs 
        integer :: dimen = 3
        double precision :: scalars(2) = (/1.d0, 1.d0/)

        double precision, dimension(:), pointer :: Cprop, Ccoefs, detJJ
        double precision, dimension(:,:,:), pointer :: Kprop, Kcoefs, invJJ
        double precision, dimension(:), allocatable :: mean

        ! Local
        integer :: nc_total

    end type thermomat

contains

    subroutine setup_geo(mat, nc, invJJ, detJJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
        integer, intent(in) :: nc
        double precision, target, intent(in) :: invJJ, detJJ
        dimension :: invJJ(mat%dimen, mat%dimen, nc), detJJ(nc)

        mat%invJJ => invJJ
        mat%detJJ => detJJ
        mat%nc_total = nc

    end subroutine setup_geo

    subroutine setupKcoefs(mat, nnz, coefs)
        implicit none
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  coefs
        dimension :: coefs(mat%dimen, mat%dimen, nnz)
        mat%Kcoefs => coefs
        mat%nc_total = nnz
    end subroutine setupKcoefs

    subroutine setupCcoefs(mat, nnz, coefs)
        implicit none
        type(thermomat), pointer :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  coefs
        dimension :: coefs(nnz)
        mat%Ccoefs => coefs
        mat%nc_total = nnz
    end subroutine setupCcoefs

    subroutine update_conductivity_coefs(mat, info)
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
        nnz  = mat%nc_total
    
        if (size(mat%Kprop, dim=3).eq.1) then 
    
            !$OMP PARALLEL PRIVATE(invJt, Kt)
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute K = invJ * prop * detJ * invJ'
                invJt = mat%invJJ(:, :, i)
                Kt = mat%detJJ(i) * matmul(invJt, mat%Kprop(:, :, 1)) 
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
                invJt = mat%invJJ(:, :, i)
                Kt = mat%detJJ(i) * matmul(invJt, mat%Kprop(:, :, i)) 
                mat%Kcoefs(:, :, i) = matmul(Kt, transpose(invJt))
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else 
            info = 0
            print*, "Error computing thermal coefficient (Conductivity)"
        end if
    
    end subroutine update_conductivity_coefs

    subroutine update_capacity_coefs(mat, info)
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
        nnz  = mat%nc_total
    
        if (size(mat%Cprop).eq.1) then 
    
            !$OMP PARALLEL
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute C = detJ  * prop
                mat%Ccoefs(i) = mat%detJJ(i) * mat%Cprop(1)
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
    
        else if (size(mat%Cprop).eq.nnz) then
    
            !$OMP PARALLEL
            nb_tasks = omp_get_num_threads()
            !$OMP DO SCHEDULE(STATIC, nnz/nb_tasks) 
            do i = 1, nnz
                ! Compute C = detJ  * prop
                mat%Ccoefs(i) = mat%detJJ(i) * mat%Cprop(i)
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
    
        else
            info = 0
            print*, "Error computing thermal coefficient (Capacity) "
        end if
    
    end subroutine update_capacity_coefs

    !! Isotropic functions

    subroutine eval_isotropic_coefs(mat, nnz, Kprop, Cprop, Kcoefs, Ccoefs)
        !! Computes the coefficients to use to calculate source vector and tangent heat matrix
        !! By the moment, it only considers isotropic materials
        !! To verify in the future

        implicit none 
        ! Input / output data
        ! -------------------
        type(thermomat), pointer :: mat
        integer :: nnz
        double precision, intent(in) :: Kprop, Cprop
        dimension :: Kprop(nnz), Cprop(nnz)

        double precision, intent(out) :: Kcoefs, Ccoefs
        dimension :: Kcoefs(mat%dimen, mat%dimen, nnz), Ccoefs(nnz)

        ! Local data
        ! ----------
        double precision :: invJJt, detJJt
        dimension :: invJJt(mat%dimen, mat%dimen)
        integer :: k

        do k = 1, nnz

            invJJt = mat%invJJ(:, :, k)
            detJJt = mat%detJJ(k)
            Kcoefs(:, :, k) = Kprop(k) * matmul(invJJt, transpose(invJJt)) * detJJt
            Ccoefs(k) = Cprop(k) * detJJt

        end do

    end subroutine eval_isotropic_coefs

    subroutine compute_mean_heat_3d(mat, nc_u, nc_v, nc_w)
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
        
        if (nc_u*nc_v*nc_w.ne.mat%nc_total) stop 'Wrong dimensions'
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
            call compute_mean_3d(3, 3, 3, mat%Kcoefs(1, 1, sample), mat%mean(1))
            call compute_mean_3d(3, 3, 3, mat%Kcoefs(2, 2, sample), mat%mean(2))
            call compute_mean_3d(3, 3, 3, mat%Kcoefs(3, 3, sample), mat%mean(3))
        end if
        if (associated(mat%Ccoefs)) then
            call compute_mean_3d(3, 3, 3, mat%Ccoefs(sample), mat%mean(4))
        end if      

    end subroutine compute_mean_heat_3d

    !! Matrix free 3d functions

    subroutine mf_wq_get_cu_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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
        
    end subroutine mf_wq_get_cu_3d

    subroutine mf_wq_get_ku_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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
        
    end subroutine mf_wq_get_ku_3d

    subroutine mf_wq_get_kcu_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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

        call mf_wq_get_cu_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_wq_get_ku_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, array_in, array_temp)

        array_out = array_out + mat%scalars(2)*array_temp
        
    end subroutine mf_wq_get_kcu_3d

end module heat_spmf
