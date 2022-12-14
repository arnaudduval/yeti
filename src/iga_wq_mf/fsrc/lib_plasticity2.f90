! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
!
! Remarks :: tensor notation is used (it is NOT Voigt notation)
! We save some memory storing only the upper triangular matrix
! For example:
! Strain e = [e11, e22, e33, e12, e13, e23]
! Stress s = [s11, s22, s33, s12, s13, s23]
! In this way, we try to avoid some misconceptions when using Voigt notation
! ==========================

subroutine tensor2array(dimen, nvoigt, matrix, array)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nvoigt
    double precision, intent(in) :: matrix
    dimension :: matrix(dimen, dimen)
    double precision, intent(out) :: array
    dimension :: array(nvoigt)

    ! Local data
    ! ----------
    integer :: i, j, k

    array = 0.d0; k = 0
    do i = 1, dimen
        k = k + 1
        array(k) = matrix(i, i)
    end do
    
    do i = 1, dimen-1
        do j = i+1, dimen
            k = k + 1
            array(k) = matrix(i, j)
        end do
    end do

end subroutine tensor2array

subroutine array2tensor(dimen, nvoigt, array, matrix)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, nvoigt
    double precision, intent(in) :: array
    dimension :: array(nvoigt)
    double precision, intent(out) :: matrix
    dimension :: matrix(dimen, dimen)
    
    ! Local data
    ! ----------
    integer :: i, j, k

    matrix = 0.d0; k = 0
    do i = 1, dimen
        k = k + 1
        matrix(i, i) = array(k)
    end do
    
    do i = 1, dimen-1
        do j = i+1, dimen
            k = k + 1
            matrix(i, j) = array(k)
            matrix(j, i) = array(k)
        end do
    end do

end subroutine array2tensor

! ==========================

subroutine compute_trace(dimen, tensor, trace)
    
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen
    double precision, intent(in) :: tensor
    dimension :: tensor(dimen, dimen)

    double precision, intent(out) :: trace

    ! Local data
    ! ----------
    integer :: i

    trace = 0.d0
    do i = 1, dimen
        trace = trace + tensor(i, i)
    end do   

end subroutine compute_trace

subroutine compute_stress_norm(dimen, tensor, norm)
    !! Returns Frobenius norm of a second-order stress-like tensor 
    !! The definition is norm = sqrt(s:s). Since s is a tensor written in Voigt notation
    !! some scale factor must be considered.

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen
    double precision, intent(in) :: tensor
    dimension :: tensor(dimen, dimen)

    double precision, intent(out) :: norm

    ! Local data
    ! ----------
    integer :: i, j

    norm = 0.0d0

    do i = 1, dimen
        do j = 1, dimen
            norm = norm + tensor(i, j)*tensor(i, j)
        end do
    end do

    norm = sqrt(norm)
    
end subroutine compute_stress_norm

subroutine compute_stress_vonmises(dimen, tensor, VM)
    !! Computes equivalent stress with Von Mises formula of a second-order stress-like tensor.
    !! Given tensor is written in Voigt notation

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen
    double precision, intent(in) :: tensor
    dimension :: tensor(dimen, dimen)

    double precision, intent(out) :: VM

    ! Local data
    ! ----------
    integer :: i
    double precision :: dev, trace
    dimension :: dev(dimen, dimen)

    call compute_trace(dimen, tensor, trace)
    dev = tensor
    do i = 1, dimen
        dev(i, i) = dev(i, i) - 1.d0/3.d0*trace
    end do
    
    call compute_stress_norm(dimen, dev, VM)
    VM = sqrt(3.d0/2.d0)*VM

end subroutine compute_stress_vonmises

module linearelastoplasticity

    implicit none
    type :: mecamat
    
        ! Inputs 
        integer :: dimen=3, nvoigt=6
        double precision :: E, H, beta, poisson, sigma_Y
        double precision, dimension(:), pointer :: detJJ
        double precision, dimension(:, :), pointer :: kwargs
        double precision, dimension(:, :, :), pointer :: invJJ
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn
        
        ! Outputs
        double precision :: lambda, mu, bulk
    
    end type mecamat

contains

    ! ======================================
    ! COMBINED ISOTROPIC/KINEMATIC HARDENING
    ! ======================================

    subroutine setupGeo(mat, nc, invJJ, detJJ)
        implicit none
        type(mecamat), pointer :: mat
        integer, intent(in) :: nc
        double precision, target, intent(in) :: invJJ, detJJ
        dimension :: invJJ(mat%dimen, mat%dimen, nc), detJJ(nc)
        integer :: i

        mat%invJJ => invJJ
        mat%detJJ => detJJ

        allocate(mat%JJjj(mat%dimen, mat%dimen, nc))
        do i = 1, nc
            mat%JJjj(:,:,i) = matmul(mat%invJJ(:,:,i), transpose(mat%invJJ(:,:,i)))
        end do
    end subroutine setupGeo

    subroutine setupScoefs(mat, nr, nc, kwargs)
        implicit none
        type(mecamat), pointer :: mat
        integer, intent(in) :: nr, nc
        double precision, target, intent(in) :: kwargs
        dimension :: kwargs(nr, nc)

        double precision :: NN(mat%dimen, mat%dimen)
        integer :: i

        mat%kwargs => kwargs(:3, :)
        allocate(mat%JJnn(mat%dimen, mat%dimen, nc))
        do i = 1, nc
            call array2tensor(mat%dimen, mat%nvoigt, kwargs(4:, i), NN)
            mat%JJnn(:,:,i) = matmul(mat%invJJ(:,:,i), NN)
        end do
    end subroutine setupScoefs

    subroutine initialize_mecamat(mat, dimen, E, H, beta, nu, sigma_Y)
        !! Creates a material using combined isotropic/kinematic hardening theory 

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: dimen
        double precision, intent(in) :: E, H, beta, nu, sigma_Y
        type(mecamat), pointer :: mat

        allocate(mat)
        mat%dimen  = dimen
        mat%nvoigt = dimen*(dimen+1)/2

        ! Compute mechanical properties
        mat%E    = E
        mat%H    = H
        mat%beta = beta
        mat%poisson = nu
        mat%sigma_Y = sigma_Y
        mat%lambda  = nu*E/((1+nu)*(1-2*nu))
        mat%mu      = E/(2*(1+nu))
        mat%bulk    = mat%lambda + 2.d0/3.d0*mat%mu

    end subroutine initialize_mecamat

    subroutine cpp_combined_hardening(mat, strain, pls, a, b, pls_new, a_new, b_new, stress, kwargs)
        !! Return closest point proyection (cpp) in combined hardening criteria

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        double precision, intent(in) :: strain, pls, a, b
        dimension :: strain(mat%nvoigt), pls(mat%nvoigt), b(mat%nvoigt)

        double precision, intent(out) ::pls_new, a_new, b_new, stress, kwargs
        dimension :: pls_new(mat%nvoigt), b_new(mat%nvoigt), stress(mat%nvoigt), kwargs(mat%nvoigt+3)
        
        ! Local data
        ! ----------
        double precision, dimension(mat%dimen, mat%dimen) :: Tstrain, Tpls, Tb, devStrain, s_trial, eta_trial, sigma, NN
        double precision :: traceStrain, norm_trial, f_trial, dgamma, c1, c2
        integer :: i

        call array2tensor(mat%dimen, mat%nvoigt, strain, Tstrain)
        call array2tensor(mat%dimen, mat%nvoigt, pls, Tpls)
        call array2tensor(mat%dimen, mat%nvoigt, b, Tb)
        call compute_trace(mat%dimen, Tstrain, traceStrain)
        devStrain = Tstrain
        do i = 1, mat%dimen
            devStrain(i, i) = devStrain(i, i) - 1.d0/3.d0*traceStrain
        end do

        ! Compute trial stress
        s_trial = 2*mat%mu*(devStrain - Tpls)

        ! Compute shifted stress
        eta_trial = s_trial - Tb

        ! Check yield condition
        call compute_stress_norm(mat%dimen, eta_trial, norm_trial)
        f_trial = norm_trial - sqrt(2.d0/3.d0)*(mat%sigma_Y + mat%beta*mat%H*a)  
        sigma = s_trial + mat%bulk*Tstrain
        kwargs = 0.d0; kwargs(1) = mat%lambda; kwargs(2) = mat%mu; 
        
        if (f_trial.le.0.d0) then ! Elastic point
            pls_new = pls; a_new = a; b_new = b  
            return
        end if

        ! Compute df/dsigma
        NN = eta_trial/norm_trial
        call tensor2array(mat%dimen, mat%nvoigt, NN, kwargs(4:))

        ! Compute plastic-strain increment        
        dgamma = 3.d0*f_trial/(6.d0*mat%mu+2.d0*mat%H)

        ! Update stress
        sigma = sigma - 2*mat%mu*dgamma*NN
        call tensor2array(mat%dimen, mat%nvoigt, sigma, stress)

        ! Update plastic strain
        Tpls = Tpls + dgamma*NN
        call tensor2array(mat%dimen, mat%nvoigt, Tpls, pls_new)
        
        ! Update internal hardening variable
        a_new  = a + sqrt(2.d0/3.d0)*dgamma
        
        ! Update back stress
        Tb  = Tb + 2.d0/3.d0*(1 - mat%beta)*mat%H*dgamma*NN
        call tensor2array(mat%dimen, mat%nvoigt, tb, b_new)

        ! Compute some coefficients
        c1 = 2.d0*mat%mu*dgamma/norm_trial
        c2 = 1.d0/(1.d0 + mat%H/(3.d0*mat%mu)) - c1
        kwargs(1) = mat%lambda + 2.d0*mat%mu*c1/3.d0
        kwargs(2) = mat%mu*(1.d0 - c1)
        kwargs(3) = -2.d0*mat%mu*c2

    end subroutine cpp_combined_hardening

    subroutine mf_wq_get_su_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        !! Computes S.u in 3D where S is stiffness matrix
        !! IN CSR FORMAT

        use heat_spmf

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
        integer :: i, j, k, l, r, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: kt1, t1, t2, t3, t4, t5, t6
        dimension :: kt1(dimen, nc_total), t1(nc_total), t2(nc_total), t3(nr_total), t4(nr_total), t5(nr_total), t6(nr_total)

        array_out = 0.d0
        do j = 1, dimen
            do l = 1, dimen
                beta = 1; beta(l) = 2
                call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in(j, :), t1) 

                do r = 1, dimen
                    kt1(r, :) = mat%kwargs(r, :)*t1*mat%detJJ
                end do
                    
                t2 = kt1(1, :)*mat%invJJ(l, j, :)
                t4 = kt1(3, :)*mat%JJnn(l, j, :)

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJJ(l, i, :)

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta = beta + (alpha - 1)*2

                        t5 = t2*mat%invJJ(k, i, :) + t3*mat%invJJ(k, j, :) + t4*mat%JJnn(k, i, :)
                        if (i.eq.j) t5 = t5 + kt1(2, :)*mat%JJjj(k, l, :)

                        call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t5, t6)
                        array_out(i, :) = array_out(i, :) + t6

                    end do

                end do

            end do
        end do

        ! array_out = 0.d0
        ! do j = 1, dimen
        !     do l = 1, dimen
        !         beta = 1; beta(l) = 2
        !         call sumproduct3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
        !                 nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
        !                 nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
        !                 nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
        !                 array_in(j, :), t1) 
        !         t1 = t1 * mat%detJJ ! diagonal scaling

        !         do k = 1, dimen
        !             alpha = 1; alpha(k) = 2
        !             zeta = beta + (alpha - 1)*2

        !             do i = 1, dimen
        !                 t2 = (mat%kwargs(1, :)*mat%invJJ(k, i, :)*mat%invJJ(l, j, :) &
        !                     + mat%kwargs(2, :)*mat%invJJ(l, i, :)*mat%invJJ(k, j, :) &
        !                     + mat%kwargs(3, :)*mat%JJnn(k, i, :)*mat%JJnn(l, j, :))! diagonal scaling
        !                 if (i.eq.j) t2 = t2 + mat%kwargs(2, :)*mat%JJjj(k, l, :)
        !                 t2 = t2*t1  

        !                 call sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
        !                         nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
        !                         nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
        !                         nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t2, t3)
        !                 array_out(i, :) = array_out(i, :) + t3

        !             end do

        !         end do

        !     end do
        ! end do

    end subroutine mf_wq_get_su_3d

end module linearelastoplasticity