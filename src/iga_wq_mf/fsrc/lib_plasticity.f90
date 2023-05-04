! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
!
! Remarks :: tensor notation is used (it is NOT Voigt notation)
! We save some memory storing only the upper triangular of a symmetric matrix
! For example:
! Strain e = [e11, e22, e33, e12, e13, e23]
! Stress s = [s11, s22, s33, s12, s13, s23]
! In this way, we try to avoid some misconceptions when using Voigt notation
! ==========================

subroutine symtensor2array(dimen, nvoigt, matrix, array)
    !! Returns the upper triangular part of a matrix

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

end subroutine symtensor2array

subroutine array2symtensor(dimen, nvoigt, array, matrix)
    !! Returns the matrix built from the upper triangular part

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

end subroutine array2symtensor

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

    call eval_trace(dimen, tensor, trace)
    dev = tensor
    do i = 1, dimen
        dev(i, i) = dev(i, i) - 1.d0/3.d0*trace
    end do
    
    call eval_frobenius_norm(dimen, dev, VM)
    VM = sqrt(3.d0/2.d0)*VM

end subroutine compute_stress_vonmises

subroutine fd_linearelasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    use solverheat
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, diag, array_in
    dimension ::    U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), U_w(nr_w, nr_w, dimen), &
                    diag(dimen, nr_total), array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv
    double precision :: temp
    dimension :: temp(nr_total) 
    integer :: i

    allocate(solv)
    do i = 1, dimen 
        temp = diag(i, :)
        call setup_preconditionerdiag(solv, nr_total, temp)
        call applyfastdiag(solv, nr_total, nr_u, nr_v, nr_w, U_u(:, :, i), U_v(:, :, i), U_w(:, :, i), &
                                array_in(i, :), array_out(i, :))
    end do
    
end subroutine fd_linearelasticity_3d

module linearelastoplasticity

    implicit none
    type :: mecamat
    
        ! Inputs 
        integer :: dimen=3, nvoigt=6
        double precision :: E, H, beta, poisson, sigma_Y
        double precision, dimension(:), pointer :: detJJ=>null()
        double precision, dimension(:, :), pointer :: kwargs=>null(), nn=>null()
        double precision, dimension(:, :, :), pointer :: invJJ=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn
        double precision, dimension(:, :), allocatable :: mean
        
        ! Outputs
        double precision :: lambda, mu, bulk

        ! Local
        logical :: isElastic = .false.
        integer :: nc_total
    
    end type mecamat

contains

    ! ======================================
    ! COMBINED ISOTROPIC/KINEMATIC HARDENING
    ! ======================================

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
        allocate(mat%mean(mat%dimen, mat%dimen))
        mat%mean    = 1.d0

    end subroutine initialize_mecamat

    subroutine setup_geo(mat, nc, invJJ, detJJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        integer, intent(in) :: nc
        double precision, target, intent(in) :: invJJ, detJJ
        dimension :: invJJ(mat%dimen, mat%dimen, nc), detJJ(nc)

        ! Local data
        ! ----------
        integer :: i

        mat%invJJ => invJJ
        mat%detJJ => detJJ
        mat%nc_total = nc

        allocate(mat%JJjj(mat%dimen, mat%dimen, nc))
        allocate(mat%JJnn(mat%dimen, mat%dimen, nc))
        mat%JJnn = 0.d0
        do i = 1, nc
            mat%JJjj(:,:,i) = matmul(mat%invJJ(:,:,i), transpose(mat%invJJ(:,:,i)))
        end do

    end subroutine setup_geo

    subroutine returnMappingAlgorithm(mat, strain, pls, a, b, pls_new, a_new, b_new, stress, kwargs)
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
        double precision, dimension(mat%dimen, mat%dimen) :: Tstrain, Tpls, Tb, e_trial, s_trial, eta_trial, sigma, NN
        double precision :: traceStrain, norm_trial, f_trial, dgamma, c1, c2
        integer :: i

        call array2symtensor(mat%dimen, mat%nvoigt, strain, Tstrain)
        call array2symtensor(mat%dimen, mat%nvoigt, pls, Tpls)
        call array2symtensor(mat%dimen, mat%nvoigt, b, Tb)
        call eval_trace(mat%dimen, Tstrain, traceStrain)
        e_trial = Tstrain
        do i = 1, mat%dimen
            e_trial(i, i) = e_trial(i, i) - 1.d0/3.d0*traceStrain
        end do

        ! Compute trial stress
        s_trial = 2*mat%mu*(e_trial - Tpls)

        ! Compute shifted stress
        eta_trial = s_trial - Tb

        ! Check yield condition
        call eval_frobenius_norm(mat%dimen, eta_trial, norm_trial)
        f_trial = norm_trial - sqrt(2.d0/3.d0)*(mat%sigma_Y + mat%beta*mat%H*a)  
        sigma = s_trial
        do i = 1, mat%dimen
            sigma(i, i) = sigma(i, i) + mat%bulk*traceStrain
        end do
        kwargs = 0.d0; kwargs(1) = mat%lambda; kwargs(2) = mat%mu; 
        
        if (f_trial.le.0.d0) then ! Elastic point
            pls_new = pls; a_new = a; b_new = b  
            call symtensor2array(mat%dimen, mat%nvoigt, sigma, stress)
            return
        end if

        ! Compute df/dsigma
        NN = eta_trial/norm_trial
        call symtensor2array(mat%dimen, mat%nvoigt, NN, kwargs(4:))

        ! Compute plastic-strain increment        
        dgamma = 3.d0*f_trial/(6.d0*mat%mu+2.d0*mat%H)

        ! Update stress
        sigma = sigma - 2*mat%mu*dgamma*NN
        call symtensor2array(mat%dimen, mat%nvoigt, sigma, stress)

        ! Update plastic strain
        Tpls = Tpls + dgamma*NN
        call symtensor2array(mat%dimen, mat%nvoigt, Tpls, pls_new)
        
        ! Update internal hardening variable
        a_new = a + sqrt(2.d0/3.d0)*dgamma
        
        ! Update back stress
        Tb = Tb + 2.d0/3.d0*(1 - mat%beta)*mat%H*dgamma*NN
        call symtensor2array(mat%dimen, mat%nvoigt, tb, b_new)

        ! Compute some coefficients
        c1 = 2.d0*mat%mu*dgamma/norm_trial
        c2 = 1.d0/(1.d0 + mat%H/(3.d0*mat%mu)) - c1
        kwargs(1) = mat%lambda + 2.d0*mat%mu*c1/3.d0
        kwargs(2) = mat%mu*(1.d0 - c1)
        kwargs(3) = -2.d0*mat%mu*c2

    end subroutine returnMappingAlgorithm

    subroutine setup_kwargs(mat, nc, kwargs)
        !! Points to data of the mechanical behavior 

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        integer, intent(in) :: nc
        double precision, target, intent(in) :: kwargs
        dimension :: kwargs(mat%nvoigt+3, nc)

        ! Local data
        ! ----------
        double precision :: NN(mat%dimen, mat%dimen)
        integer :: i

        mat%kwargs => kwargs(:3, :)
        mat%nn     => kwargs(4:, :)
        mat%JJnn   = 0.d0
        if (mat%isElastic) return
        do i = 1, nc
            call array2symtensor(mat%dimen, mat%nvoigt, mat%nn(:, i), NN)
            mat%JJnn(:,:,i) = matmul(mat%invJJ(:,:,i), NN)
        end do

    end subroutine setup_kwargs

    subroutine compute_mean_plasticity_3d(mat, nc_u, nc_v, nc_w)
        !! Computes the average of the material properties (for the moment it only considers elastic materials)

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2, samplesize = 3**3
        type(mecamat), pointer :: mat
        integer, intent(in) :: nc_u, nc_v, nc_w

        ! Local data
        ! ----------
        integer :: i, j, k, l, genPos, pos
        integer :: ind_u, ind_v, ind_w, sample
        dimension :: ind_u(3), ind_v(3), ind_w(3), sample(samplesize)
        double precision :: DD, coefs, nn, Tnn, mean
        dimension ::    DD(dimen, samplesize), coefs(dimen, dimen, samplesize), nn(nvoigt), &
                        Tnn(dimen, dimen, samplesize), mean(3)
        
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

        do i = 1, samplesize
            nn = mat%nn(sample(i), :)
            call array2symtensor(dimen, nvoigt, nn, Tnn(:,:,i))
        end do

        do i = 1, dimen
            DD = 0.d0
            do j = 1, dimen
                DD(j, :) = mat%kwargs(2, sample) + mat%kwargs(3, sample)*(Tnn(i, j, :)**2)
            end do
            DD(i, :) = DD(i, :) + mat%kwargs(1, sample) + mat%kwargs(2, sample) 
            
            do k = 1, samplesize
                l = sample(k)
                call gemm_AWB(1, 3, 3, mat%invJJ(:,:,l), 3, 3, mat%invJJ(:,:,l), DD(:, k), 3, 3, coefs(:,:,k))
            end do

            call compute_mean_3d(3, 3, 3, coefs(1, 1, :), mean(1))
            call compute_mean_3d(3, 3, 3, coefs(2, 2, :), mean(2))
            call compute_mean_3d(3, 3, 3, coefs(3, 3, :), mean(3))
            mat%mean(i, :) = mean

        end do

    end subroutine compute_mean_plasticity_3d

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
        integer :: i, j, k, l, r, alpha, beta, zeta!, info
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

                do r = 1, dimen
                    kt1(r, :) = mat%kwargs(r, :)*t1*mat%detJJ
                end do

                t2 = kt1(1, :)*mat%invJJ(l, j, :)
                t4 = kt1(3, :)*mat%JJnn(l, j, :)

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJJ(l, i, :)
                    t7 = 0.d0

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta  = beta + (alpha - 1)*2
                        t5    = t2*mat%invJJ(k, i, :) + t3*mat%invJJ(k, j, :) + t4*mat%JJnn(k, i, :)
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

    subroutine wq_forceint_3d(mat, stress, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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
            t1(:,:,i) = matmul(mat%invJJ(:,:,i), Tstress)*mat%detJJ(i)
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

    end subroutine wq_forceint_3d
    
    subroutine mf_wq_linearelasticity_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                U_u, U_v, U_w, Deigen, ndu, ndv, ndw, dod_u, dod_v, dod_w, b, &
                                nbIterPCG, threshold, isPrecond, x, resPCG)
        !! Solves elasticity problems using (Preconditioned) Bi-Conjugate gradient method
        !! This algorithm solve S x = F, where S is the stiffness matrix
        !! Moreover, it considers Dirichlet boundaries are zero
        !! IN CSR FORMAT    
    
        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3, nvoigt = dimen*(dimen+1)/2  
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        type(mecamat), pointer :: mat
        integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                        data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)
        double precision, intent(in) :: U_u, U_v, U_w, Deigen
        dimension :: U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), U_w(nr_w, nr_w, dimen), Deigen(dimen, nr_total)
        integer, intent(in) :: ndu, ndv, ndw
        integer, intent(in) :: dod_u, dod_v, dod_w
        dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)
    
        integer, intent(in) :: nbIterPCG
        double precision, intent(in) :: threshold
        logical, intent(in) :: isPrecond 
    
        double precision, intent(in) :: b
        dimension :: b(dimen, nr_total)
        
        double precision, intent(out) :: x, resPCG
        dimension :: x(dimen, nr_total), resPCG(nbIterPCG+1)
    
        ! Local data
        ! ----------
        double precision :: rsold, rsnew, alpha, omega, beta, prod, prod2, normb
        double precision, dimension(dimen, nr_total) :: r, rhat, p, s, Aptilde, Astilde
        double precision, allocatable, dimension(:, :) :: ptilde, stilde
        integer :: iter
    
        integer :: indi_T_u, indi_T_v, indi_T_w, indj_T_u, indj_T_v, indj_T_w
        dimension ::    indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1), &
                        indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)
        
        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
        call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
        call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)
    
        x = 0.d0; r = b
        call reset_dirichletbound3(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
        rhat = r; p = r
        call block_dot_product(dimen, nr_total, r, rhat, rsold)
        normb = maxval(abs(r))
        resPCG = 0.d0; resPCG(1) = 1.d0
        if (normb.lt.threshold) return
    
        if (.not.isPrecond) then ! Conjugate gradient algorithm
            do iter = 1, nbIterPCG
                
                call mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, p, Aptilde)
                call reset_dirichletbound3(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
                
                call block_dot_product(dimen, nr_total, Aptilde, rhat, prod)
                alpha = rsold/prod
                s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
    
                call mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, s, Astilde)
                call reset_dirichletbound3(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    
                call block_dot_product(dimen, nr_total, Astilde, s, prod)
                call block_dot_product(dimen, nr_total, Astilde, Astilde, prod2)
                omega = prod/prod2
                x = x + alpha*p + omega*s ! Normally x is alrady Dirichlet updated
                r = s - omega*Astilde ! Normally r is alrady Dirichlet updated
    
                resPCG(iter+1) = maxval(abs(r))/normb
                if (resPCG(iter+1).le.threshold) exit
                call block_dot_product(dimen, nr_total, r, rhat, rsnew)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                rsold = rsnew
            end do
    
        else ! Preconditioned Conjugate Gradient algorithm
    
            allocate(ptilde(dimen, nr_total), stilde(dimen, nr_total))
            do iter = 1, nbIterPCG
                call fd_linearelasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, p, ptilde)
                call reset_dirichletbound3(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    
                call mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
                call reset_dirichletbound3(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    
                call block_dot_product(dimen, nr_total, Aptilde, rhat, prod)
                alpha = rsold/prod
                s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
                
                call fd_linearelasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, s, stilde)
                call reset_dirichletbound3(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    
                call mf_wq_stiffness_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                        nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                        data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        data_W_u, data_W_v, data_W_w, stilde, Astilde)
                call reset_dirichletbound3(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    
                call block_dot_product(dimen, nr_total, Astilde, s, prod)
                call block_dot_product(dimen, nr_total, Astilde, Astilde, prod2)
    
                omega = prod/prod2
                x = x + alpha*ptilde + omega*stilde ! Normally x is alrady Dirichlet updated
                r = s - omega*Astilde ! Normally r is alrady Dirichlet updated
                
                resPCG(iter+1) = maxval(abs(r))/normb
                if (resPCG(iter+1).le.threshold) exit
                call block_dot_product(dimen, nr_total, r, rhat, rsnew)
                beta = (alpha/omega)*(rsnew/rsold)
                p = r + beta*(p - omega*Aptilde)
                rsold = rsnew
    
            end do
    
        end if
    
    end subroutine mf_wq_linearelasticity_3d
    
end module linearelastoplasticity