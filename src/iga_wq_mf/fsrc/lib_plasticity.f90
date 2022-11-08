! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
!
! Remarks :: Voigt notation is used. That is 
! Strain e = [e11, e22, e33, 2e12, 2e13, 2e23]^T
! Stress s = [s11, s22, s33, s12, s13, s23]^T
! 1 = [1, 1, 1, 0, 0, 0]^T
! I = diag(1, 1, 1, 0.5, 0.5, 0.5)
! And then, some products are:
! Double contraction st-st s:e = dot(s, e), this only works with stress and strains
! In Voigt notation if one wants to compute e:e or s:s some scale factor need to be considered.
! Double contraction C:e = matmul(C, e), this is true only with strains
! Otherwise some scale factor need to be considered.
! ==========================

subroutine clean_dirichlet_3d(nr, A, ndu, ndv, ndw, dod_u, dod_v, dod_w)
    !! Set to 0 (Dirichlet condition) the values of an array using the dod indices in each dimension
    !! A is actually a vector arranged following each dimension [Au, Av, Aw]

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, ndu, ndv, ndw
    double precision, intent(inout) :: A
    dimension :: A(3, nr)

    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    A(1, dod_u) = 0.d0 
    A(2, dod_v) = 0.d0 
    A(3, dod_w) = 0.d0 

end subroutine clean_dirichlet_3d

subroutine block_dot_product(dimen, nr, A, B, result)
    !! Computes dot product of A and B. Both are actually vectors arranged following each dimension
    !! Vector A is composed of [Au, Av, Aw] and B of [Bu, Bv, Bw]. 
    !! Dot product A.B = Au.Bu + Av.Bv + Aw.Bw 

    implicit none
    ! Input/ output data
    ! ------------------
    integer, intent(in) :: dimen, nr
    double precision, intent(in) :: A, B
    dimension :: A(dimen, nr), B(dimen, nr)

    double precision :: result

    ! Local data
    ! ----------
    integer :: i
    double precision :: rtemp

    result = 0.d0
    do i = 1, dimen 
        rtemp = dot_product(A(i, :), B(i, :))
        result = result + rtemp
    end do

end subroutine block_dot_product

subroutine compute_stress_deviatoric(dimen, ddl, tensor, dev)
    !! Returns deviatoric of a second-order tensor. That is
    !! dev(s) = s - 1/3 trace(s) 1. 
    !! s is a tensor written in Voigt notation

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor
    dimension :: tensor(ddl)

    double precision, intent(out) :: dev
    dimension :: dev(ddl)

    ! Local data
    ! ----------
    integer :: i
    double precision :: trace, one
    dimension :: one(ddl)

    one = 0.d0
    do i = 1, dimen
        one(i) = 1.d0
    end do

    trace = 0.d0
    do i = 1, dimen
        trace = trace + tensor(i)
    end do    

    dev = tensor - 1.d0/3.d0*trace*one 
    
end subroutine compute_stress_deviatoric

subroutine compute_stress_norm(dimen, ddl, tensor, norm)
    !! Returns Frobenius norm of a second-order stress-like tensor 
    !! The definition is norm = sqrt(s:s). Since s is a tensor written in Voigt notation
    !! some scale factor must be considered.

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor
    dimension :: tensor(ddl)

    double precision, intent(out) :: norm

    ! Local data
    ! ----------
    integer :: i

    norm = 0.0d0

    do i = 1, dimen
        norm = norm + tensor(i)*tensor(i)
    end do

    do i = dimen+1, ddl
        norm = norm + 2.d0*tensor(i)*tensor(i)
    end do

    norm = sqrt(norm)
    
end subroutine compute_stress_norm

subroutine fourth_order_identity(dimen, ddl, identity)
    !! Creates a fourth-order identity (in Voigt notation)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl 
    double precision, intent(out) :: identity
    dimension :: identity(ddl, ddl)

    ! Local data
    ! ----------
    integer :: i

    identity = 0.d0
    do i = 1, dimen
        identity(i, i) = 1.d0
    end do

    do i = dimen+1, ddl
        identity(i, i) = 0.5d0
    end do

end subroutine fourth_order_identity

subroutine one_kron_one(dimen, ddl, onekronone)
    !! Creates a one kron one tensor (in Voigt notation)

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(out) :: onekronone
    dimension :: onekronone(ddl, ddl)

    ! Local data
    ! ----------
    integer :: i, j

    onekronone = 0.d0
    do j = 1, dimen
        do i = 1, dimen
            onekronone(i, j) = 1.d0
        end do
    end do

end subroutine one_kron_one

subroutine create_incidence_matrix(dimen, ddl, EE)
    !! Creates incidence matrix E. E is the passage matrix from derivative to actual symetric values. 
    !! If we multiply a vector of values u_(i, j) with M matrix, one obtains the vector: 
    !! us_ij = 0.5*(u_(i,j) + u_(j,i))  

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(out) :: EE
    dimension :: EE(ddl, dimen, dimen)

    if (dimen.eq.3) then 
        EE = 0.d0
        EE(1, 1, 1) = 1.d0; EE(5, 3, 1) = 1.d0; EE(6, 2, 1) = 1.d0
        EE(2, 2, 2) = 1.d0; EE(4, 3, 2) = 1.d0; EE(6, 1, 2) = 1.d0
        EE(3, 3, 3) = 1.d0; EE(4, 2, 3) = 1.d0; EE(5, 1, 3) = 1.d0

    else if (dimen.eq.2) then
        EE = 0.d0
        EE(1, 1, 1) = 1.d0; EE(3, 2, 1) = 1.d0
        EE(2, 2, 2) = 1.d0; EE(3, 1, 2) = 1.d0

    else 
        stop 'There is no incende matrix in 1D'
    end if

end subroutine create_incidence_matrix

subroutine compute_stress_vonmises(dimen, ddl, tensor, VM)
    !! Computes equivalent stress with Von Mises formula of a second-order stress-like tensor.
    !! Given tensor is written in Voigt notation

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor
    dimension :: tensor(ddl)

    double precision, intent(out) :: VM

    ! Local data
    ! ----------
    double precision :: dev_tensor
    dimension :: dev_tensor(ddl)

    call compute_stress_deviatoric(dimen, ddl, tensor, dev_tensor)
    call compute_stress_norm(dimen, ddl, tensor, VM)
    VM = sqrt(3.d0/2.d0)*VM

end subroutine compute_stress_vonmises

module elastoplasticity

    implicit none
    type :: mecamat
        ! Inputs 
        integer :: dimen=3, ddl=6
        double precision :: young, hardening, beta, poisson, sigma_Y

        ! Outputs
        double precision :: lambda, mu, bulk
        double precision, dimension(:,:), allocatable :: Ctensor, Stensor

        ! Local
        double precision, dimension(:,:), allocatable :: Idev
    
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

        ! Local data
        ! ----------
        double precision :: lambda, mu, bulk
        double precision, dimension(:, :), allocatable :: II, Idev, onekronone, CC, SS

        allocate(mat)
        mat%dimen = dimen
        mat%ddl = dimen*(dimen+1)/2
        allocate(II(mat%ddl, mat%ddl), Idev(mat%ddl, mat%ddl), &
                onekronone(mat%ddl, mat%ddl), CC(mat%ddl, mat%ddl), SS(mat%ddl, mat%ddl))

        ! Create special tensors in Voigt notation
        call fourth_order_identity(mat%dimen, mat%ddl, II)
        call one_kron_one(mat%dimen, mat%ddl, onekronone)
        Idev = II - 1.d0/3.d0*onekronone

        ! Compute mechanical properties
        lambda = nu*E/((1+nu)*(1-2*nu))
        mu = E/(2*(1+nu))
        bulk = lambda + 2.d0/3.d0*mu
        CC = lambda*onekronone + 2*mu*II
        SS = 1.d0/(9.d0*bulk)*onekronone + 1.d0/(2.d0*mu)*(Idev)

        ! Save data computed
        allocate(mat%Ctensor(mat%ddl, mat%ddl), mat%Stensor(mat%ddl, mat%ddl), mat%Idev(mat%ddl, mat%ddl))
        mat%young = E
        mat%hardening = H
        mat%beta = beta
        mat%poisson = nu
        mat%sigma_Y = sigma_Y
        mat%lambda = lambda
        mat%mu = mu
        mat%bulk = bulk
        mat%Ctensor = CC
        mat%Stensor = SS
        mat%Idev = Idev

    end subroutine initialize_mecamat

    subroutine yield_combined_hardening(mat, sigma, alpha, ep, f, norm, NN)
        !! Computes the value of f (consistency condition) in combined hardening criteria
        
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        double precision, intent(in) :: sigma, alpha, ep
        dimension :: sigma(mat%ddl), alpha(mat%ddl)
        double precision, intent(out) :: f, norm, NN
        dimension :: NN(mat%ddl)

        ! Local data
        ! ----------
        double precision :: eta, dev_sigma
        dimension :: eta(mat%ddl), dev_sigma(mat%ddl)

        call compute_stress_deviatoric(mat%dimen, mat%ddl, sigma, dev_sigma)
        eta = dev_sigma - alpha
        call compute_stress_norm(mat%dimen, mat%ddl, eta, norm)
        f = norm - sqrt(2.d0/3.d0) * (mat%sigma_Y + (1-mat%beta)*mat%hardening*ep)   

        NN = 0.d0
        if (norm.gt.0.d0) NN = 1.d0/norm*eta

    end subroutine yield_combined_hardening

    subroutine cpp_combined_hardening(mat, deps, alpha_n0, ep_n0, sigma_n0, &
                                    alpha_n1, ep_n1, sigma_n1, Dalg)
        !! Return closest point proyection (cpp) in combined hardening criteria

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        double precision, intent(in) :: deps, alpha_n0, ep_n0, sigma_n0
        dimension :: alpha_n0(mat%ddl), deps(mat%ddl), sigma_n0(mat%ddl)

        double precision, intent(out) ::alpha_n1, ep_n1, sigma_n1, Dalg
        dimension :: alpha_n1(mat%ddl), sigma_n1(mat%ddl), Dalg(mat%ddl, mat%ddl)
        
        ! Local data
        ! ----------
        integer :: i, j
        double precision :: f, norm, dgamma, c1, c2
        double precision :: sigma_trial, N, NNT
        dimension :: sigma_trial(mat%ddl), N(mat%ddl), NNT(mat%ddl, mat%ddl)

        sigma_trial = sigma_n0 + matmul(mat%Ctensor, deps)
        call yield_combined_hardening(mat, sigma_trial, alpha_n0, ep_n0, f, norm, N)

        if (f.le.0.d0) then ! Elastic point

            sigma_n1 = sigma_trial ! stress
            ep_n1 = ep_n0 ! effective plastic strain
            alpha_n1 = alpha_n0 ! back stress
            Dalg = mat%Ctensor ! tangent matrix

        else ! Plastic point

            dgamma = f/(2.d0*mat%mu+2.d0*mat%hardening/3.d0) ! consistency parameter
            sigma_n1 = sigma_trial - 2.d0*mat%mu*dgamma*N
            ep_n1 = ep_n0 + sqrt(2.d0/3.d0)*dgamma
            alpha_n1 = alpha_n0 + 2.d0/3.d0*mat%beta*mat%hardening*dgamma*N

            ! Compute consistent tangent matrix
            c1 = 4.d0*mat%mu**2.d0/(2.d0*mat%mu+2.d0/3.d0*mat%hardening)
            c2 = 4.d0*mat%mu**2.d0*dgamma/norm
            do i = 1, mat%ddl
                do j = 1, mat%ddl
                    NNT(i, j) = N(i)*N(j)
                end do
            end do
            Dalg = mat%Ctensor - (c1-c2)*NNT - c2*mat%Idev

        end if

    end subroutine cpp_combined_hardening

    subroutine compute_meca_coefficients(dimen, ddl, nc_total, sigma, DD, invJ, detJ, coefs_Fint, coefs_Stiff)
        !! Computes the coefficients to use in internal force vector and stiffness matrix
        !! If a pseudo Newton-Raphson method is used, one must compute coef_Fint and coef_Stiff separately

        implicit none 
        ! Input / output data
        ! -------------------
        integer, intent(in) :: dimen, ddl, nc_total
        double precision, intent(in) :: sigma, DD, invJ, detJ
        dimension :: sigma(ddl, nc_total), DD(ddl, ddl, nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: coefs_Fint, coefs_Stiff
        dimension :: coefs_Fint(dimen*dimen, nc_total), coefs_Stiff(dimen*dimen, dimen*dimen, nc_total)

        ! Local data
        ! ----------
        integer :: i, j, k
        double precision :: EE, ETCE, Dij, Si
        dimension :: EE(ddl, dimen, dimen), ETCE(dimen, dimen), Dij(dimen, dimen), Si(dimen)
                    
        call create_incidence_matrix(dimen, ddl, EE)

        do k = 1, nc_total
            do i = 1, dimen
                do j = 1, dimen
                    ! Compute stiffness coefficients    
                    ETCE = matmul(transpose(EE(:,:,i)), matmul(DD(:,:,k), EE(:,:,j)))
                    Dij = matmul(invJ(:,:,k), matmul(ETCE, transpose(invJ(:,:,k))))
                    coefs_Stiff((i-1)*dimen+1:i*dimen, (j-1)*dimen+1:j*dimen, k) = Dij*detJ(k)

                end do

                ! Compute stress coefficients
                Si = matmul(invJ(:,:,k), matmul(transpose(EE(:,:,i)), sigma(:,k)))
                coefs_Fint((i-1)*dimen+1:i*dimen, k) = Si*detJ(k)

            end do
        end do

    end subroutine compute_meca_coefficients

end module elastoplasticity

subroutine mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
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
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen*dimen, dimen*dimen, nc_total)

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
    type(thermomat), pointer :: mat
    double precision :: array_temp
    dimension :: array_temp(nr_total)
    integer :: i, j
    
    array_out = 0.d0
    do i = 1, dimen
        do j = 1, dimen

            call setupKcoefs(mat, dimen, nc_total, coefs((i-1)*dimen+1:i*dimen, (j-1)*dimen+1:j*dimen, :))
            call mf_wq_get_ku_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, array_in(j, :), array_temp)
            array_out(i, :) = array_out(i, :) + array_temp    

        end do 
    end do

end subroutine mf_wq_get_su_3d

subroutine fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, eigen_diag, array_in, array_out)
    !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
    !! Applied to elasticity problems
    !! by G. Sanaglli and M. Tani
    
    use heat_solver
    implicit none
    ! Input / output  data 
    !---------------------
    integer, parameter :: dimen = 3
    integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
    double precision, intent(in) :: U_u, U_v, U_w, eigen_diag, array_in
    dimension ::    U_u(nr_u, nr_u, dimen), U_v(nr_v, nr_v, dimen), U_w(nr_w, nr_w, dimen), &
                    eigen_diag(dimen, nr_total), array_in(dimen, nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(dimen, nr_total)

    ! Local data
    ! ----------
    type(cgsolver), pointer :: solv
    double precision :: array_temp
    dimension :: array_temp(nr_total)
    integer :: i

    allocate(solv)
    do i = 1, dimen 
        call setup_eigendiag(solv, nr_total, eigen_diag(i, :))
        call fast_diagonalization(solv, nr_total, nr_u, nr_v, nr_w, U_u(:, :, i), U_v(:, :, i), U_w(:, :, i), &
                                array_in(i, :), array_temp)
        array_out(i, :) = array_temp
    end do
    
end subroutine fd_elasticity_3d

subroutine mf_wq_elasticity_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
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
    integer, parameter :: dimen = 3    
    integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(dimen*dimen, dimen*dimen, nc_total)
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

    x = 0.d0; r = b; resPCG = 0.d0; resPCG(1) = 1.d0
    call clean_dirichlet_3d(nr_total, r, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
    rhat = r; p = r
    call block_dot_product(dimen, nr_total, r, rhat, rsold)
    normb = maxval(abs(r))

    if (.not.isPrecond) then ! Conjugate gradient algorithm
        
        do iter = 1, nbIterPCG
            call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, p, Aptilde)
            call clean_dirichlet_3d(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call block_dot_product(dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated

            call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, s, Astilde)
            call clean_dirichlet_3d(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

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
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, p, ptilde)
            call clean_dirichlet_3d(nr_total, ptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
            call clean_dirichlet_3d(nr_total, Aptilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 
            
            call block_dot_product(dimen, nr_total, Aptilde, rhat, prod)
            alpha = rsold/prod
            s = r - alpha*Aptilde ! Normally s is alrady Dirichlet updated
            
            call fd_elasticity_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, Deigen, s, stilde)
            call clean_dirichlet_3d(nr_total, stilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

            call mf_wq_get_su_3d(coefs, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                    nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, stilde, Astilde)
            call clean_dirichlet_3d(nr_total, Astilde, ndu, ndv, ndw, dod_u, dod_v, dod_w) 

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

end subroutine mf_wq_elasticity_3d
