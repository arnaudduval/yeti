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

subroutine clean_dirichlet_3dim(nr, A, ndu, ndv, ndw, dod_u, dod_v, dod_w)
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

end subroutine clean_dirichlet_3dim

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
    integer, parameter :: dimen = 3
    integer, parameter :: ddl = dimen*(dimen+1)/2
    double precision, parameter :: tol1 = 1e-6, tol2 = 1e-6

    type :: mecamat
        ! Inputs 
        double precision :: young, poisson, sigma_Y, hardening, beta

        ! Outputs
        double precision :: lambda, mu, bulk
        double precision, dimension(:, :), allocatable :: Ctensor, Stensor

        ! Local
        double precision, dimension(:, :), allocatable :: Idev
    
    end type mecamat

contains

    ! ======================================
    ! COMBINED ISOTROPIC/KINEMATIC HARDENING
    ! ======================================

    subroutine initialize_mecamat(object, E, H, beta, nu, sigma_Y)
        !! Creates a material using combined isotropic/kinematic hardening theory 

        implicit none
        ! Input / output data
        ! -------------------
        double precision, intent(in) :: E, H, beta, nu, sigma_Y
        type(mecamat), pointer :: object

        ! Local data
        ! ----------
        double precision :: lambda, mu, bulk
        double precision, dimension(ddl, ddl) :: II, Idev, onekronone, CC, SS

        ! Create special tensors in Voigt notation
        call fourth_order_identity(dimen, ddl, II)
        call one_kron_one(dimen, ddl, onekronone)
        Idev = II - 1.d0/3.d0*onekronone

        ! Compute mechanical properties
        lambda = nu*E/((1+nu)*(1-2*nu))
        mu = E/(2*(1+nu))
        bulk = lambda + 2.d0/3.d0*mu
        CC = lambda*onekronone + 2*mu*II
        SS = 1.d0/(9.d0*bulk)*onekronone + 1.d0/(2.d0*mu)*(Idev)

        ! Save data computed
        allocate(object)
        allocate(object%Ctensor(ddl, ddl), object%Stensor(ddl, ddl), object%Idev(ddl, ddl))
        object%young = E
        object%hardening = H
        object%beta = beta
        object%poisson = nu
        object%sigma_Y = sigma_Y
        object%lambda = lambda
        object%mu = mu
        object%bulk = bulk
        object%Ctensor = CC
        object%Stensor = SS
        object%Idev = Idev

    end subroutine initialize_mecamat

    subroutine compute_meca_coefficients(nc_total, sigma, DD, invJ, detJ, coef_Fint, coef_Stiff)
        !! Computes the coefficients to use in internal force vector and stiffness matrix
        !! If a pseudo Newton-Raphson method is used, one must compute coef_Fint and coef_Stiff separately

        implicit none 
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nc_total
        double precision, intent(in) :: sigma, DD, invJ, detJ
        dimension :: sigma(ddl, nc_total), DD(ddl, ddl, nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: coef_Fint, coef_Stiff
        dimension :: coef_Fint(dimen*dimen, nc_total), coef_Stiff(dimen*dimen, dimen*dimen, nc_total)

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
                    coef_Stiff((i-1)*dimen+1:i*dimen, (j-1)*dimen+1:j*dimen, k) = Dij*detJ(k)

                end do

                ! Compute stress coefficients
                Si = matmul(invJ(:,:,k), matmul(transpose(EE(:,:,i)), sigma(:,k)))
                coef_Fint((i-1)*dimen+1:i*dimen, k) = Si*detJ(k)

            end do
        end do

    end subroutine compute_meca_coefficients

    subroutine yield_combined_hardening(sigma_Y, H, beta, sigma, alpha, ep, f, norm, NN)
        !! Computes the value of f (consistency condition) in combined hardening criteria
        
        implicit none
        ! Input / output data
        ! -------------------
        double precision, intent(in) :: sigma_Y, beta, H, sigma, alpha, ep
        dimension :: sigma(ddl), alpha(ddl)
        double precision, intent(out) :: f, norm, NN
        dimension :: NN(ddl)

        ! Local data
        ! ----------
        double precision :: eta, dev_sigma
        dimension :: eta(ddl), dev_sigma(ddl)

        call compute_stress_deviatoric(dimen, ddl, sigma, dev_sigma)
        eta = dev_sigma - alpha
        call compute_stress_norm(dimen, ddl, eta, norm)
        f = norm - sqrt(2.d0/3.d0) * (sigma_Y + (1-beta)*H*ep)   

        if (norm.gt.0.d0) then
            NN = 1.d0/norm*eta
        else
            NN = 0.d0
        end if

    end subroutine yield_combined_hardening

    subroutine cpp_combined_hardening(mat, deps, alpha_n0, ep_n0, sigma_n0, &
                                    alpha_n1, ep_n1, sigma_n1, Dalg)
        !! Return closest point proyection (cpp) in combined hardening criteria

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat), pointer :: mat
        double precision, intent(in) :: deps, alpha_n0, ep_n0, sigma_n0
        dimension :: alpha_n0(ddl), deps(ddl), sigma_n0(ddl)

        double precision, intent(out) ::alpha_n1, ep_n1, sigma_n1, Dalg
        dimension :: alpha_n1(ddl), sigma_n1(ddl), Dalg(ddl, ddl)
        
        ! Local data
        ! ----------
        integer :: i, j
        double precision :: f, norm, dgamma, c1, c2
        double precision :: sigma_trial, N, NNT
        dimension :: sigma_trial(ddl), N(ddl), NNT(ddl, ddl)

        sigma_trial = sigma_n0 + matmul(mat%Ctensor, deps)
        call yield_combined_hardening(mat%sigma_Y, mat%hardening, mat%beta, sigma_trial, alpha_n0, ep_n0, f, norm, N)

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
            do i = 1, ddl
                do j = 1, ddl
                    NNT(i, j) = N(i)*N(j)
                end do
            end do
            Dalg = mat%Ctensor - (c1-c2)*NNT - c2*mat%Idev

        end if

    end subroutine cpp_combined_hardening

end module elastoplasticity
