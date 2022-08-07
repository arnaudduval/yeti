! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
! ==========================

subroutine block_dot_product(dimen, nc, A, B, result)
    !! Computes dot product of A and B. Both are actually vectors arranged following each dimension
    !! Vector A is composed of [Au, Av, Aw] and B of [Bu, Bv, Bw]. 
    !! Dot product A.B = Au.Bu + Av.Bv + Aw.Bw 

    implicit none
    ! Input/ output
    ! -----------------
    integer, intent(in) :: dimen, nc
    double precision, intent(in) :: A, B
    dimension :: A(dimen, nc), B(dimen, nc)

    double precision :: result

    ! Local data
    ! ------------
    integer :: i
    double precision :: rtemp

    ! Initialize
    result = 0.d0

    do i = 1, dimen 
        rtemp = dot_product(A(i, :), B(i, :))
        result = result + rtemp
    end do

end subroutine block_dot_product

subroutine array_maps_tensor(dimen, p, i, j)
    !! Converts an array index into a symmetric second-order tensor index

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dimen, p 
    integer, intent(out) :: i, j

    if (dimen.eq.3) then 
        if (p.le.dimen) then 
            i = p; j = p;
        else
            if (p.eq.4) then
                i = 1; j = 2
            else if (p.eq.5) then
                i = 1; j = 3;
            else if (p.eq.6) then
                i = 2; j = 3;
            end if
        end if
    else if (dimen.eq.2) then 
        if (p.le.dimen) then 
            i = p; j = p;
        else
            if (p.eq.3) then
                i = 1; j = 2
            end if
        end if
    else
        print*, "Error. Only 2d or 3d"
    end if

end subroutine array_maps_tensor

subroutine tensor_maps_array(dimen, i, j, p)
    !! Converts a symmetric second-order tensor index into an array index

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dimen, i, j
    integer, intent(out) :: p

    if ((dimen.eq.3).or.(dimen.eq.2)) then 
        if (i.eq.j) then
            p = i
        else 
            p = i + j + dimen - 2
        end if
    else
        print*, "Error. Only 2d or 3d"
    end if

end subroutine tensor_maps_array

subroutine array2st(dimen, ddl, array, tensor)
    !! Returns second-order tensor from array
    !! i.e. from [t11, t22, t12], one gets [[t11, t12], [t21, t22]] (with t21 = t12)

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: array
    dimension :: array(ddl)

    double precision, intent(out) :: tensor
    dimension :: tensor(dimen, dimen)

    ! Local data
    ! -------------
    integer :: i, j, p
    double precision :: diag
    dimension :: diag(dimen, dimen)

    ! Initialize
    tensor = 0.d0
    diag = 0.d0

    do p = 1, ddl
        call array_maps_tensor(dimen, p, i, j)
        tensor(i, j) = array(p) 
    end do

    ! Find diagonal of tensor
    do i = 1, dimen
        diag(i, i) = tensor(i, i)
    end do

    ! Compute real tensor
    tensor = tensor + transpose(tensor) - diag

end subroutine array2st

subroutine stdcst(dimen, ddl, A, B, result)
    !! Returns A double contracted with B
    !! A and B are second order tensor in Voigt representation
    !! With inditial notation : result = A_ij * B_ij
    
    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: A, B
    dimension :: A(ddl), B(ddl)

    double precision, intent(out) :: result

    ! Local data
    ! ---------------
    integer :: i, j
    double precision :: At, Bt
    dimension :: At(ddl, ddl), Bt(ddl, ddl)

    ! Convert array to tensor
    call array2st(dimen, ddl, A, At)
    call array2st(dimen, ddl, B, Bt)

    ! Compute double contracted result
    result = 0.d0
    do j = 1, ddl
        do i = 1, ddl
            result = result + At(i, j)*Bt(i, j)
        end do
    end do

end subroutine stdcst

subroutine stkronst(ddl, A, B, result)
    !! Returns kron product of tensors A and B
    !! A and B are second order tensor in Voigt representation

    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: ddl
    double precision, intent(in) :: A, B
    dimension :: A(ddl), B(ddl)

    double precision, intent(out) :: result
    dimension :: result(ddl, ddl)

    ! Local data
    ! ---------------
    integer :: i, j

    do j = 1, ddl
        do i = 1, ddl
            result(i, j) = A(i)*B(j)
        end do
    end do

end subroutine stkronst

subroutine compute_deviatoric(dimen, ddl, tensor, dev)
    !! Returns deviatoric of a second-order tensor 

    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor
    dimension :: tensor(ddl)

    double precision, intent(out) :: dev
    dimension :: dev(ddl)

    ! Local data
    ! -------------
    integer :: i
    double precision :: trace, one
    dimension :: one(ddl)

    ! Compute trace of tensor
    trace = 0.d0
    do i = 1, dimen
        trace = trace + tensor(i)
    end do

    ! Compute one
    one = 0.d0
    do i = 1, dimen
        one(i) = 1.d0
    end do

    ! Definition of deviatoric
    dev = tensor - 1.d0/3.d0*trace*one 
    
end subroutine compute_deviatoric

subroutine clean_dirichlet_3d(nc, A, ndu, ndv, ndw, dod_u, dod_v, dod_w)
    !! Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
    !! A is actually a vector arranged following each dimension [Au, Av, Aw]

    implicit none
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nc, ndu, ndv, ndw
    double precision, intent(inout) :: A
    dimension :: A(3, nc)

    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    ! Clean array
    A(1, dod_u) = 0.d0 
    A(2, dod_v) = 0.d0 
    A(3, dod_w) = 0.d0 

end subroutine clean_dirichlet_3d

subroutine fourth_order_identity(ddl, identity)
    !! Creates a fourth-order identity (Voigt representation)

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: ddl 
    double precision, intent(out) :: identity
    dimension :: identity(ddl, ddl)

    ! Local data
    ! -----------
    integer :: i

    ! Create identity
    identity = 0.d0
    do i = 1, ddl
        identity(i, i) = 1.d0
    end do

end subroutine fourth_order_identity

subroutine one_kron_one(dimen, ddl, onekronone)
    !! Creates a one kron one tensor (Voigt representation)

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(out) :: onekronone
    dimension :: onekronone(ddl, ddl)

    ! Local data
    ! -----------
    integer :: i, j

    ! Create one kron one
    onekronone = 0.d0
    do j = 1, dimen
        do i = 1, dimen
            onekronone(i, j) = 1.d0
        end do
    end do

end subroutine one_kron_one

subroutine create_der2sym(dimen, ddl, MM)
    !! Creates M matrix. M is the passage matrix from derivative to actual symetric values. 
    !! If we multiply a vector of values u_(i, j) with M matrix, one obtains the vector: 
    !! us_ij = 0.5*(u_(i,j) + u_(j,i))  

    implicit none 
    ! Input / output 
    ! --------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(out) :: MM
    dimension :: MM(ddl, dimen*dimen)

    if (dimen.eq.3) then 
        MM = 0.0d0
        MM(1, 1) = 1.0d0; MM(2, 5) = 1.0d0; MM(3, 9) = 1.0d0
        MM(4, 2) = 0.5d0; MM(4, 4) = 0.5d0
        MM(5, 6) = 0.5d0; MM(5, 8) = 0.5d0
        MM(6, 3) = 0.5d0; MM(6, 7) = 0.5d0
    else if (dimen.eq.2) then
        MM = 0.0d0
        MM(1, 1) = 1.0d0; MM(2, 4) = 1.0d0
        MM(3, 2) = 0.5d0; MM(3, 3) = 0.5d0
    end if

end subroutine create_der2sym

module elastoplasticity

    implicit none
    integer, parameter :: dimen = 3! By now only consider 3D 
    integer, parameter :: ddl = dimen*(dimen+1)/2
    double precision, parameter :: tol1 = 1e-6, tol2 = 1e-6

    type :: material
        ! Inputs 
        ! ------------
        double precision :: young, poisson

        ! Outputs
        ! -----------
        double precision :: sigma_Y, hardening, beta
        double precision :: lambda, mu, bulk
        double precision, dimension(:, :), allocatable :: Ctensor, Stensor

        ! Local
        ! ---------
        double precision, dimension(:, :), allocatable :: Idev
    
    end type material

    contains

    subroutine initialize_mat(object, E, H, beta, nu, sigma_Y)
        !! Creates a material with isotropic properties only using Young's module and Poisson coefficient

        implicit none
        ! Input / output
        ! ---------------
        double precision, intent(in) :: E, H, beta, nu, sigma_Y
        type(material), pointer :: object

        ! Local data
        ! ----------------
        double precision :: lambda, mu, bulk
        double precision :: identity, onekronone, CC, SS
        dimension :: identity(ddl, ddl), onekronone(ddl, ddl), CC(ddl, ddl), SS(ddl, ddl)

        ! Compute Lame coefficients
        lambda = nu*E/((1+nu)*(1-2*nu))
        mu = E/(2*(1+nu))
        bulk = lambda + 2.d0/3.d0*mu

        ! Create special tensors in Voigt representation
        call fourth_order_identity(ddl, identity)
        call one_kron_one(dimen, ddl, onekronone)

        ! Computes C and S
        CC = lambda*onekronone + 2*mu*identity
        SS = 1.d0/(9.d0*bulk)*onekronone &
            + 1.d0/(2.d0*mu)*(identity - 1.d0/3.d0*onekronone)

        ! Save data computed
        allocate(object)
        object%young = E
        object%hardening = H
        object%beta = beta
        object%poisson = nu
        object%sigma_Y = sigma_Y
        object%lambda = lambda
        object%mu = mu
        object%bulk = bulk

        allocate(object%Ctensor(ddl, ddl), object%Stensor(ddl, ddl), object%Idev(ddl, ddl))
        object%Ctensor = CC
        object%Stensor = SS
        object%Idev = identity - 1.d0/3.d0*onekronone

    end subroutine initialize_mat

    subroutine compute_coefficients(nc_total, sigma, Dalg, invJ, detJ, coef_Fint, coef_Stiff)
        !! Computes the coefficients to use in internal force vector and stiffness matrix

        implicit none 
        ! Input / output 
        ! --------------------
        integer, intent(in) :: nc_total
        double precision, intent(in) :: sigma, Dalg, invJ, detJ
        dimension :: sigma(ddl, nc_total), Dalg(ddl, ddl, nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: coef_Fint, coef_Stiff
        dimension :: coef_Fint(dimen*dimen, nc_total), coef_Stiff(dimen*dimen, dimen*dimen, nc_total)

        ! Local data
        ! -------------
        integer :: i, j
        double precision :: MM, invJext, MM_S, MM_dSdE, MDM, coef_temp
        dimension ::    MM(ddl, dimen*dimen), invJext(dimen*dimen, dimen*dimen), MM_S(dimen*dimen), &
                        MM_dSdE(dimen*dimen, ddl), MDM(dimen*dimen, dimen*dimen), coef_temp(dimen*dimen, dimen*dimen)
                    
        ! Construct passage matrix
        call create_der2sym(dimen, ddl, MM)

        do i = 1, nc_total

            ! Compute inverse of Jacobian matrix extended
            invJext = 0.d0
            do j = 1, dimen
                invJext((j-1)*dimen+1:j*dimen, (j-1)*dimen+1:j*dimen) = invJ(:, :, i)
            end do

            ! Compute the coefficients to use in Fint
            MM_S = matmul(transpose(MM), sigma(:, i))
            coef_Fint(:, i) = matmul(transpose(invJext), MM_S) * detJ(i)

            ! Compute the coefficients to use in Stiffness matrix
            MM_dSdE = matmul(transpose(MM), Dalg(:, :, i))
            MDM = matmul(MM_dSdE, MM)
            coef_temp = matmul(transpose(invJext), MDM)
            coef_Stiff(:, :, i) = matmul(coef_temp, invJext) * detJ(i)
        
        end do

    end subroutine compute_coefficients

    ! =========================
    ! PERFECT PLASTICITY
    ! =========================

    subroutine yield_perfect_plasticity(sigma_Y, sigma, f, norm, N)
        !! Computes the value of f (consistency condition) in perfect plasticity criteria
        
        implicit none
        ! Input / output
        ! ---------------
        double precision, intent(in) ::  sigma_Y, sigma
        dimension :: sigma(ddl)
        double precision, intent(out) :: f, norm, N
        dimension :: N(ddl)

        ! Local data
        ! ----------------
        double precision :: eta
        dimension :: eta(ddl)

        ! Compute deviatoric of sigma tensor: nu_trial
        call compute_deviatoric(dimen, ddl, sigma, eta)

        ! Compute the norm sqrt(nu_trial : nu_trial)
        call stdcst(dimen, ddl, eta, eta, norm)
        norm = sqrt(norm)
        f = norm - sqrt(2.d0/3.d0) * sigma_Y   

        ! Compute gradient of the gradient of f
        if (norm.gt.0.d0) then
            N = 1.d0/norm*eta
        else
            N = 0.d0
        end if

    end subroutine yield_perfect_plasticity

    subroutine cpp_perfect_plasticity(mat, deps, ep_n0, sigma_n0, ep_n1, sigma_n1, Dalg)
        !! Return closest point proyection (cpp) in perfect plasticity criteria

        implicit none
        ! Input / output
        ! ---------------
        type(material), pointer :: mat
        double precision, intent(in) :: deps, ep_n0, sigma_n0
        dimension :: deps(ddl), sigma_n0(ddl)

        double precision, intent(out) :: ep_n1, sigma_n1, Dalg
        dimension :: sigma_n1(ddl), Dalg(ddl, ddl)
        
        ! Local data
        ! ------------
        double precision :: f, norm
        double precision :: sigma_trial, N, NNT
        dimension :: sigma_trial(ddl), N(ddl), NNT(ddl, ddl)
                       
        ! Compute elastic predictor
        sigma_trial = sigma_n0 + matmul(mat%Ctensor, deps)
        
        ! Review condition
        call yield_perfect_plasticity(mat%sigma_Y, sigma_trial, f, norm, N)

        if (f.le.0.d0) then ! Elastic point

            ! Send back values
            sigma_n1 = sigma_trial
            ep_n1 = ep_n0
            Dalg = mat%Ctensor

        else ! Plastic point

            ! Update stress
            sigma_n1 = sigma_trial - f*N

            ! Update effective plastic strain 
            ep_n1 = ep_n0 + f/(mat%mu*sqrt(6.d0))

            ! Compute consistent tangent matrix
            call stkronst(ddl, N, N, NNT)
            Dalg = mat%Ctensor - 2*mat%mu*((1-f/norm)*NNT + f/norm*mat%Idev)

        end if

    end subroutine cpp_perfect_plasticity

    ! =========================
    ! COMBINED HARDENING
    ! =========================

    subroutine yield_combined_hardening(sigma_Y, H, beta, sigma, alpha, ep, f, norm, N)
        !! Computes the value of f (consistency condition) in combined hardening criteria
        
        implicit none
        ! Input / output
        ! ---------------
        double precision, intent(in) ::  sigma_Y, beta, H, sigma, alpha, ep
        dimension :: sigma(ddl), alpha(ddl)
        double precision, intent(out) :: f, norm, N
        dimension :: N(ddl)

        ! Local data
        ! ----------------
        double precision :: eta
        dimension :: eta(ddl)

        ! Compute deviatoric of sigma tensor
        call compute_deviatoric(dimen, ddl, sigma, eta)

        ! Compute eta
        eta = eta - alpha

        ! Compute the norm sqrt(nu_trial : nu_trial)
        call stdcst(dimen, ddl, eta, eta, norm)
        norm = sqrt(norm)
        f = norm - sqrt(2.d0/3.d0) * (sigma_Y + (1-beta)*H*ep)   

        ! Compute gradient of the gradient of f
        if (norm.gt.0.d0) then
            N = 1.d0/norm*eta
        else
            N = 0.d0
        end if

    end subroutine yield_combined_hardening

    subroutine cpp_combined_hardening(mat, deps, alpha_n0, ep_n0, sigma_n0, &
                                    alpha_n1, ep_n1, sigma_n1, Dalg)
        !! Return closest point proyection (cpp) in perfect plasticity criteria

        implicit none
        ! Input / output
        ! ---------------
        type(material), pointer :: mat
        double precision, intent(in) :: deps, alpha_n0, ep_n0, sigma_n0
        dimension :: alpha_n0(ddl), deps(ddl), sigma_n0(ddl)

        double precision, intent(out) ::alpha_n1, ep_n1, sigma_n1, Dalg
        dimension :: alpha_n1(ddl), sigma_n1(ddl), Dalg(ddl, ddl)
        
        ! Local data
        ! ------------
        double precision :: f, norm, dgamma, c1, c2
        double precision :: sigma_trial, N, NNT
        dimension :: sigma_trial(ddl), N(ddl), NNT(ddl, ddl)
                       
        ! Compute elastic predictor
        sigma_trial = sigma_n0 + matmul(mat%Ctensor, deps)
        
        ! Review condition
        call yield_combined_hardening(mat%sigma_Y, mat%hardening, mat%beta, sigma_trial, alpha_n0, ep_n0, f, norm, N)

        if (f.le.0.d0) then ! Elastic point

            ! Send back values
            sigma_n1 = sigma_trial
            ep_n1 = ep_n0
            alpha_n1 = alpha_n0
            Dalg = mat%Ctensor

        else ! Plastic point

            ! Consistency parameter
            dgamma = f/(2.d0*mat%mu+2.d0*mat%hardening/3.d0)

            ! Update stress
            sigma_n1 = sigma_trial - 2.d0*mat%mu*dgamma*N

            ! Update back stress
            alpha_n1 = alpha_n0 + 2.d0/3.d0*mat%beta*mat%hardening*dgamma*N

            ! Update effective plastic strain 
            ep_n1 = ep_n0 + sqrt(2.d0/3.d0)*dgamma

            ! Compute consistent tangent matrix
            c1 = 4.d0*mat%mu**2.d0/(2.d0*mat%mu+2.d0/3.d0*mat%hardening)
            c2 = 4.d0*mat%mu**2.d0*dgamma/norm
            call stkronst(ddl, N, N, NNT)
            Dalg = mat%Ctensor - (c1-c2)*NNT - c2*mat%Idev

        end if

    end subroutine cpp_combined_hardening

end module elastoplasticity
