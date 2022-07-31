! ==========================
! module :: Elasto-plasticity 
! author :: Joaquin Cornejo
! ==========================

subroutine dot_prod_plasticity(dimen, nnz, v1, v2, result)

    implicit none
    ! Input/ output
    ! -----------------
    integer, intent(in) :: dimen, nnz
    double precision, intent(in) :: v1, v2
    dimension :: v1(dimen, nnz), v2(dimen, nnz)

    double precision :: result

    ! Local data
    ! ------------
    integer :: i
    double precision :: rtemp

    ! Initialize
    result = 0.d0

    do i = 1, dimen 
        rtemp = dot_product(v1(i, :), v2(i, :))
        result = result + rtemp
    end do

end subroutine dot_prod_plasticity

subroutine array_maps_tensor(dimen, p, i, j)
    !! Get the second-order tensor index from 1d array index

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
    !! Get the 1d array index from second-order tensor index 

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

subroutine fourth_order_identity(ddl, identity)
    !! Creates a fourth-order identity 

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
    !! Creates a one kron one tensor

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

subroutine stdcst(dimen, ddl, tensor1, tensor2, result)
    !! Returns second-order tensor double contracted with second-order tensor
    
    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor1, tensor2
    dimension :: tensor1(ddl), tensor2(ddl)

    double precision, intent(out) :: result

    ! Local data
    ! ---------------
    integer :: i, j
    double precision :: t1, t2
    dimension :: t1(ddl, ddl), t2(ddl, ddl)

    ! Convert array to tensor
    call array2st(dimen, ddl, tensor1, t1)
    call array2st(dimen, ddl, tensor2, t2)

    ! Compute double contracted result
    result = 0.d0
    do j = 1, dimen
        do i = 1, dimen
            result = result + t1(i, j)*t2(i, j)
        end do
    end do

end subroutine stdcst

subroutine stkronst(ddl, tensor1, tensor2, result)
    !! Returns kron product

    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: ddl
    double precision, intent(in) :: tensor1, tensor2
    dimension :: tensor1(ddl), tensor2(ddl)

    double precision, intent(out) :: result
    dimension :: result(ddl, ddl)

    ! Local data
    ! ---------------
    integer :: i, j

    do j = 1, ddl
        do i = 1, ddl
            result(i, j) = tensor1(i)*tensor2(j)
        end do
    end do

end subroutine stkronst

subroutine update_dirichlet_3d(nc, array, ndu, ndv, ndw, dod_u, dod_v, dod_w)
    !! Update a array using dirichlet condition
    implicit none
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: nc, ndu, ndv, ndw
    double precision, intent(inout) :: array
    dimension :: array(3, nc)

    integer, intent(in) :: dod_u, dod_v, dod_w
    dimension :: dod_u(ndu), dod_v(ndv), dod_w(ndw)

    ! Update array
    array(1, dod_u) = 0.d0 
    array(2, dod_v) = 0.d0 
    array(3, dod_w) = 0.d0 

end subroutine update_dirichlet_3d

subroutine compute_deviatoric(dimen, ddl, tensor, dev_tensor)
    !! Returns deviatoric tensor 

    implicit none
    ! Input / output data
    ! ----------------------
    integer, intent(in) :: dimen, ddl
    double precision, intent(in) :: tensor
    dimension :: tensor(ddl)

    double precision, intent(out) :: dev_tensor
    dimension :: dev_tensor(ddl)

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
    dev_tensor = tensor - 1.d0/3.d0*trace*one 
    
end subroutine compute_deviatoric

subroutine compute_MM(dimen, ddl, MM)

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

end subroutine compute_MM

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
        double precision :: lambda, mu, bulk
        double precision, dimension(:, :), allocatable :: Ctensor, Stensor
    
    end type material

    contains

    subroutine initialize_mat(object, E, nu)

        implicit none
        ! Input / output
        ! ---------------
        double precision, intent(in) :: E, nu
        type(material), pointer :: object

        ! Local data
        ! ----------------
        double precision :: lambda, mu, bulk
        double precision :: identity, onekronone
        dimension :: identity(ddl, ddl), onekronone(ddl, ddl)

        ! Compute constants
        lambda = nu*E/((1+nu)*(1-2*nu))
        mu = E/(2*(1+nu))
        bulk = lambda + 2.d0/3.d0*mu

        ! Create tensors
        call fourth_order_identity(ddl, identity)
        call one_kron_one(dimen, ddl, onekronone)

        ! Save data 
        allocate(object)
        object%young = E
        object%poisson = nu
        object%lambda = lambda
        object%mu = mu
        object%bulk = bulk

        allocate(object%Ctensor(ddl, ddl), object%Stensor(ddl, ddl))
        object%Ctensor = lambda*onekronone + 2*mu*identity
        object%Stensor = 1.d0/(9.d0*bulk)*onekronone &
                        + 1.d0/(2.d0*mu)*(identity - 1.d0/3.d0*onekronone)

    end subroutine initialize_mat

    subroutine compute_coefficients(nc_total, sigma, dSdE, invJ, detJ, coef_fint, coef_s)

        implicit none 
        ! Input / output 
        ! --------------------
        integer, intent(in) :: nc_total
        double precision, intent(in) :: sigma, dSdE, invJ, detJ
        dimension :: sigma(ddl, nc_total), dSdE(ddl, ddl, nc_total), invJ(dimen, dimen, nc_total), detJ(nc_total)

        double precision, intent(out) :: coef_fint, coef_s
        dimension :: coef_fint(dimen*dimen, nc_total), coef_s(dimen*dimen, dimen*dimen, nc_total)

        ! Local data
        ! -------------
        integer :: i, j
        double precision :: MM, invJext, MMT_S, MDM_temp, MDM, coef_s_temp
        dimension ::    MM(ddl, dimen*dimen), invJext(dimen*dimen, dimen*dimen), MMT_S(dimen*dimen), &
                        MDM_temp(dimen*dimen, ddl), MDM(dimen*dimen, dimen*dimen), coef_s_temp(dimen*dimen, dimen*dimen)
                    
        ! Construct MM
        call compute_MM(dimen, ddl, MM)

        do i = 1, nc_total

            ! Compute MM.T dot sigma
            MMT_S = matmul(transpose(MM), sigma(:, i))
    
            ! Compute coef extended
            invJext = 0.d0
            do j = 1, dimen
                invJext((j-1)*dimen+1:j*dimen, (j-1)*dimen+1:j*dimen) = invJ(:, :, i)
            end do
            
            ! Compute transpose(coef_extend) dot MMT_S
            coef_fint(:, i) = matmul(transpose(invJext), MMT_S) * detJ(i)

            ! Evaluate MM.transpose * DD * MM
            MDM_temp = matmul(transpose(MM), dSdE(:, :, i))
            MDM = matmul(MDM_temp, MM)

            ! Compute coefs_s
            coef_s_temp = matmul(transpose(invJext), MDM)
            coef_s(:, :, i) = matmul(coef_s_temp, invJext) * detJ(i)
        
        end do

    end subroutine compute_coefficients

    subroutine interpolate_strain(nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, invJ, disp_ctrlpts, strain_interp)
        !! Computes strain in 3D case (from parametric space to physical space)
        !! IN CSR FORMAT

        use tensor_methods
        implicit none 
        ! Input/ output
        ! --------------------  
        integer, intent(in) :: nr_total, nc_total, nr_u, nr_v, nr_w, nc_u, nc_v, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w
        dimension :: data_B_u(nnz_u, 2), data_B_v(nnz_v, 2), data_B_w(nnz_w, 2)
        double precision, intent(in) :: invJ, disp_ctrlpts
        dimension :: invJ(dimen, dimen, nc_total), disp_ctrlpts(dimen, nr_total)

        double precision, intent(out) :: strain_interp
        dimension :: strain_interp(ddl, nc_total)

        ! Local data
        !-----------------
        integer :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer :: i, j, k, beta
        dimension :: beta(dimen)
        double precision :: MM, invJext, result, temp
        dimension :: MM(ddl, dimen*dimen), invJext(dimen*dimen, dimen*dimen), result(dimen*dimen, nc_total), temp(dimen*dimen)

        ! Initialize
        call csr2csc(2, nr_u, nc_u, nnz_u, data_B_u, indj_u, indi_u, data_BT_u, indj_T_u, indi_T_u)
        call csr2csc(2, nr_v, nc_v, nnz_v, data_B_v, indj_v, indi_v, data_BT_v, indj_T_v, indi_T_v)
        call csr2csc(2, nr_w, nc_w, nnz_w, data_B_w, indj_w, indi_w, data_BT_w, indj_T_w, indi_T_w)

        ! Construct MM
        call compute_MM(dimen, ddl, MM)

        ! Compute displacement in physical space
        do j = 1, dimen
            do i = 1, dimen
                k = i + (j-1)*dimen
                beta = 1; beta(i) = 2
                call tensor3d_dot_vector_sp(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), &
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), &
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), &
                                disp_ctrlpts(j, :), result(k, :))
            end do
        end do

        do i = 1, nc_total

            ! Compute coef extended
            invJext = 0.d0
            do j = 1, dimen
                invJext((j-1)*dimen+1:j*dimen, (j-1)*dimen+1:j*dimen) = invJ(:, :, i)
            end do
            
            ! Compute invJext dot result
            temp = matmul(invJext, result(:, i)) 

            ! Evaluate MM dot temp
            strain_interp(:, i) = matmul(MM, temp)

        end do

    end subroutine interpolate_strain

    ! =========================
    ! PERFECT PLASTICITY
    ! =========================

    subroutine condition_perfplasticity(sigma_Y, sigma, f, grad_f, grad2_f)
        !! Returns the value of f , that is the condition in perfect plasticity
        
        implicit none
        ! Input / output
        ! ---------------
        double precision, intent(in) ::  sigma_Y, sigma
        dimension :: sigma(ddl)
        double precision, intent(out) :: f, grad_f, grad2_f
        dimension :: grad_f(ddl), grad2_f(ddl, ddl)

        ! Local data
        ! ----------------
        double precision ::  norm, dev_sigma, identity, devdev
        dimension ::    dev_sigma(ddl), devdev(ddl, ddl), identity(ddl, ddl)

        ! Compute deviatoric of sigma tensor
        call compute_deviatoric(dimen, ddl, sigma, dev_sigma)

        ! Compute the norm sqrt(sigma : sigma)
        call stdcst(dimen, ddl, dev_sigma, dev_sigma, norm)
        norm = sqrt(norm)
        f = norm - sqrt(2.d0/3.d0) * sigma_Y     

        ! Compute gradient of f
        grad_f = dev_sigma/norm

        ! Compute gradient of the gradient of f
        call fourth_order_identity(ddl, identity)
        call stkronst(ddl, dev_sigma, dev_sigma, devdev)
        grad2_f = 1.d0/norm*identity - 1.d0/(norm**3) * devdev

    end subroutine condition_perfplasticity

    subroutine cpp_perfplasticity(Ctensor, Stensor, sigma_Y, e_n1, ep_n0, ep_n1, sigma_n1, dSdE)
        !! Return closest point proyection (cpp) in perfect plasticity

        implicit none
        ! Input / output
        ! ---------------
        integer, parameter :: nbiter = 10
        double precision, intent(in) :: Ctensor, Stensor, sigma_Y
        dimension :: Ctensor(ddl, ddl), Stensor(ddl, ddl)
        double precision, intent(in) :: e_n1, ep_n0
        dimension :: e_n1(ddl), ep_n0(ddl)

        double precision, intent(out) :: ep_n1, sigma_n1, dSdE
        dimension :: ep_n1(ddl), sigma_n1(ddl), dSdE(ddl, ddl)
        
        ! Local data
        ! ------------
        integer :: iter
        double precision :: f, dgamma, d2gamma, norm, prod1, prod2
        double precision :: diff_e_ep, grad_f, grad2_f, ep_k, r_k
        dimension :: diff_e_ep(ddl),  grad_f(ddl), grad2_f(ddl, ddl), ep_k(ddl), r_k(ddl)
                       
        double precision :: xi, xi_r, xi_grad, xi_diff, diff_grad, d_ep, N, NNT
        dimension ::    xi(ddl, ddl), xi_r(ddl), xi_grad(ddl), & 
                        xi_diff(ddl), d_ep(ddl), diff_grad(ddl), N(ddl), NNT(ddl, ddl)

        ! Compute elastic predictor
        diff_e_ep = e_n1 - ep_n0
        sigma_n1 = matmul(Ctensor, diff_e_ep)
        
        ! Review condition
        call condition_perfplasticity(sigma_Y, sigma_n1, f, grad_f, grad2_f)

        if (f.le.0) then 
            ! Elastic point
            ep_n1 = ep_n0
            dSdE = Ctensor

        else 
            ! Plastic point
            dgamma = 0.d0
            ep_k = ep_n0

            ! Return-mapping iterative algorithm
            do iter = 1, nbiter
                ! Compute residuals
                diff_e_ep = e_n1 - ep_k
                sigma_n1 = matmul(Ctensor, diff_e_ep)
                call condition_perfplasticity(sigma_Y, sigma_n1, f, grad_f, grad2_f)

                r_k = ep_k - ep_n0 - dgamma*grad_f

                ! Check convergence
                call stdcst(dimen, ddl, r_k, r_k, norm)
                norm = sqrt(norm)
                if ((norm.lt.tol1).and.(f.lt.tol2)) then
                    ep_n1 = ep_k
                    exit
                else
                    ! Compute tangent moduli
                    xi = Stensor + dgamma*grad2_f   ! Actually xi is the inverse of this relation
                                                    ! but inverse is expensive. We may solve a linear system
                    ! Compute delta2 gamma
                    call solve_linear_system(ddl, ddl, xi, r_k, xi_r) 
                    call solve_linear_system(ddl, ddl, xi, grad_f, xi_grad)  
                    call stdcst(dimen, ddl, grad_f, xi_r, prod1)!!!!!!!!!!! matmul or stdcst
                    call stdcst(dimen, ddl, grad_f, xi_grad, prod2)
                    d2gamma = (f + prod1)/prod2

                    ! Compute increment
                    diff_grad = d2gamma*grad_f - r_k
                    call solve_linear_system(ddl, ddl, xi, diff_grad, xi_diff) 
                    d_ep = matmul(Stensor, xi_diff)
                    
                    ! Update plastic strains and consistency
                    dgamma = dgamma + d2gamma
                    ep_k = ep_k + d_ep
                end if
            end do

            ! Compute dSdE
            dSdE = Stensor + dgamma*grad2_f 
            call inverse_matrix(ddl, dSdE)
            
            ! Compute N
            xi_grad = matmul(dSdE, grad_f)
            call stdcst(dimen, ddl, grad_f, xi_grad, prod2)
            N = xi_grad/sqrt(prod2)
            call stkronst(ddl, N, N, NNT)

            ! Update dSdE
            dSdE = dSdE - NNT

        end if

    end subroutine cpp_perfplasticity

end module elastoplasticity
