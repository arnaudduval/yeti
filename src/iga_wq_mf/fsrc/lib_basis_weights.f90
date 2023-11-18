! ==============================================
! This module computes an accurate weighted quadrature rule using 2 methods
! For both methods, the input is the degree p and an open knot-vector of maximum regularity, that is, 
! the test functions space is S^p_r, with r = p-1.
!
! The ouput of the first method is 4 quadrature rules: 
! - For integrals of the form int N(x) f(x) dx
! - For integrals of the form int N(x) f'(x) dx
! - For integrals of the form int N'(x) f(x) dx
! - For integrals of the form int N'(x) f'(x) dx
! For that we will use functions of the target space S^{p-1}_{r-1}. 
!
! The ouput of the second method is 2 quadrature rules: 
! - For integrals of the form int N(x) f(x) dx
! - For integrals of the form int N'(x) f(x) dx
! For that we will use functions of the target space S^p_{r-1}. 
!
! The first method uses half of quadrature points of the second one (when number of elements is large).
! It seems that the second method is more adapted to structural analysis.
! ===============================================

! ------------------------------------------------------
! GLOBAL FUNCTIONS: expected to work with any knotvector
! ------------------------------------------------------

subroutine find_interpolation_kvspan(nnz, array, x, offset, span, threshold)
    !! Finds the interpolation span of the given value x. 
    !! Ex: Given the nodes {0, 0.5, 1} and x = 0.25, the interpolation span is 1

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: nnz, offset
    double precision, intent(in) :: array, x, threshold
    dimension :: array(nnz)

    integer, intent(out) :: span 

    span = 2 + offset

    do while ((span.lt.nnz-offset).and.((array(span)-x).le.threshold))
        span = span + 1
    end do

    span = span - 1 

end subroutine find_interpolation_kvspan

subroutine find_multiplicity(size_kv, knotvector, x, multiplicity, threshold)
    !! Finds the multiplicity of a given knot.
    !! Ex: Given the knot-vector {0, 0, 0, 0.5, 0.5, 1, 1, 1} and x = 0.5, the multiplicity is 2.

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_kv
    double precision, intent(in) :: knotvector, x, threshold
    dimension :: knotvector(size_kv)

    integer, intent(out) :: multiplicity 

    ! Local data
    ! ----------
    integer :: i

    multiplicity = 0

    do i = 1, size_kv
        if (abs(x-knotvector(i)).le.threshold) then 
            multiplicity = multiplicity + 1
        end if
    end do

end subroutine find_multiplicity

subroutine increase_multiplicity(repeat, degree, size_kv_in, kv_in, size_kv_out, kv_out, threshold)
    !! Computes a new knot-vector using p-refinement 
    !! Ex: Given the knot-vector = [0, 0, 0, 0.5, 1, 1, 1] and repeat = 1, then new knot vector = [0, 0, 0, 0.5, 0.5, 1, 1, 1]

    implicit none
    ! Input / output data
    ! -----------------
    integer, intent(in) :: repeat, degree, size_kv_in
    double precision, intent(in) :: kv_in, threshold
    dimension :: kv_in(size_kv_in)

    integer, intent(inout) :: size_kv_out
    double precision, intent(out) :: kv_out
    dimension :: kv_out(*)

    ! Local data
    ! -----------
    integer :: i, j, c, size_nodes, multiplicity
    double precision :: nodes
    dimension :: nodes(size_kv_in+1)

    call find_unique_array(size_kv_in, kv_in, nodes)
    size_nodes = int(nodes(size_kv_in+1))

    if (size_kv_out.le.0) then 
        size_kv_out = size_kv_in + (size_nodes - 2)*repeat

    else 
        c = 1
        do i = 1, degree + 1
            kv_out(c) = 0.d0
            c = c + 1
        end do

        do i = 2, size_nodes - 1
            call find_multiplicity(size_kv_in, kv_in, nodes(i), multiplicity, threshold)
            multiplicity = multiplicity + repeat

            do j = 1, multiplicity
                kv_out(c) = nodes(i)
                c = c + 1
            end do
        end do

        do i = 1, degree+1
            kv_out(c) = 1.d0
            c = c + 1
        end do

    end if

end subroutine increase_multiplicity

subroutine get_basis_coo(degree, size_ukv, ukv, size_kv, knotvector, nb_knots, knots, basis, indices, threshold)
    !! Finds the basis B0 and B1 for every given knot. 
    !! The algorithm computes by itself the knot-vector span of a given knot and 
    !! returns the value of the basis B0 and B1 for that knot. 
    !! Knowing that on each span, there are (degree+1) functions, it returns each time (degree+1) values.

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_ukv, size_kv, nb_knots
    double precision, intent(in) :: ukv, knotvector, knots, threshold
    dimension :: ukv(size_ukv), knotvector(size_kv), knots(nb_knots)

    integer, intent(out) :: indices
    dimension :: indices((degree+1)*nb_knots, 2)
    double precision, intent(out) :: basis
    dimension :: basis((degree+1)*nb_knots, 2)

    ! Local data
    ! ----------
    integer :: i, j, k, nbctrlpts, nbel, multiplicity
    integer :: functions_span, span
    dimension :: functions_span(degree+1), span(2)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: B0t, B1t
    dimension :: B0t(degree+1), B1t(degree+1)

    nbctrlpts = size_kv - degree - 1
    nbel      = size_ukv - 1

    ! Set table of functions 
    allocate(table_functions_span(nbel, degree+1))
    table_functions_span = 0
    do j = 1, degree+1
        table_functions_span(1, j) = j 
    end do

    do i = 2, size_ukv-1
        call find_multiplicity(size_kv, knotvector, ukv(i), multiplicity, threshold)
        table_functions_span(i, 1) = table_functions_span(i-1, 1) + multiplicity
        do j = 2, degree+1
            table_functions_span(i, j) = table_functions_span(i, 1) + j - 1
        end do
    end do

    do i = 1, nb_knots

        ! Computes B0 and B1 using a YETI function
        call find_interpolation_kvspan(size_kv, knotvector, knots(i), degree, span(1), threshold)
        call find_interpolation_kvspan(size_ukv, ukv, knots(i), 0, span(2), threshold)
        functions_span = table_functions_span(span(2), :)
        call dersbasisfuns(span(1), degree, nbctrlpts, knots(i), knotvector, B0t, B1t)

        ! Save in COO format
        do j = 1, degree+1
            k = (i - 1)*(degree + 1) + j
            basis(k, :)   = [B0t(j), B1t(j)]
            indices(k, :) = [functions_span(j), i]                                
        end do
    end do

end subroutine get_basis_coo

subroutine create_uniformmaxregular_knotvector(degree, nbel, knotvector)
    !! Gets an open uniform with maximum regularity knot-vector 

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in):: degree, nbel

    double precision, intent(out) :: knotvector 
    dimension :: knotvector(nbel+2*degree+1)

    ! Local data
    ! ----------
    integer :: i, c

    ! Create knotvector 
    knotvector = 0.d0

    c = 1
    do i = 1, degree+1
        knotvector(c) = 0.d0
        c = c + 1
    end do

    do i = 2, nbel
        knotvector(c) = dble(i - 1)/dble(nbel) 
        c = c + 1
    end do

    do i = 1, degree+1
        knotvector(c) = 1.d0
        c = c + 1
    end do
        
end subroutine create_uniformmaxregular_knotvector

! ---------------------------------------------------------------
! WQ FUNCTIONS: expected to work in weighted-quadrature approach
! ---------------------------------------------------------------

subroutine wq_solve_weights(nr_obj, Bshape_obj, nr_test, nc, BBtest, II, IIshape, weights_obj)
    !! Gets the quadrature rules of objective (or target) functions 
    !! Here the integral II = int Btest Bobj from 0 to 1
    
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_test, nr_obj, nc
    integer, intent(in) :: Bshape_obj, IIshape
    dimension :: Bshape_obj(nr_obj, 2), IIshape(nr_test, nr_obj)
    double precision, intent(in) ::  BBtest, II
    dimension :: BBtest(nr_test, nc), II(nr_test, nr_obj)

    double precision, intent(out) :: weights_obj
    dimension :: weights_obj(nr_obj, nc)

    ! Local data
    ! ----------
    integer :: i, j, Pmin, Pmax, Fmin, Fmax

    weights_obj = 0.d0  

    do i = 1, nr_obj

        Pmin = Bshape_obj(i, 1)
        Pmax = Bshape_obj(i, 2)

        Fmin = 1
        do j = 1, nr_test
            if (IIshape(j, i).gt.0) then
                Fmin = j
                exit 
            end if
        end do 

        Fmax = nr_test
        do j = nr_test, 1, -1
            if (IIshape(j, i).gt.0) then
                Fmax = j
                exit 
            end if
        end do 
        
        call solve_linear_system(Fmax-Fmin+1, Pmax-Pmin+1, BBtest(Fmin:Fmax, Pmin:Pmax), &
                                II(Fmin:Fmax, i), weights_obj(i, Pmin:Pmax))
    end do
end subroutine wq_solve_weights

module quadrature_rules

    type :: genquadrature

        integer :: degree
        double precision, dimension(:), pointer :: knotvector=>null()
        double precision :: span_threshold = 1.d-8, threshold = 1.d-12
        integer :: nbctrlpts, nbel, size_ukv, size_kv
        double precision, dimension(:), allocatable :: ukv

    end type genquadrature

    type :: gaussquadrature

        type(genquadrature) :: genquad
        integer :: order = 0
        integer :: nbqp
        double precision, dimension(:), allocatable :: quadptspos, parametricweights
        double precision, dimension(:), allocatable :: isopositions, isoweights
        
    end type gaussquadrature

    type ::  weightedquadrature

        type(genquadrature) :: genquad
        double precision, dimension(:), pointer :: quadptspos=>null()
        integer :: method = 1
        integer, dimension(:, :), pointer :: B0shape=>null(), B1shape=>null()
        integer :: nbqp

    end type weightedquadrature

contains

    subroutine init_genquad(obj, degree, knotvector)

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: degree
        double precision, dimension(:), target, intent(in) :: knotvector
        type(genquadrature) :: obj

        ! local data
        ! ----------
        double precision, allocatable, dimension(:) :: nodes

        obj%degree   = degree
        obj%knotvector => knotvector
        obj%size_kv   = size(knotvector) 
        obj%nbctrlpts = obj%size_kv - obj%degree - 1

        allocate(nodes(obj%size_kv+1))
        call find_unique_array(obj%size_kv, obj%knotvector, nodes)
        obj%size_ukv  = int(nodes(obj%size_kv+1))
        obj%nbel      = obj%size_ukv - 1
        allocate(obj%ukv(obj%size_ukv))
        obj%ukv = nodes(:obj%size_ukv)

    end subroutine init_genquad

    subroutine get_basis_simplified_coo(obj, nb_knots, knots, basis, indices)

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nb_knots
        double precision, intent(in) :: knots
        dimension :: knots(nb_knots)
        type(genquadrature) :: obj

        double precision, intent(out) :: basis
        dimension :: basis((obj%degree+1)*nb_knots, 2)

        integer, intent(out) :: indices
        dimension :: indices((obj%degree+1)*nb_knots, 2)

        call get_basis_coo(obj%degree, obj%size_ukv, obj%ukv, obj%size_kv, obj%knotvector, &
                            nb_knots, knots, basis, indices, obj%span_threshold)

    end subroutine get_basis_simplified_coo

    subroutine get_basis_simplified_dense(obj, nb_knots, knots, B0, B1)

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nb_knots
        double precision, intent(in) :: knots
        dimension :: knots(nb_knots)
        type(genquadrature) :: obj
        double precision, intent(out) :: B0, B1
        dimension :: B0(obj%nbctrlpts, nb_knots), B1(obj%nbctrlpts, nb_knots)

        ! Local data
        ! ----------
        double precision :: basis
        dimension :: basis((obj%degree+1)*nb_knots, 2)
        integer :: indices
        dimension :: indices((obj%degree+1)*nb_knots, 2)

        call get_basis_coo(obj%degree, obj%size_ukv, obj%ukv, obj%size_kv, obj%knotvector, &
                            nb_knots, knots, basis, indices, obj%span_threshold)

        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 1), size(B0, 1), size(B0, 2), B0)
        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 2), size(B1, 1), size(B1, 2), B1)

    end subroutine get_basis_simplified_dense

    ! ------------

    subroutine init_gaussquad(obj, degree, knotvector)
        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: degree
        double precision, dimension(:), target, intent(in) :: knotvector
        type(gaussquadrature) :: obj

        call init_genquad(obj%genquad, degree, knotvector)

    end subroutine init_gaussquad

    subroutine getisoinfo_gaussquad(obj)
        implicit none
        ! Input / output data
        ! -------------------
        type(gaussquadrature) :: obj
        
        ! Local data
        ! ----------
        double precision, allocatable, dimension(:, :) :: GaussPdsCoord

        if (obj%order.eq.0) obj%order = obj%genquad%degree + 1
        if (.not.allocated(obj%isopositions)) allocate(obj%isopositions((obj%order)))
        if (.not.allocated(obj%isoweights))   allocate(obj%isoweights((obj%order)))
        allocate(GaussPdsCoord(2, obj%order))
        
        call Gauss(obj%order, 1, GaussPdsCoord, 0)
        obj%isoweights   = GaussPdsCoord(1, :)
        obj%isopositions = GaussPdsCoord(2, :)

    end subroutine getisoinfo_gaussquad

    subroutine findquadpos_gaussquad(obj)
        implicit none
        ! Input / output data
        ! -------------------
        type(gaussquadrature) :: obj

        ! Local data
        ! ----------
        integer :: i, j, k, nbel
        double precision, dimension(:), allocatable :: nodes

        allocate(nodes(obj%genquad%size_ukv))
        nbel     = obj%genquad%nbel
        nodes    = obj%genquad%ukv 
        obj%nbqp = nbel*obj%order

        if (.not.allocated(obj%quadptspos)) allocate(obj%quadptspos(obj%nbqp))
        do i = 1, nbel
            do j = 1, obj%order
                k = (i - 1)*obj%order + j
                obj%quadptspos(k) = 0.5d0*(obj%isopositions(j)/dble(nbel) + nodes(i+1) + nodes(i))
            end do
        end do
        
    end subroutine findquadpos_gaussquad

    subroutine findparametricweights_gaussquad(obj)
        implicit none
        ! Input / output data
        ! -------------------
        type(gaussquadrature) :: obj

        ! Local data
        ! ----------
        integer :: i, j, k, nbel

        nbel     = obj%genquad%nbel
        obj%nbqp = nbel*obj%order

        if (.not.allocated(obj%parametricweights)) allocate(obj%parametricweights(obj%nbqp))
        do i = 1, nbel
            do j = 1, obj%order
                k = (i - 1)*obj%order + j
                obj%parametricweights(k) = 0.5d0/dble(nbel)*obj%isoweights(j)
            end do
        end do
        
    end subroutine findparametricweights_gaussquad

    ! -----------

    subroutine init_wq(obj, degree, knotvector, quadptspos, B0shape, B1shape, method)
        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: degree, method
        double precision, dimension(:), target, intent(in) :: knotvector, quadptspos
        integer, dimension(:, :), target, intent(in) :: B0shape, B1shape
        type(weightedquadrature) :: obj

        ! Local data
        ! ---------- 
        call init_genquad(obj%genquad, degree, knotvector)
        if ((method.ge.1).and.(method.le.2)) obj%method = method
        obj%nbqp       = size(quadptspos)
        obj%quadptspos => quadptspos
        obj%B0shape    => B0shape
        obj%B1shape    => B1shape
        
    end subroutine init_wq

    subroutine findbasisweightrules_md1_wq(obj, basis, weights, indices)
        !! Returns the basis and weights data at the quadrature points in WQ approach 
        !! The ouput of this first method is 4 quadrature rules: 
        !! - For integrals of the form int N(x) f(x) dx
        !! - For integrals of the form int N(x) f'(x) dx
        !! - For integrals of the form int N'(x) f(x) dx
        !! - For integrals of the form int N'(x) f'(x) dx
        !! For that we will use functions of the target space S^{p-1}_{r-1}. 
    
        implicit none 
        ! Input / output data
        ! -------------------
        type(weightedquadrature) :: obj
        double precision, intent(out) :: basis, weights
        dimension :: basis((obj%genquad%degree+1)*obj%nbqp, 2), weights((obj%genquad%degree+1)*obj%nbqp, 4)

        integer, intent(out) :: indices
        dimension :: indices((obj%genquad%degree+1)*obj%nbqp, 2)

        ! Local data
        ! ----------      
        integer :: i, j, c
        double precision :: B0wq_p0, B1wq_p0, W00, W01, W10, W11 
        dimension ::    B0wq_p0(obj%genquad%nbctrlpts,  obj%nbqp), B1wq_p0(obj%genquad%nbctrlpts,  obj%nbqp), &
                        W00(obj%genquad%nbctrlpts, obj%nbqp), W01(obj%genquad%nbctrlpts, obj%nbqp), &
                        W10(obj%genquad%nbctrlpts, obj%nbqp), W11(obj%genquad%nbctrlpts, obj%nbqp)
        
        ! For space S^p_r
        type(gaussquadrature) :: gauss_p0
        double precision, allocatable, dimension(:, :) :: B0cgg_p0, B1cgg_p0    
        integer, allocatable, dimension(:, :) :: Bcgg_p0_int
    
        ! For space S^{p-1}_{r-1}
        type(gaussquadrature) :: gauss_p1
        integer :: size_kv_p1, degree_p1
        double precision, allocatable, dimension(:) :: knotvector_p1
        double precision, allocatable, dimension(:, :) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
        integer, dimension(:, :), allocatable :: Bcgg_p1_int
    
        ! Integrals and weights
        double precision, allocatable, dimension(:, :) :: II
        integer, allocatable, dimension(:, :) :: IIshape

        ! ------------
        ! Space S^p_r
        ! ------------   
        ! Find positions and weights in IGA approach
        call init_gaussquad(gauss_p0, obj%genquad%degree, obj%genquad%knotvector)
        call getisoinfo_gaussquad(gauss_p0)
        call findquadpos_gaussquad(gauss_p0)
        call findparametricweights_gaussquad(gauss_p0)

        ! Find basis at Gauss quadrature points
        allocate(B0cgg_p0(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp), &
                B1cgg_p0(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp))
        call get_basis_simplified_dense(gauss_p0%genquad, gauss_p0%nbqp, gauss_p0%quadptspos, B0cgg_p0, B1cgg_p0)

        ! Find basis at WQ quadrature points
        call get_basis_simplified_coo(gauss_p0%genquad, obj%nbqp, obj%quadptspos, basis, indices) 
        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 1), size(B0wq_p0, 1), size(B0wq_p0, 2), B0wq_p0)
        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 2), size(B1wq_p0, 1), size(B1wq_p0, 2), B1wq_p0)

        ! -----------------------
        ! For space S^{p-1}_{r-1}
        ! -----------------------
        ! Set properties of new space
        degree_p1  = obj%genquad%degree - 1
        size_kv_p1 = obj%genquad%size_kv - 2
        allocate(knotvector_p1(size_kv_p1))
        knotvector_p1 = obj%genquad%knotvector(2:obj%genquad%size_kv-1)
        call init_gaussquad(gauss_p1, degree_p1, knotvector_p1)

        ! Find basis function values at Gauss points
        allocate(B0cgg_p1(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp), B1cgg_p1(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp))
        call get_basis_simplified_dense(gauss_p1%genquad, gauss_p0%nbqp, gauss_p0%quadptspos, B0cgg_p1, B1cgg_p1) 
        deallocate(B1cgg_p1)
    
        ! Find basis function values at WQ points
        allocate(B0wq_p1(gauss_p1%genquad%nbctrlpts, obj%nbqp), B1wq_p1(gauss_p1%genquad%nbctrlpts, obj%nbqp))
        call get_basis_simplified_dense(gauss_p1%genquad, obj%nbqp, obj%quadptspos, B0wq_p1, B1wq_p1) 
        deallocate(B1wq_p1)

        ! ---------------------
        ! Integrals and Weights
        ! ---------------------
        allocate(Bcgg_p0_int(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp))
        allocate(Bcgg_p1_int(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp))
        Bcgg_p0_int = 0; Bcgg_p1_int = 0
        where (abs(B0cgg_p0).gt.obj%genquad%threshold) Bcgg_p0_int = 1
        where (abs(B0cgg_p1).gt.obj%genquad%threshold) Bcgg_p1_int = 1
    
        ! --------------
        allocate(IIshape(gauss_p0%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        IIshape = matmul(Bcgg_p0_int, transpose(Bcgg_p0_int))
    
        ! Compute B0_p0 * W * B0_p0'
        allocate(II(gauss_p0%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        call gemm_AWB(1, size(B0cgg_p0, 1), size(B0cgg_p0, 2), B0cgg_p0, &
                    size(B0cgg_p0, 1), size(B0cgg_p0, 2), B0cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II)

        ! Compute W00
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B0shape, gauss_p0%genquad%nbctrlpts, &
                            obj%nbqp, B0wq_p0, II, IIshape, W00)
    
        ! Compute  B0_p0 * W * B1_p0'
        call gemm_AWB(1, size(B0cgg_p0, 1), size(B0cgg_p0, 2), B0cgg_p0, &
                    size(B1cgg_p0, 1), size(B1cgg_p0, 2), B1cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II)
    
        ! Compute W10
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B1shape, gauss_p0%genquad%nbctrlpts, &
                            obj%nbqp, B0wq_p0, II, IIshape, W10)
        deallocate(IIshape, II)
    
        ! --------------
        allocate(IIshape(gauss_p1%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        IIshape = matmul(Bcgg_p1_int, transpose(Bcgg_p0_int))
    
        ! Compute B0_p1 * W * B0_p0'
        allocate(II(gauss_p1%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        call gemm_AWB(1, size(B0cgg_p1, 1), size(B0cgg_p1, 2), B0cgg_p1, &
                    size(B0cgg_p0, 1), size(B0cgg_p0, 2), B0cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II)
    
        ! Compute W01
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B0shape, gauss_p1%genquad%nbctrlpts, &
                            obj%nbqp, B0wq_p1, II, IIshape, W01)
        
        ! Compute B0_p1 * W * B1_p0'
        call gemm_AWB(1, size(B0cgg_p1, 1), size(B0cgg_p1, 2), B0cgg_p1, &
                    size(B1cgg_p0, 1), size(B1cgg_p0, 2), B1cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II)

        ! Compute W11
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B1shape, gauss_p1%genquad%nbctrlpts, &
                            obj%nbqp, B0wq_p1, II, IIshape, W11) 
        deallocate(IIshape, II)

        weights = 0.d0
        do c = 1, size(indices, 1)
            i = indices(c, 1)
            j = indices(c, 2)
            if ((i.gt.0).and.(j.gt.0)) weights(c, :) = [W00(i, j), W01(i, j), W10(i, j), W11(i, j)]
        end do

    end subroutine findbasisweightrules_md1_wq

    subroutine findbasisweightrules_md2_wq(obj, basis, weights, indices)
        !! Returns the basis and weights data at the quadrature points in WQ approach 
        !! The ouput of the second method is 2 quadrature rules: 
        !! - For integrals of the form int N(x) f(x) dx
        !! - For integrals of the form int N'(x) f(x) dx
        !! For that we will use functions of the target space S^p_{r-1}. 
    
        implicit none 
        ! Input / output data
        ! -------------------
        type(weightedquadrature) :: obj
        double precision, intent(out) :: basis, weights
        dimension :: basis((obj%genquad%degree+1)*obj%nbqp, 2), weights((obj%genquad%degree+1)*obj%nbqp, 4)

        integer, intent(out) :: indices
        dimension :: indices((obj%genquad%degree+1)*obj%nbqp, 2)

        ! Local data
        ! ----------
        integer :: i, j, c
        double precision :: B0wq_p0, B1wq_p0, W00, W11 
        dimension ::    B0wq_p0(obj%genquad%nbctrlpts,  obj%nbqp), B1wq_p0(obj%genquad%nbctrlpts,  obj%nbqp), &
                        W00(obj%genquad%nbctrlpts, obj%nbqp), W11(obj%genquad%nbctrlpts, obj%nbqp)

        ! For space S^p_r
        type(gaussquadrature) :: gauss_p0
        double precision, allocatable, dimension(:, :) :: B0cgg_p0, B1cgg_p0    
        integer, allocatable, dimension(:, :) :: Bcgg_p0_int
    
        ! For space S^p_{r-1}
        type(gaussquadrature) :: gauss_p1
        integer :: size_kv_p1, degree_p1
        double precision, allocatable, dimension(:) :: knotvector_p1
        double precision, allocatable, dimension(:, :) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
        integer, dimension(:, :), allocatable :: Bcgg_p1_int
    
        ! Integrals and weights
        double precision, allocatable, dimension(:, :) :: II
        integer, allocatable, dimension(:, :) :: IIshape
    
        ! -----------
        ! Space S^p_r
        ! -----------   
        ! Find positions and weights in IGA approach
        call init_gaussquad(gauss_p0, obj%genquad%degree, obj%genquad%knotvector)
        call getisoinfo_gaussquad(gauss_p0)
        call findquadpos_gaussquad(gauss_p0)
        call findparametricweights_gaussquad(gauss_p0)

        ! Find basis at Gauss quadrature points
        allocate(B0cgg_p0(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp), &
                B1cgg_p0(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp))
        call get_basis_simplified_dense(gauss_p0%genquad, gauss_p0%nbqp, gauss_p0%quadptspos, B0cgg_p0, B1cgg_p0)
    
        ! Find basis at WQ quadrature points
        call get_basis_simplified_coo(gauss_p0%genquad, obj%nbqp, obj%quadptspos, basis, indices) 
        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 1), size(B0wq_p0, 1), size(B0wq_p0, 2), B0wq_p0)
        call coo2dense(size(basis, 1), indices(:, 1), indices(:, 2), basis(:, 2), size(B1wq_p0, 1), size(B1wq_p0, 2), B1wq_p0)

        ! ---------------
        ! Space S^p_{r-1}
        ! ---------------
        ! Set properties of new space
        degree_p1  = obj%genquad%degree 
        size_kv_p1 = -1
        call increase_multiplicity(1, degree_p1, obj%genquad%size_kv, obj%genquad%knotvector, &
                                size_kv_p1, knotvector_p1, obj%genquad%span_threshold)
        allocate(knotvector_p1(size_kv_p1))
        call increase_multiplicity(1, degree_p1, obj%genquad%size_kv, obj%genquad%knotvector, &
                                size_kv_p1, knotvector_p1, obj%genquad%span_threshold)
        call init_gaussquad(gauss_p1, degree_p1, knotvector_p1)

        ! Find basis function values at Gauss points
        allocate(B0cgg_p1(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp), B1cgg_p1(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp))
        call get_basis_simplified_dense(gauss_p1%genquad, gauss_p0%nbqp, gauss_p0%quadptspos, B0cgg_p1, B1cgg_p1) 
        deallocate(B1cgg_p1)
    
        ! Find basis function values at WQ points
        allocate(B0wq_p1(gauss_p1%genquad%nbctrlpts, obj%nbqp), B1wq_p1(gauss_p1%genquad%nbctrlpts, obj%nbqp))
        call get_basis_simplified_dense(gauss_p1%genquad, obj%nbqp, obj%quadptspos, B0wq_p1, B1wq_p1) 
        deallocate(B1wq_p1)

        ! ---------------------
        ! Integrals and Weights
        ! ---------------------
        allocate(Bcgg_p0_int(gauss_p0%genquad%nbctrlpts, gauss_p0%nbqp))
        allocate(Bcgg_p1_int(gauss_p1%genquad%nbctrlpts, gauss_p0%nbqp))
        Bcgg_p0_int = 0; Bcgg_p1_int = 0
        where (abs(B0cgg_p0).gt.obj%genquad%threshold) Bcgg_p0_int = 1
        where (abs(B0cgg_p1).gt.obj%genquad%threshold) Bcgg_p1_int = 1

        allocate(IIshape(gauss_p1%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        IIshape = matmul(Bcgg_p1_int, transpose(Bcgg_p0_int))

        ! Compute B0_p1 * W * B0_p0'
        allocate(II(gauss_p1%genquad%nbctrlpts, gauss_p0%genquad%nbctrlpts))
        call gemm_AWB(1, size(B0cgg_p1, 1),size(B0cgg_p1, 2), B0cgg_p1, &
                    size(B0cgg_p0, 1), size(B0cgg_p0, 2), B0cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II) 

        ! Compute W0
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B0shape, &
                                gauss_p1%genquad%nbctrlpts, obj%nbqp, B0wq_p1, II, IIshape, W00)

        ! Compute = B0_p1 * W * B1_p0'
        call gemm_AWB(1, size(B0cgg_p1, 1),size(B0cgg_p1, 2), B0cgg_p1, &
                    size(B0cgg_p0, 1), size(B0cgg_p0, 2), B1cgg_p0, &
                    gauss_p0%parametricweights, size(II, 1), size(II, 2), II)
                        
        ! Compute W1
        call wq_solve_weights(gauss_p0%genquad%nbctrlpts, obj%B1shape, &
                                gauss_p1%genquad%nbctrlpts, obj%nbqp, B0wq_p1, II, IIshape, W11)   
        
        weights = 0.d0
        do c = 1, size(indices, 1)
            i = indices(c, 1)
            j = indices(c, 2)
            if ((i.gt.0).and.(j.gt.0)) weights(c, :) = [W00(i, j), W00(i, j), W11(i, j), W11(i, j)]
        end do
    
    end subroutine findbasisweightrules_md2_wq

end module quadrature_rules