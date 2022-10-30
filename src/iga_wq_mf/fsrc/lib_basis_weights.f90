! ==============================================
! module: Basis and weights
! author: Joaquin Cornejo
! 
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

subroutine find_knotvector_span(degree, size_kv, knotvector, x, span, span_tol)
    !! Finds the knot-vector span of a given knot. 
    !! Ex: Given the knot-vector {0, 0, 0, 0.5, 1, 1, 1} and x = 0.25, the knot-vector span is 3.

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv 
    double precision, intent(in) :: knotvector, x, span_tol
    dimension :: knotvector(size_kv)

    integer, intent(out) :: span 

    span = degree + 2
    
    do while ((span.lt.(size_kv-degree)).and.((knotvector(span)-x).le.span_tol))
        span = span + 1
    end do
    
    span = span - 1 

end subroutine find_knotvector_span

subroutine find_parametric_span(size_kv, nodes, x, span, span_tol)
    !! Finds the parametric span of the given knot. 
    !! Ex: Given the nodes {0, 0.5, 1} and x = 0.25, the parametric span is 1

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_kv
    double precision, intent(in) :: nodes, x, span_tol
    dimension :: nodes(size_kv+1)

    integer, intent(out) :: span 

    span = 2
    
    do while ((span.lt.int(nodes(size_kv+1))).and.((nodes(span)-x).le.span_tol))
        span = span + 1
    end do
    
    span = span - 1 

end subroutine find_parametric_span

subroutine find_multiplicity(size_kv, knotvector, x, multiplicity, span_tol)
    !! Finds the multiplicity of a given knot.
    !! Ex: Given the knot-vector {0, 0, 0, 0.5, 0.5, 1, 1, 1} and x = 0.5, the multiplicity is 2.

    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_kv
    double precision, intent(in) :: knotvector, x, span_tol
    dimension :: knotvector(size_kv)

    integer, intent(out) :: multiplicity 

    ! Local data
    ! ----------
    integer :: i

    multiplicity = 0

    do i = 1, size_kv
        if (abs(x-knotvector(i)).le.span_tol) then 
            multiplicity = multiplicity + 1
        end if
    end do

end subroutine find_multiplicity

subroutine increase_multiplicity(repeat, degree, size_kv_in, kv_in, size_kv_out, kv_out, span_tol)
    !! Computes a new knot-vector using p-refinement 
    !! Ex: Given the knot-vector = [0, 0, 0, 0.5, 1, 1, 1] and repeat = 1, then new knot vector = [0, 0, 0, 0.5, 0.5, 1, 1, 1]

    implicit none
    ! Input / output data
    ! -----------------
    integer, intent(in) :: repeat, degree, size_kv_in
    double precision, intent(in) :: kv_in, span_tol
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
            call find_multiplicity(size_kv_in, kv_in, nodes(i), multiplicity, span_tol)
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

subroutine verify_regularity_uniformity(degree, size_kv, knotvector, ismaxreg, isuniform, span_tol)
    !! Verifies if a knot-vector is uniform and if it has maximum regularity. 
    !! This is important because all WQ methods works only with knot-vectors with rgeularity r = p-1. 
    !! Ex: The given knot-vector {0, 0, 0, 0.5, 1, 1, 1} works for us but {0, 0, 0, 0.5, 0.5, 1, 1, 1} does not. 
    !! Furthermore, if the knot-vector is uniform, it is possible to reduce computation time.
    !! But the latter is not mandatory for our algorithms. 

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv
    double precision, intent(in) :: knotvector, span_tol
    dimension :: knotvector(size_kv)

    logical, intent(out) :: ismaxreg, isuniform
    
    ! Local data
    ! ----------
    integer :: nbel, nb_ctrlpts, i
    double precision :: nodes
    dimension :: nodes(size_kv+1)
    double precision, dimension(:), allocatable :: unique_kv, diffknot 

    call find_unique_array(size_kv, knotvector, nodes)
    nbel = int(nodes(size_kv+1)) - 1
    nb_ctrlpts = size_kv - degree - 1

    ! Verify if knot-vector has max regularity
    ismaxreg = .false.
    if (nbel+degree.eq.nb_ctrlpts) then 
        ismaxreg = .true.
    end if

    ! Verify if knot vector is uniform
    isuniform = .true.
    allocate(unique_kv(nbel+1), diffknot(nbel+1))
    unique_kv = nodes(:nbel+1)
    call diff_array(2, size(unique_kv), unique_kv, diffknot)

    do i = 1, size(unique_kv)
        if (abs(diffknot(i)).gt.span_tol) then
            isuniform = .false.
            exit
        end if
    end do

end subroutine verify_regularity_uniformity

subroutine set_table_functions_spans(degree, size_kv, nodes, knotvector, table, span_tol)
    !! Sets a table of functions on every span given a knot-vector.
    !! Ex. Given the knot-vector {0, 0, 0, 0.5, 0.5, 1, 1, 1} (degree 2), the unique nodes are {0, 0.5, 1}. 
    !! On the first span [0, 0.5], there are the functions N1, N2, N3, and on the second span [0.5, 1], 
    !! there are the functions N3, N4, N5

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv
    double precision, intent(in) :: nodes, knotvector, span_tol
    dimension :: nodes(size_kv+1), knotvector(size_kv)

    integer, intent(out) :: table
    dimension :: table(int(nodes(size_kv+1))-1, degree+1)

    ! Local data
    ! ----------
    integer :: i, j, nbel, multiplicity

    ! Get the number of elements
    nbel = int(nodes(size_kv+1)) - 1

    ! Set table of functions 
    table = 0
    do j = 1, degree+1
        table(1, j) = j 
    end do

    do i = 2, nbel
        call find_multiplicity(size_kv, knotvector, nodes(i), multiplicity, span_tol)
        table(i, 1) = table(i-1, 1) + multiplicity
        do j = 2, degree+1
            table(i, j) = table(i, 1) + j - 1
        end do
    end do

end subroutine set_table_functions_spans

subroutine get_basis(degree, size_kv, nodes, knotvector, nb_knots, knots, B0, B1, span_tol)
    !! Finds the basis B0 and B1 for every given knot. 
    !! The algorithm computes by itself the knot-vector span of a given knot and 
    !! returns the value of the basis B0 and B1 for that knot. 
    !! Knowing that on each span, there are (degree+1) functions, it returns each time (degree+1) values.

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: nodes, knotvector, knots, span_tol
    dimension :: nodes(size_kv+1), knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: B0, B1
    dimension :: B0(size_kv-degree-1, nb_knots), B1(size_kv-degree-1, nb_knots)

    ! Local data
    ! ----------
    integer :: i, j, k, nb_ctrlpts, nbel
    integer :: functions_span, span, indices
    dimension :: functions_span(degree+1), span(2), indices((degree+1)*nb_knots, 2)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: B0t, B1t, data_B0, data_B1
    dimension :: B0t(degree+1), B1t(degree+1), data_B0((degree+1)*nb_knots), data_B1((degree+1)*nb_knots)

    nb_ctrlpts = size_kv - degree - 1
    nbel = int(nodes(size_kv+1)) - 1

    allocate(table_functions_span(nbel, degree+1))
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span, span_tol)

    do i = 1, nb_knots
        ! Computes B0 and B1 using a YETI function
        call find_knotvector_span(degree, size_kv, knotvector, knots(i), span(1), span_tol)
        call find_parametric_span(size_kv, nodes, knots(i), span(2), span_tol)
        functions_span = table_functions_span(span(2), :)

        call dersbasisfuns(span(1), degree, nb_ctrlpts, knots(i), knotvector, B0t, B1t)
        
        ! Save in COO format
        do j = 1, degree+1
            k = (i - 1)*(degree + 1) + j
            data_B0(k) = B0t(j)
            data_B1(k) = B1t(j)
            indices(k, :) = [functions_span(j), i]                                
        end do
    end do

    ! Convert COO format to dense matrix
    call coo2dense(size(data_B0), indices(:, 1), indices(:, 2), data_B0, nb_ctrlpts, nb_knots, B0)
    call coo2dense(size(data_B1), indices(:, 1), indices(:, 2), data_B1, nb_ctrlpts, nb_knots, B1)

end subroutine get_basis

subroutine get_I_csr(nr, nc, nnz_B, indi_B, indj_B, nnz_I, indi_I, indj_I)
    !! Gets the indices i and j of all non-zero values of matrix I = B * B.T 
    !! B and I in CSR format

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in):: nr, nc, nnz_B
    integer, intent(in) :: indi_B, indj_B
    dimension :: indi_B(nr+1), indj_B(nnz_B)
    
    integer, intent(inout) :: nnz_I
    integer, intent(out) :: indi_I, indj_I
    dimension :: indi_I(nr+1), indj_I(*)

    ! Local data
    ! ----------
    double precision :: MB, MI, ones
    dimension :: MB(nr, nc), MI(nr, nr), ones(nnz_B)
    double precision, allocatable, dimension(:) :: data_I

    ! Compute I = B * B.T
    ones = 1.d0 
    call csr2dense(nnz_B, indi_B, indj_B, ones, nr, nc, MB)
    MI = matmul(MB, transpose(MB))
    
    ! Convert dense matrix to CSR format
    if (nnz_I.le.0) then
        call dense2csr(nr, nr, MI, nnz_I, indi_I, indj_I, data_I)
    else
        allocate(data_I(nnz_I))
        call dense2csr(nr, nr, MI, nnz_I, indi_I, indj_I, data_I)
        deallocate(data_I)
    end if

end subroutine get_I_csr

subroutine create_uniform_knotvector(degree, nbel, nodes, knotvector)
    !! Gets an open uniform with maximum regularity knot-vector 

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in):: degree, nbel

    double precision, intent(out) :: nodes, knotvector 
    dimension :: knotvector(nbel+2*degree+1), nodes(nbel+2*degree+2)

    ! Local data
    ! ----------
    integer ::  i, c

    ! Create nodes 
    nodes = 0.d0; nodes(1) = 0.d0; nodes(nbel+1) = 1.d0

    do i = 2, nbel 
        nodes(i) = dble(i - 1)/dble(nbel) 
    end do

    nodes(size(nodes)) = nbel + 1

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
        
end subroutine create_uniform_knotvector

! --------------------------------------------------------
! IGA FUNCTIONS: expected to work in IGA-Galerkin approach
! --------------------------------------------------------

subroutine iga_get_qp_positions_weights(degree, size_kv, nodes, qp_nnz, qp_position, qp_weight)
    !! Gets quadrature points' positions and weights in IGA approach 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, qp_nnz
    double precision, intent(in) :: nodes
    dimension :: nodes(size_kv+1)

    double precision, intent(out) :: qp_position, qp_weight
    dimension :: qp_position(qp_nnz), qp_weight(qp_nnz)

    ! Local data 
    ! ----------
    double precision :: GaussPdsCoord, pg, wg
    dimension :: GaussPdsCoord(2, degree+1), pg(degree+1), wg(degree+1)
    integer :: nbel, i, j, k

    ! Find position and weight in isoparametric space
    call Gauss(size(pg), 1, GaussPdsCoord, 0)
    wg = GaussPdsCoord(1, :)
    pg = GaussPdsCoord(2, :)

    ! From isoparametric to parametric space
    nbel = int(nodes(size_kv+1)) - 1 
    do i = 1, nbel
        do j = 1, size(pg)
            k = (i - 1)*size(pg) + j
            qp_position(k) = 0.5d0*(pg(j)/dble(nbel) + nodes(i+1) + nodes(i))
            qp_weight(k) = 0.5d0/dble(nbel)*wg(j)
        end do
    end do

end subroutine iga_get_qp_positions_weights 

subroutine iga_get_B_shape(degree, size_kv, nodes, knotvector, Bshape, span_tol)
    !! Gets non zeros positions of B0 and B1 in IGA-Galerkin approach 
    
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: degree, size_kv
    double precision :: knotvector, nodes, span_tol
    dimension :: knotvector(size_kv), nodes(size_kv+1)

    integer, intent(out) :: Bshape
    dimension :: Bshape(size_kv-degree-1, 2)

    ! Local data 
    ! ----------
    integer :: nbel, nb_ctrlpts, min_span, max_span, min_knot, max_knot
    integer, allocatable, dimension(:, :) :: table_points_span, table_functions_span, table_spans_function
    integer :: i, j    

    ! Set number of elements and control points
    nbel = int(nodes(size_kv+1)) - 1 
    nb_ctrlpts = size_kv - degree - 1

    allocate(table_points_span(nbel, 2), &
            table_functions_span(nbel, degree+1), &
            table_spans_function(nb_ctrlpts, 2))

    ! Get table of points over the span
    table_points_span = 1
    table_points_span(1, 2) = degree + 1 
    
    do i = 2, nbel 
        table_points_span(i, 1) = table_points_span(i-1, 2) + 1
        table_points_span(i, 2) = table_points_span(i, 1) + degree
    end do
    
    ! Get table of functions on span 
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span, span_tol)

    ! Get table of spans for each function
    do i = 1, nb_ctrlpts
        min_span = 1
        do j = 1, nbel
            if (any(table_functions_span(j, :).eq.i)) then
                    min_span = j
                    exit 
            end if
        end do 

        max_span = nbel
        do j = nbel, 1, -1
            if (any(table_functions_span(j, :).eq.i)) then
                    max_span = j
                    exit 
            end if
        end do 

        table_spans_function(i, :) = [min_span, max_span]
    end do 
            
    ! Set shape of B0 and B1
    do i = 1, nb_ctrlpts
        min_span = table_spans_function(i, 1)
        max_span = table_spans_function(i, 2)

        min_knot = table_points_span(min_span, 1)
        max_knot = table_points_span(max_span, 2)

        Bshape(i, :) = [min_knot, max_knot]
    end do

end subroutine iga_get_B_shape

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

! ------------------------------------------------------------------------
! MODULES: In order to have a class structure to compute basis and weights
! ------------------------------------------------------------------------
module iga_basis_weights

    implicit none
    double precision, parameter :: span_tol = 1.d-8
    type :: iga
        ! Input
        integer :: degree, size_kv
        double precision, dimension(:), pointer :: knotvector
        
        ! Output
        double precision, dimension(:), pointer :: qp_position, qp_weight, data_B0, data_B1
        integer, dimension(:,:), pointer :: indices
        integer :: nnz_B, nnz_I
        
        ! Local 
        integer :: nb_ctrlpts, nb_qp, size_nodes
        integer, dimension(:,:), pointer :: Bshape
        double precision, dimension(:), pointer :: nodes

    end type iga

contains

    subroutine iga_initialize(obj, degree, size_kv, knotvector)
        !! Initialize IGA-Galerkin approach

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: degree, size_kv
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(iga), pointer :: obj

        ! Set properties
        allocate(obj)
        allocate(obj%knotvector(size_kv))
        obj%degree = degree
        obj%size_kv = size_kv
        obj%nb_ctrlpts = size_kv - degree - 1 
        obj%nnz_I = (2*obj%degree+1)*obj%nb_ctrlpts - obj%degree*(obj%degree+1)
        obj%knotvector = knotvector

        allocate(obj%nodes(size_kv+1))
        call find_unique_array(obj%size_kv, obj%knotvector, obj%nodes)
        obj%size_nodes = int(obj%nodes(size_kv+1))
        obj%nb_qp = (degree + 1)*(obj%size_nodes - 1)
        obj%nnz_B = obj%nb_qp*(obj%degree + 1) 
        
        ! Get quadrature points positions
        allocate(obj%qp_position(obj%nb_qp), obj%qp_weight(obj%nb_qp))
        call iga_get_qp_positions_weights(obj%degree, obj%size_kv, obj%nodes, & 
                                    obj%nb_qp, obj%qp_position, obj%qp_weight)

        ! Set nonzero values of B0 and B1
        allocate(obj%Bshape(obj%nb_ctrlpts, 2))
        call iga_get_B_shape(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, obj%Bshape, span_tol)  

    end subroutine iga_initialize

    subroutine iga_basis_weights_dense2coo(obj)
        !! Computes the basis and the weights in IGA-Galerkin approach
            
        implicit none 
        ! Input / output data
        ! --------------------
        type(iga), pointer :: obj
        
        ! Local data
        ! ----------
        double precision, dimension(:,:), allocatable :: B0, B1
        integer :: i, j, c

        ! Get basis and weights 
        allocate(B0(obj%nb_ctrlpts, obj%nb_qp), B1(obj%nb_ctrlpts, obj%nb_qp))
        call get_basis(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, &
                        obj%nb_qp, obj%qp_position, B0, B1, span_tol)

        ! Save values
        allocate(obj%data_B0(obj%nnz_B), obj%data_B1(obj%nnz_B), obj%indices(obj%nnz_B, 2))
        c = 0
        do i = 1, obj%nb_ctrlpts
            do j = obj%Bshape(i, 1), obj%Bshape(i, 2)
                c = c + 1
                obj%data_B0(c) = B0(i, j)
                obj%data_B1(c) = B1(i, j)
                obj%indices(c, :) = [i, j]
            end do
        end do

    end subroutine iga_basis_weights_dense2coo

end module iga_basis_weights

module wq_basis_weights

    implicit none
    integer, parameter :: r = 2
    double precision, parameter :: tol = 1.d-14, span_tol = 1.d-8

    type :: wq
        ! Input
        integer :: degree, size_kv, method
        double precision, dimension(:), pointer :: knotvector 

        ! Output
        double precision, dimension(:), pointer ::  qp_position, data_B0, data_B1, & 
                                                    data_W00, data_W01, data_W10, data_W11
        integer, dimension(:,:), pointer :: indices
        integer :: nnz_B, nnz_I

        ! Local
        integer :: maxrule, size_nodes, nb_ctrlpts, nb_qp_wq, nb_qp_cgg
        integer, dimension(:, :), pointer :: B0shape, B1shape
        double precision, dimension(:), pointer :: nodes 
        logical :: isuniform
    
    end type wq

contains

    subroutine wq_initialize(obj, degree, size_kv, knotvector, method)
        !! Initialize IGA-WQ approach

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: degree, size_kv, method
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(wq), pointer :: obj

        ! Local data
        ! -----------
        logical :: ismaxreg, isuniform

        ! Verify regularity and uniformity of knot-vector
        call verify_regularity_uniformity(degree, size_kv, knotvector, ismaxreg, isuniform, span_tol)
        if (.not.ismaxreg) stop 'Maximum regularity not respected'

        ! Set properties
        allocate(obj)
        allocate(obj%knotvector(size_kv))
        obj%degree = degree
        obj%size_kv = size_kv
        obj%nb_ctrlpts = size_kv - degree - 1
        obj%nnz_I = (2*obj%degree+1)*obj%nb_ctrlpts - obj%degree*(obj%degree+1) 
        obj%knotvector = knotvector
        obj%isuniform = isuniform

        allocate(obj%nodes(size_kv+1))
        call find_unique_array(obj%size_kv, obj%knotvector, obj%nodes)
        obj%size_nodes = int(obj%nodes(obj%size_kv+1))

        if (method.eq.1) obj%maxrule = 1
        if (method.eq.2) obj%maxrule = 2
        
        obj%method = method
        obj%nb_qp_wq = 2*(obj%degree + r) + (obj%size_nodes - 1)*(obj%maxrule + 1) - 2*obj%maxrule - 3  
        obj%nb_qp_cgg = (obj%degree + 1)*(obj%size_nodes - 1)

        allocate(obj%qp_position(obj%nb_qp_wq), &
                obj%B0shape(obj%nb_ctrlpts, 2), &
                obj%B1shape(obj%nb_ctrlpts, 2))
        obj%qp_position = 0.d0
        obj%B0shape = 0
        obj%B1shape = 0

    end subroutine wq_initialize

    subroutine wq_get_qp_positions(obj)
        !! Gets quadrature points' positions (QP) in WQ approach 
    
        implicit none
        ! Input / output data
        ! -------------------
        type(wq), pointer :: obj
        
        ! Local data
        ! ----------
        integer :: i, j, k
        double precision, allocatable, dimension(:) :: nodes, QPB, QPI
        
        allocate(QPB(obj%degree+r), QPI(2+obj%maxrule))
        allocate(nodes(obj%size_nodes))
        nodes = obj%nodes(1:obj%size_nodes)
    
        ! First span
        call linspace(nodes(1), nodes(2), size(QPB), QPB)
        do j = 1, size(QPB)
            obj%qp_position(j) = QPB(j)
        end do
    
        ! Last span 
        call linspace(nodes(size(nodes)-1), nodes(size(nodes)), size(QPB), QPB)
        do j = 1, size(QPB)
            obj%qp_position(size(obj%qp_position)-size(QPB)+j) = QPB(j)
        end do
        
        ! Inner spans
        if (size(nodes).ge.4) then
            do i = 2, size(nodes)-2
                call linspace(nodes(i), nodes(i+1), size(QPI), QPI)
                do j = 1, size(QPI) 
                    k = size(QPB) + (size(QPI) - 1)*(i - 2) + j - 1    
                    obj%qp_position(k) = QPI(j)
                end do
            end do
        end if
    
    end subroutine wq_get_qp_positions  
    
    subroutine wq_get_B_shape(obj)
        !! Gets non-zero positions of B0 and B1 in WQ approach
    
        implicit none 
        ! Input / output data 
        ! -------------------
        type(wq), pointer :: obj
    
        ! Local data 
        ! ----------
        integer :: degree, nbel, nb_ctrlpts, min_span, max_span, min_knot, max_knot
        integer, allocatable, dimension(:, :) :: table_points_span, table_functions_span, table_spans_function
        integer :: i, j, c
    
        degree = obj%degree
        nbel = obj%size_nodes - 1
        nb_ctrlpts = obj%nb_ctrlpts
    
        allocate(table_points_span(nbel, 2), &
                table_functions_span(nbel, degree+1), &
                table_spans_function(nb_ctrlpts, 2))
    
        ! Get table of points over the span
        table_points_span = 1
        table_points_span(1, 2) = degree + r 
        
        do i = 2, nbel - 1
            table_points_span(i, 1) = table_points_span(i-1, 2)
            table_points_span(i, 2) = table_points_span(i, 1) + 1 + obj%maxrule
        end do
    
        table_points_span(nbel, 1) = table_points_span(nbel-1, 2)
        table_points_span(nbel, 2) = table_points_span(nbel, 1) + degree + r - 1
    
        ! Get table of functions on span 
        call set_table_functions_spans(degree, obj%size_kv, obj%nodes, obj%knotvector, table_functions_span, span_tol)
    
        ! Get table of spans for each function
        do i = 1, nb_ctrlpts
            min_span = 1
            do j = 1, nbel
                if (any(table_functions_span(j, :).eq.i)) then
                        min_span = j
                        exit 
                end if
            end do 
    
            max_span = nbel
            do j = nbel, 1, -1
                if (any(table_functions_span(j, :).eq.i)) then
                        max_span = j
                        exit 
                end if
            end do 
    
            table_spans_function(i, :) = [min_span, max_span]
        end do 
                
        ! Set shape of B0 and B1
        do i = 1, nb_ctrlpts
            min_span = table_spans_function(i, 1)
            max_span = table_spans_function(i, 2)
    
            min_knot = table_points_span(min_span, 1)
            max_knot = table_points_span(max_span, 2)
    
            if (i.eq.1) then 
                max_knot = max_knot - 1
            else if (i.eq.nb_ctrlpts) then
                min_knot = min_knot + 1
            else
                max_knot = max_knot - 1
                min_knot = min_knot + 1
            end if
    
            obj%B0shape(i, :) = [min_knot, max_knot]
        end do
        
        do i = 1, nb_ctrlpts
            min_span = table_spans_function(i, 1)
            max_span = table_spans_function(i, 2)
    
            min_knot = table_points_span(min_span, 1)
            max_knot = table_points_span(max_span, 2)
    
            if ((i.eq.1).or.(i.eq.2)) then 
                max_knot = max_knot - 1
            else if ((i.eq.nb_ctrlpts).or.(i.eq.nb_ctrlpts-1)) then
                min_knot = min_knot + 1
            else
                max_knot = max_knot - 1
                min_knot = min_knot + 1
            end if
    
            obj%B1shape(i, :) = [min_knot, max_knot]
        end do

        ! Get number of non zero values of B
        c = 0
        do i = 1, nb_ctrlpts
            c = c + obj%B1shape(i, 2) - obj%B1shape(i, 1) + 1
        end do
        obj%nnz_B = c
    
    end subroutine wq_get_B_shape

    subroutine wq_basis_weights_method1(obj, B0, B1, W00, W01, W10, W11)
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
        type(wq), pointer :: obj
        double precision, intent(out) :: B0, B1, W00, W01, W10, W11 
        dimension ::    B0(obj%size_kv-obj%degree-1, obj%nb_qp_wq), B1(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
                        W00(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W01(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
                        W10(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W11(obj%size_kv-obj%degree-1, obj%nb_qp_wq)

        ! Local data
        ! ----------       
        ! For space S^p_r
        double precision, allocatable, dimension(:) :: qp_cgg_pos, qp_cgg_weights
        double precision, allocatable, dimension(:,:) :: B0cgg_p0, B1cgg_p0    
        integer, allocatable, dimension(:,:) :: Bcgg_p0_int
    
        ! For space S^{p-1}_{r-1}
        integer :: size_kv_p1, degree_p1, nb_ctrlpts_p1
        double precision, allocatable, dimension(:) :: knotvector_p1, nodes_p1
        double precision, allocatable, dimension(:,:) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
        integer, dimension(:,:), allocatable :: Bcgg_p1_int
    
        ! Integrals and weights
        double precision, allocatable, dimension(:,:) :: II
        integer, allocatable, dimension(:,:) :: IIshape
    
        ! ------------
        ! Space S^p_r
        ! ------------   
        ! Find positions and weights in IGA approach
        allocate(qp_cgg_pos(obj%nb_qp_cgg), qp_cgg_weights(obj%nb_qp_cgg))
        call iga_get_qp_positions_weights(obj%degree, obj%size_kv, obj%nodes, obj%nb_qp_cgg, qp_cgg_pos, qp_cgg_weights)
    
        ! Find basis at Gauss quadrature points
        allocate(B0cgg_p0(obj%nb_ctrlpts, obj%nb_qp_cgg), B1cgg_p0(obj%nb_ctrlpts, obj%nb_qp_cgg))
        call get_basis(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, obj%nb_qp_cgg, &
                        qp_cgg_pos, B0cgg_p0, B1cgg_p0, span_tol)
    
        ! Find basis at WQ quadrature points
        call get_basis(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, obj%nb_qp_wq, obj%qp_position, B0, B1, span_tol) 
    
        ! -----------------------
        ! For space S^{p-1}_{r-1}
        ! -----------------------
        ! Set properties of new space
        degree_p1 = obj%degree - 1
        size_kv_p1 = obj%size_kv - 2
        allocate(knotvector_p1(size_kv_p1), nodes_p1(size_kv_p1+1))
        knotvector_p1 = obj%knotvector(2:obj%size_kv-1)
        call find_unique_array(size_kv_p1, knotvector_p1, nodes_p1)
        nb_ctrlpts_p1 = size_kv_p1 - degree_p1 - 1 
    
        ! Find basis function values at Gauss points
        allocate(B0cgg_p1(nb_ctrlpts_p1, obj%nb_qp_cgg), B1cgg_p1(nb_ctrlpts_p1, obj%nb_qp_cgg))
        call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, obj%nb_qp_cgg, qp_cgg_pos, B0cgg_p1, B1cgg_p1, span_tol) 
        deallocate(B1cgg_p1)
    
        ! Find basis function values at WQ points
        allocate(B0wq_p1(nb_ctrlpts_p1, obj%nb_qp_wq), B1wq_p1(nb_ctrlpts_p1, obj%nb_qp_wq))
        call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, obj%nb_qp_wq, obj%qp_position, B0wq_p1, B1wq_p1, span_tol) 
        deallocate(B1wq_p1)
    
        ! ---------------------
        ! Integrals and Weights
        ! ---------------------
        allocate(Bcgg_p0_int(obj%nb_ctrlpts, obj%nb_qp_cgg))
        allocate(Bcgg_p1_int(nb_ctrlpts_p1, obj%nb_qp_cgg))
        Bcgg_p0_int = 0; Bcgg_p1_int = 0
        where (abs(B0cgg_p0).gt.tol) Bcgg_p0_int = 1
        where (abs(B0cgg_p1).gt.tol) Bcgg_p1_int = 1
    
        ! --------------
        allocate(IIshape(obj%nb_ctrlpts, obj%nb_ctrlpts))
        IIshape = matmul(Bcgg_p0_int, transpose(Bcgg_p0_int))
    
        ! Compute B0_p0 * W * B0_p0'
        allocate(II(obj%nb_ctrlpts, obj%nb_ctrlpts))
        call gemm_AWB(1, obj%nb_ctrlpts, obj%nb_qp_cgg, B0cgg_p0, obj%nb_ctrlpts, obj%nb_qp_cgg, &
                    B0cgg_p0, qp_cgg_weights, obj%nb_ctrlpts, obj%nb_ctrlpts, II)

        ! Compute W00
        call wq_solve_weights(obj%nb_ctrlpts, obj%B0shape, obj%nb_ctrlpts, obj%nb_qp_wq, B0, II, IIshape, W00)
    
        ! Compute  B0_p0 * W * B1_p0'
        call gemm_AWB(1, obj%nb_ctrlpts, obj%nb_qp_cgg, B0cgg_p0, obj%nb_ctrlpts, obj%nb_qp_cgg, & 
                    B1cgg_p0, qp_cgg_weights, obj%nb_ctrlpts, obj%nb_ctrlpts, II)
    
        ! Compute W10
        call wq_solve_weights(obj%nb_ctrlpts, obj%B1shape, obj%nb_ctrlpts, obj%nb_qp_wq, B0, II, IIshape, W10)
        deallocate(IIshape, II)
    
        ! --------------
        allocate(IIshape(nb_ctrlpts_p1, obj%nb_ctrlpts))
        IIshape = matmul(Bcgg_p1_int, transpose(Bcgg_p0_int))
    
        ! Compute B0_p1 * W * B0_p0'
        allocate(II(nb_ctrlpts_p1, obj%nb_ctrlpts))
        call gemm_AWB(1, nb_ctrlpts_p1, obj%nb_qp_cgg, B0cgg_p1, obj%nb_ctrlpts, obj%nb_qp_cgg, & 
                    B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, obj%nb_ctrlpts, II)
    
        ! Compute W01
        call wq_solve_weights(obj%nb_ctrlpts, obj%B0shape, nb_ctrlpts_p1, obj%nb_qp_wq, B0wq_p1, II, IIshape, W01)
        
        ! Compute B0_p1 * W * B1_p0'
        call gemm_AWB(1, nb_ctrlpts_p1, obj%nb_qp_cgg, B0cgg_p1, obj%nb_ctrlpts, obj%nb_qp_cgg, &
                    B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, obj%nb_ctrlpts, II)
                        
        ! Compute W11
        call wq_solve_weights(obj%nb_ctrlpts, obj%B1shape, nb_ctrlpts_p1, obj%nb_qp_wq, B0wq_p1, II, IIshape, W11)               
        deallocate(IIshape, II)
    
    end subroutine wq_basis_weights_method1

    subroutine wq_basis_weights_method2(obj, B0, B1, W00, W01, W10, W11)
        !! Returns the basis and weights data at the quadrature points in WQ approach 
        !! The ouput of the second method is 2 quadrature rules: 
        !! - For integrals of the form int N(x) f(x) dx
        !! - For integrals of the form int N'(x) f(x) dx
        !! For that we will use functions of the target space S^p_{r-1}. 
    
        implicit none 
        ! Input / output data
        ! -------------------
        type(wq), pointer :: obj
        double precision, intent(out) :: B0, B1, W00, W01, W10, W11 
        dimension ::    B0(obj%size_kv-obj%degree-1, obj%nb_qp_wq), B1(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
                        W00(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W01(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
                        W10(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W11(obj%size_kv-obj%degree-1, obj%nb_qp_wq)

        ! Local data
        ! ----------      
        ! For space S^p_r
        double precision, allocatable, dimension(:) :: qp_cgg_pos, qp_cgg_weights
        double precision, allocatable, dimension(:,:) :: B0cgg_p0, B1cgg_p0    
        integer, allocatable, dimension(:,:) :: Bcgg_p0_int
    
        ! For space S^p_{r-1}
        integer :: size_kv_p1, degree_p1, nb_ctrlpts_p1
        double precision, allocatable, dimension(:) :: knotvector_p1, nodes_p1
        double precision, allocatable, dimension(:,:) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
        integer, dimension(:,:), allocatable :: Bcgg_p1_int
    
        ! Integrals and weights
        double precision, allocatable, dimension(:,:) :: II
        integer, allocatable, dimension(:,:) :: IIshape
    
        ! -----------
        ! Space S^p_r
        ! -----------   
        ! Find positions and weights in IGA approach
        allocate(qp_cgg_pos(obj%nb_qp_cgg), qp_cgg_weights(obj%nb_qp_cgg))
        call iga_get_qp_positions_weights(obj%degree, obj%size_kv, obj%nodes, obj%nb_qp_cgg, qp_cgg_pos, qp_cgg_weights)
    
        ! Find basis at Gauss quadrature points
        allocate(B0cgg_p0(obj%nb_ctrlpts, obj%nb_qp_cgg), B1cgg_p0(obj%nb_ctrlpts, obj%nb_qp_cgg))
        call get_basis(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, obj%nb_qp_cgg, &
                        qp_cgg_pos, B0cgg_p0, B1cgg_p0, span_tol)
    
        ! Find basis at WQ quadrature points
        call get_basis(obj%degree, obj%size_kv, obj%nodes, obj%knotvector, obj%nb_qp_wq, obj%qp_position, B0, B1, span_tol) 

        ! ---------------
        ! Space S^p_{r-1}
        ! ---------------
        ! Set properties of new space
        degree_p1 = obj%degree
        size_kv_p1 = -1
        call increase_multiplicity(1, degree_p1, obj%size_kv, obj%knotvector, size_kv_p1, knotvector_p1, span_tol)
        allocate(knotvector_p1(size_kv_p1), nodes_p1(size_kv_p1+1))
        call increase_multiplicity(1, degree_p1, obj%size_kv, obj%knotvector, size_kv_p1, knotvector_p1, span_tol)
        call find_unique_array(size_kv_p1, knotvector_p1, nodes_p1)
        nb_ctrlpts_p1 = size_kv_p1 - degree_p1 - 1

        ! Find basis function values at Gauss points
        allocate(B0cgg_p1(nb_ctrlpts_p1, obj%nb_qp_cgg), B1cgg_p1(nb_ctrlpts_p1, obj%nb_qp_cgg))
        call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, obj%nb_qp_cgg, qp_cgg_pos, B0cgg_p1, B1cgg_p1, span_tol) 
        deallocate(B1cgg_p1)
    
        ! Find basis function values at WQ points
        allocate(B0wq_p1(nb_ctrlpts_p1, obj%nb_qp_wq), B1wq_p1(nb_ctrlpts_p1, obj%nb_qp_wq))
        call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, obj%nb_qp_wq, obj%qp_position, B0wq_p1, B1wq_p1, span_tol) 
        deallocate(B1wq_p1)

        ! ------------------------------------
        ! Integrals and Weights
        ! ------------------------------------
        allocate(Bcgg_p0_int(obj%nb_ctrlpts, obj%nb_qp_cgg))
        allocate(Bcgg_p1_int(nb_ctrlpts_p1, obj%nb_qp_cgg))
        Bcgg_p0_int = 0; Bcgg_p1_int = 0
        where (abs(B0cgg_p0).gt.tol) Bcgg_p0_int = 1
        where (abs(B0cgg_p1).gt.tol) Bcgg_p1_int = 1

        allocate(IIshape(nb_ctrlpts_p1, obj%nb_ctrlpts))
        IIshape = matmul(Bcgg_p1_int, transpose(Bcgg_p0_int))

        ! Compute B0_p1 * W * B0_p0'
        allocate(II(nb_ctrlpts_p1, obj%nb_ctrlpts))
        call gemm_AWB(1, nb_ctrlpts_p1, obj%nb_qp_cgg, B0cgg_p1, obj%nb_ctrlpts, obj%nb_qp_cgg, &
                        B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, obj%nb_ctrlpts, II) 

        ! Compute W0
        call wq_solve_weights(obj%nb_ctrlpts, obj%B0shape, nb_ctrlpts_p1, obj%nb_qp_wq, B0wq_p1, II, IIshape, W00)
        W01 = W00

        ! Compute = B0_p1 * W * B1_p0'
        call gemm_AWB(1, nb_ctrlpts_p1, obj%nb_qp_cgg, B0cgg_p1, obj%nb_ctrlpts, obj%nb_qp_cgg, &
                            B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, obj%nb_ctrlpts, II)
                        
        ! Compute W1
        call wq_solve_weights(obj%nb_ctrlpts, obj%B1shape, nb_ctrlpts_p1, obj%nb_qp_wq, B0wq_p1, II, IIshape, W11)   
        W10 = W11            
    
    end subroutine wq_basis_weights_method2

    subroutine wq_basis_weights_dense2coo(obj)
        !! Computes basis and weights in COO format.
        !! Object needs to be initialized
            
        implicit none 
        ! Input / output data
        ! --------------------
        type(wq), pointer :: obj
        
        ! Local data
        ! ----------
        integer :: i, j, c, nbel, method
        double precision, dimension(:,:), allocatable :: B0, B1, W00, W01, W10, W11

        type(wq), pointer :: obj_m
        integer :: nbel_m, size_kv_m
        double precision, allocatable, dimension(:) :: knotvector_m, nodes_m
        double precision, allocatable, dimension(:,:) :: B0_m, B1_m, W00_m, W01_m, W10_m, W11_m

        method = obj%method
        if (.not.any((/1, 2/).eq.method)) stop 'Only method 1 or 2 are allowed'

        nbel = obj%size_nodes - 1
        call wq_get_qp_positions(obj)
        call wq_get_B_shape(obj)
        allocate(obj%data_B0(obj%nnz_B), obj%data_B1(obj%nnz_B), &
                obj%data_W00(obj%nnz_B), obj%data_W01(obj%nnz_B), &
                obj%data_W10(obj%nnz_B), obj%data_W11(obj%nnz_B), &
                obj%indices(obj%nnz_B, 2))

        if ((nbel.le.obj%degree+3).or.(.not.obj%isuniform)) then 
            
            allocate(B0(obj%size_kv-obj%degree-1, obj%nb_qp_wq), B1(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
            W00(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W01(obj%size_kv-obj%degree-1, obj%nb_qp_wq), &
            W10(obj%size_kv-obj%degree-1, obj%nb_qp_wq), W11(obj%size_kv-obj%degree-1, obj%nb_qp_wq))

            if (method.eq.1) then
                call wq_basis_weights_method1(obj, B0, B1, W00, W01, W10, W11)
            else if (method.eq.2) then
                call wq_basis_weights_method2(obj, B0, B1, W00, W01, W10, W11)
            end if

            c = 0
            do i = 1, obj%nb_ctrlpts
                do j = obj%B1shape(i, 1), obj%B1shape(i, 2)
                    c = c + 1
                    obj%data_B0(c) = B0(i, j)
                    obj%data_W00(c) = W00(i, j)
                    obj%data_W01(c) = W01(i, j)
                    obj%data_B1(c) = B1(i, j)
                    obj%data_W10(c) = W10(i, j)
                    obj%data_W11(c) = W11(i, j)
                    obj%indices(c, :) = [i, j]
                end do
            end do
            
        else
            ! Get model 
            ! ---------
            nbel_m = obj%degree + 3
            size_kv_m = nbel_m + 2*obj%degree +  1
            
            allocate(nodes_m(size_kv_m+1), knotvector_m(size_kv_m))
            call create_uniform_knotvector(obj%degree, nbel_m, nodes_m, knotvector_m)
            call wq_initialize(obj_m, obj%degree, size_kv_m, knotvector_m, method)
            call wq_get_qp_positions(obj_m)
            call wq_get_B_shape(obj_m)

            allocate(B0_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq), B1_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq), &
                    W00_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq), W01_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq), &
                    W10_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq), W11_m(obj_m%nb_ctrlpts, obj_m%nb_qp_wq))

            if (method.eq.1) then
                call wq_basis_weights_method1(obj_m, B0_m, B1_m, W00_m, W01_m, W10_m, W11_m)
            else if (method.eq.2) then
                call wq_basis_weights_method2(obj_m, B0_m, B1_m, W00_m, W01_m, W10_m, W11_m)
            end if

            B1_m = B1_m * nbel/nbel_m
            W00_m = W00_m * nbel_m/nbel
            W01_m = W01_m * nbel_m/nbel
                            
            ! Transfer data 
            ! ----------------
            ! Set p + 1 first functions
            c = 0
            do i = 1, obj%degree + 1
                do j = obj%B1shape(i, 1), obj%B1shape(i, 2)
                    c = c + 1
                    obj%data_B0(c) = B0_m(i, j)
                    obj%data_W00(c) = W00_m(i, j)
                    obj%data_W01(c) = W01_m(i, j)
                    obj%data_B1(c) = B1_m(i, j)
                    obj%data_W10(c) = W10_m(i, j)
                    obj%data_W11(c) = W11_m(i, j)
                    obj%indices(c, :) = [i, j]
                end do
            end do

            ! Set repeated functions 
            do i = obj%degree+2, obj%nb_ctrlpts-obj%degree-1
                do j = 1, obj_m%B0shape(obj%degree+2, 2) - obj_m%B0shape(obj%degree+2, 1) + 1 
                    c = c + 1
                    obj%data_B0(c) = B0_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%data_W00(c) = W00_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%data_W01(c) = W01_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%data_B1(c) = B1_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%data_W10(c) = W10_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%data_W11(c) = W11_m(obj%degree+2, obj_m%B0shape(obj%degree+2, 1) + j - 1)
                    obj%indices(c, :) = [i, obj%B1shape(i, 1) + j - 1]
                end do
            end do

            ! Set p + 1 last functions
            do i = obj%degree + 1, 1, -1 
                do j = obj%B1shape(i, 2), obj%B1shape(i, 1), -1
                    c = c + 1
                    obj%data_B0(c) = B0_m(i, j)
                    obj%data_W00(c) = W00_m(i, j)
                    obj%data_W01(c) = W01_m(i, j)
                    obj%data_B1(c) = -B1_m(i, j)
                    obj%data_W10(c) = -W10_m(i, j)
                    obj%data_W11(c) = -W11_m(i, j)
                    obj%indices(c, :) = [obj%nb_ctrlpts - i + 1, obj%nb_qp_wq - j + 1]
                end do
            end do

        end if

    end subroutine wq_basis_weights_dense2coo

end module wq_basis_weights
