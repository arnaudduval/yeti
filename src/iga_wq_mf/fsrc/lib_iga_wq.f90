! ==========================
! module :: iga and wq methods  
! author :: Joaquin Cornejo
! Disclaimer : all these routines assume an open  
!              knot-vector and at least 2 elements
! ==========================

subroutine wq_set_properties(degree, size_kv, maxrule, & 
                            size_nodes, nb_ctrlpts, nb_qp_wq, nb_qp_cgg)
    !! Sets constants used in IGA-WQ approach 
    
    use constants_iga_wq_mf
    implicit none 
    ! Input / out data 
    ! -----------------
    integer, intent(in) :: degree, size_kv, maxrule
    integer, intent(out) :: size_nodes, nb_ctrlpts, nb_qp_wq, nb_qp_cgg

    ! Find constants values : 
    ! Number of control points
    nb_ctrlpts = size_kv - degree - 1

    ! Number of nodes
    size_nodes = nb_ctrlpts - degree + 1

    ! Number of quadrature points in WQ approach
    nb_qp_wq = 2*(degree + r) + (size_nodes - 1)*(maxrule + 1) - 2*maxrule - 3  

    ! Number of Gauss quadrature points 
    nb_qp_cgg = (degree + 1)*(size_nodes - 1)

end subroutine wq_set_properties

subroutine verify_max_regularity(degree, size_kv, knotvector, is_not_maxregular)

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv
    double precision, intent(in) :: knotvector
    dimension :: knotvector(size_kv)

    logical, intent(out) :: is_not_maxregular
    
    ! Local data
    ! --------------------
    integer :: nbel, nb_ctrlpts, i
    double precision :: nodes
    dimension :: nodes(size_kv+1)
    double precision, dimension(:), allocatable :: newnodes, diffknot 

    ! Compute nodes
    call find_unique_vector(size_kv, knotvector, nodes)
    nbel = int(nodes(size_kv+1)) - 1
    nb_ctrlpts = size_kv - degree - 1

    ! Verify
    is_not_maxregular = .false.
    if (nbel+degree.ne.nb_ctrlpts) then 
        is_not_maxregular = .true.
    end if

    ! Verify if there is maximum regularity
    allocate(diffknot(nbel+1), newnodes(nbel+1))
    newnodes = nodes(1:nbel+1)

    call diff_vector(2, size(newnodes), newnodes, diffknot)
    do i = 1, size(newnodes)
        if (abs(diffknot(i)).gt.1.d-8) then
            is_not_maxregular = .true.
            exit
        end if
    end do

end subroutine verify_max_regularity

subroutine wq_get_qp_positions(degree, size_kv, nodes, maxrule, nb_qp, qp_pos)
    !! Gets quadrature points' positions (QP) in WQ approach 

    use constants_iga_wq_mf
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv, maxrule, nb_qp 
    double precision, intent(in) :: nodes
    dimension :: nodes(size_kv+1)

    double precision, intent(out) :: qp_pos 
    dimension :: qp_pos(nb_qp)

    ! Local data
    ! ---------------
    integer :: size_nodes
    double precision, allocatable, dimension(:) :: nodest
    double precision :: QPB, QPI ! Q.P. at boundaries, Q.P. at internal spans
    dimension :: QPB(degree+r), QPI(2+maxrule)
    integer :: i, j, k

    ! Initiaalize
    size_nodes = int(nodes(size_kv+1))
    allocate(nodest(size_nodes))
    nodest = nodes(1:size_nodes)

    ! Find values for the first boundary
    call linspace(nodest(1), nodest(2), size(QPB), QPB)
    do j = 1, size(QPB)
        qp_pos(j) = QPB(j)
    end do

    ! Find values for the last boundary 
    call linspace(nodest(size(nodest)-1), nodest(size(nodest)), size(QPB), QPB)
    do j = 1, size(QPB)
        qp_pos(size(qp_pos)-size(QPB)+j) = QPB(j)
    end do

    if (size(nodest).ge.4) then
        do i = 2, size(nodest)-2
            ! Find quadrature points for inner spans
            call linspace(nodest(i), nodest(i+1), size(QPI), QPI)

            ! Assign values
            do j = 1, size(QPI) 
                k = size(QPB) + (size(QPI) - 1)*(i - 2) + j - 1    
                qp_pos(k) = QPI(j)
            end do
        end do
    end if

end subroutine wq_get_qp_positions      

subroutine iga_get_qp_positions_weights(degree, size_kv, nodes, nb_qp, qp_pos, qp_weights)
    !! Gets quadrature points' positions and weights in IGA approach 

    implicit none 
    ! Input / Output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_qp
    double precision, intent(in) :: nodes
    dimension :: nodes(size_kv+1)

    double precision, intent(out) :: qp_pos, qp_weights
    dimension :: qp_pos(nb_qp), qp_weights(nb_qp)

    ! Local data 
    ! -------------
    integer :: nbel
    double precision :: GaussPdsCoord
    dimension :: GaussPdsCoord(2, degree+1)
    double precision :: xg(degree+1), wg(degree+1)
    integer :: i, j, k

    ! Set number of elements
    nbel = int(nodes(size_kv+1)) - 1 

    ! Find position and weight in isoparametric space
    call Gauss(size(xg), 1, GaussPdsCoord, 0)

    ! Split values 
    wg = GaussPdsCoord(1, :)
    xg = GaussPdsCoord(2, :)

    ! From isoparametric to parametric space
    do i = 1, nbel
        do j = 1, size(xg)
            k = (i - 1)*size(xg) + j
            qp_pos(k) = 0.5d0*(xg(j)/dble(nbel) + nodes(i+1) + nodes(i))
            qp_weights(k) = 0.5d0/dble(nbel)*wg(j)
        end do
    end do

end subroutine iga_get_qp_positions_weights 

subroutine iga_get_B_shape(degree, size_kv, nodes, knotvector, Bshape)
    !! Gets non zeros positions of B0 and B1 in IGA approach 
    
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: degree, size_kv
    double precision :: knotvector, nodes
    dimension :: knotvector(size_kv), nodes(size_kv+1)

    integer, intent(out) :: Bshape
    dimension :: Bshape(size_kv-degree-1, 2)

    ! Local data 
    ! --------------
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
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span)

    ! Get table of spans for each function
    do i = 1, nb_ctrlpts
        ! Find min 
        min_span = 1
        do j = 1, nbel
            if (any(table_functions_span(j, :).eq.i)) then
                    min_span = j
                    exit 
            end if
        end do 

        ! Find max
        max_span = nbel
        do j = nbel, 1, -1
            if (any(table_functions_span(j, :).eq.i)) then
                    max_span = j
                    exit 
            end if
        end do 

        ! Assigning values
        table_spans_function(i, :) = [min_span, max_span]
    end do 
            
    ! Set shape of B0 and B1
    do i = 1, nb_ctrlpts
        min_span = table_spans_function(i, 1)
        max_span = table_spans_function(i, 2)

        ! For B0 and B1
        min_knot = table_points_span(min_span, 1)
        max_knot = table_points_span(max_span, 2)

        Bshape(i, :) = [min_knot, max_knot]
    end do

end subroutine iga_get_B_shape

subroutine wq_get_B0_B1_shape(degree, size_kv, nodes, knotvector, maxrule, B0shape, B1shape)
    !! Gets non-zero positions of B0 and B1 in WQ approach

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: degree, maxrule, size_kv
    double precision, intent(in) :: nodes, knotvector
    dimension :: nodes(size_kv+1), knotvector(size_kv)

    integer, intent(out) :: B0shape, B1shape
    dimension :: B0shape(size_kv-degree-1, 2), B1shape(size_kv-degree-1, 2)

    ! Local data 
    ! --------------
    integer :: nbel, nb_ctrlpts, min_span, max_span, min_knot, max_knot
    integer, allocatable, dimension(:, :) :: table_points_span, table_functions_span, table_spans_function
    integer ::  i, j

    ! Set number of elements and control points
    nbel = int(nodes(size_kv+1)) - 1
    nb_ctrlpts = size_kv - degree - 1

    allocate(table_points_span(nbel, 2), &
            table_functions_span(nbel, degree+1), &
            table_spans_function(nb_ctrlpts, 2))

    ! Get table of points over the span
    table_points_span = 1
    table_points_span(1, 2) = degree + r 
    
    do i = 2, nbel - 1
        table_points_span(i, 1) = table_points_span(i-1, 2)
        table_points_span(i, 2) = table_points_span(i, 1) + 1 + maxrule
    end do

    table_points_span(nbel, 1) = table_points_span(nbel-1, 2)
    table_points_span(nbel, 2) = table_points_span(nbel, 1) + degree + r - 1

    ! Get table of functions on span 
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span)

    ! Get table of spans for each function
    do i = 1, nb_ctrlpts
        ! Find min 
        min_span = 1
        do j = 1, nbel
            if (any(table_functions_span(j, :).eq.i)) then
                    min_span = j
                    exit 
            end if
        end do 

        ! Find max
        max_span = nbel
        do j = nbel, 1, -1
            if (any(table_functions_span(j, :).eq.i)) then
                    max_span = j
                    exit 
            end if
        end do 

        ! Assigning values
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

        B0shape(i, :) = [min_knot, max_knot]
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

        B1shape(i, :) = [min_knot, max_knot]
    end do

end subroutine wq_get_B0_B1_shape

subroutine get_I_csr(nr, nc, nnz_B, indi_B, indj_B, nnz_I, indi_I, indj_I)
    !! Gets non-zero positions for I = B1 * B1.T in WQ and IGA approach
    !! B and I in CSR format

    use constants_iga_wq_mf
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in):: nr, nc, nnz_B, nnz_I
    integer, intent(in) :: indi_B, indj_B
    dimension :: indi_B(nr+1), indj_B(nnz_B)
    
    integer, intent(out) :: indi_I, indj_I
    dimension :: indi_I(nr+1), indj_I(nnz_I)

    ! Local data
    ! -------------
    double precision :: MB, MI, ones
    dimension :: MB(nr, nc), MI(nr, nr), ones(nnz_B)

    ! Initialize matrix B
    ones = 1.d0 
    call csr2dense(nnz_B, indi_B, indj_B, ones, nr, nc, MB)

    ! Compute I = B * B.T
    MI = matmul(MB, transpose(MB))
    
    ! Get CSR format
    call dense2csr(nr, nr, MI, nnz_I, indi_I, indj_I)

end subroutine get_I_csr

subroutine wq_get_qp_weights(nr_test, BBtest, nr_trial, nc, BBtrial, II, IIshape, weights)
    !! Gets the quadrature rules
    
    implicit none
    ! Input / output data
    ! -----------------------
    integer, intent(in) :: nr_trial, nr_test, nc
    integer, intent(in) :: BBtest, IIshape
    dimension :: BBtest(nr_test, 2), IIshape(nr_trial, nr_test)
    double precision, intent(in) ::  BBtrial, II
    dimension :: BBtrial(nr_trial, nc), II(nr_trial, nr_test)

    double precision, intent(out) :: weights
    dimension :: weights(nr_test, nc)

    ! Local data
    ! ----------------
    integer :: i, j, Pmin, Pmax, Fmin, Fmax

    ! Initialize
    weights = 0.d0  

    do i = 1, nr_test

        ! Find position of points within i-function support
        Pmin = BBtest(i, 1)
        Pmax = BBtest(i, 2)

        ! Find functions which intersect i-function support
        Fmin = 1
        do j = 1, nr_trial
            if (IIshape(j, i).gt.0) then
                Fmin = j
                exit 
            end if
        end do 

        Fmax = nr_trial
        do j = nr_trial, 1, -1
            if (IIshape(j, i).gt.0) then
                Fmax = j
                exit 
            end if
        end do 
        
        ! Solve linear system
        call solve_linear_system(Fmax-Fmin+1, Pmax-Pmin+1, BBtrial(Fmin:Fmax, Pmin:Pmax), &
                                II(Fmin:Fmax, i), weights(i, Pmin:Pmax))
    end do

end subroutine wq_get_qp_weights

subroutine wq_get_basis_weights_generalized(degree, size_kv, knotvector, nb_qp_wq, maxrule, B0, B1, W00, W01, W10, W11)
    !! Returns the basis and weights data at the quadrature points in WQ approach 

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, maxrule, nb_qp_wq
    double precision :: knotvector(size_kv)
    
    double precision, intent(out) :: B0, B1, W00, W01, W10, W11 
    dimension ::    B0(size_kv-degree-1, nb_qp_wq), B1(size_kv-degree-1, nb_qp_wq), &
                    W00(size_kv-degree-1, nb_qp_wq), W01(size_kv-degree-1, nb_qp_wq), &
                    W10(size_kv-degree-1, nb_qp_wq), W11(size_kv-degree-1, nb_qp_wq)

    ! Local data
    ! -------------        
    ! For p
    integer :: size_nodes, nb_ctrlpts_p0, nb_qp_wq_p0, nb_qp_cgg_p0
    double precision, allocatable, dimension(:) :: nodes_p0, qp_cgg_pos, qp_cgg_weights, qp_wq_pos
    double precision, allocatable, dimension(:,:) :: B0cgg_p0, B1cgg_p0    
    integer, allocatable, dimension(:,:) :: B0shape, B1shape, Bcgg_p0_int

    ! For p - 1
    integer :: size_kv_p1, degree_p1, nb_ctrlpts_p1, nb_qp_wq_p1, nb_qp_cgg_p1
    double precision, allocatable, dimension(:) :: knotvector_p1, nodes_p1
    double precision, allocatable, dimension(:,:) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
    integer, dimension(:,:), allocatable :: Bcgg_p1_int

    ! Integrals and weights
    double precision, allocatable, dimension(:,:) :: II
    integer, allocatable, dimension(:,:) :: IIshape

    ! --------------------------
    ! Degree p
    ! --------------------------
    ! Setup
    allocate(nodes_p0(size_kv+1))
    call find_unique_vector(size_kv, knotvector, nodes_p0)
    call wq_set_properties(degree, size_kv, maxrule, size_nodes, nb_ctrlpts_p0, nb_qp_wq_p0, nb_qp_cgg_p0)
    allocate(B0shape(nb_ctrlpts_p0, 2), B1shape(nb_ctrlpts_p0, 2))
    call wq_get_B0_B1_shape(degree, size_kv, nodes_p0, knotvector, maxrule, B0shape, B1shape)

    ! Find positions and weights in IGA approach
    allocate(qp_cgg_pos(nb_qp_cgg_p0), qp_cgg_weights(nb_qp_cgg_p0))
    call iga_get_qp_positions_weights(degree, size_kv, nodes_p0, nb_qp_cgg_p0, qp_cgg_pos, qp_cgg_weights)

    ! Find basis at Gauss quadrature points
    allocate(B0cgg_p0(nb_ctrlpts_p0, nb_qp_cgg_p0), B1cgg_p0(nb_ctrlpts_p0, nb_qp_cgg_p0))
    call get_basis(degree, size_kv, nodes_p0, knotvector, nb_qp_cgg_p0, qp_cgg_pos, B0cgg_p0, B1cgg_p0)

    ! Find quadrature points in WQ approach
    allocate(qp_wq_pos(nb_qp_wq_p0))
    call wq_get_qp_positions(degree, size_kv, nodes_p0, maxrule, nb_qp_wq_p0, qp_wq_pos) 

    ! Find basis at WQ quadrature points
    call get_basis(degree, size_kv, nodes_p0, knotvector, nb_qp_wq_p0, qp_wq_pos, B0, B1) 

    ! --------------------------
    ! Degree p - 1 
    ! -------------------------- 
    ! Setup
    degree_p1 = degree - 1
    size_kv_p1 = size_kv-2
    allocate(knotvector_p1(size_kv_p1))
    allocate(nodes_p1(size_kv_p1+1))
    knotvector_p1 = knotvector(2:size_kv-1)
    call find_unique_vector(size_kv_p1, knotvector_p1, nodes_p1)
    call wq_set_properties(degree_p1, size_kv_p1, maxrule, size_nodes, nb_ctrlpts_p1, nb_qp_wq_p1, nb_qp_cgg_p1)

    ! Find basis function values at Gauss points
    allocate(B0cgg_p1(nb_ctrlpts_p1, nb_qp_cgg_p0), B1cgg_p1(nb_ctrlpts_p1, nb_qp_cgg_p0))
    call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, nb_qp_cgg_p0, qp_cgg_pos, B0cgg_p1, B1cgg_p1) 
    deallocate(B1cgg_p1)

    ! Find basis function values at WQ points
    allocate(B0wq_p1(nb_ctrlpts_p1, nb_qp_wq_p0), B1wq_p1(nb_ctrlpts_p1, nb_qp_wq_p0))
    call get_basis(degree_p1, size_kv_p1, nodes_p1, knotvector_p1, nb_qp_wq_p0, qp_wq_pos, B0wq_p1, B1wq_p1) 
    deallocate(B1wq_p1)

    ! ------------------------------------
    ! Integrals and Weights
    ! ------------------------------------
    ! Initialize
    allocate(Bcgg_p0_int(nb_ctrlpts_p0, nb_qp_cgg_p0))
    allocate(Bcgg_p1_int(nb_ctrlpts_p1, nb_qp_cgg_p0))
    Bcgg_p0_int = 0; Bcgg_p1_int = 0
    where (abs(B0cgg_p0).gt.tol) Bcgg_p0_int = 1
    where (abs(B0cgg_p1).gt.tol) Bcgg_p1_int = 1

    ! ===========
    allocate(IIshape(nb_ctrlpts_p0, nb_ctrlpts_p0))
    IIshape = matmul(Bcgg_p0_int, transpose(Bcgg_p0_int))

    ! Compute B0cgg_p0 * Wcgg * B0cgg_p0.T
    allocate(II(nb_ctrlpts_p0, nb_ctrlpts_p0))
    call gemm_AWB(1, nb_ctrlpts_p0, nb_qp_cgg_p0, B0cgg_p0, nb_ctrlpts_p0, nb_qp_cgg_p0, &
                B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p0, nb_ctrlpts_p0, II)
                        
    ! For W00
    call wq_get_qp_weights(nb_ctrlpts_p0, B0shape, nb_ctrlpts_p0, nb_qp_wq_p0, B0, II, IIshape, W00)

    ! ----------------------
    ! Compute  B0cgg_p0 * Wcgg * B1cgg_p0.T
    call gemm_AWB(1, nb_ctrlpts_p0, nb_qp_cgg_p0, B0cgg_p0, nb_ctrlpts_p0, nb_qp_cgg_p0, & 
                B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p0, nb_ctrlpts_p0, II)

    ! For W10
    call wq_get_qp_weights(nb_ctrlpts_p0, B1shape, nb_ctrlpts_p0, nb_qp_wq_p0, B0, II, IIshape, W10)
    deallocate(IIshape, II)

    ! ===========
    allocate(IIshape(nb_ctrlpts_p1, nb_ctrlpts_p0))
    IIshape = matmul(Bcgg_p1_int, transpose(Bcgg_p0_int))

    ! Compute B0cgg_p1 * Wcgg * B0cgg_p0.transpose
    allocate(II(nb_ctrlpts_p1, nb_ctrlpts_p0))
    call gemm_AWB(1, nb_ctrlpts_p1, nb_qp_cgg_p0, B0cgg_p1, nb_ctrlpts_p0, nb_qp_cgg_p0, & 
                B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, nb_ctrlpts_p0, II)

    ! For W01
    call wq_get_qp_weights(nb_ctrlpts_p0, B0shape, nb_ctrlpts_p1, nb_qp_wq_p0, B0wq_p1, II, IIshape, W01)
    
    ! ----------------------
    ! Compute B0cgg_p1 * Wcgg * B1cgg_p0.transpose
    call gemm_AWB(1, nb_ctrlpts_p1, nb_qp_cgg_p0, B0cgg_p1, nb_ctrlpts_p0, nb_qp_cgg_p0, &
                B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, nb_ctrlpts_p0, II)
                    
    ! For W11
    call wq_get_qp_weights(nb_ctrlpts_p0, B1shape, nb_ctrlpts_p1, nb_qp_wq_p0, B0wq_p1, II, IIshape, W11)               
    deallocate(IIshape, II)

end subroutine wq_get_basis_weights_generalized

! ===========================================================
! ===========================================================
module iga_basis_weights

    implicit none
    type :: iga
        ! Inputs :
        ! ------------
        integer :: degree, size_kv
        double precision, dimension(:), pointer :: knotvector
        
        ! Outputs :
        ! ---------
        double precision, dimension(:), pointer :: qp_pos, qp_weights, data_B0, data_B1
        integer, dimension(:,:), pointer :: data_ind
        integer :: nnz_B, nnz_I
        
        ! Local :
        ! ---------
        integer :: nb_ctrlpts, nb_qp, size_nodes
        integer, dimension(:,:), pointer :: Bshape
        double precision, dimension(:), pointer :: nodes

    end type iga

    contains

    subroutine iga_initialize(object, degree, size_kv, knotvector)
            
        implicit none
        ! Input / output
        ! ---------------
        integer, intent(in) :: degree, size_kv
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(iga), pointer :: object

        ! Set properties
        allocate(object)
        allocate(object%knotvector(size_kv))
        object%degree = degree
        object%size_kv = size_kv
        object%nb_ctrlpts = size_kv - degree - 1 
        object%knotvector = knotvector

        ! Compute unique knots
        allocate(object%nodes(size_kv+1))
        call find_unique_vector(object%size_kv, object%knotvector, object%nodes)
        object%size_nodes = int(object%nodes(size_kv+1))
        object%nb_qp = (degree + 1)*(object%size_nodes - 1)
        
        ! Get quadrature points position
        allocate(object%qp_pos(object%nb_qp), object%qp_weights(object%nb_qp))
        call iga_get_qp_positions_weights(object%degree, object%size_kv, object%nodes, & 
                                    object%nb_qp, object%qp_pos, object%qp_weights)

        ! Set non zeros values of B0 and B1
        allocate(object%Bshape(object%nb_ctrlpts, 2))
        call iga_get_B_shape(object%degree, object%size_kv, object%nodes, object%knotvector, object%Bshape)  

        ! Get necessary size of data arrays
        object%nnz_B = object%nb_qp*(object%degree + 1) 
        object%nnz_I = (2*object%degree+1)*object%nb_ctrlpts - object%degree*(object%degree+1)

    end subroutine iga_initialize

    subroutine iga_get_basis_weights(object, degree, size_kv, knotvector)
            
        implicit none 
        ! Input / output data
        ! --------------------
        integer, intent(in) :: degree, size_kv
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(iga), pointer :: object
        
        ! Local data
        ! ----------
        double precision, dimension(:,:), allocatable :: B0, B1
        integer :: i, j, count

        ! Create object
        call iga_initialize(object, degree, size_kv, knotvector)

        ! Allocate variables
        allocate(B0(object%nb_ctrlpts, object%nb_qp), B1(object%nb_ctrlpts, object%nb_qp))

        ! Get basis and weights 
        call get_basis(object%degree, object%size_kv, object%nodes, object%knotvector, &
                        object%nb_qp, object%qp_pos, B0, B1)

        ! Set size of properties
        allocate(object%data_B0(object%nnz_B), object%data_B1(object%nnz_B), object%data_ind(object%nnz_B, 2))

        ! Assign values
        count = 0
        do i = 1, object%nb_ctrlpts
            do j = object%Bshape(i, 1), object%Bshape(i, 2)
                count = count + 1
                object%data_B0(count) = B0(i, j)
                object%data_B1(count) = B1(i, j)
                object%data_ind(count, :) = [i, j]
            end do
        end do

    end subroutine iga_get_basis_weights

end module iga_basis_weights

module wq_basis_weights
   
    implicit none
    integer, parameter :: maxrule= 1
    type :: wq
        ! Inputs :
        ! ------------
        integer :: degree, size_kv
        double precision, dimension(:), pointer :: knotvector
        
        ! Outputs :
        ! ---------
        double precision, dimension(:), pointer ::  qp_pos, data_B0, data_B1, & 
                                                    data_W00, data_W01, data_W10, data_W11
        integer, dimension(:,:), pointer :: data_indices
        integer :: nnz_B, nnz_I
        
        ! Local :
        ! ---------
        integer :: nb_ctrlpts, nb_qp_wq, nb_qp_cgg, size_nodes
        integer, dimension(:, :), pointer :: Bshape
        double precision, dimension(:), pointer :: nodes
        
    end type wq

    contains

    subroutine wq_initialize(object, degree, size_kv, knotvector)
            
        implicit none
        ! Input / output
        ! ---------------
        integer, intent(in) :: degree, size_kv
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(wq), pointer :: object

        ! Local data
        ! --------------
        integer, allocatable, dimension(:,:) :: dummy
        integer :: count, i

        ! Set properties
        allocate(object)
        allocate(object%knotvector(size_kv))
        object%degree = degree
        object%size_kv = size_kv
        object%knotvector = knotvector

        ! From inputs to local data
        allocate(object%nodes(size_kv+1))
        call find_unique_vector(object%size_kv, object%knotvector, object%nodes)
        call wq_set_properties(object%degree, object%size_kv, maxrule, object%size_nodes, &
                                object%nb_ctrlpts, object%nb_qp_wq, object%nb_qp_cgg)

        allocate(object%Bshape(object%nb_ctrlpts, 2), dummy(object%nb_ctrlpts, 2))
        call wq_get_B0_B1_shape(object%degree, object%size_kv, object%nodes, & 
                    object%knotvector, maxrule, dummy, object%Bshape)  
        deallocate(dummy)
        
        ! Get quadrature points position
        allocate(object%qp_pos(object%nb_qp_wq))
        call wq_get_qp_positions(object%degree, object%size_kv, object%nodes, maxrule, object%nb_qp_wq, object%qp_pos)

        ! Get necessary size of data arrays
        count = 0
        do i = 1, object%nb_ctrlpts
            count = count + object%Bshape(i, 2) - object%Bshape(i, 1) + 1
        end do

        ! Assign values
        object%nnz_B = count
        object%nnz_I = (2*object%degree+1)*object%nb_ctrlpts - object%degree*(object%degree+1)

    end subroutine wq_initialize

    subroutine wq_get_basis_weights(object, degree, size_kv, knotvector)
            
        implicit none 
        ! Input / output data
        ! --------------------
        integer, intent(in) :: degree, size_kv
        double precision, intent(in) :: knotvector
        dimension :: knotvector(size_kv)
        type(wq), pointer :: object
        
        ! Local data
        ! ----------
        logical :: is_not_maxregular
        integer :: i, j, count
        integer :: nbel, nbel_m, nb_ctrlpts_m, nb_qp_wq_m, size_kv_m, dummy1, dummy2
        double precision, allocatable, dimension(:) :: knotvector_m, nodes_m

        double precision, dimension(:,:), allocatable :: B0, B1, W00, W01, W10, W11
        double precision, allocatable, dimension(:,:) :: B0_m, B1_m, W00_m, W01_m, W10_m, W11_m
        integer, allocatable, dimension(:,:) :: Bshape_m, dummy

        ! Create object
        call wq_initialize(object, degree, size_kv, knotvector)
        nbel = object%size_nodes - 1
        call verify_max_regularity(degree, size_kv, knotvector, is_not_maxregular)

        if ((nbel.le.degree+3).or.(is_not_maxregular)) then 

            ! Allocate variables
            allocate(B0(object%nb_ctrlpts, object%nb_qp_wq), B1(object%nb_ctrlpts, object%nb_qp_wq), &
                    W00(object%nb_ctrlpts, object%nb_qp_wq), W01(object%nb_ctrlpts, object%nb_qp_wq), &
                    W10(object%nb_ctrlpts, object%nb_qp_wq), W11(object%nb_ctrlpts, object%nb_qp_wq))

            ! Get basis and weights 
            call wq_get_basis_weights_generalized(object%degree, object%size_kv, object%knotvector, & 
                                                object%nb_qp_wq, maxrule, B0, B1, W00, W01, W10, W11)

            ! Set size of properties
            allocate(object%data_B0(object%nnz_B), object%data_B1(object%nnz_B), &
                    object%data_W00(object%nnz_B), object%data_W01(object%nnz_B), &
                    object%data_W10(object%nnz_B), object%data_W11(object%nnz_B), &
                    object%data_indices(object%nnz_B, 2))

            ! Assign values
            count = 0
            do i = 1, object%nb_ctrlpts
                do j = object%Bshape(i, 1), object%Bshape(i, 2)
                    count = count + 1
                    object%data_B0(count) = B0(i, j)
                    object%data_W00(count) = W00(i, j)
                    object%data_W01(count) = W01(i, j)
                    object%data_B1(count) = B1(i, j)
                    object%data_W10(count) = W10(i, j)
                    object%data_W11(count) = W11(i, j)
                    object%data_indices(count, :) = [i, j]
                end do
            end do
            
        else    
            ! Get model 
            ! --------------
            nbel_m = degree + 3
            size_kv_m = nbel_m + 2*degree +  1

            allocate(nodes_m(size_kv_m+1), knotvector_m(size_kv_m))
            call create_uniform_knotvector(degree, nbel_m, nodes_m, knotvector_m)
            call wq_set_properties(degree, size_kv_m, maxrule, dummy1, nb_ctrlpts_m, nb_qp_wq_m, dummy2)

            allocate(B0_m(nb_ctrlpts_m, nb_qp_wq_m), B1_m(nb_ctrlpts_m, nb_qp_wq_m), &
                    W00_m(nb_ctrlpts_m, nb_qp_wq_m), W01_m(nb_ctrlpts_m, nb_qp_wq_m), &
                    W10_m(nb_ctrlpts_m, nb_qp_wq_m), W11_m(nb_ctrlpts_m, nb_qp_wq_m))
            allocate(Bshape_m(nb_ctrlpts_m, 2), dummy(nb_ctrlpts_m, 2))
            
            ! Compute model
            call wq_get_basis_weights_generalized(degree, size_kv_m, knotvector_m, nb_qp_wq_m, maxrule, &
                                                B0_m, B1_m, W00_m, W01_m, W10_m, W11_m)
            B1_m = B1_m * nbel/nbel_m
            W00_m = W00_m * nbel_m/nbel
            W01_m = W01_m * nbel_m/nbel

            call wq_get_B0_B1_shape(degree, size_kv_m, nodes_m, knotvector_m, maxrule, dummy, Bshape_m)
            deallocate(dummy)
                            
            ! Transfer data 
            ! ----------------
            allocate(object%data_B0(object%nnz_B), object%data_B1(object%nnz_B), &
                    object%data_W00(object%nnz_B), object%data_W01(object%nnz_B), &
                    object%data_W10(object%nnz_B), object%data_W11(object%nnz_B), &
                    object%data_indices(object%nnz_B, 2))

            ! Set p + 1 first functions
            count = 0
            do i = 1, degree + 1
                do j = object%Bshape(i, 1), object%Bshape(i, 2)
                    count = count + 1
                    object%data_B0(count) = B0_m(i, j)
                    object%data_W00(count) = W00_m(i, j)
                    object%data_W01(count) = W01_m(i, j)
                    object%data_B1(count) = B1_m(i, j)
                    object%data_W10(count) = W10_m(i, j)
                    object%data_W11(count) = W11_m(i, j)
                    object%data_indices(count, :) = [i, j]
                end do
            end do

            ! Set repeated functions 
            do i = degree+2, object%nb_ctrlpts-degree-1
                do j = 1, Bshape_m(degree+2, 2) - Bshape_m(degree+2, 1) + 1 
                    count = count + 1
                    object%data_B0(count) = B0_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_W00(count) = W00_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_W01(count) = W01_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_B1(count) = B1_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_W10(count) = W10_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_W11(count) = W11_m(degree+2, Bshape_m(degree+2, 1) + j - 1)
                    object%data_indices(count, :) = [i, object%Bshape(i, 1) + j - 1]
                end do
            end do

            ! Set p + 1 last functions
            do i = degree + 1, 1, -1 
                do j = object%Bshape(i, 2), object%Bshape(i, 1), -1
                    count = count + 1
                    object%data_B0(count) = B0_m(i, j)
                    object%data_W00(count) = W00_m(i, j)
                    object%data_W01(count) = W01_m(i, j)
                    object%data_B1(count) = -B1_m(i, j)
                    object%data_W10(count) = -W10_m(i, j)
                    object%data_W11(count) = -W11_m(i, j)
                    object%data_indices(count, :) = [object%nb_ctrlpts - i + 1, object%nb_qp_wq - j + 1]
                end do
            end do

        end if

    end subroutine wq_get_basis_weights

end module wq_basis_weights
