! ==========================
! module :: iga and wq methods  
! author :: Joaquin Cornejo
! hypothesis : all these routines assume an open and uniform 
!              knot-vector and at least 2 elements
! modules :: algebra(linspace, solve_system, product_AWB),
!            bspline(set_table_functions_spans, get_parametric_nodes, 
!                   get_knotvector, get_basis)
!            GaussLegendre.f90 
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
    
    ! Assign values
    do j = 1, size(QPB)
        qp_pos(j) = QPB(j)
    end do

    ! Find values for the last boundary 
    call linspace(nodest(size(nodest)-1), nodest(size(nodest)), size(QPB), QPB)
    
    ! Assign values
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
    !! (uniform knot-vector)
    
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
    integer :: nbel, nb_ctrlpts
    integer :: min_span, max_span, min_knot, max_knot
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
        table_points_span(i, 1) = table_points_span(i-1, 1) + degree + 1
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
    !! uniform knot-vector

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
    integer :: nbel, nb_ctrlpts, nb_qp
    integer :: min_span, max_span, min_knot, max_knot
    integer, allocatable, dimension(:, :) :: table_points_span, table_functions_span, table_spans_function
    integer ::  i, j

    ! Set number of elements and control points
    nbel = int(nodes(size_kv+1)) - 1
    nb_ctrlpts = size_kv - degree - 1
    nb_qp = 2*(degree + r) + nbel*(maxrule + 1) - 2*maxrule - 3

    allocate(table_points_span(nbel, 2), &
            table_functions_span(nbel, degree+1), &
            table_spans_function(nb_ctrlpts, 2))

    ! Get table of points over the span
    table_points_span = 1
    table_points_span(1, 2) = degree + r 
    table_points_span(nbel, 1) = nb_qp + 1 - (degree + r)
    table_points_span(nbel, 2) = nb_qp
    
    do i = 2, nbel - 1
        table_points_span(i, 1) = table_points_span(i-1, 2)
        table_points_span(i, 2) = table_points_span(i, 1) + 1 + maxrule
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

        ! For B0
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

        ! For B1
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

subroutine get_I_csr(nbr, nbc, nnz_B, indi_B, indj_B, &
                    nnz_I, indi_I, indj_I)
    !! Gets non-zero positions for I = B1 * B1.T in WQ and IGA approach
    !! B and I in CSR format

    use constants_iga_wq_mf
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in):: nbr, nbc, nnz_B, nnz_I
    integer, intent(in) :: indi_B, indj_B
    dimension :: indi_B(nbr+1), indj_B(nnz_B)
    
    integer, intent(out) :: indi_I, indj_I
    dimension :: indi_I(nbr+1), indj_I(nnz_I)

    ! Local data
    ! -------------
    double precision :: MB, MI, ones
    dimension :: MB(nbr, nbc), MI(nbr, nbr), ones(nnz_I)

    ! Initialize matrix B
    ones = 1.d0 
    call csr2dense(nnz_B, indi_B, indj_B, ones, nbr, nbc, MB)

    ! Compute I = B * B.T
    MI = matmul(MB, transpose(MB))
    
    ! Get CSR format
    call dense2csr(nbr, nbr, MI, nnz_I, indi_I, indj_I)

end subroutine get_I_csr

subroutine wq_get_qp_weights(nbr, BBshape, nbr_trial, nbc, BBtrial, II, IIshape, weights)
    !! Gets the quadrature rules
    
    implicit none
    ! Input / output data
    ! -----------------------
    integer, intent(in) :: nbr_trial, nbr, nbc
    integer, intent(in) :: BBshape, IIshape
    dimension :: BBshape(nbr, 2), IIshape(nbr_trial, nbr)
    double precision, intent(in) ::  BBtrial, II
    dimension :: BBtrial(nbr_trial, nbc), II(nbr_trial, nbr)

    double precision, intent(out) :: weights
    dimension :: weights(nbr, nbc)

    ! Local data
    ! ----------------
    integer :: i, j
    integer :: Pmin, Pmax, Fmin, Fmax

    ! Initialize
    weights = 0.d0  

    do i = 1, nbr

        ! Find position of points within i-function support
        Pmin = BBshape(i, 1)
        Pmax = BBshape(i, 2)

        ! Find functions which intersect i-function support
        Fmin = 1
        do j = 1, nbr_trial
            if (IIshape(j, i).gt.0) then
                Fmin = j
                exit 
            end if
        end do 

        Fmax = nbr_trial
        do j = nbr_trial, 1, -1
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

subroutine wq_get_basis_weights_generalized(degree, size_kv, knotvector, nb_qp_wq, maxrule, &
                                            B0, B1, W00, W01, W10, W11)
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
    double precision, allocatable, dimension(:, :) :: B0cgg_p0, B1cgg_p0    
    integer, allocatable, dimension(:, :) :: B0shape, B1shape
    integer, allocatable, dimension(:, :) :: Bint_p0

    ! For p - 1
    double precision, allocatable, dimension(:) :: knotvector_p1, nodes_p1
    integer :: size_kv_p1, degree_p1, nb_ctrlpts_p1, nb_qp_wq_p1, nb_qp_cgg_p1
    double precision, allocatable, dimension(:,:) :: B0cgg_p1, B1cgg_p1, B0wq_p1, B1wq_p1
    integer, dimension(:,:), allocatable :: Bint_p1

    ! Integrals and weights
    double precision, allocatable, dimension(:,:) :: II
    integer, allocatable, dimension(:, :) :: MIint

    ! Compute unique knot vector
    allocate(nodes_p0(size_kv+1))
    call find_unique_vector(size_kv, knotvector, nodes_p0)

    ! Compute size of arrays
    call wq_set_properties(degree, size_kv, maxrule, & 
                        size_nodes, nb_ctrlpts_p0, nb_qp_wq_p0, nb_qp_cgg_p0)

    ! Get non-zero values shape
    allocate(B0shape(nb_ctrlpts_p0, 2), B1shape(nb_ctrlpts_p0, 2))
    call wq_get_B0_B1_shape(degree, size_kv, nodes_p0, knotvector, maxrule, B0shape, B1shape)

    ! --------------------------
    ! Degree p
    ! --------------------------
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
    ! Find knot-vector in WQ approach
    degree_p1 = degree - 1
    size_kv_p1 = size_kv-2
    allocate(knotvector_p1(size_kv_p1))
    knotvector_p1 = knotvector(2:size_kv-1)
    allocate(nodes_p1(size_kv_p1+1))
    call find_unique_vector(size_kv_p1, knotvector_p1, nodes_p1)
    
    ! Compute size of arrays
    call wq_set_properties(degree_p1, size_kv_p1, maxrule, & 
                        size_nodes, nb_ctrlpts_p1, nb_qp_wq_p1, nb_qp_cgg_p1)

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
    ! Initialiaze
    allocate(Bint_p0(nb_ctrlpts_p0, nb_qp_cgg_p0))
    allocate(Bint_p1(nb_ctrlpts_p1, nb_qp_cgg_p0))
    Bint_p0 = 0
    Bint_p1 = 0
    where (abs(B0cgg_p0).gt.tol)
        Bint_p0 = 1
    end where
    where (abs(B0cgg_p1).gt.tol)
        Bint_p1 = 1
    end where

    ! ----------------------
    allocate(MIint(nb_ctrlpts_p0, nb_ctrlpts_p0))
    MIint = matmul(Bint_p0, transpose(Bint_p0))

    ! I00 = B0cgg_p0 * Wcgg * B0cgg_p0.transpose
    allocate(II(nb_ctrlpts_p0, nb_ctrlpts_p0))
    call gemm_AWB(1, nb_ctrlpts_p0, nb_qp_cgg_p0, B0cgg_p0, nb_ctrlpts_p0, nb_qp_cgg_p0, &
                B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p0, nb_ctrlpts_p0, II)
                        
    ! For W00
    call wq_get_qp_weights(nb_ctrlpts_p0, B0shape, nb_ctrlpts_p0, nb_qp_wq_p0, B0, II, MIint, W00)

    ! ----------------------
    ! I10 = B0cgg_p0 * Wcgg * B1cgg_p0.transpose
    call gemm_AWB(1, nb_ctrlpts_p0, nb_qp_cgg_p0, B0cgg_p0, nb_ctrlpts_p0, nb_qp_cgg_p0,  & 
                B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p0, nb_ctrlpts_p0, II)

    ! For W10
    call wq_get_qp_weights(nb_ctrlpts_p0, B1shape, nb_ctrlpts_p0, nb_qp_wq_p0, B0, II, MIint, W10)
    deallocate(MIint, II)

    ! ----------------------
    allocate(MIint(nb_ctrlpts_p1, nb_ctrlpts_p0))
    MIint = matmul(Bint_p1, transpose(Bint_p0))

    ! I01 = B0cgg_p1 * Wcgg * B0cgg_p0.transpose
    allocate(II(nb_ctrlpts_p1, nb_ctrlpts_p0))
    call gemm_AWB(1, nb_ctrlpts_p1, nb_qp_cgg_p0, B0cgg_p1, &
    nb_ctrlpts_p0, nb_qp_cgg_p0,  B0cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, nb_ctrlpts_p0, II)

    ! For W01
    call wq_get_qp_weights(nb_ctrlpts_p0, B0shape, nb_ctrlpts_p1, nb_qp_wq_p0, B0wq_p1, II, MIint, W01)
    
    ! ----------------------
    ! I11 = B0cgg_p1 * Wcgg * B1cgg_p0.transpose
    call gemm_AWB(1, nb_ctrlpts_p1, nb_qp_cgg_p0, B0cgg_p1, nb_ctrlpts_p0, nb_qp_cgg_p0,  &
                B1cgg_p0, qp_cgg_weights, nb_ctrlpts_p1, nb_ctrlpts_p0, II)
                    
    ! For W11
    call wq_get_qp_weights(nb_ctrlpts_p0, B1shape, nb_ctrlpts_p1, nb_qp_wq_p0, B0wq_p1, II, MIint, W11)               
    deallocate(MIint, II)

end subroutine wq_get_basis_weights_generalized

subroutine wq_get_basis_weights_model(degree, nbel, nbel_m, size_kv_m, knotvector_m, nb_qp_wq_m, maxrule, &
                                    B0_m, B1_m, W00_m, W01_m, W10_m, W11_m)
    !! Gets model where nb_el = degree + 3, nb_qp_wq = 4*degree + 2*r + 1
    
    use constants_iga_wq_mf
    implicit none
    ! Input / output data 
    ! --------------------
    integer, parameter :: multiplicity = 1
    integer, intent(in) ::  degree, nbel, nbel_m, maxrule, nb_qp_wq_m, size_kv_m
    double precision :: knotvector_m
    dimension :: knotvector_m(size_kv_m)

    double precision, intent(out) :: B0_m, B1_m 
    dimension ::    B0_m(degree+2, nb_qp_wq_m), & 
                    B1_m(degree+2, nb_qp_wq_m)
    double precision, intent(out) :: W00_m, W01_m, W10_m, W11_m
    dimension ::    W00_m(degree+2, nb_qp_wq_m), &
                    W01_m(degree+2, nb_qp_wq_m), &
                    W10_m(degree+2, nb_qp_wq_m), &
                    W11_m(degree+2, nb_qp_wq_m)
    
    ! Local data 
    ! ---------------
    integer :: i, nb_ctrlpts_m
    double precision, allocatable, dimension(:,:) :: B0_mt, B1_mt, W00_mt, W01_mt, W10_mt, W11_mt

    ! Initialize
    nb_ctrlpts_m = size_kv_m - degree - 1

    allocate(B0_mt(nb_ctrlpts_m, nb_qp_wq_m), B1_mt(nb_ctrlpts_m, nb_qp_wq_m), &
            W00_mt(nb_ctrlpts_m, nb_qp_wq_m), W01_mt(nb_ctrlpts_m, nb_qp_wq_m), &
            W10_mt(nb_ctrlpts_m, nb_qp_wq_m), W11_mt(nb_ctrlpts_m, nb_qp_wq_m))

    call wq_get_basis_weights_generalized(degree, size_kv_m, knotvector_m, nb_qp_wq_m, maxrule, &
                                        B0_mt, B1_mt, W00_mt, W01_mt, W10_mt, W11_mt)

    do i = 1, degree+2
        B0_m(i, :) = B0_mt(i, :)
        B1_m(i, :) = B1_mt(i, :) * nbel / nbel_m
        W00_m(i, :) = W00_mt(i, :) * nbel_m / nbel
        W01_m(i, :) = W01_mt(i, :) * nbel_m / nbel
        W10_m(i, :) = W10_mt(i, :)
        W11_m(i, :) = W11_mt(i, :)
    end do

    deallocate(B0_mt, B1_mt, W00_mt, W01_mt, W10_mt, W11_mt)

end subroutine wq_get_basis_weights_model

! ===========================================================
! ===========================================================
module iga_basis_weights

    implicit none
    type :: iga
        ! Inputs :
        ! ------------
        integer :: degree, size_kv, size_nodes
        double precision, dimension(:), pointer :: knotvector
        
        ! Outputs :
        ! ---------
        double precision, dimension(:), pointer :: qp_pos, qp_weights
        double precision, dimension(:), pointer :: data_B0, data_B1
        integer, dimension(:, :), pointer :: data_ind
        integer :: nnz_B, nnz_I
        
        ! Local :
        ! ---------
        integer, dimension(:, :), pointer :: Bshape
        double precision, dimension(:), pointer :: nodes
        integer :: nb_ctrlpts, nb_qp

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
        object%degree = degree
        object%size_kv = size_kv
        object%nb_ctrlpts = size_kv - degree - 1 
        allocate(object%knotvector(size_kv))
        object%knotvector = knotvector

        ! Compute unique knots
        allocate(object%nodes(size_kv+1))
        call find_unique_vector(object%size_kv, object%knotvector, object%nodes)
        object%size_nodes = int(object%nodes(size_kv+1))
        object%nb_qp = (degree + 1)*(object%size_nodes - 1)
        
        ! Get quadrature points position
        allocate(object%qp_pos(object%nb_qp))
        allocate(object%qp_weights(object%nb_qp))
        call iga_get_qp_positions_weights(object%degree, object%size_kv, object%nodes, & 
                                    object%nb_qp, object%qp_pos, object%qp_weights)

        ! Set non zeros values of B0 and B1
        allocate(object%Bshape(object%nb_ctrlpts, 2))
        call iga_get_B_shape(object%degree, object%size_kv, object%nodes, &
                            object%knotvector, object%Bshape)  

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
        allocate(B0(object%nb_ctrlpts, object%nb_qp))
        allocate(B1(object%nb_ctrlpts, object%nb_qp))

        ! Get basis and weights 
        call get_basis(object%degree, object%size_kv, object%nodes, object%knotvector, &
                        object%nb_qp, object%qp_pos, B0, B1)

        ! Set size of properties
        allocate(object%data_B0(object%nnz_B))
        allocate(object%data_B1(object%nnz_B))
        allocate(object%data_ind(object%nnz_B, 2))

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
    integer, parameter :: multiplicity=1, maxrule=1
    type :: wq
        ! Inputs :
        ! ------------
        integer :: degree, size_kv, size_nodes
        double precision, dimension(:), pointer :: knotvector
        
        ! Outputs :
        ! ---------
        double precision, dimension(:), pointer :: qp_pos
        double precision, dimension(:), pointer ::  data_B0, data_B1, data_W00, & 
                                                    data_W01, data_W10, data_W11
        integer, dimension(:, :), pointer :: data_ind
        integer :: nnz_B, nnz_I
        
        ! Local :
        ! ---------
        integer, dimension(:, :), pointer :: Bshape
        double precision, dimension(:), pointer :: nodes
        integer :: nb_ctrlpts, nb_qp_wq, nb_qp_cgg

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
        integer, allocatable, dimension(:, :) :: Bshape_temp
        integer :: count, i

        ! Set properties
        allocate(object)
        object%degree = degree
        object%size_kv = size_kv
        object%nb_ctrlpts = size_kv - degree - 1 
        allocate(object%knotvector(size_kv))
        object%knotvector = knotvector

        ! Compute unique knots
        allocate(object%nodes(size_kv+1))
        call find_unique_vector(object%size_kv, object%knotvector, object%nodes)
        object%size_nodes = int(object%nodes(size_kv+1))
        
        call wq_set_properties(object%degree, object%size_kv, maxrule, &
                        object%size_nodes, object%nb_ctrlpts, object%nb_qp_wq, object%nb_qp_cgg)
        
        ! Get quadrature points position
        allocate(object%qp_pos(object%nb_qp_wq))
        call wq_get_qp_positions(object%degree, object%size_kv, object%nodes, maxrule, object%nb_qp_wq, object%qp_pos)

        ! Set B0_shape and B1_shape
        allocate(object%Bshape(object%nb_ctrlpts, 2))
        allocate(Bshape_temp(object%nb_ctrlpts, 2))
        call wq_get_B0_B1_shape(object%degree, object%size_kv, object%nodes, & 
                    object%knotvector, maxrule, Bshape_temp, object%Bshape)  
        deallocate(Bshape_temp)

        ! Get necessary size of data arrays
        count = 0
        do i = 1, object%nb_ctrlpts
            count = count + object%Bshape(i, 2) - object%Bshape(i, 1) + 1
        end do

        ! Assign values
        object%nnz_B = count
        object%nnz_I = (2*object%degree+1)*object%nb_ctrlpts -object%degree*(object%degree+1)

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
        ! Model 
        logical :: isregular
        double precision :: diffknot
        dimension :: diffknot(size_kv)
        integer :: nbel, nbel_m, nb_ctrlpts_m, nb_qp_wq_m, size_kv_m, dummy1, dummy2
        double precision, allocatable, dimension(:) :: knotvector_m, nodes_m

        ! Variables for case nb_el <= p + 3
        double precision, dimension(:,:), allocatable :: B0, B1, W00, W01, W10, W11

        ! Variables for case nb_el > p + 3
        double precision, allocatable, dimension(:,:) :: B0_model, B1_model, &
                                                            W00_model, W01_model, W10_model, W11_model
        integer, allocatable, dimension(:,:) :: B_shape_model, B_shape_model_temp

        ! Loops
        integer :: i, j, count

        ! Create object
        call wq_initialize(object, degree, size_kv, knotvector)
        nbel = object%size_nodes - 1

        ! Verify if there is maximum regularity
        isregular = .true.
        call diff_vector(2, size_kv, knotvector, diffknot)
        do i = 1, size_kv
            if (abs(diffknot(i)).gt.1.d-8) then
                isregular = .false.
                exit
            end if
        end do

        if ((nbel.le.degree+3).or.(.not.isregular)) then 

            ! Allocate variables
            allocate(B0(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(B1(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W00(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W01(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W10(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W11(object%nb_ctrlpts, object%nb_qp_wq))

            ! Get basis and weights 
            call wq_get_basis_weights_generalized(object%degree, object%size_kv, & 
            object%knotvector, object%nb_qp_wq, maxrule, B0, B1, W00, W01, W10, W11)

            ! Set size of properties
            allocate(object%data_B0(object%nnz_B))
            allocate(object%data_B1(object%nnz_B))
            allocate(object%data_W00(object%nnz_B))
            allocate(object%data_W01(object%nnz_B))
            allocate(object%data_W10(object%nnz_B))
            allocate(object%data_W11(object%nnz_B))
            allocate(object%data_ind(object%nnz_B, 2))

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

                    object%data_ind(count, :) = [i, j]
                end do
            end do
            
        else    
            ! Get model 
            ! --------------
            ! Initialize
            nbel_m = degree + 3
            size_kv_m = 2*degree + multiplicity*(nbel_m - 1) + 2

            allocate(nodes_m(size_kv_m+1), knotvector_m(size_kv_m))
            call create_knotvector(degree, nbel_m, multiplicity, nodes_m, knotvector)
            call wq_set_properties(degree, size_kv_m, maxrule, dummy1, nb_ctrlpts_m, nb_qp_wq_m, dummy2)

            ! Allocate
            allocate(B0_model(degree+2, nb_qp_wq_m))
            allocate(B1_model(degree+2, nb_qp_wq_m))
            allocate(W00_model(degree+2, nb_qp_wq_m))
            allocate(W01_model(degree+2, nb_qp_wq_m))
            allocate(W10_model(degree+2, nb_qp_wq_m))
            allocate(W11_model(degree+2, nb_qp_wq_m))
            allocate(B_shape_model(nb_ctrlpts_m, 2))
            allocate(B_shape_model_temp(nb_ctrlpts_m, 2))

            call wq_get_basis_weights_model(degree, nbel, nbel_m,  size_kv_m, knotvector_m, nb_qp_wq_m, maxrule, & 
                                        B0_model, B1_model, W00_model, W01_model, W10_model, W11_model)
            call wq_get_B0_B1_shape(degree, size_kv_m, nodes_m, knotvector_m, &
                                        maxrule, B_shape_model_temp, B_shape_model)
            deallocate(B_shape_model_temp)
                            
            ! Transfer data 
            ! ----------------
            ! Allocate 
            allocate(object%data_B0(object%nnz_B))
            allocate(object%data_B1(object%nnz_B))
            allocate(object%data_W00(object%nnz_B))
            allocate(object%data_W01(object%nnz_B))
            allocate(object%data_W10(object%nnz_B))
            allocate(object%data_W11(object%nnz_B))
            allocate(object%data_ind(object%nnz_B, 2))

            ! Set p + 1 first functions
            count = 0
            do i = 1, degree + 1
                ! For B0-type
                do j = object%Bshape(i, 1), object%Bshape(i, 2)
                    count = count + 1
                    object%data_B0(count) = B0_model(i, j)
                    object%data_W00(count) = W00_model(i, j)
                    object%data_W01(count) = W01_model(i, j)

                    object%data_B1(count) = B1_model(i, j)
                    object%data_W10(count) = W10_model(i, j)
                    object%data_W11(count) = W11_model(i, j)

                    object%data_ind(count, :) = [i, j]
                end do
            end do

            ! Set repeated functions 
            do i = degree+2, object%nb_ctrlpts-degree-1
                ! For B0-type
                do j = 1, B_shape_model(degree+2, 2) - B_shape_model(degree+2, 1) + 1 
                    count = count + 1
                    object%data_B0(count) = B0_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
                    object%data_W00(count) = W00_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
                    object%data_W01(count) = W01_model(degree+2, B_shape_model(degree+2, 1) + j - 1)

                    object%data_B1(count) = B1_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
                    object%data_W10(count) = W10_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
                    object%data_W11(count) = W11_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
                    
                    object%data_ind(count, :) = [i, object%Bshape(i, 1) + j - 1]
                end do
            end do

            ! Set p + 1 last functions
            do i = degree + 1, 1, -1 
                ! For B0-type and B1-type
                do j = object%Bshape(i, 2), object%Bshape(i, 1), -1
                    count = count + 1
                    object%data_B0(count) = B0_model(i, j)
                    object%data_W00(count) = W00_model(i, j)
                    object%data_W01(count) = W01_model(i, j)
                    
                    object%data_B1(count) = -B1_model(i, j)
                    object%data_W10(count) = -W10_model(i, j)
                    object%data_W11(count) = -W11_model(i, j)

                    object%data_ind(count, :) = [object%nb_ctrlpts - i + 1, object%nb_qp_wq - j + 1]
                end do
            end do

        end if

    end subroutine wq_get_basis_weights

end module wq_basis_weights
