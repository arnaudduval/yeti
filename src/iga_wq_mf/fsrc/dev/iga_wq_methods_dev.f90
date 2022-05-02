! ==========================
! module :: iga and wq methods  ! ALOT OF CHANGES   
! author :: Joaquin Cornejo
! hypothesis : all these routines assume an open and uniform 
!              knot-vector and at least 2 elements
! modules :: algebra(linspace, solve_system, product_AWB),
!            bspline(set_table_functions_spans, get_parametric_nodes, 
!                   get_knotvector, get_basis)
!            GaussLegendre.f90 
! ==========================

subroutine wq_set_properties(degree, nb_el, multiplicity, maxrule, nb_ctrlpts, &
                            size_kv, nb_qp_wq, nb_qp_cgg)
    !! Sets constants used in WQ approach
    
    use constants_iga_wq_mf
    implicit none 
    ! Input / out data 
    ! -----------------
    integer, intent(in) :: degree, nb_el, multiplicity, maxrule
    integer, intent(out) :: nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg

    ! Find constants values : 
    ! Number of control points
    nb_ctrlpts = degree + multiplicity*(nb_el - 1) + 1 

    ! Number of knots of knot-vector
    size_kv = degree + nb_ctrlpts + 1 

    ! Number of quadrature points in WQ approach
    nb_qp_wq = 2*(degree + r) + nb_el*(maxrule + 1) - 2*maxrule - 3  

    ! Number of Gauss quadrature points 
    nb_qp_cgg = (degree + 1)*nb_el 

end subroutine wq_set_properties

subroutine wq_get_qp_positions(degree, nb_el, maxrule, nodes, nb_qp_wq, qp_pos)
    !! Gets quadrature points' positions (QP) in WQ approach 

    use constants_iga_wq_mf
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, maxrule, nb_qp_wq 
    double precision, intent(in) :: nodes

    double precision, intent(out) :: qp_pos 
    dimension ::    qp_pos(nb_qp_wq), &
                    nodes(nb_el+1)

    ! Local data
    ! ---------------
    double precision :: QPB, QPI ! Q.P. at boundaries, Q.P. at internal spans
    dimension :: QPB(degree+r), QPI(2+maxrule)
    integer :: nb_qp_B ! nof knots at boundaries, nof Q.P. 
    integer :: i, j, k

    ! Set number of points at the boundaries
    nb_qp_B = degree + r

    ! Find values for the first boundary
    call linspace(nodes(1), nodes(2), nb_qp_B, QPB)
    
    ! Assign values
    do j = 1, nb_qp_B
        qp_pos(j) = QPB(j)
    end do

    ! Find values for the last boundary 
    call linspace(nodes(nb_el), nodes(nb_el+1), nb_qp_B, QPB)
    
    ! Assign values
    do j = 1, nb_qp_B
        qp_pos(nb_qp_wq-nb_qp_B+j) = QPB(j)
    end do

    if (nb_el .ge. 3) then
        do i = 2, nb_el-1
            ! Find quadrature points for inner spans
            call linspace(nodes(i), nodes(i+1), 2+maxrule, QPI)

            ! Assign values
            do j = 1, 2 + maxrule  
                k = nb_qp_B + (1 + maxrule)*(i - 2) + j - 1    
                qp_pos(k) = QPI(j)
            end do
        end do
    end if

end subroutine wq_get_qp_positions      

subroutine iga_get_qp_positions_weights(degree, nb_el, nodes, nb_qp_cgg, qp_pos, qp_weights)
    !! Gets quadrature points' positions and weights in IGA approach 

    implicit none 
    ! Input / Output data
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_qp_cgg
    double precision, intent(in) :: nodes
    dimension :: nodes(nb_el+1)

    double precision, intent(out) :: qp_pos, qp_weights
    dimension :: qp_pos(nb_qp_cgg), qp_weights(nb_qp_cgg)

    ! Local data 
    ! -------------
    double precision :: GaussPdsCoord
    dimension :: GaussPdsCoord(2, degree+1)
    double precision :: xg(degree+1), wg(degree+1)
    integer ::  i, j, k

    ! Find position and weight in isoparametric space
    call Gauss(degree+1, 1, GaussPdsCoord, 0)

    do i = 1, degree + 1
        wg(i) = GaussPdsCoord(1, i)
        xg(i) = GaussPdsCoord(2, i)
    end do 

    do i = 1, nb_el
        do j = 1, degree + 1
            k = (i - 1)*(degree + 1) + j
            qp_pos(k) = 0.5d0*(xg(j)/dble(nb_el) + nodes(i+1) + nodes(i))
            qp_weights(k) = 0.5d0/dble(nb_el)*wg(j)
        end do
    end do

end subroutine iga_get_qp_positions_weights 

subroutine iga_get_B_shape(degree, nb_el, nb_ctrlpts, multiplicity, B_shape)
    !! Gets non zeros positions of B0 and B1 in IGA approach 
    !! (uniform knot-vector)
    
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_ctrlpts, multiplicity

    integer, intent(out) :: B_shape
    dimension :: B_shape(nb_ctrlpts, 2)

    ! Local data 
    ! --------------
    integer :: min_span, max_span, min_knot, max_knot
    integer :: table_points_span, table_functions_span, table_spans_function
    dimension ::    table_points_span(nb_el, 2), &
                    table_functions_span(nb_el, degree+1), &
                    table_spans_function(nb_ctrlpts, 2)
    integer ::  i, j

    ! Get table of points over the span
    table_points_span = 1
    table_points_span(1, 2) = degree + 1 
    
    do i = 2, nb_el 
        table_points_span(i, 1) = table_points_span(i-1, 1) + degree + 1
        table_points_span(i, 2) = table_points_span(i, 1) + degree
    end do
    
    ! Get table of functions on span 
    call set_table_functions_spans(degree, nb_el, multiplicity, table_functions_span)

    ! Get table of spans for each function
    do i = 1, nb_ctrlpts
        ! Find min 
        min_span = 1
        jloop_min : do j = 1, nb_el
            if (any(table_functions_span(j, :).eq.i)) then
                    min_span = j
                    exit jloop_min
            end if
        end do jloop_min

        ! Find max
        max_span = nb_el
        jloop_max : do j = nb_el, 1, -1
            if (any(table_functions_span(j, :).eq.i)) then
                    max_span = j
                    exit jloop_max
            end if
        end do jloop_max

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

        B_shape(i, :) = [min_knot, max_knot]
    end do

end subroutine iga_get_B_shape

subroutine wq_get_nnz_B0_B1_shape(degree, nb_el, nb_ctrlpts, nb_qp, multiplicity, maxrule, B0_shape, B1_shape)
    !! Gets non-zero positions of B0 and B1 in WQ approach
    !! uniform knot-vector

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_ctrlpts, nb_qp, multiplicity, maxrule

    integer, intent(out) :: B0_shape, B1_shape
    dimension :: B0_shape(nb_ctrlpts, 2), B1_shape(nb_ctrlpts, 2)

    ! Local data 
    ! --------------
    integer :: min_span, max_span, min_knot, max_knot
    integer :: table_points_span, table_functions_span, table_spans_function
    dimension ::    table_points_span(nb_el, 2), &
                    table_functions_span(nb_el, degree+1), &
                    table_spans_function(nb_ctrlpts, 2)
    integer ::  i, j

    ! Get table of points over the span
    table_points_span = 1
    table_points_span(1, 2) = degree + r 
    table_points_span(nb_el, 1) = nb_qp + 1 - (degree + r)
    table_points_span(nb_el, 2) = nb_qp
    
    do i = 2, nb_el - 1
        table_points_span(i, 1) = table_points_span(i-1, 2)
        table_points_span(i, 2) = table_points_span(i, 1) + 1 + maxrule
    end do

    ! Get table of functions on span 
    call set_table_functions_spans(degree, nb_el, multiplicity, table_functions_span)

    ! Get table of spans for each function
    do i = 1, nb_ctrlpts
        ! Find min 
        min_span = 1
        jloop_min : do j = 1, nb_el
            if (any(table_functions_span(j, :).eq.i)) then
                    min_span = j
                    exit jloop_min
            end if
        end do jloop_min

        ! Find max
        max_span = nb_el
        jloop_max : do j = nb_el, 1, -1
            if (any(table_functions_span(j, :).eq.i)) then
                    max_span = j
                    exit jloop_max
            end if
        end do jloop_max

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

        B0_shape(i, :) = [min_knot, max_knot]
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

        B1_shape(i, :) = [min_knot, max_knot]
    end do

end subroutine wq_get_nnz_B0_B1_shape

subroutine wq_get_quadrature_rule_row(ifunct, nb_rows, nb_ctrlpts, nb_qp, MB, MI, Bshape, MIint, irule)
    !! Gets the quadrature rule of the i-th function
    
    implicit none
    ! Input / output data
    ! -----------------------
    integer, intent(in) :: ifunct, nb_rows, nb_ctrlpts, nb_qp
    integer, intent(in) :: Bshape, MIint
    dimension :: Bshape(nb_ctrlpts, 2), MIint(nb_rows, nb_ctrlpts)
    double precision, intent(in) ::  MB, MI
    dimension :: MB(nb_rows, nb_qp), MI(nb_rows, nb_ctrlpts)

    double precision, intent(inout) ::  irule
    dimension :: irule(nb_qp)

    ! Local data
    ! ----------------
    integer ::  Pmin, Pmax, Fmin, Fmax
    integer :: j

    ! Find position of points within i-function support
    Pmin = Bshape(ifunct, 1)
    Pmax = Bshape(ifunct, 2)

    ! Find functions which intersect i-function support
    Fmin = 1
    jloop_min : do j = 1, nb_rows
        if (MIint(j, ifunct).gt.0) then
                Fmin = j
                exit jloop_min
        end if
    end do jloop_min

    ! Find max
    Fmax = nb_rows
    jloop_max : do j = nb_rows, 1, -1
        if (MIint(j, ifunct).gt.0) then
                Fmax = j
                exit jloop_max
        end if
    end do jloop_max

    call solve_system(Fmax-Fmin+1, Pmax-Pmin+1, MB(Fmin:Fmax, Pmin:Pmax), & 
                        MI(Fmin:Fmax, ifunct), irule(Pmin:Pmax))

end subroutine wq_get_quadrature_rule_row

! Fonctions for case number of elements <= degree + 3
subroutine wq_get_quadrature_rules(nb_rows, nb_ctrlpts, nb_qp, MB, MI, Bshape, MIint, rules)
    !! Computes the weight quadrature rule given B and I

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: nb_rows, nb_ctrlpts, nb_qp
    integer, intent(in) :: Bshape, MIint
    dimension ::    Bshape(nb_ctrlpts, 2), MIint(nb_rows, nb_ctrlpts)
    double precision, intent(in) ::  MB, MI
    dimension ::    MB(nb_rows, nb_qp), MI(nb_rows, nb_ctrlpts)

    double precision, intent(out) ::  rules
    dimension :: rules(nb_ctrlpts, nb_qp)

    ! Local data
    ! -------------
    integer :: i

    ! Set values of MW
    rules = 0.0d0

    ! Compute quadrature rules
    do i = 1, nb_ctrlpts
        call wq_get_quadrature_rule_row(i, nb_rows, nb_ctrlpts, nb_qp, MB, MI, &
                                        Bshape, MIint, rules(i, :))
    end do

end subroutine wq_get_quadrature_rules

subroutine wq_get_basis_weights_generalized(degree, nb_el, nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg, maxrule, &
                                            B0wq_p0, B1wq_p0, W00, W11)
    !! Returns the basis and weights data at the quadrature points in WQ approach 

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg, maxrule
    
    double precision, intent(out) :: B0wq_p0, B1wq_p0 
    dimension ::    B0wq_p0(nb_ctrlpts, nb_qp_wq), B1wq_p0(nb_ctrlpts, nb_qp_wq)
    double precision, intent(out) :: W00, W11 
    dimension ::    W00(nb_ctrlpts, nb_qp_wq), W11(nb_ctrlpts, nb_qp_wq)

    ! Local data
    ! -------------        
    ! For p continuity p-1
    integer, parameter :: multiplicity = 1
    integer :: table_functions_span_p0
    dimension :: table_functions_span_p0(nb_el, degree+1)

    double precision :: knotvector, qp_cgg_pos, qp_cgg_weights, qp_wq_pos
    dimension ::    knotvector(size_kv), qp_cgg_pos(nb_qp_cgg), &
                    qp_cgg_weights(nb_qp_cgg), qp_wq_pos(nb_qp_wq)

    double precision :: B0cgg_p0, B1cgg_p0
    dimension ::    B0cgg_p0(nb_ctrlpts, nb_qp_cgg), & 
                    B1cgg_p0(nb_ctrlpts, nb_qp_cgg)

    integer :: Bint_p0
    dimension ::  Bint_p0(nb_ctrlpts, nb_qp_cgg)

    double precision :: nodes
    dimension :: nodes(nb_el+1)

    ! For p continuity p-2
    integer :: degree_p1, nb_el_p1, multiplicity_p1
    integer :: nb_ctrlpts_p1, nb_knots_p1, nb_qp_wq_p1, nb_qp_cgg_p1

    integer :: table_functions_span_p1
    dimension :: table_functions_span_p1(nb_el, degree+1) 
    double precision, dimension(:), allocatable :: knotvector_p1 

    double precision, dimension(:,:), allocatable :: B0cgg_p1, B1cgg_p1
    integer, dimension(:,:), allocatable :: Bint_p1
    double precision, allocatable, dimension(:,:) :: B0wq_p1, B1wq_p1

    ! Integrals
    double precision, dimension(:,:), allocatable :: I00, I11

    ! Weights
    integer :: B0shape
    dimension :: B0shape(nb_ctrlpts, 2)
    integer, allocatable, dimension(:, :) :: MIint, B1shape
    
    ! Find knots non-repeated of knot-vector
    call get_parametric_nodes(nb_el, nodes)

    ! Define degree of test functions
    degree_p1 = degree
    nb_el_p1 = nb_el
    multiplicity_p1 = 2

    ! Initialize some properties of test functions
    call wq_set_properties(degree_p1, nb_el_p1, multiplicity_p1, maxrule, nb_ctrlpts_p1, nb_knots_p1, &
                                    nb_qp_wq_p1, nb_qp_cgg_p1)

    ! Get non-zero values shape
    allocate(B1shape(nb_ctrlpts, 2))
    call wq_get_nnz_B0_B1_shape(degree, nb_el, nb_ctrlpts, nb_qp_wq, multiplicity, maxrule, B0shape, B1shape)
    deallocate(B1shape)

    ! --------------------------
    ! Degree p
    ! --------------------------

    ! Find knot-vector in WQ approach
    call get_knotvector(degree, nb_el, nodes, size_kv, knotvector, multiplicity)

    ! Find table of functions over span
    call set_table_functions_spans(degree, nb_el, multiplicity, table_functions_span_p0)

    ! Find positions and weights in IGA approach
    call iga_get_qp_positions_weights(degree, nb_el, nodes, nb_qp_cgg, qp_cgg_pos, qp_cgg_weights)

    ! Find basis at Gauss quadrature points
    call get_basis(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, nb_qp_cgg, qp_cgg_pos, & 
                    table_functions_span_p0, B0cgg_p0, B1cgg_p0)

    ! Find quadrature points in WQ approach
    call wq_get_qp_positions(degree, nb_el, maxrule, nodes, nb_qp_wq, qp_wq_pos) 

    ! Find basis at WQ quadrature points
    call get_basis(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, nb_qp_wq, qp_wq_pos, & 
                    table_functions_span_p0, B0wq_p0, B1wq_p0) 

    ! ! --------------------------
    ! ! Degree p - 1 
    ! ! -------------------------- 

    ! Find knot-vector in WQ approach
    allocate(knotvector_p1(nb_knots_p1))
    call get_knotvector(degree_p1, nb_el_p1, nodes, nb_knots_p1, knotvector_p1, multiplicity_p1)

    ! Find table of functions on span
    call set_table_functions_spans(degree_p1, nb_el_p1, multiplicity_p1, table_functions_span_p1)

    ! Find basis function values at Gauss points
    allocate(B0cgg_p1(nb_ctrlpts_p1, nb_qp_cgg), B1cgg_p1(nb_ctrlpts_p1, nb_qp_cgg))
    call get_basis(degree_p1, nb_el_p1, nodes, nb_ctrlpts_p1, nb_knots_p1, knotvector_p1, &
                    nb_qp_cgg, qp_cgg_pos, table_functions_span_p1, B0cgg_p1, B1cgg_p1) 
    deallocate(B1cgg_p1)

    ! Find basis function values at WQ points
    allocate(B0wq_p1(nb_ctrlpts_p1, nb_qp_wq), B1wq_p1(nb_ctrlpts_p1, nb_qp_wq))
    call get_basis(degree_p1, nb_el_p1, nodes, nb_ctrlpts_p1, nb_knots_p1, knotvector_p1, &
                    nb_qp_wq, qp_wq_pos, table_functions_span_p1, B0wq_p1, B1wq_p1) 
    deallocate(B1wq_p1)

    ! ------------------------------------
    ! Integrals and Weights
    ! ------------------------------------
    allocate(Bint_p1(nb_ctrlpts_p1, nb_qp_cgg))
    ! Initialiaze
    Bint_p0 = 0
    Bint_p1 = 0

    where (abs(B0cgg_p0).GT.tol)
        Bint_p0 = 1
    end where

    where (abs(B0cgg_p1).GT.tol)
        Bint_p1 = 1
    end where

    allocate(MIint(nb_ctrlpts_p1, nb_ctrlpts))
    MIint = matmul(Bint_p1, transpose(Bint_p0))

    ! ----------------------
    ! I00 = B0cgg_p1 * Wcgg * B0cgg_p0.transpose
    allocate(I00(nb_ctrlpts_p1, nb_ctrlpts))
    call product_AWB(nb_ctrlpts_p1, nb_ctrlpts, nb_qp_cgg, &
                            B0cgg_p1, qp_cgg_weights, B0cgg_p0, I00)
    
    ! For W00
    call wq_get_quadrature_rules(nb_ctrlpts_p1, nb_ctrlpts, nb_qp_wq, B0wq_p1, I00, B0shape, MIint, W00)

    ! ----------------------
    ! I11 = B0cgg_p1 * Wcgg * B1cgg_p0.transpose
    allocate(I11(nb_ctrlpts_p1, nb_ctrlpts))
    call product_AWB(nb_ctrlpts_p1, nb_ctrlpts, nb_qp_cgg, &
                            B0cgg_p1, qp_cgg_weights, B1cgg_p0, I11)
                    
    ! For W11
    call wq_get_quadrature_rules(nb_ctrlpts_p1, nb_ctrlpts, nb_qp_wq, B0wq_p1, I11, B0shape, MIint, W11)               
    deallocate(MIint)

end subroutine wq_get_basis_weights_generalized

! Fonctions for case number of elements > degree + 3
subroutine wq_get_basis_weights_model(degree_model, nb_el, maxrule, nb_qp_wp_model, B0_model, B1_model, &
                                W00_model, W11_model)
    !! Gets model where nb_el = degree + 3, nb_qp_wq = 2*(degree + r) + (degree + 3)*(maxrule + 1) - 2*maxrule - 3  
    
    use constants_iga_wq_mf
    implicit none
    ! Input / output data 
    ! --------------------
    integer, intent(in) :: degree_model, nb_el, maxrule, nb_qp_wp_model

    double precision, intent(out) :: B0_model, B1_model 
    dimension ::    B0_model(degree_model+2, nb_qp_wp_model), & 
                    B1_model(degree_model+2, nb_qp_wp_model)
    double precision, intent(out) :: W00_model, W11_model
    dimension ::    W00_model(degree_model+2, nb_qp_wp_model), &
                    W11_model(degree_model+2, nb_qp_wp_model)
    
    ! Local data 
    ! ---------------
    integer :: nb_el_m, nb_ctrlpts_m, size_kv_m, nb_qp_wq_m, nb_qp_cgg_m, multiplicity
    double precision, allocatable, dimension(:,:) :: B0_mt, B1_mt, W00_mt, W01_mt, W10_mt, W11_mt
    integer :: i
    
    ! Set properties
    multiplicity = 1
    nb_el_m = degree_model + 3
    call wq_set_properties(degree_model, nb_el_m, multiplicity, maxrule, nb_ctrlpts_m, &
                            size_kv_m, nb_qp_wq_m, nb_qp_cgg_m)

    allocate(B0_mt(nb_ctrlpts_m, nb_qp_wq_m), B1_mt(nb_ctrlpts_m, nb_qp_wq_m), &
            W00_mt(nb_ctrlpts_m, nb_qp_wq_m), W01_mt(nb_ctrlpts_m, nb_qp_wq_m), &
            W10_mt(nb_ctrlpts_m, nb_qp_wq_m), W11_mt(nb_ctrlpts_m, nb_qp_wq_m))

    call wq_get_basis_weights_generalized(degree_model, nb_el_m, nb_ctrlpts_m, size_kv_m, &
    nb_qp_wq_m, nb_qp_cgg_m, maxrule, B0_mt, B1_mt, W00_mt, W11_mt)

    do i = 1, degree_model+2
        B0_model(i, :) = B0_mt(i, :)
        B1_model(i, :) = B1_mt(i, :) * nb_el / nb_el_m
        W00_model(i, :) = W00_mt(i, :) * nb_el_m / nb_el
        W11_model(i, :) = W11_mt(i, :)
    end do

    deallocate(B0_mt, B1_mt, W00_mt, W01_mt, W10_mt, W11_mt)

end subroutine wq_get_basis_weights_model
module wq_basis_weights
   
    implicit none
    ! Parameters
    integer, parameter ::  multiplicity=1, maxrule=2
    type :: wq
        ! Inputs :
        ! ------------
        integer :: degree, nb_el, nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg
        
        ! Outputs :
        ! ---------
        integer :: nnz_B, nnz_I
        double precision, dimension(:), pointer :: qp_pos
        double precision, dimension(:), pointer ::  data_B0, data_B1, data_W00, data_W11
        integer, dimension(:, :), pointer :: data_ind
        
        ! Local :
        ! ---------
        integer, dimension(:, :), pointer :: Bshape

    end type wq

    contains

    subroutine wq_initialize(object, degree, nb_el)
            
        implicit none
        ! Input / output
        ! ---------------
        integer, intent(in) :: degree, nb_el
        type(wq), pointer :: object

        ! Local data
        ! --------------
        double precision :: nodes(nb_el+1)
        integer, allocatable, dimension(:, :) :: Bshape_temp
        integer :: count, i

        ! Set properties
        allocate(object)      
        object%degree = degree
        object%nb_el = nb_el
        call wq_set_properties(object%degree, object%nb_el, multiplicity, maxrule, object%nb_ctrlpts, &
                            object%size_kv, object%nb_qp_wq, object%nb_qp_cgg)
        
        ! Get non-repeated knots
        call get_parametric_nodes(nb_el, nodes)
        
        ! Get quadrature points position
        allocate(object%qp_pos(object%nb_qp_wq))
        call wq_get_qp_positions(degree, nb_el, maxrule, nodes, object%nb_qp_wq, object%qp_pos)

        ! Set B0_shape and B1_shape
        allocate(object%Bshape(object%nb_ctrlpts, 2))
        allocate(Bshape_temp(object%nb_ctrlpts, 2))
        call wq_get_nnz_B0_B1_shape(degree, nb_el, object%nb_ctrlpts, object%nb_qp_wq, multiplicity, maxrule, &
                                Bshape_temp, object%Bshape)  
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

    subroutine wq_get_basis_weights(object, degree, nb_el)
            
        implicit none 
        ! Input / output data
        ! --------------------
        integer, intent(in) :: degree, nb_el
        type(wq), pointer :: object
        
        ! Local data
        ! ----------
        ! Model 
        integer :: nb_el_model, nb_ctrlpts_model, & 
                nb_knots_model, nb_qp_wq_model, nb_qp_cgg_model

        ! Variables for case nb_el <= p + 3
        double precision, dimension(:,:), allocatable :: B0, B1, W00, W11

        ! Variables for case nb_el > p + 3
        double precision, dimension(:,:), allocatable :: B0_model, B1_model
        double precision, dimension(:,:), allocatable :: W00_model, W11_model
        integer, dimension(:,:), allocatable :: B_shape_model, B_shape_model_temp

        ! Loops
        integer :: i, j, count

        ! Create object
        call wq_initialize(object, degree, nb_el)

        if (nb_el.le.degree + 3) then 

            ! Allocate variables
            allocate(B0(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(B1(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W00(object%nb_ctrlpts, object%nb_qp_wq))
            allocate(W11(object%nb_ctrlpts, object%nb_qp_wq))

            ! Get basis and weights 
            call wq_get_basis_weights_generalized(degree, nb_el, object%nb_ctrlpts, object%size_kv, &
            object%nb_qp_wq, object%nb_qp_cgg, maxrule, B0, B1, W00, W11)

            ! Set size of properties
            allocate(object%data_B0(object%nnz_B))
            allocate(object%data_B1(object%nnz_B))
            allocate(object%data_W00(object%nnz_B))
            allocate(object%data_W11(object%nnz_B))
            allocate(object%data_ind(object%nnz_B, 2))

            ! Assign values
            count = 0
            do i = 1, object%nb_ctrlpts
                do j = object%Bshape(i, 1), object%Bshape(i, 2)
                    count = count + 1
                    object%data_B0(count) = B0(i, j)
                    object%data_W00(count) = W00(i, j)

                    object%data_B1(count) = B1(i, j)
                    object%data_W11(count) = W11(i, j)

                    object%data_ind(count, :) = [i, j]
                end do
            end do
            
        else    
            ! Get model 
            ! --------------
            ! Number of elements of the model
            nb_el_model = degree + 3

            ! Initialize
            call wq_set_properties(degree, nb_el_model, multiplicity, maxrule, nb_ctrlpts_model, & 
                                    nb_knots_model, nb_qp_wq_model, nb_qp_cgg_model)

            ! Allocate
            allocate(B0_model(degree+2, nb_qp_wq_model))
            allocate(B1_model(degree+2, nb_qp_wq_model))
            allocate(W00_model(degree+2, nb_qp_wq_model))
            allocate(W11_model(degree+2, nb_qp_wq_model))
            allocate(B_shape_model(nb_ctrlpts_model, 2))
            allocate(B_shape_model_temp(nb_ctrlpts_model, 2))

            call wq_get_basis_weights_model(degree, nb_el, maxrule, nb_qp_wq_model, B0_model, B1_model, & 
                                            W00_model, W11_model)
            call wq_get_nnz_B0_B1_shape(degree, nb_el_model, nb_ctrlpts_model, nb_qp_wq_model, &
                                multiplicity, maxrule, B_shape_model_temp, B_shape_model)
            deallocate(B_shape_model_temp)
                            
            ! Transfer data 
            ! ----------------
            ! Allocate 
            allocate(object%data_B0(object%nnz_B))
            allocate(object%data_B1(object%nnz_B))
            allocate(object%data_W00(object%nnz_B))
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

                    object%data_B1(count) = B1_model(i, j)
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

                    object%data_B1(count) = B1_model(degree+2, B_shape_model(degree+2, 1) + j - 1)
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
                    
                    object%data_B1(count) = -B1_model(i, j)
                    object%data_W11(count) = -W11_model(i, j)

                    object%data_ind(count, :) = [object%nb_ctrlpts - i + 1, object%nb_qp_wq - j + 1]
                end do
            end do

        end if

    end subroutine wq_get_basis_weights

end module wq_basis_weights
