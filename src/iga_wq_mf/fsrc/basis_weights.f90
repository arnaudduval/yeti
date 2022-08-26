! ==========================
! module :: Basis and weights 
! author :: Joaquin Cornejo
! 
! In this module, one can find different kind of functions to compute basis and/or weights 
! for IGA-Galerkin and IGA-WQ approaches.   
! ==========================

subroutine get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_B, indices)
    !! Gets in COO format the basis at given knots 

    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: span_tol = 1.d-8
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: data_B
    dimension :: data_B(nb_knots*(degree+1), 2)
    integer, intent(out) :: indices
    dimension :: indices((degree+1)*nb_knots, 2)

    ! Local data
    ! -------------
    integer :: i, j, k, nb_ctrlpts, nbel
    integer :: functions_span, span
    dimension :: functions_span(degree+1), span(2)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: nodes, B0t, B1t
    dimension :: nodes(size_kv+1), B0t(degree+1), B1t(degree+1)

    ! Get nodes 
    call find_unique_vector(size_kv, knotvector, nodes)

    ! Set number of control points
    nb_ctrlpts = size_kv - degree - 1

    ! Get non repeated knot vector
    nbel = int(nodes(size_kv+1)) - 1
    allocate(table_functions_span(nbel, degree+1))
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span, span_tol)

    ! Evaluate B-spline values for every knot 
    do i = 1, nb_knots
        ! Find knot-vector span
        call find_knotvector_span(degree, size_kv, knotvector, knots(i), span(1), span_tol)

        ! Find parametric span
        call find_parametric_span(size_kv, nodes, knots(i), span(2), span_tol)

        ! Find functions over the span 
        functions_span = table_functions_span(span(2), :)

        ! Evaluate B0 and B1
        call dersbasisfuns(span(1), degree, nb_ctrlpts, knots(i), knotvector, B0t, B1t)
        
        ! Assign values
        do j = 1, degree+1
            k = (i - 1)*(degree + 1) + j
            data_B(k, 1) = B0t(j)
            data_B(k, 2) = B1t(j)
            indices(k, :) = [functions_span(j), i]                                
        end do
    end do

end subroutine get_basis_generalized

subroutine get_basis_generalized_csr(degree, size_kv, knotvector, nb_knots, knots, data_B, indi, indj)
    !! Gets in CSR format the basis at given knots 

    implicit none
    ! Input / output data
    ! ---------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)

    double precision, intent(out) :: data_B
    dimension :: data_B(nb_knots*(degree+1), 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(nb_knots*(degree+1))

    ! Local data
    ! ---------------
    integer :: size_data
    integer, dimension(:, :), allocatable :: indices
    double precision, dimension(:, :), allocatable :: data_B_coo

    ! Get results in coo format
    size_data = nb_knots*(degree+1)
    allocate(data_B_coo(size_data, 2), indices(size_data, 2))
    call get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_B_coo, indices)

    ! Get CSR format
    call coo2csr(2, size_kv-degree-1, size_data, data_B_coo, indices(:, 1), indices(:, 2), data_B, indj, indi)
    deallocate(data_B_coo)

end subroutine 

! ==============================

subroutine iga_get_data(degree, size_kv, knotvector, nnz_qp, qp_pos, qp_wq, &
                        size_data, data_B, data_ind, nnz_I)
    !! Gets in COO format basis and weights in IGA approach

    use iga_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv, nnz_qp, size_data
    double precision :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_pos, qp_wq, data_B
    dimension :: qp_pos(nnz_qp), qp_wq(nnz_qp), data_B(size_data, 2)
    integer, intent(out) :: data_ind
    dimension :: data_ind(size_data, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -----------------
    type(iga), pointer :: obj

    ! Evaluate basis and weights
    call iga_initialize(obj, degree, size_kv, knotvector)
    call iga_basis_weights_dense2coo(obj)

    ! Set quadrature points
    qp_pos = obj%qp_pos
    qp_wq = obj%qp_weights

    ! Set data 
    data_B(:, 1) = obj%data_B0
    data_B(:, 2) = obj%data_B1

    ! Set indices
    data_ind = obj%data_ind

    ! Set number of non zeros of integral matrix
    nnz_I = obj%nnz_I

end subroutine iga_get_data

subroutine iga_get_data_csr(degree, size_kv, knotvector, nnz_qp, qp_pos, qp_wq, &
                            size_data, data_B, indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA approach

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv, nnz_qp, size_data
    double precision :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_pos, qp_wq, data_B
    dimension :: qp_pos(nnz_qp), qp_wq(nnz_qp), data_B(size_data, 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(size_data)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -------------
    integer :: data_ind
    dimension :: data_ind(size_data, 2)
    double precision, dimension(:,:), allocatable :: data_dummy
    
    ! Get data in COO format
    call iga_get_data(degree, size_kv, knotvector, nnz_qp, qp_pos, qp_wq, size_data, data_B, data_ind, nnz_I)

    ! Get CSR format
    allocate(data_dummy(size_data, 2))
    call coo2csr(2, size_kv-degree-1, size_data, data_B, data_ind(:,1), data_ind(:,2), data_dummy, indj, indi)
    deallocate(data_dummy)

end subroutine iga_get_data_csr

! ==============================

subroutine wq_get_size_data(degree, size_kv, knotvector, size_data, nb_qp)
    !! Gets the size of non-zeros of matrices in wq_get_data
    
    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv
    double precision :: knotvector
    dimension :: knotvector(size_kv)
    integer, intent(out) :: size_data, nb_qp

    ! Local data
    ! -----------------
    type(wq), pointer :: obj
    integer, parameter :: method = 1
    logical :: ismaxregular, isuniform

    call wq_initialize(obj, degree, size_kv, knotvector, method)
    call verify_regularity_uniformity(degree, size_kv, knotvector, ismaxregular, isuniform, span_tol)
    call wq_get_qp_positions(obj)
    call wq_get_B0_B1_shape(obj)
    size_data = obj%nnz_B
    nb_qp = obj%nb_qp_wq

end subroutine wq_get_size_data

subroutine wq_get_data(degree, size_kv, knotvector, size_data, nb_qp, qp_pos, &
                        data_B, data_W, data_ind, nnz_I)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv, size_data, nb_qp
    double precision :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_pos, data_B, data_W
    dimension :: qp_pos(nb_qp), data_B(size_data, 2), data_W(size_data, 4)
    integer, intent(out) :: data_ind
    dimension :: data_ind(size_data, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -----------------
    type(wq), pointer :: obj
    integer, parameter :: method = 1

    ! Evaluate basis and weights
    call wq_initialize(obj, degree, size_kv, knotvector, method)
    call wq_basis_weights_dense2coo(obj, degree, size_kv, knotvector, method)

    ! Set quadrature points
    qp_pos = obj%qp_pos

    ! Set data 
    data_B(:, 1) = obj%data_B0
    data_B(:, 2) = obj%data_B1
    data_W(:, 1) = obj%data_W00
    data_W(:, 2) = obj%data_W01
    data_W(:, 3) = obj%data_W10
    data_W(:, 4) = obj%data_W11

    ! Set indices
    data_ind = obj%data_indices

    ! Set number of non zeros of integral matrix
    nnz_I = obj%nnz_I

end subroutine wq_get_data

subroutine wq_get_data_csr(degree, size_kv, knotvector, size_data, nb_qp, qp_pos, &
                            data_B, data_W, indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA-WQ approach

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, size_kv, size_data, nb_qp
    double precision :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_pos, data_B, data_W
    dimension :: qp_pos(nb_qp), data_B(size_data, 2), data_W(size_data, 4)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(size_data)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -------------
    integer :: data_ind
    dimension :: data_ind(size_data, 2)
    double precision, dimension(:,:), allocatable :: data_dummy
    
    ! Get data in COO format
    call wq_get_data(degree, size_kv, knotvector, size_data, nb_qp, qp_pos, data_B, data_W, data_ind, nnz_I)

    ! Get CSR format
    allocate(data_dummy(size_data, 2))
    call coo2csr(2, size_kv-degree-1, size_data, data_B, data_ind(:,1), data_ind(:,2), data_dummy, indj, indi)
    deallocate(data_dummy)

end subroutine wq_get_data_csr
