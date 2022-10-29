! ==========================
! module :: Basis and weights 
! author :: Joaquin Cornejo
! 
! In this module, one can find different kind of functions to compute basis and/or weights 
! for IGA-Galerkin and IGA-WQ approaches.   
! ==========================

subroutine get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_basis, indices)
    !! Gets in COO format the basis at given knots 

    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: span_tol = 1.d-8
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: data_basis
    dimension :: data_basis(nb_knots*(degree+1), 2)
    integer, intent(out) :: indices
    dimension :: indices((degree+1)*nb_knots, 2)

    ! Local data
    ! ----------
    integer :: i, j, k, nb_ctrlpts, nbel
    integer :: functions_span, span
    dimension :: functions_span(degree+1), span(2)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: nodes, B0t, B1t
    dimension :: nodes(size_kv+1), B0t(degree+1), B1t(degree+1)

    call find_unique_array(size_kv, knotvector, nodes)
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
            data_basis(k, 1) = B0t(j)
            data_basis(k, 2) = B1t(j)
            indices(k, :) = [functions_span(j), i]                                
        end do
    end do

end subroutine get_basis_generalized

subroutine get_basis_generalized_csr(degree, size_kv, knotvector, nb_knots, knots, data_basis, indi, indj)
    !! Gets in CSR format the basis at given knots 

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)

    double precision, intent(out) :: data_basis
    dimension :: data_basis(nb_knots*(degree+1), 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(nb_knots*(degree+1))

    ! Local data
    ! ----------
    integer :: size_data
    integer, dimension(:, :), allocatable :: indices
    double precision, dimension(:, :), allocatable :: data_B_coo

    size_data = nb_knots*(degree+1)
    allocate(data_B_coo(size_data, 2), indices(size_data, 2))
    call get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_B_coo, indices)

    ! Convert COO to CSR format
    call coo2csr(2, size_kv-degree-1, size_data, data_B_coo, indices(:, 1), indices(:, 2), data_basis, indj, indi)
    deallocate(data_B_coo)

end subroutine get_basis_generalized_csr

subroutine iga_get_data(degree, size_kv, knotvector, nb_qp, qp_position, qp_weight, &
                        size_data, data_basis, indices, nnz_I)
    !! Gets in COO format basis and weights in IGA approach

    use iga_basis_weights
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_qp, size_data
    double precision, intent(in) :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_position, qp_weight, data_basis
    dimension :: qp_position(nb_qp), qp_weight(nb_qp), data_basis(size_data, 2)
    integer, intent(out) :: indices
    dimension :: indices(size_data, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! ----------
    type(iga), pointer :: obj

    call iga_initialize(obj, degree, size_kv, knotvector)    
    call iga_basis_weights_dense2coo(obj)

    ! Save information
    qp_position = obj%qp_position
    qp_weight = obj%qp_weight
    data_basis(:, 1) = obj%data_B0
    data_basis(:, 2) = obj%data_B1
    indices = obj%indices
    nnz_I = obj%nnz_I

end subroutine iga_get_data

subroutine iga_get_data_csr(degree, size_kv, knotvector, nb_qp, qp_position, qp_weight, &
                            size_data, data_basis, indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA approach

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_qp, size_data
    double precision, intent(in) :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_position, qp_weight, data_basis
    dimension :: qp_position(nb_qp), qp_weight(nb_qp), data_basis(size_data, 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(size_data)
    integer, intent(out) :: nnz_I

    ! Local data
    ! ----------
    integer :: data_ind
    dimension :: data_ind(size_data, 2)
    double precision, dimension(:,:), allocatable :: data_dummy
    
    call iga_get_data(degree, size_kv, knotvector, nb_qp, qp_position, qp_weight, size_data, data_basis, data_ind, nnz_I)

    ! Convert COO to CSR format
    allocate(data_dummy(size_data, 2))
    call coo2csr(2, size_kv-degree-1, size_data, data_basis, data_ind(:,1), data_ind(:,2), data_dummy, indj, indi)
    deallocate(data_dummy)

end subroutine iga_get_data_csr

subroutine wq_get_size_data(degree, size_kv, knotvector, size_data, nb_qp)
    !! Gets the size of non-zeros of matrices that will be used in wq_get_data
    
    use wq_basis_weights
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv
    double precision, intent(in) :: knotvector
    dimension :: knotvector(size_kv)
    integer, intent(out) :: size_data, nb_qp

    ! Local data
    ! ----------
    type(wq), pointer :: obj
    integer, parameter :: method = 1 ! Method 2 is possible
    logical :: ismaxregular, isuniform

    call wq_initialize(obj, degree, size_kv, knotvector, method)
    call verify_regularity_uniformity(degree, size_kv, knotvector, ismaxregular, isuniform, span_tol)
    call wq_get_qp_positions(obj)
    call wq_get_B_shape(obj)

    ! Save information
    size_data = obj%nnz_B
    nb_qp = obj%nb_qp_wq

end subroutine wq_get_size_data

subroutine wq_get_data(degree, size_kv, knotvector, size_data, nb_qp, qp_position, &
                        data_basis, data_weights, indices, nnz_I)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, size_data, nb_qp
    double precision, intent(in) :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_position, data_basis, data_weights
    dimension :: qp_position(nb_qp), data_basis(size_data, 2), data_weights(size_data, 4)
    integer, intent(out) :: indices
    dimension :: indices(size_data, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! ----------
    type(wq), pointer :: obj
    integer, parameter :: method = 1 ! Method 2 is possible

    call wq_initialize(obj, degree, size_kv, knotvector, method)
    call wq_basis_weights_dense2coo(obj)  

    ! Save information 
    qp_position = obj%qp_position
    data_basis(:, 1) = obj%data_B0
    data_basis(:, 2) = obj%data_B1
    data_weights(:, 1) = obj%data_W00
    data_weights(:, 2) = obj%data_W01
    data_weights(:, 3) = obj%data_W10
    data_weights(:, 4) = obj%data_W11
    indices = obj%indices
    nnz_I = obj%nnz_I

end subroutine wq_get_data

subroutine wq_get_data_csr(degree, size_kv, knotvector, size_data, nb_qp, qp_position, &
                            data_basis, data_weights, indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA-WQ approach

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, size_data, nb_qp
    double precision :: knotvector
    dimension :: knotvector(size_kv)

    double precision, intent(out) :: qp_position, data_basis, data_weights
    dimension :: qp_position(nb_qp), data_basis(size_data, 2), data_weights(size_data, 4)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(size_data)
    integer, intent(out) :: nnz_I

    ! Local data
    ! ----------
    integer :: data_ind
    dimension :: data_ind(size_data, 2)
    double precision, dimension(:,:), allocatable :: data_dummy
    
    call wq_get_data(degree, size_kv, knotvector, size_data, nb_qp, qp_position, data_basis, data_weights, data_ind, nnz_I)

    ! Convert COO to CSR format
    allocate(data_dummy(size_data, 2))
    call coo2csr(2, size_kv-degree-1, size_data, data_basis, data_ind(:,1), data_ind(:,2), data_dummy, indj, indi)
    deallocate(data_dummy)

end subroutine wq_get_data_csr
