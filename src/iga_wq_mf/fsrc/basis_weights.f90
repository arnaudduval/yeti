! ==========================
! module :: Basis and weights 
! author :: Joaquin Cornejo
! ==========================

subroutine get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_B, indi, indj)
    !! Gets in COO format the basis at given knots 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: data_B
    dimension :: data_B(nb_knots*(degree+1), 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(nb_knots*(degree+1)), indj(nb_knots*(degree+1))

    ! Local data
    ! -------------
    integer :: i, j, k, nb_ctrlpts, para_span, kv_span, nbel
    integer :: functions_span
    dimension :: functions_span(degree+1)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: B0temp, B1temp
    dimension :: B0temp(degree+1), B1temp(degree+1)
    double precision :: nodes(size_kv+1)

    ! Set number of control points
    nb_ctrlpts = size_kv - degree - 1

    ! Get non repeated knot vector
    call find_unique_vector(size_kv, knotvector, nodes)
    nbel = int(nodes(size_kv+1)) - 1
    allocate(table_functions_span(nbel, degree+1))
    call set_table_functions_spans(degree, size_kv, nodes, knotvector, table_functions_span)

    ! Evaluate B-spline values for every knot 
    do i = 1, nb_knots
        ! Find knot-vector span
        call find_knotvector_span(degree, size_kv, knotvector, knots(i), kv_span)

        ! Find parametric span
        call find_parametric_span(size_kv, nodes, knots(i), para_span)

        ! Find functions over the span 
        functions_span = table_functions_span(para_span, :)

        ! Evaluate B0 and B1
        call dersbasisfuns(kv_span, degree, nb_ctrlpts, knots(i), knotvector, B0temp, B1temp)
        
        ! Assign values
        do j = 1, degree+1
            k = (i - 1)*(degree + 1) + j
            data_B(k, 1) = B0temp(j)
            data_B(k, 2) = B1temp(j)
            indi(k) = functions_span(j)
            indj(k) = i                               
        end do
    end do

end subroutine get_basis_generalized

subroutine get_basis_generalized_csr(degree, size_kv, knotvector, nb_knots, knots, data_B, indi, indj)
    !! Gets in COO format the basis at given knots 

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
    integer, dimension(:), allocatable :: indi_coo, indj_coo
    double precision, dimension(:, :), allocatable :: data_B_coo

    ! Get results in coo format
    size_data = nb_knots*(degree+1)
    allocate(data_B_coo(size_data, 2), indi_coo(size_data), indj_coo(size_data))
    call get_basis_generalized(degree, size_kv, knotvector, nb_knots, knots, data_B_coo, indi_coo, indj_coo)

    ! Get CSR format
    call coo2csr(2, size_kv-degree-1, size_data, data_B_coo, indi_coo, indj_coo, data_B, indj, indi)
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
    type(iga), pointer :: object

    ! Evaluate basis and weights
    call iga_get_basis_weights(object, degree, size_kv, knotvector)

    ! Set quadrature points
    qp_pos = object%qp_pos
    qp_wq = object%qp_weights

    ! Set data 
    data_B(:, 1) = object%data_B0
    data_B(:, 2) = object%data_B1

    ! Set indices
    data_ind = object%data_ind

    ! Set number of non zeros of integral matrix
    nnz_I = object%nnz_I

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
    !! Gets the size of non-zeros in Basis matrix to use in wq_get_data
    
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
    type(wq), pointer :: object

    call wq_initialize(object, degree, size_kv, knotvector)
    size_data = object%nnz_B
    nb_qp = object%nb_qp_wq

end subroutine wq_get_size_data

subroutine wq_get_data( degree, size_kv, knotvector, size_data, nb_qp, qp_pos, &
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
    type(wq), pointer :: object

    ! Evaluate basis and weights
    call wq_get_basis_weights(object, degree, size_kv, knotvector)

    ! Set quadrature points
    qp_pos = object%qp_pos

    ! Set data 
    data_B(:, 1) = object%data_B0
    data_B(:, 2) = object%data_B1
    data_W(:, 1) = object%data_W00
    data_W(:, 2) = object%data_W01
    data_W(:, 3) = object%data_W10
    data_W(:, 4) = object%data_W11

    ! Set indices
    data_ind = object%data_indices

    ! Set number of non zeros of integral matrix
    nnz_I = object%nnz_I

end subroutine wq_get_data

subroutine wq_get_data_csr( degree, size_kv, knotvector, size_data, nb_qp, qp_pos, &
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
