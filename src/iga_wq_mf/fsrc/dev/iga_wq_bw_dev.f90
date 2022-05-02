! ==========================
! module :: Basis and weights 
! author :: Joaquin Cornejo
! modules :: iga_wq_basis_weights(iga_get_basis_weights, wq_get_basis_weights, wq_initialize)
! ==========================

subroutine wq_get_size_data(degree, nb_el, nnz_B, size_qp)
    !! Gets the size of non-zeros in Basis matrix to use in wq_get_data
    
    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el
    integer, intent(out) :: nnz_B, size_qp

    ! Local data
    ! -----------------
    type(wq), pointer :: object

    call wq_initialize(object, degree, nb_el)
    nnz_B = object%nnz_B
    size_qp = object%nb_qp_wq

end subroutine wq_get_size_data

subroutine wq_get_data( degree, nb_el, nnz_B, size_qp, qp_pos, &
                        data_B0, data_B1, data_W00, &
                        data_W11, data_ind, nnz_I)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, nnz_B, size_qp

    double precision, intent(out) :: qp_pos
    dimension :: qp_pos(size_qp)
    double precision, intent(out) :: data_B0, data_B1, data_W00, data_W11
    dimension ::    data_B0(nnz_B), data_B1(nnz_B), &
                    data_W00(nnz_B), data_W11(nnz_B)
    
    integer, intent(out) :: data_ind
    dimension :: data_ind(nnz_B, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -----------------
    type(wq), pointer :: object

    ! Evaluate basis and weights
    call wq_get_basis_weights(object, degree, nb_el)

    ! Set quadrature points
    qp_pos = object%qp_pos

    ! Set data 
    ! Hypothesis : size_data_guessed >= size_data
    data_B0 = object%data_B0
    data_B1 = object%data_B1
    data_W00 = object%data_W00
    data_W11 = object%data_W11

    ! Set indexes
    data_ind = object%data_ind

    ! Set number of non zeros of integral matrix
    nnz_I = object%nnz_I

end subroutine wq_get_data

subroutine wq_get_data_csr( degree, nb_el, nnz_B, size_qp, qp_pos, &
                            data_B0, data_B1, data_W00, &
                            data_W11, indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, nnz_B, size_qp

    double precision, intent(out) :: qp_pos
    dimension :: qp_pos(size_qp)
    double precision, intent(out) :: data_B0, data_B1, data_W00, data_W11
    dimension ::    data_B0(nnz_B), data_B1(nnz_B), &
                    data_W00(nnz_B), data_W11(nnz_B)

    integer, intent(out) :: indi, indj
    dimension :: indi(degree+nb_el+1), indj(nnz_B)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -------------
    integer :: data_ind
    dimension :: data_ind(nnz_B, 2)
    double precision, dimension(:), allocatable :: data_B0_csr
    
    call wq_get_data(degree, nb_el, nnz_B, size_qp, qp_pos, data_B0, data_B1, data_W00, &
                    data_W11, data_ind, nnz_I)

    allocate(data_B0_csr(nnz_B))
    call coo2csr(degree+nb_el, nnz_B, data_B0, data_ind(:,1), data_ind(:,2), data_B0_csr, indj, indi)
    deallocate(data_B0_csr)

end subroutine wq_get_data_csr
