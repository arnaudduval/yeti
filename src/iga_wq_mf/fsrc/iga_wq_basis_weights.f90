! ==========================
! module :: Basis and weights 
! author :: Joaquin Cornejo
! modules :: iga_wq_basis_weights(iga_get_basis_weights, wq_get_basis_weights, wq_initialize)
! ==========================

subroutine get_basis_generalized(degree, nb_el, nb_knots, knots, data_B0, data_B1, indi, indj)
    !! Gets in COO format the basis at given knots 

    implicit none

    ! Input / output data
    ! ---------------------
    !f2py intent(in) :: degree, nb_el, nb_knots
    !f2py depend(nb_knots) :: knots
    integer :: degree, nb_el, nb_knots
    double precision :: knots
    dimension :: knots(nb_knots)

    double precision, intent(out) :: data_B0, data_B1
    dimension :: data_B0(nb_knots*(degree+1)), data_B1(nb_knots*(degree+1))
    integer, intent(out) ::  indi, indj
    dimension :: indi(nb_knots*(degree+1)), indj(nb_knots*(degree+1))

    ! Local data
    ! ------------
    integer :: i, j, k
    integer :: nb_ctrlpts, size_kv
    integer:: table_functions_span
    dimension :: table_functions_span(nb_el, degree+1)
    double precision :: nodes, knotvector
    dimension :: nodes(nb_el+1), knotvector(2*degree+nb_el+1)
    integer :: table_spans
    dimension :: table_spans(nb_knots, 2)

    double precision :: B0temp, B1temp
    integer :: span, functions_span
    dimension ::    functions_span(degree+1), &
                    B0temp(degree+1), B1temp(degree+1)

    ! External data
    ! -------------
    ! Assign values
    nb_ctrlpts = degree + nb_el
    size_kv = 2*degree + nb_el + 1
                
    ! Get non-repeated knots
    call get_parametric_nodes(nb_el, nodes)

    ! Get knot-vector
    call get_knotvector(degree, nb_el, nodes, size_kv, knotvector, 1)

    ! Get table
    call set_table_functions_spans(degree, nb_el, table_functions_span)

    ! Set table of spans for every knot
    call set_table_spans(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, & 
                            nb_knots, knots, table_spans)

    ! Evaluate B-spline values for every knot 
    do i = 1, nb_knots
        ! Find parametric span in knot-vector
        span = table_spans(i, 2)

        ! Find functions over the span 
        functions_span = table_functions_span(table_spans(i, 1), :)

        ! Evaluate B0 and B1
        call dersbasisfuns(span, degree, nb_ctrlpts, knots(i), knotvector, B0temp, B1temp)
        
        ! Assign values
        do j = 1, degree+1
            k = (i - 1)*(degree + 1) + j
            data_B0(k) = B0temp(j)
            data_B1(k) = B1temp(j)
            indi(k) = functions_span(j)
            indj(k) = i
        end do
    end do

end subroutine get_basis_generalized

subroutine iga_get_data( degree, nb_el, &
                        nnz_B, qp_pos, qp_weights, &
                        data_B0, data_B1, data_ind, nnz_I)
    !! Gets in COO format basis and weights in IGA approach

    use iga_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, nnz_B

    double precision, intent(out) :: qp_pos, qp_weights
    dimension :: qp_pos(nb_el*(degree+1)), qp_weights(nb_el*(degree+1))
    double precision, intent(out) :: data_B0, data_B1
    dimension :: data_B0(nnz_B), data_B1(nnz_B)
    integer, intent(out) :: data_ind
    dimension :: data_ind(nnz_B, 2)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -----------------
    type(iga), pointer :: object

    ! Evaluate basis and weights
    call iga_get_basis_weights(object, degree, nb_el)

    ! Set quadrature points
    qp_pos = object%qp_pos
    qp_weights = object%qp_weights

    ! Set data 
    data_B0 = object%data_B0
    data_B1 = object%data_B1

    ! Set indexes
    data_ind = object%data_ind

    ! Set number of non zeros of integral matrix
    nnz_I = object%nnz_I

end subroutine iga_get_data

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
                        data_W01, data_W10, data_W11, &
                        data_ind, nnz_I)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, nnz_B, size_qp

    double precision, intent(out) :: qp_pos
    dimension :: qp_pos(size_qp)
    double precision, intent(out) :: data_B0, data_B1, data_W00, &
                                    data_W01, data_W10, data_W11
    dimension ::    data_B0(nnz_B), data_B1(nnz_B), &
                    data_W00(nnz_B), data_W01(nnz_B), &
                    data_W10(nnz_B), data_W11(nnz_B)
    
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
    data_W01 = object%data_W01
    data_W10 = object%data_W10
    data_W11 = object%data_W11

    ! Set indexes
    data_ind = object%data_ind

    ! Set number of non zeros of integral matrix
    nnz_I = object%nnz_I

end subroutine wq_get_data

subroutine wq_get_data_csr( degree, nb_el, nnz_B, size_qp, qp_pos, &
                            data_B0, data_B1, data_W00, &
                            data_W01, data_W10, data_W11, &
                            indi, indj, nnz_I)
    !! Gets in CSR format basis and weights in IGA-WQ approach

    use wq_basis_weights
    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in) :: degree, nb_el, nnz_B, size_qp

    double precision, intent(out) :: qp_pos
    dimension :: qp_pos(size_qp)
    double precision, intent(out) :: data_B0, data_B1, data_W00, &
                    data_W01, data_W10, data_W11
    dimension ::    data_B0(nnz_B), data_B1(nnz_B), &
                    data_W00(nnz_B), data_W01(nnz_B), &
                    data_W10(nnz_B), data_W11(nnz_B)

    integer, intent(out) :: indi, indj
    dimension :: indi(degree+nb_el+1), indj(nnz_B)
    integer, intent(out) :: nnz_I

    ! Local data
    ! -------------
    integer :: data_ind
    dimension :: data_ind(nnz_B, 2)
    double precision, dimension(:), allocatable :: data_B0_csr
    
    call wq_get_data(degree, nb_el, nnz_B, size_qp, qp_pos, data_B0, data_B1, data_W00, &
                    data_W01, data_W10, data_W11, data_ind, nnz_I)

    allocate(data_B0_csr(nnz_B))
    call coo2csr(degree+nb_el, nnz_B, data_B0, data_ind(:,1), data_ind(:,2), data_B0_csr, indj, indi)
    deallocate(data_B0_csr)

end subroutine wq_get_data_csr
