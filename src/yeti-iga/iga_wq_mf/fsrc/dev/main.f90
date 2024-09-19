program main 

    integer, allocatable :: a(:), b(:)

    allocate(a(3))
    a = [ 1, 2, 3 ]
    allocate(b(4))
    b = [ 1, 5, 10, 4 ]

    call move_alloc(a, b)
    print *, allocated(a), allocated(b)
    print *, b

    ! implicit none 
    ! ! --------------
    ! integer, parameter :: degree = 2, nb_el = 8
    ! integer :: nnz_B, size_qp, nnz_I
    ! double precision, allocatable, dimension(:) :: qp_pos, data_B0, data_B1, &
    !                                             data_W00, data_W11
    ! integer, allocatable, dimension(:, :) :: data_ind 

    ! integer :: nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg
    ! double precision, dimension(:), allocatable :: nodes, qp_wq, knotvector
    ! double precision, dimension(:), allocatable :: qp_cgg, w_cgg
    ! integer, dimension(:,:), allocatable :: table_funct_spans, table_spans_IGA, table_spans_WQ
    ! double precision, dimension(:,:), allocatable :: B0_cgg, B1_cgg, B0_wq, B1_wq
    ! integer, dimension(:,:), allocatable :: B_shape, B0_shape, B1_shape
    ! double precision, dimension(:,:), allocatable :: B0, B1, W00, W11
    ! integer :: i

    ! ! Eval properties
    ! call wq_set_properties(degree, nb_el, multiplicity, maxrule, nb_ctrlpts, &
    ! size_kv, nb_qp_wq, nb_qp_cgg)

    ! ! Eval nodes
    ! allocate(nodes(nb_el+1))
    ! call get_parametric_nodes(nb_el, nodes)
    ! print*, nodes
    ! print*, '--------------------'

    ! ! Eval WQ quadrature points
    ! allocate(qp_wq(nb_qp_wq))
    ! call wq_get_qp_positions(degree, nb_el, maxrule, nodes, nb_qp_wq, qp_wq)
    ! print*, qp_wq
    ! print*, '--------------------'

    ! ! Find knot-vector
    ! allocate(knotvector(size_kv))
    ! call get_knotvector(degree, nb_el, nodes, size_kv, knotvector, multiplicity)
    ! print*, knotvector
    ! print*, '--------------------'

    ! ! Find IGA quadrature points
    ! allocate(qp_cgg(nb_qp_cgg), w_cgg(nb_qp_cgg))
    ! call iga_get_qp_positions_weights(degree, nb_el, nodes, nb_qp_cgg, qp_cgg, w_cgg)
    ! print*, qp_cgg
    ! print*, w_cgg
    ! print*, '--------------------'

    ! ! Find table of functions over spans
    ! allocate(table_funct_spans(nb_el, degree+1))
    ! call set_table_functions_spans(degree, nb_el, multiplicity, table_funct_spans)
    ! do i = 1, nb_el
    !     print*, table_funct_spans(i, :)
    ! end do
    ! print*, '--------------------'

    ! ! Find table of spans of IGA points
    ! allocate(table_spans_IGA(nb_qp_cgg, 2))
    ! call set_table_spans(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, &
    ! nb_qp_cgg, qp_cgg, table_spans_IGA)
    ! do i = 1, nb_qp_cgg
    !     print*, table_spans_IGA(i, :)
    ! end do
    ! print*, '--------------------'

    ! ! Find table of spans of IGA points
    ! allocate(table_spans_WQ(nb_qp_wq, 2))
    ! call set_table_spans(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, &
    ! nb_qp_wq, qp_wq, table_spans_WQ)
    ! do i = 1, nb_qp_wq
    !     print*, table_spans_WQ(i, :)
    ! end do
    ! print*, '--------------------'

    ! ! Find basis in IGA
    ! allocate(B0_cgg(nb_ctrlpts, nb_qp_cgg), B1_cgg(nb_ctrlpts, nb_qp_cgg))
    ! call get_basis(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, nb_qp_cgg, qp_cgg, & 
    ! table_funct_spans, B0_cgg, B1_cgg)
    ! do i = 1, nb_ctrlpts
    !     print*, B0_cgg(i, :)
    !     print*, '**'
    ! end do
    ! print*, '--------------------'

    ! ! Find basis in WQ
    ! allocate(B0_wq(nb_ctrlpts, nb_qp_wq), B1_wq(nb_ctrlpts, nb_qp_wq))
    ! call get_basis(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, nb_qp_wq, qp_wq, & 
    ! table_funct_spans, B0_wq, B1_wq)
    ! do i = 1, nb_ctrlpts
    !     print*, B0_wq(i, :)
    !     print*, '**'
    ! end do
    ! print*, '--------------------'

    ! ! Find shape of B in IGA
    ! allocate(B_shape(nb_ctrlpts, 2))
    ! call iga_get_B_shape(degree, nb_el, nb_ctrlpts, multiplicity, B_shape)
    ! do i = 1, nb_ctrlpts
    !     print*, B_shape(i, :)
    ! end do
    ! print*, '--------------------'

    ! ! Find shape of B0 and B1 in WQ
    ! allocate(B0_shape(nb_ctrlpts, 2), B1_shape(nb_ctrlpts, 2))
    ! call wq_get_nnz_B0_B1_shape(degree, nb_el, nb_ctrlpts, nb_qp_wq, multiplicity, maxrule, B0_shape, B1_shape)
    ! do i = 1, nb_ctrlpts
    !     print*, B0_shape(i, :)
    ! end do
    ! print*, '--------------------'


    ! ! Testing get basis and weights
    ! call wq_set_properties(degree, nb_el, 1, 2, nb_ctrlpts, &
    ! size_kv, nb_qp_wq, nb_qp_cgg)

    ! allocate(B0(nb_ctrlpts, nb_qp_wq), B1(nb_ctrlpts, nb_qp_wq), &
    !         W00(nb_ctrlpts, nb_qp_wq), W11(nb_ctrlpts, nb_qp_wq))
    ! call wq_get_basis_weights_generalized(degree, nb_el, nb_ctrlpts, size_kv, nb_qp_wq, nb_qp_cgg, maxrule, &
    ! B0, B1, W00, W11)

    ! call wq_get_size_data(degree, nb_el, nnz_B, size_qp)

    ! allocate(qp_pos(size_qp), data_B0(nnz_B), data_B1(nnz_B), &
    !         data_W00(nnz_B), data_W11(nnz_B), data_ind(nnz_B, 2))
    ! call wq_get_data( degree, nb_el, nnz_B, size_qp, qp_pos, &
    !     data_B0, data_B1, data_W00, data_W11, data_ind, nnz_I)

end program