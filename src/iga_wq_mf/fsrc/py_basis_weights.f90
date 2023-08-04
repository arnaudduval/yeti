! ==============================================
! In this module, one can find different kind of functions to compute basis and/or weights 
! for IGA-Galerkin and IGA-WQ approaches.   
! ==============================================

subroutine gauss_quadrature_table(order, position, weight)
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: order
    double precision, intent(out) :: position, weight
    dimension :: position(order), weight(order)

    ! Local data
    ! ----------
    double precision :: GaussPdsCoord
    dimension :: GaussPdsCoord(2, order)

    ! Find position and weight in isoparametric space
    call Gauss(order, 1, GaussPdsCoord, 0)
    weight = GaussPdsCoord(1, :)
    position = GaussPdsCoord(2, :)
    
end subroutine gauss_quadrature_table

subroutine get_genbasis_coo(degree, size_kv, knotvector, nb_knots, knots, basis, indices)
    !! Gets in COO format the basis at given knots 
    
    use quadrature_rules
    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: span_tol = 1.d-8
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: basis
    dimension :: basis(nb_knots*(degree+1), 2)
    integer, intent(out) :: indices
    dimension :: indices((degree+1)*nb_knots, 2)

    ! Local data
    ! ----------
    type(genquadrature) :: obj
    call init_genquad(obj, degree, knotvector)
    call get_basis_simplified_coo(obj, nb_knots, knots, basis, indices)

end subroutine get_genbasis_coo

subroutine get_genbasis_csr(degree, size_kv, knotvector, nb_knots, knots, basis, indi, indj)
    !! Gets in CSR format the basis at given knots 

    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: knotvector, knots
    dimension :: knotvector(size_kv), knots(nb_knots)

    double precision, intent(out) :: basis
    dimension :: basis(nb_knots*(degree+1), 2)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(nb_knots*(degree+1))

    ! Local data
    ! ----------
    integer :: size_data
    integer, dimension(:, :), allocatable :: indices
    double precision, dimension(:, :), allocatable :: basis_coo

    size_data = nb_knots*(degree+1)
    allocate(basis_coo(size_data, 2), indices(size_data, 2))
    call get_genbasis_coo(degree, size_kv, knotvector, nb_knots, knots, basis_coo, indices)

    ! Convert COO to CSR format
    call coo2csr(2, size_kv-degree-1, size_data, basis_coo, indices(:, 1), indices(:, 2), basis, indj, indi)
    deallocate(basis_coo)

end subroutine get_genbasis_csr

subroutine wq_getbasisweights_coo(degree, size_kv, knotvector, nbqp, quadptspos, nbctrlpts, &
                                B0shape, B1shape, size_data, basis, weights, indices, method)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use quadrature_rules
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, size_data, nbctrlpts, nbqp, method
    double precision, intent(in) :: knotvector, quadptspos
    dimension :: knotvector(size_kv), quadptspos(nbqp)
    integer, intent(in) :: B0shape, B1shape
    dimension :: B0shape(nbctrlpts, 2), B1shape(nbctrlpts, 2)

    double precision, intent(out) :: basis, weights
    dimension :: basis(size_data, 2), weights(size_data, 4)
    integer, intent(out) :: indices
    dimension :: indices(size_data, 2)

    ! Local data
    ! ----------
    type(weightedquadrature) :: obj

    if (nbctrlpts.ne.size_kv-degree-1) stop 'Size problem in wq get basis coo'
    call init_wq(obj, degree, knotvector, quadptspos, B0shape, B1shape, method)
    if (method.eq.1) then
        call findbasisweightrules_md1_wq(obj, basis, weights, indices)
    else if (method.eq.2) then
        call findbasisweightrules_md2_wq(obj, basis, weights, indices)
    else
        stop 'Unknown method'
    end if

end subroutine wq_getbasisweights_coo

subroutine wq_getbasisweights_csr(degree, size_kv, knotvector, nbqp, quadptspos, nbctrlpts, &
                        B0shape, B1shape, size_data, basis, weights, indi, indj, method)
    !! Gets in COO format basis and weights in IGA-WQ approach

    use quadrature_rules
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, size_data, nbctrlpts, nbqp, method
    double precision, intent(in) :: knotvector, quadptspos
    dimension :: knotvector(size_kv), quadptspos(nbqp)
    integer, intent(in) :: B0shape, B1shape
    dimension :: B0shape(nbctrlpts, 2), B1shape(nbctrlpts, 2)

    double precision, intent(out), target :: basis, weights
    dimension :: basis(size_data, 2), weights(size_data, 4)
    integer, intent(out) :: indi, indj
    dimension :: indi(size_kv-degree), indj(size_data)

    ! Local data
    ! ----------
    integer :: indices(size_data, 2)
    double precision, dimension(:, :), allocatable :: dummy1, dummy2
    
    call wq_getbasisweights_coo(degree, size_kv, knotvector, nbqp, quadptspos, nbctrlpts, &
                        B0shape, B1shape, size_data, basis, weights, indices, method)

    ! Convert COO to CSR format
    allocate(dummy1(size_data, 6), dummy2(size_data, 6))
    dummy1(:, :2) = basis; dummy1(:, 3:) = weights; dummy2 = 0.d0
    call coo2csr(6, size_kv-degree-1, size_data, dummy1(1:size_data, 1:6), &
                indices(:, 1), indices(:, 2), dummy2, indj, indi)
    basis = dummy2(:, :2); weights = dummy2(:, 3:)
    deallocate(dummy1, dummy2)

end subroutine wq_getbasisweights_csr
