!! Function for evaluation of NURBS and BSplines and their derivatives
!! Will reaplace functions contained in evalnurbsfcts.f source file

subroutine evaluate_nurbs(xi, R, degrees, spans, dim, ncp_elt)
    !! Input variables
    !! ---------------
    !! xi: parametric coordinates
    !! degrees: degrees for each parametric direction
    !! spans: knot spans (should be already known)
    !! dim: patch dimension
    !! ncp_elt: number of control points per element

    !! Output variables
    !! ----------------
    !! R: NURBS functions (and their derivatives)

    !! Assume trivariate case by default

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), dimension(dim), intent(in) :: xi
    integer(c_int), dimension(dim), intent(in) :: degrees
    integer(c_int), dimension(dim), intent(in) :: spans
    integer(c_int), intent(in) :: dim, ncp_elt

    real(c_double), dimension(ncp_elt), intent(out) :: R

    real(c_double), dimension(1, degrees(1)) :: FN
    real(c_double), dimension(1, degrees(2)) :: FM
    real(c_double), dimension(1, degrees(3)) :: FL
    real(c_double), dimension(3) :: SumXi
    real(c_double) :: SumTot, SumTot_inv

    !! Initialization
    FM(1, 1) = 1.0
    FL(1, 1) = 1.0

    !! Explicitely define 3 dimension
    !! TODO implement a version with dynamic number of dimensions
    ! call gen_dersbasisfuns(spans(1), degrees(1), size(knot_vector_1),
    !                        xi(1), knot_vector_1, 0, FN)

    !! TO BE CONTINUED ....
    write(*,*) "Subroutine evaluate_nurbs is not available"
    call exit(-1)


end subroutine evaluate_nurbs