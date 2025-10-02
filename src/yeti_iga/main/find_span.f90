!! Utility routine copied outside of Fortran module nurbspatch

subroutine find_span(val, uknot, nkv, index)
    !! Find span of a given parameter value
    !! (equivalent to routine binarysearch in module nurbspatch)

    !! Input variables
    !! ---------------
    !! val: parameter value
    !! uknot: knot vector
    !! nkv: size of knotvector

    !! Output variable
    !! ---------------
    !! index: knot span

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in) :: nkv
    real(c_double), intent(in) :: val
    real(c_double), dimension(nkv), intent(in) :: uknot

    integer(c_int), intent(out):: index

    integer(c_int) :: low,high,middle

    low = 1
    high = nkv

    do while (low < high)
        middle = (low + high)/2

        if (val >= uknot(middle) .and. val < uknot(middle+1)) then
            index = middle
            return
        else if (val< uknot(middle)) then
            high = middle
        else
            low  = middle+1
        end if
    end do

    ! if not found (limit case) : return low
    index = low
end subroutine find_span


