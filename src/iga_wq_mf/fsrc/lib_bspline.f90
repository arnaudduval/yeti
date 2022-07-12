! ==========================
! module :: Bspline 
! author :: Joaquin Cornejo
! modules :: algebra(coo2matrix), 
!            dersbasisfuns.f90 
! ==========================

subroutine find_unique_vector(nnz, vec, vec_unique)
    !! It returns the non-repreated values of a vector

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nnz 
    double precision, intent(in) :: vec
    dimension :: vec(nnz)

    double precision, intent(out) :: vec_unique
    dimension :: vec_unique(nnz+1)

    ! Local data
    ! --------------
    integer :: i, num, nnz_nr
    logical, dimension(nnz) :: mask

    ! Initialize
    vec_unique = 0.d0

    ! Define mask 
    mask = .false.

    do i = 1, nnz

        ! Count the number of occurrences of this element:  
        num = count(vec(i).eq.vec)
    
        if (num.eq.1) then  
            ! There is only one, flag it:  
            mask(i) = .true.  
        else  
            !  Flag this value only if it hasn't already been flagged:  
            if (.not. any(vec(i).eq.vec .and. mask)) then
                mask(i) = .true.  
            end if
        end if
    
    end do

    ! Return only flagged elements
    nnz_nr = count(mask)
    vec_unique(1:nnz_nr) = pack(vec, mask)
    vec_unique(nnz+1) = dble(nnz_nr)

end subroutine find_unique_vector

subroutine find_knotvector_span(degree, size_kv, knotvector, x, span)
    !! Finds the span of a given knot within the knot-vector

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv 
    double precision, intent(in) :: knotvector, x 
    dimension :: knotvector(size_kv)

    integer, intent(out) :: span 

    ! Set first value of span
    span = degree + 2
    
    ! Find span
    do while ((span.lt.(size_kv-degree)) &
            .and.((knotvector(span)-x).le.span_tol))

        ! Update value
        span = span + 1
    end do
    
    ! Set result
    span = span - 1 

end subroutine find_knotvector_span

subroutine find_parametric_span(size_kv, nodes, x, span)
    !! Finds the span of the knot within the parametric space

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_kv
    double precision, intent(in) :: nodes, x 
    dimension :: nodes(size_kv+1)

    integer, intent(out) :: span 

    ! Local data
    ! -------------
    integer :: size_nodes

    ! Initialize 
    size_nodes = int(nodes(size_kv+1))

    ! Set first value of span
    span = 2
    
    ! Find span
    do while ((span.lt.size_nodes) &
            .and.((nodes(span)-x).le.span_tol))

        ! Update value 
        span = span + 1
    end do
    
    ! Set result
    span = span - 1 

end subroutine find_parametric_span

subroutine find_multiplicity(size_kv, knotvector, x, mult)

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! ------------------- 
    integer, intent(in) :: size_kv
    double precision, intent(in) :: knotvector, x 
    dimension :: knotvector(size_kv)

    integer, intent(out) :: mult 

    ! Local variables
    ! ----------------
    integer :: i

    ! Set up multiplicity
    mult = 0

    do i = 1, size_kv
        if (abs(x-knotvector(i)).le.tol) then 
            mult = mult + 1
        end if
    end do

end subroutine find_multiplicity

subroutine set_table_functions_spans(degree, size_kv, nodes, knotvector, table)
    !! Sets the table of functions on every span

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv
    double precision, intent(in) :: nodes, knotvector
    dimension :: nodes(size_kv+1), knotvector(size_kv)

    integer, intent(out) :: table
    dimension :: table(int(nodes(size_kv+1))-1, degree+1)

    ! Local data
    ! -------------
    integer :: i, j, multiplicity

    ! Initialize 
    table = 0

    ! Fist line of the table
    do j = 1, degree+1
        table(1, j) = j 
    end do

    ! Set table of functions on span 
    do i = 2, int(nodes(size_kv+1))-1
        call find_multiplicity(size_kv, knotvector, nodes(i), multiplicity)
        table(i, 1) = table(i-1, 1) + multiplicity
        do j = 2, degree+1
            table(i, j) = table(i, 1) + j - 1
        end do
    end do

end subroutine set_table_functions_spans

subroutine get_basis(degree, size_kv, nodes, knotvector, nb_knots, knots, B0, B1)
    !! Finds the basis for every given knot 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, size_kv, nb_knots
    double precision, intent(in) :: nodes, knotvector, knots
    dimension :: nodes(size_kv+1), knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) :: B0, B1
    dimension :: B0(size_kv-degree-1, nb_knots), B1(size_kv-degree-1, nb_knots)

    ! Local data
    ! -------------
    integer :: i, j, k, nb_ctrlpts, para_span, kv_span, nbel
    integer :: functions_span
    dimension :: functions_span(degree+1)
    integer, allocatable, dimension(:, :) :: table_functions_span
    double precision :: B0temp, B1temp
    dimension :: B0temp(degree+1), B1temp(degree+1)
    double precision :: data_B0, data_B1
    dimension :: data_B0((degree+1)*nb_knots), data_B1((degree+1)*nb_knots)
    integer :: ind
    dimension :: ind((degree+1)*nb_knots, 2)

    ! Set number of control points
    nb_ctrlpts = size_kv - degree - 1

    ! Get non repeated knot vector
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
            data_B0(k) = B0temp(j)
            data_B1(k) = B1temp(j)
            ind(k, :) = [functions_span(j), i]                                
        end do
    end do

    ! Matrix construction
    call coo2dense(size(data_B0), ind(:, 1), ind(:, 2), data_B0, nb_ctrlpts, nb_knots, B0)
    call coo2dense(size(data_B1), ind(:, 1), ind(:, 2), data_B1, nb_ctrlpts, nb_knots, B1)

end subroutine get_basis

subroutine create_knotvector(degree, nbel, multiplicity, nodes, knotvector)
    !! Gets an open uniform maximum regularity knot-vector

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in):: degree, nbel, multiplicity

    double precision, intent(out) :: nodes, knotvector 
    dimension :: knotvector(nbel+2*degree+1), nodes(nbel+2*degree+2)

    ! Local data
    ! -------------
    integer ::  i, j, c

    ! Create nodes
    ! =============
    ! Assign first and last values
    nodes(1) = 0.d0
    nodes(nbel+1) = 1.d0

    ! Assign values
    do i = 2, nbel 
        nodes(i) = dble(i - 1)/dble(nbel) 
    end do

    nodes(nbel+2*degree+2) = nbel + 1

    ! Create knotvector
    ! =================
    ! Set p+1 first values of knot vector 
    c = 1
    do i = 1, nbel+1
        knotvector(c) = 0.d0
        c = c + 1
    end do

    do i = 2, nbel
        do j = 1, multiplicity
            knotvector(c) = nodes(i)
            c = c + 1
        end do
    end do

    ! Set p+1 last values of knot vector 
    do i = 1, degree+1
        knotvector(c) = 1.d0
        c = c + 1
    end do
        
end subroutine create_knotvector