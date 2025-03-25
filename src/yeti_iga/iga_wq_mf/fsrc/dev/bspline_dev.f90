! ==========================
! module :: Bspline ! ADDING MULTIPLICITY
! author :: Joaquin Cornejo
! modules :: algebra(coo2matrix), 
!            dersbasisfuns.f90 
! ==========================

subroutine get_parametric_nodes(nb_el, nodes)
    !! Gets the knots non-repetited of the knot-vector
    
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_el 

    double precision, intent(out) :: nodes 
    dimension :: nodes(nb_el+1)

    ! Local data
    ! -------------
    integer ::  i

    ! Assign first and last values
    nodes(1) = 0.d0
    nodes(nb_el+1) = 1.d0

    ! Assign values
    do i = 2, nb_el 
        nodes(i) = dble(i - 1)/dble(nb_el) 
    end do

end subroutine get_parametric_nodes

subroutine find_knotvector_span(degree, nb_ctrlpts, size_kv, knotvector, x, span)
    !! Finds the span of the knot x in the knot-vector

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, nb_ctrlpts, size_kv 
    double precision, intent(in) :: knotvector, x 
    dimension ::  knotvector(size_kv)

    integer, intent(out) :: span 

    ! Set first value of span
    span = degree + 2
    
    ! Find span
    do while ((span.lt.nb_ctrlpts+1) &
            .and.((knotvector(span)-x).le.span_tol))

        ! Update value
        span = span + 1
    end do
    
    ! Set result
    span = span - 1 

end subroutine find_knotvector_span

subroutine find_parametric_span(nb_el, nodes, x, span)
    !! Finds the span of the knot x in parametric space

    use constants_iga_wq_mf
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nb_el 
    double precision, intent(in) :: nodes, x 
    dimension ::  nodes(nb_el+1)

    integer, intent(out) :: span 

    ! Set first value of span
    span = 2
    
    ! Find span
    do while ((span.lt.nb_el+1) &
            .and.((nodes(span)-x).le.span_tol))

        ! Update value 
        span = span + 1
    end do
    
    ! Set result
    span = span - 1 

end subroutine find_parametric_span

subroutine set_table_spans(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, &
                        nb_knots, knots, table)
    !! Creates a table :: 1st col: parametric span  
    !!                    2nd col: knot-vector span 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_ctrlpts, size_kv, nb_knots
    double precision, intent(in) :: nodes, knotvector, knots
    dimension :: nodes(nb_el+1), knotvector(size_kv), knots(nb_knots)

    integer, intent(out) :: table 
    dimension :: table(nb_knots, 2)

    ! Local data
    ! -------------
    integer :: i, span_1, span_2

    do i = 1, nb_knots
        ! Find parametric span
        call find_parametric_span(nb_el, nodes, knots(i), span_1)

        ! Find knot-vector span 
        call find_knotvector_span(degree, nb_ctrlpts, size_kv, knotvector, knots(i), span_2)
        
        ! Set table's i-row
        table(i, :) = [span_1, span_2] 
    end do

end subroutine set_table_spans

subroutine get_knotvector(degree, nb_el, nodes, size_kv, knotvector, multiplicity)
    !! Gets the knot-vector in IGA-WQ-Galerkin approach 
    !! Case of continuity (p-multiplicity)

    implicit none
    ! Input / output data
    ! --------------------
    integer, intent(in):: degree, nb_el, size_kv, multiplicity
    double precision, intent(in) :: nodes
    dimension :: nodes(nb_el+1)

    double precision, intent(out) :: knotvector 
    dimension :: knotvector(size_kv)

    ! Local data
    ! -------------
    integer ::  i, j, c

    ! Set p+1 first values of knot vector 
    c = 1
    do i = 1, degree+1
        knotvector(c) = 0.d0
        c = c + 1
    end do

    do i = 2, nb_el
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
        
end subroutine get_knotvector

subroutine set_table_functions_spans(degree, nb_el, multiplicity, table)
    !! Sets the table of functions on every span
    !! Case of continuity (p-multiplicity) in knot-vector

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, nb_el, multiplicity

    integer, intent(out) :: table
    dimension :: table(nb_el, degree+1)

    ! Local data
    ! -------------
    integer :: i, j

    ! Initialize 
    table = 0

    ! Fist line of the table
    do j = 1, degree + 1
        table(1, j) = j 
    end do

    ! Set table of functions on span 
    do i  = 2, nb_el
        table(i, 1) = table(i-1, 1) + multiplicity
        do j = 2, degree + 1
            table(i, j) = table(i, 1) + j - 1
        end do
    end do

end subroutine set_table_functions_spans

subroutine get_basis(degree, nb_el, nodes, nb_ctrlpts, size_kv, knotvector, nb_knots, knots, & 
                table_functions_span, B0, B1)
    !! Finds the basis for every given knot 

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: degree, nb_el, nb_ctrlpts, size_kv, nb_knots
    integer, intent(in) :: table_functions_span
    dimension :: table_functions_span(nb_el, degree+1)
    double precision, intent(in) :: nodes, knotvector, knots
    dimension :: nodes(nb_el+1), knotvector(size_kv), knots(nb_knots)
    
    double precision, intent(out) ::  B0, B1
    dimension :: B0(nb_ctrlpts, nb_knots), B1(nb_ctrlpts, nb_knots)

    ! Local data
    ! -------------
    double precision :: data_B0, data_B1
    dimension ::    data_B0((degree+1)*nb_knots), &
                    data_B1((degree+1)*nb_knots)
    integer :: indexes_B
    dimension :: indexes_B((degree+1)*nb_knots, 2)
    integer :: table_spans
    dimension :: table_spans(nb_knots, 2)
    double precision :: B0temp, B1temp
    integer :: span, functions_span
    dimension ::    functions_span(degree+1), &
                    B0temp(degree+1), B1temp(degree+1)
    integer :: i, j, k

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
            indexes_B(k, :) = [functions_span(j), i]                                
        end do
    end do

    ! Matrix construction
    call coo2matrix((degree+1)*nb_knots, indexes_B(:, 1), indexes_B(:, 2), data_B0, & 
                    nb_ctrlpts, nb_knots, B0)

    call coo2matrix((degree+1)*nb_knots, indexes_B(:, 1), indexes_B(:, 2), data_B1, & 
                    nb_ctrlpts, nb_knots, B1)

end subroutine get_basis
