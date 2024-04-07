! ---------------
! Tensor algebra 
! ---------------

module tensormode
    use omp_lib
    implicit none

    type tensoroperator
    integer :: sizelist=4
    integer, dimension(:), allocatable :: nclist
    double precision, allocatable, dimension(:) :: tensor
    end type tensoroperator

    contains

    subroutine initialize_operator(obj, sizelist, nclist, nnz, array)
        implicit none
        ! Input /  output data
        ! --------------------
        type(tensoroperator) :: obj
        integer, intent(in) :: sizelist, nnz
        integer, intent(in) :: nclist
        dimension :: nclist(sizelist)
        double precision, intent(in) :: array
        dimension :: array(nnz)

        if (sizelist.ne.obj%sizelist) stop 'Only forth-rank tensors'
        if (any(nclist.le.0)) stop 'Only positive integers'
        if (nnz.ne.product(nclist)) stop 'Size problem'
        if (.not.allocated(obj%nclist)) allocate(obj%nclist(sizelist))
        if (.not.allocated(obj%tensor)) allocate(obj%tensor(nnz))
        obj%nclist = nclist; obj%tensor = array

    end subroutine initialize_operator

    subroutine tensor_n_mode_product_dM2(obj, nr, nc, U, mode, sizelist, newnclist)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = X(nc_u, nc_v, nc_w, nc_t)
        !! Matrix U = U(nr, nc)
        !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
        !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

        implicit none
        ! Input / output data
        ! -------------------
        type(tensoroperator) :: obj
        integer, intent(in) :: nr, nc, mode, sizelist 
        integer, intent(in) :: newnclist
        dimension :: newnclist(sizelist)
        double precision, intent(in) :: U
        dimension :: U(nr, nc)

        ! Local data
        ! ----------
        double precision, allocatable, dimension(:) :: newtensor
        double precision, allocatable, dimension(:, :) :: Rt, Xt
        integer :: jj, kk, ll, genpos, nclist(obj%sizelist)

        if ((mode.lt.1).or.(mode.gt.4)) stop 'Only 1, 2, 3, 4 modes for tensor-matrix operations'
        if (.not.allocated(obj%nclist)) stop 'Unknown object'
        if (sizelist.lt.obj%sizelist) stop 'Size problem' 
        if (any(newnclist.le.0)) stop 'Only positive integers'
        if (nr.ne.newnclist(4)) stop 'Size problem'
        allocate(newtensor(product(newnclist)))
        allocate(Xt(obj%nclist(1), obj%nclist(2)), Rt(nr, obj%nclist(2)))
        nclist = obj%nclist

        !$OMP PARALLEL PRIVATE(Xt, Rt, kk, ll, genpos)
        !$OMP DO SCHEDULE(STATIC) 
        do jj = 1, obj%nclist(3)*obj%nclist(4)
            do ll = 1, obj%nclist(2)
                do kk = 1, obj%nclist(1)
                    genpos = kk + (ll-1)*nclist(1) + (jj-1)*nclist(1)*nclist(2)
                    Xt(kk, ll) = obj%tensor(genpos)
                end do
            end do
            Rt = matmul(U, Xt)
            do kk = 1, nr
                do ll = 1, obj%nclist(2)
                    genpos = ll + (jj-1)*nclist(2) + (kk-1)*nclist(2)*nclist(3)*nclist(4)
                    newtensor(genpos) = Rt(kk, ll)
                end do
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        deallocate(Xt, Rt, obj%tensor, obj%nclist)
        allocate(obj%nclist(obj%sizelist))
        obj%nclist = newnclist(1:obj%sizelist)
        allocate(obj%tensor(product(newnclist)))
        obj%tensor = newtensor

    end subroutine tensor_n_mode_product_dM2

    subroutine tensor_n_mode_product_dM(nc_u, nc_v, nc_w, nc_t, X, nr, nc, U, mode, nrR, R)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = X(nc_u, nc_v, nc_w, nc_t)
        !! Matrix U = U(nr, nc)
        !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
        !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: nc_u, nc_v, nc_w, nc_t, nr, nc, mode, nrR
        double precision, intent(in) :: X, U
        dimension :: X(nc_u*nc_v*nc_w*nc_t), U(nr, nc)

        double precision, intent(out) :: R
        dimension :: R(nrR)

        ! Local data
        ! ----------
        double precision, allocatable, dimension(:, :) :: Rt, Xt
        integer :: ju, jv, jw, jt, offset1, offset2

        R = 0.d0
        if (mode.eq.1) then 

            allocate(Xt(nc_u, nc_v), Rt(nr, nc_v))
            !$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jw = 1, nc_w
                    offset1 = (jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jw-1)*nr*nc_v+(jt-1)*nr*nc_v*nc_w
                    do jv = 1, nc_v
                        do ju = 1, nc_u
                            Xt(ju, jv) = X(ju+(jv-1)*nc_u+offset1)
                        end do
                    end do
                    Rt = matmul(U, Xt)
                    do jv = 1, nc_v
                        do ju = 1, nr
                            R(ju+(jv-1)*nr+offset2) = Rt(ju, jv)
                        end do
                    end do
                end do
            end do
            !$OMP END DO 
            !$OMP END PARALLEL
            deallocate(Xt, Rt)

        else if (mode.eq.2) then 

            allocate(Xt(nc_v, nc_u), Rt(nr, nc_u))
            !$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jw = 1, nc_w
                    offset1 = (jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jw-1)*nc_u*nr+(jt-1)*nc_u*nr*nc_w
                    do ju = 1, nc_u
                        do jv = 1, nc_v
                            Xt(jv, ju) = X(ju+(jv-1)*nc_u+offset1)
                        end do
                    end do
                    Rt = matmul(U, Xt)
                    do ju = 1, nc_u
                        do jv = 1, nr
                            R(ju+(jv-1)*nc_u+offset2) = Rt(jv, ju)
                        end do
                    end do
                end do
            end do
            !$OMP END DO 
            !$OMP END PARALLEL
            deallocate(Xt, Rt)
            
        else if (mode.eq.3) then 

            allocate(Xt(nc_w, nc_u), Rt(nr, nc_u))
            !$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jv = 1, nc_v
                    offset1 = (jv-1)*nc_u+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jv-1)*nc_u+(jt-1)*nc_u*nc_v*nr
                    do ju = 1, nc_u
                        do jw = 1, nc_w
                            Xt(jw, ju) = X(ju+(jw-1)*nc_u*nc_v+offset1)
                        end do
                    end do
                    Rt = matmul(U, Xt)
                    do ju = 1, nc_u
                        do jw = 1, nr
                            R(ju+(jw-1)*nc_u*nc_v+offset2) = Rt(jw, ju)
                        end do
                    end do
                end do
            end do
            !$OMP END DO 
            !$OMP END PARALLEL
            deallocate(Xt, Rt)

        else if (mode.eq.4) then 
            
            allocate(Xt(nc_t, nc_u), Rt(nr, nc_u))
            !$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jw = 1, nc_w
                do jv = 1, nc_v
                    offset1 = (jv-1)*nc_u+(jw-1)*nc_u*nc_v
                    offset2 = (jv-1)*nc_u+(jw-1)*nc_u*nc_v
                    do ju = 1, nc_u
                        do jt = 1, nc_t
                            Xt(jt, ju) = X(ju+(jt-1)*nc_u*nc_v*nc_w+offset1)
                        end do
                    end do
                    Rt = matmul(U, Xt)
                    do ju = 1, nc_u
                        do jt = 1, nr
                            R(ju+(jt-1)*nc_u*nc_v*nc_w+offset2) = Rt(jt, ju)
                        end do
                    end do
                end do
            end do
            !$OMP END DO 
            !$OMP END PARALLEL
            deallocate(Xt, Rt)

        end if

    end subroutine tensor_n_mode_product_dM

    subroutine tensor_n_mode_product_spM2(obj, nr, nnz, dat, indi, indj, mode, sizelist, newnclist)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = X(nc_u, nc_v, nc_w, nc_t)
        !! Matrix U = U(nr, nc)
        !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
        !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

        implicit none
        ! Input / output data
        ! -------------------
        type(tensoroperator) :: obj
        integer, intent(in) :: nr, nnz, mode, sizelist 
        integer, intent(in) :: newnclist
        dimension :: newnclist(sizelist)
        integer, intent(in) :: indi, indj
        dimension :: indi(nr+1), indj(nnz)
        double precision, intent(in) :: dat
        dimension :: dat(nnz)

        ! Local data
        ! ----------
        double precision, allocatable, dimension(:) :: newtensor
        double precision, allocatable, dimension(:, :) :: Rt, Xt
        integer :: jj, kk, ll, genpos, nclist(obj%sizelist)

        if ((mode.lt.1).or.(mode.gt.4)) stop 'Only 1, 2, 3, 4 modes for tensor-matrix operations'
        if (.not.allocated(obj%nclist)) stop 'Unknown object'
        if (sizelist.lt.obj%sizelist) stop 'Size problem' 
        if (any(newnclist.le.0)) stop 'Only positive integers'
        if (nr.ne.newnclist(4)) stop 'Size problem'
        allocate(newtensor(product(newnclist)))
        allocate(Xt(obj%nclist(1), obj%nclist(2)), Rt(nr, obj%nclist(2)))
        nclist = obj%nclist

        !$OMP PARALLEL PRIVATE(Xt, Rt, kk, ll, genpos)
        !$OMP DO SCHEDULE(STATIC) 
        do jj = 1, obj%nclist(3)*obj%nclist(4)
            do ll = 1, obj%nclist(2)
                do kk = 1, obj%nclist(1)
                    genpos = kk + (ll-1)*nclist(1) + (jj-1)*nclist(1)*nclist(2)
                    Xt(kk, ll) = obj%tensor(genpos)
                end do
            end do
            call spmat_dot_dmat(nr, nnz, indi, indj, dat, size(Xt, dim=1), size(Xt, dim=2), Xt, Rt)
            do kk = 1, nr
                do ll = 1, obj%nclist(2)
                    genpos = ll + (jj-1)*nclist(2) + (kk-1)*nclist(2)*nclist(3)*nclist(4)
                    newtensor(genpos) = Rt(kk, ll)
                end do
            end do
        end do
        !$OMP END DO 
        !$OMP END PARALLEL

        deallocate(Xt, Rt, obj%tensor, obj%nclist)
        allocate(obj%nclist(obj%sizelist))
        obj%nclist = newnclist(1:obj%sizelist)
        allocate(obj%tensor(product(newnclist)))
        obj%tensor = newtensor

    end subroutine tensor_n_mode_product_spM2

    subroutine tensor_n_mode_product_spM(nc_u, nc_v, nc_w, nc_t, X, nr, nnz, dat, indi, indj, mode, nrR, R)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = X(nc_u, nc_v, nc_w)
        !! Matrix U = U(nr, nc). Since U is in CSR format, nc is not necessary to be declared
        !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
        !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

        implicit none
        ! Input / output data 
        ! -------------------
        integer, intent(in) :: nc_u, nc_v, nc_w, nc_t, nr, nnz, mode, nrR
        integer, intent(in) :: indi, indj
        dimension :: indi(nr+1), indj(nnz)
        double precision, intent(in) :: X, dat
        dimension :: X(nc_u*nc_v*nc_w*nc_t), dat(nnz)

        double precision, intent(out) :: R
        dimension :: R(nrR)

        ! Local data
        ! ----------
        double precision, allocatable, dimension(:, :) :: Xt, Rt
        integer :: ju, jv, jw, jt, offset1, offset2

        if (mode.eq.1) then 

            allocate(Xt(nc_u, nc_v), Rt(nr, nc_v))
            !*$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !*$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jw = 1, nc_w
                    offset1 = (jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jw-1)*nr*nc_v+(jt-1)*nr*nc_v*nc_w
                    do jv = 1, nc_v
                        do ju = 1, nc_u
                            Xt(ju, jv) = X(ju+(jv-1)*nc_u+offset1)
                        end do
                    end do
                    call spmat_dot_dmat(nr, nnz, indi, indj, dat, nc_u, nc_v, Xt, Rt)
                    do jv = 1, nc_v
                        do ju = 1, nr
                            R(ju+(jv-1)*nr+offset2) = Rt(ju, jv)
                        end do
                    end do
                end do
            end do
            !*$OMP END DO
            !*$OMP END PARALLEL
            deallocate(Xt, Rt)

        else if (mode.eq.2) then 

            allocate(Xt(nc_v, nc_u), Rt(nr, nc_u))
            !*$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !*$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jw = 1, nc_w
                    offset1 = (jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jw-1)*nc_u*nr+(jt-1)*nc_u*nr*nc_w
                    do ju = 1, nc_u
                        do jv = 1, nc_v
                            Xt(jv, ju) = X(ju+(jv-1)*nc_u+offset1)
                        end do
                    end do
                    call spmat_dot_dmat(nr, nnz, indi, indj, dat, nc_v, nc_u, Xt, Rt)
                    do ju = 1, nc_u
                        do jv = 1, nr
                            R(ju+(jv-1)*nc_u+offset2) = Rt(jv, ju)
                        end do
                    end do
                end do
            end do
            !*$OMP END DO
            !*$OMP END PARALLEL
            deallocate(Xt, Rt)

        else if (mode.eq.3) then 

            allocate(Xt(nc_w, nc_u), Rt(nr, nc_u))
            !*$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !*$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jt = 1, nc_t
                do jv = 1, nc_v
                    offset1 = (jv-1)*nc_u+(jt-1)*nc_u*nc_v*nc_w
                    offset2 = (jv-1)*nc_u+(jt-1)*nc_u*nc_v*nr
                    do ju = 1, nc_u
                        do jw = 1, nc_w
                            Xt(jw, ju) = X(ju+(jw-1)*nc_u*nc_v+offset1)
                        end do
                    end do
                    call spmat_dot_dmat(nr, nnz, indi, indj, dat, nc_w, nc_u, Xt, Rt)
                    do ju = 1, nc_u
                        do jw = 1, nr
                            R(ju+(jw-1)*nc_u*nc_v+offset2) = Rt(jw, ju)
                        end do
                    end do
                end do
            end do
            !*$OMP END DO
            !*$OMP END PARALLEL
            deallocate(Xt, Rt)

        else if (mode.eq.4) then 
            
            allocate(Xt(nc_t, nc_u), Rt(nr, nc_u))
            !*$OMP PARALLEL PRIVATE(Xt, Rt, ju, jv, jw, jt, offset1, offset2)
            !*$OMP DO COLLAPSE(2) SCHEDULE(STATIC) 
            do jw = 1, nc_w
                do jv = 1, nc_v
                    offset1 = (jv-1)*nc_u+(jw-1)*nc_u*nc_v
                    offset2 = (jv-1)*nc_u+(jw-1)*nc_u*nc_v
                    do ju = 1, nc_u
                        do jt = 1, nc_t
                            Xt(jt, ju) = X(ju+(jt-1)*nc_u*nc_v*nc_w+offset1)
                        end do
                    end do
                    call spmat_dot_dmat(nr, nnz, indi, indj, dat, nc_t, nc_u, Xt, Rt)
                    do ju = 1, nc_u
                        do jt = 1, nr
                            R(ju+(jt-1)*nc_u*nc_v*nc_w+offset2) = Rt(jt, ju)
                        end do
                    end do
                end do
            end do
            !*$OMP END DO 
            !*$OMP END PARALLEL
            deallocate(Xt, Rt)

        end if

    end subroutine tensor_n_mode_product_spM

end module tensormode
    
subroutine sumfacto2d_dM(nr_u, nc_u, nr_v, nc_v, Mu, Mv, array_in, array_out)
    !! Evaluates a dot product between a tensor 2D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
    use tensormode
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v
    double precision, intent(in) :: Mu, Mv
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data 
    ! ----------
    type(tensoroperator) :: obj

    call initialize_operator(obj, 4, (/nc_u, nc_v, 1, 1/), size(array_in), array_in)
    call tensor_n_mode_product_dM2(obj, nr_u, nc_u, Mu, 1, 4, (/nc_v, 1, 1, nr_u/))
    call tensor_n_mode_product_dM2(obj, nr_v, nc_v, Mv, 2, 4, (/1, 1, nr_u, nr_v/))
    array_out = obj%tensor

    ! ! Local data 
    ! ! ----------
    ! double precision, allocatable, dimension(:) :: R1

    ! allocate(R1(nr_u*nc_v))
    ! call tensor_n_mode_product_dM(nc_u, nc_v, 1, 1, array_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    ! call tensor_n_mode_product_dM(nr_u, nc_v, 1, 1, R1, nr_v, nc_v, Mv, 2, size(array_out), array_out)
    ! deallocate(R1)

end subroutine sumfacto2d_dM

subroutine sumfacto3d_dM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, Mu, Mv, Mw, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Matrix Mw = (nb_rows_w, nb_cols_w)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
    use tensormode
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v,nr_w, nc_w
    double precision, intent(in) :: Mu, Mv, Mw
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v), Mw(nr_w, nc_w)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v*nc_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! ! Local data 
    ! ! ----------
    ! type(tensoroperator) :: obj

    ! call initialize_operator(obj, 4, (/nc_u, nc_v, nc_w, 1/), size(array_in), array_in)
    ! call tensor_n_mode_product_dM2(obj, nr_u, nc_u, Mu, 1, 4, (/nc_v, nc_w, 1, nr_u/))
    ! call tensor_n_mode_product_dM2(obj, nr_v, nc_v, Mv, 2, 4, (/nc_w, 1, nr_u, nr_v/))
    ! call tensor_n_mode_product_dM2(obj, nr_w, nc_w, Mw, 3, 4, (/1, nr_u, nr_v, nr_w/))
    ! array_out = obj%tensor

    ! Local data 
    ! ----------
    double precision, allocatable, dimension(:) :: R1, R2

    allocate(R1(nr_u*nc_v*nc_w))
    call tensor_n_mode_product_dM(nc_u, nc_v, nc_w, 1, array_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    allocate(R2(nr_u*nr_v*nc_w))
    call tensor_n_mode_product_dM(nr_u, nc_v, nc_w, 1, R1, nr_v, nc_v, Mv, 2, size(R2), R2)
    deallocate(R1)

    call tensor_n_mode_product_dM(nr_u, nr_v, nc_w, 1, R2, nr_w, nc_w, Mw, 3, size(array_out), array_out)
    deallocate(R2)  

end subroutine sumfacto3d_dM

subroutine sumfacto4d_dM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, Mu, Mv, Mw, Mt, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Matrix Mw = (nb_rows_w, nb_cols_w)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
    use tensormode
    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nr_t, nc_t
    double precision, intent(in) :: Mu, Mv, Mw, Mt
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v), Mw(nr_w, nc_w), Mt(nr_t, nc_t)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v*nc_w*nc_t)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w*nr_t)

    ! Local data 
    ! ----------
    type(tensoroperator) :: obj

    call initialize_operator(obj, 4, (/nc_u, nc_v, nc_w, nc_t/), size(array_in), array_in)
    call tensor_n_mode_product_dM2(obj, nr_u, nc_u, Mu, 1, 4, (/nc_v, nc_w, nc_t, nr_u/))
    call tensor_n_mode_product_dM2(obj, nr_v, nc_v, Mv, 2, 4, (/nc_w, nc_t, nr_u, nr_v/))
    call tensor_n_mode_product_dM2(obj, nr_w, nc_w, Mw, 3, 4, (/nc_t, nr_u, nr_v, nr_w/))
    call tensor_n_mode_product_dM2(obj, nr_t, nc_t, Mt, 4, 4, (/nr_u, nr_v, nr_w, nr_t/))
    array_out = obj%tensor

    ! ! Local data 
    ! ! ----------
    ! double precision, allocatable, dimension(:) :: R1, R2, R3

    ! allocate(R1(nr_u*nc_v*nc_w*nc_t))
    ! call tensor_n_mode_product_dM(nc_u, nc_v, nc_w, nc_t, array_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    ! allocate(R2(nr_u*nr_v*nc_w*nc_t))
    ! call tensor_n_mode_product_dM(nr_u, nc_v, nc_w, nc_t, R1, nr_v, nc_v, Mv, 2, size(R2), R2)
    ! deallocate(R1)

    ! allocate(R3(nr_u*nr_v*nr_w*nc_t))
    ! call tensor_n_mode_product_dM(nr_u, nr_v, nc_w, nc_t, R2, nr_w, nc_w, Mw, 3, size(R3), R3)
    ! deallocate(R2)   

    ! call tensor_n_mode_product_dM(nr_u, nr_v, nr_w, nc_t, R3, nr_t, nc_t, Mt, 4, size(array_out), array_out)
    ! deallocate(R3)

end subroutine sumfacto4d_dM

subroutine sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, data_u, &
                            nnz_v, indi_v, indj_v, data_v, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Vector_in = (nb_cols_u * nb_cols_v)
    use tensormode
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v)
    double precision, intent(in) :: data_u, data_v
    dimension :: data_u(nnz_u), data_v(nnz_v)
    integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                    indj_u(nnz_u), indj_v(nnz_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data
    ! __________
    type(tensoroperator) :: obj

    call initialize_operator(obj, 4, (/nc_u, nc_v, 1, 1/), size(array_in), array_in)
    call tensor_n_mode_product_spM2(obj, nr_u, nnz_u, data_u, indi_u, indj_u, 1, 4, (/nc_v, 1, 1, nr_u/))
    call tensor_n_mode_product_spM2(obj, nr_v, nnz_v, data_v, indi_v, indj_v, 2, 4, (/1, 1, nr_u, nr_v/))    
    array_out = obj%tensor

    ! ! Local data
    ! ! __________
    ! double precision, allocatable, dimension(:) :: R1

    ! allocate(R1(nr_u*nc_v))
    ! call tensor_n_mode_product_spM(nc_u, nc_v, 1, 1, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    ! call tensor_n_mode_product_spM(nr_u, nc_v, 1, 1, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(array_out), array_out)
    ! deallocate(R1)

end subroutine sumfacto2d_spM

subroutine sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, indi_u, indj_u, data_u, &
                        nnz_v, indi_v, indj_v, data_v, nnz_w, indi_w, indj_w, data_w, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector 
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Matrix Mw = (nb_rows_w, nb_cols_w)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
    use tensormode
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v*nc_w)
    double precision, intent(in) :: data_u, data_v, data_w
    dimension :: data_u(nnz_u), data_v(nnz_v), data_w(nnz_w)
    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    type(tensoroperator) :: obj

    call initialize_operator(obj, 4, (/nc_u, nc_v, nc_w, 1/), size(array_in), array_in)
    call tensor_n_mode_product_spM2(obj, nr_u, nnz_u, data_u, indi_u, indj_u, 1, 4, (/nc_v, nc_w, 1, nr_u/))
    call tensor_n_mode_product_spM2(obj, nr_v, nnz_v, data_v, indi_v, indj_v, 2, 4, (/nc_w, 1, nr_u, nr_v/))
    call tensor_n_mode_product_spM2(obj, nr_w, nnz_w, data_w, indi_w, indj_w, 3, 4, (/1, nr_u, nr_v, nr_w/))
    array_out = obj%tensor

    ! ! Local data 
    ! ! ----------
    ! double precision, allocatable, dimension(:) :: R1, R2

    ! allocate(R1(nr_u*nc_v*nc_w))
    ! call tensor_n_mode_product_spM(nc_u, nc_v, nc_w, 1, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    ! allocate(R2(nr_u*nr_v*nc_w))
    ! call tensor_n_mode_product_spM(nr_u, nc_v, nc_w, 1, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(R2), R2)
    ! deallocate(R1)

    ! call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, 1, R2, nr_w, nnz_w, data_w, indi_w, indj_w, 3, size(array_out), array_out)
    ! deallocate(R2)

end subroutine sumfacto3d_spM

subroutine sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, indi_u, indj_u, data_u, &
                        nnz_v, indi_v, indj_v, data_v, nnz_w, indi_w, indj_w, data_w, &
                        nnz_t, indi_t, indj_t, data_t, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector 
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Matrix Mw = (nb_rows_w, nb_cols_w)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
    use tensormode
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u*nc_v*nc_w*nc_t)
    double precision, intent(in) :: data_u, data_v, data_w, data_t
    dimension :: data_u(nnz_u), data_v(nnz_v), data_w(nnz_w), data_t(nnz_t)
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w*nr_t)

    ! Local data 
    ! ----------
    type(tensoroperator) :: obj

    call initialize_operator(obj, 4, (/nc_u, nc_v, nc_w, nc_t/), size(array_in), array_in)
    call tensor_n_mode_product_spM2(obj, nr_u, nnz_u, data_u, indi_u, indj_u, 1, 4, (/nc_v, nc_w, nc_t, nr_u/))
    call tensor_n_mode_product_spM2(obj, nr_v, nnz_v, data_v, indi_v, indj_v, 2, 4, (/nc_w, nc_t, nr_u, nr_v/))
    call tensor_n_mode_product_spM2(obj, nr_w, nnz_w, data_w, indi_w, indj_w, 3, 4, (/nc_t, nr_u, nr_v, nr_w/))
    call tensor_n_mode_product_spM2(obj, nr_t, nnz_t, data_t, indi_t, indj_t, 4, 4, (/nr_u, nr_v, nr_w, nr_t/))
    array_out = obj%tensor

    ! ! Local data 
    ! ! ----------
    ! double precision, allocatable, dimension(:) :: R1, R2, R3

    ! allocate(R1(nr_u*nc_v*nc_w*nc_t))
    ! call tensor_n_mode_product_spM(nc_u, nc_v, nc_w, nc_t, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    ! allocate(R2(nr_u*nr_v*nc_w*nc_t))
    ! call tensor_n_mode_product_spM(nr_u, nc_v, nc_w, nc_t, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(R2), R2)
    ! deallocate(R1)

    ! allocate(R3(nr_u*nr_v*nr_w*nc_t))
    ! call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, nc_t, R2, nr_w, nnz_w, data_w, indi_w, indj_w, 3, size(R3), R3)
    ! deallocate(R2)

    ! call tensor_n_mode_product_spM(nr_u, nr_v, nr_w, nc_t, R3, nr_t, nnz_t, data_t, indi_t, indj_t, 4, size(array_out), array_out)
    ! deallocate(R3)

end subroutine sumfacto4d_spM

! ------------------------------------------
! Sum factorization to compute the diagonals 
! ------------------------------------------

subroutine csr_get_diag_2d(coefs, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, diag)
    !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
    !! Algorithm based on sum factorization adapted to diagonal case 
    !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
    !! Indices must be in CSR format

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_B_v, data_W_u, data_W_v
    dimension ::    data_B_u(nnz_u), data_B_v(nnz_v), &
                    data_W_u(nnz_u), data_W_v(nnz_v)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v)

    ! Local data
    ! ----------
    double precision :: data_BW_u, data_BW_v
    dimension :: data_BW_u(nnz_u), data_BW_v(nnz_v)
    
    data_BW_u = data_B_u*data_W_u
    data_BW_v = data_B_v*data_W_v

    call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, data_BW_u, &
                    nnz_v, indi_v, indj_v, data_BW_v, coefs, diag)

end subroutine csr_get_diag_2d

subroutine csr_get_diag_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, diag)
    !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
    !! Algorithm based on sum factorization adapted to diagonal case 
    !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
    !! Indices must be in CSR format

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v*nc_w)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w
    dimension ::    data_B_u(nnz_u), data_B_v(nnz_v), data_B_w(nnz_w), &
                    data_W_u(nnz_u), data_W_v(nnz_v), data_W_w(nnz_w)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    ! Local data
    ! ----------
    double precision :: data_BW_u, data_BW_v, data_BW_w
    dimension :: data_BW_u(nnz_u), data_BW_v(nnz_v), data_BW_w(nnz_w)
    
    data_BW_u = data_B_u*data_W_u
    data_BW_v = data_B_v*data_W_v
    data_BW_w = data_B_w*data_W_w

    call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, indi_u, indj_u, data_BW_u, &
                    nnz_v, indi_v, indj_v, data_BW_v, nnz_w, indi_w, indj_w, data_BW_w, coefs, diag)

end subroutine csr_get_diag_3d

subroutine csr_get_diag_4d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t, &
                            data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t, diag)
    !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
    !! Algorithm based on sum factorization adapted to diagonal case 
    !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
    !! Indices must be in CSR format

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v*nc_w*nc_t)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, indi_t, indj_t
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w), &
                    indi_t(nr_t+1), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t
    dimension ::    data_B_u(nnz_u), data_B_v(nnz_v), data_B_w(nnz_w), data_B_t(nnz_t), &
                    data_W_u(nnz_u), data_W_v(nnz_v), data_W_w(nnz_w), data_W_t(nnz_t)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w*nr_t)

    ! Local data
    ! ----------
    double precision :: data_BW_u, data_BW_v, data_BW_w, data_BW_t
    dimension :: data_BW_u(nnz_u), data_BW_v(nnz_v), data_BW_w(nnz_w), data_BW_t(nnz_t)
    
    data_BW_u = data_B_u*data_W_u
    data_BW_v = data_B_v*data_W_v
    data_BW_w = data_B_w*data_W_w
    data_BW_t = data_B_t*data_W_t

    call sumfacto4d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, indi_u, indj_u, data_BW_u, &
                        nnz_v, indi_v, indj_v, data_BW_v, nnz_w, indi_w, indj_w, data_BW_w, &
                        nnz_t, indi_t, indj_t, data_BW_t, coefs, diag)

end subroutine csr_get_diag_4d

! ----------------------------
! Fast Diagonalization method
! ----------------------------

subroutine stiffmass_eigendecomposition(nr, nc, univMcoefs, univKcoefs, nnz, indi, indj, &
                                data_B, data_W, eigval, eigvec, Kdiag, Mdiag)
    !! Generalized eigen decomposition KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    double precision, intent(in) :: univMcoefs, univKcoefs
    dimension :: univMcoefs(nc), univKcoefs(nc)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B, data_W
    dimension :: data_B(nnz, 2), data_W(nnz, 4)

    double precision, intent(out) :: eigval, eigvec
    dimension :: eigval(nr), eigvec(nr, nr)
    double precision, intent(out) :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    ! Local data
    ! ----------
    integer :: i, j
    double precision :: data_Bt(nnz)
    double precision, allocatable, dimension(:, :) :: basis, weights, mass, stiff
    
    ! Masse matrix
    data_Bt = data_B(:, 1)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univMcoefs(indj(j))
        end do
    end do

    allocate(basis(nr, nc), weights(nr, nc))
    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, basis)
    call csr2dense(nnz, indi, indj, data_W(:, 1), nr, nc, weights)
    allocate(mass(nr, nr))
    mass = matmul(weights, transpose(basis))

    ! Stiffness matrix
    data_Bt = data_B(:, 2)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univKcoefs(indj(j))
        end do
    end do

    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, basis)
    call csr2dense(nnz, indi, indj, data_W(:, 4), nr, nc, weights)
    allocate(stiff(nr, nr))
    stiff = matmul(weights, transpose(basis))
    deallocate(basis, weights)

    ! Save diagonal of M and K
    do i = 1, nr
        Kdiag(i) = stiff(i, i)
        Mdiag(i) = mass(i, i)
    end do

    ! -----------------------------------
    ! Eigen decomposition KK U = MM U DD
    ! -----------------------------------
    call compute_eigdecomp_pdr(nr, stiff, mass, eigval, eigvec)
    deallocate(stiff, mass)

end subroutine stiffmass_eigendecomposition

subroutine advmass_schurdecomposition(nr, nc, univMcoefs, univAcoefs, nnz, indi, indj, &
                                data_B, data_W, UU, VV, SS, TT, Adiag, Mdiag)
    !! Generalized eigen decomposition KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr, nc, nnz
    double precision, intent(in) :: univMcoefs, univAcoefs
    dimension :: univMcoefs(nc), univAcoefs(nc)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B, data_W
    dimension :: data_B(nnz, 2), data_W(nnz, 4)

    double complex, intent(out) :: UU, VV, SS, TT
    dimension :: UU(nr, nr), VV(nr, nr), SS(nr, nr), TT(nr, nr)
    double precision, intent(out) :: Adiag, Mdiag
    dimension :: Adiag(nr), Mdiag(nr)

    ! Local data
    ! ----------
    integer :: i, j
    double precision :: data_Bt(nnz)
    double precision, allocatable, dimension(:, :) :: densebasis, denseweights, mass, adv
    
    ! Masse matrix
    data_Bt = data_B(:, 1)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univMcoefs(indj(j))
        end do
    end do

    allocate(densebasis(nr, nc), denseweights(nr, nc))
    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, densebasis)
    call csr2dense(nnz, indi, indj, data_W(:, 1), nr, nc, denseweights)
    allocate(mass(nr, nr))
    mass = matmul(denseweights, transpose(densebasis))

    ! Advention matrix
    data_Bt = data_B(:, 2)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univAcoefs(indj(j))
        end do
    end do

    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, densebasis)
    call csr2dense(nnz, indi, indj, data_W(:, 2), nr, nc, denseweights)
    allocate(adv(nr, nr))
    adv = matmul(denseweights, transpose(densebasis))
    deallocate(densebasis, denseweights)

    ! Save diagonal of M and K
    do i = 1, nr
        Adiag(i) = adv(i, i)
        Mdiag(i) = mass(i, i)
    end do

    ! -----------------------------------
    ! Schur decomposition
    ! -----------------------------------
    call compute_schurdecomp_gc(nr, adv, mass, UU, VV, SS, TT)
    deallocate(adv, mass)

end subroutine advmass_schurdecomposition

subroutine find_parametric_diag_2d(nr_u, nr_v, Mu, Mv, Ku, Kv, coefs, diag)
    !! Computes the diagonal given by cu (Mw x Mv x Ku) + cv (Mw x Kv x Mu) + cw (Kw x Mv x Mu)
                        
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nr_v
    double precision, intent(in) :: Mu, Mv, Ku, Kv, coefs
    dimension :: Mu(nr_u), Mv(nr_v), &
                Ku(nr_u), Kv(nr_v), coefs(2)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v)

    diag = 0.d0

    ! Mv x Ku
    call kronvec2d(nr_v, Mv, nr_u, Ku, diag, coefs(1))

    ! Kv x Mu
    call kronvec2d(nr_v, Kv, nr_u, Mu, diag, coefs(2))

end subroutine find_parametric_diag_2d

subroutine find_parametric_diag_3d(nr_u, nr_v, nr_w, Mu, Mv, Mw, Ku, Kv, Kw, coefs, diag)
    !! Computes the diagonal given by cu (Mw x Mv x Ku) + cv (Mw x Kv x Mu) + cw (Kw x Mv x Mu)
                        
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nr_v, nr_w
    double precision, intent(in) :: Mu, Mv, Mw, Ku, Kv, Kw, coefs
    dimension :: Mu(nr_u), Mv(nr_v), Mw(nr_w), &
                Ku(nr_u), Kv(nr_v), Kw(nr_w), coefs(3)

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    diag = 0.d0

    ! Mw x Mv x Ku
    call kronvec3d(nr_w, Mw, nr_v, Mv, nr_u, Ku, diag, coefs(1))

    ! Mw x Kv x Mu
    call kronvec3d(nr_w, Mw, nr_v, Kv, nr_u, Mu, diag, coefs(2))

    ! Kw x Mv x Mu
    call kronvec3d(nr_w, Kw, nr_v, Mv, nr_u, Mu, diag, coefs(3))
    
end subroutine find_parametric_diag_3d
