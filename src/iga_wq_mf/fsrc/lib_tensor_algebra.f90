! ---------------
! Tensor algebra 
! ---------------

subroutine tensor_n_mode_product_dM(nc_u, nc_v, nc_w, nc_t, X, nr, nc, U, mode, nrR, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = X(nc_u, nc_v, nc_w, nc_t)
    !! Matrix U = U(nr, nc)
    !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
    !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

    use omp_lib
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
    double precision, allocatable, dimension(:, :) :: Rt
    double precision, allocatable, dimension(:, :) :: Xt
    integer :: ju, jv, jw, jt, i, nb_tasks

    R = 0.d0
    if (mode.eq.1) then 

        allocate(Xt(nc_u, nc_v), Rt(nr, nc_v))
        
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(X)/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        Xt(ju, jv) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w)
                    end do
                end do
                Rt = matmul(U, Xt)
                do jv = 1, nc_v
                    do i = 1, nr
                        R(i+(jv-1)*nr+(jw-1)*nr*nc_v+(jt-1)*nr*nc_v*nc_w) = Rt(i, jv)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        deallocate(Xt)
        deallocate(Rt)

    else if (mode.eq.2) then 

        allocate(Xt(nc_v, nc_u*nc_w*nc_t))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(X)/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do ju = 1, nc_u
                    do jv = 1, nc_v
                        Xt(jv, ju+(jw-1)*nc_u+(jt-1)*nc_u*nc_w) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        allocate(Rt(nr, nc_u*nc_w*nc_t))
        Rt = matmul(U, Xt)
        deallocate(Xt)

        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(R)/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do ju = 1, nc_u
                    do i = 1, nr
                        R(ju+(i-1)*nc_u+(jw-1)*nc_u*nr+(jt-1)*nc_u*nr*nc_w) = Rt(i, ju+(jw-1)*nc_u+(jt-1)*nc_u*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(Rt)
        
    else if (mode.eq.3) then 

        allocate(Xt(nc_w, nc_u*nc_v*nc_t))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(X)/nb_tasks) 
        do jt = 1, nc_t
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do jw = 1, nc_w
                        Xt(jw, ju+(jv-1)*nc_u+(jt-1)*nc_u*nc_v) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        allocate(Rt(nr, nc_u*nc_v*nc_t))
        Rt = matmul(U, Xt)
        deallocate(Xt)

        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(R)/nb_tasks) 
        do jt = 1, nc_t
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do i = 1, nr
                        R(ju+(jv-1)*nc_u+(i-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nr) = Rt(i, ju+(jv-1)*nc_u+(jt-1)*nc_u*nc_v)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(Rt)

    else if (mode.eq.4) then 
        
        allocate(Xt(nc_t, nc_u*nc_v*nc_w))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(X)/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do jt = 1, nc_t
                        Xt(jt, ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        allocate(Rt(nr, nc_u*nc_v*nc_w))
        Rt = matmul(U, Xt)
        deallocate(Xt)

        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(4) SCHEDULE(STATIC, size(R)/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do i = 1, nr
                        R(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(i-1)*nc_u*nc_v*nc_w) = Rt(i, ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(Rt)
        
    end if

end subroutine tensor_n_mode_product_dM

subroutine tensor_n_mode_product_spM(s_in, X, nr, nc, nnz, U, indi, indj, mode, s_out, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = X(nc_u, nc_v, nc_w)
    !! Matrix U = U(nr, nc). Since U is in CSR format, nc is not necessary to be declared
    !! Tensor R = R(nu, nv, nw, nt) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw*nt
    !! Ex: if n=1, then nc = nc_u and dim(R) = [nr, nc_v, nc_w, nc_t]

    use omp_lib
    implicit none
    ! Input / output data 
    ! -------------------
    integer, parameter :: dimen = 4
    integer, intent(in) :: nr, nc, nnz, mode, s_in(dimen), s_out(dimen)
    double precision, intent(in) :: X, U
    dimension :: X(s_in(1), s_in(2), s_in(3), s_in(4)), U(nnz)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)

    double precision, intent(out) :: R
    dimension :: R(s_out(1), s_out(2), s_out(3), s_out(4))

    ! Local data
    ! ----------
    integer :: ju, jv, jw, jt, nb_tasks, nb_oper
    double precision :: tmp_out(nr), tmp_in(s_in(mode))

    nb_oper = product(s_in)/s_in(mode)

    if (mode.eq.1) then 

        !$OMP PARALLEL PRIVATE(tmp_in, tmp_out)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_oper/nb_tasks) 
        do jt = 1, s_in(4)
            do jw = 1, s_in(3)
                do jv = 1, s_in(2)
                    tmp_in = X(:, jv, jw, jt)
                    call spmat_dot_dvec(nr, nc, nnz, indi, indj, U, tmp_in, tmp_out)
                    R(:, jv, jw, jt) = tmp_out
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.2) then 

        !$OMP PARALLEL PRIVATE(tmp_in, tmp_out)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_oper/nb_tasks) 
        do jt = 1, s_in(4)
            do jw = 1, s_in(3)
                do ju = 1, s_in(1)
                    call spmat_dot_dvec(nr, nc, nnz, indi, indj, U, X(ju, :, jw, jt), tmp_out)
                    R(ju, :, jw, jt) = tmp_out
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.3) then 

        !$OMP PARALLEL PRIVATE(tmp_in, tmp_out)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_oper/nb_tasks) 
        do jt = 1, s_in(4)
            do jv = 1, s_in(2)
                do ju = 1, s_in(1)
                    tmp_in = X(ju, jv, :, jt)
                    call spmat_dot_dvec(nr, nc, nnz, indi, indj, U, tmp_in, tmp_out)
                    R(ju, jv, :, jt) = tmp_out
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.4) then 

        !$OMP PARALLEL PRIVATE(tmp_in, tmp_out)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_oper/nb_tasks) 
        do jw = 1, s_in(3)
            do jv = 1, s_in(2)
                do ju = 1, s_in(1)
                    tmp_in = X(ju, jv, jw, :)
                    call spmat_dot_dvec(nr, nc, nnz, indi, indj, U, tmp_in, tmp_out)
                    R(ju, jv, jw, :) = tmp_out
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
    end if

end subroutine tensor_n_mode_product_spM

subroutine sumfacto2d_dM(nr_u, nc_u, nr_v, nc_v, Mu, Mv, array_in, array_out)
    !! Evaluates a dot product between a tensor 2D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v
    double precision, intent(in) :: Mu, Mv
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, 1, 1)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, 1, 1)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, 1, 1))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_dM(s_in, array_in, nr_u, nc_u, Mu, 1, s_out, R1)

    s_in = shape(R1); s_out = shape(array_out)
    call tensor_n_mode_product_dM(s_in, R1, nr_v, nc_v, Mv, 2, s_out, array_out)
    deallocate(R1)

end subroutine sumfacto2d_dM

subroutine sumfacto3d_dM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, Mu, Mv, Mw, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Matrix Mw = (nb_rows_w, nb_cols_w)
    !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v,nr_w, nc_w
    double precision, intent(in) :: Mu, Mv, Mw
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v), Mw(nr_w, nc_w)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, nc_w, 1)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, nr_w, 1)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1, R2
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, nc_w, 1))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_dM(s_in, array_in, nr_u, nc_u, Mu, 1, s_out, R1)

    allocate(R2(nr_u, nr_v, nc_w, 1))
    s_in = shape(R1); s_out = shape(R2)
    call tensor_n_mode_product_dM(s_in, R1, nr_v, nc_v, Mv, 2, s_out, R2)
    deallocate(R1)

    s_in = shape(R2); s_out = shape(array_out)
    call tensor_n_mode_product_dM(s_in, R2, nr_w, nc_w, Mw, 3, s_out, array_out)
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

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nr_t, nc_t
    double precision, intent(in) :: Mu, Mv, Mw, Mt
    dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v), Mw(nr_w, nc_w), Mt(nr_t, nc_t)
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, nc_w, nc_t)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, nr_w, nr_t)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1, R2, R3
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, nc_w, nc_t))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_dM(s_in, array_in, nr_u, nc_u, Mu, 1, s_out, R1)

    allocate(R2(nr_u, nr_v, nc_w, nc_t))
    s_in = shape(R1); s_out = shape(R2)
    call tensor_n_mode_product_dM(s_in, R1, nr_v, nc_v, Mv, 2, s_out, R2)
    deallocate(R1)

    allocate(R3(nr_u, nr_v, nr_w, nc_t))
    s_in = shape(R2); s_out = shape(R3)
    call tensor_n_mode_product_dM(s_in, R2, nr_w, nc_w, Mw, 3, s_out, R3)
    deallocate(R2)
    
    s_in = shape(R3); s_out = shape(array_out)
    call tensor_n_mode_product_dM(s_in, R3, nr_t, nc_t, Mt, 4, s_out, array_out)
    deallocate(R3)

end subroutine sumfacto4d_dM

subroutine sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, data_u, &
                            nnz_v, indi_v, indj_v, data_v, array_in, array_out)
    !! Evaluates a dot product between a tensor 3D and a vector using sum factorization
    !! Based on "Preconditioners for IGA" by Montardini
    !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
    !! Matrix Mu = (nb_rows_u, nb_cols_u)
    !! Matrix Mv = (nb_rows_v, nb_cols_v)
    !! Vector_in = (nb_cols_u * nb_cols_v)

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, 1, 1)
    double precision, intent(in) :: data_u, data_v
    dimension :: data_u(nnz_u), data_v(nnz_v)
    integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                    indj_u(nnz_u), indj_v(nnz_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, 1, 1)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, 1, 1))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_spM(s_in, array_in, nr_u, nc_u, nnz_u, data_u, indi_u, indj_u, 1, s_out, R1)

    s_in = shape(R1); s_out = shape(array_out)
    call tensor_n_mode_product_spM(s_in, R1, nr_v, nc_v, nnz_v, data_v, indi_v, indj_v, 2, s_out, array_out)
    deallocate(R1)

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

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nnz_u, nnz_v, nnz_w
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, nc_w, 1)
    double precision, intent(in) :: data_u, data_v, data_w
    dimension :: data_u(nnz_u), data_v(nnz_v), data_w(nnz_w)
    integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, nr_w, 1)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1, R2
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, nc_w, 1))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_spM(s_in, array_in, nr_u, nc_u, nnz_u, data_u, indi_u, indj_u, 1, s_out, R1)

    allocate(R2(nr_u, nr_v, nc_w, 1))
    s_in = shape(R1); s_out = shape(R2)
    call tensor_n_mode_product_spM(s_in, R1, nr_v, nc_v, nnz_v, data_v, indi_v, indj_v, 2, s_out, R2)
    deallocate(R1)

    s_in = shape(R2); s_out = shape(array_out)
    call tensor_n_mode_product_spM(s_in, R2, nr_w, nc_w, nnz_w, data_w, indi_w, indj_w, 3, s_out, array_out)
    deallocate(R2)

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

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
    double precision, intent(in) :: array_in
    dimension :: array_in(nc_u, nc_v, nc_w, nc_t)
    double precision, intent(in) :: data_u, data_v, data_w, data_t
    dimension :: data_u(nnz_u), data_v(nnz_v), data_w(nnz_w), data_t(nnz_t)
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                    indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u, nr_v, nr_w, nr_t)

    ! Local data 
    ! ----------
    integer, parameter :: dim = 4
    double precision, allocatable, dimension(:, :, :, :) :: R1, R2, R3
    integer :: s_in(dim), s_out(dim)

    allocate(R1(nr_u, nc_v, nc_w, nc_t))
    s_in = shape(array_in); s_out = shape(R1)
    call tensor_n_mode_product_spM(s_in, array_in, nr_u, nc_u, nnz_u, data_u, indi_u, indj_u, 1, s_out, R1)

    allocate(R2(nr_u, nr_v, nc_w, nc_t))
    s_in = shape(R1); s_out = shape(R2)
    call tensor_n_mode_product_spM(s_in, R1, nr_v, nc_v, nnz_v, data_v, indi_v, indj_v, 2, s_out, R2)
    deallocate(R1)

    allocate(R3(nr_u, nr_v, nr_w, nc_t))
    s_in = shape(R2); s_out = shape(R3)
    call tensor_n_mode_product_spM(s_in, R2, nr_w, nc_w, nnz_w, data_w, indi_w, indj_w, 3, s_out, R3)
    deallocate(R2)

    s_in = shape(R3); s_out = shape(array_out)
    call tensor_n_mode_product_spM(s_in, R3, nr_t, nc_t, nnz_t, data_t, indi_t, indj_t, 4, s_out, array_out)
    deallocate(R3)

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

    use omp_lib
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

    use omp_lib
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

! ----------------------------
! Fast Diagonalization method
! ----------------------------

subroutine eigen_decomposition(nr, nc, univMcoefs, univKcoefs, nnz, indi, indj, &
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
    double precision, parameter :: penalty = 100.0d0
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
    double precision, allocatable, dimension(:, :) :: BB0, WW0, MM, BB1, WW1, KK
    
    ! Masse matrix
    data_Bt = data_B(:, 1)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univMcoefs(indj(j))
        end do
    end do

    allocate(BB0(nr, nc))
    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, BB0)
    allocate(WW0(nr, nc))
    call csr2dense(nnz, indi, indj, data_W(:, 1), nr, nc, WW0)
    allocate(MM(nr, nr))
    MM = matmul(WW0, transpose(BB0))
    deallocate(BB0, WW0)

    ! Stiffness matrix
    data_Bt = data_B(:, 2)
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_Bt(j) = data_Bt(j)*univKcoefs(indj(j))
        end do
    end do

    allocate(BB1(nr, nc))
    call csr2dense(nnz, indi, indj, data_Bt, nr, nc, BB1)
    allocate(WW1(nr, nc))
    call csr2dense(nnz, indi, indj, data_W(:, 4), nr, nc, WW1)
    allocate(KK(nr, nr))

    KK = matmul(WW1, transpose(BB1))
    deallocate(BB1, WW1)

    ! Save diagonal of M and K
    do i = 1, nr
        Kdiag(i) = KK(i, i)
        Mdiag(i) = MM(i, i)
    end do

    ! -----------------------------------
    ! Eigen decomposition KK U = MM U DD
    ! -----------------------------------
    call compute_eigdecomp_pdr(nr, KK, MM, eigval, eigvec)
    deallocate(KK, MM)

end subroutine eigen_decomposition

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
