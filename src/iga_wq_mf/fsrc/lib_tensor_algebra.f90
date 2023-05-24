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
    double precision, allocatable, dimension(:,:) :: Xt, Rt
    integer :: ju, jv, jw, jt, i, nb_tasks

    if (mode.eq.1) then 

        allocate(Xt(nc_u, nc_v*nc_w*nc_t))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(X)/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        Xt(ju, jv+(jw-1)*nc_v+(jt-1)*nc_v*nc_w) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v+(jt-1)*nc_u*nc_v*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        allocate(Rt(nr, nc_v*nc_w*nc_t))
        Rt = matmul(U, Xt)
        deallocate(Xt)

        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(R)/nb_tasks) 
        do jt =1, nc_t
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do i = 1, nr
                        R(i+(jv-1)*nr+(jw-1)*nr*nc_v+(jt-1)*nr*nc_v*nc_w) = Rt(i, jv+(jw-1)*nc_v+(jt-1)*nc_v*nc_w)
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
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

subroutine tensor_n_mode_product_spM(nc_u, nc_v, nc_w, nc_t, X, nrU, nnzU, dataU, indi, indj, mode, nrR, R)
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
    integer, intent(in) :: nc_u, nc_v, nc_w, nc_t, nrU, nnzU, mode, nrR
    integer, intent(in) :: indi, indj
    dimension :: indi(nrU+1), indj(nnzU)
    double precision, intent(in) :: X, dataU
    dimension :: X(nc_u*nc_v*nc_w*nc_t), dataU(nnzU)

    double precision, intent(out) :: R
    dimension :: R(nrR)

    ! Local data
    ! ----------
    integer :: ju, jv, jw, jt, i, jX, jR, k, nb_tasks
    double precision :: s

    if (mode.eq.1) then 
        !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nc_t*nc_w*nc_v/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do jv = 1, nc_v
                    jX = (jv-1)*nc_u + (jw-1)*nc_u*nc_v + (jt-1)*nc_u*nc_v*nc_w
                    jR = (jv-1)*nrU + (jw-1)*nrU*nc_v + (jt-1)*nrU*nc_v*nc_w
                    do i = 1, nrU
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + dataU(k) * X(indj(k) + jX)
                        end do
                        R(i + jR) = s
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.2) then 
        !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nc_t*nc_w*nc_u/nb_tasks) 
        do jt = 1, nc_t
            do jw = 1, nc_w
                do ju = 1, nc_u
                    jX = ju + (jw-1)*nc_u*nc_v + (jt-1)*nc_u*nc_v*nc_w
                    jR = ju + (jw-1)*nc_u*nrU + (jt-1)*nc_u*nrU*nc_w
                    do i = 1, nrU
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + dataU(k) * X((indj(k)-1)*nc_u + jX)
                        end do
                        R((i-1)*nc_u + jR) = s
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.3) then 
        !$OMP PARALLEL PRIVATE(jX, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nc_t*nc_v*nc_u/nb_tasks) 
        do jt = 1, nc_t
            do jv = 1, nc_v
                do ju = 1, nc_u
                    jX = ju + (jv-1)*nc_u + (jt-1)*nc_u*nc_v*nc_w
                    jR = ju + (jv-1)*nc_u + (jt-1)*nc_u*nc_v*nrU
                    do i = 1, nrU
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + dataU(k) * X((indj(k)-1)*nc_u*nc_v + jX)
                        end do
                        R((i-1)*nc_u*nc_v + jR) = s
                    end do
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.4) then 
        !$OMP PARALLEL PRIVATE(jX, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nc_w*nc_v*nc_u/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                do ju = 1, nc_u
                    jX = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                    jR = jX
                    do i = 1, nrU
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + dataU(k) * X((indj(k)-1)*nc_u*nc_v*nc_w + jX)
                        end do
                        R((i-1)*nc_u*nc_v*nc_w + jR) = s
                    end do
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
    dimension :: array_in(nc_u*nc_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data 
    ! ----------
    double precision, allocatable, dimension(:) :: R1

    allocate(R1(nr_u*nc_v))
    call tensor_n_mode_product_dM(nc_u, nc_v, 1, 1, array_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    call tensor_n_mode_product_dM(nr_u, nc_v, 1, 1, R1, nr_v, nc_v, Mv, 2, size(array_out), array_out)
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
    dimension :: array_in(nc_u*nc_v*nc_w)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v*nr_w)

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
    double precision, allocatable, dimension(:) :: R1, R2, R3

    allocate(R1(nr_u*nc_v*nc_w*nc_t))
    call tensor_n_mode_product_dM(nc_u, nc_v, nc_w, nc_t, array_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    allocate(R2(nr_u*nr_v*nc_w*nc_t))
    call tensor_n_mode_product_dM(nr_u, nc_v, nc_w, nc_t, R1, nr_v, nc_v, Mv, 2, size(R2), R2)
    deallocate(R1)

    allocate(R3(nr_u*nr_v*nr_w*nc_t))
    call tensor_n_mode_product_dM(nr_u, nr_v, nc_w, nc_t, R2, nr_w, nc_w, Mw, 3, size(R3), R3)
    deallocate(R2)   

    call tensor_n_mode_product_dM(nr_u, nr_v, nr_w, nc_t, R3, nr_t, nc_t, Mt, 4, size(array_out), array_out)
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
    dimension :: array_in(nc_u*nc_v)
    double precision, intent(in) :: data_u, data_v
    dimension :: data_u(nnz_u), data_v(nnz_v)
    integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
    dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                    indj_u(nnz_u), indj_v(nnz_v)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_u*nr_v)

    ! Local data 
    ! ----------
    double precision, allocatable, dimension(:) :: R1

    allocate(R1(nr_u*nc_v))
    call tensor_n_mode_product_spM(nc_u, nc_v, 1, 1, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    call tensor_n_mode_product_spM(nr_u, nc_v, 1, 1, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(array_out), array_out)
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
    double precision, allocatable, dimension(:) :: R1, R2

    allocate(R1(nr_u*nc_v*nc_w))
    call tensor_n_mode_product_spM(nc_u, nc_v, nc_w, 1, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    allocate(R2(nr_u*nr_v*nc_w))
    call tensor_n_mode_product_spM(nr_u, nc_v, nc_w, 1, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(R2), R2)
    deallocate(R1)

    call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, 1, R2, nr_w, nnz_w, data_w, indi_w, indj_w, 3, size(array_out), array_out)
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
    double precision, allocatable, dimension(:) :: R1, R2, R3

    allocate(R1(nr_u*nc_v*nc_w*nc_t))
    call tensor_n_mode_product_spM(nc_u, nc_v, nc_w, nc_t, array_in, nr_u, nnz_u, data_u, indi_u, indj_u, 1, size(R1), R1)

    allocate(R2(nr_u*nr_v*nc_w*nc_t))
    call tensor_n_mode_product_spM(nr_u, nc_v, nc_w, nc_t, R1, nr_v, nnz_v, data_v, indi_v, indj_v, 2, size(R2), R2)
    deallocate(R1)

    allocate(R3(nr_u*nr_v*nr_w*nc_t))
    call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, nc_t, R2, nr_w, nnz_w, data_w, indi_w, indj_w, 3, size(R3), R3)
    deallocate(R2)

    call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, nc_t, R3, nr_t, nnz_t, data_t, indi_t, indj_t, 4, size(array_out), array_out)
    deallocate(R3)

end subroutine sumfacto4d_spM

! -------------------------------------
! Sum factorization to compute matrices 
! -------------------------------------

subroutine csr_get_row_2d(coefs, nr_u, nc_u, nr_v, nc_v, row_u, row_v, &
                        nnz_row_u, nnz_row_v, i_nnz_u, i_nnz_v, &
                        nnz_col_u, nnz_col_v, j_nnz_u, j_nnz_v, &
                        basis_u, basis_v, weights_u, weights_v, data_row)
    !! Computes a row of a matrix using row_u and row_v

    implicit none 
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, row_u, row_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v)
    integer, intent(in) :: nnz_row_u, nnz_row_v, nnz_col_u, nnz_col_v
    integer, intent(in) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
    dimension ::    i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), &
                    j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v)
    double precision, intent(in) :: basis_u, basis_v, weights_u, weights_v
    dimension ::    basis_u(nr_u, nc_u), weights_u(nr_u, nc_u), &   
                    basis_v(nr_v, nc_v), weights_v(nr_v, nc_v)                         

    double precision, intent(out) :: data_row
    dimension :: data_row(nnz_row_u*nnz_row_v)

    ! Local data 
    ! ----------
    integer :: iu, iv, ju, jv, posu, posv, posCoef, genPosC
    double precision, allocatable, dimension(:) :: Ci0, Wt
    double precision, allocatable, dimension(:,:) :: BW_u, BW_v, Bt

    ! Set values of C
    allocate(Ci0(nnz_col_u*nnz_col_v))
    do jv = 1, nnz_col_v
        do ju = 1, nnz_col_u
            genPosC = ju + (jv-1)*nnz_col_u 
            posu = j_nnz_u(ju)
            posv = j_nnz_v(jv)
            posCoef = posu + (posv-1)*nc_u 
            Ci0(genPosC) = coefs(posCoef)                    
        end do
    end do

    ! Set values of BW in u-dimension
    allocate(BW_u(nnz_row_u, nnz_col_u), Bt(nnz_row_u, nnz_col_u), Wt(nnz_col_u))
    Bt = basis_u(i_nnz_u, j_nnz_u)
    Wt = weights_u(row_u, j_nnz_u)
    
    forall (iu = 1 : size(Bt, 1), ju = 1 : size(Bt, 2)) 
        BW_u(iu, ju) = Bt(iu, ju) * Wt(ju)
    end forall
    deallocate(Bt, Wt)

    ! Set values of BW in v-dimension
    allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
    Bt = basis_v(i_nnz_v, j_nnz_v)
    Wt = weights_v(row_v, j_nnz_v)

    forall (iv = 1 : size(Bt, 1), jv = 1 : size(Bt, 2)) 
        BW_v(iv, jv) = Bt(iv, jv) * Wt(jv)
    end forall
    deallocate(Bt, Wt)

    ! Compute non zero values of the row
    call sumfacto2d_dM(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, BW_u, BW_v, Ci0, data_row)

end subroutine csr_get_row_2d 

subroutine csr_get_matrix_2d(coefs, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, &
                            data_B_u, data_B_v, data_W_u, data_W_v, &
                            nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                            indi_result, nnz_result, data_result)
    !! Computes a matrix in 2D case (Wv . Bv) x (Wu . Bu)
    !! x: kronecker product and .: inner product
    !! Indices must be in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, nnz_I_u, nnz_I_v
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
    dimension ::    data_B_u(nnz_u), data_W_u(nnz_u), &
                    data_B_v(nnz_v), data_W_v(nnz_v)
    integer, intent(in) :: indi_I_u, indi_I_v, indj_I_u, indj_I_v 
    dimension ::    indi_I_u(nr_u+1), indi_I_v(nr_v+1), &
                    indj_I_u(nnz_I_u), indj_I_v(nnz_I_v)
    integer, intent(in) :: indi_result
    dimension :: indi_result(nr_u*nr_v+1) 
    integer, intent(in) :: nnz_result

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_result)

    ! Local data 
    ! ----------
    integer :: iu, iv, ju, jv, k, genPos, offset, nb_tasks
    double precision, allocatable, dimension(:,:) :: B_u, B_v, W_u, W_v
    
    integer :: nnz_row_u, nnz_row_v, nnz_col_u, nnz_col_v, nnz_gen_row
    integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
    double precision, allocatable, dimension(:) :: data_row

    allocate(B_u(nr_u, nc_u), B_v(nr_v, nc_v), W_u(nr_u, nc_u), W_v(nr_v, nc_v))
    call csr2dense(nnz_u, indi_u, indj_u, data_B_u, nr_u, nc_u, B_u)
    call csr2dense(nnz_v, indi_v, indj_v, data_B_v, nr_v, nc_v, B_v)
    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W_u)
    call csr2dense(nnz_v, indi_v, indj_v, data_W_v, nr_v, nc_v, W_v)

    !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,offset,j_nnz_u,ju,j_nnz_v,jv,nnz_row_u) &
    !$OMP PRIVATE(nnz_row_v,i_nnz_u,i_nnz_v,data_row,nnz_gen_row,genPos,k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nr_u*nr_v/nb_tasks)
    do iv = 1, nr_v
        do iu = 1, nr_u

            ! Set values of rows
            nnz_row_u = indi_I_u(iu+1) - indi_I_u(iu)
            allocate(i_nnz_u(nnz_row_u))
            offset = indi_I_u(iu)
            do ju = 1, nnz_row_u
                i_nnz_u(ju) = indj_I_u(ju+offset-1)
            end do
            
            nnz_row_v = indi_I_v(iv+1) - indi_I_v(iv)
            allocate(i_nnz_v(nnz_row_v))
            offset = indi_I_v(iv)
            do jv = 1, nnz_row_v
                i_nnz_v(jv) = indj_I_v(jv+offset-1)
            end do

            ! Set values of columns
            nnz_col_u = indi_u(iu+1) - indi_u(iu)
            allocate(j_nnz_u(nnz_col_u))
            offset = indi_u(iu)
            do ju = 1, nnz_col_u
                j_nnz_u(ju) = indj_u(ju+offset-1)
            end do

            nnz_col_v = indi_v(iv+1) - indi_v(iv)
            allocate(j_nnz_v(nnz_col_v))
            offset = indi_v(iv)
            do jv = 1, nnz_col_v
                j_nnz_v(jv) = indj_v(jv+offset-1)
            end do

            ! Get data of row
            nnz_gen_row = nnz_row_u*nnz_row_v 
            allocate(data_row(nnz_gen_row))
            call csr_get_row_2d(coefs, nr_u, nc_u, nr_v, nc_v, iu, iv, &
                            nnz_row_u, nnz_row_v, i_nnz_u, i_nnz_v, &
                            nnz_col_u, nnz_col_v, j_nnz_u, j_nnz_v, &
                            B_u, B_v, W_u, W_v, data_row)
            deallocate(i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v)

            ! Get offset in result 
            genPos = iu + (iv - 1)*nr_u 
            offset = indi_result(genPos)
        
            ! Get result
            do k = 1, nnz_gen_row
                data_result(offset+k-1) = data_row(k)
            end do   

            deallocate(data_row)
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    deallocate(B_u, B_v, W_u, W_v)

end subroutine csr_get_matrix_2d 

subroutine csr_get_row_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, row_u, row_v, row_w, &
                        nnz_row_u, nnz_row_v, nnz_row_w, i_nnz_u, i_nnz_v, i_nnz_w, &
                        nnz_col_u, nnz_col_v, nnz_col_w, j_nnz_u, j_nnz_v, j_nnz_w, &
                        basis_u, basis_v, basis_w, weights_u, weights_v, weights_w, data_row)
    !! Computes a row of a matrix using row_u, row_v and row_w

    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, row_u, row_v, row_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v*nc_w)
    integer, intent(in) :: nnz_row_u, nnz_row_v, nnz_row_w, nnz_col_u, nnz_col_v, nnz_col_w
    integer, intent(in) :: i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w
    dimension ::    i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), i_nnz_w(nnz_row_w), &
                    j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v), j_nnz_w(nnz_col_w)
    double precision, intent(in) :: basis_u, basis_v, basis_w, weights_u, weights_v, weights_w
    dimension ::    basis_u(nr_u, nc_u), weights_u(nr_u, nc_u), &   
                    basis_v(nr_v, nc_v), weights_v(nr_v, nc_v), &
                    basis_w(nr_w, nc_w), weights_w(nr_w, nc_w)

    double precision, intent(out) :: data_row
    dimension :: data_row(nnz_row_u*nnz_row_v*nnz_row_w)

    ! Local data 
    ! ----------
    integer :: iu, iv, iw, ju, jv, jw, posu, posv, posw, posCoef, genPosC
    double precision, allocatable, dimension(:) :: Ci0, Wt
    double precision, allocatable, dimension(:,:) :: BW_u, BW_v, BW_w, Bt

    ! Set values of C
    allocate(Ci0(nnz_col_u*nnz_col_v*nnz_col_w))
    do jw = 1, nnz_col_w
        do jv = 1, nnz_col_v
            do ju = 1, nnz_col_u
                genPosC = ju + (jv-1)*nnz_col_u + (jw-1)*nnz_col_u*nnz_col_v
                posu = j_nnz_u(ju)
                posv = j_nnz_v(jv)
                posw = j_nnz_w(jw)
                posCoef = posu + (posv-1)*nc_u + (posw-1)*nc_u*nc_v
                Ci0(genPosC) = coefs(posCoef)                    
            end do
        end do
    end do

    ! Set values of BW in u-dimension
    allocate(BW_u(nnz_row_u,nnz_col_u), Bt(nnz_row_u,nnz_col_u), Wt(nnz_col_u))
    Bt = basis_u(i_nnz_u, j_nnz_u)
    Wt = weights_u(row_u, j_nnz_u)
    
    forall (iu = 1 : size(Bt, 1), ju = 1 : size(Bt, 2)) 
        BW_u(iu, ju) = Bt(iu, ju) * Wt(ju)
    end forall
    deallocate(Bt, Wt)

    ! Set values of BW in v-dimension
    allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
    Bt = basis_v(i_nnz_v, j_nnz_v)
    Wt = weights_v(row_v, j_nnz_v)
    
    forall (iv = 1 : size(Bt, 1), jv = 1 : size(Bt, 2)) 
        BW_v(iv, jv) = Bt(iv, jv) * Wt(jv)
    end forall
    deallocate(Bt, Wt)

    ! Set values of BW in w-dimension
    allocate(BW_w(nnz_row_w,nnz_col_w), Bt(nnz_row_w,nnz_col_w), Wt(nnz_col_w))
    Bt = basis_w(i_nnz_w, j_nnz_w)
    Wt = weights_w(row_w, j_nnz_w)

    forall (iw = 1 : size(Bt, 1), jw = 1 : size(Bt, 2)) 
        BW_w(iw, jw) = Bt(iw, jw) * Wt(jw)
    end forall
    deallocate(Bt, Wt)

    ! Compute non zero values of the row
    call sumfacto3d_dM(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, &
                            nnz_row_w, nnz_col_w, BW_u, BW_v, BW_w, Ci0, data_row)

end subroutine csr_get_row_3d 

subroutine csr_get_matrix_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                            nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                            indi_result, nnz_result, data_result)
    !! Computes a matrix in 3D case (Ww . Bw) x (Wv . Bv) x (Wu . Bu)
    !! x: kronecker product and .: inner product
    !! Indices must be in CSR format

    use omp_lib
    implicit none 
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, nnz_I_u, nnz_I_v, nnz_I_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v*nc_w)
    integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
    dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                    indi_v(nr_v+1), indj_v(nnz_v), &
                    indi_w(nr_w+1), indj_w(nnz_w)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
    dimension ::    data_B_u(nnz_u), data_W_u(nnz_u), &
                    data_B_v(nnz_v), data_W_v(nnz_v), &
                    data_B_w(nnz_w), data_W_w(nnz_w)
    integer, intent(in) :: indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w 
    dimension ::    indi_I_u(nr_u+1), indi_I_v(nr_v+1), indi_I_w(nr_w+1), &
                    indj_I_u(nnz_I_u), indj_I_v(nnz_I_v), indj_I_w(nnz_I_w)
    integer, intent(in) :: indi_result
    dimension :: indi_result(nr_u*nr_v*nr_w+1) 
    integer, intent(in) :: nnz_result

    double precision, intent(out) :: data_result
    dimension :: data_result(nnz_result)

    ! Local data 
    ! ----------
    integer :: iu, iv, iw, ju, jv, jw, k, genPos, offset, nb_tasks
    double precision, allocatable, dimension(:,:) :: B_u, B_v, B_w, W_u, W_v, W_w
    
    integer :: nnz_row_u, nnz_row_v, nnz_row_w, nnz_col_u, nnz_col_v, nnz_col_w, nnz_gen_row
    integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w
    double precision, allocatable :: data_row(:)

    allocate(B_u(nr_u, nc_u), B_v(nr_v, nc_v), B_w(nr_w, nc_w), &
            W_u(nr_u, nc_u), W_v(nr_v, nc_v), W_w(nr_w, nc_w))
    call csr2dense(nnz_u, indi_u, indj_u, data_B_u, nr_u, nc_u, B_u)
    call csr2dense(nnz_v, indi_v, indj_v, data_B_v, nr_v, nc_v, B_v)
    call csr2dense(nnz_w, indi_w, indj_w, data_B_w, nr_w, nc_w, B_w)
    call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W_u)
    call csr2dense(nnz_v, indi_v, indj_v, data_W_v, nr_v, nc_v, W_v)
    call csr2dense(nnz_w, indi_w, indj_w, data_W_w, nr_w, nc_w, W_w)

    !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,nnz_col_w,offset,j_nnz_u,ju,j_nnz_v,jv,j_nnz_w,jw,nnz_row_u) &
    !$OMP PRIVATE(nnz_row_v,nnz_row_w,i_nnz_u,i_nnz_v,i_nnz_w,data_row,nnz_gen_row,genPos,k)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_u*nr_v*nr_w/nb_tasks)
    do iw = 1, nr_w
        do iv = 1, nr_v
            do iu = 1, nr_u
                
                ! Set values of rows
                nnz_row_u = indi_I_u(iu+1) - indi_I_u(iu)
                allocate(i_nnz_u(nnz_row_u))
                offset = indi_I_u(iu)
                do ju = 1, nnz_row_u
                    i_nnz_u(ju) = indj_I_u(ju+offset-1)
                end do

                nnz_row_v = indi_I_v(iv+1) - indi_I_v(iv)
                allocate(i_nnz_v(nnz_row_v))
                offset = indi_I_v(iv)
                do jv = 1, nnz_row_v
                    i_nnz_v(jv) = indj_I_v(jv+offset-1)
                end do

                nnz_row_w = indi_I_w(iw+1) - indi_I_w(iw)
                allocate(i_nnz_w(nnz_row_w))
                offset = indi_I_w(iw)
                do jw = 1, nnz_row_w
                    i_nnz_w(jw) = indj_I_w(jw+offset-1)
                end do

                ! Set values of columns
                nnz_col_u = indi_u(iu+1) - indi_u(iu)
                allocate(j_nnz_u(nnz_col_u))
                offset = indi_u(iu)
                do ju = 1, nnz_col_u
                    j_nnz_u(ju) = indj_u(ju+offset-1)
                end do

                nnz_col_v = indi_v(iv+1) - indi_v(iv)
                allocate(j_nnz_v(nnz_col_v))
                offset = indi_v(iv)
                do jv = 1, nnz_col_v
                    j_nnz_v(jv) = indj_v(jv+offset-1)
                end do

                nnz_col_w = indi_w(iw+1) - indi_w(iw)
                allocate(j_nnz_w(nnz_col_w))
                offset = indi_w(iw)
                do jw = 1, nnz_col_w
                    j_nnz_w(jw) = indj_w(jw+offset-1)
                end do

                ! Get data of row
                nnz_gen_row = nnz_row_u*nnz_row_v*nnz_row_w
                allocate(data_row(nnz_gen_row))
                call csr_get_row_3d(coefs, nr_u, nc_u,  nr_v, nc_v, nr_w, nc_w, iu, iv, iw, &
                                nnz_row_u, nnz_row_v, nnz_row_w, i_nnz_u, i_nnz_v, i_nnz_w, &
                                nnz_col_u, nnz_col_v, nnz_col_w, j_nnz_u, j_nnz_v, j_nnz_w, &
                                B_u, B_v, B_w, W_u, W_v, W_w, data_row)
                deallocate(i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w)

                ! Get offset in result 
                genPos = iu + (iv - 1)*nr_u + (iw - 1)*nr_u*nr_v
                offset = indi_result(genPos)
            
                ! Get result
                do k = 1, nnz_gen_row
                    data_result(offset+k-1) = data_row(k)
                end do   

                deallocate(data_row)
            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    deallocate(B_u, B_v, B_w, W_u, W_v, W_w)

end subroutine csr_get_matrix_3d 

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

subroutine eigen_decomposition(nr, nc, Mcoefs, Kcoefs, nnz, indi, indj, &
                                data_B0, data_W0, data_B1, data_W1, robin_condition, &
                                eigenvalues, eigenvectors, Kdiag, Mdiag)
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
    double precision, intent(in) :: Mcoefs, Kcoefs
    dimension :: Mcoefs(nc), Kcoefs(nc)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)
    double precision, intent(in) :: data_B0, data_W0, data_B1, data_W1
    dimension :: data_B0(nnz), data_W0(nnz), data_B1(nnz), data_W1(nnz)
    integer, intent(in) :: robin_condition
    dimension :: robin_condition(2)
            
    double precision, intent(out) :: eigenvalues, eigenvectors
    dimension :: eigenvalues(nr), eigenvectors(nr, nr)
    double precision, intent(out) :: Kdiag, Mdiag
    dimension :: Kdiag(nr), Mdiag(nr)

    ! Local data
    ! ----------
    integer :: i, j
    double precision, allocatable, dimension(:) :: data_B0t, data_B1t
    double precision, allocatable, dimension(:,:) :: BB0, WW0, MM, BB1, WW1, KK
    
    ! Masse matrix
    allocate(data_B0t(nnz))
    data_B0t = data_B0
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_B0t(j) = data_B0t(j)*Mcoefs(indj(j))
        end do
    end do

    allocate(BB0(nr, nc))
    call csr2dense(nnz, indi, indj, data_B0t, nr, nc, BB0)
    deallocate(data_B0t)
    allocate(WW0(nr, nc))
    call csr2dense(nnz, indi, indj, data_W0, nr, nc, WW0)
    allocate(MM(nr, nr))
    MM = matmul(WW0, transpose(BB0))
    deallocate(BB0, WW0)

    ! Stiffness matrix
    allocate(data_B1t(nnz))
    data_B1t = data_B1
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_B1t(j) = data_B1t(j)*Kcoefs(indj(j))
        end do
    end do

    allocate(BB1(nr, nc))
    call csr2dense(nnz, indi, indj, data_B1t, nr, nc, BB1)
    deallocate(data_B1t)
    allocate(WW1(nr, nc))
    call csr2dense(nnz, indi, indj, data_W1, nr, nc, WW1)
    allocate(KK(nr, nr))
    KK = matmul(WW1, transpose(BB1))
    deallocate(BB1, WW1)

    ! Modify K to avoid singular matrix (using Robin boundary condition)
    if (robin_condition(1).eq.1) then 
        KK(1, 1) = penalty*KK(1, 1)
    end if

    if (robin_condition(2).eq.1) then 
        KK(nr, nr) = penalty*KK(nr, nr)
    end if

    ! Save diagonal of M and K
    do i = 1, nr
        Kdiag(i) = KK(i, i)
        Mdiag(i) = MM(i, i)
    end do

    ! -----------------------------------
    ! Eigen decomposition KK U = MM U DD
    ! -----------------------------------
    call compute_geneigs(nr, KK, MM, eigenvalues, eigenvectors)
    deallocate(KK, MM)

end subroutine eigen_decomposition

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
