! =========================================
! module :: sum product (Tensor operations)
! author :: Joaquin Cornejo
! =========================================
module tensor_methods

    contains

    ! ----------------------------------------------------
    ! Tensor algebra 
    ! ----------------------------------------------------

    subroutine tensor_n_mode_product(nc_u, nc_v, nc_w, X, nr, nc, U, n, nu, nv, nw, R)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = (nc_u, nc_v, nc_w)
        !! Matrix U = (nr, nc)
        !! Tensor R = (nu, nv, nw) (It depends on 'n')
        !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u
        !! THESE FUNCTIONS IS NOT OPTIMIZED AND IT COULD BE BETTER

        use omp_lib
        implicit none
        ! Input / output data 
        ! -------------------- 
        integer, intent(in) :: nc_u, nc_v, nc_w, nr, nc, n, nu, nv, nw
        double precision, intent(in) :: X, U
        dimension :: X(nc_u*nc_v*nc_w), U(nr, nc)

        double precision, intent(out) :: R
        dimension :: R(nu*nv*nw)

        ! Local data
        ! ---------------
        double precision, allocatable, dimension(:, :) :: Xt, Rt
        integer :: ju, jv, jw, i, nb_tasks

        if (n.eq.1) then 

            allocate(Xt(nc_u, nc_v*nc_w))
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nc_u/nb_tasks) 
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        Xt(ju, jv+(jw-1)*nc_v) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

            allocate(Rt(nr, nc_v*nc_w))
            Rt = matmul(U, Xt)
            deallocate(Xt)

            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nr/nb_tasks) 
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do i = 1, nr
                        R(i+(jv-1)*nr+(jw-1)*nr*nc_v) = Rt(i, jv+(jw-1)*nc_v)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            deallocate(Rt)

        else if (n.eq.2) then 

            allocate(Xt(nc_v, nc_u*nc_w))
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nc_u/nb_tasks) 
            do jw = 1, nc_w
                do ju = 1, nc_u
                    do jv = 1, nc_v
                        Xt(jv, ju+(jw-1)*nc_u) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

            allocate(Rt(nr, nc_u*nc_w))
            Rt = matmul(U, Xt)
            deallocate(Xt)

            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nr*nc_u/nb_tasks) 
            do jw = 1, nc_w
                do ju = 1, nc_u
                    do i = 1, nr
                        R(ju+(i-1)*nc_u+(jw-1)*nc_u*nr) = Rt(i, ju+(jw-1)*nc_u)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            deallocate(Rt)
            
        else if (n.eq.3) then 

            allocate(Xt(nc_w, nc_u*nc_v))
            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nc_u/nb_tasks) 
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do jw = 1, nc_w
                        Xt(jw, ju+(jv-1)*nc_u) = X(ju+(jv-1)*nc_u+(jw-1)*nc_u*nc_v)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

            allocate(Rt(nr, nc_u*nc_v))
            Rt = matmul(U, Xt)
            deallocate(Xt)

            !$OMP PARALLEL 
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nr*nc_v*nc_u/nb_tasks) 
            do jv = 1, nc_v
                do ju = 1, nc_u
                    do i = 1, nr
                        R(ju+(jv-1)*nc_u+(i-1)*nc_u*nc_v) = Rt(i, ju+(jv-1)*nc_u)
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            deallocate(Rt)
            
        end if

    end subroutine tensor_n_mode_product

    subroutine tensor_n_mode_product_sp(nc_u, nc_v, nc_w, X, nr, nc, nnz, U, indi, indj, n, nu, nv, nw, R)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = (nc_u, nc_v, nc_w)
        !! Matrix U = (nr, nc)
        !! Tensor R = (nu, nv, nw) (It depends on 'n')
        !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u
        !! U is in CSR format

        use omp_lib
        implicit none
        ! Input / output data 
        ! -------------------- 
        integer, intent(in) :: nc_u, nc_v, nc_w, nr, nc, n, nu, nv, nw, nnz
        double precision, intent(in) :: X, U
        dimension :: X(nc_u*nc_v*nc_w), U(nnz)
        integer, intent(in) :: indi, indj
        dimension :: indi(nr+1), indj(nnz)

        double precision, intent(out) ::  R
        dimension :: R(nu*nv*nw)

        ! Local data
        ! ---------------
        integer :: ju, jv, jw, i, jX, jR, k, dummy, nb_tasks
        double precision :: s

        ! Initialize
        dummy = nc

        if (n.eq.1) then 
            !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_w*nc_v/nb_tasks) 
            do jw = 1, nc_w
                do jv = 1, nc_v
                    jX = (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                    jR = (jv-1)*nr + (jw-1)*nr*nc_v
                    do i = 1, nr
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + U(k) * X(indj(k) + jX)
                        end do
                        R(i + jR) = s
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        else if (n.eq.2) then 
            !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_w*nc_u/nb_tasks) 
            do jw = 1, nc_w
                do ju = 1, nc_u
                    jX = ju + (jw-1)*nc_u*nc_v
                    jR = ju + (jw-1)*nc_u*nr
                    do i = 1, nr
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + U(k) * X((indj(k)-1)*nc_u + jX)
                        end do
                        R((i-1)*nc_u + jR) = s
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        else if (n.eq.3) then 
            !$OMP PARALLEL PRIVATE(jX, i, s, k)
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_v*nc_u/nb_tasks) 
            do jv = 1, nc_v
                do ju = 1, nc_u
                    ! In this case jR = jX
                    jX = ju + (jv-1)*nc_u
                    do i = 1, nr
                        s = 0.d0
                        do k = indi(i), indi(i+1)-1
                            s = s + U(k) * X((indj(k)-1)*nc_u*nc_v + jX)
                        end do
                        R((i-1)*nc_u*nc_v + jX) = s
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            
        end if

    end subroutine tensor_n_mode_product_sp

    subroutine tensor2d_dot_vector(nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, &
                                    Mu, Mv, vector_in, vector_out)
        !! Evaluates a dot product between a tensor 2D and a vector 
        !! Based on "Preconditioners for IGA" by Montardini
        !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v)
        double precision, intent(in) :: Mu, Mv
        dimension :: Mu(nb_rows_u, nb_cols_u), Mv(nb_rows_v, nb_cols_v)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1
    
        ! First product
        allocate(R1(nb_rows_u*nb_cols_v))
        call tensor_n_mode_product(nb_cols_u, nb_cols_v, 1, vector_in, &
        nb_rows_u, nb_cols_u, Mu, 1, nb_rows_u, nb_cols_v, 1, R1)

        ! Second product
        call tensor_n_mode_product(nb_rows_u, nb_cols_v, 1, R1, &
        nb_rows_v, nb_cols_v, Mv, 2, nb_rows_u, nb_rows_v, 1, vector_out)
        deallocate(R1)

    end subroutine tensor2d_dot_vector

    subroutine tensor3d_dot_vector(nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, &
                                    nb_rows_w, nb_cols_w, &
                                    Mu, Mv, Mw, vector_in, vector_out)
        !! Evaluates a dot product between a tensor 3D and a vector 
        !! Based on "Preconditioners for IGA" by Montardini
        !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Matrix Mw = (nb_rows_w, nb_cols_w)
        !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v,nb_rows_w, nb_cols_w
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v*nb_cols_w)
        double precision, intent(in) :: Mu, Mv, Mw
        dimension :: Mu(nb_rows_u, nb_cols_u), Mv(nb_rows_v, nb_cols_v), Mw(nb_rows_w, nb_cols_w)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data 
        ! -------------
        integer :: u, v, w, nb_tasks
        double precision, allocatable, dimension(:, :) :: X, R

        ! First product
        allocate(X(nb_cols_u, nb_cols_v*nb_cols_w))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(X)/nb_tasks)
        do w = 1, nb_cols_w
            do v = 1, nb_cols_v
                do u = 1, nb_cols_u
                    X(u, v+(w-1)*nb_cols_v) = vector_in(u+(v-1)*nb_cols_u+(w-1)*nb_cols_u*nb_cols_v)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        allocate(R(nb_rows_u, nb_cols_v*nb_cols_w))
        R = matmul(Mu, X)
        deallocate(X)

        ! Second product
        allocate(X(nb_cols_v, nb_rows_u*nb_cols_w))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(R)/nb_tasks)
        do w = 1, nb_cols_w
            do u = 1, nb_rows_u
                do v = 1, nb_cols_v
                    X(v, u+(w-1)*nb_rows_u) = R(u, v+(w-1)*nb_cols_v)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(R)

        allocate(R(nb_rows_v, nb_rows_u*nb_cols_w))
        R = matmul(Mv, X)
        deallocate(X)

        ! Third product
        allocate(X(nb_cols_w, nb_rows_u*nb_rows_v))
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(R)/nb_tasks)
        do w = 1, nb_cols_w 
            do u = 1, nb_rows_u
                do v = 1, nb_rows_v
                    X(w, u+(v-1)*nb_rows_u) = R(v, u+(w-1)*nb_rows_u)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(R)
        
        allocate(R(nb_rows_w, nb_rows_u*nb_rows_v))
        R = matmul(Mw, X)
        deallocate(X)

        ! Re-arrange output
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size(R)/nb_tasks)
        do v = 1, nb_rows_v
            do u = 1, nb_rows_u
                do w = 1, nb_rows_w
                    vector_out(u+(v-1)*nb_rows_u+(w-1)*nb_rows_u*nb_rows_v) = R(w, u+(v-1)*nb_rows_u)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(R)

    end subroutine tensor3d_dot_vector

    subroutine tensor2d_dot_vector_sp(nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, &
                                    size_data_u, indi_u, indj_u, data_u, &
                                    size_data_v, indi_v, indj_v, data_v, &
                                    vector_in, vector_out)
        !! Evaluates a dot product between a tensor 3D and a vector 
        !! Based on "Preconditioners for IGA" by Montardini
        !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Vector_in = (nb_cols_u * nb_cols_v)

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, &
                                nb_cols_v, size_data_u, size_data_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v)
        double precision, intent(in) :: data_u, data_v
        dimension :: data_u(size_data_u), data_v(size_data_v)
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), &
                        indj_u(size_data_u), indj_v(size_data_v)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1

        ! First product
        allocate(R1(nb_rows_u*nb_cols_v))
        call tensor_n_mode_product_sp(nb_cols_u, nb_cols_v, 1, vector_in, &
        nb_rows_u, nb_cols_u, size_data_u, data_u, indi_u, indj_u, 1, nb_rows_u, nb_cols_v, 1, R1)

        ! Second product
        call tensor_n_mode_product_sp(nb_rows_u, nb_cols_v, 1, R1, &
        nb_rows_v, nb_cols_v, size_data_v, data_v, indi_v, indj_v, 2, nb_rows_u, nb_rows_v, 1, vector_out)
        deallocate(R1)

    end subroutine tensor2d_dot_vector_sp

    subroutine tensor3d_dot_vector_sp(nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                    size_data_u, indi_u, indj_u, data_u, &
                                    size_data_v, indi_v, indj_v, data_v, &
                                    size_data_w, indi_w, indj_w, data_w, &
                                    vector_in, vector_out)
        !! Evaluates a dot product between a tensor 3D and a vector 
        !! Based on "Preconditioners for IGA" by Montardini
        !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Matrix Mw = (nb_rows_w, nb_cols_w)
        !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, &
                                nb_cols_v,nb_rows_w, nb_cols_w, size_data_u, size_data_v, size_data_w
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v*nb_cols_w)
        double precision, intent(in) :: data_u, data_v, data_w
        dimension :: data_u(size_data_u), data_v(size_data_v), data_w(size_data_w)
        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), indi_w(nb_rows_w+1), &
                        indj_u(size_data_u), indj_v(size_data_v), indj_w(size_data_w)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1, R2
        integer :: u, v, w, pX, pR1, pR2, k, nb_tasks
        double precision :: s

        ! First product
        allocate(R1(nb_rows_u*nb_cols_v*nb_cols_w))
        !$OMP PARALLEL PRIVATE(pX, pR1, u, v, w, k, s)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nb_cols_v*nb_cols_w/nb_tasks)
        do w = 1, nb_cols_w
            do v = 1, nb_cols_v
                pX = (v-1)*nb_cols_u + (w-1)*nb_cols_u*nb_cols_v
                pR1 = v + (w-1)*nb_rows_u*nb_cols_v
                do u = 1, nb_rows_u
                    s = 0.d0
                    do k = indi_u(u), indi_u(u+1)-1
                        s = s + data_u(k) * vector_in(indj_u(k) + pX)
                    end do
                    R1((u-1)*nb_cols_v + pR1) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        ! Second product
        allocate(R2(nb_rows_u*nb_rows_v*nb_cols_w))
        !$OMP PARALLEL PRIVATE(pR1, pR2, u, v, w, k, s)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nb_rows_u*nb_cols_w/nb_tasks)
        do w = 1, nb_cols_w
            do u = 1, nb_rows_u
                pR1 = (u-1)*nb_cols_v + (w-1)*nb_rows_u*nb_cols_v
                pR2 = w + (u-1)*nb_cols_w
                do v = 1, nb_rows_v
                    s = 0.d0
                    do k = indi_v(v), indi_v(v+1)-1
                        s = s + data_v(k) * R1(indj_v(k) + pR1)
                    end do
                    R2((v-1)*nb_rows_u*nb_cols_w + pR2) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(R1)

        !$OMP PARALLEL PRIVATE(pX, pR2, u, v, w, k, s)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nb_rows_u*nb_rows_v/nb_tasks)
        ! Third product
        do v = 1, nb_rows_v
            do u = 1, nb_rows_u
                pR2 = (u-1)*nb_cols_w + (v-1)*nb_rows_u*nb_cols_w
                pX = u + (v-1)*nb_rows_u
                do w = 1, nb_rows_w
                    s = 0.d0
                    do k = indi_w(w), indi_w(w+1)-1
                        s = s + data_w(k) * R2(indj_w(k) + pR2)
                    end do
                    vector_out((w-1)*nb_rows_u*nb_rows_v + pX) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        deallocate(R2)

    end subroutine tensor3d_dot_vector_sp

    ! ----------------------------------------------------
    ! Sum product to compute matrices 
    ! ----------------------------------------------------

    subroutine csr_get_row_2d(coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            row_u, row_v, &
                            nnz_row_u, nnz_row_v, &
                            i_nnz_u, i_nnz_v, &
                            nnz_col_u, nnz_col_v, &
                            j_nnz_u, j_nnz_v, &
                            B_u, B_v, W_u, W_v, &
                            data_row)
        !! Computes a row of a matrix constructed with row_u and row_v

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v)
        integer, intent(in) :: row_u, row_v
        integer, intent(in) :: nnz_row_u, nnz_row_v
        integer, intent(in) :: nnz_col_u, nnz_col_v
        integer, intent(in) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
        dimension :: i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), &
                    j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v)
        double precision, intent(in) :: B_u, B_v, W_u, W_v
        dimension ::    B_u(nb_rows_u, nb_cols_u), &   
                        B_v(nb_rows_v, nb_cols_v), &
                        W_u(nb_rows_u, nb_cols_u), &
                        W_v(nb_rows_v, nb_cols_v)

        double precision, intent(out) :: data_row
        dimension :: data_row(nnz_row_u*nnz_row_v)

        ! Local data 
        ! ----------------- 
        ! Create Ci
        integer :: pos_coef
        double precision, allocatable :: Ci0(:)

        ! Create Bl and Wl
        integer :: i, j
        double precision, allocatable, dimension(:,:) :: BW_u, BW_v, Bt
        double precision, allocatable, dimension(:) :: Wt

        ! Loops
        integer :: genPosC
        integer :: ju, jv, posu, posv

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v))

        ! Set values of C
        do jv = 1, nnz_col_v
            do ju = 1, nnz_col_u
                posu = j_nnz_u(ju)
                posv = j_nnz_v(jv)
                pos_coef = posu + (posv-1)*nb_cols_u 
                genPosC = ju + (jv-1)*nnz_col_u 
                Ci0(genPosC) = coefs(pos_coef)                    
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u,nnz_col_u), Bt(nnz_row_u,nnz_col_u), Wt(nnz_col_u))
        Bt = B_u(i_nnz_u, j_nnz_u)
        Wt = W_u(row_u, j_nnz_u)
        
        forall (i = 1 : size(Bt, 1), j = 1 : size(Bt, 2)) 
            BW_u(i, j) = Bt(i, j) * Wt(j)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
        Bt = B_v(i_nnz_v, j_nnz_v)
        Wt = W_v(row_v, j_nnz_v)

        forall (i = 1 : size(Bt, 1), j = 1 : size(Bt, 2)) 
            BW_v(i, j) = Bt(i, j) * Wt(j)
        end forall
        deallocate(Bt, Wt)

        call tensor2d_dot_vector(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, &
                                BW_u, BW_v, Ci0, data_row)

    end subroutine csr_get_row_2d 

    subroutine csr_get_matrix_2d(coefs, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                size_data_u, size_data_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, &
                                size_data_I_u, size_data_I_v, &
                                indi_I_u, indi_I_v, & 
                                indj_I_u, indj_I_v, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 2D case (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v)
        integer, intent(in) :: size_data_u, size_data_v
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
        dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                        indi_v(nb_rows_v+1), indj_v(size_data_v)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
        dimension ::    data_B_u(size_data_u), data_W_u(size_data_u), &
                        data_B_v(size_data_v), data_W_v(size_data_v)
        integer, intent(in) :: size_data_I_u, size_data_I_v
        integer, intent(in) :: indi_I_u, indi_I_v, &
                                indj_I_u, indj_I_v 
        dimension ::    indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), &
                        indj_I_u(size_data_I_u), indj_I_v(size_data_I_v)
        integer, intent(in) :: indi_result
        dimension :: indi_result(nb_rows_u*nb_rows_v+1) 
        integer, intent(in) :: size_data_result

        double precision, intent(out) :: data_result
        dimension :: data_result(size_data_result)

        ! Local data 
        ! -----------------  
        integer :: genPosResult, result_offset, offset, j, nb_tasks
        double precision, allocatable, dimension(:,:) :: B_u, B_v, W_u, W_v
        integer :: iu, iv, ju, jv
        integer :: nnz_col_u, nnz_col_v
        integer, allocatable, dimension(:) :: j_nnz_u, j_nnz_v
        integer :: nnz_row_u, nnz_row_v
        integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v
        integer :: size_data_row
        double precision, allocatable :: data_row(:)

        ! ====================================================
        ! Initialize
        allocate(B_u(nb_rows_u, nb_cols_u), &   
                B_v(nb_rows_v, nb_cols_v), &
                W_u(nb_rows_u, nb_cols_u), &
                W_v(nb_rows_v, nb_cols_v))
        call csr2matrix(size_data_u, indi_u, indj_u, data_B_u, nb_rows_u, nb_cols_u, B_u)
        call csr2matrix(size_data_v, indi_v, indj_v, data_B_v, nb_rows_v, nb_cols_v, B_v)
        call csr2matrix(size_data_u, indi_u, indj_u, data_W_u, nb_rows_u, nb_cols_u, W_u)
        call csr2matrix(size_data_v, indi_v, indj_v, data_W_v, nb_rows_v, nb_cols_v, W_v)
        ! ====================================================

        !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,offset,j_nnz_u,ju,j_nnz_v,jv,nnz_row_u) &
        !$OMP PRIVATE(nnz_row_v,i_nnz_u,i_nnz_v,data_row,size_data_row,genPosResult,result_offset,j)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nb_rows_u*nb_rows_v/nb_tasks)
        do iv = 1, nb_rows_v
            do iu = 1, nb_rows_u
                
                ! FOR COLUMNS
                ! Number of nonzeros  
                nnz_col_u = indi_u(iu+1) - indi_u(iu)
                nnz_col_v = indi_v(iv+1) - indi_v(iv)

                ! Set values
                allocate(j_nnz_u(nnz_col_u))
                offset = indi_u(iu)
                do ju = 1, nnz_col_u
                    j_nnz_u(ju) = indj_u(ju+offset-1)
                end do

                allocate(j_nnz_v(nnz_col_v))
                offset = indi_v(iv)
                do jv = 1, nnz_col_v
                    j_nnz_v(jv) = indj_v(jv+offset-1)
                end do

                ! FOR ROWS
                ! Number of nonzeros 
                nnz_row_u = indi_I_u(iu+1) - indi_I_u(iu)
                nnz_row_v = indi_I_v(iv+1) - indi_I_v(iv)

                ! Set values
                allocate(i_nnz_u(nnz_row_u))
                offset = indi_I_u(iu)
                do ju = 1, nnz_row_u
                    i_nnz_u(ju) = indj_I_u(ju+offset-1)
                end do

                allocate(i_nnz_v(nnz_row_v))
                offset = indi_I_v(iv)
                do jv = 1, nnz_row_v
                    i_nnz_v(jv) = indj_I_v(jv+offset-1)
                end do

                size_data_row = nnz_row_u * nnz_row_v 
                allocate(data_row(size_data_row))

                call csr_get_row_2d(coefs, nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, iu, iv, &
                                nnz_row_u, nnz_row_v, i_nnz_u, i_nnz_v, &
                                nnz_col_u, nnz_col_v, j_nnz_u, j_nnz_v, &
                                B_u, B_v, W_u, W_v, data_row)
                deallocate(i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v)

                ! Get offset in result 
                genPosResult = iu + (iv-1)*nb_rows_u 
                result_offset = indi_result(genPosResult)
            
                ! Get result
                do j = 1, size_data_row
                    data_result(result_offset+j-1) = data_row(j)
                end do   

                deallocate(data_row)
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        deallocate(B_u, B_v, W_u, W_v)

    end subroutine csr_get_matrix_2d 

    subroutine csr_get_row_3d(coefs, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            nb_rows_w, nb_cols_w, &
                            row_u, row_v, row_w, &
                            nnz_row_u, nnz_row_v, nnz_row_w, &
                            i_nnz_u, i_nnz_v, i_nnz_w, &
                            nnz_col_u, nnz_col_v, nnz_col_w, &
                            j_nnz_u, j_nnz_v, j_nnz_w, &
                            B_u, B_v, B_w, W_u, W_v, W_w, &
                            data_row)
        !! Computes a row of a matrix constructed with row_u, row_v and row_w

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v*nb_cols_w)
        integer, intent(in) :: row_u, row_v, row_w
        integer, intent(in) :: nnz_row_u, nnz_row_v, nnz_row_w
        integer, intent(in) :: nnz_col_u, nnz_col_v, nnz_col_w
        integer, intent(in) :: i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w
        dimension :: i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), i_nnz_w(nnz_row_w), &
                    j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v), j_nnz_w(nnz_col_w)
        double precision, intent(in) :: B_u, B_v, B_w, W_u, W_v, W_w
        dimension ::    B_u(nb_rows_u, nb_cols_u), &   
                        B_v(nb_rows_v, nb_cols_v), &
                        B_w(nb_rows_w, nb_cols_w), &
                        W_u(nb_rows_u, nb_cols_u), &
                        W_v(nb_rows_v, nb_cols_v), &
                        W_w(nb_rows_w, nb_cols_w)

        double precision, intent(out) :: data_row
        dimension :: data_row(nnz_row_u*nnz_row_v*nnz_row_w)

        ! Local data 
        ! ----------------- 
        ! Create Ci
        integer :: genPosCoef
        double precision, allocatable :: Ci0(:)

        ! Create Bl and Wl
        integer :: i, j
        double precision, allocatable, dimension(:,:) :: BW_u, BW_v, BW_w, Bt
        double precision, allocatable, dimension(:) :: Wt

        ! Loops
        integer :: genPosC
        integer :: ju, jv, jw, posu, posv, posw

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v*nnz_col_w))

        ! Set values of C
        do jw = 1, nnz_col_w
            do jv = 1, nnz_col_v
                do ju = 1, nnz_col_u
                    posu = j_nnz_u(ju)
                    posv = j_nnz_v(jv)
                    posw = j_nnz_w(jw)
                    genPosCoef = posu + (posv-1)*nb_cols_u + (posw-1)*nb_cols_u*nb_cols_v
                    genPosC = ju + (jv-1)*nnz_col_u + (jw-1)*nnz_col_u*nnz_col_v
                    Ci0(genPosC) = coefs(genPosCoef)                    
                end do
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u,nnz_col_u), Bt(nnz_row_u,nnz_col_u), Wt(nnz_col_u))
        Bt = B_u(i_nnz_u, j_nnz_u)
        Wt = W_u(row_u, j_nnz_u)
        
        forall (i = 1 : size(Bt, 1), j = 1 : size(Bt, 2)) 
            BW_u(i, j) = Bt(i, j) * Wt(j)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
        Bt = B_v(i_nnz_v, j_nnz_v)
        Wt = W_v(row_v, j_nnz_v)
        
        forall (i = 1 : size(Bt, 1), j = 1 : size(Bt, 2)) 
            BW_v(i, j) = Bt(i, j) * Wt(j)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_w(nnz_row_w,nnz_col_w), Bt(nnz_row_w,nnz_col_w), Wt(nnz_col_w))
        Bt = B_w(i_nnz_w, j_nnz_w)
        Wt = W_w(row_w, j_nnz_w)

        forall (i = 1 : size(Bt, 1), j = 1 : size(Bt, 2)) 
            BW_w(i, j) = Bt(i, j) * Wt(j)
        end forall
        deallocate(Bt, Wt)

        call tensor3d_dot_vector(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, &
                                nnz_row_w, nnz_col_w, BW_u, BW_v, BW_w, Ci0, data_row)

    end subroutine csr_get_row_3d 

    subroutine csr_get_matrix_3d(coefs, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, &
                                indi_I_u, indi_I_v, indi_I_w, & 
                                indj_I_u, indj_I_v, indj_I_w, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 3D case (Ww . Bw) x (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v*nb_cols_w)
        integer, intent(in) :: size_data_u, size_data_v, size_data_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                        indi_v(nb_rows_v+1), indj_v(size_data_v), &
                        indi_w(nb_rows_w+1), indj_w(size_data_w)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
        dimension ::    data_B_u(size_data_u), data_W_u(size_data_u), &
                        data_B_v(size_data_v), data_W_v(size_data_v), &
                        data_B_w(size_data_w), data_W_w(size_data_w)
        integer, intent(in) :: size_data_I_u, size_data_I_v, size_data_I_w
        integer, intent(in) :: indi_I_u, indi_I_v, indi_I_w, &
                                indj_I_u, indj_I_v, indj_I_w 
        dimension ::    indi_I_u(nb_rows_u+1), indi_I_v(nb_rows_v+1), indi_I_w(nb_rows_w+1), &
                        indj_I_u(size_data_I_u), indj_I_v(size_data_I_v), indj_I_w(size_data_I_w)
        integer, intent(in) :: indi_result
        dimension :: indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1) 
        integer, intent(in) :: size_data_result

        double precision, intent(out) :: data_result
        dimension :: data_result(size_data_result)

        ! Local data 
        ! -----------------  
        integer :: genPosResult, result_offset, offset, j, nb_tasks
        double precision, allocatable, dimension(:,:) :: B_u, B_v, B_w, W_u, W_v, W_w
        integer :: iu, iv, iw, ju, jv, jw
        integer :: nnz_col_u, nnz_col_v, nnz_col_w
        integer, allocatable, dimension(:) :: j_nnz_u, j_nnz_v, j_nnz_w
        integer :: nnz_row_u, nnz_row_v, nnz_row_w
        integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v, i_nnz_w
        integer :: size_data_row
        double precision, allocatable :: data_row(:)

        ! ====================================================
        ! Initialize
        allocate(B_u(nb_rows_u, nb_cols_u), &   
                B_v(nb_rows_v, nb_cols_v), &
                B_w(nb_rows_w, nb_cols_w), &
                W_u(nb_rows_u, nb_cols_u), &
                W_v(nb_rows_v, nb_cols_v), &
                W_w(nb_rows_w, nb_cols_w))
        call csr2matrix(size_data_u, indi_u, indj_u, data_B_u, nb_rows_u, nb_cols_u, B_u)
        call csr2matrix(size_data_v, indi_v, indj_v, data_B_v, nb_rows_v, nb_cols_v, B_v)
        call csr2matrix(size_data_w, indi_w, indj_w, data_B_w, nb_rows_w, nb_cols_w, B_w)
        call csr2matrix(size_data_u, indi_u, indj_u, data_W_u, nb_rows_u, nb_cols_u, W_u)
        call csr2matrix(size_data_v, indi_v, indj_v, data_W_v, nb_rows_v, nb_cols_v, W_v)
        call csr2matrix(size_data_w, indi_w, indj_w, data_W_w, nb_rows_w, nb_cols_w, W_w)
        ! ====================================================

        !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,nnz_col_w,offset,j_nnz_u,ju,j_nnz_v,jv,j_nnz_w,jw,nnz_row_u) &
        !$OMP PRIVATE(nnz_row_v,nnz_row_w,i_nnz_u,i_nnz_v,i_nnz_w,data_row,size_data_row,genPosResult,result_offset,j)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nb_rows_u*nb_rows_v*nb_rows_w/nb_tasks)
        do iw = 1, nb_rows_w
            do iv = 1, nb_rows_v
                do iu = 1, nb_rows_u
                    
                    ! FOR COLUMNS
                    ! Number of nonzeros  
                    nnz_col_u = indi_u(iu+1) - indi_u(iu)
                    nnz_col_v = indi_v(iv+1) - indi_v(iv)
                    nnz_col_w = indi_w(iw+1) - indi_w(iw)

                    ! Set values
                    allocate(j_nnz_u(nnz_col_u))
                    offset = indi_u(iu)
                    do ju = 1, nnz_col_u
                        j_nnz_u(ju) = indj_u(ju+offset-1)
                    end do

                    allocate(j_nnz_v(nnz_col_v))
                    offset = indi_v(iv)
                    do jv = 1, nnz_col_v
                        j_nnz_v(jv) = indj_v(jv+offset-1)
                    end do

                    allocate(j_nnz_w(nnz_col_w))
                    offset = indi_w(iw)
                    do jw = 1, nnz_col_w
                        j_nnz_w(jw) = indj_w(jw+offset-1)
                    end do

                    ! FOR ROWS
                    ! Number of nonzeros 
                    nnz_row_u = indi_I_u(iu+1) - indi_I_u(iu)
                    nnz_row_v = indi_I_v(iv+1) - indi_I_v(iv)
                    nnz_row_w = indi_I_w(iw+1) - indi_I_w(iw)

                    ! Set values
                    allocate(i_nnz_u(nnz_row_u))
                    offset = indi_I_u(iu)
                    do ju = 1, nnz_row_u
                        i_nnz_u(ju) = indj_I_u(ju+offset-1)
                    end do

                    allocate(i_nnz_v(nnz_row_v))
                    offset = indi_I_v(iv)
                    do jv = 1, nnz_row_v
                        i_nnz_v(jv) = indj_I_v(jv+offset-1)
                    end do

                    allocate(i_nnz_w(nnz_row_w))
                    offset = indi_I_w(iw)
                    do jw = 1, nnz_row_w
                        i_nnz_w(jw) = indj_I_w(jw+offset-1)
                    end do

                    size_data_row = nnz_row_u * nnz_row_v * nnz_row_w
                    allocate(data_row(size_data_row))

                    call csr_get_row_3d(coefs, nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, iu, iv, iw, &
                                    nnz_row_u, nnz_row_v, nnz_row_w, i_nnz_u, i_nnz_v, i_nnz_w, &
                                    nnz_col_u, nnz_col_v, nnz_col_w, j_nnz_u, j_nnz_v, j_nnz_w, &
                                    B_u, B_v, B_w, W_u, W_v, W_w, data_row)
                    deallocate(i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w)

                    ! Get offset in result 
                    genPosResult = iu + (iv-1)*nb_rows_u + (iw-1)*nb_rows_u*nb_rows_v
                    result_offset = indi_result(genPosResult)
                
                    ! Get result
                    do j = 1, size_data_row
                        data_result(result_offset+j-1) = data_row(j)
                    end do   

                    deallocate(data_row)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        deallocate(B_u, B_v, B_w, W_u, W_v, W_w)

    end subroutine csr_get_matrix_3d 

    ! ----------------------------------------------------
    ! Functions for Fast Diagonalization method
    ! ----------------------------------------------------
    ! "Fast Diagonalization" 

    subroutine eigen_decomposition(nb_rows, nb_cols, &
                                    Mcoef, Kcoef, size_data, indi, indj, &
                                    data_B0, data_W00, data_B1, data_W11, &
                                    Method, D, U, Kdiag, Mdiag)
        !! Eigen decomposition generalized KU = MUD
        !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
        !! M: mass matrix, M = int B0 B0 dx = W00 * B0
        !! U: eigenvectors matrix
        !! D: diagonal of eigenvalues
        !! IN CSR FORMAT
        
        implicit none 
        ! Input / output 
        ! -------------------
        integer, intent(in) :: nb_rows, nb_cols
        double precision, dimension(*), intent(in) :: Mcoef, Kcoef
        integer, intent(in) :: size_data
        integer, intent(in) :: indi, indj
        dimension :: indi(nb_rows+1), indj(size_data)
        double precision, intent(in) :: data_B0, data_W00, data_B1, data_W11
        dimension ::    data_B0(size_data), data_W00(size_data), &
                        data_B1(size_data), data_W11(size_data)
        character(len=10), intent(in) :: Method
                
        double precision, intent(out) :: D, U
        dimension :: D(nb_rows), U(nb_rows, nb_rows)
        double precision, intent(out) :: Kdiag, Mdiag
        dimension :: Kdiag(nb_rows), Mdiag(nb_rows)

        ! Local data
        ! ----------------
        double precision, allocatable, dimension(:) :: data_B0_new, data_B1_new
        double precision, allocatable, dimension(:,:) :: BB0, WW0, MM, BB1, WW1, KK
        integer :: i, j

        ! Use Lapack
        integer, parameter :: itype = 1 ! Type A x = B x D where in this proble A =  KK and B = MM
        character, parameter :: jobz = 'V' ! Computes eigen values and vectors
        character, parameter :: uplo = 'U' ! Consider upper triangle of A and B
        integer :: N, lda, ldb ! Size of A and B

        integer :: lwork, liwork
        double precision, allocatable, dimension(:) :: W, work, iwork
        double precision :: dummy(1)
        integer :: idum(1)
        integer :: info
        
        ! Initialize Masse matrix
        allocate(data_B0_new(size_data))
        data_B0_new = data_B0
        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            do i = 1, nb_rows
                do j = indi(i), indi(i+1)-1
                    data_B0_new(j) = data_B0_new(j)*Mcoef(indj(j))
                end do
            end do
        end if

        allocate(BB0(nb_rows, nb_cols))
        call csr2matrix(size_data, indi, indj, data_B0_new, nb_rows, nb_cols, BB0)
        deallocate(data_B0_new)
        allocate(WW0(nb_rows, nb_cols))
        call csr2matrix(size_data, indi, indj, data_W00, nb_rows, nb_cols, WW0)
        allocate(MM(nb_rows, nb_rows))
        MM = matmul(WW0, transpose(BB0))
        deallocate(BB0, WW0)

        ! Initialize Stiffness matrix
        allocate(data_B1_new(size_data))
        data_B1_new = data_B1
        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then
            do i = 1, nb_rows
                do j = indi(i), indi(i+1)-1
                    data_B1_new(j) = data_B1_new(j)*Kcoef(indj(j))
                end do
            end do
        end if

        allocate(BB1(nb_rows, nb_cols))
        call csr2matrix(size_data, indi, indj, data_B1_new, nb_rows, nb_cols, BB1)
        deallocate(data_B1_new)
        allocate(WW1(nb_rows, nb_cols))
        call csr2matrix(size_data, indi, indj, data_W11, nb_rows, nb_cols, WW1)
        allocate(KK(nb_rows, nb_rows))
        KK = matmul(WW1, transpose(BB1))
        deallocate(BB1, WW1)

        ! Select diagonal of M and K
        do j = 1, nb_rows
            Kdiag(j) = KK(j, j)
            Mdiag(j) = MM(j, j)
        end do

        ! ====================================================
        ! Eigen decomposition KK U = MM U DD
        N = nb_rows
        lda = nb_rows
        ldb = nb_rows

        ! Set eigen vectors and eigen values
        allocate(W(N))

        ! Use routine workspace query to get optimal workspace.
        lwork = -1
        liwork = -1
        call dsygvd(itype, jobz, uplo, N, KK, lda, MM, ldb, W, dummy, lwork, idum, liwork, info)

        ! Make sure that there is enough workspace 
        lwork = max(1+(6+2*n)*n, nint(dummy(1)))
        liwork = max(3+5*n, idum(1))
        allocate (work(lwork), iwork(liwork))

        ! Solve
        call dsygvd(itype, jobz, uplo, N, KK, lda, MM, ldb, W, work, lwork, iwork, liwork, info)

        ! Get values
        U = KK
        D = W
        deallocate(KK, MM, W, work, iwork)

    end subroutine eigen_decomposition

    subroutine fast_diagonalization_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                                        U_u, U_v, U_w, diagonal, array_in, array_out)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w
        double precision, intent(in) :: U_u, U_v, U_w, diagonal
        dimension ::    U_u(nb_rows_u, nb_rows_u), &
                        U_v(nb_rows_v, nb_rows_v), &
                        U_w(nb_rows_w, nb_rows_w), &
                        diagonal(nb_rows_total)

        double precision, intent(in) :: array_in
        dimension :: array_in(nb_rows_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nb_rows_total)

        ! Local data
        ! -------------
        integer :: i, nb_tasks
        double precision :: array_temp_1
        dimension :: array_temp_1(nb_rows_total)

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nb_rows_u, nb_rows_u, nb_rows_v, nb_rows_v, nb_rows_w, nb_rows_w, &
                            transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp_1)

        ! ---------------------------------
        ! Second part 
        ! ---------------------------------
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nb_rows_total/nb_tasks)
        do i = 1, nb_rows_total
            array_temp_1(i) = array_temp_1(i)/diagonal(i)
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        ! ----------------------------------
        ! Third part
        ! ---------------------------------
        call tensor3d_dot_vector(nb_rows_u, nb_rows_u, nb_rows_v, nb_rows_v, nb_rows_w, nb_rows_w, &
                                U_u, U_v, U_w, array_temp_1, array_out)
        
    end subroutine fast_diagonalization_3d

    subroutine precond_interp_3d(nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w, &
                                        U_u, U_v, U_w, array_in, array_out)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nb_rows_total, nb_rows_u, nb_rows_v, nb_rows_w
        double precision, intent(in) :: U_u, U_v, U_w
        dimension ::    U_u(nb_rows_u, nb_rows_u), &
                        U_v(nb_rows_v, nb_rows_v), &
                        U_w(nb_rows_w, nb_rows_w)

        double precision, intent(in) :: array_in
        dimension :: array_in(nb_rows_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nb_rows_total)

        ! Local data
        ! -------------
        double precision :: array_temp_1
        dimension :: array_temp_1(nb_rows_total)

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nb_rows_u, nb_rows_u, nb_rows_v, nb_rows_v, nb_rows_w, nb_rows_w, &
                            transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp_1)

        ! ----------------------------------
        ! Second part
        ! ---------------------------------
        call tensor3d_dot_vector(nb_rows_u, nb_rows_u, nb_rows_v, nb_rows_v, nb_rows_w, nb_rows_w, &
                                U_u, U_v, U_w, array_temp_1, array_out)

    end subroutine precond_interp_3d

    ! For improving fast diagonalisation (TD, TDS, JM and JMS)

    subroutine tensor_decomposition_2d(nb_cols_total, nb_cols_u, nb_cols_v, CC, &
                                        M_u, M_v, K_u, K_v)
        !! Tensor decomposition to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini
        
        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_cols_total, nb_cols_u, nb_cols_v
        double precision, intent(in) :: CC
        dimension :: CC(2, 2, nb_cols_total)

        double precision, intent(inout) :: M_u, M_v, K_u, K_v
        dimension ::    M_u(nb_cols_u), M_v(nb_cols_v), &
                        K_u(nb_cols_u), K_v(nb_cols_v)

        ! Local data
        ! ---------------
        double precision :: Vscript(nb_cols_u, nb_cols_v), & 
                            Wscript(nb_cols_u, nb_cols_v)
        ! Loop
        integer :: k, l, i1, i2
        integer :: genpos
        double precision :: vmin, vmax
        double precision :: UU(2), WW(2)

        do k = 1, 2
            ! Set Vscript
            Vscript = 0.d0
            do i2 = 1, nb_cols_v
                do i1 = 1, nb_cols_u
                    genpos = i1 + (i2-1)*nb_cols_u 
                    UU = [M_u(i1), M_v(i2)]
                    Vscript(i1, i2) = CC(k, k, genpos)*UU(k)/(UU(1)*UU(2))
                end do
            end do

            ! Update w
            if (k.eq.1) then 
                do i1 = 1, nb_cols_u
                    vmin = minval(Vscript(i1, :))
                    vmax = maxval(Vscript(i1, :))
                    K_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_cols_v
                    vmin = minval(Vscript(:, i2))
                    vmax = maxval(Vscript(:, i2))
                    K_v(i2) = sqrt(vmin*vmax)
                end do
            end if
        end do

        do k = 1, 2
            ! Initialize
            Wscript = 0.d0
            do l = 1, 2
                if (k.ne.l) then 
                    ! Set Wscript
                    do i2 = 1, nb_cols_v
                        do i1 = 1, nb_cols_u
                            genpos = i1 + (i2-1)*nb_cols_u 
                            UU = [M_u(i1), M_v(i2)]
                            WW = [K_u(i1), K_v(i2)]
                            Wscript(i1, i2) = CC(k, k, genpos)*UU(k)*UU(l)&
                                                        /(UU(1)*UU(2)*WW(k))
                        end do
                    end do
                end if
            end do

            ! Update u
            if (k.eq.1) then 
                do i1 = 1, nb_cols_u
                    vmin = minval(Wscript(i1, :))
                    vmax = maxval(Wscript(i1, :))
                    M_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_cols_v
                    vmin = minval(Wscript(:, i2))
                    vmax = maxval(Wscript(:, i2))
                    M_v(i2) = sqrt(vmin*vmax)
                end do
            end if
        end do
                                    
    end subroutine tensor_decomposition_2d

    subroutine tensor_decomposition_3d(nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w, CC, &
                                        M_u, M_v, M_w, K_u, K_v, K_w)
        !! Tensor decomposition to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_cols_total, nb_cols_u, nb_cols_v, nb_cols_w
        double precision, intent(in) :: CC
        dimension :: CC(3, 3, nb_cols_total)

        double precision, intent(inout) :: M_u, M_v, M_w, K_u, K_v, K_w
        dimension ::    M_u(nb_cols_u), M_v(nb_cols_v), M_w(nb_cols_w), &
                        K_u(nb_cols_u), K_v(nb_cols_v), K_w(nb_cols_w)

        ! Local data
        ! ---------------
        double precision :: Vscript(nb_cols_u, nb_cols_v, nb_cols_w), & 
                            Wscript(2, nb_cols_u, nb_cols_v, nb_cols_w), &
                            Mscript(nb_cols_u, nb_cols_v, nb_cols_w), &
                            Nscript(nb_cols_u, nb_cols_v, nb_cols_w)

        ! Loop
        integer :: k, l, i1, i2, i3
        integer :: genpos
        double precision :: vmin, vmax
        integer :: cont
        double precision :: UU(3), WW(3), WWlk(2)

        do k = 1, 3
            ! Set Vscript
            Vscript = 0.d0
            do i3 = 1, nb_cols_w
                do i2 = 1, nb_cols_v
                    do i1 = 1, nb_cols_u
                        genpos = i1 + (i2-1)*nb_cols_u + (i3-1)*nb_cols_u*nb_cols_v
                        UU = [M_u(i1), M_v(i2), M_w(i3)]
                        Vscript(i1, i2, i3) = CC(k, k, genpos)*UU(k)/(UU(1)*UU(2)*UU(3))
                    end do
                end do
            end do

            ! Update w
            if (k.eq.1) then 
                do i1 = 1, nb_cols_u
                    vmin = minval(Vscript(i1, :, :))
                    vmax = maxval(Vscript(i1, :, :))
                    K_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_cols_v
                    vmin = minval(Vscript(:, i2, :))
                    vmax = maxval(Vscript(:, i2, :))
                    K_v(i2) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do i3 = 1, nb_cols_v
                    vmin = minval(Vscript(:, :, i3))
                    vmax = maxval(Vscript(:, :, i3))
                    K_w(i3) = sqrt(vmin*vmax)
                end do
            end if
        end do

        do k = 1, 3
            ! Initialize
            Wscript = 0.d0
            Nscript = 0.d0
            Mscript = 0.d0
            cont = 0
            do l = 1, 3
                if (k.ne.l) then 
                    cont = cont + 1
                    ! Set Wscript
                    do i3 = 1, nb_cols_w
                        do i2 = 1, nb_cols_v
                            do i1 = 1, nb_cols_u
                                genpos = i1 + (i2-1)*nb_cols_u + (i3-1)*nb_cols_u*nb_cols_v
                                UU = [M_u(i1), M_v(i2), M_w(i3)]
                                WW = [K_u(i1), K_v(i2), K_w(i3)]
                                Wscript(cont, i1, i2, i3) = CC(k, k, genpos)*UU(k)*UU(l)&
                                                            /(UU(1)*UU(2)*UU(3)*WW(k))
                            end do
                        end do
                    end do
                end if
            end do

            ! Compute Nscript and Mscript
            do i3 = 1, nb_cols_w
                do i2 = 1, nb_cols_v
                    do i1 = 1, nb_cols_u
                        WWlk = Wscript(:, i1, i2, i3)
                        Nscript(i1, i2, i3) = minval(WWlk)
                        Mscript(i1, i2, i3) = maxval(WWlk)
                    end do
                end do
            end do

            ! Update u
            if (k.eq.1) then 
                do i1 = 1, nb_cols_u
                    vmin = minval(Nscript(i1, :, :))
                    vmax = maxval(Mscript(i1, :, :))
                    M_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_cols_v
                    vmin = minval(Nscript(:, i2, :))
                    vmax = maxval(Mscript(:, i2, :))
                    M_v(i2) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do i3 = 1, nb_cols_v
                    vmin = minval(Nscript(:, :, i3))
                    vmax = maxval(Mscript(:, :, i3))
                    M_w(i3) = sqrt(vmin*vmax)
                end do
            end if
        end do

    end subroutine tensor_decomposition_3d

    subroutine jacobien_mean_3d(nb_cols_u, nb_cols_v, nb_cols_w, size_J, JJ, size_K, KK, &
                                L1, L2, L3, lambda1, lambda2, lambda3)
        
        use omp_lib
        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_cols_u, nb_cols_v, nb_cols_w, size_J, size_K
        double precision, intent(in) :: JJ, KK
        dimension :: JJ(3, 3, size_J), KK(3, 3, size_K)
    
        double precision, intent(out) :: L1, L2, L3
        double precision, intent(out) :: lambda1, lambda2, lambda3
        integer :: step
    
        ! Local data
        ! --------------
        integer :: i, j, k, nb_pts, nb_pts_temp, nb_tasks, pos
        double precision, dimension(3,3) :: dist, lambda, MatrixT, Q
        double precision :: sq

        ! Define Step
        step = min((nb_cols_u-1)/2, (nb_cols_v-1)/2, (nb_cols_w-1)/2)
        
        ! Count number of quadrature points
        nb_pts = 1
        nb_pts_temp = 0
        do i = 1, nb_cols_u, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
        nb_pts_temp = 0

        do j = 1, nb_cols_v, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
        nb_pts_temp = 0

        do k = 1, nb_cols_w, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
    
        ! Compute distance
        !--------------------------
        L1 = 0.d0
        L2 = 0.d0
        L3 = 0.d0
    
        !$OMP PARALLEL PRIVATE(pos, MatrixT, dist) REDUCTION(+:L1, L2, L3)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size_J/(nb_tasks*step*step*step))
        do k = 1, nb_cols_w, step
            do j = 1, nb_cols_v, step
                do i = 1, nb_cols_u, step
                    pos = i + (j-1)*nb_cols_u + (k-1)*nb_cols_u*nb_cols_v
                    MatrixT = JJ(:, :, pos)
                    call polar_decomposition(MatrixT, Q, dist, 1, 1)
    
                    ! Find mean of diagonal of jacobien
                    L1 = L1 + dist(1, 1)/nb_pts
                    L2 = L2 + dist(2, 2)/nb_pts
                    L3 = L3 + dist(3, 3)/nb_pts
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

        ! Dimension normalized
        sq = sqrt(L1**2 + L2**2 + L3**2)
        L1 = L1/sq
        L2 = L2/sq
        L3 = L3/sq

        ! Compute conductivity
        !--------------------------
        if (size_K.eq.1) then
            call polar_decomposition(KK(:, :, 1), Q, lambda, 1, 1)

            lambda1 = lambda(1, 1)
            lambda2 = lambda(2, 2)
            lambda3 = lambda(3, 3)

        else if (size_K.eq.size_J) then

            lambda1 = 0.d0
            lambda2 = 0.d0
            lambda3 = 0.d0

            !$OMP PARALLEL PRIVATE(pos, MatrixT, Q, lambda) REDUCTION(+:lambda1, lambda2, lambda3)
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, size_J/(nb_tasks*step*step*step))
            do k = 1, nb_cols_w, step
                do j = 1, nb_cols_v, step
                    do i = 1, nb_cols_u, step
                        pos = i + (j-1)*nb_cols_u + (k-1)*nb_cols_u*nb_cols_v
                        MatrixT = KK(:, :, pos)
                        call polar_decomposition(MatrixT, Q, lambda, 1, 1)
        
                        ! Find mean of diagonal of jacobien
                        lambda1 = lambda1 + lambda(1, 1)/nb_pts
                        lambda2 = lambda2 + lambda(2, 2)/nb_pts
                        lambda3 = lambda3 + lambda(3, 3)/nb_pts
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 

            ! Conductivity normalized
            sq = sqrt(lambda1**2 + lambda2**2 + lambda3**2)
            lambda1 = lambda1/sq
            lambda2 = lambda2/sq
            lambda3 = lambda3/sq

        end if
        
    end subroutine jacobien_mean_3d

    ! For scaling (TDS and JMS)
    
    subroutine find_parametric_diag_3d(nb_rows_u, nb_rows_v, nb_rows_w, c_u, c_v, c_w, &
                                Mdiag_u, Mdiag_v, Mdiag_w, &
                                Kdiag_u, Kdiag_v, Kdiag_w, diag)
        !! Find the diagonal of the preconditioner "fast diagonalization"
                            
        implicit none
        ! Input / output data
        ! -------------------------
        integer, intent(in) :: nb_rows_u, nb_rows_v, nb_rows_w
        double precision, intent(in) :: c_u, c_v, c_w
        double precision, intent(in) :: Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w
        dimension :: Mdiag_u(nb_rows_u), Mdiag_v(nb_rows_v), Mdiag_w(nb_rows_w), &
                    Kdiag_u(nb_rows_u), Kdiag_v(nb_rows_v), Kdiag_w(nb_rows_w)

        double precision, intent(out) :: diag
        dimension :: diag(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Initialize
        diag = 0.d0

        ! Find K3 M2 M1
        call kron_product_3vec(nb_rows_w, Kdiag_w, nb_rows_v, Mdiag_v, nb_rows_u, Mdiag_u, diag, c_w)

        ! Find M3 K2 M1
        call kron_product_3vec(nb_rows_w, Mdiag_w, nb_rows_v, Kdiag_v, nb_rows_u, Mdiag_u, diag, c_v)

        ! Find M3 M2 K1
        call kron_product_3vec(nb_rows_w, Mdiag_w, nb_rows_v, Mdiag_v, nb_rows_u, Kdiag_u, diag, c_u)

    end subroutine find_parametric_diag_3d

    subroutine find_physical_diag_3d(coefs, nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w, &
                                        size_data_u, size_data_v, size_data_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, &
                                        data_W_u, data_W_v, data_W_w, diag)
        !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
        !! Algotihm based on sum factorization adapted to diagonal case 
        !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! -------------------
        integer, intent(in) :: nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v*nb_cols_w)
        integer, intent(in) :: size_data_u, size_data_v, size_data_w
        integer, intent(in) :: indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                        indi_v(nb_rows_v+1), indj_v(size_data_v), &
                        indi_w(nb_rows_w+1), indj_w(size_data_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w
        dimension ::    data_B_u(size_data_u), data_B_v(size_data_v), data_B_w(size_data_w), &
                        data_W_u(size_data_u), data_W_v(size_data_v), data_W_w(size_data_w)

        double precision, intent(inout) :: diag
        dimension :: diag(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data
        ! ------------------
        integer :: offset, nb_tasks
        integer :: iu, iv, iw, ju, jv, jw, Cpos, Ipos
        double precision :: sum1, sum2, sum3

        integer :: nnz_u, nnz_v, nnz_w
        integer, allocatable, dimension(:) :: indj_nnz_u, indj_nnz_v, indj_nnz_w
        double precision, allocatable, dimension(:) :: data_nnz_B_u, data_nnz_B_v, data_nnz_B_w, &
                                                        data_nnz_W_u, data_nnz_W_v, data_nnz_W_w        

        !$OMP PARALLEL PRIVATE(ju,jv,jw,nnz_u,nnz_v,nnz_w,offset,indj_nnz_u,data_nnz_B_u,data_nnz_W_u) &
        !$OMP PRIVATE(indj_nnz_v,data_nnz_B_v,data_nnz_W_v,indj_nnz_w,data_nnz_B_w,data_nnz_W_w,sum1,sum2,sum3,Cpos,Ipos)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nb_rows_u*nb_rows_v*nb_rows_w /nb_tasks) 
        do iw = 1, nb_rows_w
            do iv = 1, nb_rows_v
                do iu = 1, nb_rows_u

                    ! Number of nonzeros
                    nnz_u = indi_u(iu+1) - indi_u(iu)
                    nnz_v = indi_v(iv+1) - indi_v(iv)
                    nnz_w = indi_w(iw+1) - indi_w(iw)

                    ! Set values
                    allocate(indj_nnz_u(nnz_u), data_nnz_B_u(nnz_u), data_nnz_W_u(nnz_u))
                    offset = indi_u(iu)
                    do ju = 1, nnz_u
                        indj_nnz_u(ju) = indj_u(ju+offset-1)
                        data_nnz_B_u(ju) = data_B_u(ju+offset-1)
                        data_nnz_W_u(ju) = data_W_u(ju+offset-1)
                    end do

                    allocate(indj_nnz_v(nnz_v), data_nnz_B_v(nnz_v), data_nnz_W_v(nnz_v))
                    offset = indi_v(iv)
                    do jv = 1, nnz_v
                        indj_nnz_v(jv) = indj_v(jv+offset-1)
                        data_nnz_B_v(jv) = data_B_v(jv+offset-1)
                        data_nnz_W_v(jv) = data_W_v(jv+offset-1)
                    end do

                    allocate(indj_nnz_w(nnz_w), data_nnz_B_w(nnz_w), data_nnz_W_w(nnz_w))
                    offset = indi_w(iw)
                    do jw = 1, nnz_w
                        indj_nnz_w(jw) = indj_w(jw+offset-1)
                        data_nnz_B_w(jw) = data_B_w(jw+offset-1)
                        data_nnz_W_w(jw) = data_W_w(jw+offset-1)
                    end do

                    sum3 = 0.d0
                    do jw = 1, nnz_w
                        sum2 = 0.d0
                        do jv = 1, nnz_v
                            sum1 = 0.d0
                            do ju = 1, nnz_u
                                Cpos = indj_nnz_u(ju) + (indj_nnz_v(jv)-1)*nb_cols_u + (indj_nnz_w(jw)-1)*nb_cols_u*nb_cols_v
                                sum1 = sum1 + data_nnz_W_u(ju)*data_nnz_B_u(ju)*coefs(Cpos)
                            end do
                            sum2 = sum2 + data_nnz_W_v(jv)*data_nnz_B_v(jv)*sum1
                        end do
                        sum3 = sum3 + data_nnz_W_w(jw)*data_nnz_B_w(jw)*sum2
                    end do 

                    ! General position
                    Ipos = iu + (iv-1)*nb_rows_u + (iw-1)*nb_rows_u*nb_rows_w
                    
                    ! Update diagonal
                    diag(Ipos) = diag(Ipos) + sum3

                    deallocate(indj_nnz_u, data_nnz_B_u, data_nnz_W_u)
                    deallocate(indj_nnz_v, data_nnz_B_v, data_nnz_W_v)
                    deallocate(indj_nnz_w, data_nnz_B_w, data_nnz_W_w)

                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine find_physical_diag_3d

    subroutine scaling_FastDiag(size_vector, diag_par, diag_phy, vector)
        !! Scaling in fast diagonalization
    
        use omp_lib
        implicit none
        ! Input / output data
        ! -------------------
        integer, intent(in) :: size_vector
        double precision, intent(in) :: diag_par, diag_phy
        dimension :: diag_par(size_vector), diag_phy(size_vector)
    
        double precision, intent(inout) :: vector
        dimension :: vector(size_vector)
    
        ! Local data
        ! -------------
        integer :: i, nb_tasks
    
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, size_vector/nb_tasks)
        do i = 1, size_vector
            vector(i) = sqrt(diag_par(i)/diag_phy(i)) * vector(i) 
        end do  
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 
    
    end subroutine scaling_FastDiag

end module tensor_methods