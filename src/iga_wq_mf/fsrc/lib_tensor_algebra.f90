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
        double precision, allocatable, dimension(:,:) :: Xt, Rt
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

    subroutine tensor2d_dot_vector(nr_u, nc_u, nr_v, nc_v, &
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
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nc_u*nc_v)
        double precision, intent(in) :: Mu, Mv
        dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nr_u*nr_v)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1
    
        ! First product
        allocate(R1(nr_u*nc_v))
        call tensor_n_mode_product(nc_u, nc_v, 1, vector_in, nr_u, nc_u, Mu, 1, nr_u, nc_v, 1, R1)

        ! Second product
        call tensor_n_mode_product(nr_u, nc_v, 1, R1, nr_v, nc_v, Mv, 2, nr_u, nr_v, 1, vector_out)
        deallocate(R1)

    end subroutine tensor2d_dot_vector

    subroutine tensor3d_dot_vector(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
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
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v,nr_w, nc_w
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nc_u*nc_v*nc_w)
        double precision, intent(in) :: Mu, Mv, Mw
        dimension :: Mu(nr_u, nc_u), Mv(nr_v, nc_v), Mw(nr_w, nc_w)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nr_u*nr_v*nr_w)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1, R2

        ! First product
        allocate(R1(nr_u*nc_v*nc_w))
        call tensor_n_mode_product(nc_u, nc_v, nc_w, vector_in, nr_u, nc_u, Mu, 1, nr_u, nc_v, nc_w, R1)

        ! Second product
        allocate(R2(nr_u*nr_v*nc_w))
        call tensor_n_mode_product(nr_u, nc_v, nc_w, R1, nr_v, nc_v, Mv, 2, nr_u, nr_v, nc_w, R2)
        deallocate(R1)

        ! Third product
        call tensor_n_mode_product(nr_u, nr_v, nc_w, R2, nr_w, nc_w, Mw, 3, nr_u, nr_v, nr_w, vector_out)
        deallocate(R2)

    end subroutine tensor3d_dot_vector

    subroutine tensor2d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, &
                                    nnz_u, indi_u, indj_u, data_u, &
                                    nnz_v, indi_v, indj_v, data_v, &
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
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nc_u*nc_v)
        double precision, intent(in) :: data_u, data_v
        dimension :: data_u(nnz_u), data_v(nnz_v)
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nr_u*nr_v)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1

        ! First product
        allocate(R1(nr_u*nc_v))
        call tensor_n_mode_product_sp(nc_u, nc_v, 1, vector_in, nr_u, nc_u, nnz_u, &
                                    data_u, indi_u, indj_u, 1, nr_u, nc_v, 1, R1)

        ! Second product
        call tensor_n_mode_product_sp(nr_u, nc_v, 1, R1, nr_v, nc_v, nnz_v, &
                                    data_v, indi_v, indj_v, 2, nr_u, nr_v, 1, vector_out)
        deallocate(R1)

    end subroutine tensor2d_dot_vector_sp

    subroutine tensor3d_dot_vector_sp(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                    nnz_u, indi_u, indj_u, data_u, &
                                    nnz_v, indi_v, indj_v, data_v, &
                                    nnz_w, indi_w, indj_w, data_w, &
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
        integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v,nr_w, nc_w, nnz_u, nnz_v, nnz_w
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nc_u*nc_v*nc_w)
        double precision, intent(in) :: data_u, data_v, data_w
        dimension :: data_u(nnz_u), data_v(nnz_v), data_w(nnz_w)
        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nr_u*nr_v*nr_w)

        ! Local data 
        ! -------------
        double precision, allocatable, dimension(:) :: R1, R2

        ! First product
        allocate(R1(nr_u*nc_v*nc_w))
        call tensor_n_mode_product_sp(nc_u, nc_v, nc_w, vector_in, nr_u, nc_u, nnz_u, &
                                    data_u, indi_u, indj_u, 1, nr_u, nc_v, nc_w, R1)

        ! Second product
        allocate(R2(nr_u*nr_v*nc_w))
        call tensor_n_mode_product_sp(nr_u, nc_v, nc_w, R1, nr_v, nc_v, nnz_v, &
                                    data_v, indi_v, indj_v, 2, nr_u, nr_v, nc_w, R2)
        deallocate(R1)

        ! Third product
        call tensor_n_mode_product_sp(nr_u, nr_v, nc_w, R2, nr_w, nc_w, nnz_w, &
                                    data_w, indi_w, indj_w, 3, nr_u, nr_v, nr_w, vector_out)
        deallocate(R2)

    end subroutine tensor3d_dot_vector_sp

    ! ----------------------------------------------------
    ! Sum product to compute matrices 
    ! ----------------------------------------------------

    subroutine csr_get_row_2d(coefs, nr_u, nc_u, nr_v, nc_v, row_u, row_v, &
                            nnz_row_u, nnz_row_v, i_nnz_u, i_nnz_v, &
                            nnz_col_u, nnz_col_v, j_nnz_u, j_nnz_v, &
                            B_u, B_v, W_u, W_v, data_row)
        !! Computes a row of a matrix constructed with row_u and row_v

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v, row_u, row_v
        double precision, intent(in) :: coefs
        dimension :: coefs(nc_u*nc_v)
        integer, intent(in) :: nnz_row_u, nnz_row_v, nnz_col_u, nnz_col_v
        integer, intent(in) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
        dimension ::    i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), &
                        j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v)
        double precision, intent(in) :: B_u, B_v, W_u, W_v
        dimension ::    B_u(nr_u, nc_u), W_u(nr_u, nc_u), &   
                        B_v(nr_v, nc_v), W_v(nr_v, nc_v)                         

        double precision, intent(out) :: data_row
        dimension :: data_row(nnz_row_u*nnz_row_v)

        ! Local data 
        ! ----------------- 
        integer :: iu, iv, ju, jv, posu, posv, posCoef, genPosC
        double precision, allocatable, dimension(:) :: Ci0, Wt
        double precision, allocatable, dimension(:,:) :: BW_u, BW_v, Bt

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v))

        ! Set values of C
        do jv = 1, nnz_col_v
            do ju = 1, nnz_col_u
                posu = j_nnz_u(ju)
                posv = j_nnz_v(jv)
                posCoef = posu + (posv-1)*nc_u 
                genPosC = ju + (jv-1)*nnz_col_u 
                Ci0(genPosC) = coefs(posCoef)                    
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u, nnz_col_u), Bt(nnz_row_u, nnz_col_u), Wt(nnz_col_u))
        Bt = B_u(i_nnz_u, j_nnz_u)
        Wt = W_u(row_u, j_nnz_u)
        
        forall (iu = 1 : size(Bt, 1), ju = 1 : size(Bt, 2)) 
            BW_u(iu, ju) = Bt(iu, ju) * Wt(ju)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
        Bt = B_v(i_nnz_v, j_nnz_v)
        Wt = W_v(row_v, j_nnz_v)

        forall (iv = 1 : size(Bt, 1), jv = 1 : size(Bt, 2)) 
            BW_v(iv, jv) = Bt(iv, jv) * Wt(jv)
        end forall
        deallocate(Bt, Wt)

        call tensor2d_dot_vector(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, BW_u, BW_v, Ci0, data_row)

    end subroutine csr_get_row_2d 

    subroutine csr_get_matrix_2d(coefs, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v, &
                                nnz_I_u, nnz_I_v, indi_I_u, indi_I_v, indj_I_u, indj_I_v, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 2D case (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, nnz_I_u, nnz_I_v
        double precision, intent(in) :: coefs
        dimension :: coefs(nc_u*nc_v)
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
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
        integer, intent(in) :: size_data_result

        double precision, intent(out) :: data_result
        dimension :: data_result(size_data_result)

        ! Local data 
        ! -----------------  
        integer :: iu, iv, ju, jv, k, nb_tasks, genPos, offset
        double precision, allocatable, dimension(:,:) :: B_u, B_v, W_u, W_v
        
        integer :: nnz_row_u, nnz_row_v, nnz_col_u, nnz_col_v, nnz_gen_row
        integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
        double precision, allocatable :: data_row(:)

        ! ====================================================
        ! Initialize
        allocate(B_u(nr_u, nc_u), B_v(nr_v, nc_v), W_u(nr_u, nc_u), W_v(nr_v, nc_v))
        call csr2dense(nnz_u, indi_u, indj_u, data_B_u, nr_u, nc_u, B_u)
        call csr2dense(nnz_v, indi_v, indj_v, data_B_v, nr_v, nc_v, B_v)
        call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W_u)
        call csr2dense(nnz_v, indi_v, indj_v, data_W_v, nr_v, nc_v, W_v)
        ! ====================================================

        !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,offset,j_nnz_u,ju,j_nnz_v,jv,nnz_row_u) &
        !$OMP PRIVATE(nnz_row_v,i_nnz_u,i_nnz_v,data_row,nnz_gen_row,genPos,k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nr_u*nr_v/nb_tasks)
        do iv = 1, nr_v
            do iu = 1, nr_u

                ! FOR ROWS
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
                
                ! FOR COLUMNS
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

                nnz_gen_row = nnz_row_u * nnz_row_v 
                allocate(data_row(nnz_gen_row))

                call csr_get_row_2d(coefs, nr_u, nc_u, nr_v, nc_v, iu, iv, &
                                nnz_row_u, nnz_row_v, i_nnz_u, i_nnz_v, &
                                nnz_col_u, nnz_col_v, j_nnz_u, j_nnz_v, &
                                B_u, B_v, W_u, W_v, data_row)
                deallocate(i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v)

                ! Get offset in result 
                genPos = iu + (iv-1)*nr_u 
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
                            B_u, B_v, B_w, W_u, W_v, W_w, data_row)
        !! Computes a row of a matrix constructed with row_u, row_v and row_w

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, row_u, row_v, row_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nc_u*nc_v*nc_w)
        integer, intent(in) :: nnz_row_u, nnz_row_v, nnz_row_w, nnz_col_u, nnz_col_v, nnz_col_w
        integer, intent(in) :: i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w
        dimension ::    i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), i_nnz_w(nnz_row_w), &
                        j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v), j_nnz_w(nnz_col_w)
        double precision, intent(in) :: B_u, B_v, B_w, W_u, W_v, W_w
        dimension ::    B_u(nr_u, nc_u), W_u(nr_u, nc_u), &   
                        B_v(nr_v, nc_v), W_v(nr_v, nc_v), &
                        B_w(nr_w, nc_w), W_w(nr_w, nc_w)

        double precision, intent(out) :: data_row
        dimension :: data_row(nnz_row_u*nnz_row_v*nnz_row_w)

        ! Local data 
        ! ----------------- 
        integer :: iu, iv, iw, ju, jv, jw, posu, posv, posw, posCoef, genPosC
        double precision, allocatable, dimension(:) :: Ci0, Wt
        double precision, allocatable, dimension(:,:) :: BW_u, BW_v, BW_w, Bt

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v*nnz_col_w))

        ! Set values of C
        do jw = 1, nnz_col_w
            do jv = 1, nnz_col_v
                do ju = 1, nnz_col_u
                    posu = j_nnz_u(ju)
                    posv = j_nnz_v(jv)
                    posw = j_nnz_w(jw)
                    posCoef = posu + (posv-1)*nc_u + (posw-1)*nc_u*nc_v
                    genPosC = ju + (jv-1)*nnz_col_u + (jw-1)*nnz_col_u*nnz_col_v
                    Ci0(genPosC) = coefs(posCoef)                    
                end do
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u,nnz_col_u), Bt(nnz_row_u,nnz_col_u), Wt(nnz_col_u))
        Bt = B_u(i_nnz_u, j_nnz_u)
        Wt = W_u(row_u, j_nnz_u)
        
        forall (iu = 1 : size(Bt, 1), ju = 1 : size(Bt, 2)) 
            BW_u(iu, ju) = Bt(iu, ju) * Wt(ju)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_v(nnz_row_v,nnz_col_v), Bt(nnz_row_v,nnz_col_v), Wt(nnz_col_v))
        Bt = B_v(i_nnz_v, j_nnz_v)
        Wt = W_v(row_v, j_nnz_v)
        
        forall (iv = 1 : size(Bt, 1), jv = 1 : size(Bt, 2)) 
            BW_v(iv, jv) = Bt(iv, jv) * Wt(jv)
        end forall
        deallocate(Bt, Wt)

        allocate(BW_w(nnz_row_w,nnz_col_w), Bt(nnz_row_w,nnz_col_w), Wt(nnz_col_w))
        Bt = B_w(i_nnz_w, j_nnz_w)
        Wt = W_w(row_w, j_nnz_w)

        forall (iw = 1 : size(Bt, 1), jw = 1 : size(Bt, 2)) 
            BW_w(iw, jw) = Bt(iw, jw) * Wt(jw)
        end forall
        deallocate(Bt, Wt)

        call tensor3d_dot_vector(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, &
                                nnz_row_w, nnz_col_w, BW_u, BW_v, BW_w, Ci0, data_row)

    end subroutine csr_get_row_3d 

    subroutine csr_get_matrix_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, &
                                nnz_I_u, nnz_I_v, nnz_I_w, indi_I_u, indi_I_v, indi_I_w, indj_I_u, indj_I_v, indj_I_w, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 3D case (Ww . Bw) x (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, nnz_I_u, nnz_I_v, nnz_I_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nc_u*nc_v*nc_w)
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
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
        integer, intent(in) :: size_data_result

        double precision, intent(out) :: data_result
        dimension :: data_result(size_data_result)

        ! Local data 
        ! -----------------  
        integer :: iu, iv, iw, ju, jv, jw, k, nb_tasks, genPos, offset
        double precision, allocatable, dimension(:,:) :: B_u, B_v, B_w, W_u, W_v, W_w
        
        integer :: nnz_row_u, nnz_row_v, nnz_row_w, nnz_col_u, nnz_col_v, nnz_col_w, nnz_gen_row
        integer, allocatable, dimension(:) :: i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w
        double precision, allocatable :: data_row(:)

        ! ====================================================
        ! Initialize
        allocate(B_u(nr_u, nc_u), B_v(nr_v, nc_v), B_w(nr_w, nc_w), &
                W_u(nr_u, nc_u), W_v(nr_v, nc_v), W_w(nr_w, nc_w))
        call csr2dense(nnz_u, indi_u, indj_u, data_B_u, nr_u, nc_u, B_u)
        call csr2dense(nnz_v, indi_v, indj_v, data_B_v, nr_v, nc_v, B_v)
        call csr2dense(nnz_w, indi_w, indj_w, data_B_w, nr_w, nc_w, B_w)
        call csr2dense(nnz_u, indi_u, indj_u, data_W_u, nr_u, nc_u, W_u)
        call csr2dense(nnz_v, indi_v, indj_v, data_W_v, nr_v, nc_v, W_v)
        call csr2dense(nnz_w, indi_w, indj_w, data_W_w, nr_w, nc_w, W_w)
        ! ====================================================

        !$OMP PARALLEL PRIVATE(nnz_col_u,nnz_col_v,nnz_col_w,offset,j_nnz_u,ju,j_nnz_v,jv,j_nnz_w,jw,nnz_row_u) &
        !$OMP PRIVATE(nnz_row_v,nnz_row_w,i_nnz_u,i_nnz_v,i_nnz_w,data_row,nnz_gen_row,genPos,k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_u*nr_v*nr_w/nb_tasks)
        do iw = 1, nr_w
            do iv = 1, nr_v
                do iu = 1, nr_u

                    ! FOR ROWS
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
                    
                    ! FOR COLUMNS
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

                    nnz_gen_row = nnz_row_u * nnz_row_v * nnz_row_w
                    allocate(data_row(nnz_gen_row))

                    call csr_get_row_3d(coefs, nr_u, nc_u,  nr_v, nc_v, nr_w, nc_w, iu, iv, iw, &
                                    nnz_row_u, nnz_row_v, nnz_row_w, i_nnz_u, i_nnz_v, i_nnz_w, &
                                    nnz_col_u, nnz_col_v, nnz_col_w, j_nnz_u, j_nnz_v, j_nnz_w, &
                                    B_u, B_v, B_w, W_u, W_v, W_w, data_row)
                    deallocate(i_nnz_u, i_nnz_v, i_nnz_w, j_nnz_u, j_nnz_v, j_nnz_w)

                    ! Get offset in result 
                    genPos = iu + (iv-1)*nr_u + (iw-1)*nr_u*nr_v
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

    ! ----------------------------------------------------
    ! Functions for Fast Diagonalization method
    ! ----------------------------------------------------
    ! "Fast Diagonalization" 

    subroutine eigen_decomposition(nr, nc, Mcoef, Kcoef, nnz, indi, indj, &
                                    data_B0, data_W0, data_B1, data_W1, Method, dorobin, &
                                    eigenvalues, eigenvectors, Kdiag, Mdiag)
        !! Eigen decomposition generalized KU = MUD
        !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
        !! M: mass matrix, M = int B0 B0 dx = W00 * B0
        !! U: eigenvectors matrix
        !! D: diagonal of eigenvalues
        !! IN CSR FORMAT
        
        implicit none 
        ! Input / output 
        ! -------------------
        integer, intent(in) :: nr, nc, nnz
        double precision, dimension(*), intent(in) :: Mcoef, Kcoef
        integer, intent(in) :: indi, indj
        dimension :: indi(nr+1), indj(nnz)
        double precision, intent(in) :: data_B0, data_W0, data_B1, data_W1
        dimension :: data_B0(nnz), data_W0(nnz), data_B1(nnz), data_W1(nnz)
        character(len=10), intent(in) :: Method
        integer, intent(in) :: dorobin
        dimension :: dorobin(2)
                
        double precision, intent(out) :: eigenvalues, eigenvectors
        dimension :: eigenvalues(nr), eigenvectors(nr, nr)
        double precision, intent(out) :: Kdiag, Mdiag
        dimension :: Kdiag(nr), Mdiag(nr)

        ! Local data
        ! ----------------
        integer :: i, j, info, lwork, liwork
        double precision, allocatable, dimension(:) :: data_B0n, data_B1n
        double precision, allocatable, dimension(:,:) :: BB0, WW0, MM, BB1, WW1, KK
        double precision, allocatable, dimension(:) :: W, work, iwork
        double precision :: dummy(1)
        integer :: idum(1)
        
        ! Initialize Masse matrix
        allocate(data_B0n(nnz))
        data_B0n = data_B0
        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then 
            do i = 1, nr
                do j = indi(i), indi(i+1)-1
                    data_B0n(j) = data_B0n(j)*Mcoef(indj(j))
                end do
            end do
        end if

        allocate(BB0(nr, nc))
        call csr2dense(nnz, indi, indj, data_B0n, nr, nc, BB0)
        deallocate(data_B0n)
        allocate(WW0(nr, nc))
        call csr2dense(nnz, indi, indj, data_W0, nr, nc, WW0)
        allocate(MM(nr, nr))
        MM = matmul(WW0, transpose(BB0))
        deallocate(BB0, WW0)

        ! Initialize Stiffness matrix
        allocate(data_B1n(nnz))
        data_B1n = data_B1
        if ((Method.eq.'TDS').or.(Method.eq.'TDC')) then
            do i = 1, nr
                do j = indi(i), indi(i+1)-1
                    data_B1n(j) = data_B1n(j)*Kcoef(indj(j))
                end do
            end do
        end if

        allocate(BB1(nr, nc))
        call csr2dense(nnz, indi, indj, data_B1n, nr, nc, BB1)
        deallocate(data_B1n)
        allocate(WW1(nr, nc))
        call csr2dense(nnz, indi, indj, data_W1, nr, nc, WW1)
        allocate(KK(nr, nr))
        KK = matmul(WW1, transpose(BB1))
        deallocate(BB1, WW1)

        ! Select diagonal of M and K
        do j = 1, nr
            Kdiag(j) = KK(j, j)
            Mdiag(j) = MM(j, j)
        end do

        ! Modify K to avoid singular matrix (We consider a Robin boundary condition)
        if (dorobin(1).eq.1) then 
            KK(1, 1) = 1000 * KK(1,1)
        end if

        if (dorobin(2).eq.1) then 
            KK(nr,nr) = 1000 * KK(nr, nr)
        end if

        ! -----------------------------------
        ! Eigen decomposition KK U = MM U DD
        ! -----------------------------------
        ! Use routine workspace query to get optimal workspace.
        allocate(W(nr))
        call dsygvd(1, 'V', 'U', nr, KK, nr, MM, nr, W, dummy, -1, idum, -1, info)

        ! Make sure that there is enough workspace 
        lwork = max(1+(6+2*nr)*nr, nint(dummy(1)))
        liwork = max(3+5*nr, idum(1))
        allocate (work(lwork), iwork(liwork))

        ! Solve
        call dsygvd(1, 'V', 'U', nr, KK, nr, MM, nr, W, work, lwork, iwork, liwork, info)

        ! Get values
        eigenvectors = KK
        eigenvalues = W
        deallocate(KK, MM, W, work, iwork)

    end subroutine eigen_decomposition

    subroutine fast_diag_steady_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, &
                                        diagonal, array_in, array_out)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
        double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
        dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                        diagonal(nr_total), array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! -------------
        integer :: i, nb_tasks
        double precision :: array_temp
        dimension :: array_temp(nr_total)

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

        ! ---------------------------------
        ! Second part 
        ! ---------------------------------
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
        do i = 1, nr_total
            array_temp(i) = array_temp(i)/diagonal(i)
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        ! ----------------------------------
        ! Third part
        ! ----------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, array_out)
        
    end subroutine fast_diag_steady_3d

    subroutine fast_diag_interp_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, array_in, array_out)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
        double precision, intent(in) :: U_u, U_v, U_w, array_in
        dimension :: U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! -------------
        double precision :: array_temp
        dimension :: array_temp(nr_total)

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

        ! ----------------------------------
        ! Second part
        ! ----------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, array_out)

    end subroutine fast_diag_interp_3d

    subroutine fast_diag_transient_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, &
                                        diagonal, array_in, dt, nu, array_out)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
        double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
        dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
                        diagonal(nr_total), array_in(nr_total)
        double precision :: dt, nu

        double precision, intent(out) :: array_out
        dimension :: array_out(nr_total)

        ! Local data
        ! -------------
        integer :: i, nb_tasks
        double precision :: array_temp, diagonal_new
        dimension :: array_temp(nr_total), diagonal_new(nr_total)

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, &
                            transpose(U_u), transpose(U_v), transpose(U_w), array_in, array_temp)

        ! ---------------------------------
        ! Second part 
        ! ---------------------------------
        diagonal_new = 1.d0/dt
        if (nu.gt.0.d0) then 
            diagonal_new = diagonal_new + nu*diagonal
        end if
        !$OMP PARALLEL 
        nb_tasks = omp_get_num_threads()
        !$OMP DO SCHEDULE(STATIC, nr_total/nb_tasks)
        do i = 1, nr_total
            array_temp(i) = array_temp(i)/diagonal_new(i)
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        ! ----------------------------------
        ! Third part
        ! ----------------------------------
        call tensor3d_dot_vector(nr_u, nr_u, nr_v, nr_v, nr_w, nr_w, U_u, U_v, U_w, array_temp, array_out)
    
    end subroutine fast_diag_transient_3d

    subroutine fast_diag_static_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, &
                                        diagonal, array_in, array_out)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        use omp_lib
        implicit none
        ! Input / output  data 
        !---------------------
        integer, parameter :: d = 3
        integer, intent(in) :: nr_total, nr_u, nr_v, nr_w
        double precision, intent(in) :: U_u, U_v, U_w, diagonal, array_in
        dimension ::    U_u(nr_u, nr_u, d), U_v(nr_v, nr_v, d), U_w(nr_w, nr_w, d), &
                        diagonal(d, nr_total), array_in(d, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(d, nr_total)

        ! Local data
        ! -------------
        integer :: i
        double precision :: array_temp
        dimension :: array_temp(nr_total)

        do i = 1, d 
            call fast_diag_steady_3d(nr_total, nr_u, nr_v, nr_w, U_u(:, :, i), U_v(:, :, i), U_w(:, :, i), &
                                    diagonal(i, :), array_in(i, :), array_temp)
            array_out(i, :) = array_temp
        end do
        
    end subroutine fast_diag_static_3d

    ! For improving fast diagonalisation (TD, TDS, JM and JMS)

    subroutine tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, CC, &
                                        M_u, M_v, M_w, K_u, K_v, K_w)
        !! Tensor decomposition to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nc_total, nc_u, nc_v, nc_w
        double precision, intent(in) :: CC
        dimension :: CC(3, 3, nc_total)

        double precision, intent(inout) :: M_u, M_v, M_w, K_u, K_v, K_w
        dimension :: M_u(nc_u), M_v(nc_v), M_w(nc_w), K_u(nc_u), K_v(nc_v), K_w(nc_w)

        ! Local data
        ! ---------------
        double precision :: Vscript, Wscript, Mscript, Nscript
        dimension ::    Vscript(nc_u, nc_v, nc_w), Wscript(2, nc_u, nc_v, nc_w), &
                        Mscript(nc_u, nc_v, nc_w), Nscript(nc_u, nc_v, nc_w)
        integer :: ju, jv, jw, k, l, genPos, cont
        double precision :: vmin, vmax
        double precision :: UU(3), WW(3), WWlk(2)

        do k = 1, 3
            ! Set Vscript
            Vscript = 0.d0
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                        UU = [M_u(ju), M_v(jv), M_w(jw)]
                        Vscript(ju, jv, jw) = CC(k, k, genPos)*UU(k)/(UU(1)*UU(2)*UU(3))
                    end do
                end do
            end do

            ! Update w
            if (k.eq.1) then 
                do ju = 1, nc_u
                    vmin = minval(Vscript(ju, :, :))
                    vmax = maxval(Vscript(ju, :, :))
                    K_u(ju) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do jv = 1, nc_v
                    vmin = minval(Vscript(:, jv, :))
                    vmax = maxval(Vscript(:, jv, :))
                    K_v(jv) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do jw = 1, nc_v
                    vmin = minval(Vscript(:, :, jw))
                    vmax = maxval(Vscript(:, :, jw))
                    K_w(jw) = sqrt(vmin*vmax)
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
                    do jw = 1, nc_w
                        do jv = 1, nc_v
                            do ju = 1, nc_u
                                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                                UU = [M_u(ju), M_v(jv), M_w(jw)]
                                WW = [K_u(ju), K_v(jv), K_w(jw)]
                                Wscript(cont, ju, jv, jw) = CC(k, k, genPos)*UU(k)*UU(l)&
                                                            /(UU(1)*UU(2)*UU(3)*WW(k))
                            end do
                        end do
                    end do
                end if
            end do

            ! Compute Nscript and Mscript
            do jw = 1, nc_w
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        WWlk = Wscript(:, ju, jv, jw)
                        Nscript(ju, jv, jw) = minval(WWlk)
                        Mscript(ju, jv, jw) = maxval(WWlk)
                    end do
                end do
            end do

            ! Update u
            if (k.eq.1) then 
                do ju = 1, nc_u
                    vmin = minval(Nscript(ju, :, :))
                    vmax = maxval(Mscript(ju, :, :))
                    M_u(ju) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do jv = 1, nc_v
                    vmin = minval(Nscript(:, jv, :))
                    vmax = maxval(Mscript(:, jv, :))
                    M_v(jv) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do jw = 1, nc_v
                    vmin = minval(Nscript(:, :, jw))
                    vmax = maxval(Mscript(:, :, jw))
                    M_w(jw) = sqrt(vmin*vmax)
                end do
            end if
        end do

    end subroutine tensor_decomposition_3d

    subroutine jacobien_mean_3d(nc_u, nc_v, nc_w, nnz_J, JJ, nnz_K, KK, &
                                L1, L2, L3, lambda1, lambda2, lambda3)
        
        use omp_lib
        implicit none
        ! Input /  output data
        ! -----------------------
        integer :: step
        integer, intent(in) :: nc_u, nc_v, nc_w, nnz_J, nnz_K
        double precision, intent(in) :: JJ, KK
        dimension :: JJ(3, 3, nnz_J), KK(3, 3, nnz_K)
    
        double precision, intent(out) :: L1, L2, L3, lambda1, lambda2, lambda3
        
        ! Local data
        ! --------------
        integer :: i, j, k, nb_tasks, pos, nb_pts, nb_pts_temp
        double precision, dimension(3,3) :: dist, MatrixT, Q
        double precision, dimension(3) :: dummy
        double precision :: sq

        ! Define step
        step = min((nc_u-1)/2, (nc_v-1)/2, (nc_w-1)/2)
        
        ! Count number of quadrature points
        nb_pts = 1
        nb_pts_temp = 0
        do i = 1, nc_u, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
        nb_pts_temp = 0

        do j = 1, nc_v, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
        nb_pts_temp = 0

        do k = 1, nc_w, step
            nb_pts_temp = nb_pts_temp + 1
        end do
        nb_pts = nb_pts * nb_pts_temp
    
        ! Compute distance
        !--------------------------
        L1 = 0.d0; L2 = 0.d0; L3 = 0.d0
    
        !$OMP PARALLEL PRIVATE(pos, MatrixT, dist, dummy) REDUCTION(+:L1, L2, L3)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nnz_J/(nb_tasks*step*step*step))
        do k = 1, nc_w, step
            do j = 1, nc_v, step
                do i = 1, nc_u, step
                    pos = i + (j-1)*nc_u + (k-1)*nc_u*nc_v
                    MatrixT = JJ(:, :, pos)
                    call polar_decomposition(MatrixT, Q, dist, dummy, 1, 1)
    
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
        L1 = L1/sq; L2 = L2/sq; L3 = L3/sq

        ! Compute conductivity
        !--------------------------
        if (nnz_K.eq.1) then

            lambda1 = KK(1, 1, 1)
            lambda2 = KK(2, 2, 1)
            lambda3 = KK(3, 3, 1)

        else if (nnz_K.eq.nnz_J) then

            lambda1 = 0.d0; lambda2 = 0.d0; lambda3 = 0.d0

            !$OMP PARALLEL PRIVATE(pos) REDUCTION(+:lambda1, lambda2, lambda3)
            nb_tasks = omp_get_num_threads()
            !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nnz_J/(nb_tasks*step*step*step))
            do k = 1, nc_w, step
                do j = 1, nc_v, step
                    do i = 1, nc_u, step
                        pos = i + (j-1)*nc_u + (k-1)*nc_u*nc_v
        
                        ! Find mean of diagonal 
                        lambda1 = lambda1 + KK(1, 1, pos)/nb_pts
                        lambda2 = lambda2 + KK(2, 2, pos)/nb_pts
                        lambda3 = lambda3 + KK(3, 3, pos)/nb_pts
                    end do
                end do
            end do
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL 
            
        end if

        ! Conductivity normalized
        sq = sqrt(lambda1**2 + lambda2**2 + lambda3**2)
        lambda1 = lambda1/sq; lambda2 = lambda2/sq; lambda3 = lambda3/sq
        
    end subroutine jacobien_mean_3d

    ! For scaling (TDS and JMS)
    
    subroutine find_parametric_diag_3d(nr_u, nr_v, nr_w, c_u, c_v, c_w, &
                                        Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w, diag)
        !! Find the diagonal of the preconditioner "fast diagonalization"
                            
        implicit none
        ! Input / output data
        ! -------------------------
        integer, intent(in) :: nr_u, nr_v, nr_w
        double precision, intent(in) :: c_u, c_v, c_w
        double precision, intent(in) :: Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
        dimension :: Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), &
                    Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w)

        double precision, intent(out) :: diag
        dimension :: diag(nr_u*nr_v*nr_w)

        ! Initialize
        diag = 0.d0

        ! Find K3 M2 M1
        call kron_product_3vec(nr_w, Kdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, diag, c_w)

        ! Find M3 K2 M1
        call kron_product_3vec(nr_w, Mdiag_w, nr_v, Kdiag_v, nr_u, Mdiag_u, diag, c_v)

        ! Find M3 M2 K1
        call kron_product_3vec(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Kdiag_u, diag, c_u)

    end subroutine find_parametric_diag_3d

    subroutine find_physical_diag_3d(coefs, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w, diag)
        !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
        !! Algotihm based on sum factorization adapted to diagonal case 
        !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
        !! Indices must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
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
        ! ------------------
        integer :: iu, iv, iw, ju, jv, jw, nb_tasks, Cpos, Ipos, offset
        double precision :: sum1, sum2, sum3

        integer :: nnz_row_u, nnz_row_v, nnz_row_w
        integer, allocatable, dimension(:) :: indj_nnz_u, indj_nnz_v, indj_nnz_w
        double precision, allocatable, dimension(:) :: data_nnz_B_u, data_nnz_B_v, data_nnz_B_w, &
                                                        data_nnz_W_u, data_nnz_W_v, data_nnz_W_w        

        !$OMP PARALLEL PRIVATE(ju,jv,jw,nnz_row_u,nnz_row_v,nnz_row_w,offset,indj_nnz_u,data_nnz_B_u,data_nnz_W_u) &
        !$OMP PRIVATE(indj_nnz_v,data_nnz_B_v,data_nnz_W_v,indj_nnz_w,data_nnz_B_w,data_nnz_W_w,sum1,sum2,sum3,Cpos,Ipos)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_u*nr_v*nr_w /nb_tasks) 
        do iw = 1, nr_w
            do iv = 1, nr_v
                do iu = 1, nr_u

                    ! Number of nonzeros
                    nnz_row_u = indi_u(iu+1) - indi_u(iu)
                    nnz_row_v = indi_v(iv+1) - indi_v(iv)
                    nnz_row_w = indi_w(iw+1) - indi_w(iw)

                    ! Set values
                    allocate(indj_nnz_u(nnz_row_u), data_nnz_B_u(nnz_row_u), data_nnz_W_u(nnz_row_u))
                    offset = indi_u(iu)
                    do ju = 1, nnz_row_u
                        indj_nnz_u(ju) = indj_u(ju+offset-1)
                        data_nnz_B_u(ju) = data_B_u(ju+offset-1)
                        data_nnz_W_u(ju) = data_W_u(ju+offset-1)
                    end do

                    allocate(indj_nnz_v(nnz_row_v), data_nnz_B_v(nnz_row_v), data_nnz_W_v(nnz_row_v))
                    offset = indi_v(iv)
                    do jv = 1, nnz_row_v
                        indj_nnz_v(jv) = indj_v(jv+offset-1)
                        data_nnz_B_v(jv) = data_B_v(jv+offset-1)
                        data_nnz_W_v(jv) = data_W_v(jv+offset-1)
                    end do

                    allocate(indj_nnz_w(nnz_row_w), data_nnz_B_w(nnz_row_w), data_nnz_W_w(nnz_row_w))
                    offset = indi_w(iw)
                    do jw = 1, nnz_row_w
                        indj_nnz_w(jw) = indj_w(jw+offset-1)
                        data_nnz_B_w(jw) = data_B_w(jw+offset-1)
                        data_nnz_W_w(jw) = data_W_w(jw+offset-1)
                    end do

                    sum3 = 0.d0
                    do jw = 1, nnz_row_w
                        sum2 = 0.d0
                        do jv = 1, nnz_row_v
                            sum1 = 0.d0
                            do ju = 1, nnz_row_u
                                Cpos = indj_nnz_u(ju) + (indj_nnz_v(jv)-1)*nc_u + (indj_nnz_w(jw)-1)*nc_u*nc_v
                                sum1 = sum1 + data_nnz_W_u(ju)*data_nnz_B_u(ju)*coefs(Cpos)
                            end do
                            sum2 = sum2 + data_nnz_W_v(jv)*data_nnz_B_v(jv)*sum1
                        end do
                        sum3 = sum3 + data_nnz_W_w(jw)*data_nnz_B_w(jw)*sum2
                    end do 

                    ! General position
                    Ipos = iu + (iv-1)*nr_u + (iw-1)*nr_u*nr_w
                    
                    ! Update diagonal
                    diag(Ipos) = sum3

                    deallocate(indj_nnz_u, data_nnz_B_u, data_nnz_W_u)
                    deallocate(indj_nnz_v, data_nnz_B_v, data_nnz_W_v)
                    deallocate(indj_nnz_w, data_nnz_B_w, data_nnz_W_w)

                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine find_physical_diag_3d

end module tensor_methods