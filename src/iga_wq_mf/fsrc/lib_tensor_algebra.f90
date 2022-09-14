! =========================
! module :: tensor algebra
! author :: Joaquin Cornejo
!
! This modules is intended to have all the tensor algebra functions necessary for the good
! performance of the rest of functions.
! =========================

! ---------------
! Tensor algebra 
! ---------------

subroutine tensor_n_mode_product_dM(nc_u, nc_v, nc_w, X, nr, nc, U, mode, nrR, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = (nc_u, nc_v, nc_w)
    !! Matrix U = (nr, nc)
    !! Tensor R = (nu, nv, nw) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw
    !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u

    use omp_lib
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nc_u, nc_v, nc_w, nr, nc, mode, nrR
    double precision, intent(in) :: X, U
    dimension :: X(nc_u*nc_v*nc_w), U(nr, nc)

    double precision, intent(out) :: R
    dimension :: R(nrR)

    ! Local data
    ! ----------
    double precision, allocatable, dimension(:,:) :: Xt, Rt
    integer :: ju, jv, jw, i, nb_tasks

    if (mode.eq.1) then 

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

    else if (mode.eq.2) then 

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
        
    else if (mode.eq.3) then 

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

end subroutine tensor_n_mode_product_dM

subroutine tensor_n_mode_product_spM(nc_u, nc_v, nc_w, X, nrU, nnzU, dataU, indi, indj, mode, nrR, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = (nc_u, nc_v, nc_w)
    !! Matrix U = (nr, nc). Since U is in CSR format, nc is not necessary to be declared
    !! Tensor R = (nu, nv, nw) (It depends on 'mode'). It is mandatory that nrR = nu*nv*nw
    !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u

    use omp_lib
    implicit none
    ! Input / output data 
    ! -------------------
    integer, intent(in) :: nc_u, nc_v, nc_w, nrU, nnzU, mode, nrR
    integer, intent(in) :: indi, indj
    dimension :: indi(nrU+1), indj(nnzU)
    double precision, intent(in) :: X, dataU
    dimension :: X(nc_u*nc_v*nc_w), dataU(nnzU)

    double precision, intent(out) :: R
    dimension :: R(nrR)

    ! Local data
    ! ----------
    integer :: ju, jv, jw, i, jX, jR, k, nb_tasks
    double precision :: s

    if (mode.eq.1) then 
        !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_w*nc_v/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                jX = (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                jR = (jv-1)*nrU + (jw-1)*nrU*nc_v
                do i = 1, nrU
                    s = 0.d0
                    do k = indi(i), indi(i+1)-1
                        s = s + dataU(k) * X(indj(k) + jX)
                    end do
                    R(i + jR) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.2) then 
        !$OMP PARALLEL PRIVATE(jX, jR, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_w*nc_u/nb_tasks) 
        do jw = 1, nc_w
            do ju = 1, nc_u
                jX = ju + (jw-1)*nc_u*nc_v
                jR = ju + (jw-1)*nc_u*nrU
                do i = 1, nrU
                    s = 0.d0
                    do k = indi(i), indi(i+1)-1
                        s = s + dataU(k) * X((indj(k)-1)*nc_u + jX)
                    end do
                    R((i-1)*nc_u + jR) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    else if (mode.eq.3) then 
        !$OMP PARALLEL PRIVATE(jX, i, s, k)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, nc_v*nc_u/nb_tasks) 
        do jv = 1, nc_v
            do ju = 1, nc_u
                jX = ju + (jv-1)*nc_u
                jR = jX
                do i = 1, nrU
                    s = 0.d0
                    do k = indi(i), indi(i+1)-1
                        s = s + dataU(k) * X((indj(k)-1)*nc_u*nc_v + jX)
                    end do
                    R((i-1)*nc_u*nc_v + jR) = s
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
    end if

end subroutine tensor_n_mode_product_spM

subroutine sumproduct2d_dM(nr_u, nc_u, nr_v, nc_v, Mu, Mv, vector_in, vector_out)
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
    double precision, intent(in) :: vector_in
    dimension :: vector_in(nc_u*nc_v)

    double precision, intent(out) :: vector_out
    dimension :: vector_out(nr_u*nr_v)

    ! Local data 
    ! ----------
    double precision, allocatable, dimension(:) :: R1

    ! First product
    allocate(R1(nr_u*nc_v))
    call tensor_n_mode_product_dM(nc_u, nc_v, 1, vector_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    ! Second product
    call tensor_n_mode_product_dM(nr_u, nc_v, 1, R1, nr_v, nc_v, Mv, 2, size(vector_out), vector_out)
    deallocate(R1)

end subroutine sumproduct2d_dM

subroutine sumproduct3d_dM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, Mu, Mv, Mw, vector_in, vector_out)
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
    double precision, intent(in) :: vector_in
    dimension :: vector_in(nc_u*nc_v*nc_w)

    double precision, intent(out) :: vector_out
    dimension :: vector_out(nr_u*nr_v*nr_w)

    ! Local data 
    ! ----------
    double precision, allocatable, dimension(:) :: R1, R2

    ! First product
    allocate(R1(nr_u*nc_v*nc_w))
    call tensor_n_mode_product_dM(nc_u, nc_v, nc_w, vector_in, nr_u, nc_u, Mu, 1, size(R1), R1)

    ! Second product
    allocate(R2(nr_u*nr_v*nc_w))
    call tensor_n_mode_product_dM(nr_u, nc_v, nc_w, R1, nr_v, nc_v, Mv, 2, size(R2), R2)
    deallocate(R1)

    ! Third product
    call tensor_n_mode_product_dM(nr_u, nr_v, nc_w, R2, nr_w, nc_w, Mw, 3, size(vector_out), vector_out)
    deallocate(R2)

end subroutine sumproduct3d_dM

subroutine sumproduct2d_spM(nr_u, nc_u, nr_v, nc_v, nnz_u, indi_u, indj_u, data_u, &
                            nnz_v, indi_v, indj_v, data_v, vector_in, vector_out)
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
    ! ----------
    double precision, allocatable, dimension(:) :: R1

    ! First product
    allocate(R1(nr_u*nc_v))
    call tensor_n_mode_product_spM(nc_u, nc_v, 1, vector_in, nr_u, nnz_u, &
                                data_u, indi_u, indj_u, 1, size(R1), R1)

    ! Second product
    call tensor_n_mode_product_spM(nr_u, nc_v, 1, R1, nr_v, nnz_v, &
                                data_v, indi_v, indj_v, 2, size(vector_out), vector_out)
    deallocate(R1)

end subroutine sumproduct2d_spM

subroutine sumproduct3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, indi_u, indj_u, data_u, &
                        nnz_v, indi_v, indj_v, data_v, nnz_w, indi_w, indj_w, data_w, &
                        vector_in, vector_out)
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
    ! ----------
    double precision, allocatable, dimension(:) :: R1, R2

    ! First product
    allocate(R1(nr_u*nc_v*nc_w))
    call tensor_n_mode_product_spM(nc_u, nc_v, nc_w, vector_in, nr_u, nnz_u, &
                                data_u, indi_u, indj_u, 1, size(R1), R1)

    ! Second product
    allocate(R2(nr_u*nr_v*nc_w))
    call tensor_n_mode_product_spM(nr_u, nc_v, nc_w, R1, nr_v, nnz_v, &
                                data_v, indi_v, indj_v, 2, size(R2), R2)
    deallocate(R1)

    ! Third product
    call tensor_n_mode_product_spM(nr_u, nr_v, nc_w, R2, nr_w, nnz_w, &
                                data_w, indi_w, indj_w, 3, size(vector_out), vector_out)
    deallocate(R2)

end subroutine sumproduct3d_spM

! --------------------------------
! Sum product to compute matrices 
! --------------------------------

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

    ! Initiliaze
    allocate(Ci0(nnz_col_u*nnz_col_v))

    ! Set values of C
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
    call sumproduct2d_dM(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, BW_u, BW_v, Ci0, data_row)

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

    ! Initialize
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

    ! Initiliaze
    allocate(Ci0(nnz_col_u*nnz_col_v*nnz_col_w))

    ! Set values of C
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
    call sumproduct3d_dM(nnz_row_u, nnz_col_u, nnz_row_v, nnz_col_v, &
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

    ! Initialize
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
    !! Algotihm based on sum factorization adapted to diagonal case 
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
    integer :: iu, iv, iw, ju, jv, jw, posCoef, genPos, offset, nb_tasks
    double precision :: sum_u, sum_v, sum_w

    integer :: nnz_col_u, nnz_col_v, nnz_col_w
    integer, allocatable, dimension(:) :: j_nnz_u, j_nnz_v, j_nnz_w
    double precision, allocatable, dimension(:) :: B_nnz_u, B_nnz_v, B_nnz_w, W_nnz_u, W_nnz_v, W_nnz_w        

    !$OMP PARALLEL PRIVATE(ju,jv,jw,nnz_col_u,nnz_col_v,nnz_col_w,offset,j_nnz_u,B_nnz_u,W_nnz_u) &
    !$OMP PRIVATE(j_nnz_v,B_nnz_v,W_nnz_v,j_nnz_w,B_nnz_w,W_nnz_w,sum_u,sum_v,sum_w,posCoef,genPos)
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(DYNAMIC, nr_u*nr_v*nr_w /nb_tasks) 
    do iw = 1, nr_w
        do iv = 1, nr_v
            do iu = 1, nr_u

                ! Set values 
                nnz_col_u = indi_u(iu+1) - indi_u(iu)
                allocate(j_nnz_u(nnz_col_u), B_nnz_u(nnz_col_u), W_nnz_u(nnz_col_u))
                offset = indi_u(iu)
                do ju = 1, nnz_col_u
                    j_nnz_u(ju) = indj_u(ju+offset-1)
                    B_nnz_u(ju) = data_B_u(ju+offset-1)
                    W_nnz_u(ju) = data_W_u(ju+offset-1)
                end do

                nnz_col_v = indi_v(iv+1) - indi_v(iv)
                allocate(j_nnz_v(nnz_col_v), B_nnz_v(nnz_col_v), W_nnz_v(nnz_col_v))
                offset = indi_v(iv)
                do jv = 1, nnz_col_v
                    j_nnz_v(jv) = indj_v(jv+offset-1)
                    B_nnz_v(jv) = data_B_v(jv+offset-1)
                    W_nnz_v(jv) = data_W_v(jv+offset-1)
                end do

                nnz_col_w = indi_w(iw+1) - indi_w(iw)
                allocate(j_nnz_w(nnz_col_w), B_nnz_w(nnz_col_w), W_nnz_w(nnz_col_w))
                offset = indi_w(iw)
                do jw = 1, nnz_col_w
                    j_nnz_w(jw) = indj_w(jw+offset-1)
                    B_nnz_w(jw) = data_B_w(jw+offset-1)
                    W_nnz_w(jw) = data_W_w(jw+offset-1)
                end do

                ! Sum factorization
                sum_w = 0.d0
                do jw = 1, nnz_col_w
                    sum_v = 0.d0
                    do jv = 1, nnz_col_v
                        sum_u = 0.d0
                        do ju = 1, nnz_col_u
                            posCoef = j_nnz_u(ju) + (j_nnz_v(jv)-1)*nc_u + (j_nnz_w(jw)-1)*nc_u*nc_v
                            sum_u = sum_u + W_nnz_u(ju)*B_nnz_u(ju)*coefs(posCoef)
                        end do
                        sum_v = sum_v + W_nnz_v(jv)*B_nnz_v(jv)*sum_u
                    end do
                    sum_w = sum_w + W_nnz_w(jw)*B_nnz_w(jw)*sum_v
                end do 

                ! Update diagonal
                genPos = iu + (iv - 1)*nr_u + (iw - 1)*nr_u*nr_w
                diag(genPos) = sum_w

                deallocate(j_nnz_u, B_nnz_u, W_nnz_u)
                deallocate(j_nnz_v, B_nnz_v, W_nnz_v)
                deallocate(j_nnz_w, B_nnz_w, W_nnz_w)

            end do
        end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

end subroutine csr_get_diag_3d

! ----------------------------
! Fast Diagonalization method
! ----------------------------
subroutine eigen_decomposition(nr, nc, Mcoef, Kcoef, nnz, indi, indj, &
                                data_B0, data_W0, data_B1, data_W1, robin_condition, &
                                eigenvalues, eigenvectors, Kdiag, Mdiag)
    !! Eigen decomposition generalized KU = MUD
    !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
    !! M: mass matrix, M = int B0 B0 dx = W00 * B0
    !! U: eigenvectors matrix
    !! D: diagonal of eigenvalues
    !! IN CSR FORMAT
    
    implicit none 
    ! Input / output data
    ! -------------------
    double precision, parameter :: penalty = 1001
    integer, intent(in) :: nr, nc, nnz
    double precision, intent(in) :: Mcoef, Kcoef
    dimension :: Mcoef(nc), Kcoef(nc)
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
    integer :: i, j, info, lwork, liwork, idum(1)
    double precision, allocatable, dimension(:) :: data_B0t, data_B1t
    double precision, allocatable, dimension(:,:) :: BB0, WW0, MM, BB1, WW1, KK
    double precision, allocatable, dimension(:) :: work, iwork
    double precision :: dummy(1)
    
    ! Initialize masse matrix
    allocate(data_B0t(nnz))
    data_B0t = data_B0
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_B0t(j) = data_B0t(j)*Mcoef(indj(j))
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

    ! Initialize stiffness matrix
    allocate(data_B1t(nnz))
    data_B1t = data_B1
    do i = 1, nr
        do j = indi(i), indi(i+1)-1
            data_B1t(j) = data_B1t(j)*Kcoef(indj(j))
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

    ! Modify K to avoid singular matrix (Like a Robin boundary condition)
    if (robin_condition(1).eq.1) then 
        KK(1, 1) = penalty*KK(1,1)
    end if

    if (robin_condition(2).eq.1) then 
        KK(nr,nr) = penalty*KK(nr, nr)
    end if

    ! Save diagonal of M and K
    do j = 1, nr
        Kdiag(j) = KK(j, j)
        Mdiag(j) = MM(j, j)
    end do

    ! -----------------------------------
    ! Eigen decomposition KK U = MM U DD
    ! -----------------------------------
    ! Use routine workspace query to get optimal workspace
    call dsygvd(1, 'V', 'L', nr, KK, nr, MM, nr, eigenvalues, dummy, -1, idum, -1, info)

    ! Make sure that there is enough workspace 
    lwork = max(1+(6+2*nr)*nr, nint(dummy(1)))
    liwork = max(3+5*nr, idum(1))
    allocate (work(lwork), iwork(liwork))

    ! Get eigen decomposition
    eigenvectors = KK
    call dsygvd(1, 'V', 'L', nr, eigenvectors, nr, MM, nr, eigenvalues, work, lwork, iwork, liwork, info)
    deallocate(KK, MM, work, iwork)

end subroutine eigen_decomposition

subroutine find_parametric_diag_3d(nr_u, nr_v, nr_w, Mdiag_u, Mdiag_v, Mdiag_w, &
                                    Kdiag_u, Kdiag_v, Kdiag_w, c_u, c_v, c_w, diag)
    !! Find the diagonal of the preconditioner in fast diagonalization method
                        
    implicit none
    ! Input / output data
    ! -------------------
    integer, intent(in) :: nr_u, nr_v, nr_w
    double precision, intent(in) :: Mdiag_u, Mdiag_v, Mdiag_w, Kdiag_u, Kdiag_v, Kdiag_w
    dimension :: Mdiag_u(nr_u), Mdiag_v(nr_v), Mdiag_w(nr_w), &
                Kdiag_u(nr_u), Kdiag_v(nr_v), Kdiag_w(nr_w)
    double precision :: c_u, c_v, c_w

    double precision, intent(out) :: diag
    dimension :: diag(nr_u*nr_v*nr_w)

    ! Initialize
    diag = 0.d0

    ! Find M3 x M2 x K1
    call kron_product_3vec(nr_w, Mdiag_w, nr_v, Mdiag_v, nr_u, Kdiag_u, diag, c_u)

    ! Find M3 x K2 x M1
    call kron_product_3vec(nr_w, Mdiag_w, nr_v, Kdiag_v, nr_u, Mdiag_u, diag, c_v)

    ! Find K3 x M2 x M1
    call kron_product_3vec(nr_w, Kdiag_w, nr_v, Mdiag_v, nr_u, Mdiag_u, diag, c_w)
    
end subroutine find_parametric_diag_3d

subroutine tensor_decomposition_3d(nc_total, nc_u, nc_v, nc_w, CC, &
                                    M_u, M_v, M_w, K_u, K_v, K_w)
    !! Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
    !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

    implicit none
    ! Input /  output data
    ! --------------------
    integer, intent(in) :: nc_total, nc_u, nc_v, nc_w
    double precision, intent(in) :: CC
    dimension :: CC(3, 3, nc_total)

    double precision, intent(inout) :: M_u, M_v, M_w, K_u, K_v, K_w
    dimension :: M_u(nc_u), M_v(nc_v), M_w(nc_w), K_u(nc_u), K_v(nc_v), K_w(nc_w)

    ! Local data
    ! ----------
    double precision :: Vscript, Wscript, Mscript, Nscript
    dimension ::    Vscript(nc_u, nc_v, nc_w), Wscript(2, nc_u, nc_v, nc_w), &
                    Mscript(nc_u, nc_v, nc_w), Nscript(nc_u, nc_v, nc_w)
    integer :: ju, jv, jw, k, l, genPos, c
    double precision :: vmin, vmax
    double precision :: UU(3), WW(3), WWlk(2)

    do k = 1, 3
        ! Set Vscript
        do jw = 1, nc_w
            do jv = 1, nc_v
                do ju = 1, nc_u
                    genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                    UU = [M_u(ju), M_v(jv), M_w(jw)] 
                    Vscript(ju, jv, jw) = CC(k, k, genPos)*UU(k)/product(UU)
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

        c = 0
        do l = 1, 3
            if (k.ne.l) then 
                c = c + 1
                ! Set Wscript
                do jw = 1, nc_w
                    do jv = 1, nc_v
                        do ju = 1, nc_u
                            genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                            UU = [M_u(ju), M_v(jv), M_w(jw)]
                            WW = [K_u(ju), K_v(jv), K_w(jw)]
                            Wscript(c, ju, jv, jw) = CC(k, k, genPos)*UU(k)*UU(l)&
                                                        /(product(UU)*WW(k))
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

subroutine jacobien_mean_3d(nc_u, nc_v, nc_w, nnz, JJ, Lu, Lv, Lw)
    !! This method aims to calculate the average "deformation" that 
    !! the hypercube undergoes through the topological transformation F 
    !! (from parametric to physical space) 
    !! We suppose that this transformation is a composition of a rotation and a stretching, then 
    !! one can apply polar decomposition to the jacobian matrix of F.
    !! THIS FUNCTION IS DEPRECATED
    
    implicit none
    ! Input /  output data
    ! --------------------
    integer :: step
    integer, intent(in) :: nc_u, nc_v, nc_w, nnz
    double precision, intent(in) :: JJ
    dimension :: JJ(3, 3, nnz)

    double precision, intent(out) :: Lu, Lv, Lw
    
    ! Local data
    ! ----------
    integer :: i, j, k, l, nb_pts, nb_pts_temp
    double precision, dimension(3,3) :: dist, JJtemp, Q
    double precision, dimension(3) :: dummy

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

    ! Initialize    
    Lu = 0.d0; Lv = 0.d0; Lw = 0.d0
    do k = 1, nc_w, step
        do j = 1, nc_v, step
            do i = 1, nc_u, step
                l = i + (j-1)*nc_u + (k-1)*nc_u*nc_v
                JJtemp = JJ(:, :, l)
                call polar_decomposition(JJtemp, Q, dist, dummy, .true., .true.)

                ! Find mean of diagonal of jacobien
                Lu = Lu + dist(1, 1)/nb_pts
                Lv = Lv + dist(2, 2)/nb_pts
                Lw = Lw + dist(3, 3)/nb_pts
            end do
        end do
    end do

end subroutine jacobien_mean_3d

subroutine conductivity_mean_3d(nc_u, nc_v, nc_w, nnz, KK, lamb_u, lamb_v, lamb_w)
    !! This function is DEPRECATED. 

    implicit none
    ! Input /  output data
    ! --------------------
    integer :: step
    integer, intent(in) :: nc_u, nc_v, nc_w, nnz
    double precision, intent(in) :: KK
    dimension :: KK(3, 3, nnz)

    double precision, intent(out) :: lamb_u, lamb_v, lamb_w
    
    ! Local data
    ! ----------
    integer :: i, j, k, pos, nb_pts, nb_pts_temp

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

    if (nnz.eq.1) then

        lamb_u = KK(1, 1, 1)
        lamb_v = KK(2, 2, 1)
        lamb_w = KK(3, 3, 1)

    else 

        ! Initialize
        lamb_u = 0.d0; lamb_v = 0.d0; lamb_w = 0.d0

        do k = 1, nc_w, step
            do j = 1, nc_v, step
                do i = 1, nc_u, step
                    pos = i + (j-1)*nc_u + (k-1)*nc_u*nc_v
    
                    ! Find mean of diagonal 
                    lamb_u = lamb_u + KK(1, 1, pos)/nb_pts
                    lamb_v = lamb_v + KK(2, 2, pos)/nb_pts
                    lamb_w = lamb_w + KK(3, 3, pos)/nb_pts
                end do
            end do
        end do
        
    end if

end subroutine conductivity_mean_3d

subroutine compute_mean_3d(nc_u, nc_v, nc_w, coefs, integral)
    !! Compute the "mean" of a property. This is a factor that improves fast diagonalization method 

    implicit none
    ! Input /  output data
    ! --------------------
    integer, intent(in) :: nc_u, nc_v, nc_w
    double precision, intent(in) :: coefs
    dimension :: coefs(nc_u*nc_v*nc_w)

    double precision, intent(out) :: integral
    
    ! Local data
    ! ----------
    integer :: i, j, k, genPos, pos
    integer :: ind_u, ind_v, ind_w
    dimension :: ind_u(3), ind_v(3), ind_w(3)
    double precision :: coefs_temp
    dimension :: coefs_temp(3, 3, 3)

    ! Initialize
    pos = int((nc_u+1)/2); ind_u = (/1, pos, nc_u/)
    pos = int((nc_v+1)/2); ind_v = (/1, pos, nc_v/)
    pos = int((nc_w+1)/2); ind_w = (/1, pos, nc_w/)

    ! Create new coefficients
    coefs_temp = 0.d0
    do k = 1, 3
        do j = 1, 3
            do i = 1, 3
                genPos = ind_u(i) + (ind_v(j) - 1)*nc_u + (ind_w(k) - 1)*nc_u*nc_v
                coefs_temp(i, j, k) = coefs(genPos)
            end do
        end do
    end do

    ! Compute integral
    call trapezoidal_rule_3d(3, 3, 3, coefs_temp, integral)

end subroutine compute_mean_3d