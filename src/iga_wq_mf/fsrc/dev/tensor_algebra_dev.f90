! =========================================
! module :: sum product (Tensor operations) ! ALOT OF CHANGES
! author :: Joaquin Cornejo
! =========================================
module tensor_methods

    contains

    subroutine tensor_n_mode_product(I1, I2, I3, X, I, J, U, n, Ir1, Ir2, Ir3, R)
        !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
        !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
        !! Tensor X = (I1, I2, I3)
        !! Matrix U = (I, J)
        !! Tensor R = (Ir1, Ir2, Ir3) (It depends on 'n')

        implicit none
        ! Input / output data 
        ! -------------------- 
        integer, intent(in) :: I1, I2, I3, I, J, n, Ir1, Ir2, Ir3
        double precision, intent(in) :: X, U
        dimension :: X(I1*I2*I3), U(I*J)

        double precision, intent(out) ::  R(Ir1*Ir2*Ir3)

        ! Local data
        ! ---------------
        integer :: genPosX, genPosU, genPosR
        integer :: ii1, ii2, ii3, ii
        double precision :: sum

        if (n.eq.1) then 
            do ii3 = 1, I3
                do ii2 = 1, I2
                    do ii = 1, I
                        sum = 0.0d0
                        do ii1 = 1, I1
                            genPosX = ii1 + (ii2-1)*I1 + (ii3-1)*I1*I2
                            genPosU = ii + (ii1-1)*I
                            sum = sum + X(genPosX) * U(genPosU)
                        end do
                        genPosR = ii + (ii2-1)*I + (ii3-1)*I*I2
                        R(genPosR) = sum
                    end do
                end do
            end do
        else if (n.eq.2) then 
            do ii3 = 1, I3
                do ii = 1, I
                    do ii1 = 1, I1
                        sum = 0.0d0
                        do ii2 = 1, I2
                            genPosX = ii1 + (ii2-1)*I1 + (ii3-1)*I1*I2
                            genPosU = ii + (ii2-1) * I
                            sum = sum + X(genPosX) * U(genPosU)
                        end do
                        genPosR = ii1 + (ii-1)*I1 + (ii3-1)*I1*I
                        R(genPosR) = sum
                    end do
                end do
            end do
        else if (n.eq.3) then 
            do ii = 1, I
                do ii2 = 1, I2
                    do ii1 = 1, I1
                        sum = 0.0d0
                        do ii3 = 1, I3
                            genPosX = ii1 + (ii2-1)*I1 + (ii3-1)*I1*I2
                            genPosU = ii + (ii3-1) * I
                            sum = sum + X(genPosX) * U(genPosU)
                        end do
                        genPosR = ii1 + (ii2-1)*I1 + (ii-1)*I1*I2
                        R(genPosR) = sum
                    end do
                end do
            end do
        end if

    end subroutine tensor_n_mode_product

    subroutine rankone2d_dot_vector(size_u, size_v, Vu, Vv, X0, result)
        !! Evaluates a dot product between a 2D rank-one tensor and a vector 
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
        !! result = (Vv x Vu) . X0 (x = tensor prod, . = dot product)
        !! Vector Vu = (size_u)
        !! Vector Vv = (size_v)
        !! Vector X0 = (size_u * size_v)

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  size_u, size_v
        double precision, intent(in) :: Vu, Vv
        dimension :: Vu(size_u), Vv(size_v)
        double precision, intent(in) :: X0
        dimension :: X0(size_u*size_v)

        double precision, intent(out) :: result

        ! Local data 
        ! -------------
        integer :: ju, jv, posGen
        double precision :: sum

        ! Initialize
        result = 0.0d0 
        do ju = 1, size_u
            sum = 0.0d0
            do jv = 1, size_v
                posGen = ju + (jv-1)*size_u
                sum = sum + X0(posGen) * Vv(jv)
            end do
            result = result + sum * Vu(ju)
        end do

    end subroutine rankone2d_dot_vector

    subroutine rankone3d_dot_vector(size_u, size_v, size_w, &
                                Vu, Vv, Vw, X0, result)
        !! Evaluates a dot product between a 3D rank-one tensor and a vector 
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
        !! result = (Vw x Vv x Vu) . X0 (x = tensor prod, . = dot product)
        !! Vector Vu = (size_u)
        !! Vector Vv = (size_v)
        !! Vector Vw = (size_w)
        !! Vector X0 = (size_u * size_v * size_w)

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  size_u, size_v, size_w
        double precision, intent(in) :: X0
        dimension :: X0(size_u*size_v*size_w)
        double precision, intent(in) :: Vu, Vv, Vw
        dimension :: Vu(size_u), Vv(size_v), Vw(size_w)

        double precision, intent(out) :: result

        ! Local data 
        ! -------------
        integer :: ju, jv, jw, posGen
        double precision :: sum1, sum2

        ! Initialize
        result = 0.0d0 
        do ju = 1, size_u
            sum2 = 0.0d0
            do jv = 1, size_v
                sum1 = 0.0d0
                do jw = 1, size_w
                    posGen = ju + (jv-1)*size_u + (jw-1)*size_u*size_v
                    sum1 = sum1 + X0(posGen) * Vw(jw)
                end do
                sum2 = sum2 + sum1 * Vv(jv)
            end do
            result = result + sum2 * Vu(ju)
        end do

    end subroutine rankone3d_dot_vector

    subroutine tensor2d_dot_vector( nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, &
                                    Mu, Mv, vector_in, vector_out)
        !! Evaluates a dot product between a tensor 2D and a vector
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
        !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Vector_in = (nb_cols_u * nb_cols_v)

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v)
        double precision, intent(in) :: Mu, Mv
        dimension ::    Mu(nb_rows_u, nb_cols_u), &
                        Mv(nb_rows_v, nb_cols_v)

        double precision, intent(inout) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v)

        ! Local data 
        ! -------------
        integer :: genPos_out
        double precision :: sum
        integer :: iu, iv, nb_tasks

        !$OMP PARALLEL PRIVATE(genPos_out, sum)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nb_rows_u*nb_rows_v/nb_tasks) 
        do iv = 1, nb_rows_v
            do iu = 1, nb_rows_u
                ! General position
                genPos_out = iu + (iv-1)*nb_rows_u 

                call rankone2d_dot_vector(nb_cols_u, nb_cols_v, Mu(iu, :), Mv(iv, :), vector_in, sum)

                ! Update vector
                vector_out(genPos_out) = sum
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    end subroutine tensor2d_dot_vector

    subroutine tensor3d_dot_vector(nb_rows_u, nb_cols_u, &
                                    nb_rows_v, nb_cols_v, &
                                    nb_rows_w, nb_cols_w, &
                                    Mu, Mv, Mw, vector_in, vector_out)
        !! Evaluates a dot product between a tensor 3D and a vector 
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
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
        dimension ::    Mu(nb_rows_u, nb_cols_u), &
                        Mv(nb_rows_v, nb_cols_v), &
                        Mw(nb_rows_w, nb_cols_w)

        double precision, intent(out) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data 
        ! -------------
        integer :: genPos_out
        double precision :: sum
        integer :: iu, iv, iw, nb_tasks

        !$OMP PARALLEL PRIVATE(genPos_out, sum) 
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u * nb_rows_v * nb_rows_w /nb_tasks) 
        do iw = 1, nb_rows_w
            do iv = 1, nb_rows_v
                do iu = 1, nb_rows_u
                    ! General position
                    genPos_out = iu + (iv-1)*nb_rows_u + (iw-1)*nb_rows_u*nb_rows_v

                    call rankone3d_dot_vector(nb_cols_u, nb_cols_v, nb_cols_w, &
                                        Mu(iu, :), Mv(iv, :), Mw(iw, :), vector_in, sum)

                    ! Update vector
                    vector_out(genPos_out) = sum
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    end subroutine tensor3d_dot_vector

    subroutine tensor2d_sparsedot_vector(nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        size_data_u, indi_u, indj_u, data_u, &
                                        size_data_v, indi_v, indj_v, data_v, &
                                        vector_in, vector_out)
        !! Evaluates a dot product between a tensor 2D and a vector 
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
        !! Vector_out = (Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Vector_in = (nb_cols_u * nb_cols_v)
        !! Indexes must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v)
        integer, intent(in) :: size_data_u, size_data_v
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), &
                        indj_u(size_data_u), indj_v(size_data_v)
        double precision, intent(in) :: data_u, data_v
        dimension ::    data_u(size_data_u), data_v(size_data_v)

        double precision, intent(inout) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v)

        ! Local data 
        ! -------------
        ! Loops
        integer :: offset, iu, iv, ju, jv, nb_tasks
        integer :: genPos_in, genPos_out, genPos_tensor

        ! Eval product
        double precision :: sum

        ! Select row
        integer :: nnz_u, nnz_v
        integer, allocatable, dimension(:) :: indj_nnz_u, indj_nnz_v
        double precision, allocatable, dimension(:) :: data_nnz_u, data_nnz_v, tensor
        integer :: dummy_var

        ! Initiliaze
        dummy_var = nb_cols_v

        !$OMP PARALLEL PRIVATE(nnz_u,nnz_v,offset,ju,jv,indj_nnz_u,data_nnz_u,indj_nnz_v,data_nnz_v) &
        !$OMP PRIVATE(tensor,genPos_tensor,genPos_in,genPos_out,sum)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nb_rows_u * nb_rows_v /nb_tasks) 
        do iv = 1, nb_rows_v
            do iu = 1, nb_rows_u
                ! General position
                genPos_out = iu + (iv-1)*nb_rows_u 

                ! Number of nonzeros
                nnz_u = indi_u(iu+1) - indi_u(iu)
                nnz_v = indi_v(iv+1) - indi_v(iv)

                ! Set values of row in Mu
                allocate(indj_nnz_u(nnz_u), data_nnz_u(nnz_u))
                indj_nnz_u = 0
                data_nnz_u = 0.0d0
                offset = indi_u(iu)
                do ju = 1, nnz_u
                    indj_nnz_u(ju) = indj_u(ju+offset-1)
                    data_nnz_u(ju) = data_u(ju+offset-1)
                end do

                ! Set values of row in Mv
                allocate(indj_nnz_v(nnz_v), data_nnz_v(nnz_v))
                indj_nnz_v = 0
                data_nnz_v = 0.0d0
                offset = indi_v(iv)
                do jv = 1, nnz_v
                    indj_nnz_v(jv) = indj_v(jv+offset-1)
                    data_nnz_v(jv) = data_v(jv+offset-1)
                end do

                allocate(tensor(nnz_u*nnz_v))
                tensor = 0.0d0
                do jv = 1, nnz_v
                    do ju = 1, nnz_u
                        genPos_tensor = ju + (jv-1)*nnz_u 
                        genPos_in = indj_nnz_u(ju) + (indj_nnz_v(jv)-1)*nb_cols_u 
                        tensor(genPos_tensor) = vector_in(genPos_in)
                    end do
                end do

                call rankone2d_dot_vector(nnz_u, nnz_v, data_nnz_u, data_nnz_v, tensor, sum)

                ! Update vector
                vector_out(genPos_out) = vector_out(genPos_out) + sum

                deallocate(indj_nnz_u, indj_nnz_v)
                deallocate(data_nnz_u, data_nnz_v)
                deallocate(tensor)
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    end subroutine tensor2d_sparsedot_vector

    subroutine tensor3d_sparsedot_vector(nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        nb_rows_w, nb_cols_w, &
                                        size_data_u, indi_u, indj_u, data_u, &
                                        size_data_v, indi_v, indj_v, data_v, &
                                        size_data_w, indi_w, indj_w, data_w, &
                                        vector_in, vector_out)
        !! Evaluates a dot product between a tensor 3D and a vector 
        !! Based on "Matrix-free weighted quadrature for a computationally efficient" by Sangalli and Tani
        !! Vector_out = (Mw x Mv x Mu) . Vector_in (x = tensor prod, . = dot product)
        !! Matrix Mu = (nb_rows_u, nb_cols_u)
        !! Matrix Mv = (nb_rows_v, nb_cols_v)
        !! Matrix Mw = (nb_rows_w, nb_cols_w)
        !! Vector_in = (nb_cols_u * nb_cols_v * nb_cols_w)
        !! Indexes must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: vector_in
        dimension :: vector_in(nb_cols_u*nb_cols_v*nb_cols_w)
        integer, intent(in) :: size_data_u, size_data_v, size_data_w
        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nb_rows_u+1), indi_v(nb_rows_v+1), indi_w(nb_rows_w+1), &
                        indj_u(size_data_u), indj_v(size_data_v), indj_w(size_data_w)
        double precision, intent(in) :: data_u, data_v, data_w
        dimension ::    data_u(size_data_u), data_v(size_data_v), data_w(size_data_w)

        double precision, intent(inout) :: vector_out
        dimension :: vector_out(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data 
        ! -------------
        ! Loops
        integer :: offset, iu, iv, iw, ju, jv, jw, nb_tasks

        ! Get X where vec(X) = vector
        integer :: genPos_in, genPos_out, genPos_tensor

        ! Eval product
        double precision :: sum

        ! Select row
        integer :: nnz_u, nnz_v, nnz_w
        integer, allocatable, dimension(:) :: indj_nnz_u, indj_nnz_v, indj_nnz_w
        double precision, allocatable, dimension(:) :: data_nnz_u, data_nnz_v, data_nnz_w, tensor
        integer :: dummy_var

        ! Initiliaze
        dummy_var = nb_cols_w

        !$OMP PARALLEL PRIVATE(nnz_u,nnz_v,nnz_w,offset,ju,jv,jw,indj_nnz_u,data_nnz_u,indj_nnz_v) &
        !$OMP PRIVATE(data_nnz_v,indj_nnz_w,data_nnz_w,tensor,genPos_tensor,genPos_in,genPos_out,sum)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u * nb_rows_v * nb_rows_w /nb_tasks) 
        do iw = 1, nb_rows_w
            do iv = 1, nb_rows_v
                do iu = 1, nb_rows_u

                    ! General position
                    genPos_out = iu + (iv-1)*nb_rows_u + (iw-1)*nb_rows_u*nb_rows_v

                    ! Number of nonzeros
                    nnz_u = indi_u(iu+1) - indi_u(iu)
                    nnz_v = indi_v(iv+1) - indi_v(iv)
                    nnz_w = indi_w(iw+1) - indi_w(iw)

                    ! Set values
                    allocate(indj_nnz_u(nnz_u), data_nnz_u(nnz_u))
                    offset = indi_u(iu)
                    do ju = 1, nnz_u
                        indj_nnz_u(ju) = indj_u(ju+offset-1)
                        data_nnz_u(ju) = data_u(ju+offset-1)
                    end do

                    allocate(indj_nnz_v(nnz_v), data_nnz_v(nnz_v))
                    offset = indi_v(iv)
                    do jv = 1, nnz_v
                        indj_nnz_v(jv) = indj_v(jv+offset-1)
                        data_nnz_v(jv) = data_v(jv+offset-1)
                    end do

                    allocate(indj_nnz_w(nnz_w), data_nnz_w(nnz_w))
                    offset = indi_w(iw)
                    do jw = 1, nnz_w
                        indj_nnz_w(jw) = indj_w(jw+offset-1)
                        data_nnz_w(jw) = data_w(jw+offset-1)
                    end do

                    allocate(tensor(nnz_u*nnz_v*nnz_w))
                    tensor = 0.0d0
                    do jw = 1, nnz_w
                        do jv = 1, nnz_v
                            do ju = 1, nnz_u
                                genPos_tensor = ju + (jv-1)*nnz_u + (jw-1)*nnz_u*nnz_v
                                genPos_in = indj_nnz_u(ju) + (indj_nnz_v(jv)-1)*nb_cols_u + (indj_nnz_w(jw)-1)*nb_cols_u*nb_cols_v
                                tensor(genPos_tensor) = vector_in(genPos_in)
                            end do
                        end do
                    end do

                    call rankone3d_dot_vector(nnz_u, nnz_v, nnz_w, data_nnz_u, data_nnz_v, data_nnz_w, tensor, sum)

                    ! Update vector
                    vector_out(genPos_out) = vector_out(genPos_out) + sum

                    deallocate(indj_nnz_u, indj_nnz_v, indj_nnz_w)
                    deallocate(data_nnz_u, data_nnz_v, data_nnz_w)
                    deallocate(tensor)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 

    end subroutine tensor3d_sparsedot_vector

    ! ----------------------------------------------------
    ! Sum product to evaluate matrices 
    ! ----------------------------------------------------
    subroutine csr_get_row_2d(coefficients, &
                            nb_rows_u, nb_cols_u, &
                            nb_rows_v, nb_cols_v, &
                            row_u, row_v, &
                            nnz_row_u, nnz_row_v, &
                            i_nnz_u, i_nnz_v, &
                            nnz_col_u, nnz_col_v, &
                            j_nnz_u, j_nnz_v, &
                            B_u, B_v, W_u, W_v, &
                            data_row)
        !! Computes a row of a matrix constructed with row_u, row_v 

        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: coefficients
        dimension :: coefficients(nb_cols_u*nb_cols_v)
        integer, intent(in) :: row_u, row_v
        integer, intent(in) :: nnz_row_u, nnz_row_v
        integer, intent(in) :: nnz_col_u, nnz_col_v
        integer, intent(in) :: i_nnz_u, i_nnz_v, j_nnz_u, j_nnz_v
        dimension :: i_nnz_u(nnz_row_u), i_nnz_v(nnz_row_v), j_nnz_u(nnz_col_u), j_nnz_v(nnz_col_v)
        double precision, intent(in) :: B_u, B_v, W_u, W_v
        dimension ::    B_u(nb_rows_u, nb_cols_u), &   
                        B_v(nb_rows_v, nb_cols_v), &
                        W_u(nb_rows_u, nb_cols_u), &
                        W_v(nb_rows_v, nb_cols_v)

        double precision, intent(out) :: data_row(nnz_row_u*nnz_row_v)

        ! Local data 
        ! ----------------- 
        ! Create Ci
        integer :: pos_coef
        double precision, allocatable :: Ci0(:), Ci1(:)
        double precision, allocatable :: BW_u(:), BW_v(:)

        ! Loops
        integer :: genPosBWl, genPosC
        integer :: iu, iv, ju, jv, posu, posv

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v))

        ! Set values of C
        do jv = 1, nnz_col_v
            do ju = 1, nnz_col_u
                posu = j_nnz_u(ju)
                posv = j_nnz_v(jv)
                pos_coef = posu + (posv-1)*nb_cols_u 
                genPosC = ju + (jv-1)*nnz_col_u 
                Ci0(genPosC) = coefficients(pos_coef)                    
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u*nnz_col_u))
        allocate(BW_v(nnz_row_v*nnz_col_v))

        do ju = 1, nnz_col_u
            do iu = 1, nnz_row_u
                genPosBWl = iu + (ju-1)*nnz_row_u
                BW_u(genPosBWl) = B_u(i_nnz_u(iu), j_nnz_u(ju)) * W_u(row_u, j_nnz_u(ju))
            end do
        end do

        do jv = 1, nnz_col_v
            do iv = 1, nnz_row_v
                genPosBWl = iv + (jv-1)*nnz_row_v
                BW_v(genPosBWl) = B_v(i_nnz_v(iv), j_nnz_v(jv)) * W_v(row_v, j_nnz_v(jv))
            end do
        end do

        ! Evaluate tensor product
        allocate(Ci1(nnz_col_u*nnz_row_v))
        call tensor_n_mode_product(nnz_col_u, nnz_col_v, 1, Ci0, & 
                    nnz_row_v, nnz_col_v, BW_v, 2, &
                    nnz_col_u, nnz_row_v, 1, Ci1)
        deallocate(Ci0, BW_v)

        call tensor_n_mode_product(nnz_col_u, nnz_row_v, 1, Ci1, & 
                    nnz_row_u, nnz_col_u, BW_u, 1, &
                    nnz_row_u, nnz_row_v, 1, data_row)
        deallocate(Ci1, BW_u)

        deallocate(Ci0, Ci1)
        deallocate(BW_u, BW_v)

    end subroutine csr_get_row_2d 

    subroutine csr_get_matrix_2d(coefficients, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                size_data_u, size_data_v, &
                                indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_W_u, &
                                data_B_v, data_W_v, &
                                size_data_I_u, size_data_I_v, &
                                indi_I_u, indi_I_v, & 
                                indj_I_u, indj_I_v, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 2D case (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indexes must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v
        double precision, intent(in) :: coefficients
        dimension :: coefficients(nb_cols_u*nb_cols_v)
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

        integer, intent(in) :: indi_result(nb_rows_u*nb_rows_v+1) 
        integer, intent(in) :: size_data_result
        double precision, intent(inout) :: data_result(size_data_result)

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
        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nb_rows_u * nb_rows_v /nb_tasks)
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

                call csr_get_row_2d(coefficients, nb_rows_u, nb_cols_u, &
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
                    data_result(result_offset+j-1) = data_result(result_offset+j-1) + data_row(j)
                end do   

                deallocate(data_row)
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        deallocate(B_u, B_v, W_u, W_v)

    end subroutine csr_get_matrix_2d 

    subroutine csr_get_row_3d(coefficients, &
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
        double precision, intent(in) :: coefficients
        dimension :: coefficients(nb_cols_u*nb_cols_v*nb_cols_w)
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

        double precision, intent(out) :: data_row(nnz_row_u*nnz_row_v*nnz_row_w)

        ! Local data 
        ! ----------------- 
        ! Create Ci
        integer :: pos_coef
        double precision, allocatable :: Ci0(:), Ci1(:), Ci2(:)

        ! Create Bl and Wl
        double precision, allocatable :: BW_u(:), BW_v(:), BW_w(:)

        ! Loops
        integer :: genPosBWl, genPosC
        integer :: iu, iv, iw, ju, jv, jw, posu, posv, posw

        ! Initiliaze
        allocate(Ci0(nnz_col_u*nnz_col_v*nnz_col_w))

        ! Set values of C
        do jw = 1, nnz_col_w
            do jv = 1, nnz_col_v
                do ju = 1, nnz_col_u
                    posu = j_nnz_u(ju)
                    posv = j_nnz_v(jv)
                    posw = j_nnz_w(jw)
                    pos_coef = posu + (posv-1)*nb_cols_u + (posw-1)*nb_cols_u*nb_cols_v
                    genPosC = ju + (jv-1)*nnz_col_u + (jw-1)*nnz_col_u*nnz_col_v
                    Ci0(genPosC) = coefficients(pos_coef)                    
                end do
            end do
        end do

        ! Set values of BW
        allocate(BW_u(nnz_row_u*nnz_col_u))
        allocate(BW_v(nnz_row_v*nnz_col_v))
        allocate(BW_w(nnz_row_w*nnz_col_w))

        do ju = 1, nnz_col_u
            do iu = 1, nnz_row_u
                genPosBWl = iu + (ju-1)*nnz_row_u
                BW_u(genPosBWl) = B_u(i_nnz_u(iu), j_nnz_u(ju)) * W_u(row_u, j_nnz_u(ju))
            end do
        end do

        do jv = 1, nnz_col_v
            do iv = 1, nnz_row_v
                genPosBWl = iv + (jv-1)*nnz_row_v
                BW_v(genPosBWl) = B_v(i_nnz_v(iv), j_nnz_v(jv)) * W_v(row_v, j_nnz_v(jv))
            end do
        end do

        do jw = 1, nnz_col_w
            do iw = 1, nnz_row_w
                genPosBWl = iw + (jw-1)*nnz_row_w
                BW_w(genPosBWl) = B_w(i_nnz_w(iw), j_nnz_w(jw)) * W_w(row_w, j_nnz_w(jw))
            end do
        end do

        ! Evaluate tensor product
        allocate(Ci1(nnz_col_u*nnz_col_v*nnz_row_w))
        call tensor_n_mode_product(nnz_col_u, nnz_col_v, nnz_col_w, Ci0, & 
                                nnz_row_w, nnz_col_w, BW_w, 3, &
                                nnz_col_u, nnz_col_v, nnz_row_w, Ci1)

        allocate(Ci2(nnz_col_u*nnz_row_v*nnz_row_w))
        call tensor_n_mode_product(nnz_col_u, nnz_col_v, nnz_row_w, Ci1, & 
                                nnz_row_v, nnz_col_v, BW_v, 2, &
                                nnz_col_u, nnz_row_v, nnz_row_w, Ci2)

        call tensor_n_mode_product(nnz_col_u, nnz_row_v, nnz_row_w, Ci2, & 
                                nnz_row_u, nnz_col_u, BW_u, 1, &
                                nnz_row_u, nnz_row_v, nnz_row_w, data_row)

        deallocate(Ci0, Ci1, Ci2)
        deallocate(BW_u, BW_v, BW_w)

    end subroutine csr_get_row_3d 

    subroutine csr_get_matrix_3d(coefficients, &
                                nb_rows_u, nb_cols_u, &
                                nb_rows_v, nb_cols_v, &
                                nb_rows_w, nb_cols_w, &
                                size_data_u, size_data_v, size_data_w, &
                                indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_B_u, data_W_u, &
                                data_B_v, data_W_v, &
                                data_B_w, data_W_w, &
                                size_data_I_u, size_data_I_v, size_data_I_w, &
                                indi_I_u, indi_I_v, indi_I_w, & 
                                indj_I_u, indj_I_v, indj_I_w, &
                                indi_result, size_data_result, data_result)
        !! Computes a matrix in 3D case (Ww . Bw) x (Wv . Bv) x (Wu . Bu)
        !! x: kronecker product and .: inner product
        !! Indexes must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! ------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: coefficients
        dimension :: coefficients(nb_cols_u*nb_cols_v*nb_cols_w)
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
        integer, intent(in) :: indi_result(nb_rows_u*nb_rows_v*nb_rows_w+1) 
        integer, intent(in) :: size_data_result

        double precision, intent(inout) :: data_result(size_data_result)

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
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u * nb_rows_v * nb_rows_w /nb_tasks)
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

                    call csr_get_row_3d(coefficients, nb_rows_u, nb_cols_u, &
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
                        data_result(result_offset+j-1) = data_result(result_offset+j-1) + data_row(j)
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
    ! Jacobi diagonal
    ! ----------------------------------------------------
    subroutine find_diagonal_fd_3d(nb_rows_u, nb_rows_v, nb_rows_w, L1, L2, L3, &
                                Mdiag_u, Mdiag_v, Mdiag_w, &
                                Kdiag_u, Kdiag_v, Kdiag_w, diag)
        implicit none

        ! Input / output data
        ! -------------------------
        integer, intent(in) :: nb_rows_u, nb_rows_v, nb_rows_w
        double precision, intent(in) :: L1, L2, L3
        double precision, intent(in) :: Mdiag_u, Mdiag_v, Mdiag_w, &
                                        Kdiag_u, Kdiag_v, Kdiag_w
        dimension :: Mdiag_u(nb_rows_u), Mdiag_v(nb_rows_v), Mdiag_w(nb_rows_w), &
                    Kdiag_u(nb_rows_u), Kdiag_v(nb_rows_v), Kdiag_w(nb_rows_w)

        double precision, intent(out) :: diag
        dimension :: diag(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data
        ! -----------------------
        double precision, dimension(:), allocatable :: diag_temp1, diag_temp2

        ! Initialize
        diag = 0.0d0

        ! Find K3 M2 M1
        allocate(diag_temp1(nb_rows_u*nb_rows_v))
        allocate(diag_temp2(nb_rows_u*nb_rows_v*nb_rows_w))
        call vector_kron_vector(nb_rows_v, Mdiag_v, nb_rows_u, Mdiag_u, diag_temp1)
        call vector_kron_vector(nb_rows_w, Kdiag_w, size(diag_temp1), diag_temp1, diag_temp2)
        diag = diag + L1*L2/L3*diag_temp2

        ! Find M3 K2 M1
        call vector_kron_vector(nb_rows_v, Kdiag_v, nb_rows_u, Mdiag_u, diag_temp1)
        call vector_kron_vector(nb_rows_w, Mdiag_w, size(diag_temp1), diag_temp1, diag_temp2)
        diag = diag + L3*L1/L2*diag_temp2

        ! Find M3 M2 K1
        call vector_kron_vector(nb_rows_v, Mdiag_v, nb_rows_u, Kdiag_u, diag_temp1)
        call vector_kron_vector(nb_rows_w, Mdiag_w, size(diag_temp1), diag_temp1, diag_temp2)
        diag = diag + L3*L1/L2*diag_temp2

        deallocate(diag_temp1, diag_temp2)

    end subroutine find_diagonal_fd_3d

    subroutine find_jacobi_diagonal_3d(coefs, &
                                        nb_rows_u, nb_cols_u, &
                                        nb_rows_v, nb_cols_v, &
                                        nb_rows_w, nb_cols_w, &
                                        size_data_u, size_data_v, size_data_w, &
                                        indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                        data_B_u, data_B_v, data_B_w, &
                                        data_W_u, data_W_v, data_W_w, &
                                        jacobi_diag)
        !! Find diagonal without constructing all the matrix (WQ-IGA Analysis)
        !! Algotihm based on sum factorization adapted to diagonal case 
        !! See more in "Efficient matrix computation for tensor-product isogeometric analysis" by G. Sanaglli et al.
        !! Indexes must be in CSR format

        use omp_lib
        implicit none 
        ! Input / output 
        ! -------------------
        integer, intent(in) ::  nb_rows_u, nb_cols_u, nb_rows_v, nb_cols_v, nb_rows_w, nb_cols_w
        double precision, intent(in) :: coefs
        dimension :: coefs(nb_cols_u*nb_cols_v*nb_cols_w)
        integer, intent(in) ::  size_data_u, size_data_v, size_data_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nb_rows_u+1), indj_u(size_data_u), &
                        indi_v(nb_rows_v+1), indj_v(size_data_v), &
                        indi_w(nb_rows_w+1), indj_w(size_data_w)
        double precision, intent(in) :: data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w
        dimension ::    data_B_u(size_data_u), data_B_v(size_data_v), data_B_w(size_data_w), &
                        data_W_u(size_data_u), data_W_v(size_data_v), data_W_w(size_data_w)

        double precision, intent(out) :: jacobi_diag
        dimension :: jacobi_diag(nb_rows_u*nb_rows_v*nb_rows_w)

        ! Local data
        ! ------------------
        integer :: offset, nb_tasks
        integer :: iu, iv, iw, ju, jv, jw, Cpos, Ipos
        double precision :: sum1, sum2, sum3

        integer :: nnz_u, nnz_v, nnz_w
        integer, allocatable, dimension(:) :: indj_nnz_u, indj_nnz_v, indj_nnz_w
        double precision, allocatable, dimension(:) :: data_nnz_B_u, data_nnz_B_v, data_nnz_B_w, &
                                                        data_nnz_W_u, data_nnz_W_v, data_nnz_W_w        
        
        ! Initialize
        jacobi_diag = 0.0d0

        !$OMP PARALLEL PRIVATE(ju,jv,jw,nnz_u,nnz_v,nnz_w,offset,indj_nnz_u,data_nnz_B_u,data_nnz_W_u) &
        !$OMP PRIVATE(indj_nnz_v,data_nnz_B_v,data_nnz_W_v,indj_nnz_w,data_nnz_B_w,data_nnz_W_w,sum1,sum2,sum3,Cpos,Ipos)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u * nb_rows_v * nb_rows_w /nb_tasks) 
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

                    sum3 = 0.0d0
                    do jw = 1, nnz_w
                        sum2 = 0.0d0
                        do jv = 1, nnz_v
                            sum1 = 0.0d0
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
                    jacobi_diag(Ipos) = sum3

                    deallocate(indj_nnz_u, data_nnz_B_u, data_nnz_W_u)
                    deallocate(indj_nnz_v, data_nnz_B_v, data_nnz_W_v)
                    deallocate(indj_nnz_w, data_nnz_B_w, data_nnz_W_w)

                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine find_jacobi_diagonal_3d

    ! ----------------------------------------------------
    ! Sum product for fast diagonalization
    ! ----------------------------------------------------
    subroutine eigen_decomposition(nb_rows, nb_cols, &
                                    Mcoef, Kcoef, size_data, indi, indj, &
                                    data_B0, data_W00, data_B1, data_W11, &
                                    kmod1, kmod2, eigenvalues, eigenvectors, Kdiag, Mdiag)
        !! Eigen decomposition generalized KU = MUD
        !! K: stiffness matrix, K = int B1 B1 dx = W11 * B1
        !! M: mass matrix, M = int B0 B0 dx = W00 * B0
        !! U: eigenvectors matrix
        !! D: diagonal of eigenvalues
        
        implicit none 
        ! Input / output 
        ! -------------------
        integer, intent(in) ::  nb_rows, nb_cols
        double precision, intent(in) :: Mcoef, Kcoef
        dimension :: Mcoef(nb_cols), Kcoef(nb_cols)
        integer, intent(in) ::  size_data
        integer, intent(in) ::  indi, indj
        dimension :: indi(size_data), indj(size_data)
        double precision, intent(in) :: data_B0, data_W00, data_B1, data_W11
        dimension ::    data_B0(size_data), data_W00(size_data), &
                        data_B1(size_data), data_W11(size_data)
        integer, intent(in) :: kmod1, kmod2
                
        double precision, intent(out) :: eigenvalues, eigenvectors
        dimension :: eigenvalues(nb_rows), eigenvectors(nb_rows, nb_rows)
        double precision, intent(out) :: Kdiag, Mdiag
        dimension :: Kdiag(nb_rows), Mdiag(nb_rows)

        ! Local data
        ! ----------------
        double precision, allocatable, dimension(:,:) :: BB0, WW0, MM, BB1, WW1, KK
        integer :: i, j

        ! Use Lapack
        external :: dlamch, dsygvx
        integer, parameter  :: ITYPE = 1 ! Type A x = B x D where in this proble A =  KK and B = MM
        character, parameter :: JOBZ = 'V' ! Computes eigen values and vectors
        character, parameter :: UPLO = 'U' ! Consider upper triangle of A and B
        integer :: N, LDA, LDB ! Size of A and B

        double precision, allocatable :: W(:) ! Array of eigen values
        integer :: LWORk, LIWORK
        double precision, allocatable :: WORK(:), IWORK(:)
        double precision :: dummy(1)
        integer :: idum(1)
        integer :: INFO
        
        ! ====================================================
        ! Initialize Masse matrix
        allocate(BB0(nb_rows, nb_cols), &   
                WW0(nb_rows, nb_cols))
        call coo2matrix(size_data, indi, indj, data_B0, nb_rows, nb_cols, BB0)
        call coo2matrix(size_data, indi, indj, data_W00, nb_rows, nb_cols, WW0)
        allocate(MM(nb_rows, nb_rows))
        do j = 1, nb_cols
            do i = 1, nb_rows
                BB0(i, j) = BB0(i, j) * Mcoef(j)
            end do
        end do
        MM = matmul(WW0, transpose(BB0))
        deallocate(BB0, WW0)

        ! Initialize Stiffness matrix
        allocate(BB1(nb_rows, nb_cols), &   
                WW1(nb_rows, nb_cols))
        call coo2matrix(size_data, indi, indj, data_B1, nb_rows, nb_cols, BB1)
        call coo2matrix(size_data, indi, indj, data_W11, nb_rows, nb_cols, WW1)
        allocate(KK(nb_rows, nb_rows))
        do j = 1, nb_cols
            do i = 1, nb_rows
                BB1(i, j) = BB1(i, j) * Kcoef(j)
            end do
        end do
        KK = matmul(WW1, transpose(BB1))
        deallocate(BB1, WW1)

        ! Modify K to avoid singular matrix (We consider a Robin boundary condition)
        if (kmod1.eq.1) then 
            KK(1, 1) = 100 * KK(1, 1)
        end if

        if (kmod2.eq.1) then 
            KK(nb_rows, nb_rows) = 100 * KK(nb_rows, nb_rows)
        end if

        ! Select diagonal of M and K
        do j = 1, nb_rows
            Kdiag(j) = KK(j, j)
            Mdiag(j) = MM(j, j)
        end do

        ! ====================================================
        ! Eigen decomposition KK U = MM U DD
        N = nb_rows
        LDA = nb_rows
        LDB = nb_rows

        ! Set eigen vectors and eigen values
        allocate(W(N))

        ! Use routine workspace query to get optimal workspace.
        LWORk = -1
        LIWORK = -1
        call dsygvd(ITYPE, JOBZ, UPLO, N, KK, LDA, MM, LDB, W, dummy, LWORk, idum, LIWORK, INFO)

        ! Make sure that there is enough workspace 
        LWORK = max(1+(6+2*n)*n, nint(dummy(1)))
        LIWORK = max(3+5*n, idum(1))
        allocate (WORK(LWORk), IWORK(LIWORK))

        ! Solve
        call dsygvd(ITYPE, JOBZ, UPLO, N, KK, LDA, MM, LDB, W, WORK, LWORk, IWORK, LIWORK, INFO)

        ! Get values
        eigenvectors = KK
        eigenvalues = W

        deallocate(KK, MM, W, WORK, IWORK)
    end subroutine eigen_decomposition

    subroutine fast_diagonalization_2d(nb_ctrlpts_1, nb_ctrlpts_2, &
                                        U_1, D_1, U_2, D_2, array_input, array_output)
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani

        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nb_ctrlpts_1, nb_ctrlpts_2
        double precision, intent(in) :: U_1, D_1, U_2, D_2
        dimension :: U_1(nb_ctrlpts_1, nb_ctrlpts_1), D_1(nb_ctrlpts_1), &
                    U_2(nb_ctrlpts_2, nb_ctrlpts_2), D_2(nb_ctrlpts_2)

        double precision, intent(in) :: array_input
        dimension :: array_input(nb_ctrlpts_1*nb_ctrlpts_2)

        double precision, intent(out) :: array_output
        dimension :: array_output(nb_ctrlpts_1*nb_ctrlpts_2)

        ! Local data
        ! -------------
        integer :: i
        double precision, allocatable, dimension(:) :: array_temp_1, diagonal_temp, diagonal
        double precision, allocatable, dimension(:) :: Ident1, Ident2

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        allocate(array_temp_1(nb_ctrlpts_1*nb_ctrlpts_2))
        array_temp_1 = 0.0d0

        call tensor2d_dot_vector(nb_ctrlpts_1, nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_2, &
                                    transpose(U_1), transpose(U_2), array_input, array_temp_1)

        ! ---------------------------------
        ! Second part 
        ! ---------------------------------
        ! Define identities
        allocate(Ident1(nb_ctrlpts_1), Ident2(nb_ctrlpts_2))
        Ident1 = 1.0d0
        Ident2 = 1.0d0

        allocate(diagonal(nb_ctrlpts_1*nb_ctrlpts_2))
        allocate(diagonal_temp(nb_ctrlpts_1*nb_ctrlpts_2))
        diagonal = 0.0d0

        call vector_kron_vector(nb_ctrlpts_2, Ident2, nb_ctrlpts_1, D_1, diagonal_temp)
        diagonal = diagonal + diagonal_temp

        call vector_kron_vector(nb_ctrlpts_2, D_2, nb_ctrlpts_1, Ident1, diagonal_temp)
        diagonal = diagonal + diagonal_temp

        deallocate(Ident1, Ident2, diagonal_temp)

        do i = 1, size(array_temp_1)
            array_temp_1(i) = array_temp_1(i) / diagonal(i)
        end do

        ! ----------------------------------
        ! Third part
        ! ---------------------------------
        array_output = 0.0d0
        call tensor2d_dot_vector(nb_ctrlpts_1, nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_2, &
                                    U_1, U_2, array_temp_1, array_output)
        deallocate(array_temp_1)

    end subroutine fast_diagonalization_2d

    subroutine fast_diagonalization_3d(nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_3, &
                                        U_1, D_1, U_2, D_2, U_3, D_3, L1, L2, L3, array_input, array_output)
        
        !! Fast diagonalization based on "Isogeometric preconditionners based on fast solvers for the Sylvester equations"
        !! by G. Sanaglli and M. Tani
        
        implicit none
        ! Input / output  data 
        !---------------------
        integer, intent(in) :: nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_3
        double precision, intent(in) :: U_1, D_1, U_2, D_2, U_3, D_3
        dimension :: U_1(nb_ctrlpts_1, nb_ctrlpts_1), D_1(nb_ctrlpts_1), &
                    U_2(nb_ctrlpts_2, nb_ctrlpts_2), D_2(nb_ctrlpts_2), &
                    U_3(nb_ctrlpts_3, nb_ctrlpts_3), D_3(nb_ctrlpts_3)
        double precision, intent(in) :: L1, L2, L3

        double precision, intent(in) :: array_input
        dimension :: array_input(nb_ctrlpts_1*nb_ctrlpts_2*nb_ctrlpts_3)

        double precision, intent(out) :: array_output
        dimension :: array_output(nb_ctrlpts_1*nb_ctrlpts_2*nb_ctrlpts_3)

        ! Local data
        ! -------------
        integer :: i
        double precision :: array_temp_1
        dimension :: array_temp_1(nb_ctrlpts_1*nb_ctrlpts_2*nb_ctrlpts_3)
        double precision, allocatable, dimension(:) :: diagonal_temp1, diagonal_temp2, diagonal
        double precision, allocatable, dimension(:) :: Ident1, Ident2, Ident3

        ! ---------------------------------
        ! First part 
        ! ---------------------------------
        call tensor3d_dot_vector(nb_ctrlpts_1, nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_2, nb_ctrlpts_3, nb_ctrlpts_3, &
                            transpose(U_1), transpose(U_2), transpose(U_3), array_input, array_temp_1)

        ! ---------------------------------
        ! Second part 
        ! ---------------------------------
        ! Define identities
        allocate(Ident1(nb_ctrlpts_1), Ident2(nb_ctrlpts_2), Ident3(nb_ctrlpts_3))
        Ident1 = 1.0d0
        Ident2 = 1.0d0
        Ident3 = 1.0d0

        allocate(diagonal(nb_ctrlpts_1*nb_ctrlpts_2*nb_ctrlpts_3))
        allocate(diagonal_temp1(nb_ctrlpts_1*nb_ctrlpts_2))
        allocate(diagonal_temp2(nb_ctrlpts_1*nb_ctrlpts_2*nb_ctrlpts_3))
        diagonal = 0.0d0

        call vector_kron_vector(nb_ctrlpts_2, Ident2, nb_ctrlpts_1, D_1, diagonal_temp1)
        call vector_kron_vector(nb_ctrlpts_3, Ident3, size(diagonal_temp1), diagonal_temp1, diagonal_temp2)
        diagonal = diagonal + L2*L3/L1*diagonal_temp2

        call vector_kron_vector(nb_ctrlpts_2, D_2, nb_ctrlpts_1, Ident1, diagonal_temp1)
        call vector_kron_vector(nb_ctrlpts_3, Ident3, size(diagonal_temp1), diagonal_temp1, diagonal_temp2)
        diagonal = diagonal + L1*L3/L2*diagonal_temp2

        call vector_kron_vector(nb_ctrlpts_2, Ident2, nb_ctrlpts_1, Ident1, diagonal_temp1)
        call vector_kron_vector(nb_ctrlpts_3, D_3, size(diagonal_temp1), diagonal_temp1, diagonal_temp2)
        diagonal = diagonal + L1*L2/L3*diagonal_temp2

        deallocate(Ident1, Ident2, diagonal_temp1, diagonal_temp2)

        do i = 1, size(array_temp_1)
            array_temp_1(i) = array_temp_1(i) / diagonal(i)
        end do

        ! ----------------------------------
        ! Third part
        ! ---------------------------------
        array_output = 0.0d0
        call tensor3d_dot_vector(nb_ctrlpts_1, nb_ctrlpts_1, nb_ctrlpts_2, nb_ctrlpts_2, nb_ctrlpts_3, nb_ctrlpts_3, &
                                U_1, U_2, U_3, array_temp_1, array_output)

    end subroutine fast_diagonalization_3d

    ! ----------------------------------------------------
    ! Tensor decomposition
    ! ----------------------------------------------------
    subroutine diagonal_decomposition_2d(nb_qp_u, nb_qp_v, CC, &
                                        U_u, U_v, W_u, W_v)
        !! Tensor decomposition to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini
        
        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_qp_u, nb_qp_v
        double precision, intent(in) :: CC
        dimension :: CC(2, nb_qp_u*nb_qp_v)

        double precision, intent(inout) :: U_u, U_v, W_u, W_v
        dimension ::    U_u(nb_qp_u), U_v(nb_qp_v), &
                        W_u(nb_qp_u), W_v(nb_qp_v)

        ! Local data
        ! ---------------
        double precision :: Vscript(nb_qp_u, nb_qp_v), & 
                            Wscript(nb_qp_u, nb_qp_v)
        ! Loop
        integer :: k, l, i1, i2
        integer :: genpos
        double precision :: vmin, vmax
        double precision :: UU(2), WW(2)

        do k = 1, 2
            ! Set Vscript
            Vscript = 0.0d0
            do i2 = 1, nb_qp_v
                do i1 = 1, nb_qp_u
                    genpos = i1 + (i2-1)*nb_qp_u 
                    UU = [U_u(i1), U_v(i2)]
                    Vscript(i1, i2) = CC(k, genpos)*UU(k)/(UU(1)*UU(2))
                end do
            end do

            ! Update w
            if (k.eq.1) then 
                do i1 = 1, nb_qp_u
                    vmin = minval(Vscript(i1, :))
                    vmax = maxval(Vscript(i1, :))
                    W_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_qp_v
                    vmin = minval(Vscript(:, i2))
                    vmax = maxval(Vscript(:, i2))
                    W_v(i2) = sqrt(vmin*vmax)
                end do
            end if
        end do

        do k = 1, 2
            ! Initialize
            Wscript = 0.0d0
            do l = 1, 2
                if (k.ne.l) then 
                    ! Set Wscript
                    do i2 = 1, nb_qp_v
                        do i1 = 1, nb_qp_u
                            genpos = i1 + (i2-1)*nb_qp_u 
                            UU = [U_u(i1), U_v(i2)]
                            WW = [W_u(i1), W_v(i2)]
                            Wscript(i1, i2) = CC(k, genpos)*UU(k)*UU(l)&
                                                        /(UU(1)*UU(2)*WW(k))
                        end do
                    end do
                end if
            end do

            ! Update u
            if (k.eq.1) then 
                do i1 = 1, nb_qp_u
                    vmin = minval(Wscript(i1, :))
                    vmax = maxval(Wscript(i1, :))
                    U_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_qp_v
                    vmin = minval(Wscript(:, i2))
                    vmax = maxval(Wscript(:, i2))
                    U_v(i2) = sqrt(vmin*vmax)
                end do
            end if
        end do
                                    
    end subroutine diagonal_decomposition_2d

    subroutine diagonal_decomposition_3d(nb_qp_u, nb_qp_v, nb_qp_w, CC, &
                                        U_u, U_v, U_w, W_u, W_v, W_w)
        !! Tensor decomposition to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_qp_u, nb_qp_v, nb_qp_w
        double precision, intent(in) :: CC
        dimension :: CC(3, nb_qp_u*nb_qp_v*nb_qp_w)

        double precision, intent(inout) :: U_u, U_v, U_w, W_u, W_v, W_w
        dimension ::    U_u(nb_qp_u), U_v(nb_qp_v), U_w(nb_qp_w), &
                        W_u(nb_qp_u), W_v(nb_qp_v), W_w(nb_qp_w)

        ! Local data
        ! ---------------
        double precision :: Vscript(nb_qp_u, nb_qp_v, nb_qp_w), & 
                            Wscript(2, nb_qp_u, nb_qp_v, nb_qp_w), &
                            Mscript(nb_qp_u, nb_qp_v, nb_qp_w), &
                            Nscript(nb_qp_u, nb_qp_v, nb_qp_w)

        ! Loop
        integer :: k, l, i1, i2, i3
        integer :: genpos
        double precision :: vmin, vmax
        integer :: cont
        double precision :: UU(3), WW(3), WWlk(2)

        do k = 1, 3
            ! Set Vscript
            Vscript = 0.0d0
            do i3 = 1, nb_qp_w
                do i2 = 1, nb_qp_v
                    do i1 = 1, nb_qp_u
                        genpos = i1 + (i2-1)*nb_qp_u + (i3-1)*nb_qp_u*nb_qp_v
                        UU = [U_u(i1), U_v(i2), U_w(i3)]
                        Vscript(i1, i2, i3) = CC(k, genpos)*UU(k)/(UU(1)*UU(2)*UU(3))
                    end do
                end do
            end do

            ! Update w
            if (k.eq.1) then 
                do i1 = 1, nb_qp_u
                    vmin = minval(Vscript(i1, :, :))
                    vmax = maxval(Vscript(i1, :, :))
                    W_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_qp_v
                    vmin = minval(Vscript(:, i2, :))
                    vmax = maxval(Vscript(:, i2, :))
                    W_v(i2) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do i3 = 1, nb_qp_v
                    vmin = minval(Vscript(:, :, i3))
                    vmax = maxval(Vscript(:, :, i3))
                    W_w(i3) = sqrt(vmin*vmax)
                end do
            end if
        end do

        do k = 1, 3
            ! Initialize
            Wscript = 0.0d0
            Nscript = 0.0d0
            Mscript = 0.0d0
            cont = 0
            do l = 1, 3
                if (k.ne.l) then 
                    cont = cont + 1
                    ! Set Wscript
                    do i3 = 1, nb_qp_w
                        do i2 = 1, nb_qp_v
                            do i1 = 1, nb_qp_u
                                genpos = i1 + (i2-1)*nb_qp_u + (i3-1)*nb_qp_u*nb_qp_v
                                UU = [U_u(i1), U_v(i2), U_w(i3)]
                                WW = [W_u(i1), W_v(i2), W_w(i3)]
                                Wscript(cont, i1, i2, i3) = CC(k, genpos)*UU(k)*UU(l)&
                                                            /(UU(1)*UU(2)*UU(3)*WW(k))
                            end do
                        end do
                    end do
                end if
            end do

            ! Compute Nscript and Mscript
            do i3 = 1, nb_qp_w
                do i2 = 1, nb_qp_v
                    do i1 = 1, nb_qp_u
                        WWlk = Wscript(:, i1, i2, i3)
                        Nscript(i1, i2, i3) = minval(WWlk)
                        Mscript(i1, i2, i3) = maxval(WWlk)
                    end do
                end do
            end do

            ! Update u
            if (k.eq.1) then 
                do i1 = 1, nb_qp_u
                    vmin = minval(Nscript(i1, :, :))
                    vmax = maxval(Mscript(i1, :, :))
                    U_u(i1) = sqrt(vmin*vmax)
                end do

            else if (k.eq.2) then
                do i2 = 1, nb_qp_v
                    vmin = minval(Nscript(:, i2, :))
                    vmax = maxval(Mscript(:, i2, :))
                    U_v(i2) = sqrt(vmin*vmax)
                end do

            else if (k.eq.3) then 
                do i3 = 1, nb_qp_v
                    vmin = minval(Nscript(:, :, i3))
                    vmax = maxval(Mscript(:, :, i3))
                    U_w(i3) = sqrt(vmin*vmax)
                end do
            end if
        end do

    end subroutine diagonal_decomposition_3d

    subroutine jacobien_mean_3d(nb_qp_u, nb_qp_v, nb_qp_w, Jacob, &
                                L1, L2, L3)

        implicit none
        ! Input /  output data
        ! -----------------------
        integer, intent(in) :: nb_qp_u, nb_qp_v, nb_qp_w
        double precision, intent(in) :: Jacob
        dimension :: Jacob(3, 3, nb_qp_u*nb_qp_v*nb_qp_w)

        double precision, intent(inout) :: L1, L2, L3

        ! Local data
        ! --------------
        double precision :: LNS

        ! Find mean of diagonal of jacobien
        L1 = sum(Jacob(1, 1, :))/(nb_qp_u*nb_qp_v*nb_qp_w)
        L2 = sum(Jacob(2, 2, :))/(nb_qp_u*nb_qp_v*nb_qp_w)
        L3 = sum(Jacob(3, 3, :))/(nb_qp_u*nb_qp_v*nb_qp_w)

        LNS = sqrt(L1**2 + L2**2 + L3**2)
        L1 = L1/LNS
        L2 = L2/LNS
        L3 = L3/LNS

    end subroutine jacobien_mean_3d

end module tensor_methods