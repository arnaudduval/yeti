subroutine tensor_n_mode_product(nc_u, nc_v, nc_w, X, nr, nc, U, n, nu, nv, nw, R)
    !! Evaluates tensor n-mode product with a matrix (R = X x_n U) (x_n: tensor n-mode product) 
    !! Based on "Tensor Decompositions and Applications" by Tamara Kolda and Brett Bader
    !! Tensor X = (nc_u, nc_v, nc_w)
    !! Matrix U = (nr, nc)
    !! Tensor R = (nu, nv, nw) (It depends on 'n')
    !! Ex: if n=1, R(nr, nc_v, nc_w) and nc=nc_u

    use omp_lib
    implicit none
    ! Input / output data 
    ! -------------------- 
    integer, intent(in) :: nc_u, nc_v, nc_w, nr, nc, n, nu, nv, nw
    double precision, intent(in) :: X, U
    dimension :: X(nc_u*nc_v*nc_w), U(nr*nc)

    double precision, intent(out) ::  R
    dimension :: R(nu*nv*nw)

    ! Local data
    ! ---------------
    integer :: genPosX, genPosU, genPosR, nb_tasks
    integer :: ju, jv, jw, i
    double precision :: sum

    ! Initialize
    R = 0.d0

    if (n.eq.1) then 
        !$OMP PARALLEL PRIVATE(sum, ju, genPosX, genPosU, genPosR)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nr/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                do i = 1, nr
                    sum = 0.d0
                    do ju = 1, nc_u
                        genPosX = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                        genPosU = i + (ju-1)*nr
                        sum = sum + X(genPosX)*U(genPosU)
                    end do
                    genPosR = i + (jv-1)*nr + (jw-1)*nr*nc_v
                    R(genPosR) = sum
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 
    else if (n.eq.2) then 
        !$OMP PARALLEL PRIVATE(sum, jv, genPosX, genPosU, genPosR)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nr*nc_u/nb_tasks) 
        do jw = 1, nc_w
            do ju = 1, nc_u
                do i = 1, nr
                    sum = 0.d0
                    do jv = 1, nc_v
                        genPosX = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                        genPosU = i + (jv-1)*nr
                        sum = sum + X(genPosX)*U(genPosU)
                    end do
                    genPosR = ju + (i-1)*nc_u + (jw-1)*nc_u*nr
                    R(genPosR) = sum
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL 
    else if (n.eq.3) then 
        !$OMP PARALLEL PRIVATE(sum, jw, genPosX, genPosU, genPosR)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nr*nc_v*nc_u/nb_tasks)
        do jv = 1, nc_v
            do ju = 1, nc_u
                do i = 1, nr
                    sum = 0.d0
                    do jw = 1, nc_w
                        genPosX = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                        genPosU = i + (jw-1)*nr
                        sum = sum + X(genPosX)*U(genPosU)
                    end do
                    genPosR = ju + (jv-1)*nc_u + (i-1)*nc_u*nc_v
                    R(genPosR) = sum
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    end if

end subroutine tensor_n_mode_product

subroutine tensor_n_mode_product_sp(nc_u, nc_v, nc_w, X, nr, nc, nnz, data_u, indi, indj, n, nu, nv, nw, R)
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
    double precision, intent(in) :: X, data_u
    dimension :: X(nc_u*nc_v*nc_w), data_u(nnz)
    integer, intent(in) :: indi, indj
    dimension :: indi(nr+1), indj(nnz)

    double precision, intent(out) ::  R
    dimension :: R(nu*nv*nw)

    ! Local data
    ! ---------------
    integer :: genPosX, genPosR
    integer :: ju, jv, jw, i, dummy, nb_tasks
    double precision, allocatable, dimension(:) :: data_nnz
    integer, allocatable, dimension(:) :: j_nnz
    double precision :: sum

    ! Initialize
    R = 0.d0
    dummy = nc

    if (n.eq.1) then 
        !$OMP PARALLEL PRIVATE(sum, ju, genPosX, genPosR, j_nnz, data_nnz)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nc_v*nr/nb_tasks) 
        do jw = 1, nc_w
            do jv = 1, nc_v
                do i = 1, nr
                    ! Define non zeros values of U
                    allocate(j_nnz(indi(i+1)-indi(i)), data_nnz(indi(i+1)-indi(i)))
                    j_nnz = indj(indi(i):indi(i+1)-1)
                    data_nnz = data_u(indi(i):indi(i+1)-1)
                    ! Get sum
                    sum = 0.d0
                    do ju = 1, size(j_nnz)
                        genPosX = j_nnz(ju) + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                        sum = sum + X(genPosX)*data_nnz(ju)
                    end do
                    genPosR = i + (jv-1)*nr + (jw-1)*nr*nc_v
                    R(genPosR) = sum
                    deallocate(j_nnz, data_nnz)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    else if (n.eq.2) then 
        !$OMP PARALLEL PRIVATE(sum, jv, genPosX, genPosR, j_nnz, data_nnz)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nc_w*nr*nc_u/nb_tasks) 
        do jw = 1, nc_w
            do ju = 1, nc_u
                do i = 1, nr
                    ! Define non zeros values of U
                    allocate(j_nnz(indi(i+1)-indi(i)), data_nnz(indi(i+1)-indi(i)))
                    j_nnz = indj(indi(i):indi(i+1)-1)
                    data_nnz = data_u(indi(i):indi(i+1)-1)
                    ! Get sum
                    sum = 0.d0
                    do jv = 1, size(j_nnz)
                        genPosX = ju + (j_nnz(jv)-1)*nc_u + (jw-1)*nc_u*nc_v
                        sum = sum + X(genPosX)*data_nnz(jv)
                    end do
                    genPosR = ju + (i-1)*nc_u + (jw-1)*nc_u*nr
                    R(genPosR) = sum
                    deallocate(j_nnz, data_nnz)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    else if (n.eq.3) then 
        !$OMP PARALLEL PRIVATE(sum, jw, genPosX, genPosR, j_nnz, data_nnz)
        nb_tasks = omp_get_num_threads()
        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nr*nc_v*nc_u/nb_tasks)
        do jv = 1, nc_v
            do ju = 1, nc_u
                do i = 1, nr
                    ! Define non zeros values of U
                    allocate(j_nnz(indi(i+1)-indi(i)), data_nnz(indi(i+1)-indi(i)))
                    j_nnz = indj(indi(i):indi(i+1)-1)
                    data_nnz = data_u(indi(i):indi(i+1)-1)
                    ! Get sum
                    sum = 0.d0
                    do jw = 1, size(j_nnz)
                        genPosX = ju + (jv-1)*nc_u + (j_nnz(jw)-1)*nc_u*nc_v
                        sum = sum + X(genPosX)*data_nnz(jw)
                    end do
                    genPosR = ju + (jv-1)*nc_u + (i-1)*nc_u*nc_v
                    R(genPosR) = sum
                    deallocate(j_nnz, data_nnz)
                end do
            end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    end if

end subroutine tensor_n_mode_product_sp

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
    result = 0.d0 
    do ju = 1, size_u
        sum = 0.d0
        do jv = 1, size_v
            posGen = ju + (jv-1)*size_u
            sum = sum + X0(posGen)*Vv(jv)
        end do
        result = result + sum*Vu(ju)
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
    result = 0.d0 
    do ju = 1, size_u
        sum2 = 0.d0
        do jv = 1, size_v
            sum1 = 0.d0
            do jw = 1, size_w
                posGen = ju + (jv-1)*size_u + (jw-1)*size_u*size_v
                sum1 = sum1 + X0(posGen)*Vw(jw)
            end do
            sum2 = sum2 + sum1*Vv(jv)
        end do
        result = result + sum2*Vu(ju)
    end do

end subroutine rankone3d_dot_vector

subroutine sumfact2d_dot_vector(nb_rows_u, nb_cols_u, &
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
    dimension :: Mu(nb_rows_u, nb_cols_u), Mv(nb_rows_v, nb_cols_v)

    double precision, intent(inout) :: vector_out
    dimension :: vector_out(nb_rows_u*nb_rows_v)

    ! Local data 
    ! -------------
    double precision :: sum
    integer :: iu, iv, nb_tasks, genPos_out

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

end subroutine sumfact2d_dot_vector

subroutine sumfact3d_dot_vector(nb_rows_u, nb_cols_u, &
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
    dimension :: Mu(nb_rows_u, nb_cols_u), Mv(nb_rows_v, nb_cols_v), Mw(nb_rows_w, nb_cols_w)

    double precision, intent(out) :: vector_out
    dimension :: vector_out(nb_rows_u*nb_rows_v*nb_rows_w)

    ! Local data 
    ! -------------
    double precision :: sum
    integer :: iu, iv, iw, nb_tasks, genPos_out

    !$OMP PARALLEL PRIVATE(genPos_out, sum) 
    nb_tasks = omp_get_num_threads()
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u*nb_rows_v*nb_rows_w/nb_tasks) 
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

end subroutine sumfact3d_dot_vector

subroutine sumfact2d_dot_vector_sp(nb_rows_u, nb_cols_u, &
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
    !! Indices must be in CSR format

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
    !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, nb_rows_u*nb_rows_v/nb_tasks) 
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
            data_nnz_u = 0.d0
            offset = indi_u(iu)
            do ju = 1, nnz_u
                indj_nnz_u(ju) = indj_u(ju+offset-1)
                data_nnz_u(ju) = data_u(ju+offset-1)
            end do

            ! Set values of row in Mv
            allocate(indj_nnz_v(nnz_v), data_nnz_v(nnz_v))
            indj_nnz_v = 0
            data_nnz_v = 0.d0
            offset = indi_v(iv)
            do jv = 1, nnz_v
                indj_nnz_v(jv) = indj_v(jv+offset-1)
                data_nnz_v(jv) = data_v(jv+offset-1)
            end do

            allocate(tensor(nnz_u*nnz_v))
            tensor = 0.d0
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

end subroutine sumfact2d_dot_vector_sp

subroutine sumfact3d_dot_vector_sp(nb_rows_u, nb_cols_u, &
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
    !! Indices must be in CSR format

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
    !$OMP DO COLLAPSE(3) SCHEDULE(STATIC, nb_rows_u*nb_rows_v*nb_rows_w/nb_tasks) 
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
                tensor = 0.d0
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

end subroutine sumfact3d_dot_vector_sp