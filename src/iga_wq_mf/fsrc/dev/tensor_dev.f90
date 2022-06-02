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