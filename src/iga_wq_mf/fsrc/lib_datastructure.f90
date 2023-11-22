module datastructure 

    implicit none
    type :: structure
        logical :: isspacetime = .false.
        integer :: dimen, nr_total, nc_total
        integer, allocatable, dimension(:) :: nrows, ncols, nnzs, dof, dod
        integer, allocatable, dimension(:, :) :: indi, indj, indi_T, indj_T
        double precision, allocatable, dimension(:) :: diageigval_sp
        double precision, allocatable, dimension(:, :) :: eigval_sp_dir
        double precision, allocatable, dimension(:, :, :) :: bw, bw_T, eigvec_sp_dir
        double precision, allocatable, dimension(:, :) :: univrightcoefs, univleftcoefs
        double complex, allocatable, dimension(:, :) :: Lschur_tm, Rschur_tm, Luptr_tm, Ruptr_tm
    end type structure

contains

    subroutine init_1datastructure(datstruct, nr_u, nc_u, nnz_u, indi_u, indj_u, data_B_u, data_W_u)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 1
        type(structure) :: datstruct
        integer, intent(in) :: nr_u, nc_u, nnz_u
        integer, intent(in) :: indi_u, indj_u
        dimension :: indi_u(nr_u+1), indj_u(nnz_u)
        double precision, intent(in) :: data_B_u, data_W_u
        dimension :: data_B_u(nnz_u, 2), data_W_u(nnz_u, 4)

        datstruct%dimen = dimen
        allocate(datstruct%nrows(dimen), datstruct%ncols(dimen), datstruct%nnzs(dimen))
        datstruct%nrows = (/nr_u/)
        datstruct%ncols = (/nc_u/)
        datstruct%nnzs  = (/nnz_u/)

        allocate(datstruct%indi(dimen, maxval(datstruct%nrows)+1), & 
                datstruct%indj(dimen, maxval(datstruct%nnzs)), &
                datstruct%bw(dimen, maxval(datstruct%nnzs), 6))

        datstruct%indi = 0
        datstruct%indi(1, 1:nr_u+1) = indi_u

        datstruct%indj = 0
        datstruct%indj(1, 1:nnz_u) = indj_u

        datstruct%bw = 0.d0
        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%nr_total = product(datstruct%nrows)
        datstruct%nc_total = product(datstruct%ncols)
        
    end subroutine init_1datastructure

    subroutine init_2datastructure(datstruct, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 2
        type(structure) :: datstruct
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

        datstruct%dimen = dimen
        allocate(datstruct%nrows(dimen), datstruct%ncols(dimen), datstruct%nnzs(dimen))
        datstruct%nrows = (/nr_u, nr_v/)
        datstruct%ncols = (/nc_u, nc_v/)
        datstruct%nnzs  = (/nnz_u, nnz_v/)

        allocate(datstruct%indi(dimen, maxval(datstruct%nrows)+1), & 
                datstruct%indj(dimen, maxval(datstruct%nnzs)), &
                datstruct%bw(dimen, maxval(datstruct%nnzs), 6))

        datstruct%indi = 0
        datstruct%indi(1, 1:nr_u+1) = indi_u
        datstruct%indi(2, 1:nr_v+1) = indi_v

        datstruct%indj = 0
        datstruct%indj(1, 1:nnz_u) = indj_u
        datstruct%indj(2, 1:nnz_v) = indj_v

        datstruct%bw = 0.d0
        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(2, 1:nnz_v, 1:2) = data_B_v
        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%bw(2, 1:nnz_v, 3:6) = data_W_v
        datstruct%nr_total = product(datstruct%nrows)
        datstruct%nc_total = product(datstruct%ncols)
        
    end subroutine init_2datastructure

    subroutine init_3datastructure(datstruct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 3
        type(structure) :: datstruct
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                        data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

        datstruct%dimen = dimen
        allocate(datstruct%nrows(dimen), datstruct%ncols(dimen), datstruct%nnzs(dimen))
        datstruct%nrows = (/nr_u, nr_v, nr_w/)
        datstruct%ncols = (/nc_u, nc_v, nc_w/)
        datstruct%nnzs  = (/nnz_u, nnz_v, nnz_w/)

        allocate(datstruct%indi(dimen, maxval(datstruct%nrows)+1), & 
                datstruct%indj(dimen, maxval(datstruct%nnzs)), &
                datstruct%bw(dimen, maxval(datstruct%nnzs), 6))

        datstruct%indi = 0
        datstruct%indi(1, 1:nr_u+1) = indi_u
        datstruct%indi(2, 1:nr_v+1) = indi_v
        datstruct%indi(3, 1:nr_w+1) = indi_w

        datstruct%indj = 0
        datstruct%indj(1, 1:nnz_u) = indj_u
        datstruct%indj(2, 1:nnz_v) = indj_v
        datstruct%indj(3, 1:nnz_w) = indj_w

        datstruct%bw = 0.d0
        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(2, 1:nnz_v, 1:2) = data_B_v
        datstruct%bw(3, 1:nnz_w, 1:2) = data_B_w
        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%bw(2, 1:nnz_v, 3:6) = data_W_v
        datstruct%bw(3, 1:nnz_w, 3:6) = data_W_w
        datstruct%nr_total = product(datstruct%nrows)
        datstruct%nc_total = product(datstruct%ncols)

    end subroutine init_3datastructure

    subroutine init_4datastructure(datstruct, nr_u, nc_u, nr_v, nr_t, nc_v, nr_w, nc_w, nc_t, &
            nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 4
        type(structure) :: datstruct
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                        data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                        data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)

        datstruct%dimen = dimen
        allocate(datstruct%nrows(dimen), datstruct%ncols(dimen), datstruct%nnzs(dimen))
        datstruct%nrows = (/nr_u, nr_v, nr_w, nr_t/)
        datstruct%ncols = (/nc_u, nc_v, nc_w, nc_t/)
        datstruct%nnzs  = (/nnz_u, nnz_v, nnz_w, nnz_t/)

        allocate(datstruct%indi(dimen, maxval(datstruct%nrows)+1), & 
                datstruct%indj(dimen, maxval(datstruct%nnzs)), &
                datstruct%bw(dimen, maxval(datstruct%nnzs), 6))

        datstruct%indi = 0
        datstruct%indi(1, 1:nr_u+1) = indi_u
        datstruct%indi(2, 1:nr_v+1) = indi_v
        datstruct%indi(3, 1:nr_w+1) = indi_w
        datstruct%indi(4, 1:nr_t+1) = indi_t

        datstruct%indj = 0
        datstruct%indj(1, 1:nnz_u) = indj_u
        datstruct%indj(2, 1:nnz_v) = indj_v
        datstruct%indj(3, 1:nnz_w) = indj_w
        datstruct%indj(4, 1:nnz_t) = indj_t

        datstruct%bw = 0.d0
        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(2, 1:nnz_v, 1:2) = data_B_v
        datstruct%bw(3, 1:nnz_w, 1:2) = data_B_w
        datstruct%bw(4, 1:nnz_t, 1:2) = data_B_t
        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%bw(2, 1:nnz_v, 3:6) = data_W_v
        datstruct%bw(3, 1:nnz_w, 3:6) = data_W_w
        datstruct%bw(4, 1:nnz_t, 3:6) = data_W_t
        datstruct%nr_total = product(datstruct%nrows)
        datstruct%nc_total = product(datstruct%ncols)

    end subroutine init_4datastructure

    subroutine setup_univariatecoefs(datstruct, nr, nc, univrightcoefs, univleftcoefs)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct
        integer, intent(in) :: nr, nc
        double precision, intent(in) :: univrightcoefs, univleftcoefs
        dimension :: univrightcoefs(nr, nc), univleftcoefs(nr, nc)

        allocate(datstruct%univrightcoefs(nr, nc), datstruct%univleftcoefs(nr, nc))
        datstruct%univrightcoefs = univrightcoefs
        datstruct%univleftcoefs  = univleftcoefs
    end subroutine setup_univariatecoefs

    subroutine get_innernodes__(datstruct, dimen, table)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct
        integer, intent(in) :: dimen
        logical, intent(in) :: table
        dimension :: table(dimen, 2)

        ! Local data
        ! ----------
        logical :: mask(datstruct%nr_total)
        integer :: ndof, ndod, c, i, j, k, l
        integer, dimension(dimen) :: inf, sup, tmp

        if (dimen.ne.datstruct%dimen) stop 'Not the same dimension'

        inf = 1; sup = datstruct%nrows
        do i = 1, dimen
            if (table(i, 1)) inf(i) = inf(i) + 1
            if (table(i, 2)) sup(i) = sup(i) - 1  
        end do

        tmp  = sup - inf + 1
        ndof = product(tmp)
        allocate(datstruct%dof(ndof)); datstruct%dof = 0
        c = 1

        if (dimen.eq.1) then
            
            do i = inf(1), sup(1)
                datstruct%dof(c) = i
                c = c + 1 
            end do

        else if (dimen.eq.2) then

            do j = inf(2), sup(2)
                do i = inf(1), sup(1)
                    datstruct%dof(c) = i + (j-1)*datstruct%nrows(1)
                    c = c + 1 
                end do
            end do

        else if (dimen.eq.3) then

            do k = inf(3), sup(3)
                do j = inf(2), sup(2)
                    do i = inf(1), sup(1)
                        datstruct%dof(c) = i + (j-1)*datstruct%nrows(1) + (k-1)*datstruct%nrows(1)*datstruct%nrows(2)
                        c = c + 1 
                    end do
                end do
            end do

        else if (dimen.eq.4) then
            
            do l = inf(4), sup(4)
                do k = inf(3), sup(3)
                    do j = inf(2), sup(2)
                        do i = inf(1), sup(1)
                            datstruct%dof(c) = i + (j-1)*datstruct%nrows(1) + (k-1)*datstruct%nrows(1)*datstruct%nrows(2) &
                                                + (l-1)*datstruct%nrows(1)*datstruct%nrows(2)*datstruct%nrows(3)
                            c = c + 1 
                        end do
                    end do
                end do
            end do

        else
            stop 'Try 2, 3 or 4 dimensions'
        end if

        mask = .true.
        mask(datstruct%dof) = .false.
        ndod = count(mask)
        allocate(datstruct%dod(ndod))
        if (ndod.ge.1) datstruct%dod = pack([(i, i = 1, size(mask))], mask)

    end subroutine get_innernodes__

    subroutine update_datastructure(datstruct, dimen, table)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct
        integer, intent(in) :: dimen
        logical, intent(in) :: table
        dimension :: table(dimen, 2)

        ! Local data
        ! ----------
        integer :: i, nr2er, nr, nnz, nr_new_list(dimen), nnz_new_list(dimen)
        integer, allocatable, dimension(:) :: rows2er
        integer, dimension(:, :), allocatable :: indi, indj
        double precision, dimension(:, :, :), allocatable :: bw
        integer, dimension(:), allocatable :: indi_out, indj_out
        double precision, dimension(:, :), allocatable :: bw_out

        if (dimen.ne.datstruct%dimen) stop 'Not the same dimension'
        call get_innernodes__(datstruct, dimen, table)

        nr_new_list  = datstruct%nrows 
        nnz_new_list = 0
        do i = 1, dimen
            nr  = datstruct%nrows(i)
            nnz = datstruct%nnzs(i)

            if (all(table(i, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, nr/)
            else if (all(table(i, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(i, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = nr
            else 
                nnz_new_list(i) = nnz
                cycle
            end if

            nr2er = size(rows2er)
            call erase_rows_csr(nr2er, rows2er, 6, nr, nnz, datstruct%indi(i, 1:nr+1), &
                            datstruct%indj(i, 1:nnz), datstruct%bw(i, 1:nnz, :), nnz_new_list(i), &
                            indi_out, indj_out, bw_out)
            deallocate(rows2er)
            nr_new_list(i) = nr - nr2er

        end do

        allocate(indi(dimen, maxval(datstruct%nrows)+1), indj(dimen, maxval(datstruct%nnzs)), &
                bw(dimen, maxval(datstruct%nnzs), 6))
        indi = datstruct%indi; indj = datstruct%indj; bw = datstruct%bw
        deallocate(datstruct%indi, datstruct%indj, datstruct%bw)
        allocate(datstruct%indi(dimen, maxval(nr_new_list)+1), &
                datstruct%indj(dimen, maxval(nnz_new_list)), &
                datstruct%bw(dimen, maxval(nnz_new_list), 6))
        datstruct%indi = 0; datstruct%indj = 0; datstruct%bw = 0.d0

        do i = 1, dimen
            nr  = datstruct%nrows(i)
            nnz = datstruct%nnzs(i)

            if (all(table(i, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, nr/)
            else if (all(table(i, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(i, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = nr
            else 
                datstruct%indi(i, 1:nr+1) = indi(i, 1:nr+1)
                datstruct%indj(i, 1:nnz)  = indj(i, 1:nnz)
                datstruct%bw(i, 1:nnz, :) = bw(i, 1:nnz, :)
                cycle
            end if
            
            nr2er = size(rows2er)
            allocate(indi_out(nr_new_list(i)+1), indj_out(nnz_new_list(i)), bw_out(nnz_new_list(i), 6))
            call erase_rows_csr(nr2er, rows2er, 6, nr, nnz, indi(i, 1:nr+1), indj(i, 1:nnz), &
                bw(i, 1:nnz, :), nnz_new_list(i), indi_out, indj_out, bw_out)
            datstruct%indi(i, 1:nr_new_list(i)+1) = indi_out
            datstruct%indj(i, 1:nnz_new_list(i))  = indj_out
            datstruct%bw(i, 1:nnz_new_list(i), :) = bw_out
            deallocate(rows2er, indi_out, indj_out, bw_out)
        end do

        datstruct%nrows = nr_new_list
        datstruct%nnzs  = nnz_new_list
        datstruct%nr_total = product(datstruct%nrows)
        datstruct%nc_total = product(datstruct%ncols)

    end subroutine update_datastructure

    subroutine getcsr2csc(datstruct)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct

        ! Local data
        ! ----------
        integer :: i, nr, nc, nnz
        integer, dimension(:), allocatable :: indi, indj, indi_T, indj_T
        double precision, dimension(:, :), allocatable :: bw, bw_T

        allocate(datstruct%indi_T(datstruct%dimen, maxval(datstruct%ncols)+1), &
                datstruct%indj_T(datstruct%dimen, maxval(datstruct%nnzs)), &
                datstruct%bw_T(datstruct%dimen, maxval(datstruct%nnzs), 6))

        do i = 1, datstruct%dimen
            nr  = datstruct%nrows(i)
            nc  = datstruct%ncols(i)
            nnz = datstruct%nnzs(i)
            allocate(indi(nr+1), indi_T(nc+1), &
                    indj(nnz), indj_T(nnz), &
                    bw(nnz, 6), bw_T(nnz, 6))
            indi = datstruct%indi(i, 1:nr+1)
            indj = datstruct%indj(i, 1:nnz)
            bw   = datstruct%bw(i, 1:nnz, :)
            call csr2csc(6, nr, nc, nnz, bw, indj, indi, bw_T, indj_T, indi_T)
            datstruct%indi_T(i, 1:nc+1) = indi_T
            datstruct%indj_T(i, 1:nnz)  = indj_T
            datstruct%bw_T(i, 1:nnz, :) = bw_T
            deallocate(indi, indj, indi_T, indj_T, bw, bw_T)
        end do

    end subroutine getcsr2csc

    subroutine space_eigendecomposition(datstruct, size_in, mean)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct
        integer :: size_in
        double precision, intent(in) :: mean
        dimension :: mean(size_in)
        
        ! Local data
        ! ----------
        integer :: i, nr, nc, nnz, ncols, nrows, dimen_sp
        integer, dimension(:), allocatable :: indi, indj
        double precision, dimension(:), allocatable :: ones, eigvalues, massdiag, stiffdiag
        double precision, dimension(:, :), allocatable :: bw, eigvectors
        
        if (datstruct%isspacetime) then 
            dimen_sp = datstruct%dimen - 1
        else
            dimen_sp = datstruct%dimen
        end if
        if (size_in.lt.dimen_sp) stop 'Size problem'  
        nrows = maxval(datstruct%nrows); ncols = maxval(datstruct%ncols)
        allocate(datstruct%eigval_sp_dir(dimen_sp, nrows), &
                datstruct%eigvec_sp_dir(dimen_sp, nrows, nrows), &
                ones(ncols))
        datstruct%eigval_sp_dir = 0.d0; datstruct%eigvec_sp_dir = 0.d0; ones = 1.d0

        ! Eigen decomposition
        do i = 1, dimen_sp
            nr = datstruct%nrows(i); nc = datstruct%ncols(i); nnz = datstruct%nnzs(i)
            allocate(indi(nr+1), indj(nnz), bw(nnz, 6), &
            eigvalues(nr), eigvectors(nr, nr), stiffdiag(nr), massdiag(nr))
            indi = datstruct%indi(i, 1:nr+1)
            indj = datstruct%indj(i, 1:nnz)
            bw   = datstruct%bw(i, 1:nnz, :)

            if (allocated(datstruct%univrightcoefs).and.allocated(datstruct%univleftcoefs)) then 
                call stiffmass_eigendecomposition(nr, nc, datstruct%univrightcoefs(i, 1:nc), &
                                        datstruct%univleftcoefs(i, 1:nc), nnz, indi, indj, &
                                        bw(:, 1:2), bw(:, 3:6), eigvalues, eigvectors, stiffdiag, massdiag)
            else
                call stiffmass_eigendecomposition(nr, nc, ones(1:nc), ones(1:nc), nnz, indi, indj, &
                                        bw(:, 1:2), bw(:, 3:6), eigvalues, eigvectors, stiffdiag, massdiag)
            end if
            datstruct%eigval_sp_dir(i, 1:nr) = eigvalues
            datstruct%eigvec_sp_dir(i, 1:nr, 1:nr) = eigvectors
            deallocate(indi, indj, bw, eigvalues, eigvectors, massdiag, stiffdiag)
        end do
        deallocate(ones)

        allocate(ones(maxval(datstruct%nrows(:dimen_sp))), &
                datstruct%diageigval_sp(product(datstruct%nrows(:dimen_sp))))
        ones = 1.d0; datstruct%diageigval_sp = 0.d0

        if (dimen_sp.eq.2) then 
            call find_parametric_diag_2d(datstruct%nrows(1), datstruct%nrows(2), ones(1:datstruct%nrows(1)), &
                                    ones(1:datstruct%nrows(2)), datstruct%eigval_sp_dir(1, 1:datstruct%nrows(1)), &
                                    datstruct%eigval_sp_dir(2, 1:datstruct%nrows(2)), mean(1:dimen_sp), &
                                    datstruct%diageigval_sp)
        else if (dimen_sp.eq.3) then
            call find_parametric_diag_3d(datstruct%nrows(1), datstruct%nrows(2), datstruct%nrows(3), &
                                    ones(1:datstruct%nrows(1)), ones(1:datstruct%nrows(2)), &
                                    ones(1:datstruct%nrows(3)), datstruct%eigval_sp_dir(1, 1:datstruct%nrows(1)), &
                                    datstruct%eigval_sp_dir(2, 1:datstruct%nrows(2)), &
                                    datstruct%eigval_sp_dir(3, 1:datstruct%nrows(3)), mean(1:dimen_sp), &
                                    datstruct%diageigval_sp)
        else
            stop 'Until not coded'
        end if

    end subroutine space_eigendecomposition

    subroutine time_schurdecomposition(datstruct)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure) :: datstruct

        ! Local data
        ! ----------
        integer :: nr, nc, nnz, ncols, dimen_tm
        integer, dimension(:), allocatable :: indi, indj
        double precision, dimension(:), allocatable :: Mdiag, Adiag
        double precision, dimension(:, :), allocatable :: bw
        double precision, dimension(:), allocatable :: univrightcoefs, univleftcoefs
        
        ncols = maxval(datstruct%ncols)
        allocate(univrightcoefs(ncols), univleftcoefs(ncols))
        univrightcoefs = 1.d0; univleftcoefs = 1.d0

        if (.not.(datstruct%isspacetime)) stop 'Only for space-time formulation'
        dimen_tm = datstruct%dimen
        nr = datstruct%nrows(dimen_tm); nc = datstruct%ncols(dimen_tm); nnz = datstruct%nnzs(dimen_tm)
        allocate(indi(nr+1), indj(nnz), bw(nnz, 6))
        indi = datstruct%indi(dimen_tm, 1:nr+1)
        indj = datstruct%indj(dimen_tm, 1:nnz)
        bw   = datstruct%bw(dimen_tm, 1:nnz, :)

        allocate(datstruct%Lschur_tm(nr, nr), datstruct%Rschur_tm(nr, nr), &
                datstruct%Luptr_tm(nr, nr), datstruct%Ruptr_tm(nr, nr), &
                Adiag(nr), Mdiag(nr))

        if (allocated(datstruct%univrightcoefs).and.allocated(datstruct%univleftcoefs)) then 
            call advmass_schurdecomposition(nr, nc, datstruct%univrightcoefs(dimen_tm, 1:nc), &
                                    datstruct%univleftcoefs(dimen_tm, 1:nc), nnz, indi, indj, &
                                    bw(:, 1:2), bw(:, 3:6), datstruct%Lschur_tm, datstruct%Rschur_tm,&
                                    datstruct%Luptr_tm, datstruct%Ruptr_tm, Adiag, Mdiag)
        else
            call advmass_schurdecomposition(nr, nc, univrightcoefs(1:nc), univleftcoefs(1:nc), nnz, indi, indj, &
                                    bw(:, 1:2), bw(:, 3:6), datstruct%Lschur_tm, datstruct%Rschur_tm,&
                                    datstruct%Luptr_tm, datstruct%Ruptr_tm, Adiag, Mdiag)
        end if

    end subroutine time_schurdecomposition

    subroutine set2zero(datstruct, nc_total, array_inout)
        implicit none 
        ! Input / output data
        ! --------------------
        integer, intent(in) :: nc_total
        type(structure) :: datstruct
        double precision, intent(inout) :: array_inout(nc_total)

        array_inout(datstruct%dod) = 0.d0

    end subroutine set2zero

end module datastructure