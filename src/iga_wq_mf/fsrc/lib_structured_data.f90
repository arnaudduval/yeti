module structured_data 

    implicit none
    type :: basis_data
        integer :: dimen, nr_total, nc_total
        integer, allocatable, dimension(:) :: nrows, ncols, nnzs
        integer, allocatable, dimension(:, :) :: indi, indj, indiT, indjT
        double precision, allocatable, dimension(:, :, :) :: data_bw, data_bwT
        double precision, allocatable, dimension(:, :, :, :) :: BTdense, Wdense
    end type basis_data

    type :: reduced_system
        logical :: isspacetime = .false., isalreadyreduced=.false.
        type(basis_data) :: basisdata
        integer, allocatable, dimension(:) :: dof, dod
        double precision, allocatable, dimension(:) :: meancoefs
        double precision, allocatable, dimension(:, :) :: univrightcoefs, univleftcoefs
        double precision, allocatable, dimension(:, :) :: eigval_sp_dir
        double precision, allocatable, dimension(:, :, :) :: eigvec_sp_dir
        double complex, allocatable, dimension(:, :) :: Lschur_tm, Rschur_tm, Luptr_tm, Ruptr_tm
        double precision, allocatable, dimension(:) :: diageigval_sp
    end type reduced_system

contains

    subroutine init_1basisdata(basisdata, nr_u, nc_u, nnz_u, indi_u, indj_u, data_B_u, data_W_u)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 1
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_u, nc_u, nnz_u
        integer, intent(in) :: indi_u, indj_u
        dimension :: indi_u(nr_u+1), indj_u(nnz_u)
        double precision, intent(in) :: data_B_u, data_W_u
        dimension :: data_B_u(nnz_u, 2), data_W_u(nnz_u, 4)

        basisdata%dimen = dimen
        allocate(basisdata%nrows(dimen), basisdata%ncols(dimen), basisdata%nnzs(dimen))
        basisdata%nrows = (/nr_u/)
        basisdata%ncols = (/nc_u/)
        basisdata%nnzs  = (/nnz_u/)

        allocate(basisdata%indi(dimen, maxval(basisdata%nrows)+1), & 
                basisdata%indj(dimen, maxval(basisdata%nnzs)), &
                basisdata%data_bw(dimen, maxval(basisdata%nnzs), 6))

        basisdata%indi = 0
        basisdata%indi(1, 1:nr_u+1) = indi_u

        basisdata%indj = 0
        basisdata%indj(1, 1:nnz_u) = indj_u

        basisdata%data_bw = 0.d0
        basisdata%data_bw(1, 1:nnz_u, 1:2) = data_B_u
        basisdata%data_bw(1, 1:nnz_u, 3:6) = data_W_u
        basisdata%nr_total = product(basisdata%nrows)
        basisdata%nc_total = product(basisdata%ncols)
        
    end subroutine init_1basisdata

    subroutine init_2basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 2
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) :: indi_u, indi_v, indj_u, indj_v
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), &
                        indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4)

        basisdata%dimen = dimen
        allocate(basisdata%nrows(dimen), basisdata%ncols(dimen), basisdata%nnzs(dimen))
        basisdata%nrows = (/nr_u, nr_v/)
        basisdata%ncols = (/nc_u, nc_v/)
        basisdata%nnzs  = (/nnz_u, nnz_v/)

        allocate(basisdata%indi(dimen, maxval(basisdata%nrows)+1), & 
                basisdata%indj(dimen, maxval(basisdata%nnzs)), &
                basisdata%data_bw(dimen, maxval(basisdata%nnzs), 6))

        basisdata%indi = 0
        basisdata%indi(1, 1:nr_u+1) = indi_u
        basisdata%indi(2, 1:nr_v+1) = indi_v

        basisdata%indj = 0
        basisdata%indj(1, 1:nnz_u) = indj_u
        basisdata%indj(2, 1:nnz_v) = indj_v

        basisdata%data_bw = 0.d0
        basisdata%data_bw(1, 1:nnz_u, 1:2) = data_B_u
        basisdata%data_bw(2, 1:nnz_v, 1:2) = data_B_v
        basisdata%data_bw(1, 1:nnz_u, 3:6) = data_W_u
        basisdata%data_bw(2, 1:nnz_v, 3:6) = data_W_v
        basisdata%nr_total = product(basisdata%nrows)
        basisdata%nc_total = product(basisdata%ncols)
        
    end subroutine init_2basisdata

    subroutine init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 3
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) :: indi_u, indi_v, indi_w, indj_u, indj_v, indj_w
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                        data_B_w(nnz_w, 2), data_W_w(nnz_w, 4)

        basisdata%dimen = dimen
        allocate(basisdata%nrows(dimen), basisdata%ncols(dimen), basisdata%nnzs(dimen))
        basisdata%nrows = (/nr_u, nr_v, nr_w/)
        basisdata%ncols = (/nc_u, nc_v, nc_w/)
        basisdata%nnzs  = (/nnz_u, nnz_v, nnz_w/)

        allocate(basisdata%indi(dimen, maxval(basisdata%nrows)+1), & 
                basisdata%indj(dimen, maxval(basisdata%nnzs)), &
                basisdata%data_bw(dimen, maxval(basisdata%nnzs), 6))

        basisdata%indi = 0
        basisdata%indi(1, 1:nr_u+1) = indi_u
        basisdata%indi(2, 1:nr_v+1) = indi_v
        basisdata%indi(3, 1:nr_w+1) = indi_w

        basisdata%indj = 0
        basisdata%indj(1, 1:nnz_u) = indj_u
        basisdata%indj(2, 1:nnz_v) = indj_v
        basisdata%indj(3, 1:nnz_w) = indj_w

        basisdata%data_bw = 0.d0
        basisdata%data_bw(1, 1:nnz_u, 1:2) = data_B_u
        basisdata%data_bw(2, 1:nnz_v, 1:2) = data_B_v
        basisdata%data_bw(3, 1:nnz_w, 1:2) = data_B_w
        basisdata%data_bw(1, 1:nnz_u, 3:6) = data_W_u
        basisdata%data_bw(2, 1:nnz_v, 3:6) = data_W_v
        basisdata%data_bw(3, 1:nnz_w, 3:6) = data_W_w
        basisdata%nr_total = product(basisdata%nrows)
        basisdata%nc_total = product(basisdata%ncols)

    end subroutine init_3basisdata

    subroutine init_4basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
            nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, data_W_t)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 4
        type(basis_data) :: basisdata
        integer, intent(in) :: nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t
        integer, intent(in) :: indi_u, indi_v, indi_w, indi_t, indj_u, indj_v, indj_w, indj_t
        dimension ::    indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1), &
                        indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
        double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
        dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                        data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                        data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                        data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)

        basisdata%dimen = dimen
        allocate(basisdata%nrows(dimen), basisdata%ncols(dimen), basisdata%nnzs(dimen))
        basisdata%nrows = (/nr_u, nr_v, nr_w, nr_t/)
        basisdata%ncols = (/nc_u, nc_v, nc_w, nc_t/)
        basisdata%nnzs  = (/nnz_u, nnz_v, nnz_w, nnz_t/)

        allocate(basisdata%indi(dimen, maxval(basisdata%nrows)+1), & 
                basisdata%indj(dimen, maxval(basisdata%nnzs)), &
                basisdata%data_bw(dimen, maxval(basisdata%nnzs), 6))

        basisdata%indi = 0
        basisdata%indi(1, 1:nr_u+1) = indi_u
        basisdata%indi(2, 1:nr_v+1) = indi_v
        basisdata%indi(3, 1:nr_w+1) = indi_w
        basisdata%indi(4, 1:nr_t+1) = indi_t

        basisdata%indj = 0
        basisdata%indj(1, 1:nnz_u) = indj_u
        basisdata%indj(2, 1:nnz_v) = indj_v
        basisdata%indj(3, 1:nnz_w) = indj_w
        basisdata%indj(4, 1:nnz_t) = indj_t

        basisdata%data_bw = 0.d0
        basisdata%data_bw(1, 1:nnz_u, 1:2) = data_B_u
        basisdata%data_bw(2, 1:nnz_v, 1:2) = data_B_v
        basisdata%data_bw(3, 1:nnz_w, 1:2) = data_B_w
        basisdata%data_bw(4, 1:nnz_t, 1:2) = data_B_t
        basisdata%data_bw(1, 1:nnz_u, 3:6) = data_W_u
        basisdata%data_bw(2, 1:nnz_v, 3:6) = data_W_v
        basisdata%data_bw(3, 1:nnz_w, 3:6) = data_W_w
        basisdata%data_bw(4, 1:nnz_t, 3:6) = data_W_t
        basisdata%nr_total = product(basisdata%nrows)
        basisdata%nc_total = product(basisdata%ncols)

    end subroutine init_4basisdata

    subroutine clearbasisdata(basisdata)
        implicit none
        ! Input / output data
        ! --------------------
        type(basis_data) :: basisdata
        basisdata%dimen = 1
        basisdata%nr_total = 1
        basisdata%nc_total = 1
        if (allocated(basisdata%nrows)) deallocate(basisdata%nrows)
        if (allocated(basisdata%ncols)) deallocate(basisdata%ncols)
        if (allocated(basisdata%nnzs)) deallocate(basisdata%nnzs)
        if (allocated(basisdata%indi)) deallocate(basisdata%indi)
        if (allocated(basisdata%indj)) deallocate(basisdata%indj)
        if (allocated(basisdata%data_bw)) deallocate(basisdata%data_bw)
        if (allocated(basisdata%indiT)) deallocate(basisdata%indiT)
        if (allocated(basisdata%indjT)) deallocate(basisdata%indjT)
        if (allocated(basisdata%data_bwT)) deallocate(basisdata%data_bwT)
        if (allocated(basisdata%Wdense)) deallocate(basisdata%Wdense)
        if (allocated(basisdata%BTdense)) deallocate(basisdata%BTdense)
    end subroutine clearbasisdata

    subroutine copybasisdata(basisdata_in, basisdata_inout)
        implicit none
        ! Input / output data
        ! --------------------
        type(basis_data), intent(in) :: basisdata_in
        type(basis_data), intent(inout) :: basisdata_inout
        ! Local data
        ! ----------
        integer :: dimen, size2, size3
        dimen = basisdata_in%dimen
        call clearbasisdata(basisdata_inout)
        basisdata_inout%dimen = basisdata_in%dimen
        basisdata_inout%nr_total = basisdata_in%nr_total
        basisdata_inout%nc_total = basisdata_in%nc_total
        allocate(basisdata_inout%nrows(dimen)); basisdata_inout%nrows=basisdata_in%nrows
        allocate(basisdata_inout%ncols(dimen)); basisdata_inout%ncols=basisdata_in%ncols
        allocate(basisdata_inout%nnzs(dimen)); basisdata_inout%nnzs=basisdata_in%nnzs

        allocate(basisdata_inout%indi(dimen, maxval(basisdata_in%nrows)+1)) 
        basisdata_inout%indi=basisdata_in%indi
        allocate(basisdata_inout%indj(dimen, maxval(basisdata_in%nnzs))) 
        basisdata_inout%indj=basisdata_in%indj
        allocate(basisdata_inout%data_bw(dimen, maxval(basisdata_in%nnzs), 6)) 
        basisdata_inout%data_bw=basisdata_in%data_bw

        if (allocated(basisdata_in%indiT)) then
            allocate(basisdata_inout%indiT(dimen, maxval(basisdata_in%ncols)+1)) 
            basisdata_inout%indiT=basisdata_in%indiT
        end if
        if (allocated(basisdata_in%indjT)) then
            allocate(basisdata_inout%indjT(dimen, maxval(basisdata_in%nnzs)))
            basisdata_inout%indjT=basisdata_in%indjT
        end if
        if (allocated(basisdata_in%data_bwT)) then
            allocate(basisdata_inout%data_bwT(dimen, maxval(basisdata_in%nnzs), 6)) 
            basisdata_inout%data_bwT=basisdata_in%data_bwT
        end if

        if (allocated(basisdata_in%Wdense)) then
            size2 = size(basisdata_in%Wdense, dim=2)
            size3 = size(basisdata_in%Wdense, dim=3)
            allocate(basisdata_inout%Wdense(dimen, size2, size3, 4))
            basisdata_inout%Wdense=basisdata_in%Wdense
        end if
        if (allocated(basisdata_in%BTdense)) then
            size2 = size(basisdata_in%BTdense, dim=2)
            size3 = size(basisdata_in%BTdense, dim=3)
            allocate(basisdata_inout%BTdense(dimen, size2, size3, 2)) 
            basisdata_inout%BTdense=basisdata_in%BTdense
        end if
    end subroutine copybasisdata

    subroutine getcsrc2dense(basisdata)
        implicit none 
        ! Input / output data
        ! --------------------
        type(basis_data) :: basisdata

        ! Local data
        ! ----------
        integer :: it, j, nr, nc, nnz
        integer, dimension(:), allocatable :: indi, indj
        double precision, dimension(:, :), allocatable :: data_basis
        double precision, dimension(:, :, :), allocatable :: basisdense

        if (allocated(basisdata%BTdense)) deallocate(basisdata%BTdense)
        if (allocated(basisdata%Wdense)) deallocate(basisdata%Wdense)
        allocate(basisdata%BTdense(basisdata%dimen, maxval(basisdata%ncols), maxval(basisdata%nrows), 2), &
                basisdata%Wdense(basisdata%dimen, maxval(basisdata%nrows), maxval(basisdata%ncols), 4))
        basisdata%BTdense = 0.d0; basisdata%Wdense = 0.d0
        do it = 1, basisdata%dimen
            nr  = basisdata%nrows(it)
            nc  = basisdata%ncols(it)
            nnz = basisdata%nnzs(it)
            allocate(indi(nr+1), indj(nnz), data_basis(nnz, 6), basisdense(nr, nc, 6))
            indi = basisdata%indi(it, 1:nr+1)
            indj = basisdata%indj(it, 1:nnz)
            data_basis = basisdata%data_bw(it, 1:nnz, :)
            call multicsr2dense(6, nnz, indi, indj, data_basis, nr, nc, basisdense)
            do j = 1, 2
                basisdata%BTdense(it, 1:nc, 1:nr, j) = transpose(basisdense(:, :, j))
            end do
            do j = 1, 4
                basisdata%Wdense(it, 1:nr, 1:nc, j) = basisdense(:, :, 2+j)
            end do
            deallocate(indi, indj, data_basis, basisdense)
        end do
    end subroutine getcsrc2dense

    subroutine getcsr2csc(basisdata)
        implicit none 
        ! Input / output data
        ! --------------------
        type(basis_data) :: basisdata

        ! Local data
        ! ----------
        integer :: i, nr, nc, nnz
        integer, dimension(:), allocatable :: indi, indj, indiT, indjT
        double precision, dimension(:, :), allocatable :: bw, bwT

        if (allocated(basisdata%indiT)) deallocate(basisdata%indiT)
        if (allocated(basisdata%indjT)) deallocate(basisdata%indjT)
        if (allocated(basisdata%data_bwT)) deallocate(basisdata%data_bwT)
        allocate(basisdata%indiT(basisdata%dimen, maxval(basisdata%ncols)+1), &
                basisdata%indjT(basisdata%dimen, maxval(basisdata%nnzs)), &
                basisdata%data_bwT(basisdata%dimen, maxval(basisdata%nnzs), 6))
        basisdata%indiT = 0.d0; basisdata%indjT = 0.d0; basisdata%data_bwT = 0.d0 
        do i = 1, basisdata%dimen
            nr  = basisdata%nrows(i)
            nc  = basisdata%ncols(i)
            nnz = basisdata%nnzs(i)
            allocate(indi(nr+1), indiT(nc+1), &
                    indj(nnz), indjT(nnz), &
                    bw(nnz, 6), bwT(nnz, 6))
            indi = basisdata%indi(i, 1:nr+1)
            indj = basisdata%indj(i, 1:nnz)
            bw   = basisdata%data_bw(i, 1:nnz, :)
            call multicsr2csc(6, nr, nc, nnz, bw, indj, indi, bwT, indjT, indiT)
            basisdata%indiT(i, 1:nc+1) = indiT
            basisdata%indjT(i, 1:nnz)  = indjT
            basisdata%data_bwT(i, 1:nnz, :) = bwT
            deallocate(indi, indj, indiT, indjT, bw, bwT)
        end do

    end subroutine getcsr2csc

    subroutine get_innernodes__(basisdata, dimen, table, ndof, ndod, dof, dod)
        implicit none 
        ! Input / output data
        ! --------------------
        type(basis_data) :: basisdata
        integer, intent(in) :: dimen
        logical, intent(in) :: table
        dimension :: table(dimen, 2)
        integer, intent(inout) :: ndof, ndod
        integer, dimension(*), intent(out) :: dof, dod

        ! Local data
        ! ----------
        logical :: mask(basisdata%nr_total)
        integer :: c, i, j, k, l
        integer, dimension(dimen) :: inf, sup, tmp

        if (dimen.ne.basisdata%dimen) stop 'Not the same dimension'
        inf = 1; sup = basisdata%nrows
        do i = 1, dimen
            if (table(i, 1)) inf(i) = inf(i) + 1
            if (table(i, 2)) sup(i) = sup(i) - 1  
        end do
        
        if ((ndof.lt.0).or.(ndod.lt.0)) then
            tmp  = sup - inf + 1
            ndof = product(tmp)
            ndod = basisdata%nr_total - ndof
            return
        end if

        c = 1
        if (dimen.eq.1) then
            
            do i = inf(1), sup(1)
                dof(c) = i
                c = c + 1 
            end do

        else if (dimen.eq.2) then

            do j = inf(2), sup(2)
                do i = inf(1), sup(1)
                    dof(c) = i + (j-1)*basisdata%nrows(1)
                    c = c + 1 
                end do
            end do

        else if (dimen.eq.3) then

            do k = inf(3), sup(3)
                do j = inf(2), sup(2)
                    do i = inf(1), sup(1)
                        dof(c) = i + (j-1)*basisdata%nrows(1) + (k-1)*basisdata%nrows(1)*basisdata%nrows(2)
                        c = c + 1 
                    end do
                end do
            end do

        else if (dimen.eq.4) then
            
            do l = inf(4), sup(4)
                do k = inf(3), sup(3)
                    do j = inf(2), sup(2)
                        do i = inf(1), sup(1)
                            dof(c) = i + (j-1)*basisdata%nrows(1) + (k-1)*basisdata%nrows(1)*basisdata%nrows(2) &
                                                + (l-1)*basisdata%nrows(1)*basisdata%nrows(2)*basisdata%nrows(3)
                            c = c + 1 
                        end do
                    end do
                end do
            end do

        else
            stop 'Try 1, 2, 3 or 4 dimensions'
        end if

        mask = .true.
        mask(dof(1:ndof)) = .false.
        if (ndod.ge.1) dod(1:ndod) = pack([(i, i = 1, size(mask))], mask)

    end subroutine get_innernodes__

    !! -------------

    subroutine update_reducedsystem(redsyst, dimen, table)
        implicit none 
        ! Input / output data
        ! --------------------
        type(reduced_system) :: redsyst
        integer, intent(in) :: dimen
        logical, intent(in) :: table
        dimension :: table(dimen, 2)

        ! Local data
        ! ----------
        integer :: ndof, ndod, it, row2er_it, nr_it, nnz_it
        integer, dimension(dimen) :: nr_new_list, nnz_new_list
        integer, allocatable, dimension(:) :: dof, dod, rows2er 
        integer, dimension(:, :), allocatable :: indi_in, indj_in
        double precision, dimension(:, :, :), allocatable :: basis_in
        integer, dimension(:), allocatable :: indi_out, indj_out
        double precision, dimension(:, :), allocatable :: basis_out

        if (dimen.ne.redsyst%basisdata%dimen) stop 'Not the same dimension'
        if (redsyst%isalreadyreduced) stop 'Data basis has already been reduced'

        ndof = -1; ndod = -1
        call get_innernodes__(redsyst%basisdata, dimen, table, ndof, ndod, dof, dod)
        allocate(dof(ndof), dod(ndod))
        call get_innernodes__(redsyst%basisdata, dimen, table, ndof, ndod, dof, dod)
        allocate(redsyst%dof(ndof), redsyst%dod(ndod))
        if (ndof.gt.0) redsyst%dof = dof; deallocate(dof)
        if (ndod.gt.0) redsyst%dod = dod; deallocate(dod)

        nr_new_list  = redsyst%basisdata%nrows 
        nnz_new_list = 0
        do it = 1, dimen
            nr_it  = redsyst%basisdata%nrows(it)
            nnz_it = redsyst%basisdata%nnzs(it)

            if (all(table(it, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, nr_it/)
            else if (all(table(it, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(it, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = nr_it
            else 
                nnz_new_list(it) = nnz_it
                cycle
            end if

            row2er_it = size(rows2er)
            call erase_rows_csr(row2er_it, rows2er, 6, nr_it, nnz_it, redsyst%basisdata%indi(it, 1:nr_it+1), &
                            redsyst%basisdata%indj(it, 1:nnz_it), redsyst%basisdata%data_bw(it, 1:nnz_it, :), &
                            nnz_new_list(it), indi_out, indj_out, basis_out)
            deallocate(rows2er)
            nr_new_list(it) = nr_it - row2er_it
        end do

        allocate(indi_in(dimen, maxval(redsyst%basisdata%nrows)+1), &
                indj_in(dimen, maxval(redsyst%basisdata%nnzs)), &
                basis_in(dimen, maxval(redsyst%basisdata%nnzs), 6))
        indi_in = redsyst%basisdata%indi; indj_in = redsyst%basisdata%indj; basis_in = redsyst%basisdata%data_bw
        deallocate(redsyst%basisdata%indi, redsyst%basisdata%indj, redsyst%basisdata%data_bw)
        allocate(redsyst%basisdata%indi(dimen, maxval(nr_new_list)+1), &
                redsyst%basisdata%indj(dimen, maxval(nnz_new_list)), &
                redsyst%basisdata%data_bw(dimen, maxval(nnz_new_list), 6))
        redsyst%basisdata%indi = 0; redsyst%basisdata%indj = 0; redsyst%basisdata%data_bw = 0.d0

        do it = 1, dimen
            nr_it  = redsyst%basisdata%nrows(it)
            nnz_it = redsyst%basisdata%nnzs(it)

            if (all(table(it, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, nr_it/)
            else if (all(table(it, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(it, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = nr_it
            else 
                redsyst%basisdata%indi(it, 1:nr_it+1) = indi_in(it, 1:nr_it+1)
                redsyst%basisdata%indj(it, 1:nnz_it)  = indj_in(it, 1:nnz_it)
                redsyst%basisdata%data_bw(it, 1:nnz_it, :) = basis_in(it, 1:nnz_it, :)
                cycle
            end if

            row2er_it = size(rows2er)
            allocate(indi_out(nr_new_list(it)+1), indj_out(nnz_new_list(it)), basis_out(nnz_new_list(it), 6))
            call erase_rows_csr(row2er_it, rows2er, 6, nr_it, nnz_it, indi_in(it, 1:nr_it+1), indj_in(it, 1:nnz_it), &
                basis_in(it, 1:nnz_it, :), nnz_new_list(it), indi_out, indj_out, basis_out)
            redsyst%basisdata%indi(it, 1:nr_new_list(it)+1) = indi_out
            redsyst%basisdata%indj(it, 1:nnz_new_list(it))  = indj_out
            redsyst%basisdata%data_bw(it, 1:nnz_new_list(it), :) = basis_out
            deallocate(rows2er, indi_out, indj_out, basis_out)
        end do

        redsyst%basisdata%nrows = nr_new_list
        redsyst%basisdata%nnzs  = nnz_new_list
        redsyst%basisdata%nr_total = product(nr_new_list)
        redsyst%isalreadyreduced = .true.

    end subroutine update_reducedsystem

    subroutine setup_univariatecoefs(redsyst, nr, nc, univrightcoefs, univleftcoefs)
        implicit none 
        ! Input / output data
        ! --------------------
        type(reduced_system) :: redsyst
        integer, intent(in) :: nr, nc
        double precision, intent(in) :: univrightcoefs, univleftcoefs
        dimension :: univrightcoefs(nr, nc), univleftcoefs(nr, nc)
        if (nr.lt.redsyst%basisdata%dimen-1) stop 'Size problem'
        if (allocated(redsyst%univleftcoefs)) deallocate(redsyst%univleftcoefs)
        if (allocated(redsyst%univrightcoefs)) deallocate(redsyst%univrightcoefs)
        allocate(redsyst%univrightcoefs(nr, nc), redsyst%univleftcoefs(nr, nc))
        redsyst%univrightcoefs = univrightcoefs
        redsyst%univleftcoefs  = univleftcoefs
    end subroutine setup_univariatecoefs

    subroutine setup_meancoefs(redsyst, size_array, array)
        implicit none 
        ! Input / output data
        ! --------------------
        type(reduced_system) :: redsyst
        integer, intent(in) :: size_array
        double precision, intent(in) :: array
        dimension :: array(size_array)
        if (size_array.lt.redsyst%basisdata%dimen-1) stop 'Size problem'
        if (allocated(redsyst%meancoefs)) deallocate(redsyst%meancoefs)
        allocate(redsyst%meancoefs(size_array))
        redsyst%meancoefs = array
    end subroutine setup_meancoefs

    subroutine space_eigendecomposition(redsyst)
        implicit none 
        ! Input / output data
        ! --------------------
        type(reduced_system) :: redsyst
        
        ! Local data
        ! ----------
        integer :: i, nr, nc, nnz, ncols, nrows, dimen_sp
        integer, dimension(:), allocatable :: indi, indj
        double precision, dimension(:), allocatable :: ones, means, eigvalues, massdiag, stiffdiag
        double precision, dimension(:, :), allocatable :: basis, eigvectors
        
        dimen_sp = redsyst%basisdata%dimen
        if (redsyst%isspacetime) dimen_sp = redsyst%basisdata%dimen - 1

        nrows = maxval(redsyst%basisdata%nrows); ncols = maxval(redsyst%basisdata%ncols)
        allocate(redsyst%eigval_sp_dir(dimen_sp, nrows), &
                redsyst%eigvec_sp_dir(dimen_sp, nrows, nrows), ones(ncols))
        redsyst%eigval_sp_dir = 0.d0; redsyst%eigvec_sp_dir = 0.d0; ones = 1.d0

        ! Eigen decomposition
        do i = 1, dimen_sp
            nr = redsyst%basisdata%nrows(i); nc = redsyst%basisdata%ncols(i); nnz = redsyst%basisdata%nnzs(i)
            allocate(indi(nr+1), indj(nnz), basis(nnz, 6), &
            eigvalues(nr), eigvectors(nr, nr), stiffdiag(nr), massdiag(nr))
            indi = redsyst%basisdata%indi(i, 1:nr+1)
            indj = redsyst%basisdata%indj(i, 1:nnz)
            basis = redsyst%basisdata%data_bw(i, 1:nnz, :)

            if (allocated(redsyst%univrightcoefs).and.allocated(redsyst%univleftcoefs)) then 
                call stiffmass_eigendecomposition(nr, nc, redsyst%univrightcoefs(i, 1:nc), &
                                        redsyst%univleftcoefs(i, 1:nc), nnz, indi, indj, &
                                        basis(:, 1:2), basis(:, 3:6), eigvalues, eigvectors, stiffdiag, massdiag)
            else
                call stiffmass_eigendecomposition(nr, nc, ones(1:nc), ones(1:nc), nnz, indi, indj, &
                                        basis(:, 1:2), basis(:, 3:6), eigvalues, eigvectors, stiffdiag, massdiag)
            end if
            redsyst%eigval_sp_dir(i, 1:nr) = eigvalues
            redsyst%eigvec_sp_dir(i, 1:nr, 1:nr) = eigvectors
            deallocate(indi, indj, basis, eigvalues, eigvectors, massdiag, stiffdiag)
        end do
        deallocate(ones)

        allocate(ones(maxval(redsyst%basisdata%nrows(:dimen_sp))), &
                redsyst%diageigval_sp(product(redsyst%basisdata%nrows(:dimen_sp))))
        ones = 1.d0; redsyst%diageigval_sp = 0.d0

        allocate(means(dimen_sp)); means = 1.d0
        if (allocated(redsyst%meancoefs)) means = redsyst%meancoefs(1:dimen_sp)
        if (dimen_sp.eq.2) then 
            call find_parametric_diag_2d(redsyst%basisdata%nrows(1), redsyst%basisdata%nrows(2), &
                                    ones(1:redsyst%basisdata%nrows(1)), ones(1:redsyst%basisdata%nrows(2)), &
                                    redsyst%eigval_sp_dir(1, 1:redsyst%basisdata%nrows(1)), &
                                    redsyst%eigval_sp_dir(2, 1:redsyst%basisdata%nrows(2)), &
                                    means(1:dimen_sp), redsyst%diageigval_sp)
        else if (dimen_sp.eq.3) then
            call find_parametric_diag_3d(redsyst%basisdata%nrows(1), redsyst%basisdata%nrows(2), &
                                    redsyst%basisdata%nrows(3), ones(1:redsyst%basisdata%nrows(1)), &
                                    ones(1:redsyst%basisdata%nrows(2)), ones(1:redsyst%basisdata%nrows(3)), &
                                    redsyst%eigval_sp_dir(1, 1:redsyst%basisdata%nrows(1)), &
                                    redsyst%eigval_sp_dir(2, 1:redsyst%basisdata%nrows(2)), &
                                    redsyst%eigval_sp_dir(3, 1:redsyst%basisdata%nrows(3)), &
                                    means(1:dimen_sp), redsyst%diageigval_sp)
        else
            stop 'Try 2 or 3 dimensions'
        end if

    end subroutine space_eigendecomposition

    subroutine time_schurdecomposition(redsyst)
        implicit none 
        ! Input / output data
        ! --------------------
        type(reduced_system) :: redsyst

        ! Local data
        ! ----------
        integer :: nr, nc, nnz, pos_tm
        integer, dimension(:), allocatable :: indi, indj
        double precision, dimension(:), allocatable :: Mdiag, Adiag
        double precision, dimension(:, :), allocatable :: basis
        double precision, dimension(:), allocatable :: univrightcoefs, univleftcoefs
        
        if (.not.(redsyst%isspacetime)) stop 'Only for space-time formulation'
        pos_tm = redsyst%basisdata%dimen
        nr = redsyst%basisdata%nrows(pos_tm); nc = redsyst%basisdata%ncols(pos_tm); nnz = redsyst%basisdata%nnzs(pos_tm)
        
        allocate(univrightcoefs(nc), univleftcoefs(nc), indi(nr+1), indj(nnz), basis(nnz, 6))
        univrightcoefs = 1.d0; univleftcoefs = 1.d0
        indi = redsyst%basisdata%indi(pos_tm, 1:nr+1)
        indj = redsyst%basisdata%indj(pos_tm, 1:nnz)
        basis = redsyst%basisdata%data_bw(pos_tm, 1:nnz, :)

        allocate(redsyst%Lschur_tm(nr, nr), redsyst%Rschur_tm(nr, nr), &
                redsyst%Luptr_tm(nr, nr), redsyst%Ruptr_tm(nr, nr), &
                Adiag(nr), Mdiag(nr))

        if (allocated(redsyst%univrightcoefs).and.allocated(redsyst%univleftcoefs)) then 
            if (size(redsyst%univrightcoefs, dim=1).lt.pos_tm) stop 'Size problem'
            if (size(redsyst%univleftcoefs, dim=1).lt.pos_tm) stop 'Size problem'
            call advmass_schurdecomposition(nr, nc, redsyst%univrightcoefs(pos_tm, 1:nc), &
                                    redsyst%univleftcoefs(pos_tm, 1:nc), nnz, indi, indj, &
                                    basis(:, 1:2), basis(:, 3:6), redsyst%Lschur_tm, redsyst%Rschur_tm,&
                                    redsyst%Luptr_tm, redsyst%Ruptr_tm, Adiag, Mdiag)
        else
            call advmass_schurdecomposition(nr, nc, univrightcoefs(1:nc), univleftcoefs(1:nc), nnz, indi, indj, &
                                    basis(:, 1:2), basis(:, 3:6), redsyst%Lschur_tm, redsyst%Rschur_tm,&
                                    redsyst%Luptr_tm, redsyst%Ruptr_tm, Adiag, Mdiag)
        end if

    end subroutine time_schurdecomposition

    subroutine set2zero(redsyst, size_inout, array_inout)
        implicit none 
        ! Input / output data
        ! --------------------
        integer, intent(in) :: size_inout
        type(reduced_system) :: redsyst
        double precision, intent(inout) :: array_inout(size_inout)
        array_inout(redsyst%dod) = 0.d0
    end subroutine set2zero

end module structured_data