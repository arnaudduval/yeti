module datastructure 

    implicit none
    type :: structure
        integer :: dimen
        integer, allocatable, dimension(:) :: nrows, ncols, nnzs
        integer, allocatable, dimension(:, :) :: indi, indj, indi_T, indj_T
        double precision, allocatable, dimension(:, :, :) :: bw, bw_T

    end type structure

contains

    subroutine init_2datastructure(datstruct, nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, nnz_v, indi_u, indj_u, indi_v, indj_v, &
                                data_B_u, data_B_v, data_W_u, data_W_v)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 2
        type(structure), allocatable :: datstruct
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

        datstruct%indi(1, 1:nr_u+1) = indi_u
        datstruct%indi(2, 1:nr_v+1) = indi_v

        datstruct%indj(1, 1:nnz_u) = indj_u
        datstruct%indj(2, 1:nnz_v) = indj_v

        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(2, 1:nnz_v, 1:2) = data_B_v

        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%bw(2, 1:nnz_v, 3:6) = data_W_v
        
    end subroutine init_2datastructure

    subroutine init_3datastructure(datstruct, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
            nnz_u, nnz_v, nnz_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
            data_B_u, data_B_v, data_B_w, data_W_u, data_W_v, data_W_w)

        implicit none 
        ! Input / output data
        ! --------------------
        integer, parameter :: dimen = 3
        type(structure), allocatable :: datstruct
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

        datstruct%indi(1, 1:nr_u+1) = indi_u
        datstruct%indi(2, 1:nr_v+1) = indi_v
        datstruct%indi(3, 1:nr_w+1) = indi_w

        datstruct%indj(1, 1:nnz_u) = indj_u
        datstruct%indj(2, 1:nnz_v) = indj_v
        datstruct%indj(3, 1:nnz_w) = indj_w

        datstruct%bw(1, 1:nnz_u, 1:2) = data_B_u
        datstruct%bw(2, 1:nnz_v, 1:2) = data_B_v
        datstruct%bw(3, 1:nnz_w, 1:2) = data_B_w

        datstruct%bw(1, 1:nnz_u, 3:6) = data_W_u
        datstruct%bw(2, 1:nnz_v, 3:6) = data_W_v
        datstruct%bw(3, 1:nnz_w, 3:6) = data_W_w
        
    end subroutine init_3datastructure

    subroutine update_datastructure(datstruct, dimen, table)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure), allocatable :: datstruct
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

        if (dimen.ne.datstruct%dimen) stop 'table not well defined.'
        nr_new_list  = datstruct%nrows 
        nnz_new_list = 0
        do i = 1, dimen
            nr  = datstruct%nrows(i)
            nnz = datstruct%nnzs(i)

            if (all(table(i, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, datstruct%nrows(i)/)
            else if (all(table(i, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(i, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = datstruct%nrows(i)
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

        do i = 1, dimen
            nr  = datstruct%nrows(i)
            nnz = datstruct%nnzs(i)

            if (all(table(i, :).eqv.(/.true., .true./))) then
                allocate(rows2er(2))
                rows2er = (/1, datstruct%nrows(i)/)
            else if (all(table(i, :).eqv.(/.true., .false./))) then
                allocate(rows2er(1))
                rows2er = 1
            else if (all(table(i, :).eqv.(/.false., .true./))) then
                allocate(rows2er(1))
                rows2er = datstruct%nrows(i)
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
            deallocate(indi_out, indj_out, bw_out)
        end do

        datstruct%nrows = nr_new_list
        datstruct%nnzs  = nnz_new_list

    end subroutine update_datastructure

    subroutine getcsr2csc(datstruct)
        implicit none 
        ! Input / output data
        ! --------------------
        type(structure), allocatable :: datstruct

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

end module datastructure