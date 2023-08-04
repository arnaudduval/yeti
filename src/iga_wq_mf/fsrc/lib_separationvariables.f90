module separatevariables

    implicit none

    type sepoperator
        integer :: dimen, ncols, maxit = 2
        integer, dimension(:), allocatable :: nclist
        logical, dimension(:), allocatable :: update
        double precision, dimension(:, :), allocatable :: Mcoefs, Kcoefs

    end type sepoperator

contains

    subroutine initialize_operator(obj, dimen, nclist, update)
        
        implicit none
        ! Input /  output data
        ! --------------------
        type(sepoperator) :: obj
        integer :: dimen, nclist
        logical :: update
        dimension :: nclist(dimen), update(dimen)

        obj%dimen  = dimen
        allocate(obj%nclist(dimen), obj%update(dimen))
        obj%nclist = nclist
        obj%update = update
        obj%ncols  = product(obj%nclist)

        allocate(obj%Mcoefs(obj%dimen, maxval(obj%nclist)))
        allocate(obj%Kcoefs(obj%dimen, maxval(obj%nclist)))
        obj%Mcoefs = 1.d0; obj%Kcoefs = 1.d0

    end subroutine initialize_operator

    subroutine separatevariables_2d(obj, CC)
        !! Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! --------------------
        integer, parameter :: dimen = 2
        type(sepoperator) :: obj
        double precision, intent(in) :: CC
        dimension :: CC(dimen, dimen, obj%ncols)

        ! Local data
        ! ----------
        integer :: nc_u, nc_v
        double precision, dimension(:, :), allocatable :: Vscript, Mscript, Nscript
        double precision, dimension(:, :, :), allocatable :: Wscript
        integer :: ju, jv, k, l, genPos, c, iter
        double precision :: vmin, vmax
        double precision :: UU(dimen), WW(dimen), WWlk(dimen-1)     
        
        if (obj%dimen.ne.dimen) stop 'Dimension problem'
        nc_u = obj%nclist(1); nc_v = obj%nclist(2)

        do iter = 1, obj%maxit
            allocate(Vscript(nc_u, nc_v))
            do k = 1, dimen
                do jv = 1, nc_v
                    do ju = 1, nc_u
                        genPos = ju + (jv-1)*nc_u
                        UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv)] 
                        Vscript(ju, jv) = CC(k, k, genPos)*UU(k)/product(UU)
                    end do
                end do

                ! Update K
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Vscript(ju, :))
                        vmax = maxval(Vscript(ju, :))
                        obj%Kcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Vscript(:, jv))
                        vmax = maxval(Vscript(:, jv))
                        obj%Kcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

            end do
            deallocate(Vscript)

            do k = 1, dimen
                allocate(Wscript(dimen-1, nc_u, nc_v))
                c = 0
                do l = 1, dimen
                    if (l.ne.k) then 
                        c = c + 1
                        do jv = 1, nc_v
                            do ju = 1, nc_u
                                genPos = ju + (jv-1)*nc_u
                                UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv)]
                                WW = [obj%Kcoefs(1, ju), obj%Kcoefs(2, jv)]
                                Wscript(c, ju, jv) = CC(k, k, genPos)*UU(k)*UU(l)&
                                                            /(product(UU)*WW(k))
                            end do
                        end do
                    end if
                end do

                ! Compute Nscript and Mscript
                allocate(Mscript(nc_u, nc_v), Nscript(nc_u, nc_v))
                    do jv = 1, nc_v
                        do ju = 1, nc_u
                            WWlk = Wscript(:, ju, jv)
                            Nscript(ju, jv) = minval(WWlk)
                            Mscript(ju, jv) = maxval(WWlk)
                        end do
                    end do
                deallocate(Wscript)

                ! Update M
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Nscript(ju, :))
                        vmax = maxval(Mscript(ju, :))
                        obj%Mcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Nscript(:, jv))
                        vmax = maxval(Mscript(:, jv))
                        obj%Mcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

                deallocate(Mscript, Nscript)
            end do
        end do

    end subroutine separatevariables_2d

    subroutine separatevariables_3d(obj, CC)
        !! Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! --------------------
        integer, parameter :: dimen = 3
        type(sepoperator) :: obj
        double precision, intent(in) :: CC
        dimension :: CC(dimen, dimen, obj%ncols)

        ! Local data
        ! ----------
        integer :: nc_u, nc_v, nc_w
        double precision, dimension(:, :, :), allocatable :: Vscript, Mscript, Nscript
        double precision, dimension(:, :, :, :), allocatable :: Wscript
        integer :: ju, jv, jw, k, l, genPos, c, iter
        double precision :: vmin, vmax
        double precision :: UU(dimen), WW(dimen), WWlk(dimen-1)     
        
        if (obj%dimen.ne.dimen) stop 'Dimension problem'
        nc_u = obj%nclist(1); nc_v = obj%nclist(2); nc_w = obj%nclist(3)

        do iter = 1, obj%maxit
            allocate(Vscript(nc_u, nc_v, nc_w))
            do k = 1, dimen
                do jw = 1, nc_w
                    do jv = 1, nc_v
                        do ju = 1, nc_u
                            genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                            UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv), obj%Mcoefs(3, jw)] 
                            Vscript(ju, jv, jw) = CC(k, k, genPos)*UU(k)/product(UU)
                        end do
                    end do
                end do

                ! Update K
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Vscript(ju, :, :))
                        vmax = maxval(Vscript(ju, :, :))
                        obj%Kcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Vscript(:, jv, :))
                        vmax = maxval(Vscript(:, jv, :))
                        obj%Kcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.3).and.(obj%update(3))) then 
                    do jw = 1, nc_w
                        vmin = minval(Vscript(:, :, jw))
                        vmax = maxval(Vscript(:, :, jw))
                        obj%Kcoefs(k, jw) = sqrt(vmin*vmax)
                    end do
                end if

            end do
            deallocate(Vscript)

            do k = 1, dimen
                allocate(Wscript(dimen-1, nc_u, nc_v, nc_w))
                c = 0
                do l = 1, dimen
                    if (l.ne.k) then 
                        c = c + 1
                        do jw = 1, nc_w
                            do jv = 1, nc_v
                                do ju = 1, nc_u
                                    genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v
                                    UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv), obj%Mcoefs(3, jw)]
                                    WW = [obj%Kcoefs(1, ju), obj%Kcoefs(2, jv), obj%Kcoefs(3, jw)]
                                    Wscript(c, ju, jv, jw) = CC(k, k, genPos)*UU(k)*UU(l)&
                                                                /(product(UU)*WW(k))
                                end do
                            end do
                        end do
                    end if
                end do

                ! Compute Nscript and Mscript
                allocate(Mscript(nc_u, nc_v, nc_w), Nscript(nc_u, nc_v, nc_w))
                    do jw = 1, nc_w
                        do jv = 1, nc_v
                            do ju = 1, nc_u
                                WWlk = Wscript(:, ju, jv, jw)
                                Nscript(ju, jv, jw) = minval(WWlk)
                                Mscript(ju, jv, jw) = maxval(WWlk)
                            end do
                        end do
                    end do
                deallocate(Wscript)

                ! Update M
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Nscript(ju, :, :))
                        vmax = maxval(Mscript(ju, :, :))
                        obj%Mcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Nscript(:, jv, :))
                        vmax = maxval(Mscript(:, jv, :))
                        obj%Mcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.3).and.(obj%update(3))) then 
                    do jw = 1, nc_w
                        vmin = minval(Nscript(:, :, jw))
                        vmax = maxval(Mscript(:, :, jw))
                        obj%Mcoefs(k, jw) = sqrt(vmin*vmax)
                    end do
                end if

                deallocate(Mscript, Nscript)
            end do
        end do

    end subroutine separatevariables_3d

    subroutine separatevariables_4d(obj, CC)
        !! Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
        !! Based on "Preconditioners for Isogemetric Analysis" by M. Montardini

        implicit none
        ! Input /  output data
        ! --------------------
        integer, parameter :: dimen = 4
        type(sepoperator) :: obj
        double precision, intent(in) :: CC
        dimension :: CC(dimen, dimen, obj%ncols)

        ! Local data
        ! ----------
        integer :: nc_u, nc_v, nc_w, nc_t
        double precision, dimension(:, :, :, :), allocatable :: Vscript, Mscript, Nscript
        double precision, dimension(:, :, :, :, :), allocatable :: Wscript
        integer :: ju, jv, jw, jt, k, l, genPos, c, iter
        double precision :: vmin, vmax
        double precision :: UU(dimen), WW(dimen), WWlk(dimen-1)     
        
        if (obj%dimen.ne.dimen) stop 'Dimension problem'
        nc_u = obj%nclist(1); nc_v = obj%nclist(2); nc_w = obj%nclist(3); nc_t = obj%nclist(4)

        do iter = 1, obj%maxit
            allocate(Vscript(nc_u, nc_v, nc_w, nc_t))
            do k = 1, dimen
                do jt = 1, nc_t
                    do jw = 1, nc_w
                        do jv = 1, nc_v
                            do ju = 1, nc_u
                                genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v + (jt-1)*nc_u*nc_v*nc_w
                                UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv), obj%Mcoefs(3, jw), obj%Mcoefs(4, jt)] 
                                Vscript(ju, jv, jw, jt) = CC(k, k, genPos)*UU(k)/product(UU)
                            end do
                        end do
                    end do
                end do

                ! Update K
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Vscript(ju, :, :, :))
                        vmax = maxval(Vscript(ju, :, :, :))
                        obj%Kcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Vscript(:, jv, :, :))
                        vmax = maxval(Vscript(:, jv, :, :))
                        obj%Kcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.3).and.(obj%update(3))) then 
                    do jw = 1, nc_w
                        vmin = minval(Vscript(:, :, jw, :))
                        vmax = maxval(Vscript(:, :, jw, :))
                        obj%Kcoefs(k, jw) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.4).and.(obj%update(4))) then 
                    do jt = 1, nc_t
                        vmin = minval(Vscript(:, :, :, jt))
                        vmax = maxval(Vscript(:, :, :, jt))
                        obj%Kcoefs(k, jt) = sqrt(vmin*vmax)
                    end do
                end if

            end do
            deallocate(Vscript)

            do k = 1, dimen
                allocate(Wscript(dimen-1, nc_u, nc_v, nc_w, nc_t))
                c = 0
                do l = 1, dimen
                    if (k.ne.l) then 
                        c = c + 1
                        do jt = 1, nc_t
                            do jw = 1, nc_w
                                do jv = 1, nc_v
                                    do ju = 1, nc_u
                                        genPos = ju + (jv-1)*nc_u + (jw-1)*nc_u*nc_v + (jt-1)*nc_u*nc_v*nc_w
                                        UU = [obj%Mcoefs(1, ju), obj%Mcoefs(2, jv), obj%Mcoefs(3, jw), obj%Mcoefs(4, jt)]
                                        WW = [obj%Kcoefs(1, ju), obj%Kcoefs(2, jv), obj%Kcoefs(3, jw), obj%Kcoefs(4, jt)]
                                        Wscript(c, ju, jv, jw, jt) = CC(k, k, genPos)*UU(k)*UU(l)&
                                                                    /(product(UU)*WW(k))
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do

                ! Compute Nscript and Mscript
                allocate(Mscript(nc_u, nc_v, nc_w, nc_t), Nscript(nc_u, nc_v, nc_w, nc_t))
                do jt = 1, nc_t
                    do jw = 1, nc_w
                        do jv = 1, nc_v
                            do ju = 1, nc_u
                                WWlk = Wscript(:, ju, jv, jw, jt)
                                Nscript(ju, jv, jw, jt) = minval(WWlk)
                                Mscript(ju, jv, jw, jt) = maxval(WWlk)
                            end do
                        end do
                    end do
                end do
                deallocate(Wscript)

                ! Update M
                if ((k.eq.1).and.(obj%update(1))) then 
                    do ju = 1, nc_u
                        vmin = minval(Nscript(ju, :, :, :))
                        vmax = maxval(Mscript(ju, :, :, :))
                        obj%Mcoefs(k, ju) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.2).and.(obj%update(2))) then
                    do jv = 1, nc_v
                        vmin = minval(Nscript(:, jv, :, :))
                        vmax = maxval(Mscript(:, jv, :, :))
                        obj%Mcoefs(k, jv) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.3).and.(obj%update(3))) then 
                    do jw = 1, nc_w
                        vmin = minval(Nscript(:, :, jw, :))
                        vmax = maxval(Mscript(:, :, jw, :))
                        obj%Mcoefs(k, jw) = sqrt(vmin*vmax)
                    end do
                end if

                if ((k.eq.4).and.(obj%update(4))) then 
                    do jt = 1, nc_t
                        vmin = minval(Nscript(:, :, :, jt))
                        vmax = maxval(Mscript(:, :, :, jt))
                        obj%Mcoefs(k, jt) = sqrt(vmin*vmax)
                    end do
                end if
                deallocate(Mscript, Nscript)
            end do
        end do

    end subroutine separatevariables_4d

end module separatevariables