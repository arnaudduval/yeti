!! Copyright 2010 Thomas Elguedj
!! Copyright 2011 Florian Maurin
!! Copyright 2016-2018 Thibaut Hirschler
!! Copyright 2021 Arnaud Duval

!! This file is part of Yeti.
!!
!! Yeti is free software: you can redistribute it and/or modify it under the terms
!! of the GNU Lesser General Public License as published by the Free Software
!! Foundation, either version 3 of the License, or (at your option) any later version.
!!
!! Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!! PURPOSE. See the GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License along
!! with Yeti. If not, see <https://www.gnu.org/licenses/>


!! This file contains the basis functions and derivatives routine.
!! The algorithm is taken from Piegl, Les. "The NURBS Book". Springer-Verlag:
!! Berlin 1995; p. 72-73.


!! ...  The routine consumes a knot index, parameter value, and a knot
!! ...  vector and returns a vector containing all nonzero 1D b-spline shape
!! ...  functions evaluated at that parameter as well as their derivatives.


subroutine dersbasisfuns(ni, Jp, NKv, u, u_knotl, ders1, ders2)
    implicit none

    !! Input parameters
    !! ----------------
    !! knot span, degree of curve, number of control points, counters
    !! parameter value, vector of knots, derivative matrix

    integer, intent(in) :: Jp, NKv, ni
    double precision, intent(in) :: u, u_knotl
    double precision, intent(out) :: ders1, ders2


    !! Local variables
    !! ---------------
    double precision :: Aleft, right, ders, Andu
    integer :: nders
    dimension :: Aleft(Jp + 1)
    dimension :: right(Jp + 1)
    dimension :: a(2, Jp + 1)
    dimension :: u_knotl(NKv)
    dimension :: ders(2, Jp + 1)
    dimension :: Andu(Jp + 1, Jp + 1)
    dimension :: ders1(Jp + 1), ders2(Jp + 1)
    integer :: j, k, l, Kr, Kp, j1, j2, ls1, ls2
    double precision :: u_knot1, saved, temp, A, D

    nders = 1

    Andu(1, 1) = 1.d0
    do j = 1, Jp
        Aleft(j + 1) = u - u_knotl(ni + 1 - j)
        right(j + 1) = u_knotl(ni + j) - u
        saved = 0.d0
        do l = 0, j - 1
            Andu(j + 1, l + 1) = right(l + 2) + Aleft(j - l + 1)
            temp = Andu(l + 1, j) / Andu(j + 1, l + 1)
            Andu(l + 1, j + 1) = saved + right(l + 2) * temp
            saved = Aleft(j - l + 1) * temp
        end do
        Andu(j + 1, j + 1) = saved
    end do

    !! load basis functions
    do j = 0, Jp
        ders(1, j + 1) = Andu(j + 1, Jp + 1)
    end do

    !! compute derivatives
    do l = 0, Jp ! loop over function index
        ls1 = 0
        ls2 = 1                 ! alternate rows in array a
        a(1, 1) = 1.d0
        ! loop to compute kth derivative
        do k = 1, nders
            d = 0.d0
            kr = l - k
            kp = Jp - k
            if (l >= k) then
                a(ls2 + 1, 1) = a(ls1 + 1, 1) / Andu(kp + 2, kr + 1)
                d = a(ls2 + 1, 1) * Andu(kr + 1, kp + 1)
            end if
            if (kr >= - 1) then
                j1 = 1
            else
                j1 = - kr
            end if
            if ((l - 1) <= kp) then
                j2 = k - 1
            else
                j2 = Jp - l
            end if
            do j = j1, j2
                a(ls2 + 1, j + 1) = (a(ls1 + 1, j + 1) - a(ls1 + 1, j)) / Andu(kp + 2, kr + j + 1)
                d = d + a(ls2 + 1, j + 1) * Andu(kr + j + 1, kp + 1)
            end do
            if (l <= kp) then
                a(ls2 + 1, k + 1) = - a(ls1 + 1, k) / Andu(kp + 2, l + 1)
                d = d + a(ls2 + 1, k + 1) * Andu(l + 1, kp + 1)
            end if
            ders(k + 1, l + 1) = d
            j = ls1
            ls1 = ls2
            ls2 = j              ! switch rows
        end do
    end do

    !! Multiply through by the correct factors
    l = Jp
    do k = 1, nders
        do j = 0, Jp
            ders(k + 1, j + 1) = ders(k + 1, j + 1) * l
        end do
        l = l * (Jp - k)
    end do

    do j = 1, Jp + 1
        ders1(j) = ders(1, j)
        ders2(j) = ders(2, j)
    end do

end subroutine dersbasisfuns


subroutine dersbasisfuns2(ni, Jp, NKv, u, u_knotl, ders1, ders2, ders3)
    implicit none

    !! Input variables
    !! ---------------
    !! knot span, degree of curve, number of control points, counters
    !! parameter value, vector of knots, derivative matrix

    integer, intent(in) :: Jp, NKv, ni
    double precision, intent(in) :: u, u_knotl
    double precision, intent(out) :: ders1, ders2, ders3

    !! Local variables
    !! ---------------
    double precision :: Aleft, right, ders, Andu
    integer ::  nders

    dimension Aleft(Jp + 1)
    dimension right(Jp + 1)
    dimension a(2, Jp + 1)
    dimension u_knotl(NKv)
    dimension ders(3, Jp + 1)
    dimension Andu(Jp + 1, Jp + 1)
    dimension ders1(Jp + 1), ders2(Jp + 1), ders3(Jp + 1)

    integer j, k, l, Kr, Kp, j1, j2, ls1, ls2
    double precision saved, temp, A, D

    nders = 2

    Andu(1, 1) = 1.d0
    do j = 1, Jp
        Aleft(j + 1) = u - u_knotl(ni + 1 - j)
        right(j + 1) = u_knotl(ni + j) - u
        saved = 0.d0
        do l = 0, j - 1
            Andu(j + 1, l + 1) = right(l + 2) + Aleft(j - l + 1)
            temp = Andu(l + 1, j) / Andu(j + 1, l + 1)
            Andu(l + 1, j + 1) = saved + right(l + 2) * temp
            saved = Aleft(j - l + 1) * temp
        end do
        Andu(j + 1, j + 1) = saved
    end do

    !! Load basis functions
    do j = 0, Jp
        ders(1, j + 1) = Andu(j + 1, Jp + 1)
    end do
    !! Compute derivatives
    do l = 0, Jp               ! loop over function index
        ls1 = 0
        ls2 = 1                ! alternate rows in array a
        a(1, 1) = 1.d0
        !! loop to compute kth derivative
        do k = 1, nders
            d = 0.d0
            kr = l - k
            kp = Jp - k
            if (l >= k) then
                a(ls2 + 1, 1) = a(ls1 + 1, 1) / Andu(kp + 2, kr + 1)
                d = a(ls2 + 1, 1) * Andu(kr + 1, kp + 1)
            end if
            if (kr >= - 1) then
                j1 = 1
            else
                j1 = - kr
            end if
            if ((l - 1) <= kp) then
                j2 = k - 1
            else
                j2 = Jp - l
            end if
            do j = j1, j2
                a(ls2 + 1, j + 1) = (a(ls1 + 1, j + 1) - a(ls1 + 1, j)) / Andu(kp + 2, kr + j + 1)
                d = d + a(ls2 + 1, j + 1) * Andu(kr + j + 1, kp + 1)
            end do
            if (l <= kp) then
                a(ls2 + 1, k + 1) = - a(ls1 + 1, k) / Andu(kp + 2, l + 1)
                d = d + a(ls2 + 1, k + 1) * Andu(l + 1, kp + 1)
            end if
            ders(k + 1, l + 1) = d
            j = ls1
            ls1 = ls2
            ls2 = j             ! switch rows
        end do
    end do

    !! Multiply through by the correct factors
    l = Jp
    do k = 1, nders
        do j = 0, Jp
            ders(k + 1, j + 1) = ders(k + 1, j + 1) * l
        end do
        l = l * (Jp - k)
    end do

    do j = 1, Jp + 1
        ders1(j) = ders(1, j)
        ders2(j) = ders(2, j)
        ders3(j) = ders(3, j)
    end do

end subroutine dersbasisfuns2

subroutine dersbasisfuns3(ni, Jp, NKv, u, u_knotl, ders1, ders2, ders3, ders4)
    implicit none

    !! Input variables
    !! ---------------
    !! knot span, degree of curve, number of control points, counters
    !! parameter value, vector of knots, derivative matrix

    integer, intent(in) :: Jp, NKv, ni
    double precision, intent(in) :: u, u_knotl
    double precision, intent(out) :: ders1, ders2, ders3, ders4

    !! Local variables
    !! ---------------
    double precision :: Aleft, right, ders, Andu
    integer ::  nders

    dimension Aleft(Jp + 1)
    dimension right(Jp + 1)
    dimension a(2, Jp + 1)
    dimension u_knotl(NKv)
    dimension ders(4, Jp + 1)
    dimension Andu(Jp + 1, Jp + 1)
    dimension ders1(Jp + 1), ders2(Jp + 1), ders3(Jp + 1), ders4(Jp + 1)

    integer j, k, l, Kr, Kp, j1, j2, ls1, ls2
    double precision saved, temp, A, D

    nders = 3
    ders(:,:) = 0.D0

    Andu(1, 1) = 1.d0
    do j = 1, Jp
        Aleft(j + 1) = u - u_knotl(ni + 1 - j)
        right(j + 1) = u_knotl(ni + j) - u
        saved = 0.d0
        do l = 0, j - 1
            Andu(j + 1, l + 1) = right(l + 2) + Aleft(j - l + 1)
            temp = Andu(l + 1, j) / Andu(j + 1, l + 1)
            Andu(l + 1, j + 1) = saved + right(l + 2) * temp
            saved = Aleft(j - l + 1) * temp
        end do
        Andu(j + 1, j + 1) = saved
    end do

    !! Load basis functions
    do j = 0, Jp
        ders(1, j + 1) = Andu(j + 1, Jp + 1)
    end do
    !! Compute derivatives
    do l = 0, Jp               ! loop over function index
        ls1 = 0
        ls2 = 1                ! alternate rows in array a
        a(1, 1) = 1.d0
        !! loop to compute kth derivative
        do k = 1, nders
            d = 0.d0
            kr = l - k
            kp = Jp - k
            if (l >= k) then
                a(ls2 + 1, 1) = a(ls1 + 1, 1) / Andu(kp + 2, kr + 1)
                d = a(ls2 + 1, 1) * Andu(kr + 1, kp + 1)
            end if
            if (kr >= - 1) then
                j1 = 1
            else
                j1 = - kr
            end if
            if ((l - 1) <= kp) then
                j2 = k - 1
            else
                j2 = Jp - l
            end if
            do j = j1, j2
                a(ls2 + 1, j + 1) = (a(ls1 + 1, j + 1) - a(ls1 + 1, j)) / Andu(kp + 2, kr + j + 1)
                d = d + a(ls2 + 1, j + 1) * Andu(kr + j + 1, kp + 1)
            end do
            if (l <= kp) then
                a(ls2 + 1, k + 1) = - a(ls1 + 1, k) / Andu(kp + 2, l + 1)
                d = d + a(ls2 + 1, k + 1) * Andu(l + 1, kp + 1)
            end if
            ders(k + 1, l + 1) = d
            j = ls1
            ls1 = ls2
            ls2 = j             ! switch rows
        end do
    end do

    !! Multiply through by the correct factors
    l = Jp
    do k = 1, nders
        do j = 0, Jp
            ders(k + 1, j + 1) = ders(k + 1, j + 1) * l
        end do
        l = l * (Jp - k)
    end do

    do j = 1, Jp + 1
        ders1(j) = ders(1, j)
        ders2(j) = ders(2, j)
        ders3(j) = ders(3, j)
        ders4(j) = ders(4, j)
    end do

end subroutine dersbasisfuns3


subroutine gen_dersbasisfuns(ni, Jp, NKv, u, u_knotl, nders, ders)
    !! Generic function for computation basis function and derivatives up
    !! to a given order
    !!!! TODO : reecrire les autres fonctions specifiques pour appeler la fonction generique
    implicit none

    !! Input variables
    !! ---------------
    !! ni : knot span
    !! Jp : degree of curve
    !! Nkv : knot vector size
    !! u : parameter value
    !! u_knotl : knot vector
    !! nders : derivation order
    !! ders : derivatives matrix

    integer, intent(in) :: Jp, NKv, ni
    double precision, intent(in) :: u, u_knotl
    integer, intent(in) ::  nders
    double precision, intent(out) :: ders

    !! Local variables
    !! ---------------
    double precision :: Aleft, right, Andu


    dimension Aleft(Jp + 1)
    dimension right(Jp + 1)
    dimension a(2, Jp + 1)
    dimension u_knotl(NKv)
    dimension ders(nders + 1, Jp + 1)
    dimension Andu(Jp + 1, Jp + 1)


    integer j, k, l, Kr, Kp, j1, j2, ls1, ls2
    double precision saved, temp, A, D

    Andu(1, 1) = 1.d0
    do j = 1, Jp
        Aleft(j + 1) = u - u_knotl(ni + 1 - j)
        right(j + 1) = u_knotl(ni + j) - u
        saved = 0.d0
        do l = 0, j - 1
            Andu(j + 1, l + 1) = right(l + 2) + Aleft(j - l + 1)
            temp = Andu(l + 1, j) / Andu(j + 1, l + 1)
            Andu(l + 1, j + 1) = saved + right(l + 2) * temp
            saved = Aleft(j - l + 1) * temp
        end do
        Andu(j + 1, j + 1) = saved
    end do

    !! Load basis functions
    do j = 0, Jp
        ders(1, j + 1) = Andu(j + 1, Jp + 1)
    end do
    !! Compute derivatives
    do l = 0, Jp               ! loop over function index
        ls1 = 0
        ls2 = 1                ! alternate rows in array a
        a(1, 1) = 1.d0
        !! loop to compute kth derivative
        do k = 1, nders
            d = 0.d0
            kr = l - k
            kp = Jp - k
            if (l >= k) then
                a(ls2 + 1, 1) = a(ls1 + 1, 1) / Andu(kp + 2, kr + 1)
                d = a(ls2 + 1, 1) * Andu(kr + 1, kp + 1)
            end if
            if (kr >= - 1) then
                j1 = 1
            else
                j1 = - kr
            end if
            if ((l - 1) <= kp) then
                j2 = k - 1
            else
                j2 = Jp - l
            end if
            do j = j1, j2
                a(ls2 + 1, j + 1) = (a(ls1 + 1, j + 1) - a(ls1 + 1, j)) / Andu(kp + 2, kr + j + 1)
                d = d + a(ls2 + 1, j + 1) * Andu(kr + j + 1, kp + 1)
            end do
            if (l <= kp) then
                a(ls2 + 1, k + 1) = - a(ls1 + 1, k) / Andu(kp + 2, l + 1)
                d = d + a(ls2 + 1, k + 1) * Andu(l + 1, kp + 1)
            end if
            ders(k + 1, l + 1) = d
            j = ls1
            ls1 = ls2
            ls2 = j             ! switch rows
        end do
    end do

    !! Multiply through by the correct factors
    l = Jp
    do k = 1, nders
        do j = 0, Jp
            ders(k + 1, j + 1) = ders(k + 1, j + 1) * l
        end do
        l = l * (Jp - k)
    end do

end subroutine gen_dersbasisfuns
