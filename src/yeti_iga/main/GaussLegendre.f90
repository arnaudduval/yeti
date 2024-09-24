!! Copyright 2022 Arnaud Duval

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

!! Computes abscissas and weights for Gauss-Legendre quadrature

!! Parameters (input):
!! norder (integer) : order of the rule (between 1 and 9)
!! dim (integer) : dimension of the integration space
!! iFace (integer) : index of the face if integration is made in a subspace
!! iFace should be set to 0 for integration on full space
!! Parameters (output):
!! GaussWtCoord (double precision, dimension(dim+1, norder)) : weights (index 1) and abscissas of the rule

!! TODO : rewrite to return separate variables for weight and coordinates
subroutine Gauss (norder, dim, GaussWtCoord, iFace)
    implicit none

    integer :: norder, dim, iFace
    double precision :: GaussWtCoord
    dimension GaussWtCoord(dim + 1, norder ** dim)

    double precision zero, one
    parameter (zero = 0.d0, one = 1.d0)

    double precision xi, weight

    dimension xi(norder)
    dimension weight(norder)

    integer NumPtInt, i, j, k

    select case (norder)
        case (1)
        xi(1) =   0.0D+00

        weight(1) = 2.0D+00
        case (2)
        xi(1) = - 0.577350269189625764509148780502D+00
        xi(2) =   0.577350269189625764509148780502D+00

        weight(1) = 1.0D+00
        weight(2) = 1.0D+00
        case (3)
        xi(1) = - 0.774596669241483377035853079956D+00
        xi(2) =   0.0D+00
        xi(3) =   0.774596669241483377035853079956D+00

        weight(1) = 5.0D+00 / 9.0D+00
        weight(2) = 8.0D+00 / 9.0D+00
        weight(3) = 5.0D+00 / 9.0D+00
        case (4)
        xi(1) = - 0.861136311594052575223946488893D+00
        xi(2) = - 0.339981043584856264802665759103D+00
        xi(3) =   0.339981043584856264802665759103D+00
        xi(4) =   0.861136311594052575223946488893D+00

        weight(1) = 0.347854845137453857373063949222D+00
        weight(2) = 0.652145154862546142626936050778D+00
        weight(3) = 0.652145154862546142626936050778D+00
        weight(4) = 0.347854845137453857373063949222D+00
        case (5)
        xi(1) = - 0.906179845938663992797626878299D+00
        xi(2) = - 0.538469310105683091036314420700D+00
        xi(3) =   0.0D+00
        xi(4) =   0.538469310105683091036314420700D+00
        xi(5) =   0.906179845938663992797626878299D+00

        weight(1) = 0.236926885056189087514264040720D+00
        weight(2) = 0.478628670499366468041291514836D+00
        weight(3) = 0.568888888888888888888888888889D+00
        weight(4) = 0.478628670499366468041291514836D+00
        weight(5) = 0.236926885056189087514264040720D+00
        case (6)
        xi(1) = - 0.932469514203152027812301554494D+00
        xi(2) = - 0.661209386466264513661399595020D+00
        xi(3) = - 0.238619186083196908630501721681D+00
        xi(4) =   0.238619186083196908630501721681D+00
        xi(5) =   0.661209386466264513661399595020D+00
        xi(6) =   0.932469514203152027812301554494D+00

        weight(1) = 0.171324492379170345040296142173D+00
        weight(2) = 0.360761573048138607569833513838D+00
        weight(3) = 0.467913934572691047389870343990D+00
        weight(4) = 0.467913934572691047389870343990D+00
        weight(5) = 0.360761573048138607569833513838D+00
        weight(6) = 0.171324492379170345040296142173D+00
        case (7)
        xi(1) = - 0.949107912342758524526189684048D+00
        xi(2) = - 0.741531185599394439863864773281D+00
        xi(3) = - 0.405845151377397166906606412077D+00
        xi(4) =   0.0D+00
        xi(5) =   0.405845151377397166906606412077D+00
        xi(6) =   0.741531185599394439863864773281D+00
        xi(7) =   0.949107912342758524526189684048D+00

        weight(1) = 0.129484966168869693270611432679D+00
        weight(2) = 0.279705391489276667901467771424D+00
        weight(3) = 0.381830050505118944950369775489D+00
        weight(4) = 0.417959183673469387755102040816D+00
        weight(5) = 0.381830050505118944950369775489D+00
        weight(6) = 0.279705391489276667901467771424D+00
        weight(7) = 0.129484966168869693270611432679D+00
        case (8)
        xi(1) = - 0.960289856497536231683560868569D+00
        xi(2) = - 0.796666477413626739591553936476D+00
        xi(3) = - 0.525532409916328985817739049189D+00
        xi(4) = - 0.183434642495649804939476142360D+00
        xi(5) =   0.183434642495649804939476142360D+00
        xi(6) =   0.525532409916328985817739049189D+00
        xi(7) =   0.796666477413626739591553936476D+00
        xi(8) =   0.960289856497536231683560868569D+00

        weight(1) = 0.101228536290376259152531354310D+00
        weight(2) = 0.222381034453374470544355994426D+00
        weight(3) = 0.313706645877887287337962201987D+00
        weight(4) = 0.362683783378361982965150449277D+00
        weight(5) = 0.362683783378361982965150449277D+00
        weight(6) = 0.313706645877887287337962201987D+00
        weight(7) = 0.222381034453374470544355994426D+00
        weight(8) = 0.101228536290376259152531354310D+00
        case (9)
        xi(1) = - 0.968160239507626089835576202904D+00
        xi(2) = - 0.836031107326635794299429788070D+00
        xi(3) = - 0.613371432700590397308702039341D+00
        xi(4) = - 0.324253423403808929038538014643D+00
        xi(5) =   0.0D+00
        xi(6) =   0.324253423403808929038538014643D+00
        xi(7) =   0.613371432700590397308702039341D+00
        xi(8) =   0.836031107326635794299429788070D+00
        xi(9) =   0.968160239507626089835576202904D+00

        weight(1) = 0.812743883615744119718921581105D-01
        weight(2) = 0.180648160694857404058472031243D+00
        weight(3) = 0.260610696402935462318742869419D+00
        weight(4) = 0.312347077040002840068630406584D+00
        weight(5) = 0.330239355001259763164525069287D+00
        weight(6) = 0.312347077040002840068630406584D+00
        weight(7) = 0.260610696402935462318742869419D+00
        weight(8) = 0.180648160694857404058472031243D+00
        weight(9) = 0.812743883615744119718921581105D-01
    end select

    GaussWtCoord(:,:) = zero

    select case (iFace + dim * 10)
        case (10)
        NumPtInt = 0
        do i = 1, norder
            NumPtInt = NumPtInt + 1
            GaussWtCoord(1, NumPtInt) = weight(i)
            GaussWtCoord(2, NumPtInt) = xi(i)
        end do
        case (20)
        NumPtInt = 0
        do i = 1, norder
            do j = 1, norder
                NumPtInt = NumPtInt + 1
                GaussWtCoord(1, NumPtInt) = weight(i) * weight(j)
                GaussWtCoord(2, NumPtInt) = xi(i)
                GaussWtCoord(3, NumPtInt) = xi(j)
            end do
        end do
        case (30)
        NumPtInt = 0
        do i = 1, norder
            do j = 1, norder
                do k = 1, norder
                    NumPtInt = NumPtInt + 1
                    GaussWtCoord(1, NumPtInt) = weight(i) * weight(j) * weight(k)
                    GaussWtCoord(2, NumPtInt) = xi(i)
                    GaussWtCoord(3, NumPtInt) = xi(j)
                    GaussWtCoord(4, NumPtInt) = xi(k)
                end do
            end do
        end do

        !! Creation de la matrice GaussFace qui correspond aux points en
        !! surface de la geometrie
        ! 2D cases
        case (21)                 ! Cote 1
        do i = 1, norder
            GaussWtCoord(1, i) = weight(i)
            GaussWtCoord(2, i) = - one
            GaussWtCoord(3, i) = xi(i)
        end do
        case (22)                 ! Cote 2
        do i = 1, norder
            GaussWtCoord(1, i) = weight(i)
            GaussWtCoord(2, i) = one
            GaussWtCoord(3, i) = xi(i)
        end do
        case (23)                 ! Cote 3
        do i = 1, norder
            GaussWtCoord(1, i) = weight(i)
            GaussWtCoord(2, i) = xi(i)
            GaussWtCoord(3, i) = - one
        end do
        case (24)                 ! Cote 4
        do i = 1, norder
            GaussWtCoord(1, i) = weight(i)
            GaussWtCoord(2, i) = xi(i)
            GaussWtCoord(3, i) = one
        end do
        ! 3D cases
        case (31)                 ! Face 1
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = - one
                GaussWtCoord(3, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(4, i + norder * (j - 1)) = xi(j)
            end do
        end do
        case (32)                 ! Face 2
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = one
                GaussWtCoord(3, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(4, i + norder * (j - 1)) = xi(j)
            end do
        end do
        case (33)                 ! Face 3
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(3, i + norder * (j - 1)) = - one
                GaussWtCoord(4, i + norder * (j - 1)) = xi(j)
            end do
        end do
        case (34)                 ! Face 4
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(3, i + norder * (j - 1)) = one
                GaussWtCoord(4, i + norder * (j - 1)) = xi(j)
            end do
        end do
        case (35)                 ! Face 5
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = xi(j)
                GaussWtCoord(3, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(4, i + norder * (j - 1)) = - one
            end do
        end do
        case (36)                 ! Face 6
        do j = 1, norder
            do i = 1, norder
                GaussWtCoord(1, i + norder * (j - 1)) = weight(i) * weight(j)
                GaussWtCoord(2, i + norder * (j - 1)) = xi(j)
                GaussWtCoord(3, i + norder * (j - 1)) = xi(i)
                GaussWtCoord(4, i + norder * (j - 1)) = one
            end do
        end do
    end select


    return
end subroutine Gauss

! !! Routine for Python binding defining number of integration points explicitely
! subroutine GaussPts(norder,dim,GaussWtCoord,iFace,n_gps)
! implicit none

! integer, intent(in) :: norder, dim, iFace,n_gps
! double precision, intent(out) :: GaussWtCoord
! dimension GaussWtCoord(dim+1,n_gps)

! call Gauss(norder, dim, GaussWtCoord, iFace)

! end subroutine GaussPts
