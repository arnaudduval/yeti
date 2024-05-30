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

!! Compute shape functions and their 1st derivative at given point
!! Given point has coordinates in isoparametric element [-1 1]

subroutine shap(dRdx, R, DetJac, COORDS, PtGauss, MCRD)

    use parameters
    use nurbspatch

    implicit none

    !! Input arguments
    !! ---------------
    integer :: MCRD

    double precision :: COORDS, PtGauss
    dimension COORDS(MCRD, nnode_patch), PtGauss(MCRD)

    !! Output variables
    !! ----------------
    double precision :: R, dRdx, DetJac
    dimension dRdx(MCRD, nnode_patch), R(nnode_patch)


    !! Local variables
    !! ---------------
    double precision :: FN, FM, FL, dNdxi, dMdEta, dLdZeta, dxdxi, dxidx,    &
        &     dXidtildexi, AJmat, xi, dRdxi, Detdxdxi, SumTot,  &
        &     SumTot_inv, SumXi
    dimension FN(Jpqr_patch(1) + 1), dNdxi(Jpqr_patch(1) + 1),   &
        &     FM(Jpqr_patch(2) + 1), dMdEta(Jpqr_patch(2) + 1),  &
        &     FL(Jpqr_patch(3) + 1), dLdZeta(Jpqr_patch(3) + 1), &
        &     dxdxi(MCRD, MCRD), dxidx(MCRD, MCRD), dXidtildexi(MCRD, MCRD),     &
        &     AJmat(MCRD, MCRD), xi(dim_patch), dRdxi(nnode_patch, 3), SumXi(3)

    Integer ::  i, j, k, Ni, NumLoc, Na, Nb, Nc
    dimension Ni(dim_patch)


    !! Initialization

    !! 1D and 2D cases
    FM(1)      = one
    dMdEta(1)  = zero
    FL(1)      = one
    dLdZeta(1) = zero

    !! Knot spans of current element
    do i = 1, dim_patch
        Ni(i) = Nijk_patch(i, current_elem)
    end do

    dxdxi(:,:) = zero
    dRdx(:,:)  = zero
    dXidtildexi(:,:) = zero
    AJmat(:,:) = zero

    !! Compute parametric coordinates of the given point
    do i = 1, dim_patch
        xi(i) = ((Ukv_elem(2, i) - Ukv_elem(1, i)) * PtGauss(i)  &
            &+  (Ukv_elem(2, i) + Ukv_elem(1, i)) ) * 0.5d0
    end do

    !! Compute univariate B-Spline function
    call dersbasisfuns(Ni(1), Jpqr_patch(1), Nkv_patch(1), xi(1),  &
        &   Ukv1_patch(:), FN, dNdxi)
    call dersbasisfuns(Ni(2), Jpqr_patch(2), Nkv_patch(2), xi(2),  &
        &   Ukv2_patch(:), FM, dMdeta)
    if (dim_patch == 3) then
        call dersbasisfuns(Ni(3), Jpqr_patch(3), Nkv_patch(3), xi(3),  &
            &   Ukv3_patch(:), FL, dLdZeta)
    end if

    !! Build numerators and denominators
    NumLoc   = 0
    SumTot   = zero
    SumXi(:) = zero
    do k = 0, Jpqr_patch(3)
        do j = 0, Jpqr_patch(2)
            do i = 0, Jpqr_patch(1)
                NumLoc = NumLoc + 1
                R(NumLoc) = FN(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j) &
                    &   * FL(Jpqr_patch(3) + 1 - k) * Weight_elem(NumLoc)

                SumTot = SumTot + R(NumLoc)

                dRdxi(NumLoc, 1) =   &
                    &   dNdxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)  &
                    &   * FL(Jpqr_patch(3) + 1 - k) * Weight_elem(NumLoc)
                SumXi(1) = SumXi(1) + dRdxi(NumLoc, 1)

                dRdxi(NumLoc, 2) =   &
                    &   FN(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j) &
                    &   * FL(Jpqr_patch(3) + 1 - k) * Weight_elem(NumLoc)
                SumXi(2) = SumXi(2) + dRdxi(NumLoc, 2)

                dRdxi(NumLoc, 3) =   &
                    &   FN(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j) &
                    &   * dLdZeta(Jpqr_patch(3) + 1 - k) * Weight_elem(NumLoc)
                SumXi(3) = SumXi(3) + dRdxi(NumLoc, 3)
            end do
        end do
    end do

    !! Divide by denominator to complete definition of fct and deriv.
    SumTot_inv = one / SumTot
    do NumLoc = 1, nnode_patch
        R(NumLoc) = R(NumLoc) * SumTot_inv
        do i = 1, MCRD
            dRdxi(NumLoc, i) &
                &   = (dRdxi(NumLoc, i) - R(NumLoc) * SumXi(i)) * SumTot_inv
        end do
    end do


    !! Gradient of mapping from parameter space to physical space
    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dxdxi(Na, Nb) = dxdxi(Na, Nb) &
                    &   + COORDS(Na, NumLoc) * dRdxi(NumLoc, Nb)
            end do
        end do
    end do

    !! Compute inverse of gradient
    call MatrixInv(dxidx, dxdxi, Detdxdxi, MCRD)

    !! Compute derivatives of basis functions with respect to physical coordinates
    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dRdx(Na, NumLoc) = dRdx(Na, NumLoc)   &
                    &   + dRdxi(NumLoc, Nb) * dxidx(Nb, Na)
            end do
        end do
    end do

    !! Gradient of mapping from parent element to parameter space
    do i = 1, dim_patch
        dXidtildexi(i, i) = 0.5d0 * ( Ukv_elem(2, i) - Ukv_elem(1, i) )
    end do

    do Na = 1, MCRD
        do Nb = 1, MCRD
            do Nc = 1, MCRD
                AJmat(Na, Nb) = AJmat(Na, Nb) + dxdxi(Na, Nc) * dXidtildexi(Nc, Nb)
            end do
        end do
    end do


    if (dim_patch == 2) then
        DetJac = AJmat(1, 1) * AJmat(2, 2) - AJmat(2, 1) * AJmat(1, 2)
    else
        DetJac = AJmat(1, 1) * AJmat(2, 2) * AJmat(3, 3)   &
            &   + AJmat(1, 2) * AJmat(2, 3) * AJmat(3, 1)  &
            &   + AJmat(2, 1) * AJmat(3, 2) * AJmat(1, 3)  &
            &   - AJmat(1, 3) * AJmat(2, 2) * AJmat(3, 1)  &
            &   - AJmat(1, 2) * AJmat(2, 1) * AJmat(3, 3)  &
            &   - AJmat(2, 3) * AJmat(3, 2) * AJmat(1, 1)
    end if

end subroutine shap
