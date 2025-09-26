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

subroutine shap_omp(patch, element, dRdx,R,DetJac,COORDS,PtGauss,MCRD)
    use m_nurbspatch, only : nurbselement, nurbspatch

    implicit none

    !! Input arguments
    !! ---------------
    type(nurbspatch), intent(in) :: patch
    type(nurbselement), intent(in) :: element

    integer :: MCRD

    double precision :: COORDS,PtGauss
    dimension COORDS(MCRD, patch%nnode_patch),PtGauss(MCRD)

    !! Output variables
    !! ----------------
    double precision :: R,dRdx,DetJac
    dimension dRdx(MCRD,patch%nnode_patch),R(patch%nnode_patch)


    !! Local variables
    !! ---------------
    double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, dxdxi,dxidx,    &
        &     dXidtildexi, AJmat, xi, dRdxi, Detdxdxi, SumTot,  &
        &     SumTot_inv, SumXi
    dimension FN(patch%Jpqr_patch(1)+1),dNdxi(patch%Jpqr_patch(1)+1),   &
        &     FM(patch%Jpqr_patch(2)+1),dMdEta(patch%Jpqr_patch(2)+1),  &
        &     FL(patch%Jpqr_patch(3)+1),dLdZeta(patch%Jpqr_patch(3)+1), &
        &     dxdxi(MCRD,MCRD),dxidx(MCRD,MCRD),dXidtildexi(MCRD,MCRD),     &
        &     AJmat(MCRD,MCRD),xi(patch%dim_patch),dRdxi(patch%nnode_patch,3),SumXi(3)

    integer ::  i,j,k, Ni, NumLoc, Na,Nb,Nc
    dimension Ni(patch%dim_patch)

    !! Initialization

    !! 1D and 2D cases
    FM(1)      = 1.0
    dMdEta(1)  = 0.0
    FL(1)      = 1.0
    dLdZeta(1) = 0.0

    !! Knot spans of current element
    do i = 1, patch%dim_patch
        Ni(i) = patch%Nijk_patch(i, element%current_elem)
    enddo

    dxdxi(:,:) = 0.0
    dRdx(:,:)  = 0.0
    dXidtildexi(:,:) = 0.0
    AJmat(:,:) = 0.0

    !! Compute parametric coordinates of the given point
    do i = 1, patch%dim_patch
        xi(i)= ((element%Ukv_elem(2,i) - element%Ukv_elem(1,i))*PtGauss(i)  &
            &+  (element%Ukv_elem(2,i) + element%Ukv_elem(1,i)) ) * 0.5d0
    enddo

    !! Compute univariate B-Spline function
    call dersbasisfuns(Ni(1),patch%Jpqr_patch(1),patch%Nkv_patch(1),xi(1),  &
        &   patch%Ukv1_patch(:),FN,dNdxi)
    call dersbasisfuns(Ni(2),patch%Jpqr_patch(2),patch%Nkv_patch(2),xi(2),  &
        &   patch%Ukv2_patch(:),FM,dMdeta)
    if (patch%dim_patch .eq. 3) then
        call dersbasisfuns(Ni(3),patch%Jpqr_patch(3),patch%Nkv_patch(3),xi(3),  &
            &   patch%Ukv3_patch(:),FL,dLdZeta)
    endif

    !! Build numerators and denominators
    NumLoc   = 0
    SumTot   = 0.0
    SumXi(:) = 0.0
    do k = 0, patch%Jpqr_patch(3)
        do j = 0, patch%Jpqr_patch(2)
            do i = 0, patch%Jpqr_patch(1)
                NumLoc = NumLoc+1
                R(NumLoc) = FN(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j) &
                    &   *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)

                SumTot = SumTot + R(NumLoc)

                dRdxi(NumLoc,1) =   &
                    &   dNdxi(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j)  &
                    &   *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)

                dRdxi(NumLoc,2) =   &
                    &   FN(patch%Jpqr_patch(1)+1-i)*dMdEta(patch%Jpqr_patch(2)+1-j) &
                    &   *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)

                dRdxi(NumLoc,3) =   &
                    &   FN(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j) &
                    &   *dLdZeta(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
            enddo
        enddo
    enddo

    !! Divide by denominator to complete definition of fct and deriv.
    SumTot_inv = 1.0/SumTot
    do NumLoc = 1,patch%nnode_patch
        R(NumLoc) = R(NumLoc)*SumTot_inv
        do i = 1,MCRD
            dRdxi(NumLoc,i) &
                &   = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
        enddo
    enddo


    !! Gradient of mapping from parameter space to physical space
    do NumLoc= 1,patch%nnode_patch
        do Na = 1,MCRD
            do Nb = 1,MCRD
                dxdxi(Na,Nb) = dxdxi(Na,Nb) &
                    &   + COORDS(Na,NumLoc)*dRdxi(NumLoc,Nb)
            enddo
        enddo
    enddo

    !! Compute inverse of gradient
    call MatrixInv(dxidx, dxdxi,Detdxdxi, MCRD)

    !! Compute derivatives of basis functions with respect to physical coordinates
    do NumLoc= 1,patch%nnode_patch
        do Na = 1,MCRD
            do Nb = 1,MCRD
                dRdx(Na,NumLoc) = dRdx(Na,NumLoc)   &
                    &   + dRdxi(NumLoc,Nb)*dxidx(Nb,Na)
            enddo
        enddo
    enddo

    !! Gradient of mapping from parent element to parameter space
    do i = 1,patch%dim_patch
        dXidtildexi(i,i) = 0.5d0*( element%Ukv_elem(2,i) - element%Ukv_elem(1,i) )
    enddo

    do Na = 1,MCRD
        do Nb = 1,MCRD
            do Nc = 1,MCRD
                AJmat(Na,Nb) = AJmat(Na,Nb)+dxdxi(Na,Nc)*dXidtildexi(Nc,Nb)
            enddo
        enddo
    enddo


    if (patch%dim_patch==2) then
        DetJac = AJmat(1,1)*AJmat(2,2) - AJmat(2,1)*AJmat(1,2)
    else
        DetJac = AJmat(1,1)*AJmat(2,2)*AJmat(3,3)   &
            &   + AJmat(1,2)*AJmat(2,3)*AJmat(3,1)  &
            &   + AJmat(2,1)*AJmat(3,2)*AJmat(1,3)  &
            &   - AJmat(1,3)*AJmat(2,2)*AJmat(3,1)  &
            &   - AJmat(1,2)*AJmat(2,1)*AJmat(3,3)  &
            &   - AJmat(2,3)*AJmat(3,2)*AJmat(1,1)
    endif

end subroutine shap_omp
