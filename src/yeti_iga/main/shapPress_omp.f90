!! Copyright 2011 Florian Maurin
!! Copyright 2016-2018 Thibaut Hirschler
!! Copyright 2020 Arnaud Duval

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

subroutine shapPress_omp (patch, element, R, Vect, det_jac, COORDS, PtGauss, mcrd, k_numface, k_typedload)
    use m_nurbspatch, only : nurbselement, nurbspatch

    implicit none

    !! Input arguments
    !! ---------------
    type(nurbspatch), intent(in) :: patch
    type(nurbselement), intent(in) :: element
    integer :: mcrd, k_numface, k_typedload

    double precision, dimension(mcrd, patch%nnode_patch), intent(in) :: COORDS
    double precision, dimension(mcrd), intent(in) :: PtGauss

    !! Output variables
    !! ----------------
    double precision, dimension(patch%nnode_patch) :: R
    double precision, dimension(mcrd) :: Vect
    double precision :: det_jac

    !! Local variables
    !! ---------------
    double precision, dimension(patch%Jpqr_patch(1)+1) :: FN, dNdxi
    double precision, dimension(patch%Jpqr_patch(2)+1) :: FM, dMdEta
    double precision, dimension(patch%Jpqr_patch(3)+1) :: FL, dLdZeta
    double precision, dimension(mcrd, mcrd) :: dxdxi, dXidtildexi, AJmat
    double precision, dimension(patch%dim_patch) :: xi
    double precision, dimension(patch%nnode_patch, 3) :: dRdxi
    double precision, dimension(3) :: SumXi
    double precision :: Detdxdxi, SumTot, SumTot_inv

    integer ::  i,j,k, NumLoc, Na,Nb,Nc
    integer, dimension(patch%dim_patch) :: Ni

    !! Initialisation

    !! 1D and 2D cases
    FM(1)      = 1.0
    dMdEta(1)  = 0.0
    FL(1)      = 1.0
    dLdZeta(1) = 0.0

    !! Intervals of current element
    do i = 1, patch%dim_patch
        !! TODO verify if this information is not already stored in element object
        Ni(i) = patch%Nijk_patch(i, element%current_elem)
    enddo

    !! Initialize matrices
    dxdxi(:,:) = 0.0
    dXidtildexi(:,:) = 0.0
    AJmat(:,:) = 0.0

    !! Compute parametric coordinates form parents elements
    do i = 1, patch%dim_patch
        xi(i)= ((element%Ukv_elem(2,i) - element%Ukv_elem(1,i))*PtGauss(i)      &
            &   +  (element%Ukv_elem(2,i) + element%Ukv_elem(1,i)) ) * 0.5d0
    enddo

    !! Compute univariate B-spline function
    call dersbasisfuns(Ni(1),patch%Jpqr_patch(1),patch%Nkv_patch(1),xi(1),  &
            &       patch%Ukv1_patch(:),FN,dNdxi)
    call dersbasisfuns(Ni(2),patch%Jpqr_patch(2),patch%Nkv_patch(2),xi(2),  &
            &       patch%Ukv2_patch(:),FM,dMdeta)
    if (patch%dim_patch==3) then
        call dersbasisfuns(Ni(3),patch%Jpqr_patch(3),patch%Nkv_patch(3),xi(3),  &
                &       patch%Ukv3_patch(:),FL,dLdZeta)
    endif

    !! Build numerators and denominators
    NumLoc   = 0
    SumTot   = 0.0
    SumXi(:) = 0.0
    do k = 0,patch%Jpqr_patch(3)
        do j = 0,patch%Jpqr_patch(2)
            do i = 0,patch%Jpqr_patch(1)
                NumLoc = NumLoc+1
                R(NumLoc) = FN(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j)     &
                    &       *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)

                SumTot = SumTot + R(NumLoc)

                dRdxi(NumLoc,1) = dNdxi(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j)      &
                            &   *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)

                dRdxi(NumLoc,2)= FN(patch%Jpqr_patch(1)+1-i)*dMdEta(patch%Jpqr_patch(2)+1-j)        &
                            &   *FL(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)

                dRdxi(NumLoc,3) = FN(patch%Jpqr_patch(1)+1-i)*FM(patch%Jpqr_patch(2)+1-j)           &
                            &   *dLdZeta(patch%Jpqr_patch(3)+1-k)*element%weight_elem(NumLoc)
                SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
            enddo
        enddo
    enddo

    !! Divide by denominator to complete definition of fct and deriv.
    SumTot_inv = 1.0/SumTot
    do NumLoc = 1,patch%nnode_patch
        R(NumLoc) = R(NumLoc)*SumTot_inv
        do i = 1, mcrd
            dRdxi(NumLoc,i) = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
        enddo
    enddo


    !! Gradient of mapping from parameter space to physical space
    do NumLoc= 1, patch%nnode_patch
        do Na = 1,MCRD
            do Nb = 1,MCRD
                dxdxi(Na,Nb) = dxdxi(Na,Nb) + COORDS(Na,NumLoc)*dRdxi(NumLoc,Nb)
            enddo
        enddo
    enddo

    !! Gradient of mapping from parent element to parameter space
    do i = 1, patch%dim_patch
        dXidtildexi(i,i) = 0.5d0*( element%Ukv_elem(2,i) - element%Ukv_elem(1,i) )
    enddo

    do Na = 1,MCRD
        do Nb = 1,MCRD
            do Nc = 1,MCRD
               AJmat(Na,Nb)=AJmat(Na,Nb)+dxdxi(Na,Nc)*dXidtildexi(Nc,Nb)
            enddo
        enddo
    enddo


    !! Determination of elementary surface as a function of face number
    select case (k_numface + mcrd*10)
        case(21,22)         !! sides 1 and 2
            det_jac = sqrt(AJmat(1,2)**2.0 + AJmat(2,2)**2.0)
        case(23,24)         !! sides 3 and 4
            det_jac = sqrt(AJmat(1,1)**2.0 + AJmat(2,1)**2.0)
        case(31,32)         !! Faces 1 and 2
            call SurfElem(AJmat(:,2),AJmat(:,3), det_jac)
        case(33,34)         !! Faces 3 and 4
            call SurfElem(AJmat(:,3),AJmat(:,1), det_jac)
        case(35,36)         !! Faces 5 and 6
            call SurfElem(AJmat(:,1),AJmat(:,2), det_jac)
    endselect

    !! Determination of normalized vector
    select case (k_typedload + mcrd*10)
        case(21)                    !! X direction
            Vect(1) = 1.0
            Vect(2) = 0.0
        case(22)                    !! Y direction
            Vect(1) = 0.0
            Vect(2) = 1.0
        case(31)                    !! X direction
            Vect(1) = 1.0
            Vect(2) = 0.0
            Vect(3) = 0.0
        case(32)                    !! Y direction
            Vect(1) = 0.0
            Vect(2) = 1.0
            Vect(3) = 0.0
        case(33)                    !! Z direction
            Vect(1) = 0.0
            Vect(2) = 0.0
            Vect(3) = 1.0
        case(20)                    !! Normal pressure in 2D
            select case (k_numface)
                case(1)                 !! Side 1
                    Vect(1) = AJmat(2,2)/det_jac
                    Vect(2) =-AJmat(1,2)/det_jac
                case(2)                 !! Side 2
                    Vect(1) =-AJmat(2,2)/det_jac
                    Vect(2) = AJmat(1,2)/det_jac
                case(3)                 !! Side 3
                    Vect(1) =-AJmat(2,1)/det_jac
                    Vect(2) = AJmat(1,1)/det_jac
                case(4)                 !! Side 4
                    Vect(1) = AJmat(2,1)/det_jac
                    Vect(2) =-AJmat(1,1)/det_jac
            end select
        case(29)                    !! Tangent pressure in 2D
            select case (k_numface)
            case(1)                     !! Side 1
                Vect(1)= AJmat(2,2)/det_jac
                Vect(2)= AJmat(1,2)/det_jac
            case(2)                     !! Side 2
                Vect(1)=-AJmat(2,2)/det_jac
                Vect(2)=-AJmat(1,2)/det_jac
            case(3)                     !! Side 3
                Vect(1)= AJmat(2,1)/det_jac
                Vect(2)= AJmat(1,1)/det_jac
            case(4)                     !! Side 4
                Vect(1)=-AJmat(2,1)/det_jac
                Vect(2)=-AJmat(1,1)/det_jac
            end select
        case(30)                    !! Normal pressure in 3D
            select case (k_numface)
                case(1)                 !! Face 1
                    call VectNormNorm(Vect,AJmat(:,2),AJmat(:,3),det_jac)
                case(2)                 !! Face 2
                    call VectNormNorm(Vect,AJmat(:,3),AJmat(:,2),det_jac)
                case(3)                 !! Face 3
                    call VectNormNorm(Vect,AJmat(:,3),AJmat(:,1),det_jac)
                case(4)                 !! Face 4
                    call VectNormNorm(Vect,AJmat(:,1),AJmat(:,3),det_jac)
                case(5)                 !! Face 5
                    call VectNormNorm(Vect,AJmat(:,1),AJmat(:,2),det_jac)
                case(6)                 !! Face 6
                    call VectNormNorm(Vect,AJmat(:,2),AJmat(:,1),det_jac)
            end select
        case(34)                    !! Non-uniform normal pressure
            select case (k_numface)
                case(1)                 !! Face 1
                    call VectNormNorm(Vect,AJmat(:,2),AJmat(:,3),det_jac)
                case(2)                 !! Face 2
                    call VectNormNorm(Vect,AJmat(:,3),AJmat(:,2),det_jac)
                case(3)                 !! Face 3
                    call VectNormNorm(Vect,AJmat(:,3),AJmat(:,1),det_jac)
                case(4)                 !! Face 4
                    call VectNormNorm(Vect,AJmat(:,1),AJmat(:,3),det_jac)
                case(5)                 !! Face 5
                    call VectNormNorm(Vect,AJmat(:,1),AJmat(:,2),det_jac)
                case(6)                 !! Face 6
                    call VectNormNorm(Vect,AJmat(:,2),AJmat(:,1),det_jac)
            end select
     end select

end subroutine shapPress_omp



