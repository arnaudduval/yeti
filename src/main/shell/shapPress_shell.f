!! Copyright 2017-2019 Thibaut Hirschler

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

c     --
C     Calcul des quantites necessaire pour l'actualisation au point de
c     gauss du second membre pour les coques types KL
c
c     Retourne toutes les fonctions non nulles pour le point de gauss 
c     considere, le vecteur directeur de l'effort, ainsi que le jacobien
c     pour l'integration
c     --
      
C     ******************************************************************
      
      subroutine shapPress3(Vect,R,DetJac, KTypeDload,KNumFace,COORDS,
     1     PtGauss,MCRD)
      
      use parameters
      use nurbspatch      
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: KTypeDload,KNumFace,MCRD
      Double precision, intent(in) :: PtGauss,COORDS
      dimension PtGauss(2), COORDS(MCRD,nnode_patch)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: Vect,R,DetJac
      dimension Vect(MCRD),R(nnode_patch)
      
      
c     Local variables :
c     ---------------
      
!     Nurbs functions and derivatives
      Double precision :: xi, FN, FM, dNdxi, dMdEta, SumTot, SumXi, 
     &     dXidtildexi, dRdxi(nnode_patch,2)
      dimension xi(2), FN(Jpqr_patch(1)+1), FM(Jpqr_patch(2)+1), 
     &     dNdxi(Jpqr_patch(1)+1),dMdEta(Jpqr_patch(2)+1), SumXi(2),
     &     dXidtildexi(2)
      
!     Curvilinear quantities
      Double precision :: A1,A2,A3, norm1,norm2,normV
      dimension A1(MCRD),A2(MCRD),A3(MCRD)
      
      Integer :: i,j,Ni,NumLoc
      dimension Ni(2)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Extraction des Ni,Nj,et Nk de la matrice Nijk(Nknock,3)       
      do i = 1,2
         Ni(i) = Nijk_patch(i,current_elem)
      enddo
      
!     Initialisation des differentes matrices
      dXidtildexi = zero
      R = zero
      dRdxi = zero
      
c     Fin Initialisation ...............................................
c
c
c
c
c     Debut Analyse ....................................................
      
c     Calculate parametrique coordinates from parents element
      do i = 1,dim_patch
         xi(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &        +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
      enddo
      
c     Calculate univariate B-spline function
      CALL dersbasisfuns(Ni(1),Jpqr_patch(1),Nkv_patch(1),xi(1),
     &     Ukv1_patch(:),FN,dNdxi)
      CALL dersbasisfuns(Ni(2),Jpqr_patch(2),Nkv_patch(2),xi(2),
     &     Ukv2_patch(:),FM,dMdeta)
      
      
c     Build numerators and denominators
      NumLoc = 0
      SumTot = zero
      SumXi  = zero
      do j = 0,Jpqr_patch(2)
         do i = 0,Jpqr_patch(1)
            NumLoc = NumLoc+1
!     Nurbs functions
            R(NumLoc) = FN(Jpqr_patch(1)+1-i)*FM(Jpqr_patch(2)+1-j)
     &           *Weight_elem(NumLoc)
            SumTot = SumTot + R(NumLoc)
            
!     Nurbs first derivatives
            dRdxi(NumLoc,1) = dNdxi(Jpqr_patch(1)+1-i)
     &           *FM(Jpqr_patch(2)+1-j)*Weight_elem(NumLoc)
            SumXi(1) = SumXi(1)+dRdxi(NumLoc,1)
c     
            dRdxi(NumLoc,2) = FN(Jpqr_patch(1)+1-i)
     &           *dMdEta(Jpqr_patch(2)+1-j)*Weight_elem(NumLoc)
            SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)
         enddo
      enddo
      
      
!     Complete definition of function and 1st derivatives
      do NumLoc = 1,nnode_patch
         R(NumLoc) = R(NumLoc)/SumTot
         do i = 1,2
            dRdxi(NumLoc,i) = 
     &           (dRdxi(NumLoc,i) - R(NumLoc)*SumXi(i)) / SumTot
         enddo
      enddo
            
      
c     Gradient of mapping from parent element to parameter space
      do i = 1,2  
         dXidtildexi(i) = 0.5d0*( Ukv_elem(2,i) - Ukv_elem(1,i) )
      enddo
      ! DetJac = dXidtildexi(1,1)*dXidtildexi(2,2)
      
      
c     Covariant basis vectors
      A1(:) = zero
      A2(:) = zero
      Do numLoc = 1,nnode_patch
         A1(:) = A1(:) + dRdxi(numLoc,1)*COORDS(:,numLoc)
         A2(:) = A2(:) + dRdxi(numLoc,2)*COORDS(:,numLoc)
      Enddo
      call cross(A1(:), A2(:), A3(:))
      call norm(A3(:), 3, normV)
      A3(:) = A3(:)/normV

c     Compute jacobian from parent element to physical space
      SELECT CASE (KNumFace)
      case(1,2)
         norm2  = sqrt(A2(1)*A2(1) + A2(2)*A2(2) + A2(3)*A2(3))
         DetJac = dXidtildexi(2)*norm2
      case(3,4)
         norm1  = sqrt(A1(1)*A1(1) + A1(2)*A1(2) + A1(3)*A1(3))
         DetJac = dXidtildexi(1)*norm1
      END SELECT
      
      
c     Compute load direction
      SELECT CASE (KTypeDload)
!     -
!     Pression normale a la surface chargee
      case(0)
         SELECT CASE (KNumFace)
         case(1)
            call cross(A2(:), A3(:), Vect(:))
         case(2)
            call cross(A3(:), A2(:), Vect(:))
         case(3)
            call cross(A3(:), A1(:), Vect(:))
         case(4)
            call cross(A1(:), A3(:), Vect(:))
         END SELECT
!     -
!     Pression tangentielle dans le plan median
      case(9)
         SELECT CASE (KNumFace)
         case(1,2)
            Vect(:) = A2(:)
         case(3,4)
            Vect(:) = A1(:)
         END SELECT
!     -
!     Pression tangentielle hors plan
      case(8)
         Vect(:) = A3(:)
!     -
!     Pression dans direction x, y, ou z
      case(1:3)
         Vect(:) = zero
         Vect(KTypeDload) = one
      END SELECT
      
c     Normalisation si necessaire
      SELECT CASE (KTypeDload)
      case(0,9,8)
         normV = Vect(1)*Vect(1) + Vect(2)*Vect(2) + Vect(3)*Vect(3)
         Vect(:) = Vect(:)/sqrt(normV)
      END SELECT
      
      
      
      End Subroutine shapPress3
