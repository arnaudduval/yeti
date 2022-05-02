!! Copyright 2017-2018 Thibaut Hirschler

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
C     Calcul des fonctions de base Nurbs pour le cas des surfaces
c     uniquement !
c
c     Retourne toutes les fonctions non nulles pour le point de gauss 
c     considere, ainsi que les derivees premieres et secondes.
c     --
      
C     ******************************************************************
      
      subroutine nurbsbasis(R,dRdxi,ddRddxi,DetJac,PtGauss)
      
      use parameters
      use nurbspatch
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------      
      Double precision, intent(in) :: PtGauss
      dimension PtGauss(2)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: R, dRdxi, ddRddxi, DetJac
      dimension R(nnode_patch),dRdxi(nnode_patch,2),
     &     ddRddxi(nnode_patch,3)
      
      
c     Local variables :
c     ---------------
      
!     Nurbs functions and derivatives
      Double precision :: xi, FN, FM, dNdxi, dMdEta, ddNddxi, ddMddEta,
     &     SumTot, SumXi, SumXixi, dXidtildexi
      dimension xi(2), FN(Jpqr_patch(1)+1), FM(Jpqr_patch(2)+1), 
     &     dNdxi(Jpqr_patch(1)+1), dMdEta(Jpqr_patch(2)+1),
     &     ddNddxi(Jpqr_patch(1)+1), ddMddEta(Jpqr_patch(2)+1),
     &     SumXi(2), SumXixi(3), dXidtildexi(2,2)
      
      Integer :: i, j, Ni, NumLoc
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
      dXidtildexi(:,:) = zero
      R(:)         = zero
      dRdxi(:,:)   = zero
      ddRddxi(:,:) = zero
      
c     Fin Initialisation ...............................................
c
c
c
c
c     Debut Analyse ....................................................
      
c     Calculate parametrique coordinates from parents element
      Do i = 1,dim_patch
         xi(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &        +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
      End do
      
c     Calculate univariate B-spline function
      CALL dersbasisfuns2(Ni(1),Jpqr_patch(1),Nkv_patch(1),xi(1),
     &     Ukv1_patch(:),FN,dNdxi,ddNddxi)
      CALL dersbasisfuns2(Ni(2),Jpqr_patch(2),Nkv_patch(2),xi(2),
     &     Ukv2_patch(:),FM,dMdeta,ddMddEta)
      
      
c     Build numerators and denominators
      NumLoc = 0
      SumTot = zero
      SumXi  = zero
      SumXixi= zero
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
            
!     Nurbs second derivatives
            ddRddxi(NumLoc,1) = ddNddxi(Jpqr_patch(1)+1-i)
     &           *FM(Jpqr_patch(2)+1-j)*Weight_elem(NumLoc)
            SumXixi(1) = SumXixi(1) + ddRddxi(NumLoc,1)
c     
            ddRddxi(NumLoc,2) = FN(Jpqr_patch(1)+1-i)
     &           *ddMddEta(Jpqr_patch(2)+1-j)*Weight_elem(NumLoc)
            SumXixi(2) = SumXixi(2) + ddRddxi(NumLoc,2)
c     
            ddRddxi(NumLoc,3) = dNdxi(Jpqr_patch(1)+1-i)
     &           *dMdEta(Jpqr_patch(2)+1-j)*Weight_elem(NumLoc)
            SumXixi(3) = SumXixi(3) + ddRddxi(NumLoc,3)
         enddo
      enddo
      
      
!     Divide by denominator to complete definition of function and 1st derivatives
      do NumLoc = 1,nnode_patch
         R(NumLoc) = R(NumLoc)/SumTot
         do i = 1,2
            dRdxi(NumLoc,i) = 
     &           (dRdxi(NumLoc,i) - R(NumLoc)*SumXi(i)) / SumTot
         enddo
      enddo
      
!     Complete definition of 2nd derivatives
      do NumLoc = 1,nnode_patch
         ddRddxi(NumLoc,1) = ( ddRddxi(NumLoc,1) - R(NumLoc)*SumXixi(1)
     &        - two*dRdxi(NumLoc,1)*SumXi(1)
     &        )/SumTot
         ddRddxi(NumLoc,2) = ( ddRddxi(NumLoc,2) - R(NumLoc)*SumXixi(2)
     &        - two*dRdxi(NumLoc,2)*SumXi(2)
     &        )/SumTot
         ddRddxi(NumLoc,3) = ( ddRddxi(NumLoc,3) - R(NumLoc)*SumXixi(3)
     &        - dRdxi(NumLoc,1)*SumXi(2) - dRdxi(NumLoc,2)*SumXi(1)
     &        )/SumTot
      enddo
      
      
      
c     Gradient of mapping from parent element to parameter space
      do i = 1,2
         dXidtildexi(i,i) = 0.5d0*( Ukv_elem(2,i) - Ukv_elem(1,i) )
      enddo
      DetJac = dXidtildexi(1,1)*dXidtildexi(2,2)
      
      
      End Subroutine nurbsbasis
