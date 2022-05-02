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
C     Calcul des fonctions de base B-spline pour le cas des surfaces
c     uniquement !
c
c     Retourne toutes les fonctions non nulles pour le point de gauss 
c     considere, ainsi que les derivees premieres et secondes.
c     --
      
C     ******************************************************************
      
      subroutine bsplinebasis(R,dRdxi,ddRddxi,DetJac,NNODE,PtGaus,JELEM,
     1     NumPatch)
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, JELEM, NumPatch
      
      Double precision, intent(in) :: PtGaus
      dimension PtGaus(2)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: R, dRdxi, ddRddxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
      
c     Local variables :
c     ---------------
      
!     Parameters and COMMON variables
      Double precision :: zero, one, two
      parameter (zero=0.0D0, one=1.0D0, two=2.0D0)
      Integer, parameter :: Lelement=10000, LNODE=10000, Lknot=10000,
     &     Lpatch=100
      
      Common /BndStripParameter/ Ukv_bs, Nkv_bs, Jpqr_bs, Nijk_bs
      Double precision :: Ukv_bs
      Integer :: Nkv_bs, Jpqr_bs, Nijk_bs
      dimension Ukv_bs(Lknot,3,Lpatch), Nkv_bs(3,Lpatch), Jpqr_bs(3),
     &     Nijk_bs(Lelement,3)
      
!     Nurbs functions and derivatives
      Double precision :: xi, FN, FM, dNdxi, dMdEta, ddNddxi, ddMddEta,
     &     dXidtildexi
      dimension xi(2), FN(Jpqr_bs(1)+1), FM(Jpqr_bs(2)+1),
     &     dNdxi(Jpqr_bs(1)+1),dMdEta(Jpqr_bs(2)+1), 
     &     ddNddxi(Jpqr_bs(1)+1),ddMddEta(Jpqr_bs(2)+1),dXidtildexi(2,2)
      
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
         Ni(i) = Nijk_bs(JELEM,i)
      enddo
      
!     Initialisation des differentes matrices
      dXidtildexi = zero
      R = zero
      dRdxi = zero
      ddxddxi = zero
      
c     Fin Initialisation ...............................................
c
c
c
c
c     Debut Analyse ....................................................
      
c     Calculate parametrique coordinates from parents element
      do i = 1,2
         xi(i) = 
     &        ((Ukv_bs(Ni(i)+1,i,NumPatch)-Ukv_bs(Ni(i),i,NumPatch))
     &        * PtGaus(i)
     &        + (Ukv_bs(Ni(i)+1,i,NumPatch)+Ukv_bs(Ni(i),i,NumPatch)))
     &        / two
      enddo
      
c     Calculate univariate B-spline function
      call dersbasisfuns2(Ni(1),Jpqr_bs(1),Nkv_bs(1,NumPatch),xi(1),
     &     Ukv_bs(:Nkv_bs(1,NumPatch),1,NumPatch),FN,dNdxi,ddNddxi)
      call dersbasisfuns2(Ni(2),Jpqr_bs(2),NKv_bs(2,NumPatch),xi(2),
     &     Ukv_bs(:NKv_bs(2,NumPatch),2,NumPatch),FM,dMdEta,ddMddEta)
      
      
c     Build numerators and denominators
      NumLoc = 0
      do j = 0,Jpqr_bs(2)
         do i = 0,Jpqr_bs(1)
            NumLoc = NumLoc+1
!     B-spline functions
            R(NumLoc) = FN(Jpqr_bs(1)+1-i)*FM(Jpqr_bs(2)+1-j)
            
!     B-spline first derivatives
            dRdxi(NumLoc,1) = dNdxi(Jpqr_bs(1)+1-i)*FM(Jpqr_bs(2)+1-j)
            dRdxi(NumLoc,2) = FN(Jpqr_bs(1)+1-i)*dMdEta(Jpqr_bs(2)+1-j)
            
!     B-spline second derivatives
            ddRddxi(NumLoc,1) = 
     &           ddNddxi(Jpqr_bs(1)+1-i)*FM(Jpqr_bs(2)+1-j)
            ddRddxi(NumLoc,2) = 
     &           FN(Jpqr_bs(1)+1-i)*ddMddEta(Jpqr_bs(2)+1-j)
            ddRddxi(NumLoc,3) = 
     &           dNdxi(Jpqr_bs(1)+1-i)*dMdEta(Jpqr_bs(2)+1-j)
         enddo
      enddo
      
      
c     Gradient of mapping from parent element to parameter space
      do i = 1,2  
         dXidtildexi(i,i) = 
     &        (Ukv_bs(Ni(i)+1,i,NumPatch)-Ukv_bs(Ni(i),i,NumPatch))/two
      enddo
      DetJac = dXidtildexi(1,1)*dXidtildexi(2,2)
      
      
      return
      
      End Subroutine bsplinebasis
