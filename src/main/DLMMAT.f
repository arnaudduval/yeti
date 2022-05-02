!! Copyright 2016-2018 Thibaut Hirschler

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

!>    Routine calculant la "Diagonally-Lumped Mapping Matrix" elementaire
!!    Retourne un vecteur comprenant les termes diagonaux
      

C     ******************************************************************

      Subroutine DLMMAT(MCRD,NNODE,JELEM,NBINT,NumPatch,COORDS,ELT_TYPE,
     1     DLMMATRX)
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: MCRD, NNODE, JELEM, NumPatch, NBINT
      Double precision, intent(in) :: COORDS
      dimension COORDS(MCRD,NNODE)
      
      Character(len=*), intent(in) :: ELT_TYPE
      
      
c     Output variables : matrice de masse elementaire
c     ----------------
      Double precision, intent(out) :: DLMMATRX
      dimension DLMMATRX(NNODE)
      
      
      
c     Local variables :
c     ---------------
      
      Integer :: n,k,numLoc
      Double precision :: dvol
      
!     Parameters
      Double precision zero, one, two
      parameter (zero=0.d0, one=1.d0,two=2.d0)
      
!     for gauss.f
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(MCRD+1,NBINT)
      
!     for shap.f
      Integer :: NDOFEL
      Double precision :: R, dRdx, DetJac, PtGauss
      dimension R(NNODE), dRdx(MCRD,NNODE), PtGauss(MCRD)

!     for nurbsbasis.f
      Double precision :: A1, A2, dRdxi, ddRddxi, Area
      dimension A1(MCRD), A2(MCRD), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Initialize DLMMATRX to zero
      DLMMATRX(:) = zero

!     Number of integration point in each direction
      if (MCRD==2 .OR. ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         NbPtInt = int(NBINT**(1.0/2.0))
      else if (MCRD==3) then
         NbPtInt = int(NBINT**(1.0/3.0))
      endif
      
!     Number of DOF by element
      NDOFEL = MCRD*NNODE

!     Defining Gauss points coordinates and Gauss weights
      If (ELT_TYPE=='U1') then
         call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)
      elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         call Gauss(NbPtInt, 2, GaussPdsCoord(:MCRD,:), 0)
      Endif
      
c     Fin initialisaton ................................................
c     
c     
c     
c     Fin calcul .......................................................
      
!     Loop on integration points
      Do n = 1,NBINT
!     Computing Nurbs basis functions at integration points
         If (ELT_TYPE=='U1') then
            call shap(COORDS,dRdx,NNODE,R,GaussPdsCoord(2:,n),DetJac,
     &           NDOFEL,MCRD,JELEM,NumPatch)
            dvol = GaussPdsCoord(1,n)*DetJac
            
         elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
            call nurbsbasis(R,dRdxi,ddRddxi,DetJac,NNODE,
     &           GaussPdsCoord(2:MCRD,n),JELEM,NumPatch)
            
            A1(:) = zero; A2(:) = zero
            Do numLoc = 1,NNODE
               A1(:) = A1(:) + dRdxi(numLoc,1)*COORDS(:,numLoc)
               A2(:) = A2(:) + dRdxi(numLoc,2)*COORDS(:,numLoc)
            Enddo
            call SurfElem(A1(:), A2(:), Area)
            dvol = GaussPdsCoord(1,n)*DetJac*Area
            
         Endif
         
         
!     Assembling DLMMATRIX
         Do k = 1,NNODE
            DLMMATRX(k) = DLMMATRX(k) + R(k)*dvol
         Enddo  
      Enddo
!     End loop on integration points
      
c     Fin calcul .......................................................
      
      End Subroutine DLMMAT
