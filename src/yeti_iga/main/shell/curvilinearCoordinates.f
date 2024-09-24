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

C     --
c     Calcul des quantites necessaires a la construction du systeme
c     lineaire pour les coques en comportment elastique lineaire.
c     
c     Retourne les vecteurs covariants pour le point de gauss considere,
c     leurs derivees, et les coefficients metriques contravariants.
c     D'autres quantites sont calcules et peuvent etre choisies d'etre
c     retournees si besoin.
c     --
      
      
C     ******************************************************************
      
      subroutine curvilinear(AI, dAI1dxi, dAI2dxi, AAE,
     1     R,dRdxi,ddRddxi,MCRD,NNODE,COORDS)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: MCRD, NNODE
      Double precision, intent(in) :: R, dRdxi, ddRddxi, COORDS
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3),
     &     COORDS(MCRD,NNODE)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: AI, dAI1dxi, dAI2dxi, AAE
      dimension AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3), AAE(2,2)
      
      
c     Local variables :
c     ---------------
      
!     Parameters and COMMON variables
      Double precision :: zero, one, two
      parameter (zero=0.0D0, one=1.0D0, two=2.0D0)
      
!     Curvilinear Coordinates
      Double precision :: AAI, vect, scal, det, norm
      dimension AAI(2,2), vect(3)
      
      Integer :: i,j,numLoc
      
C     Fin declaration des variables ....................................
      
      
      
      
c     Initialisation
      AI = zero
      AAI = zero
      AAE = zero
      dAI1dxi = zero
      dAI2dxi = zero
      
c     Covariant basis vectors
      Do numLoc = 1,NNODE
         AI(1,:) = AI(1,:) + dRdxi(numLoc,1)*COORDS(:,numLoc)
         AI(2,:) = AI(2,:) + dRdxi(numLoc,2)*COORDS(:,numLoc)
      Enddo
      call cross(AI(1,:), AI(2,:), vect)
      norm = sqrt( vect(1)**two + vect(2)**two + vect(3)**two  )
      AI(3,:) = vect(:)/norm
      
c     Derivatives of Covariant basis vectors
      Do numLoc = 1,NNODE
         dAI1dxi(1,:) = dAI1dxi(1,:)+ddRddxi(numLoc,1)*COORDS(:,numLoc)
         dAI1dxi(2,:) = dAI1dxi(2,:)+ddRddxi(numLoc,3)*COORDS(:,numLoc)
         dAI2dxi(1,:) = dAI2dxi(1,:)+ddRddxi(numLoc,3)*COORDS(:,numLoc)
         dAI2dxi(2,:) = dAI2dxi(2,:)+ddRddxi(numLoc,2)*COORDS(:,numLoc)
      Enddo
      

c     Covariant Mertic coefficients
      AAI = zero
      Do i = 1,2
         Do j = 1,2
            call dot(AI(i,:), AI(j,:), scal)
            AAI(i,j) = scal
         Enddo
      Enddo
      
c     Contravariant Mertic coefficients
      call MatrixInv(AAE, AAI, det, 2)
      
      
      End subroutine curvilinear
