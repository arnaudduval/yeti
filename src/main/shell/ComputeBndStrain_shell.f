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
      
C     ******************************************************************
      
      Subroutine uStrainBnd_shell(sol,NNODE,MCRD,AI,dAI1dxi,dAI2dxi,
     1     dRdxi,ddRddxi,StrainBnd)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD
      Double precision, intent(in) :: sol, AI, dAI1dxi, dAI2dxi, dRdxi,
     &     ddRddxi
      dimension sol(MCRD,NNODE), AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3),
     &     dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: StrainBnd
      dimension StrainBnd(MCRD)
      
c     Local variables :
c     ---------------
      Double precision :: BoJ, dNjdxi, ddNjddxi, B1, B2, B3, Area,
     &     dA1d1_A2, dA2d2_A2, dA1d2_A2, A1_dA1d1, A1_dA2d2, A1_dA1d2,
     &     A2_A3, A3_A1, A3dA1d1, A3dA2d2, A3dA1d2, StrainBnd_loc
      dimension BoJ(MCRD,MCRD), dNjdxi(2), ddNjddxi(3), B1(3), B2(3),
     &     B3(3), dA1d1_A2(3), dA2d2_A2(3), dA1d2_A2(3), A1_dA1d1(3),
     &     A1_dA2d2(3), A1_dA1d2(3), A2_A3(3), A3_A1(3),
     &     StrainBnd_loc(MCRD)
      
      Integer :: nodj, jdim
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Calcul deformation flexion .......................................
c     
c     Initialisation
      StrainBnd(:) = 0.0D0
      
c     Calculs preliminaire
      call SurfElem(AI(1,:), AI(2,:), Area)
      call cross(dAI1dxi(1,:), AI(2,:), dA1d1_A2(:))
      call cross(dAI2dxi(2,:), AI(2,:), dA2d2_A2(:))
      call cross(dAI1dxi(2,:), AI(2,:), dA1d2_A2(:))
      call cross(AI(1,:), dAI1dxi(1,:), A1_dA1d1(:))
      call cross(AI(1,:), dAI2dxi(2,:), A1_dA2d2(:))
      call cross(AI(1,:), dAI1dxi(2,:), A1_dA1d2(:))
      call cross(AI(2,:), AI(3,:), A2_A3(:))
      call cross(AI(3,:), AI(1,:), A3_A1(:))
      call dot(AI(3,:), dAI1dxi(1,:), A3dA1d1)
      call dot(AI(3,:), dAI2dxi(2,:), A3dA2d2)
      call dot(AI(3,:), dAI1dxi(2,:), A3dA1d2)
      
      
c     Boucle element
      do nodj = 1,NNODE
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
            ddNjddxi(jdim) = ddRddxi(nodj,jdim)
         enddo
         ddNjddxi(3) = ddRddxi(nodj,3)
         
c     Calcul matrice BoJ

c     Calcul matrice BoJ
         B1(:) = 0.0d0; B2(:) = 0.0d0; B3(:) = 0.0d0
         B1(:) = -ddNjddxi(1)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d1_A2(:) + dNjdxi(2)*A1_dA1d1(:)
     &        + A3dA1d1*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B2(:) = -ddNjddxi(2)*AI(3,:) +
     &        ( dNjdxi(1)*dA2d2_A2(:) + dNjdxi(2)*A1_dA2d2(:)
     &        + A3dA2d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B3(:) = -ddNjddxi(3)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d2_A2(:) + dNjdxi(2)*A1_dA1d2(:)
     &        + A3dA1d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         BoJ(1,:) = B1(:)
         BoJ(2,:) = B2(:)
         BoJ(3,:) = 2.0D0*B3(:)
         
c     Mise a jour champ deformations
         call MulVect(BoJ, sol(:,nodj), StrainBnd_loc, MCRD, MCRD)
         StrainBnd(:) = StrainBnd(:) + StrainBnd_loc(:)
         
      enddo
      
      End subroutine  uStrainBnd_shell
      
