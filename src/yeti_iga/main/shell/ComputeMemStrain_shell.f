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
      
      Subroutine uStrainMem_shell(sol,NNODE,MCRD,AI,dRdxi,StrainMem)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD
      Double precision, intent(in) :: sol, AI, dRdxi
      dimension sol(MCRD,NNODE), AI(3,3), dRdxi(NNODE,2)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: StrainMem
      dimension StrainMem(MCRD)
      
c     Local variables :
c     ---------------
      Double precision :: BoJ, dNjdxi, StrainMem_loc
      dimension BoJ(MCRD,MCRD), dNjdxi(2), StrainMem_loc(MCRD)
      
      Integer :: nodj, jdim
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Calcul deformation membrane ......................................
c     
c     Initialisation
      StrainMem(:) = 0.0D0

c     Boucle element
      do nodj = 1,NNODE
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
         enddo  
         
c     Calcul matrice BoJ
         BoJ(:,:) = 0.0D0
         BoJ(1,:) = dNjdxi(1)*AI(1,:)
         BoJ(2,:) = dNjdxi(2)*AI(2,:)
         BoJ(3,:) = dNjdxi(2)*AI(1,:) + dNjdxi(1)*AI(2,:)
         
c     Mise a jour champ deformations
         call MulVect(BoJ, sol(:,nodj), StrainMem_loc, MCRD, MCRD)
         StrainMem(:) = StrainMem(:) + StrainMem_loc(:)
         
      enddo
      
      End subroutine  uStrainMem_shell
      
