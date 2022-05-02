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
      
      Subroutine apply_dispBC(K,F,nb_cp,MCRD,nb_bc,bc_target,
     1     bc_target_nbelem,bc_values,nb_dof_tot,U)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_cp, MCRD, nb_bc, bc_target,
     &     bc_target_nbelem, nb_dof_tot
      dimension bc_target(SUM(bc_target_nbelem)),bc_target_nbelem(nb_bc)
      
      Double precision, intent(in) :: bc_values
      dimension bc_values(2,nb_bc)
      
c     Output variables :
c     ----------------
      Double precision, intent(inout) :: K,F,U
      dimension K(nb_cp*MCRD,nb_cp*MCRD),F(nb_cp*MCRD),U(nb_cp*MCRD)
      
c     Local variables :
c     ---------------
      Integer :: i,j,k1,sctr,node,count
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Applications des conditions aux limites ..........................
      
!     Loop on boundary conditions
      count = 0
      Do i  = 1,nb_bc
         If (bc_values(i,1) <= 3) then
            Do j = 1,bc_target_nbelem(i)
               node = bc_target(count+j)
               sctr = (node-1)*MCRD
      End subroutine apply_dispBC
