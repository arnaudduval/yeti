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

c     Include iga subroutines
c     include "./operateurs.f"
      
      
C     ******************************************************************
      
      Subroutine reconstruction(sol,U,U_inv,COORDS,ind_dof_free,
     1     bc_target,bc_target_nbelem,bc_values,COUPLG_flag,nb_bc,MCRD,
     2     nb_cp,nb_dof_free)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_cp,MCRD,nb_dof_free,ind_dof_free,nb_bc,
     &     bc_target,bc_target_nbelem
      dimension ind_dof_free(nb_cp*MCRD), bc_target(:),
     &     bc_target_nbelem(nb_bc)
      
      Double precision, intent(in) :: U_inv,COORDS,bc_values
      dimension U_inv(nb_dof_free), COORDS(3,nb_cp), bc_values(2,nb_bc)
      
      Logical, intent(in) :: COUPLG_flag
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: sol,U
      dimension sol(nb_cp,MCRD), U(nb_cp*MCRD)
      
c     Local variables :
c     ---------------
      Integer :: ddl,i,num_bc,num_cp,direction, ddl_slave,
     &     ddl_master1, ddl_master2, num_cp_slave, num_cp_master1,
     &     num_cp_master2, count
      Double precision :: vectu, vectv, norm1, norm2, c
      Dimension vectu(3), vectv(3)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Application des couplages ........................................
      
c     Reconstruction disp pour ddl maitre
      Do ddl = 1,nb_dof_free
         U(ind_dof_free(ddl)) = U_inv(ddl)
      Enddo
      
      count = 0
      Do num_bc = 1,nb_bc
         If (bc_values(1,num_bc) <= 3) then
            Do num_cp = 1,bc_target_nbelem(num_bc)
               ddl = (bc_target(count+num_cp) - 1)*MCRD
     &              + int(bc_values(1,num_bc))
               U(ddl) = bc_values(2,num_bc)
            Enddo
         Endif
         count = count+bc_target_nbelem(num_bc)
      Enddo
      
c     Reconstruction disp pour ddl slave
c      If (COUPLG_flag) then
c         Do num_bc = 1,nb_bc
c            If (6.<bc_values(num_bc,1) .and.
c     &           bc_values(num_bc,1)<10.) then
c               direction = int(bc_values(num_bc,1)) - 6
c               Do num_cp = 1,bc_target_nbelem(num_bc)
c                  ddl_slave = (bc_target(num_bc, num_cp) - 1)*MCRD
c     &                 + direction
c                  ddl_master1 = (bc_target(num_bc+200, num_cp) - 1)*MCRD
c     &                 + direction
c                  
c                  U(ddl_slave) = U(ddl_master1) 
c               Enddo
c            
c            
c            Elseif (bc_values(num_bc,1) == 10.) then
c               Do num_cp = 1,bc_target_nbelem(num_bc)
c                  
c                  num_cp_slave   = bc_target(num_bc, num_cp)
c                  num_cp_master1 = bc_target(num_bc+300, num_cp)
c                  num_cp_master2 = bc_target(num_bc+400, num_cp)
c               
c                  vectu = COORDS(num_cp_slave,:)
c     &                 - COORDS(num_cp_master1,:)
c                  vectv = COORDS(num_cp_master1,:)
c     &                 - COORDS(num_cp_master2,:)
!                 call dot(vectu,vectv, norm1)
!                 call dot(vectu,vectu, norm2)
c                  norm1=vectu(1)**2.D0 + vectu(2)**2.D0 + vectu(3)**2.D0
c                  norm2=vectv(1)**2.D0 + vectv(2)**2.D0 + vectv(3)**2.D0
c                  
c                  c = norm2/norm1
c                  
c                  Do direction = 1,3
c                     ddl_slave   = (num_cp_slave   - 1)*MCRD + direction
c                     ddl_master1 = (num_cp_master1 - 1)*MCRD + direction
c                     ddl_master2 = (num_cp_master2 - 1)*MCRD + direction
c                     
c                     U(ddl_slave) = (c+1.0D0) * U(ddl_master1) 
c     &                    - c*U(ddl_master2)
c                  Enddo
c               Enddo
c            Endif
c
c         Enddo
c      Endif
      
c     Construction solution finale
      Do i = 1,MCRD
         Do num_cp = 1,nb_cp
            sol(num_cp,i) = U( (num_cp-1)*MCRD + i )
         Enddo
      Enddo
      
      End subroutine reconstruction
