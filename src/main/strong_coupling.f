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
      
      subroutine coupling(K,F,nb_cp,MCRD,nb_bc,bc_target,
     1     bc_target_nbelem,bc_values,COORDS)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_cp,MCRD,nb_bc,bc_target,
     &     bc_target_nbelem
      dimension bc_target(1000,1000), bc_target_nbelem(1000)
      
      Double precision, intent(in) :: bc_values, COORDS
      dimension bc_values(1000,2), COORDS(nb_cp,3)
      
c     Output variables :
c     ----------------
      Double precision, intent(inout) :: K,F
      dimension K(nb_cp*MCRD,nb_cp*MCRD),F(nb_cp*MCRD)
      
c     Local variables :
c     ---------------
      Integer :: nb_dof_tot, num_bc, num_cp, direction, ddl_slave,
     &     ddl_master1, ddl_master2, num_cp_slave, num_cp_master1,
     &     num_cp_master2
      Double precision :: vectu, vectv, norm1, norm2, c
      Dimension vectu(3), vectv(3)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Application des couplages ........................................
      
c     Initialisation
      nb_dof_tot = MCRD*nb_cp
c     Construction matrice de couplage
      Do num_bc = 1,nb_bc
         
         If (6.<bc_values(num_bc,1) .and. bc_values(num_bc,1)<10.) then
!     Cas condition egalite maitre/escave (pour symmetrie par ex.)
            direction = int(bc_values(num_bc,1)) - 6
            Do num_cp = 1,bc_target_nbelem(num_bc)
               
               ddl_slave = (bc_target(num_bc, num_cp) - 1)*MCRD
     &              + direction
               ddl_master1 = (bc_target(num_bc+200, num_cp) - 1)*MCRD
     &              + direction
               
c     Modification second membre
               F(ddl_master1) = F(ddl_master1) + F(ddl_slave)
               F(ddl_slave)   = 0.0D0
               
c     Modification matrice de rigidite
               K(:,ddl_master1) = K(:,ddl_master1) + K(:,ddl_slave)
               K(:,ddl_slave)   = 0.0D0
               
               K(ddl_master1,:) = K(ddl_master1,:) + K(ddl_slave,:)
               K(ddl_slave,:)   = 0.0D0
            Enddo
            
            
            
            
         Elseif (bc_values(num_bc,1) == 10.) then
!     Cas condition de couplage fort entre patches tangents
            Do num_cp = 1,bc_target_nbelem(num_bc)
               
               num_cp_slave   = bc_target(num_bc, num_cp)
               num_cp_master1 = bc_target(num_bc+300, num_cp)
               num_cp_master2 = bc_target(num_bc+400, num_cp)
               
               vectu = COORDS(num_cp_slave,:)-COORDS(num_cp_master1,:)
               vectv = COORDS(num_cp_master1,:)-COORDS(num_cp_master2,:)
               call dot(vectu,vectv, norm1)
               call dot(vectu,vectu, norm2)
               c = norm2/norm1
               
               Do direction = 1,3
                  ddl_slave   = (num_cp_slave   - 1)*MCRD + direction
                  ddl_master1 = (num_cp_master1 - 1)*MCRD + direction
                  ddl_master2 = (num_cp_master2 - 1)*MCRD + direction
                  
c     Modification second membre
                  F(ddl_master2) = F(ddl_master2) - c*F(ddl_slave)
                  F(ddl_master1) = F(ddl_master1) + c*F(ddl_slave)
     &                 + F(ddl_slave)
                  F(ddl_slave)   = 0.0D0
                  
c     Modification matrice de rigidite
                  K(:,ddl_master2) = K(:,ddl_master2) - c*K(:,ddl_slave)
                  K(:,ddl_master1) = K(:,ddl_master1) + c*K(:,ddl_slave)
     &                 + K(:,ddl_slave)
                  K(:,ddl_slave)   = 0.0D0
                  
                  K(ddl_master2,:) = K(ddl_master2,:) - c*K(ddl_slave,:)                  
                  K(ddl_master1,:) = K(ddl_master1,:) + c*K(ddl_slave,:)
     &                 + K(ddl_slave,:)
                  K(ddl_slave,:)   = 0.0D0
               Enddo
            Enddo
            
         End if
      Enddo
      
      
      
      End subroutine coupling
