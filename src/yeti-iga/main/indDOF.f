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
      
      subroutine getindDOF(MCRD,nb_bc,bc_target,bc_target_nbelem,
     1     bc_values,bc_target_size,nb_dof_tot,nb_dof_bloq,nb_dof_free,
     2     ind_dof_bloq,ind_dof_free, COUPLG_flag, BNDSTRIP_flag)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: MCRD, nb_bc, bc_target, bc_target_nbelem,
     &     bc_target_size, nb_dof_tot
      dimension bc_target(bc_target_size),bc_target_nbelem(nb_bc)
      
      Double precision, intent(in) :: bc_values
      dimension bc_values(2,nb_bc)

      
c     Output variables :
c     ----------------
      Integer, intent(out) :: nb_dof_bloq, nb_dof_free, ind_dof_bloq,
     &     ind_dof_free
      dimension ind_dof_bloq(nb_dof_tot), ind_dof_free(nb_dof_tot)
      
      Logical, intent(out) :: COUPLG_flag, BNDSTRIP_flag
      
      
c     Local variables :
c     ---------------
      Integer :: num_bc, num_cp, direction, ddl, ind_dof_test,
     &     comp_dof_bloq, comp_dof_free, count
      dimension ind_dof_test(nb_dof_tot)
            
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Determination des ddl bloques ....................................
      
!     Construction vecteur test sur les ddl bloques
      COUPLG_flag   = .FALSE.
      BNDSTRIP_flag = .FALSE.
      ind_dof_test(:) = 0
      count = 0
      Do num_bc = 1,nb_bc
         If (bc_values(1,num_bc) == 11.) then
!     Cas bending strip : non ajout de ddl bloque
            BNDSTRIP_flag = .TRUE.
            
         Elseif (bc_values(1,num_bc) == 10.) then
!     Cas couplage fort entre patches tangents
            COUPLG_flag   = .TRUE.
            Do num_cp = 1,bc_target_nbelem(num_bc)
               ddl = (bc_target(count+num_cp) - 1)*MCRD + 1
               ind_dof_test(ddl:ddl+2) = 1
            Enddo
            
         Else
!     Autre cas de conditions aux limites
            direction = int(bc_values(1,num_bc))
            If (6.<bc_values(1,num_bc).and.bc_values(1,num_bc)<10.) then
               direction = direction - 6
               COUPLG_flag   = .TRUE.
            Endif
            
            Do num_cp = 1,bc_target_nbelem(num_bc)
               ddl = (bc_target(count+num_cp) - 1)*MCRD + direction
               ind_dof_test(ddl) = 1
            Enddo
         Endif
         count = count+bc_target_nbelem(num_bc)
      Enddo
      
!     Construction des vecteurs d'indice des ddl bloques et libres
      ind_dof_bloq(:) = 0
      ind_dof_free(:) = 0
      comp_dof_bloq = 0
      comp_dof_free = 0
      Do ddl = 1,nb_dof_tot
         If (ind_dof_test(ddl) == 1) then
            comp_dof_bloq = comp_dof_bloq + 1
            ind_dof_bloq(comp_dof_bloq) = ddl
         Else
            comp_dof_free = comp_dof_free + 1
            ind_dof_free(comp_dof_free) = ddl
         Endif
      Enddo
      nb_dof_bloq = comp_dof_bloq
      nb_dof_free = comp_dof_free

      End subroutine getindDOF
