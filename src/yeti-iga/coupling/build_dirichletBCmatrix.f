!! Copyright 2016-2019 Thibaut Hirschler

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
c
c     Creation d'une matrice pour l'imposition des conditions aux 
c     limites de Dirichlet par multiplicateurs de lagrange.
c     Renvoie egalement le second membre (composante nun nulle si
c     deplacement impose par exemple)
      
      subroutine cplg_dirichlet(Cdata,Ccol,Crow,uBC,
     1     nb_data,MCRD,nb_bc,bc_target,bc_target_nbelem,bc_values,
     2     bc_target_size,nb_dof_tot)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_data,MCRD, nb_bc, bc_target, 
     &     bc_target_nbelem,bc_target_size, nb_dof_tot
      dimension bc_target(bc_target_size),bc_target_nbelem(nb_bc)
      
      Double precision, intent(in) :: bc_values
      dimension bc_values(2,nb_bc)
      
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: uBC,Cdata
      Integer,          intent(out) :: Crow,Ccol
      dimension uBC(nb_data),Cdata(nb_data),Crow(nb_data),Ccol(nb_data)
      
c     Local variables :
c     ---------------
      Integer :: num_bc,num_cp,direction,ddl,count,ind_dof_test
      dimension ind_dof_test(nb_dof_tot)
      Double precision :: value_pre_dof
      dimension value_pre_dof(nb_dof_tot)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Determination des ddl bloques ....................................
      
!     Construction vecteur test sur les ddl bloques
      ind_dof_test(:) = 0
      value_pre_dof(:)= 0.d0
      count = 0
      Do num_bc = 1,nb_bc
         direction = int(bc_values(1,num_bc))
         If (direction<4) then
         Do num_cp = 1,bc_target_nbelem(num_bc)
            ddl = (bc_target(count+num_cp) - 1)*MCRD + direction
            ind_dof_test(ddl) = 1
            value_pre_dof(ddl)= bc_values(2,num_bc)
         Enddo
         Endif
         count = count+bc_target_nbelem(num_bc)
      Enddo
      
!     Construction des vecteurs d'indice des ddl bloques et libres
      count = 1
      Do ddl = 1,nb_dof_tot
         If (ind_dof_test(ddl) == 1) then
            uBC(  count) = value_pre_dof(ddl)
            Crow( count) = ddl   -1 ! python indices start at 0
            Ccol( count) = count -1 ! python indices start at 0
            Cdata(count) = 1.d0
            count = count+1
         Endif
      Enddo
      
      End subroutine cplg_dirichlet
