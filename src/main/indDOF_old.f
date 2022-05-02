!! Copyright 2016-2017 Thibaut Hirschler

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


      subroutine getindDOF(MCRD,nb_bc,bc_target,
     1 bc_target_nbelem,bc_values,nb_dof_tot,nb_dof_bloq,
     2 nb_dof_free,ind_dof_bloq,ind_dof_free)

      ! Input arguments
      integer, intent(in) :: MCRD, nb_bc, bc_target
      integer, intent(in) :: nb_dof_tot,nb_dof_bloq, nb_dof_free
      integer, intent(in) :: bc_target_nbelem
      double precision, intent(in) :: bc_values
      
      dimension bc_target(1000,1000)
      dimension bc_values(1000,2)
      dimension bc_target_nbelem(1000)
      
      
      ! Output variables
      integer, intent(out) :: ind_dof_bloq, ind_dof_free 
      dimension ind_dof_bloq(nb_dof_bloq)
      dimension ind_dof_free(nb_dof_free)


      !Subroutines variables
      integer i,j,compt,test,sctr
      integer nb_dofb,node
      

      compt = 1
      do i=1,nb_bc
         nb_dofb = bc_target_nbelem(i) ! nb of DOF blocked with this BC
         do j=1,nb_dofb         ! loop of DOF blocked
            node = bc_target(i,j) 
            if (bc_values(i,1)==1) then
               if (MCRD==2) then
                  sctr = MCRD*node - 1
               else if (MCRD==3) then
                  sctr = MCRD*node -2
               endif

            else if (bc_values(i,1)==2) then 
               if (MCRD==2)  then
                  sctr = MCRD*node
               else if (MCRD==3) then
                  sctr = MCRD*node - 1
               endif
            else if (bc_values(i,1)==3) then 
               sctr = MCRD*node 
            endif  
            ind_dof_bloq(compt) = sctr
            compt = compt + 1
         enddo
      enddo

      compt = 0
      test = 0
      ! compute ind_dof_free
      do i=1,nb_dof_tot
         test=0
         do j=1,nb_dof_bloq
            if (i==ind_dof_bloq(j)) then
               test=1
            endif
         enddo
         if (test==0) then
            compt=compt+1
            ind_dof_free(compt) = i
         endif
      enddo

      end subroutine getindDOF
