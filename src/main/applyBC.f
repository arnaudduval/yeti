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

      
C     ******************************************************************
      
      subroutine apply_BC(K,F,nb_cp,MCRD,nb_bc,bc_target,
     1     bc_target_nbelem,bc_values,nb_dof_tot,nb_dof_bloq,
     2     nb_dof_free,ind_dof_bloq,ind_dof_free,U)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_cp, MCRD, nb_bc, bc_target,
     &     bc_target_nbelem, nb_dof_tot,nb_dof_bloq, nb_dof_free
      dimension bc_target(1000,1000), bc_target_nbelem(1000)
      
      Double precision, intent(in) :: bc_values
      dimension bc_values(1000,2)
      
c     Output variables :
c     ----------------
      Double precision, intent(inout) :: K,F,U
      Integer, intent(out) :: ind_dof_bloq,ind_dof_free
      dimension K(nb_cp*MCRD,nb_cp*MCRD),F(nb_cp*MCRD),U(nb_cp*MCRD),
     &     ind_dof_bloq(nb_dof_bloq),ind_dof_free(nb_dof_free)
      
c     Local variables :
c     ---------------
      Integer :: i,j,k1,compt,test,sctr,nb_dofb,node

      Integer :: slave, master1, master2, dof1, dof2, numPC, dof
      Double precision :: coef
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Applications des conditions aux limites ..........................
      
      compt = 1
!     Loop on boundary conditions
      do i = 1,nb_bc
         nb_dofb = bc_target_nbelem(i)
!     Loop of DOF blocked
         do j = 1,nb_dofb
            node = bc_target(i,j)
            
c     --
c     Boundary Condition on displacement
            if (bc_values(i,1) <= 3) then
!     Find indice of blocked DOF
               if (bc_values(i,1)==1) then
                  if (MCRD==2) then
                     sctr = MCRD*node - 1
                  else if (MCRD==3) then
                     sctr = MCRD*node - 2
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
               
               U(sctr) = bc_values(i,2) ! imposed displacement in vector U
               
!     Delete DOF in matrix K and take account of BC in F
               do k1 = 1,nb_dof_tot
                  F(k1) = F(k1) - bc_values(i,2)*K(k1,sctr)
                  K(k1,sctr) = 0.0D0
                  K(sctr,k1) = 0.0D0
               enddo
               
               
c     --
c     Patch coupling : G1 continuity
            else if (bc_values(i,1)==10) then
               slave = node
               master1 = bc_target(i+500,j)
               master2 = bc_target(i+600,j)
               
               sctr = (slave-1)*MCRD + 1
               dof1 = (master1-1)*MCRD + 1
               dof2 = (master2-1)*MCRD + 1

               coef = 1.0D0 ! /!\ temporaire -> calcul norm points de controle
               
!     Modification second membre
               F(dof1:dof1+MCRD-1) = F(dof1:dof1+MCRD-1)
     &              + (1.0D0+coef) * F(sctr:sctr+MCRD-1)
               F(dof2:dof2+MCRD-1) = F(dof2:dof2+MCRD-1)
     &              - coef * F(sctr:sctr+MCRD-1)
               F(sctr:sctr+MCRD-1) = 0.0D0
               
!     Modification matrice rigidite
               do numPC = 1,nb_cp
                  dof = (numPC-1)*MCRD + 1
                  if (.NOT. (numPC==master1 .OR. numPC==master2)) then
                     K(dof:dof+MCRD-1, dof1:dof1+MCRD-1)
     &                    = K(dof:dof+MCRD-1, dof1:dof1+MCRD-1)
     &                    + (1.0D0 + coef)
     &                    *   K(dof:dof+MCRD-1, sctr:sctr+MCRD-1)
                     K(dof1:dof1+MCRD-1, dof:dof+MCRD-1)
     &                    = K(dof:dof+MCRD-1, dof1:dof1+MCRD-1)
c     
                     K(dof:dof+MCRD-1, dof2:dof2+MCRD-1)
     &                    = K(dof:dof+MCRD-1, dof2:dof2+MCRD-1)
     &                    - coef * K(dof:dof+MCRD-1, sctr:sctr+MCRD-1)
                     K(dof2:dof2+MCRD-1, dof:dof+MCRD-1)
     &                    = K(dof:dof+MCRD-1, dof2:dof2+MCRD-1)                     
                  end if
               end do
c     
               K(dof2:dof2+MCRD-1, dof2:dof2+MCRD-1)
     &              = K(dof2:dof2+MCRD-1, dof2:dof2+MCRD-1)
     &              + coef**2.D0 * K(sctr:sctr+MCRD-1, sctr:sctr+MCRD-1)
c     
               K(dof2:dof2+MCRD-1, dof1:dof1+MCRD-1)
     &              = K(dof2:dof2+MCRD-1, dof1:dof1+MCRD-1)
     &              - coef*K(sctr:sctr+MCRD-1, dof1:dof1+MCRD-1)
     &              - coef*(coef+1.0D0)
     &              *   K(sctr:sctr+MCRD-1, sctr:sctr+MCRD-1)
c     
               K(dof1:dof1+MCRD-1, dof2:dof2+MCRD-1)
     &              = K(dof1:dof1+MCRD-1, dof2:dof2+MCRD-1)
     &              - coef*K(dof1:dof1+MCRD-1, sctr:sctr+MCRD-1)
     &              - coef*(coef+1.0D0)
     &              *   K(sctr:sctr+MCRD-1, sctr:sctr+MCRD-1)
c     
               K(dof1:dof1+MCRD-1, dof1:dof1+MCRD-1)
     &              = K(dof1:dof1+MCRD-1, dof1:dof1+MCRD-1)
     &              + (coef+1.0D0)
     &              *   K(dof1:dof1+MCRD-1, sctr:sctr+MCRD-1)
     &              + (coef+1.0D0)
     &              *   K(sctr:sctr+MCRD-1, dof1:dof1+MCRD-1)
     &              + (coef+1.0D0)**2.D0
     &              *   K(sctr:sctr+MCRD-1, sctr:sctr+MCRD-1)
c
c
               K(sctr:sctr+MCRD-1,:) = 0.0D0
               K(:,sctr:sctr+MCRD-1) = 0.0D0
               call SetMatrixToIdentity(
     &              K(sctr:sctr+MCRD-1, sctr:sctr+MCRD-1), MCRD)
               
            end if
         enddo
      enddo
      
      
!     Compute ind_dof_free
      compt = 0
      test = 0
      do i = 1,nb_dof_tot
         test = 0
         do j = 1,nb_dof_bloq
            if (i==ind_dof_bloq(j)) then
               test=1
            endif
         enddo
         if (test==0) then
            compt = compt+1
            ind_dof_free(compt) = i
         endif
      enddo


      end subroutine apply_BC
