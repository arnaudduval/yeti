!! Copyright 2018 Thibaut Hirschler

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
      
C     Calcul des deformations, contraintes et deplacements au niveau des
c     frontiere delimitant un element nurbs. Cela permettra de creer le
c     fichier VTU pour la visualisation.
      subroutine compute_svars_Q1(COORDS,sol,svars,nsvint,Output_FLAG,
     1     nb_vertice,nb_REF,MCRD,NNODE,MATERIAL_PROPERTIES,TENSOR)
      
      use parameters

      implicit none

C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: COORDS,sol,MATERIAL_PROPERTIES
      dimension COORDS(MCRD,NNODE)
      dimension sol(MCRD,NNODE)
      dimension MATERIAL_PROPERTIES(2)
      Character(len=*), intent(in) :: TENSOR
      
      Integer,intent(in) :: nsvint,MCRD,NNODE,nb_vertice,nb_REF
      dimension nb_REF(3)
      
      Logical, intent(in) :: Output_FLAG
      dimension Output_FLAG(3)
      
C     Output variables :
c     ----------------
      Double precision, intent(inout) :: svars
      dimension svars(nsvint*nb_vertice)
      
C     Local variables :
c     ---------------
      Double precision :: vertice, R, dRdx, DetJac, stran, stress, 
     1     svarsip,coords_ip,u_ip,ddsdde,dNidx,dNidy,dNidz
      dimension vertice(MCRD,nb_vertice),R(NNODE),dRdx(MCRD,NNODE),
     &     stran(2*MCRD),stress(2*MCRD),svarsip(nsvint),coords_ip(3),
     &     u_ip(3),ddsdde(2*MCRD,2*MCRD)
      
      Integer :: NDOFEL,n,k1,ntens,nb_xi,nb_eta,nb_zeta,i_xi,
     &     i_eta,i_zeta,offset,i
      
      
C     ------------------------------------------------------------------
      
C     Initialization :
c     --------------
      
      NDOFEL = NNODE*MCRD  
      
c     Get material behaviour tensor
      call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde) 
      
c     Defining element bounds : coords in parent space
      if (MCRD==2) then
         vertice(:,1) = (/-one, -one/)
         vertice(:,2) = (/ one, -one/)
         vertice(:,3) = (/ one,  one/)
         vertice(:,4) = (/-one,  one/)
      else if (MCRD==3) then
         vertice(:,1) = (/-one, -one, -one/)
         vertice(:,2) = (/-one,  one, -one/)
         vertice(:,3) = (/-one, -one,  one/)
         vertice(:,4) = (/-one,  one,  one/)
         vertice(:,5) = (/ one, -one, -one/)
         vertice(:,6) = (/ one,  one, -one/)
         vertice(:,7) = (/ one, -one,  one/)
         vertice(:,8) = (/ one,  one,  one/)
      endif
      
      
      
C     ------------------------------------------------------------------
      
C     Compute disp., stress, strain :
c     -----------------------------
      
      R = zero
      dRdx = zero
      svars = zero

c     do n = 1,nb_vertice
      nb_xi  = 2**max(nb_REF(1)-1,0)+1
      nb_eta = 2**max(nb_REF(2)-1,0)+1
      nb_zeta= 2**max(nb_REF(3)-1,0)+1
      if (MCRD==2) nb_zeta= 1
      
      Do i_zeta= 1,nb_zeta
      Do i_eta = 1,nb_eta
      Do i_xi  = 1,nb_xi

         n = (i_zeta-1)*nb_eta*nb_xi + (i_eta-1)*nb_xi + i_xi
         vertice(1,n) = two/dble(nb_xi -1)*dble(i_xi -1) - one
         vertice(2,n) = two/dble(nb_eta-1)*dble(i_eta-1) - one
         If (MCRD==3) then
            vertice(3,n) = two/dble(nb_zeta-1)*dble(i_zeta-1) - one
         Endif

         call shap (dRdx,R,DetJac,COORDS,vertice(1:,n),MCRD)
         
c     Get integration points coordinates and displacements
         coords_ip = zero
         u_ip = zero
         do k1 = 1,NNODE
            coords_ip(:MCRD) = coords_ip(:MCRD)+R(k1)*coords(:,k1)
            if (Output_FLAG(1)) then
               u_ip(:MCRD) = u_ip(:MCRD) + R(k1)*sol(:,k1)
            endif
         enddo
         
c     Get strain
         ntens = 2*MCRD
         stran = zero
         If (Output_FLAG(2) .OR. Output_FLAG(3)) then
            Do k1 = 1,NNODE
               dNidx = dRdx(1,k1)
               dNidy = dRdx(2,k1)
               
               stran(1) = stran(1) + dNidx*sol(1,k1)
               stran(2) = stran(2) + dNidy*sol(2,k1)
               stran(4) = stran(4) + dNidy*sol(1,k1)+dNidx*sol(2,k1)
               If (MCRD==3) then
                  dNidz = dRdx(3,k1)
                  
                  stran(3) = stran(3) + dNidz*sol(3,k1)
                  stran(5) = stran(5) + dNidz*sol(1,k1)+dNidx*sol(3,k1)
                  stran(6) = stran(6) + dNidy*sol(3,k1)+dNidz*sol(2,k1)
               Endif
            Enddo

c     Get stress
            If (Output_FLAG(2)) then
               call MulVect(ddsdde, stran, stress, ntens, ntens)
            Endif
         Endif
         
c     Sum up all variables into svarsip
         svarsip = zero
         offset = 1
         
         svarsip(offset:offset+2) = coords_ip(:)
         offset = offset + 3
         
         If (Output_FLAG(1)) then
            svarsip(offset:offset+2) = u_ip(:)
            offset = offset + 3
         Endif
         
         If (Output_FLAG(2)) then
            svarsip(offset:offset+ntens-1) = stress(:)
            offset = offset + ntens
         Endif
         
         If (Output_FLAG(3)) then
            svarsip(offset:offset+ntens-1) = stran(:)
         Endif
         
         
c     Update global variable : all variables at each intergration point
         do i = 1,nsvint
            svars(nsvint*(n-1)+i) = svarsip(i)
         enddo
         
      Enddo
      Enddo
      Enddo

C     ------------------------------------------------------------------
      
      End subroutine compute_svars_Q1
