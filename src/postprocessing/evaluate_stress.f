!! Copyright 2020 Thibaut Hirschler

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

      
      subroutine evalStress(stress,NumPatch,Xi,
     1     SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,ntens,NNODE,
     3     nb_patch,nb_elem,nb_cp,MCRD)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in) :: ntens,MCRD,NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(:),
     &     JPROPS(nb_patch),NNODE(nb_patch),IEN(:),
     &     nb_elem_patch(nb_patch)
      
!     Analysis solution      
      Double precision, intent(in) :: SOL
      dimension SOL(MCRD,nb_cp)
      
!     Infos
      Integer,          intent(in) :: NumPatch
      Double precision, intent(in) :: Xi
      dimension Xi(3)
      
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress
      dimension stress(ntens)
      
      
C     Local variables :
c     ---------------
      Integer :: sctr
      dimension sctr(NNODE(numPatch))
      Double precision :: R,dRdxi,ddRddxi
      dimension R(NNODE(numPatch)), dRdxi(NNODE(numPatch),3), 
     &     ddRddxi(NNODE(numPatch),6)
      
      ! displacement field
      Double precision :: Uelem,U,dUdxi,ddUddxi
      dimension Uelem(MCRD,NNODE(numPatch)),U(3),dUdxi(3,3),ddUddxi(3,6)
      
      ! curvilinear quantities
      Double precision :: COORDSelem,AI,dAIdxi,AAI,AAE, AIxAJ, Area, Det
      dimension COORDSelem(3,NNODE(numPatch)), AI(3,3), dAIdxi(3,3), 
     &     AAI(3,3), AAE(3,3), AIxAJ(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(3,2)
      Double precision :: matH,E,nu,lambda,mu,h
      dimension matH(3,3)

      ! strain/stress
      Double precision :: strain,stressCurv,eI,Pmtx,e_AI,normV,
     &     coef1,coef2,temp
      dimension strain(ntens),stressCurv(ntens),eI(3,3),Pmtx(3,3),
     &     e_AI(2,2)
      
      ! external functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Extract infos
      
      CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &     weight,nb_elem_patch)
      CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &     NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
      
      
      E = MATERIAL_PROPERTIES(1,NumPatch)
      nu= MATERIAL_PROPERTIES(2,NumPatch)
      

      CALL updateElementNumber(Xi)
      sctr(:) = IEN_patch(:,current_elem)

      COORDSelem(:,:) = zero
      Uelem(:,:)      = zero
      Do cp = 1,NNODE(numPatch)
         COORDSelem(:,cp) = COORDS3D(:,sctr(cp))
         Uelem(:,cp)      = sol(:,sctr(cp))
      Enddo

      ! Initialize
      stress(:) = zero
      
      ! Compute
      If (ELT_TYPE_patch == 'U3') then
         ! KL shell
         h = PROPS(2)
         mu     = E/two/(one+nu)
         lambda = E*nu/(one+nu)/(one-two*nu)
         lambda = two*lambda*mu/(lambda+two*mu)
         
         ! displacement field
         call evalnurbs_w2ndDerv(Xi,R,dRdxi,ddRddxi)
         U(:) = zero
         Do cp = 1,NNODE(numPatch)
            U(:) = U(:) + R(cp)*Uelem(:,cp)
         Enddo
         
         dUdxi(:,:) = zero
         Do i=1,2
            Do cp = 1,NNODE(numPatch)
               dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         ddUddxi(:,:) = zero
         Do i=1,3
            Do cp = 1,NNODE(numPatch)
               ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         
         ! curvilinear 
         AI(:,:) = zero
         Do i = 1,2
            Do cp = 1,NNODE(numPatch)
               AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
            Enddo
         Enddo
         call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
         call norm( AIxAJ(:,3),3, Area)
         AI(:,3) = AIxAJ(:,3)/Area
         call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
         call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
                  
         dAIdxi(:,:) = zero
         Do i = 1,3
            Do cp = 1,NNODE(numPatch)
               dAIdxi(:,i) = dAIdxi(:,i) +ddRddxi(cp,i)*COORDSelem(:,cp)
            Enddo
         Enddo

         ! material
         AAI(:,:) = zero
         AAE(:,:) = zero
         AAI(3,3) = one
         AAE(3,3) = one
         Do j = 1,2
            Do i = 1,2
               call dot(AI(:,i), AI(:,j), AAI(i,j))
            Enddo
         Enddo
         call MatrixInv(AAE(:2,:2), AAI(:2,:2), det, 2)
         
         voigt(:,1) = (/ 1,2,1 /)
         voigt(:,2) = (/ 1,2,2 /)
         matH(:,:) = zero
         Do kl = 1,3
            k=voigt(kl,1); l=voigt(kl,2)
            Do ij = 1,kl
               i=voigt(ij,1); j=voigt(ij,2)
               matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)
     &              + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
            Enddo
         Enddo
         Do kl = 2,3
            matH(kl:,kl-1) = matH(kl-1,kl:)
         Enddo
         matH(:,:) = h*matH(:,:)
         
         
         ! membrane part
         strain(:) = zero
         
         Do ij = 1,3
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUdxi(:,i), strain(ij))
            Else
               call dot(AI(:,i),dUdxi(:,j), coef1)
               call dot(AI(:,j),dUdxi(:,i), coef2)
               strain(ij) = coef1 + coef2
            Endif
         Enddo
         call MulVect(matH,strain(:3),stressCurv(:3),3,3)
         
         ! bending part
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
         Do ij = 1,3            
            strain(ij+3) =
     &           - SUM( ddUddxi(:,ij)*AI(:,3) )
     &           + one/Area*( 
     &             ScaTripleProd(dUdxi(:,1), dAIdxi(:,ij), AI(:,2))
     &             + ScaTripleProd(dUdxi(:,2), AI(:,1), dAIdxi(:,ij)) )
     &           + SUM(AI(:,3)*dAIdxi(:,ij))*temp
         Enddo
         strain(6) = two*strain(6)
         call MulVect(matH,strain(4:6),stressCurv(4:6),3,3)

         ! change of basis
         eI(:,:) = zero
         call norm(AI(:,1), 3, normV)
         eI(:,1) = AI(:,1)/normV
         eI(:,3) = AI(:,3)
         call cross(eI(:,3),eI(:,1),eI(:,2))
         
         e_AI(:,:) = zero
         do j = 1,2
            do i = 1,2
               call dot(eI(:,i), AI(:,j), e_AI(i,j))
            enddo
         enddo
         Pmtx(:,:) = zero
         Pmtx(1,1) = e_AI(1,1)*e_AI(1,1)
         Pmtx(1,2) = e_AI(1,2)*e_AI(1,2)
         Pmtx(1,3) = e_AI(1,1)*e_AI(1,2)*two
         Pmtx(2,1) = e_AI(2,1)*e_AI(2,1)
         Pmtx(2,2) = e_AI(2,2)*e_AI(2,2)
         Pmtx(2,3) = e_AI(2,1)*e_AI(2,2)*two
         Pmtx(3,1) = e_AI(1,1)*e_AI(2,1)
         Pmtx(3,2) = e_AI(1,2)*e_AI(2,2)
         Pmtx(3,3) = e_AI(1,1)*e_AI(2,2)+e_AI(1,2)*e_AI(2,1)

         call MulVect(Pmtx, stressCurv(1:3), stress(1:3), 3, 3)
         call MulVect(Pmtx, stressCurv(4:6), stress(4:6), 3, 3)
         
         
      Endif
      
      CALL finalizeNurbsPatch()
      
C     ------------------------------------------------------------------
      
      end subroutine evalStress
      
