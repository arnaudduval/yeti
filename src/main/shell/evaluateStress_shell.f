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

      subroutine evalStress3(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE, stress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(3,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress
      dimension stress(6)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi,ddRddxi
      dimension R(NNODE), dRdxi(NNODE,3), ddRddxi(NNODE,6)
      
      ! displacement field
      Double precision :: U,dUdxi,ddUddxi
      dimension U(3),dUdxi(3,2),ddUddxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,dAIdxi,AAI,AAE, AIxAJ, Area, Det
      dimension AI(3,3), dAIdxi(3,3), AAI(3,3), AAE(3,3), AIxAJ(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(3,2)
      Double precision :: matH,E,nu,lambda,mu,h
      dimension matH(3,3)

      ! strain/stress
      Double precision :: strain,stressCurv,eI,Pmtx,e_AI,normV,
     &     coef1,coef2,temp
      dimension strain(6),stressCurv(6),eI(3,3),Pmtx(3,3),
     &     e_AI(2,2)
      
      ! external functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      h = PROPS(2)
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      lambda = two*lambda*mu/(lambda+two*mu)

      stress(:) = zero
      
c     Compute
!     displacement field
      call evalnurbs_w2ndDerv(Xi,R,dRdxi,ddRddxi)
      U(:) = zero
      Do cp = 1,NNODE
         U(:) = U(:) + R(cp)*Uelem(:,cp)
      Enddo
      
      dUdxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      ddUddxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         ddUddxi(:,3) = ddUddxi(:,3) + ddRddxi(cp,4)*Uelem(:,cp)
      Enddo
      
!     curvilinear 
      AI(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
      call norm( AIxAJ(:,3),3, Area)
      AI(:,3) = AIxAJ(:,3)/Area
      call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
      call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
      
      dAIdxi(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            dAIdxi(:,i) = dAIdxi(:,i) +ddRddxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         dAIdxi(:,3) = dAIdxi(:,3) + ddRddxi(cp,4)*COORDSelem(:,cp)
      Enddo
      
!     material
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
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
         Enddo
      Enddo
      Do kl = 2,3
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      matH(:,:) = h*matH(:,:)
      
      
!     membrane part
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
      
!     bending part
      matH(:,:) = h**two/12.0D0 * matH(:,:)
      temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
      Do ij = 1,3            
         strain(ij+3) =
     &        - SUM( ddUddxi(:,ij)*AI(:,3) )
     &        + one/Area*( 
     &        ScaTripleProd(dUdxi(:,1), dAIdxi(:,ij), AI(:,2))
     &        + ScaTripleProd(dUdxi(:,2), AI(:,1), dAIdxi(:,ij)) )
     &        + SUM(AI(:,3)*dAIdxi(:,ij))*temp
      Enddo
      strain(6) = two*strain(6)
      call MulVect(matH,strain(4:6),stressCurv(4:6),3,3)

!     change of basis
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
      
      
C     ------------------------------------------------------------------
      
      end subroutine evalStress3
      




































      
      

c     --
c     Derivative of the stress w.r.t. the degrees of freedom.
c     --
      
      subroutine adjointRHSstress3(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE, stress, gradUstress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(3,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress,gradUstress
      dimension stress(6),gradUstress(6,3,NNODE)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi,ddRddxi
      dimension R(NNODE), dRdxi(NNODE,3), ddRddxi(NNODE,6)
      
      ! displacement field
      Double precision :: U,dUdxi,ddUddxi
      dimension U(3),dUdxi(3,2),ddUddxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,dAIdxi,AAI,AAE, AIxAJ, Area, Det, A3dAIdxi,
     &     dAIdxi_A2, A1_dAIdxi
      dimension AI(3,3), dAIdxi(3,3), AAI(3,3), AAE(3,3), AIxAJ(3,3),
     &     A3dAIdxi(3), dAIdxi_A2(3,3), A1_dAIdxi(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(3,2)
      Double precision :: matH,matPH,E,nu,lambda,mu,h
      dimension matH(3,3),matPH(3,3)

      ! strain/stress
      Double precision :: strain,stressCurv,eI,Pmtx,e_AI,normV,
     &     coef1,coef2,temp,vect,matB
      dimension strain(6),stressCurv(6),eI(3,3),Pmtx(3,3),
     &     e_AI(2,2),vect(3),matB(3,3)
      
      ! external functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      h = PROPS(2)
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      lambda = two*lambda*mu/(lambda+two*mu)
      
      stress(:) = zero
      gradUstress(:,:,:) = zero
      
c     Compute
!     displacement field
      call evalnurbs_w2ndDerv(Xi,R,dRdxi,ddRddxi)
      U(:) = zero
      Do cp = 1,NNODE
         U(:) = U(:) + R(cp)*Uelem(:,cp)
      Enddo
      
      dUdxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      ddUddxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         ddUddxi(:,3) = ddUddxi(:,3) + ddRddxi(cp,4)*Uelem(:,cp)
      Enddo
      
      
!     curvilinear 
      AI(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
      call norm( AIxAJ(:,3),3, Area)
      AI(:,3) = AIxAJ(:,3)/Area
      call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
      call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
      
      dAIdxi(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            dAIdxi(:,i) = dAIdxi(:,i) +ddRddxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         dAIdxi(:,3) = dAIdxi(:,3) + ddRddxi(cp,4)*COORDSelem(:,cp)
      Enddo

!     material
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
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
         Enddo
      Enddo
      Do kl = 2,3
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      matH(:,:) = h*matH(:,:)
      
!     change of basis
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

      matPH(:,:) = zero
      call MulMat(Pmtx,matH,matPH,3,3,3)

!     membrane part
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
      call MulVect(matPH,strain(:3),stress(:3),3,3)

      
!     derivative wrt the dof (membrane part)
      Do cp = 1,NNODE
         matB(:,:) = zero
         matB(1,:) = AI(:,1)*dRdxi(cp,1)
         matB(2,:) = AI(:,2)*dRdxi(cp,2)
         matB(3,:) = AI(:,1)*dRdxi(cp,2) + AI(:,2)*dRdxi(cp,1)
         call MulMat(matPH,matB,gradUstress(1:3,:,cp),3,3,3)
      Enddo
      
      
!     bending part
      Do ij = 1,3
         call dot(AI(:,3),dAIdxi(:,ij),A3dAIdxi(ij))
         call cross(dAIdxi(:,ij),AI(:,2),dAIdxi_A2(:,ij))
         call cross(AI(:,1),dAIdxi(:,ij),A1_dAIdxi(:,ij))
      Enddo
      
      matPH(:,:) = h**two/12.0D0 * matPH(:,:)
      temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
      Do ij = 1,3            
         strain(ij+3) =
     &        - SUM( ddUddxi(:,ij)*AI(:,3) )
     &        + one/Area*( 
     &        DotProd(dUdxi(:,1), dAIdxi_A2(:,ij))
     &        + DotProd(dUdxi(:,2), A1_dAIdxi(:,ij)) )
     &        + A3dAIdxi(ij)*temp
      Enddo
      strain(6) = two*strain(6)
      call MulVect(matPH,strain(4:6),stress(4:6),3,3)
      
      
!     derivative wrt the dof (bending part)      
      Do cp = 1,NNODE
         matB(:,:) = zero
         vect(:) = dRdxi(cp,1)*AIxAJ(:,1)+dRdxi(cp,2)*AIxAJ(:,2)
         Do ij = 1,3
            i = ij
            if (ij==3) i = 4
            matB(ij,:) = -ddRddxi(cp,i)*AI(:,3) +
     &        ( dRdxi(cp,1)*dAIdxi_A2(:,ij)+dRdxi(cp,2)*A1_dAIdxi(:,ij)
     &          + A3dAIdxi(ij)*vect(:)
     &        )/Area
         Enddo
         matB(3,:) = two*matB(3,:)
         call MulMat(matPH,matB,gradUstress(4:6,:,cp),3,3,3)
      Enddo
      
      
C     ------------------------------------------------------------------
      
      end subroutine adjointRHSstress3







































      
      

c     --
c     Derivative of the stress w.r.t. the degrees of freedom.
c     --
      
      subroutine partialDervStress3(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE, stress, gradQstress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(3,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress,gradQstress
      dimension stress(6),gradQstress(6,3,NNODE)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi,ddRddxi
      dimension R(NNODE), dRdxi(NNODE,3), ddRddxi(NNODE,6)
      
      ! displacement field
      Double precision :: U,dUdxi,ddUddxi
      dimension U(3),dUdxi(3,2),ddUddxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,dAIdxi,AAI,AAE, AIxAJ, Area, Det, A3dAIdxi,
     &     dAIdxi_A2, A1_dAIdxi, normAI
      dimension AI(3,3), dAIdxi(3,3), AAI(3,3), AAE(3,3), AIxAJ(3,3),
     &     A3dAIdxi(3), dAIdxi_A2(3,3), A1_dAIdxi(3,3), normAI(3)
      
      ! material behavior
      Integer :: voigt,ntens
      dimension voigt(3,2)
      Double precision :: matH,matPH,E,nu,lambda,mu,h
      dimension matH(3,3),matPH(3,3)

      ! strain/stress
      Double precision :: strain,stressCurv,eI,Pmtx,e_AI,normV,
     &     coef1,coef2,temp,vect,matB,coef
      dimension strain(6),stressCurv(6),eI(3,3),Pmtx(3,3),
     &     e_AI(2,2),vect(3),matB(3,3)
      
!     For derivatives (membrane)
      Double precision :: dAAIdP,dAAEdP,dJdP,dEdP,dCdP,C,dNdP,dNNdP
      dimension dAAIdP(3,2,2),dAAEdP(3,2,2),dJdP(3),dEdP(3,3),dCdP(3),
     &     dNdP(3,3),dNNdP(3,3)

!     For derivatives (bending)
      Double precision :: dKdP,dMdP,dMMdP,vdA3dP,vectsave
      dimension dKdP(3,3),dMdP(3,3),dMMdP(3,3),vdA3dP(3),vectsave(3)

!     For derivatives (basis change)
      Double precision :: deIdP_AI
      dimension deIdP_AI(3,2,2)

      ! external functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      INTERFACE
         FUNCTION CROSSPROD(u,v)
         IMPLICIT NONE
         double precision, dimension(3) :: CROSSPROD
         double precision, dimension(3), intent(in) :: u,v
         END FUNCTION CROSSPROD
      END INTERFACE 
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      h = PROPS(2)
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      lambda = two*lambda*mu/(lambda+two*mu)
      
      ntens = 3
      stress(:)     = zero
      stresscurv(:) = zero
      gradQstress(:,:,:) = zero
      
c     Compute
!     displacement field
      call evalnurbs_w2ndDerv(Xi,R,dRdxi,ddRddxi)
      U(:) = zero
      Do cp = 1,NNODE
         U(:) = U(:) + R(cp)*Uelem(:,cp)
      Enddo
      
      dUdxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      ddUddxi(:,:) = zero
      Do i=1,2
         Do cp = 1,NNODE
            ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         ddUddxi(:,3) = ddUddxi(:,3) + ddRddxi(cp,4)*Uelem(:,cp)
      Enddo
      
      
!     curvilinear 
      AI(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
      call norm( AIxAJ(:,3),3, Area)
      AI(:,3) = AIxAJ(:,3)/Area
      call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
      call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
      
      dAIdxi(:,:) = zero
      Do i = 1,2
         Do cp = 1,NNODE
            dAIdxi(:,i) = dAIdxi(:,i) +ddRddxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      Do cp = 1,NNODE
         dAIdxi(:,3) = dAIdxi(:,3) + ddRddxi(cp,4)*COORDSelem(:,cp)
      Enddo

!     material
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
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
         Enddo
      Enddo
      Do kl = 2,3
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      matH(:,:) = h*matH(:,:)
      
!     change of basis
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

      matPH(:,:) = zero
      call MulMat(Pmtx,matH,matPH,3,3,3)
      
!     membrane part
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
      call MulVect( matH,strain(:3),stresscurv(:3),3,3)
      call MulVect(matPH,strain(:3),stress(:3),3,3)
      
!     bending part
      Do ij = 1,3
         call dot(AI(:,3),dAIdxi(:,ij),A3dAIdxi(ij))
         call cross(dAIdxi(:,ij),AI(:,2),dAIdxi_A2(:,ij))
         call cross(AI(:,1),dAIdxi(:,ij),A1_dAIdxi(:,ij))
      Enddo
      
      matH(:,:)  = h**two/12.0D0 * matH(:,:)
      matPH(:,:) = h**two/12.0D0 * matPH(:,:)
      temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
      Do ij = 1,3            
         strain(ij+3) =
     &        - SUM( ddUddxi(:,ij)*AI(:,3) )
     &        + one/Area*( 
     &        DotProd(dUdxi(:,1), dAIdxi_A2(:,ij))
     &        + DotProd(dUdxi(:,2), A1_dAIdxi(:,ij)) )
     &        + A3dAIdxi(ij)*temp
      Enddo
      strain(6) = two*strain(6)
      call MulVect( matH,strain(4:6),stresscurv(4:6),3,3)
      call MulVect(matPH,strain(4:6),stress(4:6),3,3)
      
      
      
!     derivative wrt the dof
      Do cp = 1,NNODE
         
         ! 1. Derivatives of the jacobian
         dJdP(:) = zero
         dJdP(:) = AIxAJ(:,1)*dRdxi(cp,1) + AIxAJ(:,2)*dRdxi(cp,2)
         
         ! 2. Derivative of the membrane stress
         
         ! - derivatives of covariant metrics
         dAAIdP(:,:,:) = zero
         Do j = 1,2
            Do i = 1,2
               dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)
     &              + AI(:,i)*dRdxi(cp,j)
            Enddo
         Enddo
         ! - derivatives of contravariant metrics
         dAAEdP(:,:,:) = zero
         Do j = 1,2
            Do i = 1,2
               Do l = 1,2
                  Do k = 1,2
                     dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                    - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                  Enddo
               Enddo
            Enddo
         Enddo
         ! - first subterm (material dot derivative strain)
         dEdP(:,:) = zero
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               dEdP(:,ij) = dUdxi(:,i)*dRdxi(cp,i)
            Else
               dEdP(:,ij) =dUdxi(:,i)*dRdxi(cp,j)+dUdxi(:,j)*dRdxi(cp,i)
            Endif
         Enddo
         
         ! - second subterm (derivative material tensor)
         dNdP(:,:) = zero
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               
               dCdP(:) =
     &              lambda*dAAEdP(:,i,j)*AAE(k,l)
     &              + lambda*AAE(i,j)*dAAEdP(:,k,l)
     &              + mu*dAAEdP(:,i,k)*AAE(j,l)
     &              + mu*AAE(i,k)*dAAEdP(:,j,l)
     &              + mu*dAAEdP(:,i,l)*AAE(j,k)
     &              + mu*AAE(i,l)*dAAEdP(:,j,k)
               
               C = lambda*AAE(i,j)*AAE(k,l)
     &              + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
               
               dNdP(:,kl) = dNdP(:,kl) 
     &              + h*dCdP(:)*strain(ij) + h*C*dEdP(:,ij)
            Enddo
         Enddo
         
         ! - change of basis
         dNNdP(:,:) = zero
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            Do i = 1,2
               Do j =1,2
                  ij = i
                  if (i/=j) ij = 3
                  dNNdP(:,kl) = dNNdP(:,kl) 
     &                 + e_AI(k,i)*e_AI(l,j)*dNdP(:,ij)
               Enddo
            Enddo
         Enddo
         
         ! - derivatives of the change of basis
         deIdP_AI(:,:,:) = zero
         
         normAI(:) = zero
         call norm(AI(:,1),3,normAI(1))
         call norm(AI(:,2),3,normAI(2))
         normAI(3) = one
         
         Do j = 1,2
            deIdP_AI(:,1,j) = one/normAI(1)*dRdxi(cp,1)*AI(:,j) 
     &        - one/normAI(1)*dRdxi(cp,1)*eI(:,1)*e_AI(1,j)
         Enddo
         
         Do j = 1,2
            vect(:) = CROSSPROD(eI(:,1),AI(:,j))
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            
            vect(:) = CROSSPROD(AI(:,j),eI(:,3))
            temp    = DOTPROD(vect(:),eI(:,1))
            deIdP_AI(:,2,j) = vdA3dP(:)
     &           + one/normAI(1)*dRdxi(cp,1)*vect(:)
     &           - one/normAI(1)*dRdxi(cp,1)*eI(:,1)*temp
         Enddo
         
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            Do i = 1,2
               Do j =1,2
                  ij = i
                  if (i/=j) ij = 3
                  dNNdP(:,kl) = dNNdP(:,kl) 
     &                 + (deIdP_AI(:,k,i)+eI(:,k)*dRdxi(cp,i))*e_AI(l,j)
     &                   *stressCurv(ij)
     &                 + e_AI(k,i)*(eI(:,l)*dRdxi(cp,j)+deIdP_AI(:,l,j))
     &                   *stressCurv(ij)
               Enddo
            Enddo
         Enddo
         
         Do kl = 1,ntens
            gradQstress(kl,:,cp) = dNNdP(:,kl)
         Enddo
         
         
         
         
         ! 3. Derivative of the bending stress
         temp = SUM(dUdxi(:,1)*AIxAJ(:,1))+SUM(dUdxi(:,2)*AIxAJ(:,2))
         
         vect(:)   = CROSSPROD(dUdxi(:,1),AI(:,2)) 
     &        + CROSSPROD(AI(:,1),dUdxi(:,2))
         vdA3dP(:) = 
     &        -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &        - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &        - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
         vdA3dP(:) = vdA3dP(:)/Area
         vectsave(:) = vdA3dP(:)
         
         coef = one
         Do ij = 1,ntens
            i = ij
            if (ij==3) i = 4

            ! 1st term
            vect(:)   = ddUddxi(:,ij)
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            dKdP(:,ij) = - vdA3dP(:)
            
            ! 2nd term
            vect(:) =  
     &           -dJdP(:)/Area/Area*(
     &           ScaTripleProd(dUdxi(:,1),dAIdxi(:,ij),AI(:,2))
     &           +ScaTripleProd(dUdxi(:,2),AI(:,1),dAIdxi(:,ij)))
     &           +one/Area*(
     &           CROSSPROD(dUdxi(:,1),dAIdxi(:,ij))*dRdxi(cp,2)
     &           -CROSSPROD(dUdxi(:,1),AI(:,2))*ddRddxi(cp,i))
     &           +one/Area*(
     &           CROSSPROD(dAIdxi(:,ij),dUdxi(:,2))*dRdxi(cp,1)
     &           -CROSSPROD(AI(:,1),dUdxi(:,2))*ddRddxi(cp,i))
            dKdP(:,ij) = dKdP(:,ij) + vect(:)
            
            ! 3rd term
            vect(:)   = dAIdxi(:,ij)
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            
            vect(:) = 
     &           (-dJdP(:)/Area/Area*DOTPROD(AI(:,3),dAIdxi(:,ij))
     &           + vdA3dP(:)/Area + AI(:,3)*ddRddxi(cp,i)/Area)*temp
     &           +DOTPROD(AI(:,3),dAIdxi(:,ij))/Area*(
     &           CROSSPROD(AI(:,3),dUdxi(:,1))*dRdxi(cp,2)
     &           + CROSSPROD(dUdxi(:,2),AI(:,3))*dRdxi(cp,1)
     &           + vectsave(:) )
            dKdP(:,ij) = dKdP(:,ij) + vect(:)
         Enddo
         dKdP(:,3) = two*dKdP(:,3)
         
         ! - second subterm (derivative material tensor)
         dMdP(:,:) = zero
         coef = h*h*h/12.d0
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               
               dCdP(:) =
     &              lambda*dAAEdP(:,i,j)*AAE(k,l)
     &              + lambda*AAE(i,j)*dAAEdP(:,k,l)
     &              + mu*dAAEdP(:,i,k)*AAE(j,l)
     &              + mu*AAE(i,k)*dAAEdP(:,j,l)
     &              + mu*dAAEdP(:,i,l)*AAE(j,k)
     &              + mu*AAE(i,l)*dAAEdP(:,j,k)
               
               C = lambda*AAE(i,j)*AAE(k,l)
     &              + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
               
               dMdP(:,kl) = dMdP(:,kl) 
     &              + coef*dCdP(:)*strain(ij+3) + coef*C*dKdP(:,ij)
            Enddo
         Enddo
         
         ! - change of basis
         dMMdP(:,:) = zero
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            Do i = 1,2
               Do j =1,2
                  ij = i
                  if (i/=j) ij = 3
                  dMMdP(:,kl) = dMMdP(:,kl) 
     &                 + e_AI(k,i)*e_AI(l,j)*dMdP(:,ij)
               Enddo
            Enddo
         Enddo
         
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            Do i = 1,2
               Do j =1,2
                  ij = i
                  if (i/=j) ij = 3
                  dMMdP(:,kl) = dMMdP(:,kl) 
     &                 + (deIdP_AI(:,k,i)+eI(:,k)*dRdxi(cp,i))*e_AI(l,j)
     &                   *stressCurv(ij+3)
     &                 + e_AI(k,i)*(eI(:,l)*dRdxi(cp,j)+deIdP_AI(:,l,j))
     &                   *stressCurv(ij+3)
               Enddo
            Enddo
         Enddo
         
         Do kl = 1,ntens
            gradQstress(3+kl,:,cp) = dMMdP(:,kl)
         Enddo
         
      Enddo
      
      
C     ------------------------------------------------------------------
      
      end subroutine partialDervStress3
