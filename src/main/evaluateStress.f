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

      
      subroutine evalStress1(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE,MCRD, stress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,MCRD,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(MCRD,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress
      dimension stress(6)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi
      dimension R(NNODE), dRdxi(NNODE,3)
      
      ! displacement field
      Double precision :: U,dUdxi
      dimension U(3),dUdxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,AAI,AAE, Area, Det
      dimension AI(3,3), AAI(3,3), AAE(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(6,2)
      Double precision :: matH,E,nu,lambda,mu
      dimension matH(2*MCRD,2*MCRD)

      ! strain/stress
      Integer          :: ntens
      Double precision :: strain,stressCurv,eI,e_AI,normV,coef1,coef2,
     &     temp
      dimension strain(6),stressCurv(6),eI(3,3),e_AI(3,3)
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)

      stress(:) = zero
      ntens = 2*MCRD

c     Compute
!     displacement field
      call evalnurbs(Xi,R,dRdxi)     
      dUdxi(:,:) = zero
      Do i=1,MCRD
         Do cp = 1,NNODE
            dUdxi(:MCRD,i) = dUdxi(:MCRD,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      
!     curvilinear quantities
      AI(:,:) = zero
      Do i = 1,MCRD
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      
      If (MCRD==2) then
         call cross(AI(:,1),AI(:,2), AI(:,3))
         call norm( AI(:,3),3, Area)
         AI(:,3) = AI(:,3)/Area
      Endif
      
      AAI(:,:) = zero
      AAE(:,:) = zero
      AAI(3,3) = one
      AAE(3,3) = one
      Do i = 1,MCRD
         Do j = 1,MCRD
            call dot(AI(:,i), AI(:,j), AAI(i,j))
         Enddo
      Enddo
      call MatrixInv(AAE(:MCRD,:MCRD), AAI(:MCRD,:MCRD), det, MCRD)
      
      
c     Computing material matrix
      voigt(:,1) = (/ 1,2,3,1,1,2 /)
      voigt(:,2) = (/ 1,2,3,2,3,3 /)
      
      matH(:,:) = zero
      Do kl = 1,ntens
         k=voigt(kl,1); l=voigt(kl,2)
         Do ij = 1,kl
            i=voigt(ij,1); j=voigt(ij,2)
            
            matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
            
         Enddo
      Enddo
      Do kl = 2,ntens
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      
c     Computing state strain and stress
      strain(:) = zero
      stressCurv(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         If (i==j) then
            call dot(AI(:,i),dUdxi(:,i), strain(ij))
         Else
            call dot(AI(:,i),dUdxi(:,j), coef1)
            call dot(AI(:,j),dUdxi(:,i), coef2)
            strain(ij) = coef1 + coef2
         Endif
      Enddo
      call MulVect(matH,strain(:ntens),stressCurv(:ntens),ntens,ntens)
      
      
!     change of basis
      eI(:,:) = zero
      eI(:,1) = (/ one,zero,zero /)
      eI(:,2) = (/ zero,one,zero /)
      eI(:,3) = (/ zero,zero,one /)
      
      e_AI(:,:) = zero
      do j = 1,MCRD
         do i = 1,MCRD
            call dot(eI(:,i), AI(:,j), e_AI(i,j))
         enddo
      enddo
      
      stress(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         Do kl = 1,MCRD
            stress(ij) = stress(ij) 
     &           + e_AI(i,kl)*e_AI(j,kl)*stressCurv(kl)
         Enddo
         Do kl = 4,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            stress(ij) = stress(ij) 
     &           + e_AI(i,k)*e_AI(j,l)*stressCurv(kl) 
     &           + e_AI(i,l)*e_AI(j,k)*stressCurv(kl)
         Enddo
      Enddo
      
      
C     ------------------------------------------------------------------
      
      end subroutine evalStress1
      




































      
      

c     --
c     Derivative of the stress w.r.t. the degrees of freedom.
c     --
      
      subroutine adjointRHSstress1(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE,MCRD, stress,gradUstress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,MCRD,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(MCRD,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress,gradUstress
      dimension stress(6),gradUstress(6,MCRD,NNODE)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi
      dimension R(NNODE), dRdxi(NNODE,3)
      
      ! displacement field
      Double precision :: dUdxi
      dimension dUdxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,AAI,AAE, Area, Det
      dimension AI(3,3), AAI(3,3), AAE(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(6,2)
      Double precision :: matH,E,nu,lambda,mu
      dimension matH(2*MCRD,2*MCRD)
      
      ! strain/stress
      Integer          :: ntens
      Double precision :: strain,stressCurv,eI,e_AI,normV,coef1,coef2,
     &     temp,vect,matB,gradUstressCurv
      dimension strain(6),stressCurv(6),eI(3,3),e_AI(3,3),
     &     matB(2*MCRD,3),gradUstressCurv(6,MCRD,NNODE)
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)
      
      ntens = 2*MCRD
      stress(:) = zero
      gradUstress(:,:,:) = zero
      
c     Compute
!     displacement field
      call evalnurbs(Xi,R,dRdxi)
      dUdxi(:,:) = zero
      Do i=1,MCRD
         Do cp = 1,NNODE
            dUdxi(:MCRD,i) = dUdxi(:MCRD,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      
!     curvilinear quantities
      AI(:,:) = zero
      Do i = 1,MCRD
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      
      If (MCRD==2) then
         call cross(AI(:,1),AI(:,2), AI(:,3))
         call norm( AI(:,3),3, Area)
         AI(:,3) = AI(:,3)/Area
      Endif
      
      AAI(:,:) = zero
      AAE(:,:) = zero
      AAI(3,3) = one
      AAE(3,3) = one
      Do i = 1,MCRD
         Do j = 1,MCRD
            call dot(AI(:,i), AI(:,j), AAI(i,j))
         Enddo
      Enddo
      call MatrixInv(AAE(:MCRD,:MCRD), AAI(:MCRD,:MCRD), det, MCRD)
      
      
c     Computing material matrix
      voigt(:,1) = (/ 1,2,3,1,1,2 /)
      voigt(:,2) = (/ 1,2,3,2,3,3 /)
      
      matH(:,:) = zero
      Do kl = 1,ntens
         k=voigt(kl,1); l=voigt(kl,2)
         Do ij = 1,kl
            i=voigt(ij,1); j=voigt(ij,2)
            
            matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
            
         Enddo
      Enddo
      Do kl = 2,ntens
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      
      
c     Computing state strain and stress
      strain(:) = zero
      stressCurv(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         If (i==j) then
            call dot(AI(:,i),dUdxi(:,i), strain(ij))
         Else
            call dot(AI(:,i),dUdxi(:,j), coef1)
            call dot(AI(:,j),dUdxi(:,i), coef2)
            strain(ij) = coef1 + coef2
         Endif
      Enddo
      call MulVect(matH,strain(:ntens),stressCurv(:ntens),ntens,ntens)
      
      
!     derivative wrt the dof
      gradUstressCurv(:,:,:) = zero
      Do cp = 1,NNODE
         matB(:,:) = zero
         Do i = 1,MCRD
            matB(i,:) = AI(:,i)*dRdxi(cp,i)
         Enddo
         Do ij = 4,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            matB(ij,:) = AI(:,i)*dRdxi(cp,j) + AI(:,j)*dRdxi(cp,i)
         Enddo
         call MulMat(matH,matB(:,:MCRD),gradUstressCurv(:ntens,:,cp),
     &        ntens,MCRD,ntens)
      Enddo
      
!     change of basis
      eI(:,:) = zero
      eI(:,1) = (/ one,zero,zero /)
      eI(:,2) = (/ zero,one,zero /)
      eI(:,3) = (/ zero,zero,one /)
      
      e_AI(:,:) = zero
      do j = 1,MCRD
         do i = 1,MCRD
            call dot(eI(:,i), AI(:,j), e_AI(i,j))
         enddo
      enddo
      
      stress(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         Do kl = 1,MCRD
            stress(ij) = stress(ij) 
     &           + e_AI(i,kl)*e_AI(j,kl)*stressCurv(kl)
         Enddo
         Do kl = 4,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            stress(ij) = stress(ij) 
     &           + e_AI(i,k)*e_AI(j,l)*stressCurv(kl) 
     &           + e_AI(i,l)*e_AI(j,k)*stressCurv(kl)
         Enddo
      Enddo
      
      gradUstress(:,:,:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         Do kl = 1,MCRD
            gradUstress(ij,:,:) = gradUstress(ij,:,:) 
     &           + e_AI(i,kl)*e_AI(j,kl)*gradUstressCurv(kl,:,:)
         Enddo
         Do kl = 4,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            gradUstress(ij,:,:) = gradUstress(ij,:,:) 
     &           + e_AI(i,k)*e_AI(j,l)*gradUstressCurv(kl,:,:)
     &           + e_AI(i,l)*e_AI(j,k)*gradUstressCurv(kl,:,:)
         Enddo
      Enddo


C     ------------------------------------------------------------------
      
      end subroutine adjointRHSstress1







































      
      

c     --
c     Derivative of the stress w.r.t. the degrees of freedom.
c     --
      
      subroutine partialDervStress1(XI,
     &     COORDSelem,Uelem,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     &     NNODE,MCRD, stress, gradQstress)
      
      use parameters
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI
      dimension XI(3)
      
      Integer, intent(in) :: NNODE,MCRD,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDSelem, Uelem, 
     &     MATERIAL_PROPERTIES, PROPS
      dimension COORDSelem(3,NNODE),Uelem(MCRD,NNODE),
     &     MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
C     Output variables :
c     ----------------
      Double precision, intent(out) :: stress,gradQstress
      dimension stress(6),gradQstress(6,3,NNODE)
      
C     Local variables :
c     ---------------
      ! basis functions
      Double precision :: R,dRdxi
      dimension R(NNODE), dRdxi(NNODE,3)
      
      ! displacement field
      Double precision :: dUdxi
      dimension dUdxi(3,3)
      
      ! curvilinear quantities
      Double precision :: AI,AIxAJ,AAI,AAE, Area, Det
      dimension AI(3,3),AIxAJ(3,3),AAI(3,3),AAE(3,3)
      
      ! material behavior
      Integer :: voigt
      dimension voigt(6,2)
      Double precision :: matH,E,nu,lambda,mu
      dimension matH(2*MCRD,2*MCRD)
      
      ! strain/stress
      Integer          :: ntens
      Double precision :: strain,stressCurv,eI,e_AI,normV,coef1,coef2,
     &     temp,vect
      dimension strain(6),stressCurv(6),eI(3,3),e_AI(3,3)
      
      ! For derivatives
      Double precision :: dAAIdP,dAAEdP,dJdP,dEdP,dCdP,C,dSdP,
     &     gradQstressCurv
      dimension dAAIdP(3,3,3),dAAEdP(3,3,3),dJdP(3),dEdP(3,2*MCRD),
     &     dCdP(3),dSdP(3,2*MCRD),gradQstressCurv(6,3,NNODE)
      
      ! other
      Integer :: i,j,k,l,ij,kl,cp
      
      
C     ------------------------------------------------------------------
      
      ! Initialize
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)
      
      ntens = 2*MCRD
      stress(:)     = zero
      stresscurv(:) = zero
      gradQstress(:,:,:)     = zero
      gradQstressCurv(:,:,:) = zero
      
c     Compute
!     displacement field
      call evalnurbs(Xi,R,dRdxi)
      dUdxi(:,:) = zero
      Do i=1,MCRD
         Do cp = 1,NNODE
            dUdxi(:MCRD,i) = dUdxi(:MCRD,i) + dRdxi(cp,i)*Uelem(:,cp)
         Enddo
      Enddo
      
!     curvilinear quantities
      AI(:,:) = zero
      Do i = 1,MCRD
         Do cp = 1,NNODE
            AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDSelem(:,cp)
         Enddo
      Enddo
      
      call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
      If (MCRD==2) then
         call norm( AIxAJ(:,3),3, Area)
         AI(:,3) = AIxAJ(:,3)/Area
      Else
         call dot(AIxAJ(:,3),AI(:,3),Area)
      Endif
      call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
      call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
      
      AAI(:,:) = zero
      AAE(:,:) = zero
      AAI(3,3) = one
      AAE(3,3) = one
      Do i = 1,MCRD
         Do j = 1,MCRD
            call dot(AI(:,i), AI(:,j), AAI(i,j))
         Enddo
      Enddo
      call MatrixInv(AAE(:MCRD,:MCRD), AAI(:MCRD,:MCRD), det, MCRD)
      
      
c     Computing material matrix
      voigt(:,1) = (/ 1,2,3,1,1,2 /)
      voigt(:,2) = (/ 1,2,3,2,3,3 /)
      
      matH(:,:) = zero
      Do kl = 1,ntens
         k=voigt(kl,1); l=voigt(kl,2)
         Do ij = 1,kl
            i=voigt(ij,1); j=voigt(ij,2)
            
            matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)
     &           + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
            
         Enddo
      Enddo
      Do kl = 2,ntens
         matH(kl:,kl-1) = matH(kl-1,kl:)
      Enddo
      
      
c     Computing state strain and stress
      strain(:) = zero
      stressCurv(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         If (i==j) then
            call dot(AI(:,i),dUdxi(:,i), strain(ij))
         Else
            call dot(AI(:,i),dUdxi(:,j), coef1)
            call dot(AI(:,j),dUdxi(:,i), coef2)
            strain(ij) = coef1 + coef2
         Endif
      Enddo
      call MulVect(matH,strain(:ntens),stressCurv(:ntens),ntens,ntens)
            
      
!     derivative wrt the control points
      Do cp = 1,NNODE
         
         ! 1. Derivatives of the jacobian
         dJdP(:) = zero
         Do i = 1,MCRD
            dJdP(:) = dJdP(:) + AIxAJ(:,i)*dRdxi(cp,i)
         Enddo
         
         
         ! 2. Derivative of the membrane stress
         !    - derivatives of covariant metrics
         dAAIdP(:,:,:) = zero
         Do j = 1,MCRD
            Do i = 1,MCRD
               dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)
     &              + AI(:,i)*dRdxi(cp,j)
            Enddo
         Enddo
         !    - derivatives of contravariant metrics
         dAAEdP(:,:,:) = zero
         Do j = 1,MCRD
            Do i = 1,MCRD
               Do l = 1,MCRD
                  Do k = 1,MCRD
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
         dSdP(:,:) = zero
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
               
               dSdP(:,kl) = dSdP(:,kl) 
     &              + dCdP(:)*strain(ij) + C*dEdP(:,ij)
            Enddo
         Enddo
         
         Do kl = 1,ntens
            gradQstressCurv(kl,:,cp) = dSdP(:,kl)
         Enddo
         
      Enddo
      

!     change of basis
      eI(:,:) = zero
      eI(:,1) = (/ one,zero,zero /)
      eI(:,2) = (/ zero,one,zero /)
      eI(:,3) = (/ zero,zero,one /)
      
      e_AI(:,:) = zero
      do j = 1,MCRD
         do i = 1,MCRD
            call dot(eI(:,i), AI(:,j), e_AI(i,j))
         enddo
      enddo
      
      stress(:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
         Do kl = 1,MCRD
            stress(ij) = stress(ij) 
     &           + e_AI(i,kl)*e_AI(j,kl)*stressCurv(kl)
         Enddo
         Do kl = 4,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            stress(ij) = stress(ij) 
     &           + e_AI(i,k)*e_AI(j,l)*stressCurv(kl) 
     &           + e_AI(i,l)*e_AI(j,k)*stressCurv(kl)
         Enddo
      Enddo
      
      gradQstress(:,:,:) = zero
      Do ij = 1,ntens
         i=voigt(ij,1); j=voigt(ij,2)
      
         ! - link with the derivatives of the stress in covariant basis
         Do kl = 1,MCRD
            gradQstress(ij,:,:) = gradQstress(ij,:,:) 
     &           + e_AI(i,kl)*e_AI(j,kl)*gradQstressCurv(kl,:,:)
         Enddo
         Do kl = 4,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            gradQstress(ij,:,:) = gradQstress(ij,:,:) 
     &           + e_AI(i,k)*e_AI(j,l)*gradQstressCurv(kl,:,:)
     &           + e_AI(i,l)*e_AI(j,k)*gradQstressCurv(kl,:,:)
         Enddo
      Enddo
      
      ! - derivatibes of the change of basis
      Do cp = 1,NNODE
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            Do kl = 1,MCRD
               gradQstress(ij,:,cp) = gradQstress(ij,:,cp) 
     &              + eI(:,i)*dRdxi(cp,kl)*e_AI(j,kl)*stressCurv(kl)
     &              + e_AI(i,kl)*eI(:,j)*dRdxi(cp,kl)*stressCurv(kl)
            Enddo
            Do kl = 4,ntens
               k=voigt(kl,1); l=voigt(kl,2)
               gradQstress(ij,:,cp) = gradQstress(ij,:,cp)
     &              + eI(:,i)*dRdxi(cp,k)*e_AI(j,l)*stressCurv(kl)
     &              + e_AI(i,k)*eI(:,j)*dRdxi(cp,l)*stressCurv(kl)
     &              + eI(:,i)*dRdxi(cp,l)*e_AI(j,k)*stressCurv(kl)
     &              + e_AI(i,l)*eI(:,j)*dRdxi(cp,k)*stressCurv(kl)
            Enddo
         Enddo
      Enddo
      
C     ------------------------------------------------------------------
      
      end subroutine partialDervStress1
