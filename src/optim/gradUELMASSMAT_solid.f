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

      
c     --
c     Construction de la sentitivite de la matrice de rigite et de la 
c     matrice de masse  pour la formulation solide classique.
c     --
      
      SUBROUTINE gradUELMASSMAT1(activeDir,VECTelem,VALS,nb_frq,
     1     NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,TENSOR,
     2     MATERIAL_PROPERTIES,Rho,PROPS,JPROPS,gradV_elem,massElem)
      
      use parameters
      use nurbspatch
      
      implicit none
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES,Rho,
     &     PROPS
      dimension COORDS(3,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
      Double precision, intent(in) :: VECTelem,VALS
      Integer,          intent(in) :: nb_frq,activeDir
      dimension VECTelem(nb_frq,3,NNODE),VALS(nb_frq),activeDir(3)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: gradV_elem,massElem
      dimension gradV_elem(nb_frq,3,NNODE),massElem(nb_frq)
      
      
c     Local variables :
c     ---------------
      
!     For gauss points
      Integer :: NbPtInt, n
      Double precision :: GaussPdsCoord,PtGauss
      dimension GaussPdsCoord(MCRD+1,NBINT),PtGauss(MCRD+1)
      
!     For nurbs basis functions
      Double precision :: XI, R, dRdxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,3), XI(3)
      
!     For curvilinear coordinate objects
      Double precision :: AI,AIxAJ,AAI,AAE,det,Area
      dimension AI(3,3),AIxAJ(3,3),AAI(3,3),AAE(3,3)
      
!     For material matrix
      Integer :: voigt
      dimension voigt(6,2)
      Double precision :: E,nu,lambda,mu, matH,coef
      dimension matH(2*MCRD,2*MCRD)
      
!     For disp/strain/stress fields (state and adjoint)
      Double precision :: U,dUdxi,strain,stress,work,coef1,coef2
      dimension U(3,nb_frq),dUdxi(3,3,nb_frq),strain(2*MCRD,nb_frq),
     &     stress(2*MCRD,nb_frq),work(nb_frq)
      
!     For derivatives
      Double precision :: dAAIdP,dAAEdP,dJdP,dEdP_S,dCdP,dSdP_E
      dimension dAAIdP(3,3,3),dAAEdP(3,3,3),dJdP(3),dEdP_S(3,nb_frq),
     &     dCdP(3),dSdP_E(3,nb_frq)
      
!     For loops
      Integer ntens
      Integer k1,k2,i,j,k,l,iV,kk,ll,ij,kl,cp
      Double precision :: temp,temp1,temp2
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens   = 2*MCRD          ! size of stiffness tensor
      NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
      if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and weights
      call Gauss(NbPtInt,MCRD,GaussPdsCoord,0)
      
c     Stiffness matrix and force vector are initialized to zero
      gradV_elem(:,:,:)  = zero
      massElem(:) = zero
      
c     Material behaviour
      E  = MATERIAL_PROPERTIES(1)
      nu = MATERIAL_PROPERTIES(2)
      lambda = E*nu/(one+nu)/(one-two*nu)
      mu     = E/two/(one+nu)
      if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)
      
c
c     ..................................................................
c
C     Computation :
c     -----------
      
c     Loop on integration points on main surface
      do n = 1,NBINT
         
c     PRELIMINARY QUANTITES
c     Computing NURBS basis functions and derivatives
         R(:)       = zero
         dRdxi(:,:) = zero
         
         XI(:)      = zero
         PtGauss(:) = GaussPdsCoord(2:,n)
         DetJac     = GaussPdsCoord( 1,n)
         Do i = 1,MCRD
            XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &            +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
            DetJac  = DetJac * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
         End do
         call evalnurbs(XI,R,dRdxi)
         
c     Computing Covariant basis vectors
         AI(:,:) = zero
         Do i = 1,MCRD
            Do cp = 1,NNODE
               AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDS(:,cp)
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
         
         

c     Computing metrics
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
     &              + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
               
            Enddo
         Enddo
         Do kl = 2,ntens
            matH(kl:,kl-1) = matH(kl-1,kl:)
         Enddo
         
c     Computing disp and adjoint derivative
         U(:,:) = zero
         Do iV = 1,nb_frq
            Do cp = 1,NNODE
               U(:,iV) = U(:,iV) + R(cp)*VECTelem(iV,:,cp)
            Enddo
         Enddo
         
         dUdxi(:,:,:) = zero
         Do iV = 1,nb_frq
         Do i = 1,MCRD
            Do cp = 1,NNODE
               dUdxi(:,i,iV) 
     &              = dUdxi(:,i,iV) + dRdxi(cp,i)*VECTelem(iV,:,cp)
            Enddo
         Enddo
         Enddo
         
c     Computing strain and stress
         strain(:,:) = zero
         stress(:,:) = zero
         Do iV = 1,nb_frq
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUdxi(:,i,iV), strain(ij,iV))
            Else
               call dot(AI(:,i),dUdxi(:,j,iV), coef1)
               call dot(AI(:,j),dUdxi(:,i,iV), coef2)
               strain(ij,iV) = coef1 + coef2
            Endif
         Enddo 
         call MulVect(matH,strain(:,iV),stress(:,iV),ntens,ntens)       
         Enddo
         
c     Computing local work
         work(:)  = zero
         Do ij = 1,ntens
            work(:) = work(:) + strain(ij,:)*stress(ij,:)
         Enddo
         
c     Compute mass
         Do iV = 1,nb_frq
            massElem(iV) = massElem(iV) 
     &           + Rho*SUM(U(:,iV)*U(:,iV))*Area*detJac
         Enddo
         
c     --
c     Derivatives
         Do cp = 1,NNODE
         
            ! 1. derivatives of the jacobian
            dJdP(:) = zero
            Do i = 1,MCRD
               dJdP(:) = dJdP(:) + AIxAJ(:,i)*dRdxi(cp,i)
            Enddo
            Do iV = 1,nb_frq
               gradV_elem(iV,:,cp) = gradV_elem(iV,:,cp)
     &              + work(iV)*dJdP(:)*detJac
     &              - VALS(iV)*Rho*SUM(U(:,iV)*U(:,iV))*dJdP(:)*detJac
            Enddo
            
            
            ! 2. derivatives of the strain 
            !    (with dble prod. by stress)
            dEdP_S(:,:) = zero
            Do iV = 1,nb_frq
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               If (i==j) then
                  dEdP_S(:,iV) = dEdP_S(:,iV) 
     &                 + stress(ij,iV)*dUdxi(:,i,iV)*dRdxi(cp,i)
               Else
                  dEdP_S(:,iV) = dEdP_S(:,iV) 
     &                 + stress(ij,iV)*dUdxi(:,i,iV)*dRdxi(cp,j)
     &                 + stress(ij,iV)*dUdxi(:,j,iV)*dRdxi(cp,i)
               Endif
            Enddo
            gradV_elem(iV,:,cp) = gradV_elem(iV,:,cp) 
     &           + dEdP_S(:,iV)*Area*detJac
            Enddo
            
            
            ! 3. derivatives of the stress 
            !    (with dble prod. by strain)
            dSdP_E(:,:) = zero
            !    - derivatives of covariant metrics
            dAAIdP(:,:,:) = zero
            Do j = 1,MCRD
               Do i = 1,MCRD
                  dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)
     &                 + AI(:,i)*dRdxi(cp,j)
               Enddo
            Enddo
            !    - derivatives of contravariant metrics
            dAAEdP(:,:,:) = zero
            Do j = 1,MCRD
               Do i = 1,MCRD
                  Do l = 1,MCRD
                     Do k = 1,MCRD
                        dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                       - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
            !    - first subterm (material dot derivative strain)
            dSdP_E(:,:) = dEdP_S(:,:)
            
            !    - second subterm (derivative material tensor)
            Do kl = 1,ntens
               k=voigt(kl,1); l=voigt(kl,2)
               
               Do ij = 1,ntens
                  i=voigt(ij,1); j=voigt(ij,2)
                  
                  dCdP(:) =
     &                   lambda*dAAEdP(:,i,j)*AAE(k,l)
     &                 + lambda*AAE(i,j)*dAAEdP(:,k,l)
     &                 + mu*dAAEdP(:,i,k)*AAE(j,l)
     &                 + mu*AAE(i,k)*dAAEdP(:,j,l)
     &                 + mu*dAAEdP(:,i,l)*AAE(j,k)
     &                 + mu*AAE(i,l)*dAAEdP(:,j,k)
                  
                  Do iV = 1,nb_frq
                     dSdP_E(:,iV) = dSdP_E(:,iV) 
     &                    + dCdP(:)*strain(ij,iV)*strain(kl,iV)
                  Enddo
               Enddo
            Enddo

            Do iV = 1,nb_frq
               gradV_elem(iV,:,cp) = gradV_elem(iV,:,cp) 
     &              + dSdP_E(:,iV)*Area*detJac
            Enddo
         Enddo
         
      Enddo
c     End of the loop on integration points on main surf
c     
c     ..................................................................
c     
      End SUBROUTINE gradUELMASSMAT1



      
