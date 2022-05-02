!! Copyright 2018-2020 Thibaut Hirschler

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

      
      
      SUBROUTINE gradUELMAT3adj(activeDir,Uelem,UAelem,
     1     NADJ,NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,
     2     TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,indDLoad,
     3     load_target_nbelem,JDLType,ADLMAG,computeWint,computeWext,
     4     gradWint_elem,gradWext_elem)
      
      use parameters
      
      implicit none
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NADJ,NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(3,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
      Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
      Double precision, intent(in) :: ADLMAG
      dimension ADLMAG(nb_load),indDLoad(SUM(load_target_nbelem)),
     &     load_target_nbelem(nb_load),JDLType(nb_load)
      
      Double precision, intent(in) :: Uelem,UAelem
      Integer, intent(in)          :: activeDir
      dimension Uelem(3,NNODE),UAelem(3,NNODE,NADJ),activeDir(3)

      Logical, intent(in) :: computeWint,computeWext
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: gradWint_elem,gradWext_elem
      dimension gradWint_elem(NADJ,3,NNODE),gradWext_elem(NADJ,3,NNODE)
      
      
c     Local variables :
c     ---------------
      
!     For gauss points
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(3,NBINT)
      
!     For nurbs basis functions
      Double precision :: R, dRdxi, ddRddxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI, dAIdxi,AAI,AAE, AIxAJ, Area, Det
      dimension AI(3,3), dAIdxi(3,3), AAI(3,3), AAE(3,3),AIxAJ(3,3)

!     For material matrix
      Integer :: voigt
      dimension voigt(3,2)
      Double precision :: E, nu, lambda, mu, h, matH, coef
      dimension matH(3,3)
      
!     For disp/strain/stress fields (state and adjoint)
      Double precision :: dUdxi,dUAdxi,ddUddxi,ddUAddxi,
     &     strainMem,stressMem,strainMemAdj,stressMemAdj,
     &     strainBnd,stressBnd,strainBndAdj,stressBndAdj,
     &     work,coef1,coef2,UA
      dimension dUdxi(3,2),dUAdxi(3,2,NADJ),ddUddxi(3,3),
     &     ddUAddxi(3,3,NADJ),strainMem(3),stressMem(3),
     &     strainMemAdj(3,nadj),stressMemAdj(3,nadj),strainBnd(3),
     &     stressBnd(3),strainBndAdj(3,nadj),stressBndAdj(3,nadj),
     &     work(nadj),UA(3,nadj)

!     For derivatives (membrane)
      Double precision :: dAAIdP,dAAEdP,dJdP,dEAdP_N,dEdP_NA,dCdP,
     &     dNdP_EA
      dimension dAAIdP(3,2,2),dAAEdP(3,2,2),dJdP(3),
     &     dEAdP_N(3,nadj),dEdP_NA(3,nadj),dCdP(3),dNdP_EA(3,nadj)

!     For derivatives (bending)
      Double precision :: dKAdP_M,dKdP_MA,dMdP_KA,temp,vect,vdA3dP,
     &     vectsave
      dimension dKAdP_M(3,nadj),dKdP_MA(3,nadj),dMdP_KA(3,nadj),vect(3),
     &     vdA3dP(3),vectsave(3)
      
!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      Double precision :: VectNorm,normV,UA_V,UA_dVdP
      dimension VectNorm(3),UA_V(nadj),UA_dVdP(3)
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,j,k,l,iA,ij,kl,kk,ll,cp, KTypeDload, KNumFace
      Double precision :: temp1,temp2
      
!     External functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      INTERFACE
         FUNCTION CROSSPROD(u,v)
         IMPLICIT NONE
         double precision, dimension(3) :: CROSSPROD
         double precision, dimension(3), intent(in) :: u,v
         END FUNCTION CROSSPROD
      END INTERFACE 
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
c     Stiffness matrix and force vector are initialized to zero
      gradWint_elem  = zero
      gradWext_elem  = zero
      
c     Material behaviour
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      lambda = two*lambda*mu/(lambda+two*mu)
      
c     Loads
      kk = 0
      nb_load_bnd = 0
      nb_load_srf = 0
      ind_load_loc(:) = 0
      Do i = 1,nb_load
         If (JDLTYPE(i)/10>0 .AND. 
     &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i)) == JELEM))then
            If (JDLType(i)/10 < 5) then
               nb_load_bnd = nb_load_bnd + 1
               ind_load_loc(nb_load_bnd) = i
            Else
               nb_load_srf = nb_load_srf + 1
               ind_load_loc(nb_load+1-nb_load_srf) = i
            Endif
         Endif
         kk = kk + load_target_nbelem(i)
      Enddo
      

c
c     ..................................................................
c
C     Computation :
c     -----------
      
c     Loop on integration points on main surface
      Do n = 1,NBINT
         
c     PRELIMINARY QUANTITES
c     Computing NURBS basis functions and derivatives
         R      = zero
         dRdxi  = zero
         ddRddxi= zero
         DetJac = zero
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         DetJac = DetJac * GaussPdsCoord(1,n)
         
c     Computing Curvilinear Coordinate objects
c     - Covariant basis vectors and derivatives
         AI(:,:) = zero
         Do i = 1,2
            Do cp = 1,NNODE
               AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDS(:,cp)
            Enddo
         Enddo
         
         call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
         call norm( AIxAJ(:,3),3, Area)
         AI(:,3) = AIxAJ(:,3)/Area
         call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
         call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
         
         dAIdxi(:,:) = zero
         Do i = 1,ntens
            Do cp = 1,NNODE
               dAIdxi(:,i) = dAIdxi(:,i) + ddRddxi(cp,i)*COORDS(:,cp)
            Enddo
         Enddo
c     - Computing metrics
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
         
c     Computing material matrix
         voigt(:,1) = (/ 1,2,1 /)
         voigt(:,2) = (/ 1,2,2 /)
         
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
         matH(:,:) = h*matH(:,:)
         
c     Computing disp and adjoint derivatives
c     - 1st derivatives
         dUdxi(:,:) = zero
         Do i = 1,2
            Do cp = 1,NNODE
               dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         
         dUAdxi(:,:,:) = zero
         Do iA = 1,nadj
         Do i = 1,2
            Do cp = 1,NNODE
               dUAdxi(:,i,iA) 
     &              = dUAdxi(:,i,iA) + dRdxi(cp,i)*UAelem(:,cp,iA)
            Enddo
         Enddo
         Enddo
c     - 2nd derivatives
         ddUddxi(:,:) = zero
         Do i = 1,3
            Do cp = 1,NNODE
               ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         
         ddUAddxi(:,:,:) = zero
         Do iA = 1,nadj
         Do i = 1,3
            Do cp = 1,NNODE
               ddUAddxi(:,i,iA) 
     &              = ddUAddxi(:,i,iA) + ddRddxi(cp,i)*UAelem(:,cp,iA)
            Enddo
         Enddo
         Enddo


         If (computeWint) then
c     Computing strain and stress (state and adjoint)
c     - Membrane
         strainMem(:) = zero
         stressMem(:) = zero
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUdxi(:,i), strainMem(ij))
            Else
               call dot(AI(:,i),dUdxi(:,j), coef1)
               call dot(AI(:,j),dUdxi(:,i), coef2)
               strainMem(ij) = coef1 + coef2
            Endif
         Enddo
         call MulVect(matH,strainMem,stressMem,ntens,ntens)
         
         strainMemAdj(:,:) = zero
         stressMemAdj(:,:) = zero
         Do iA = 1,nadj
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUAdxi(:,i,iA), strainMemAdj(ij,iA))
            Else
               call dot(AI(:,i),dUAdxi(:,j,iA), coef1)
               call dot(AI(:,j),dUAdxi(:,i,iA), coef2)
               strainMemAdj(ij,iA) = coef1 + coef2
            Endif
         Enddo 
         call MulVect(matH,strainMemAdj(:,iA),stressMemAdj(:,iA),ntens,
     &        ntens)       
         Enddo
         
c     - Bending
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         
         strainBnd(:) = zero
         stressBnd(:) = zero
         temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
         Do ij = 1,ntens            
            strainBnd(ij) =
     &           - SUM( ddUddxi(:,ij)*AI(:,3) )
     &           + one/Area*( 
     &             ScaTripleProd(dUdxi(:,1), dAIdxi(:,ij), AI(:,2))
     &             + ScaTripleProd(dUdxi(:,2), AI(:,1), dAIdxi(:,ij)) )
     &           + SUM(AI(:,3)*dAIdxi(:,ij))*temp
         Enddo
         strainBnd(3) = two*strainBnd(3)
         call MulVect(matH,strainBnd,stressBnd,ntens,ntens)
         
         strainBndAdj(:,:) = zero
         stressBndAdj(:,:) = zero
         Do iA = 1,nadj
         temp = SUM(dUAdxi(:,1,iA)*AIxAJ(:,1)+dUAdxi(:,2,iA)*AIxAJ(:,2))
     &           /Area
         Do ij = 1,ntens
            strainBndAdj(ij,iA) =
     &           - SUM( ddUAddxi(:,ij,iA)*AI(:,3) )
     &           + one/Area*( 
     &             ScaTripleProd(dUAdxi(:,1,iA), dAIdxi(:,ij), AI(:,2))
     &           + ScaTripleProd(dUAdxi(:,2,iA), AI(:,1), dAIdxi(:,ij)))
     &           + SUM(AI(:,3)*dAIdxi(:,ij))*temp
         Enddo 
         strainBndAdj(3,iA) = two*strainBndAdj(3,iA)
         call MulVect(matH,strainBndAdj(:,iA),stressBndAdj(:,iA),ntens,
     &        ntens)       
         Enddo
         
         
c     Computing local work
         work(:)  = zero
         Do ij = 1,ntens
            work(:) = work(:) 
     &           + strainMemAdj(ij,:)*stressMem(ij) 
     &           + strainBndAdj(ij,:)*stressBnd(ij)
         Enddo

         Endif  ! test computeWint is True
         
c     --
c     Derivatives
         If (computeWint) then
         Do cp = 1,NNODE
            
            ! 1. derivatives of the jacobian
            dJdP(:) = zero
            dJdP(:) = AIxAJ(:,1)*dRdxi(cp,1) + AIxAJ(:,2)*dRdxi(cp,2)
            Do iA = 1,nadj
               gradWint_elem(iA,:,cp)
     &              = gradWint_elem(iA,:,cp) - work(iA)*dJdP(:)*detJac
            Enddo
            
            
            ! 2. derivatives of the adjoint membrane strain 
            !    (with dble prod. by membrane stress)
            dEAdP_N(:,:) = zero
            Do iA = 1,nadj
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               If (i==j) then
                  dEAdP_N(:,iA) = dEAdP_N(:,iA) 
     &                 + stressMem(ij)*dUAdxi(:,i,iA)*dRdxi(cp,i)
               Else
                  dEAdP_N(:,iA) = dEAdP_N(:,iA) 
     &                 + stressMem(ij)*dUAdxi(:,i,iA)*dRdxi(cp,j)
     &                 + stressMem(ij)*dUAdxi(:,j,iA)*dRdxi(cp,i)
               Endif
            Enddo
            gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &           - dEAdP_N(:,iA)*Area*detJac
            Enddo
            
            
            ! 3. derivatives of the membrane stress 
            !    (with dble prod. by adjoint membrane strain)
            dNdP_EA(:,:) = zero
            !    - derivatives of covariant metrics
            dAAIdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)
     &                 + AI(:,i)*dRdxi(cp,j)
               Enddo
            Enddo
            !    - derivatives of contravariant metrics
            dAAEdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  Do l = 1,2
                     Do k = 1,2
                        dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                       - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
            !    - first subterm (material dot derivative strain)
            dEdP_NA(:,:) = zero
            Do iA = 1,nadj
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               If (i==j) then
                  dEdP_NA(:,iA) = dEdP_NA(:,iA) 
     &                 + stressMemAdj(ij,iA)*dUdxi(:,i)*dRdxi(cp,i)
               Else
                  dEdP_NA(:,iA) = dEdP_NA(:,iA) 
     &                 + stressMemAdj(ij,iA)*dUdxi(:,i)*dRdxi(cp,j)
     &                 + stressMemAdj(ij,iA)*dUdxi(:,j)*dRdxi(cp,i)
               Endif
            Enddo
            Enddo
            dNdP_EA(:,:) = dEdP_NA(:,:)
            
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
                  
                  Do iA = 1,nadj
                     dNdP_EA(:,iA) = dNdP_EA(:,iA) 
     &                    + h*dCdP(:)*strainMem(ij)*strainMemAdj(kl,iA)
                  Enddo
               Enddo
            Enddo

            Do iA = 1,nadj
               gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &              - dNdP_EA(:,iA)*Area*detJac
            Enddo
            
            
            ! 4. derivatives of the adjoint bending strain 
            !    (with dble prod. by bending stress)
            dKAdP_M(:,:) = zero
            Do iA = 1,nadj
            temp = SUM(dUAdxi(:,1,iA)*AIxAJ(:,1))
     &              + SUM(dUAdxi(:,2,iA)*AIxAJ(:,2))
            
            vect(:)   = CROSSPROD(dUAdxi(:,1,iA),AI(:,2)) 
     &                + CROSSPROD(AI(:,1),dUAdxi(:,2,iA))
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            vectsave(:) = vdA3dP(:)
            
            coef = one
            Do ij = 1,ntens
               
               if (ij==3) coef=two

               ! 1st term
               vect(:)   = ddUAddxi(:,ij,iA)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               dKAdP_M(:,iA) =dKAdP_M(:,iA)-vdA3dP(:)*stressBnd(ij)*coef
               
               ! 2nd term
               vect(:) =  
     &              -dJdP(:)/Area/Area*(
     &               ScaTripleProd(dUAdxi(:,1,iA),dAIdxi(:,ij),AI(:,2))
     &              +ScaTripleProd(dUAdxi(:,2,iA),AI(:,1),dAIdxi(:,ij)))
     &              +one/Area*(
     &                CROSSPROD(dUAdxi(:,1,iA),dAIdxi(:,ij))*dRdxi(cp,2)
     &               -CROSSPROD(dUAdxi(:,1,iA),AI(:,2))*ddRddxi(cp,ij))
     &              +one/Area*(
     &                CROSSPROD(dAIdxi(:,ij),dUAdxi(:,2,iA))*dRdxi(cp,1)
     &               -CROSSPROD(AI(:,1),dUAdxi(:,2,iA))*ddRddxi(cp,ij))
               
               dKAdP_M(:,iA) = dKAdP_M(:,iA) +vect(:)*stressBnd(ij)*coef
               
               ! 3rd term
               vect(:)   = dAIdxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               vect(:) = 
     &             -DOTPROD(AI(:,3),dAIdxi(:,ij))*dJdP(:)/Area/Area*temp
     &              +(vdA3dP(:)/Area + AI(:,3)*ddRddxi(cp,ij)/Area)*temp
     &              +DOTPROD(AI(:,3),dAIdxi(:,ij))/Area*(
     &                  CROSSPROD(AI(:,3),dUAdxi(:,1,iA))*dRdxi(cp,2)
     &                + CROSSPROD(dUAdxi(:,2,iA),AI(:,3))*dRdxi(cp,1)
     &                + vectsave(:) )
               
               dKAdP_M(:,iA) = dKAdP_M(:,iA) +vect(:)*stressBnd(ij)*coef
               
               
            Enddo
            gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &           - dKAdP_M(:,iA)*Area*detJac
            Enddo
            
            
            ! 5. derivatives of the bending stress 
            !    (with dble prod. by adjoint bending strain)
            dMdP_KA(:,:) = zero
            !    - derivatives of covariant metrics
            dAAIdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)
     &                 + AI(:,i)*dRdxi(cp,j)
               Enddo
            Enddo
            !    - derivatives of contravariant metrics
            dAAEdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  Do l = 1,2
                     Do k = 1,2
                        dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                       - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
            
            !    - first subterm (material dot derivative strain)
            dKdP_MA(:,:) = zero
            temp = SUM(dUdxi(:,1)*AIxAJ(:,1))+SUM(dUdxi(:,2)*AIxAJ(:,2))
            
            vect(:)   = CROSSPROD(dUdxi(:,1),AI(:,2)) 
     &                + CROSSPROD(AI(:,1),dUdxi(:,2))
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            vectsave(:) = vdA3dP(:)
            
            coef = one
            Do ij = 1,ntens
               
               if (ij==3) coef=two

               ! 1st term
               vect(:)   = ddUddxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA)
     &                 - vdA3dP(:)*stressBndAdj(ij,iA)*coef
               Enddo
               
               ! 2nd term
               vect(:) =  
     &              -dJdP(:)/Area/Area*(
     &               ScaTripleProd(dUdxi(:,1),dAIdxi(:,ij),AI(:,2))
     &              +ScaTripleProd(dUdxi(:,2),AI(:,1),dAIdxi(:,ij)))
     &              +one/Area*(
     &                CROSSPROD(dUdxi(:,1),dAIdxi(:,ij))*dRdxi(cp,2)
     &               -CROSSPROD(dUdxi(:,1),AI(:,2))*ddRddxi(cp,ij))
     &              +one/Area*(
     &                CROSSPROD(dAIdxi(:,ij),dUdxi(:,2))*dRdxi(cp,1)
     &               -CROSSPROD(AI(:,1),dUdxi(:,2))*ddRddxi(cp,ij))
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA) 
     &                 + vect(:)*stressBndAdj(ij,iA)*coef
               Enddo
               
               ! 3rd term
               vect(:)   = dAIdxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               vect(:) = 
     &              (-dJdP(:)/Area/Area*DOTPROD(AI(:,3),dAIdxi(:,ij))
     &              + vdA3dP(:)/Area + AI(:,3)*ddRddxi(cp,ij)/Area)*temp
     &              +DOTPROD(AI(:,3),dAIdxi(:,ij))/Area*(
     &                  CROSSPROD(AI(:,3),dUdxi(:,1))*dRdxi(cp,2)
     &                + CROSSPROD(dUdxi(:,2),AI(:,3))*dRdxi(cp,1)
     &                + vectsave(:) )
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA) 
     &                 + vect(:)*stressBndAdj(ij,iA)*coef
               Enddo
            Enddo
            dMdP_KA(:,:) = dKdP_MA(:,:)
            
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
                  
                  Do iA = 1,nadj
                     dMdP_KA(:,iA) = dMdP_KA(:,iA) 
     &                    + h*h*h/12.0d0*dCdP(:)
     &                      *strainBnd(ij)*strainBndAdj(kl,iA)
                  Enddo
               Enddo
            Enddo
            
            Do iA = 1,nadj
               gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &              - dMdP_KA(:,iA)*Area*detJac
            Enddo
            
         Enddo
         Endif ! test computeWint is True
         
         
         
c     Load vector
         If (computeWext) then
         If (nb_load_srf>0) then
            UA(:,:) = zero
            Do iA = 1,nadj
               Do cp = 1,NNODE
                  UA(:,iA) = UA(:,iA) + R(cp)*UAelem(:,cp,iA)
               Enddo
            Enddo
         Endif
         
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            
            coef = one
            if (KTypeDload == 0) then
               if (KNumFace == 5) then
                  VectNorm(:) = AI(:,3)
               else
                  VectNorm(:) =-AI(:,3)
                  coef = -one
               endif
            elseif (KTypeDload == 1) then
               VectNorm(:) = (/ one, zero, zero/)
            elseif (KTypeDload == 2) then
               VectNorm(:) = (/ zero, one, zero/)
            elseif (KTypeDload == 3) then
               VectNorm(:) = (/ zero, zero, one/)
            elseif (KTypeDload == 4) then
               print*,'Warning: incomplete sensitivity in the load term'
               ! --> the derivative needs to be completed
               call norm(AI(1,:), MCRD, normV)
               VectNorm(:) = AI(:,1)/normV
            elseif (KTypeDload == 5) then
               print*,'Warning: incomplete sensitivity in the load term'
               ! --> the derivative needs to be completed
               call norm(AI(2,:), MCRD, normV)
               VectNorm(:) = AI(:,2)/normV
            elseif (KTypeDload == 6) then
               ! constant vertical load
               VectNorm(:) = (/ zero, zero, AI(3,3) /)
            endif
            
            UA_V(:) = zero
            Do iA = 1,nadj
               UA_V(iA) = DOTPROD(UA(:,iA),VectNorm(:))
            Enddo
            UA_V(:) = UA_V(:)*ADLMAG(numLoad)*detJac
            
            Do cp = 1,NNODE
               dJdP(:) = zero
               dJdP(:) = AIxAJ(:,1)*dRdxi(cp,1) + AIxAJ(:,2)*dRdxi(cp,2)

               Do iA = 1,nadj
                  UA_dVdP(:) = zero
                  if (KTypeDload == 0) then
                     vect(:)   = UA(:,iA)
                     vdA3dP(:) = 
     &                    -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &                    - CROSSPROD(vect(:),AI(:,2))*dRdxi(cp,1)
     &                    - CROSSPROD(AI(:,1),vect(:))*dRdxi(cp,2)
                     vdA3dP(:) = vdA3dP(:)/Area
                     UA_dVdP(:) = vdA3dP(:)*ADLMAG(numLoad)*detJac*coef
                  elseif (KTypeDload == 6) then
                     vect(:) = (/ zero, zero, one/AI(3,3) /)
                     dJdP(:) = CROSSPROD(AI(:,2),vect(:))*dRdxi(cp,1)
     &                       + CROSSPROD(vect(:),AI(:,1))*dRdxi(cp,2)
                  Endif
                  
                  gradWext_elem(iA,:,cp) = gradWext_elem(iA,:,cp) 
     &                 + UA_V(iA)*dJdP(:)
     &                 + UA_dVdP(:)*Area
               Enddo
            Enddo
            
         Enddo
         Endif ! test computeWext is True
         
c     End load traitement


         
      Enddo
c     End of the loop on integration points on main surf
      
      
      
c     Build of element stiffness matrix done.
c
c     ..................................................................
c
c     Loop on load : find boundary loads
      
      
      End SUBROUTINE gradUELMAT3adj
