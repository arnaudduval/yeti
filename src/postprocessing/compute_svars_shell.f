!! Copyright 2018-2019 Thibaut Hirschler

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
      subroutine compute_svars_Q1_shell(COORDS,sol,svars,nsvint,
     1     Output_FLAG,nb_vertice,nb_REF,MCRD,NNODE,MATERIAL_PROPERTIES,
     2     TENSOR,PROPS,JPROPS)
      
      use parameters

      Implicit None
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: COORDS,sol,MATERIAL_PROPERTIES,
     &     PROPS
      dimension COORDS(MCRD,NNODE)
      dimension sol(MCRD,NNODE)
      dimension MATERIAL_PROPERTIES(2)
      dimension PROPS(JPROPS)
      Character(len=*) , intent(in) :: TENSOR
      
      Integer, intent(in) :: nsvint,MCRD,NNODE,nb_vertice,JPROPS,nb_REF
      dimension nb_REF(3)
      
      Logical, intent(in) :: Output_FLAG
      dimension Output_FLAG(3)
      
C     Output variables :
c     ----------------
      Double precision, intent(inout) :: svars
      dimension svars(nsvint*nb_vertice)
      
C     Local variables :
c     ---------------
!     For gauss points
      Double precision :: vertice
      dimension vertice(2,nb_vertice)
      
!     For nurbs basis functions
      Double precision :: R, dRdxi, ddRddxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI, dAI1dxi, dAI2dxi, AAE, eI, normV, AE
      dimension AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3), AAE(2,2), eI(3,3),
     &     AE(3,3)

!     For material matrix
      Double precision :: E, nu, h, matH, coef
      dimension matH(3,3)
      
!     Other
      Double precision :: stran, svarsip, coords_ip,u_ip,stranC,stressC,
     &     e_AE, Pmtx
      dimension stran(2*MCRD), svarsip(nsvint),coords_ip(3), u_ip(3),
     &     stranC(2*MCRD), stressC(2*MCRD), e_AE(2,2), Pmtx(3,3)
      
      Integer :: n,k1,ntens,i,j,nb_xi,nb_eta,i_xi,i_eta,offset
      
!     For test rotation field
      Double precision :: vect,phi1,phi2,Psi,w_ip,temp,dAI3dxi,dAb3dxi,
     &     tensB,d1,d2,d3_1,d3_2,rot, vectT,vectN,a3,alpha
      dimension vect(3),Psi(3),w_ip(3),dAI3dxi(2,3),dAb3dxi(2,3),
     &     tensB(2,2),rot(3), vectT(3),vectN(3),a3(3)

C     ------------------------------------------------------------------
      
C     Initialization :
c     --------------
      
      ntens = 6
      
c     Get material behaviour tensor
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      coef = E/(one-nu**two)
      matH = zero
      
      
c     Defining element bounds : coords in parent space
      vertice(:,1) = (/-one, -one/)
      vertice(:,2) = (/ one, -one/)
      vertice(:,3) = (/ one,  one/)
      vertice(:,4) = (/-one,  one/)
      
      
      
C     ------------------------------------------------------------------
      
C     Compute disp., stress, strain :
c     -----------------------------
      
      R = zero
      dRdxi  = zero
      ddRddxi= zero
      DetJac = zero
      svars  = zero
!      do n = 1,nb_vertice
      
      nb_xi  = 2**max(nb_REF(1)-1,0)+1
      nb_eta = 2**max(nb_REF(2)-1,0)+1

      Do i_eta = 1,nb_eta
      Do i_xi  = 1,nb_xi
         
         n = (i_eta-1)*nb_xi + i_xi
         vertice(1,n) = two/dble(nb_xi -1)*dble(i_xi -1) - one
         vertice(2,n) = two/dble(nb_eta-1)*dble(i_eta-1) - one
         
c     Computing NURBS basis functions and derivatives
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,vertice(1:,n))
         
         
c     Get intergration points coordinates and displacements
         coords_ip = zero
         u_ip = zero
         do k1 = 1,NNODE
            coords_ip(:MCRD) = coords_ip(:MCRD)+R(k1)*coords(:,k1)
            if (Output_FLAG(1)) then
               u_ip(:MCRD) = u_ip(:MCRD) + R(k1)*sol(:,k1)
            endif
         enddo
         
         

         
         If (.False.) then
!     test calcul du champ de rotation
            AI = zero
            Do k1 = 1,NNODE
               AI(:,1) = AI(:,1) + dRdxi(k1,1)*COORDS(:,k1)
               AI(:,2) = AI(:,2) + dRdxi(k1,2)*COORDS(:,k1)
            Enddo
            call cross(AI(:,1), AI(:,2), vect)
            normV = sqrt( vect(1)**two + vect(2)**two + vect(3)**two )
            AI(:,3) = vect(:)/normV

            phi1 = zero
            phi2 = zero
            Do k1 = 1,NNODE
               call dot(sol(:,k1), AI(:,3), temp)
               phi1 = phi1 + dRdxi(k1,2)*temp
               phi2 = phi2 - dRdxi(k1,1)*temp
            Enddo
            phi1 = phi1/normV
            phi2 = phi2/normV
            
            Psi(:) = phi1*AI(:,1) + phi2*AI(:,2)
            call cross(Psi(:),AI(:,3), w_ip(:))
            
            
            a3(:) = AI(:,3) + w_ip(:)
            vectT(:) = (/ zero,one,zero /)
            call norm(vectT(:), 3, temp)
            vectT(:) = vectT(:)/temp

            call cross(AI(:,3), vectT(:), vectN(:))
            call dot(a3(:),vectN(:), temp)
            temp  = MIN(temp, one)
            temp  = MAX(temp,-one)
            alpha = ASIN(-temp)

!     ** test calcul rotation totale pour couplage
            call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,R,dRdxi,ddRddxi,
     &           MCRD,NNODE,COORDS)

c     Get A3 derivatives
            call cross(dAI1dxi(1,:),AI(2,:), dAb3dxi(1,:))
            call cross(AI(1,:),dAI2dxi(1,:), vect)
            dAb3dxi(1,:) = dAb3dxi(1,:) + vect(:)
            call cross(dAI1dxi(2,:),AI(2,:), dAb3dxi(2,:))
            call cross(AI(1,:),dAI2dxi(2,:), vect)
            dAb3dxi(2,:) = dAb3dxi(2,:) + vect(:)
            
            call dot(AI(3,:),dAb3dxi(1,:), temp)
            dAI3dxi(1,:) = (dAb3dxi(1,:) - AI(3,:)*temp)/normV
            call dot(AI(3,:),dAb3dxi(2,:), temp)
            dAI3dxi(2,:) = (dAb3dxi(2,:) - AI(3,:)*temp)/normV
            
c     Get covariant basis
            AE(:,:) = zero
            AE(1,:) = AAE(1,1)*AI(1,:) + AAE(1,2)*AI(2,:)
            AE(2,:) = AAE(2,1)*AI(1,:) + AAE(2,2)*AI(2,:)
            AE(3,:) = AI(3,:)
            
c     Build tensB
            call dot(dAI1dxi(1,:),AI(3,:), tensB(1,1))
            call dot(dAI1dxi(2,:),AI(3,:), tensB(1,2))
            call dot(dAI2dxi(1,:),AI(3,:), tensB(2,1))
            call dot(dAI2dxi(2,:),AI(3,:), tensB(2,2))
            
c     Get disp
            call dot(u_ip(:),AE(1,:), d1)
            call dot(u_ip(:),AE(2,:), d2)
            
            call dot(u_ip(:),dAI3dxi(1,:), d3_1)
            vect = zero
            Do k1 = 1,NNODE
               vect(:) = vect(:) + dRdxi(k1,1)*sol(:,k1)
            Enddo
            call dot(vect,AI(3,:),temp)
            d3_1 = d3_1 + temp
            
            call dot(u_ip(:),dAI3dxi(2,:), d3_2)
            vect = zero
            Do k1 = 1,NNODE
               vect(:) = vect(:) + dRdxi(k1,2)*sol(:,k1)
            Enddo
            call dot(vect,AI(3,:),temp)
            d3_2 = d3_2 + temp
            
            !d3_1 = zero
            !d3_2 = zero

c     Get tot rotation
            rot(:) = (d3_2 + tensB(2,1)*d1 + tensB(2,2)*d2)*AI(1,:)
     &             - (d3_1 + tensB(1,1)*d1 + tensB(1,2)*d2)*AI(2,:)
            rot(:) = rot(:)/normV
         Endif
         
         


c     Computing Curvilinear Coordinate objects
         If (Output_FLAG(2) .OR. Output_FLAG(3)) then
            call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,R,dRdxi,ddRddxi,
     &           MCRD,NNODE,COORDS)
         
      
c     Computing local Cartesian basis
            eI(:,:) = zero
            call norm(AI(1,:), 3, normV)
            eI(1,:) = AI(1,:)/normV
            
            call dot(AI(2,:), eI(1,:), normV)
            eI(2,:) = AI(2,:) - normV*eI(1,:)
            call norm(eI(2,:), 3, normV)
            eI(2,:) = eI(2,:)/normV
            
            eI(3,:) = AI(3,:)
            
            
c     Get covariant basis
            AE(:,:) = zero
            AE(1,:) = AAE(1,1)*AI(1,:) + AAE(1,2)*AI(2,:)
            AE(2,:) = AAE(2,1)*AI(1,:) + AAE(2,2)*AI(2,:)
            AE(3,:) = AI(3,:)
            
            
c     Computing strain in curvilinear system
            call uStrainMem_shell(sol,NNODE,MCRD,AI,dRdxi,stran(1:3))
            call uStrainBnd_shell(sol,NNODE,MCRD,AI,dAI1dxi,dAI2dxi,
     &           dRdxi,ddRddxi,stran(4:6))
            
            
c     Get strain in local cartesian basis
            e_AE(:,:) = zero
            do i = 1,2
               do j = 1,2
                  call dot(eI(i,:), AE(j,:), e_AE(i,j))
               enddo
            enddo
            
            Pmtx(:,:) = zero
            Pmtx(1,1) = e_AE(1,1)*e_AE(1,1)
            Pmtx(1,2) = e_AE(1,2)*e_AE(1,2)
            Pmtx(1,3) = e_AE(1,1)*e_AE(1,2)
            Pmtx(2,1) = e_AE(2,1)*e_AE(2,1)
            Pmtx(2,2) = e_AE(2,2)*e_AE(2,2)
            Pmtx(2,3) = e_AE(2,1)*e_AE(2,2)
            Pmtx(3,1) = e_AE(1,1)*e_AE(2,1)*two
            Pmtx(3,2) = e_AE(1,2)*e_AE(2,2)*two
            Pmtx(3,3) =(e_AE(1,1)*e_AE(2,2)+e_AE(1,2)*e_AE(2,1))
            call MulVect(Pmtx, stran(1:3), stranC(1:3), 3, 3)
            call MulVect(Pmtx, stran(4:6), stranC(4:6), 3, 3)
            
c     Get stress
            If (Output_FLAG(2)) then
               matH(:,:) = zero
               matH(1,1) = one; matH(1,2) = nu
               matH(2,1) = nu;  matH(2,2) = one
               matH(3,3) = 0.5D0*(one-nu)
               matH(:,:) = coef*matH(:,:)
               
               stressC(:) = zero
               matH(:,:) = h*matH(:,:)
               call MulVect(matH, stranC(1:3), stressC(1:3), 3, 3)
               
               matH(:,:) = h**two/12.0D0 * matH(:,:)
               call MulVect(matH, stranC(4:6), stressC(4:6), 3, 3)
            Endif
            
         Endif
         
c     Sum up all variables into svarsip
         svarsip = zero
         offset = 1
         
         svarsip(offset:offset+2) = coords_ip(:)
         offset = offset + 3
         
         If (Output_FLAG(1)) then
            svarsip(offset:offset+2) =  u_ip(:) !w_ip(:)
            offset = offset + 3
         Endif
         
         If (Output_FLAG(2)) then
            !temp = sqrt(rot(1)**two+rot(2)**two+rot(3)**two)
            !svarsip(offset:offset+2) = rot(:)
            !offset = offset + 3
            !svarsip(offset)   = temp
            !svarsip(offset+1) = alpha
            !svarsip(offset+2) = one
            !offset = offset + 3
            svarsip(offset:offset+ntens-1) = stressC(:)
            offset = offset + ntens
         Endif
         
         If (Output_FLAG(3)) then
            svarsip(offset:offset+ntens-1) = stranC(:)
         Endif
         
         
                  
c     Update global variable : all variables at each intergration point
         do i = 1,nsvint
            svars(nsvint*(n-1)+i) = svarsip(i)
         enddo
      Enddo
      Enddo
      
C     ------------------------------------------------------------------
      
      End subroutine compute_svars_Q1_shell
