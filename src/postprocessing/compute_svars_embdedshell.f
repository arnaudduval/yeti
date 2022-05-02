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

            
      
C     ******************************************************************
      
C     Calcul des deformations, contraintes et deplacements au niveau des
c     frontiere delimitant un element nurbs. Cela permettra de creer le
c     fichier VTU pour la visualisation.
      subroutine compute_svars_Q1_embdedshell(COORDS,COORDSall,sol,
     1     svars,nsvint,Output_FLAG,nb_vertice,nb_REF,MCRD,NNODE,
     2     NNODEmap,nb_cp,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS)
      
      use parameters
      use embeddedMapping

      Implicit None
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
      Double precision, intent(in) :: COORDS,COORDSall,sol,
     &     MATERIAL_PROPERTIES,PROPS
      dimension COORDS(MCRD,NNODE)
      dimension COORDSall(3,nb_cp)
      dimension sol(MCRD,NNODE)
      dimension MATERIAL_PROPERTIES(2)
      dimension PROPS(JPROPS)
      Character(len=*) , intent(in) :: TENSOR
      
      Integer, intent(in) :: nsvint,MCRD,NNODE,NNODEmap,nb_cp,
     &     nb_vertice,JPROPS,nb_REF
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
      
!     Embedded Surface
      ! - nurbs basis functions
      Double precision :: Rs, dRdTheta, ddRddTheta, DetJac
      dimension Rs(NNODE), dRdTheta(NNODE,2), ddRddTheta(NNODE,3)
      ! - curvilinear coordinate objects
      Double precision :: XI,BI,dBI
      dimension XI(3),BI(3,2),dBI(3,3)
      
!     Mapping
      ! - nurbs basis functions
      Double precision :: R, dRdxi, ddRddxi
      dimension R(NNODEmap), dRdxi(NNODEmap,3), ddRddxi(NNODEmap,6)
      ! - element infos
      Double precision :: COORDSmap
      dimension COORDSmap(MCRD,NNODEmap)
      Integer          :: sctr_map
      dimension sctr_map(NNODEmap)
      
!     Composition Mapping+Surface
      Double precision :: dRRdTheta, ddRRddTheta, AI, dAI1dTheta,
     &     dAI2dTheta, AAE
      dimension dRRdTheta(NNODEmap,2),ddRRddTheta(NNODEmap,3),AI(3,3),
     &     dAI1dTheta(2,3), dAI2dTheta(2,3), AAE(2,2)

!     Cartesian basis vectors
      Double precision :: eI, normV, AE
      dimension eI(3,3),AE(3,3)

!     For material matrix
      Double precision :: E, nu, h, matH, coef
      dimension matH(3,3)
      
!     Other
      Double precision :: stran, svarsip, coords_ip,u_ip,stranC,stressC,
     &     e_AE, Pmtx
      dimension stran(2*MCRD), svarsip(nsvint),coords_ip(3), u_ip(3),
     &     stranC(2*MCRD), stressC(2*MCRD), e_AE(2,2), Pmtx(3,3)
      
      Integer :: isave
      Integer :: n,k1,ntens,i,j,kk,nb_xi,nb_eta,i_xi,i_eta,offset
      

C     ------------------------------------------------------------------
      
C     Initialization :
c     --------------
      
      ntens = 6
      
c     Get material behaviour tensor
      h = PROPS(3)
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

      isave = 0
      
      svars  = zero
      nb_xi  = 2**max(nb_REF(1)-1,0)+1
      nb_eta = 2**max(nb_REF(2)-1,0)+1

      Do i_eta = 1,nb_eta
      Do i_xi  = 1,nb_xi
         
         n = (i_eta-1)*nb_xi + i_xi
         vertice(1,n) = two/dble(nb_xi -1)*dble(i_xi -1) - one
         vertice(2,n) = two/dble(nb_eta-1)*dble(i_eta-1) - one
         
c     Computing NURBS basis functions and derivatives
         DetJac = zero
         Rs(:)  = zero
         dRdTheta(:,:)   = zero
         ddRddTheta(:,:) = zero
         call nurbsbasis(Rs,dRdTheta,ddRddTheta,DetJac,vertice(1:,n))
         
c     Find mapping parametric position
         XI(:)    = zero
         BI(:,:)  = zero
         dBI(:,:) = zero
         Do k1 = 1,NNODE
            XI(:)    =  XI(:)   +           Rs(k1)*COORDS(:,k1)
            BI(:,1)  =  BI(:,1) +   dRdTheta(k1,1)*COORDS(:,k1)
            BI(:,2)  =  BI(:,2) +   dRdTheta(k1,2)*COORDS(:,k1)
            dBI(:,1) = dBI(:,1) + ddRddTheta(k1,1)*COORDS(:,k1)
            dBI(:,2) = dBI(:,2) + ddRddTheta(k1,2)*COORDS(:,k1)
            dBI(:,3) = dBI(:,3) + ddRddTheta(k1,3)*COORDS(:,k1)
         Enddo
         
c     Computing NURBS basis functions and derivatives of the mapping
         ! get active element number
         call updateMapElementNumber(XI(:))
         
         call evalnurbs_mapping_w2ndDerv(XI(:),R(:),dRdxi(:,:),
     &        ddRddxi(:,:))
         
         ! extract COORDS
         If (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:,current_map_elem)
            
            Do k1 = 1,NNODEmap
               COORDSmap(:,k1) = COORDSall(:,sctr_map(k1))
            Enddo
            
            isave = current_map_elem
         Endif
         
c     Get intergration points coordinates and displacements
         coords_ip = zero
         Do k1 = 1,NNODEmap
            coords_ip(:MCRD) = coords_ip(:MCRD)+R(k1)*COORDSmap(:,k1)
         Enddo
         
         u_ip = zero
         if (Output_FLAG(1)) then
            Do k1 = 1,NNODE   
               u_ip(:MCRD) = u_ip(:MCRD) + Rs(k1)*sol(:,k1)
            Enddo
         endif
         
         
         
         
         
c     Computing Curvilinear Coordinate objects
         If (Output_FLAG(2) .OR. Output_FLAG(3)) then

!     - Composition of the basis functions
            dRRdTheta(:,:) = zero
            Do i = 1,3
               dRRdTheta(:,1) = dRRdTheta(:,1) + BI(i,1)*dRdxi(:,i)
               dRRdTheta(:,2) = dRRdTheta(:,2) + BI(i,2)*dRdxi(:,i)
            Enddo
            
            ddRRddTheta(:,:) = zero
            Do i = 1,3
               ddRRddTheta(:,1) = ddRRddTheta(:,1)
     &              + dBI(i,1)*dRdxi(:,i) + BI(i,1)*BI(i,1)*ddRddxi(:,i)
               ddRRddTheta(:,2) = ddRRddTheta(:,2) 
     &              + dBI(i,2)*dRdxi(:,i) + BI(i,2)*BI(i,2)*ddRddxi(:,i)
               ddRRddTheta(:,3) = ddRRddTheta(:,3) 
     &              + dBI(i,3)*dRdxi(:,i) + BI(i,1)*BI(i,2)*ddRddxi(:,i)
               Do j = i+1,3
                  kk = i+j+1
                  ddRRddTheta(:,1) = ddRRddTheta(:,1)
     &                 + two*BI(i,1)*BI(j,1)*ddRddxi(:,kk)
                  ddRRddTheta(:,2) = ddRRddTheta(:,2)
     &                 + two*BI(i,2)*BI(j,2)*ddRddxi(:,kk)
                  ddRRddTheta(:,3) = ddRRddTheta(:,3)
     &                 + two*BI(i,1)*BI(j,2)*ddRddxi(:,kk)
               Enddo
            Enddo
            
!     - curvilinear vectors and derivatives
            call curvilinear(AI,dAI1dTheta,dAI2dTheta,AAE,R,dRRdTheta,
     &        ddRRddTheta,MCRD,NNODEmap,COORDSmap)
            
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
            call uStrainMem_shell(sol,NNODE,MCRD,AI,dRdTheta,stran(1:3))
            call uStrainBnd_shell(sol,NNODE,MCRD,AI,dAI1dTheta,
     &           dAI2dTheta,dRdTheta,ddRddTheta,stran(4:6))
            
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
            Pmtx(3,3) =(e_AE(1,1)*e_AE(2,2)+e_AE(1,2)*e_AE(2,1)) !*two
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
            svarsip(offset:offset+2) =  u_ip(:)
            offset = offset + 3
         Endif
         
         If (Output_FLAG(2)) then
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
      
      End subroutine compute_svars_Q1_embdedshell
