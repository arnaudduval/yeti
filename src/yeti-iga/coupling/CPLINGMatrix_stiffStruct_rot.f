!! Copyright 2016-2019 Thibaut Hirschler

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
      
C     ASSEMBLAGE DE LA MATRICE DE COUPLAGE PAR MULTIPLICATEUR DE 
C     LAGRANGE
c     --> cas des structures raidies
      
      subroutine build_CplgMatrix_rot(CMatrix, TABinterfaceDef,NBINT,
     1     nb_stiff,nb_cp_l,ind_dof_free,nb_dof_free,nb_dof_tot,
     2     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,
     3     MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,nb_patch,
     4     nb_elem,nb_cp,MCRD)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
!     Geometry NURBS
      Integer,          intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer,          intent(in) :: MCRD,NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(:),
     &     JPROPS(nb_patch),NNODE(nb_patch),IEN(:),
     &     nb_elem_patch(nb_patch)
      
!     Multiplicateur de lagrange
      Integer, intent(in) :: nb_cp_l,NBINT,nb_stiff,TABinterfaceDef
      dimension NBINT(nb_patch), TABinterfaceDef(2,nb_stiff)

!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)
      
      
!     Output arguments :
!     ----------------
      Double precision, intent(out):: CMatrix
      dimension CMatrix(nb_dof_free,nb_cp_l)
      
      
C     Local variables :
c     ---------------
      
      ! Gauss info
      Integer          :: NbPtInt,n
      Double precision :: xi_bar,Xi,GaussPdsCoord,PtGauss
      dimension xi_bar(3),Xi(3),GaussPdsCoord(2,MAXVAL(NBINT))
      
      ! Nurbs functions
      Double precision :: COORDS_elem,COORDS_l, R,dRdxi, R_l,dRdxi_l,
     &     AI,B1,normV,dvol,DetJac
      dimension COORDS_elem(MCRD,MAXVAL(NNODE)), 
     &     COORDS_l(MCRD,MAXVAL(NNODE)), 
     &     R(MAXVAL(NNODE)),  dRdxi(MAXVAL(NNODE),3),
     &     R_l(MAXVAL(NNODE)),dRdxi_l(MAXVAL(NNODE),3),
     &     AI(3,3),B1(3)
      
      ! for rotation field
      Double precision :: vectT,vectN,A3_A1,A2_A3, coef1,coef2, normA
      dimension vectT(3),vectN(3),A3_A1(3),A2_A3(3)

      ! for assembly
      Double precision :: temp
      Integer          :: numI,numC,numS,numP, tab,offset,sctr,sctr_l,
     &     num_elem,num_elem_l,numel,num_elemDir,nb_elemDir,nbel,
     &     k1,ddl1,ddl2,num_l,numCP, count
      dimension tab(nb_dof_tot),sctr(MAXVAL(NNODE)),
     &     sctr_l(MAXVAL(NNODE)),num_elemDir(3),nb_elemDir(3)

      ! for loops
      Integer :: i
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
      
      CMatrix(:,:) = zero
      
!     Create tab to link dof_tot <--> dof_free
!     tab(i) = | 0  if dof is blocked
!              | ii if dof is free (ii: index for the reduced matrix)
      tab(:) = 0
      Do i = 1,nb_dof_free
         tab(ind_dof_free(i)) = i
      Enddo
      
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      count = 0
      Do numI = 1,nb_stiff
         
!     Get patch numbers
         numC = TABinterfaceDef(1,numI) ! curve
         numS = TABinterfaceDef(2,numI) ! stiffener
         CALL extractNurbsPatchMechInfos(numC,IEN,PROPS,JPROPS,NNODE,
     &        nb_elem_patch,ELT_TYPE,TENSOR)
         numP = int(PROPS_patch(2))     ! panel

!     Initialize stiffener element numbering
         numel= 1
         CALL extractNurbsPatchGeoInfos(numS, Nkv,Jpqr,Nijk,Ukv,weight,
     &        nb_elem_patch)
         nbel = 0
         Do i = Jpqr_patch(1)+1,Nkv_patch(1)-Jpqr_patch(1)-1
            If (Ukv1_patch(i)<Ukv1_patch(i+1)) then
               nbel = nbel+1
            Endif
         Enddo
         
!     Initialize lagrangian ddl
         sctr_l(:) = 0
         n = NNODE(numC) + 1 + count
         Do i = 1,NNODE(numC)
            sctr_l(i) = n-i
         Enddo
         count = count + Nkv(1,numC) - ( Jpqr(1,numC) + 1 )
         
c     Get gauss rule
         GaussPdsCoord(:,:) = zero
         NbPtInt = NBINT(numC)
         call Gauss(NbPtInt,1,GaussPdsCoord(:2,:NbPtInt),0)
         
c     Loop on elements
         Do num_elem_l = 1,nb_elem_patch(numC)
            
!     Extract curve infos
            CALL extractNurbsPatchGeoInfos(numC, Nkv,Jpqr,Nijk,Ukv,
     &           weight,nb_elem_patch)
            CALL extractNurbsPatchMechInfos(numC,IEN,PROPS,JPROPS,NNODE,
     &           nb_elem_patch,ELT_TYPE,TENSOR)
            
            COORDS_l(:,:) = zero
            Do i = 1,nnode_patch
               COORDS_l(:,i) = COORDS3D(:,IEN_patch(i,num_elem_l))
            Enddo
            offset = Nijk_patch(1,num_elem_l) - Jpqr_patch(1) - 1
            
c     Loop on gauss points .............................................
            Do n = 1,NbPtInt
               
c     Step1: Evaluate on curve
!     - extract infos
               CALL extractNurbsPatchGeoInfos(numC, Nkv,Jpqr,Nijk,Ukv,
     &              weight,nb_elem_patch)
               CALL extractNurbsPatchMechInfos(numC,IEN,PROPS,JPROPS,
     &              NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
               CALL extractNurbsElementInfos(num_elem_l)
!     - get gauss cordinates
               PtGauss   = GaussPdsCoord(2,n)
               xi_bar(:) = zero
               xi_bar(1) = ((Ukv_elem(2,1) - Ukv_elem(1,1)) * PtGauss 
     &                   +  (Ukv_elem(2,1) + Ukv_elem(1,1)) ) * 0.5d0
               DetJac = 0.5d0*( Ukv_elem(2,1) - Ukv_elem(1,1) )
!     - Calculate univariate B-spline function
               call evalnurbs(xi_bar, R_l(:nnode_patch),
     &              dRdxi_l(:nnode_patch,:))
               
!     Compute directional vector and parametric position
               B1(:) = zero
               Xi(:) = zero
               Do i = 1,nnode_patch
                  B1(:2) = B1(:2) + dRdxi_l(i,1)*COORDS_l(:2,i)
                  Xi(:2) = Xi(:2) + R_l(i)*COORDS_l(:2,i)
               Enddo
               
               
c     Step2: Evaluate on panel
c     - get panel infos
               CALL extractNurbsPatchGeoInfos(numP, Nkv,Jpqr,Nijk,Ukv,
     &              weight,nb_elem_patch)
               CALL extractNurbsPatchMechInfos(numP,IEN,PROPS,JPROPS,
     &              NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
               
c     - get position of the points (return element number)
               nb_elemDir  = 0
               num_elemDir = 1
               ! xi
               Do i = Jpqr_patch(1)+1,Nkv_patch(1)-Jpqr_patch(1)-1
                  If (Ukv1_patch(i)<Ukv1_patch(i+1)) then
                     nb_elemDir(1) = nb_elemDir(1)+1
                     If (Xi(1)>=Ukv1_patch(i)) then
                        num_elemDir(1) = nb_elemDir(1)
                     Endif
                  Endif
               Enddo
              ! eta
               Do i = Jpqr_patch(2)+1,Nkv_patch(2)-Jpqr_patch(2)-1
                  If (Ukv2_patch(i)<Ukv2_patch(i+1)) then
                     nb_elemDir(2) = nb_elemDir(2)+1
                     If (Xi(2)>=Ukv2_patch(i)) then
                        num_elemDir(2) = nb_elemDir(2)
                     Endif
                  Endif
               Enddo
               num_elem =(num_elemDir(2)-1)*nb_elemDir(1)+num_elemDir(1)
               
c     - get element infos
               sctr(:nnode_patch) = IEN_patch(:,num_elem)
               Do i = 1,nnode_patch
                  COORDS_elem(:,i)=COORDS3D(:MCRD,IEN_patch(i,num_elem))
               Enddo
               CALL extractNurbsElementInfos(num_elem)
               
c     - Calculate univariate B-spline function
               call evalnurbs(Xi,R(:nnode_patch),dRdxi(:nnode_patch,:))
               
c     - compute Jacobien
               AI(:,:) = zero
               Do i = 1,nnode_patch
                  AI(:,1) = AI(:,1) + dRdxi(i,1)*COORDS_elem(:,i)
                  AI(:,2) = AI(:,2) + dRdxi(i,2)*COORDS_elem(:,i)
               Enddo
               call norm(B1(1)*AI(:,1) + B1(2)*AI(:,2), 3, normV)
               dvol = GaussPdsCoord(1,n)*normV*DetJac
               
c     - compute rotation field
               vectT(:) = (B1(1)*AI(:,1) + B1(2)*AI(:,2))/normV
               
               call cross(AI(:,1),AI(:,2),AI(:,3))
               call norm(AI(:,3),3,normA)
               AI(:,3) = AI(:,3)/normA
               
               call cross(AI(:,3),AI(:,1), A3_A1)
               call cross(AI(:,2),AI(:,3), A2_A3)
               call cross(AI(:,3),vectT(:),vectN(:))
               
               call dot(A3_A1(:),vectN(:), coef1)
               call dot(A2_A3(:),vectN(:), coef2)
               

c     Update coupling matrix (panel side)
               Do num_l = 1,NNODE(numC)
                  temp  = R_l(num_l)*dvol/normA
                  ddl2  = sctr_l(num_l) + offset
                  Do numCP = 1,NNODE(numP)
                     k1    = sctr(numCP)
                     Do i  = 1,MCRD
                        ddl1  = tab((k1-1)*MCRD+i)
                        If (ddl1>0) then
                           CMatrix(ddl1,ddl2) = CMatrix(ddl1,ddl2)
     &                          - dRdxi(numCP,2)*AI(i,3)*coef1*temp
     &                          - dRdxi(numCP,1)*AI(i,3)*coef2*temp
                        Endif
                     Enddo
                  Enddo
               Enddo
               
               
c     Step3: Evaluate on stiffener
c     - get panel infos
               CALL extractNurbsPatchGeoInfos(numS, Nkv,Jpqr,Nijk,Ukv,
     &              weight,nb_elem_patch)
               CALL extractNurbsPatchMechInfos(numS,IEN,PROPS,JPROPS,
     &              NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
               CALL extractNurbsElementInfos(numel)
c     - update element number
               Do while ((xi_bar(1) > Ukv_elem(2,1)).AND.(numel<nbel+1))
                  numel = numel+1
                  CALL extractNurbsElementInfos(numel)
               Enddo
               sctr(:nnode_patch) = IEN_patch(:,numel)
               Do i = 1,nnode_patch
                  COORDS_elem(:,i) = COORDS3D(:MCRD,sctr(i))
               Enddo
c     - Calculate univariate B-spline function
               Xi(:) = (/ xi_bar(1), Ukv2_patch(1), zero /)
               call evalnurbs(Xi,R(:nnode_patch),dRdxi(:nnode_patch,:))
c     - compute rotation field
               AI(:,:) = zero
               Do i = 1,nnode_patch
                  AI(:,1) = AI(:,1) + dRdxi(i,1)*COORDS_elem(:,i)
                  AI(:,2) = AI(:,2) + dRdxi(i,2)*COORDS_elem(:,i)
               Enddo
               call cross(AI(:,1),AI(:,2),AI(:,3))
               call norm(AI(:,3),3,normA)
               AI(:,3) = AI(:,3)/normA

               call cross(AI(:,3),AI(:,1), A3_A1)
               call cross(AI(:,2),AI(:,3), A2_A3)
               call cross(AI(:,3),vectT(:),vectN(:))
               
               call dot(A3_A1(:),vectN(:), coef1)
               call dot(A2_A3(:),vectN(:), coef2)
               
c     Update coupling matrix (stiffener side)
               Do num_l = 1,NNODE(numC)
                  temp  = R_l(num_l)*dvol/normA
                  ddl2  = sctr_l(num_l) + offset
                  Do numCP = 1,NNODE(numS)
                     k1    = sctr(numCP)
                     Do i  = 1,MCRD
                        ddl1  = tab((k1-1)*MCRD+i)
                        If (ddl1>0) then
                           CMatrix(ddl1,ddl2) = CMatrix(ddl1,ddl2)
     &                          + dRdxi(numCP,2)*AI(i,3)*coef1*temp
     &                          + dRdxi(numCP,1)*AI(i,3)*coef2*temp
                        Endif
                     Enddo
                  Enddo
               Enddo
               
            Enddo ! int points
         Enddo ! curve elts
         
         CALL finalizeNurbsPatch()
           
      Enddo ! stiffener
      
      
C     ------------------------------------------------------------------
      
      end subroutine build_CplgMatrix_rot
