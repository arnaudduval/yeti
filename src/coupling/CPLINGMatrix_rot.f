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
           
      
      
C     ******************************************************************
      
      
      Subroutine build_CplgMatrix_rot(CMatrx, ind_dof_free,nb_dof_free,
     1     nb_dof_tot,nb_gauss_tot,IGJelem,GaussCoordsGlo,
     2     COORDS3D,  IEN,  nb_elem_patch,  Nkv,  Ukv,  Nijk,  weight,
     3     Jpqr,  ELT_TYPE,  PROPS,  JPROPS,  TENSOR,  NNODE,  nb_cp,
     4     nb_patch,  nb_elem,  MCRD,
     5     COORDS3D_l,IEN_l,nb_elem_patch_l,Nkv_l,Ukv_l,Nijk_l,weight_l,
     6     Jpqr_l,ELT_TYPE_l,PROPS_l,JPROPS_l,TENSOR_l,NNODE_l,nb_cp_l,
     7     nb_patch_l,nb_elem_l,NBINT,dim_interface)
      
      use parameters
      use nurbspatch      

      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     **
!     Geometry NURBS
      Integer,          intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,nb_elem_patch,
     &     IEN,JPROPS
      dimension NNODE(nb_patch),nb_elem_patch(nb_patch),IEN(:),PROPS(:),
     &     JPROPS(nb_patch)
      
!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)
      
      
!     **
!     Multiplicateur de Lagrange
      Integer,          intent(in) :: nb_cp_l
      Double precision, intent(in) :: COORDS3D_l
      dimension COORDS3D_l(3,nb_cp_l)
      
      Double precision, intent(in) :: Ukv_l, weight_l
      Integer, intent(in) :: Nkv_l, Jpqr_l, Nijk_l
      dimension Nkv_l(3,nb_patch_l), Jpqr_l(3,nb_patch_l), 
     &     Nijk_l(3,nb_elem_l),Ukv_l(:),weight_l(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR_l, ELT_TYPE_l
      Double precision, intent(in) :: PROPS_l
      Integer, intent(in) :: NNODE_l,nb_patch_l,nb_elem_l,
     &     nb_elem_patch_l,IEN_l,JPROPS_l
      dimension NNODE_l(nb_patch_l),nb_elem_patch_l(nb_patch_l),
     &     IEN_l(:),PROPS_l(:),JPROPS_l(nb_patch_l)
      
      
!     Gaussian points
      Integer, intent(in) :: nb_gauss_tot,IGJelem,dim_interface,NBINT
      dimension NBINT(nb_patch_l)
      Double precision, intent(in) :: GaussCoordsGlo
      dimension IGJelem(nb_gauss_tot), GaussCoordsGlo(3,nb_gauss_tot)
            
      
c     Output variables : operateur de couplage
c     ----------------
      Double precision, intent(out) :: CMatrx
      dimension CMatrx(nb_dof_free,nb_cp_l)
!     + nb_elem_l)
      
      
      
      
c     Local Variables :
c     ---------------
      
      ! gauss info
      Integer          :: num_gauss,NbPtInt,IGNumPatch
      Double precision :: PtGauss,xi,GaussPdsCoord
      dimension PtGauss(dim_interface),xi(3),
     &     GaussPdsCoord(dim_interface+1,MAXVAL(NBINT)),
     &     IGNumPatch(nb_gauss_tot)
      
      ! computing Nurbs basis functions
      Double precision :: R_l,dRdxi_l,R,dRdxi
      dimension R_l(MAXVAL(NNODE_l)),dRdxi_l(MAXVAL(NNODE_l),3),
     &     R(MAXVAL(NNODE)),dRdxi(MAXVAL(NNODE),3)
      
      ! jacobian
      Double precision :: A1,A2,v,COORDS_l,DetJac,norm,dvol
      dimension A1(3),A2(3),v(3),COORDS_l(3,MAXVAL(NNODE_l))
      
      ! for rotation field
      Double precision :: A3,COORDS_elem,A3_A1,A2_A3,vectT,vectN,coef1,
     &     coef2,coef3,coef4
      dimension A3(3),COORDS_elem(3,MAXVAL(NNODE)),A3_A1(3),A2_A3(3),
     &     vectT(3),vectN(3)

      ! for assembly
      Integer :: num_l,numCP,ddl1,ddl2,k1, sctr,sctr_l,tab
      Integer :: numPatch_l,numPatch,num_elem_l,num_elem
      dimension sctr(MAXVAL(NNODE)),sctr_l(MAXVAL(NNODE_l)),
     &     tab(nb_dof_tot)
      
      ! for local loops
      Integer          :: i,n,JELEM_l,JELEM
      Double precision :: temp
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
      
!     Initialize coupling matrix
      CMatrx(:,:) = zero
      
!     Create tab to link dof_tot <--> dof_free
!     tab(i) = | 0  if dof is blocked
!              | ii if dof is free (ii: index for the reduced matrix)
      tab(:) = 0
      Do i = 1,nb_dof_free
         tab(ind_dof_free(i)) = i
      Enddo
      
!     get IGNumPatch
      IGNumPatch(:) = 1
      JELEM = nb_elem_patch(1)
      Do numPatch = 2,nb_patch
         JELEM = JELEM + nb_elem_patch(numPatch)
         Do n = 1,nb_gauss_tot
            If (IGJelem(n)<JELEM) then
               IGNumPatch(n) = numPatch
            Endif
         Enddo
      Enddo
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      num_gauss = 0
      JELEM_l   = 0
c     Loop on patches
      Do NumPatch_l = 1,nb_patch_l
         
!     Gauss points
         NbPtInt=int(
     &        float(NBINT(NumPatch_l))**(one/float(dim_interface)))
         if (NbPtInt**dim_interface<NBINT(NumPatch_l)) NbPtInt=NbPtInt+1
         call Gauss(NbPtInt,dim_interface,GaussPdsCoord,0)
         
         
c     Loop on elements
         Do num_elem_l = 1,nb_elem_patch_l(NumPatch_l)
            JELEM_l = JELEM_l + 1
            
            CALL extractNurbsPatchGeoInfos(NumPatch_l, Nkv_l,Jpqr_l,
     &           Nijk_l,Ukv_l,weight_l,nb_elem_patch_l)
            CALL extractNurbsPatchMechInfos(NumPatch_l,IEN_l,PROPS_l,
     &           JPROPS_l,NNODE_l,nb_elem_patch_l,ELT_TYPE_l,TENSOR_l)
            
            
            COORDS_l(:,:) = zero
            Do i = 1,NNODE_l(NumPatch_l)
               COORDS_l(:,i) = COORDS3D_l(:,IEN_patch(i,num_elem_l))
            Enddo
            CALL extractNurbsElementInfos(num_elem_l)
            
c     Loop on gauss points .............................................
            Do n = 1,NBINT(numPatch_l)
               
c     
c     ...
c     
               
C     Compute non-zero nurbs basis functions for the Lagrange Multiplier
c     - exctract infos
               CALL extractNurbsPatchGeoInfos(NumPatch_l, Nkv_l,Jpqr_l,
     &              Nijk_l,Ukv_l,weight_l,nb_elem_patch_l)
               CALL extractNurbsPatchMechInfos(NumPatch_l,IEN_l,PROPS_l,
     &              JPROPS_l,NNODE_l,nb_elem_patch_l,ELT_TYPE_l,
     &              TENSOR_l)
               CALL extractNurbsElementInfos(num_elem_l)
               sctr_l(:nnode_patch) = IEN_patch(:,num_elem_l)
!     - get gauss coordinates
               PtGauss(:) = GaussPdsCoord(2:,n)
               xi(:) = zero
               Do i  = 1,dim_interface
                  xi(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
               Enddo
!     - Calculate univariate B-spline function
               call evalnurbs(xi,R_l(:nnode_patch),
     &              dRdxi_l(:nnode_patch,:))
               
!     Compute Jacobien
               A1(:)  = zero
               A2(:)  = zero
               DetJac = zero
               v(:)   = zero
               If (dim_interface == 1) then
!     - interface 1d
                  DetJac = 0.5d0*( Ukv_elem(2,1) - Ukv_elem(1,1) )
                  Do num_l = 1,NNODE_l(NumPatch_l) 
                     A1(:)  = A1(:) + dRdxi_l(num_l,1)*COORDS_l(:,num_l)
                  Enddo
                  v(:) = A1(:)
               Elseif (dim_interface == 2) then
!     - interface 2d
                  DetJac = ( Ukv_elem(2,1) - Ukv_elem(1,1) )*0.5d0
     &                   * ( Ukv_elem(2,2) - Ukv_elem(1,2) )*0.5d0
                  Do num_l = 1,NNODE_l(NumPatch_l)
                     A1(:)  = A1(:) + dRdxi_l(num_l,1)*COORDS_l(:,num_l)
                     A2(:)  = A2(:) + dRdxi_l(num_l,2)*COORDS_l(:,num_l)
                  Enddo
                  call cross(A1(:), A2(:), v)
               Endif
               norm = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) )
               dvol = GaussPdsCoord(1,n)*norm*DetJac
               
               vectT(:) = v(:)/norm
c     
c     ...
c     
               
               
c     Extract element infos of the geometry
               num_gauss = num_gauss+1
               xi(:)     = GaussCoordsGlo(:,num_gauss)
               JELEM     = IGJelem(num_gauss)
               numPatch  = IGNumPatch(num_gauss)
               
               CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,
     &              Nijk,Ukv,weight,nb_elem_patch)
               CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,
     &              JPROPS,NNODE,nb_elem_patch,ELT_TYPE,
     &              TENSOR)
               num_elem = JELEM
               Do i = 1,NumPatch-1
                  num_elem = num_elem - nb_elem_patch(i)
               Enddo
               CALL extractNurbsElementInfos(num_elem)
               sctr(:nnode_patch) = IEN_patch(:,num_elem)
               
               COORDS_elem(:,:) = zero
               Do i = 1,NNODE(numPatch)
                  COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               Enddo
                              
!     Compute non-zero nurbs basis functions for the Geometry
               call evalnurbs(xi,R(:nnode_patch),dRdxi(:nnode_patch,:))
               
c     Compute rotation field
               A1 = zero
               A2 = zero
               Do numCP = 1,NNODE(numPatch)
                  A1(:)  = A1(:) + dRdxi(numCP,1)*COORDS_elem(:,numCP)
                  A2(:)  = A2(:) + dRdxi(numCP,2)*COORDS_elem(:,numCP)
               Enddo
               call cross(A1(:), A2(:), A3(:))
               norm = sqrt( A3(1)*A3(1) + A3(2)*A3(2) + A3(3)*A3(3) )
               A3(:) = A3(:)/norm
               
               call cross(A3(:),A1(:), A3_A1)
               call cross(A2(:),A3(:), A2_A3)
               call cross(A3(:),vectT(:),vectN(:))
               
               call dot(A3_A1(:),vectN(:), coef1)
               call dot(A2_A3(:),vectN(:), coef2)
               
               call dot(A3_A1(:),vectT(:), coef3)
               call dot(A2_A3(:),vectT(:), coef4)
 
c     
c     ...
c     
               
               
c     UPDATE coupling matrix
c     - rotation along the interface direction
               Do num_l = 1,NNODE_l(numPatch_l)
                  temp  = R_l(num_l)*dvol/norm
                  ddl2  = sctr_l(num_l)
                  Do numCP = 1,NNODE(numPatch)
                     k1    = sctr(numCP)
                     Do i  = 1,MCRD
                        ddl1  = tab((k1-1)*MCRD+i)
                        If (ddl1>0) then
                           CMatrx(ddl1,ddl2) = CMatrx(ddl1,ddl2)
     &                          - dRdxi(numCP,2)*A3(i)*coef1*temp
     &                          - dRdxi(numCP,1)*A3(i)*coef2*temp
                        Endif
                     Enddo
                  Enddo
               Enddo

c     - rotation along the normal to the interface direction
               If (.False.) then
               temp = dvol/norm
               ddl2 = JELEM_l + nb_cp_l
               Do numCP = 1,NNODE(numPatch)
                  k1    = sctr(numCP)
                  Do i  = 1,MCRD
                     ddl1  = tab((k1-1)*MCRD+i)
                     If (ddl1>0) then
                        CMatrx(ddl1,ddl2) = CMatrx(ddl1,ddl2)
     &                       - dRdxi(numCP,2)*A3(i)*coef3*temp
     &                       - dRdxi(numCP,1)*A3(i)*coef4*temp
                     Endif
                  Enddo
               Enddo
               Endif
               
            Enddo
!     End loop gauss points ............................................
         Enddo
!     End loop on Lagrange mult. elements
         
         CALL finalizeNurbsPatch()
         
         
      Enddo
!     End loop on Lagrange mult. patches
      
      
      End subroutine build_CplgMatrix_rot
