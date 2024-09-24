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
      
      
      subroutine gradComp_semiAN(gradC, SOL,epsilon,
     1     activeElement,activeDir,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,
     2     Nijk,weight,Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,
     3     TENSOR,indDLoad,JDLType,ADLMAG,load_target_nbelem,MCRD,NBINT,
     4     nb_load,nb_patch,nb_elem,nnode,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
      
      Implicit none
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer, intent(in)          :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in)          :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in)          :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN,nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),
     &     PROPS(:),
     &     NNODE(nb_patch),
     &     IEN(:),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),
     &     NBINT(nb_patch)
      
      
!     Loads
      Double precision, intent(in) :: ADLMAG
      Integer, intent(in)          :: nb_load,indDLoad,JDLType,
     &     load_target_nbelem
      dimension ADLMAG(nb_load),
     &     indDLoad(:),
     &     JDLType(nb_load),
     &     load_target_nbelem(nb_load)
      
!     INFOS
      Integer, intent(in)          :: activeElement, activeDir
      dimension activeElement(nb_elem), activeDir(3)
      
!     Solution
      Double precision, intent(in) :: epsilon,SOL
      dimension SOL(3,nb_cp)
      
c     Output variables : system lineaire a resoudre
c     ----------------
      Double precision, intent(out):: gradC
      dimension gradC(3,nb_cp)
            
      
c     Local variables :
c     ---------------
      
!     UELMAT.f
      Integer :: NDOFEL
      Double precision :: COORDS_elem, RHS0, AMATRX0, MAT_patch
      dimension COORDS_elem(MCRD,MAXVAL(NNODE)),
     &     RHS0(MCRD*MAXVAL(NNODE)),
     &     AMATRX0(MCRD,MCRD,MAXVAL(NNODE)*(MAXVAL(NNODE)+1)/2),
     &     MAT_patch(2)
      
!     FD approx
      Integer          :: dir
      Double precision :: gradC_elem,RHS,AMATRX,coef,savecp,gC
      dimension gradC_elem(3,MAXVAL(NNODE)),
     &     RHS(MCRD*MAXVAL(NNODE)),
     &     AMATRX(MCRD,MCRD,MAXVAL(NNODE)*(MAXVAL(NNODE)+1)/2)
      
!     Solution
      Double precision :: U_elem
      dimension U_elem(3,MAXVAL(NNODE))

!     Assembly
      Integer :: num_elem,numcp,i,j,kk,JELEM,Numpatch,sctr,nnodeSum,cpi,
     &     cpj,cpij
      dimension sctr(MAXVAL(NNODE))
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Compute gradient .................................................
      

      gradC(:,:) = zero
      
      
      
      JELEM = 0
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         If (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif
         
         NDOFEL = nnode_patch*MCRD
         nnodeSum = nnode_patch*(nnode_patch+1)/2
         
c     Loop on element
         Do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            
            IF (activeElement(JELEM)==1) then

c     Get element infos
            CALL extractNurbsElementInfos(num_elem)
            sctr(:nnode_patch)  = IEN_patch(:,num_elem)
            
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               U_elem(:,i)      = SOL(:,sctr(i))
            Enddo
            
            
c     Compute initial elementary matrix and load vector
            RHS0    = zero
            AMATRX0 = zero
            MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
            if (ELT_TYPE_patch == 'U1') then
               ! 'Element classique solide'
               call UELMAT_byCP(NDOFEL,MCRD,nnode_patch,JELEM,
     1              NBINT(NumPatch),COORDS_elem(:,:nnode_patch),
     2              TENSOR_patch,MAT_patch,nb_load,indDLoad,
     3              load_target_nbelem,JDLType,ADLMAG,RHS0(:NDOFEL),
     4              AMATRX0(:,:,:nnodeSum))
            elseif (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               call UELMAT3_byCP(NDOFEL,MCRD,nnode_patch,JELEM,
     1              NBINT(NumPatch),COORDS_elem,TENSOR_patch,MAT_patch,
     2              PROPS_patch,JPROPS_patch,nb_load,indDLoad,
     3              load_target_nbelem,JDLType,ADLMAG,RHS0(:NDOFEL),
     4              AMATRX0(:,:,:nnodeSum))
            elseif (ELT_TYPE_patch == 'U30') then
               ! 'Element coque immerge'
               call UELMAT30(NDOFEL,MCRD,nnode_patch,nnode_map,nb_cp,
     1              JELEM,NBINT(NumPatch),COORDS_elem,COORDS3D,
     2              TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,
     3              nb_load,indDLoad,load_target_nbelem,JDLType,ADLMAG,
     4              RHS0(:NDOFEL),AMATRX0(:,:,:nnodeSum))
            endif
            
            
c     Finite Difference Approximation
            gradC_elem(:,:) = zero
            Do numcp = 1,nnode_patch
               Do dir = 1,3
               If (activeDir(dir) == 1) then
                  
               savecp = COORDS_elem(dir,numcp)
               COORDS_elem(dir,numcp) = savecp + epsilon
               
               RHS(:) = zero
               AMATRX(:,:,:) = zero
               if (ELT_TYPE_patch == 'U1') then
                  ! 'Element classique solide'
                  call UELMAT_byCP(NDOFEL,MCRD,nnode_patch,JELEM,
     1                 NBINT(NumPatch),COORDS_elem(:,:nnode_patch),
     2                 TENSOR_patch,MAT_patch,nb_load,indDLoad,
     3                 load_target_nbelem,JDLType,ADLMAG,RHS(:NDOFEL),
     4                 AMATRX(:,:,:nnodeSum))
               elseif (ELT_TYPE_patch == 'U3') then
                  ! 'Element coque'
                  call UELMAT3_byCP(NDOFEL,MCRD,nnode_patch,JELEM,
     1               NBINT(NumPatch),COORDS_elem,TENSOR_patch,MAT_patch,
     2                 PROPS_patch,JPROPS_patch,nb_load,indDLoad,
     3                 load_target_nbelem,JDLType,ADLMAG,RHS(:NDOFEL),
     4                 AMATRX(:,:,:nnodeSum))
               elseif (ELT_TYPE_patch == 'U30') then
                  ! 'Element coque immerge'
                  call UELMAT30(NDOFEL,MCRD,nnode_patch,nnode_map,nb_cp,
     1                 JELEM,NBINT(NumPatch),COORDS_elem,COORDS3D,
     2                 TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,
     3                 nb_load,indDLoad,load_target_nbelem,JDLType,
     4                 ADLMAG,RHS(:NDOFEL),AMATRX(:,:,:nnodeSum))
                  
                  COORDS_elem(dir,numcp) = savecp - epsilon
                  call UELMAT30(NDOFEL,MCRD,nnode_patch,nnode_map,nb_cp,
     1                 JELEM,NBINT(NumPatch),COORDS_elem,COORDS3D,
     2                 TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,
     3                 nb_load,indDLoad,load_target_nbelem,JDLType,
     4                 ADLMAG,RHS0(:NDOFEL),AMATRX0(:,:,:nnodeSum))
               endif
               

               gC = zero
               
               cpij = 0
               j = 0
               Do cpj = 1,nnode_patch
                  Do cpi = 1,cpj
                     coef = 1.0d0
                     if (cpi==cpj) coef = 0.5d0
                     cpij = cpij + 1
                     
                     Do i = 1,3
                        gC = gC 
     &                       - coef*U_elem(i,cpj)*SUM(U_elem(:,cpi)*(
     &                       AMATRX(:,i,cpij) - AMATRX0(:,i,cpij)))
                     Enddo
                  Enddo
                  
                  Do i = 1,3
                     j = j +1
                     gC = gC + (RHS(j)-RHS0(j))* U_elem(i,cpj)
                  Enddo
                  
               Enddo
                              
               gradC_elem(dir,numcp) = gC/epsilon/two
                              
               COORDS_elem(dir,numcp) = savecp
               
               Endif ! active dir
               Enddo ! end loop on direction
            Enddo ! end loop on cp
            
c     Update global grad
            Do numcp = 1,nnode_patch
               gradC(:,sctr(numcp)) = gradC(:,sctr(numcp))
     &              + gradC_elem(:,numcp)
            Enddo
            
            Endif ! active elem
            
         Enddo ! end loop on element
         
         call deallocateMappingData()
         call finalizeNurbsPatch()
         
      Enddo ! end loop on patch

c     Fin calcul .......................................................
      
      end subroutine gradComp_semiAN
