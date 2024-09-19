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

      
      
      
C     ******************************************************************
      
      
      subroutine gradVibration_AN(gradV, VECT, VALS, nb_frq,
     1     activeElement,activeDir,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,
     2     Nijk,weight,Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,
     3     RHO,TENSOR,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp)
      
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
      Double precision, intent(in) :: MATERIAL_PROPERTIES,RHO,PROPS
      Integer, intent(in)          :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN,nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),
     &     RHO(nb_patch),
     &     PROPS(:),
     &     NNODE(nb_patch),
     &     IEN(:),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),
     &     NBINT(nb_patch)
      
      
!     INFOS
      Integer, intent(in)          :: activeElement, activeDir
      dimension activeElement(nb_elem), activeDir(3)
      
!     Solution (eigenvectors and eigenvalues)
      Integer,          intent(in) :: nb_frq
      Double precision, intent(in) :: VECT,VALS
      dimension VECT(nb_frq,3,nb_cp),VALS(nb_frq)
      
c     Output variables :
c     ----------------
      Double precision, intent(out):: gradV
      dimension gradV(nb_frq,3,nb_cp)
      
      
c     Local variables :
c     ---------------
      
!     gradUELMASSMAT.f
      Integer :: NDOFEL,dir
      Double precision :: COORDS_elem,gradV_elem,MAT_patch,massElem
      dimension COORDS_elem(3,MAXVAL(NNODE)),MAT_patch(2),
     &     gradV_elem(nb_frq,3,MAXVAL(NNODE)),massElem(nb_frq)
      
!     Solution
      Double precision :: U_elem
      dimension U_elem(nb_frq,3,MAXVAL(NNODE))
      
!     Assembly
      Double precision :: massTot
      dimension massTot(nb_frq)
      Integer :: num_elem,numcp,i,JELEM,Numpatch,sctr,iV
      dimension sctr(MAXVAL(NNODE))
      Integer :: activeElementMap
      dimension activeElementMap(MAXVAL(nb_elem_patch))
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Compute gradient .................................................
      
      gradV(:,:,:) = zero
      massTot(:)   = zero
      
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
c     Loop on element
         Do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            
c     Get element infos
            CALL extractNurbsElementInfos(num_elem)
            sctr(:nnode_patch)  = IEN_patch(:,num_elem)
            
            U_elem(:,:,:) = zero
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               U_elem(:,:,i)    = VECT(:,:,sctr(i))
            Enddo
            
            gradV_elem(:,:,:) = zero
            massElem(:) = zero
            
c     Compute initial elementary matrix and load vector
            MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
            If (ELT_TYPE_patch == 'U1') then
               ! 'Element classique solide'
               If (activeElement(JELEM)==1) then
               call gradUELMASSMAT1(activeDir,
     1              U_elem(:,:,:nnode_patch),VALS(:),nb_frq,NDOFEL,
     1              MCRD,nnode_patch,JELEM,NBINT(NumPatch),COORDS_elem,
     2              TENSOR_patch,MAT_patch,RHO(NumPatch),PROPS_patch,
     4              JPROPS_patch,gradV_elem(:,:,:nnode_patch),
     5              massElem(:))
               Endif
               
            elseif (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               If (activeElement(JELEM)==1) then
                  print*,'Element type U3 not available.'
               Endif
               
            elseif (ELT_TYPE_patch == 'U30') then
               ! 'Element coque'
               If (activeElement(JELEM)==1) then
                  print*,'Element type U30 not available.'
               Endif
            Endif
            
            Do numcp = 1,nnode_patch
               gradV(:,:,sctr(numcp)) = gradV(:,:,sctr(numcp))
     &              + gradV_elem(:,:,numcp)
            Enddo
            massTot(:) = massTot(:) + massElem(:)
            
         Enddo ! end loop on element
         
         call deallocateMappingData()
         call finalizeNurbsPatch()

      Enddo ! end loop on patch

c     Scaling
      Do iV = 1,nb_frq
         gradV(iV,:,:) = gradV(iV,:,:)/massTot(iV)
      Enddo
      
c     Fin calcul .......................................................
      
      end subroutine gradVibration_AN
