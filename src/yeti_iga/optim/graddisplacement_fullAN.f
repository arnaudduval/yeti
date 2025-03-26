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


      subroutine gradDisp_AN(gradD, SOL,ADJ,nadj,
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
      Integer,          intent(in) :: nadj
      Double precision, intent(in) :: SOL,ADJ
      dimension SOL(MCRD,nb_cp),ADJ(nadj,MCRD,nb_cp)

c     Output variables :
c     ----------------
      Double precision, intent(out):: gradD
      dimension gradD(nadj,3,nb_cp)


c     Local variables :
c     ---------------

!     gradUELMAT.f
      Integer :: NDOFEL,dir
      Double precision :: COORDS_elem,gradD_elem,MAT_patch
      dimension COORDS_elem(3,MAXVAL(NNODE)),MAT_patch(2),
     &     gradD_elem(nadj,3,MAXVAL(NNODE))

!     Solution
      Double precision :: U_elem,UA_elem
      dimension U_elem(3,MAXVAL(NNODE)),UA_elem(3,MAXVAL(NNODE),nadj)

!     Assembly
      Integer :: num_elem,numcp,n,i,j,iA,kk,ll,JELEM,Numpatch,sctr
      dimension sctr(MAXVAL(NNODE))
      Integer :: activeElementMap
      dimension activeElementMap(MAXVAL(nb_elem_patch))

C     Fin declaration des variables ....................................
c
c
c
c
c     Compute gradient .................................................

      gradD(:,:,:) = zero

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

            U_elem(:,:) = zero
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               U_elem(:MCRD,i)      = SOL(:,sctr(i))
            Enddo

            UA_elem(:,:,:) = zero
            Do i = 1,nnode_patch
            Do iA = 1,nadj
               UA_elem(:MCRD,i,iA) = ADJ(iA,:,sctr(i))
            Enddo
            Enddo

            gradD_elem(:,:,:) = zero

c     Compute initial elementary matrix and load vector
            MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
            If (ELT_TYPE_patch == 'U1') then
               ! 'Element classique solide'
               If (activeElement(JELEM)==1) then
               call gradUELMAT1adj(activeDir,U_elem(:,:nnode_patch),
     &                 UA_elem(:,:nnode_patch,:),NADJ,NDOFEL,MCRD,
     &                 nnode_patch,JELEM,NBINT(NumPatch),COORDS_elem,
     &                 TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,
     &                 nb_load,indDLoad,load_target_nbelem,JDLType,
     &                 ADLMAG,gradD_elem(:,:,:nnode_patch))

               Do numcp = 1,nnode_patch
                  gradD(:,:,sctr(numcp)) = gradD(:,:,sctr(numcp))
     &                 + gradD_elem(:,:,numcp)
               Enddo
               Endif

            elseif (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               If (activeElement(JELEM)==1) then
               call gradUELMAT3adj(activeDir,U_elem(:,:nnode_patch),
     &                 UA_elem(:,:nnode_patch,:),NADJ,NDOFEL,MCRD,
     &                 nnode_patch,JELEM,NBINT(NumPatch),COORDS_elem,
     &                 TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,
     &                 nb_load,indDLoad,load_target_nbelem,JDLType,
     &                 ADLMAG,gradD_elem(:,:,:nnode_patch))

               Do numcp = 1,nnode_patch
                  gradD(:,:,sctr(numcp)) = gradD(:,:,sctr(numcp))
     &                 + gradD_elem(:,:,numcp)
               Enddo
               Endif

            elseif (ELT_TYPE_patch == 'U30') then
               ! 'Element coque immerge'
               If (activeElement(JELEM)==1) then
               print*,'Elt type U30 not available yet. (graddisplacement_fullAN)'
               Endif


               kk = int(PROPS_patch(2))
               n  = nb_elem_patch(kk)
               j = 0
               Do ll = 1,kk
                  i = j+1
                  j = j+nb_elem_patch(ll)
               Enddo
               activeElementMap(:n) = activeElement(i:j)

               IF (SUM(activeElementMap(:n))>0) then

               Endif

            Endif

         Enddo ! end loop on element

         call deallocateMappingData()
         call finalizeNurbsPatch()

      Enddo ! end loop on patch

c     Fin calcul .......................................................

      end subroutine gradDisp_AN
