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


subroutine gradLinElastWork_AN(     &
    &      gradWint,gradWext, computeWint,computeWext,SOL,ADJ,nadj,     &
    &      activeElement,activeDir,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,  &
    &     Nijk,weight,Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,  &
    &     RHO,TENSOR,indDLoad,JDLType,ADLMAG,load_target_nbelem,       &
    &     load_additionalInfos,MCRD,NBINT,nb_load,nb_patch,nb_elem,    &
    &     nnode,nb_cp)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------

    !! NURBS geometry
    integer, intent(in)          :: nb_cp
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3,nb_cp)

    double precision, intent(in) :: Ukv, weight
    integer, intent(in)          :: Nkv, Jpqr, Nijk
    dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),   &
        &     Ukv(:),weight(:)

    !! Patches and Elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: MATERIAL_PROPERTIES,RHO,PROPS
    integer, intent(in)          :: MCRD,NNODE,nb_patch,nb_elem,NBINT,  &
        &                           IEN,nb_elem_patch, JPROPS
    dimension MATERIAL_PROPERTIES(2,nb_patch),    &
        &     RHO(nb_patch),                &
        &     PROPS(:),                     &
        &     NNODE(nb_patch),              &
        &     IEN(:),                       &
        &     nb_elem_patch(nb_patch),      &
        &     JPROPS(nb_patch),             &
        &     NBINT(nb_patch)

    !! Loads
    double precision, intent(in) :: ADLMAG,load_additionalInfos
    integer, intent(in)          :: nb_load,indDLoad,JDLType,       &
        &                           load_target_nbelem
    dimension ADLMAG(nb_load),              &
        &     load_additionalInfos(:),      &
        &     indDLoad(:),                  &
        &     JDLType(nb_load),             &
        &     load_target_nbelem(nb_load)

    !! INFOS
    logical, intent(in)          :: computeWint,computeWext
    integer, intent(in)          :: activeElement, activeDir
    dimension activeElement(nb_elem), activeDir(3)

    !! Solution
    integer,          intent(in) :: nadj
    double precision, intent(in) :: SOL,ADJ
    dimension SOL(MCRD,nb_cp),ADJ(nadj,MCRD,nb_cp)

    !! Output variables
    !! ----------------
    double precision, intent(out):: gradWint,gradWext
    dimension gradWint(nadj,3,nb_cp),gradWext(nadj,3,nb_cp)

    !! Local variables
    !! ---------------

    !! gradUELMAT
    integer :: NDOFEL,dir
    double precision :: COORDS_elem,MAT_patch,      &
        &               gradWint_elem,gradWext_elem
    dimension COORDS_elem(3,MAXVAL(NNODE)),MAT_patch(2),    &
        &     gradWint_elem(nadj,3,MAXVAL(NNODE)),          &
        &     gradWext_elem(nadj,3,MAXVAL(NNODE))

    !! Solution
    double precision :: U_elem,UA_elem
    dimension U_elem(3,MAXVAL(NNODE)),UA_elem(3,MAXVAL(NNODE),nadj)

    !! Assembly
    integer :: num_elem,numcp,n,i,j,iA,kk,ll,JELEM,Numpatch,sctr
    dimension sctr(MAXVAL(NNODE))
    integer :: activeElementMap
    dimension activeElementMap(MAXVAL(nb_elem_patch))


    !! Compute gradient

    gradWint(:,:,:) = zero
    gradWext(:,:,:) = zero

    JELEM = 0
    do NumPatch = 1,nb_patch
        call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,     &
            &        weight,nb_elem_patch)
        call extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,      &
            &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

        if ((ELT_TYPE_patch == 'U30')           &
                &    .or.(ELT_TYPE_patch == 'U10')) then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,     &
                &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)
        endif


        NDOFEL = nnode_patch*MCRD

        !! Loop on element
        do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1

            !! Get element infos
            call extractNurbsElementInfos(num_elem)
            sctr(:nnode_patch)  = IEN_patch(:,num_elem)

            U_elem(:,:) = zero
            do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               U_elem(:MCRD,i)      = SOL(:,sctr(i))
            enddo

            UA_elem(:,:,:) = zero
            do i = 1,nnode_patch
                do iA = 1,nadj
                    UA_elem(:MCRD,i,iA) = ADJ(iA,:,sctr(i))
                enddo
            enddo

            gradWint_elem(:,:,:) = zero
            gradWext_elem(:,:,:) = zero

            !! Compute initial elementary matrix and load vector
            MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
            if (ELT_TYPE_patch == 'U1') then
                !! Solid element
                if (activeElement(JELEM)==1) then
                    call gradUELMAT1adj(activeDir,U_elem(:,:nnode_patch),                   &
                        &               UA_elem(:,:nnode_patch,:),NADJ,NDOFEL,MCRD,         &
                        &               nnode_patch,JELEM,NBINT(NumPatch),COORDS_elem,      &
                        &               TENSOR_patch,MAT_patch,RHO(NumPatch),PROPS_patch,   &
                        &               JPROPS_patch,nb_load,indDLoad,load_target_nbelem,   &
                        &               JDLType,ADLMAG,load_additionalInfos,                &
                        &               size(load_additionalInfos),                         &
                        &               computeWint,computeWext,                            &
                        &               gradWint_elem(:,:,:nnode_patch),                    &
                        &               gradWext_elem(:,:,:nnode_patch))

                    do numcp = 1,nnode_patch
                        if (computeWint) then
                            gradWint(:,:,sctr(numcp))=gradWint(:,:,sctr(numcp))     &
                                &       + gradWint_elem(:,:,numcp)
                        endif
                        if (computeWext) then
                            gradWext(:,:,sctr(numcp))=gradWext(:,:,sctr(numcp))     &
                                &       + gradWext_elem(:,:,numcp)
                        endif
                    enddo
                endif

            elseif (ELT_TYPE_patch == 'U3') then
                !! Shell element
                if (activeElement(JELEM)==1) then
                    call gradUELMAT3adj(activeDir,U_elem(:,:nnode_patch),                   &
                        &               UA_elem(:,:nnode_patch,:),NADJ,NDOFEL,MCRD,         &
                        &               nnode_patch,JELEM,NBINT(NumPatch),COORDS_elem,      &
                        &               TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,    &
                        &               nb_load,indDLoad,load_target_nbelem,JDLType,        &
                        &               ADLMAG,computeWint,computeWext,                     &
                        &               gradWint_elem(:,:,:nnode_patch),                    &
                        &                 gradWext_elem(:,:,:nnode_patch))

                    do numcp = 1,nnode_patch
                        if (computeWint) then
                            gradWint(:,:,sctr(numcp))=gradWint(:,:,sctr(numcp))     &
                                &               + gradWint_elem(:,:,numcp)
                        endif
                        if (computeWext) then
                            gradWext(:,:,sctr(numcp))=gradWext(:,:,sctr(numcp))     &
                                &               + gradWext_elem(:,:,numcp)
                        endif
                    enddo
                endif

            elseif(ELT_TYPE_patch == 'U10') then
                !! Embedded solid element
                if (activeElement(jelem) == 1) then
                    ! kk = int(PROPS_patch(2))
                    ! call gradUELMAT10adj(U_elem(:,:nnode_patch),                            &
                    !     &                UA_elem(:,:nnode_patch,:), nadj, mcrd,             &
                    !     &                nnode_patch,                                       &
                    !     &                nnode(kk), nb_cp, jelem, nbint(numpatch),          &
                    !     &                COORDS_elem,                                       &
                    !     &                COORDS3D, TENSOR_patch, MAT_patch,                 &
                    !     &                RHO(numpatch),                                     &
                    !     &                nb_load, indDload, load_target_nbelem,             &
                    !     &                JDLType,                                           &
                    !     &                ADLMAG, load_additionalInfos,                      &
                    !     &                size(load_additionalInfos),                        &
                    !     &                computeWint, computeWext,                          &
                    !     &                gradWint_elem(:,:,:nnode_patch),                   &
                    !     &                gradWext_elem(:,:,:nnode_patch))

                    ! do numcp = 1,nnode_patch
                    !     if (computeWint) then
                    !         gradWint(:,:,sctr(numcp)) =             &
                    !             &       gradWint(:,:,sctr(numcp))   &
                    !             &       + gradWint_elem(:,:,numcp)
                    !     endif
                    !     if (computeWext) then
                    !         gradWext(:,:,sctr(numcp)) =             &
                    !             &       gradWext(:,:,sctr(numcp))   &
                    !             &       + gradWext_elem(:,:,numcp)
                    !     endif
                    ! enddo

                    !! New version with global computation to take into account derivatives
                    !! w.r.t hulle control points
                    kk = int(PROPS_patch(2))
                    n = nb_elem_patch(kk)
                    j = 0
                    do ll = 1, kk
                        i = j+1
                        j = j + nb_elem_patch(ll)
                    enddo
                    activeElementMap(:n) = activeElement(i:j)
                    ! if (sum(activeElementMap(:n)) > 0) then
                        call gradUELMAT10adj(activeElementMap, n,                               &
                            &                sctr(:nnode_patch),  &
                            &                U_elem(:,:nnode_patch),                            &
                            &                UA_elem(:,:nnode_patch,:), nadj, mcrd,             &
                            &                nnode_patch,                                       &
                            &                nnode(kk), nb_cp, jelem, nbint(numpatch),          &
                            &                COORDS_elem,                                       &
                            &                COORDS3D, TENSOR_patch, MAT_patch,                 &
                            &                RHO(numpatch),                                     &
                            &                nb_load, indDload, load_target_nbelem,             &
                            &                JDLType,                                           &
                            &                ADLMAG, load_additionalInfos,                      &
                            &                size(load_additionalInfos),                        &
                            &                computeWint, computeWext,                          &
                            &                gradWint(:,:,:),                   &
                            &                gradWext(:,:,:))

                    ! endif
                endif

            elseif (ELT_TYPE_patch == 'U30') then
                !! Embedded shell element
                !! Disable warning
                ! if (activeElement(JELEM)==1) then
                !     print*,'Elt type U30 not available yet. (gradtotalwork_fullAN)'
                ! endif

                kk = int(PROPS_patch(2))
                n  = nb_elem_patch(kk)
                j = 0
                do ll = 1,kk
                    i = j+1
                    j = j+nb_elem_patch(ll)
                enddo
                activeElementMap(:n) = activeElement(i:j)
                if (sum(activeElementMap(:n))>0) then
                    call gradUELMAT30gloAdj(activeDir,activeElementMap,n,                       &
                        &                   U_elem(:,:nnode_patch),UA_elem(:,:nnode_patch,:),   &
                        &                   NADJ,NDOFEL,MCRD,nnode_patch,NNODE(kk),nb_cp,       &
                        &                   JELEM,NBINT(NumPatch),COORDS_elem,COORDS3D,         &
                        &                   TENSOR_patch,MAT_patch,PROPS_patch,JPROPS_patch,    &
                        &                   nb_load,indDLoad,load_target_nbelem,JDLType,        &
                        &                   ADLMAG,computeWint,computeWext,                     &
                        &                   gradWint(:,:,:),gradWext(:,:,:))
                endif
            endif

        enddo !! end loop on element

        call deallocateMappingData()
        call finalizeNurbsPatch()

    enddo ! end loop on patch

end subroutine gradLinElastWork_AN
