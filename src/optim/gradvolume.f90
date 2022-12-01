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

subroutine computeGradVolume(gradV, listpatch,      &
            & COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE, &
            & PROPS,JPROPS,TENSOR,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp)
      
    use parameters
    use nurbspatch
    use embeddedMapping
      
    implicit None
     
    !! Input arguments
    !! ---------------
      
    !! NURBS geometry
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3,nb_cp)

    double precision, intent(in) :: Ukv, weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),   &
     &     Ukv(:),weight(:)
      
    !! Patches and Elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: PROPS
    integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,       &
     &     nb_elem_patch,IEN,JPROPS
    dimension PROPS(:), JPROPS(nb_patch), nb_elem_patch(nb_patch),  &
     &     IEN(:), NNODE(nb_patch), NBINT(nb_patch)
      
    !! Other infos
    integer, intent(in) :: listpatch
    dimension listpatch(nb_patch)
      
      
    !! Output variables
    !! ----------------
    double precision, intent(out) :: gradV
    dimension gradV(3,nb_cp)
      
    !! Local Variables :
    !! ---------------
            
    !! NURBS basis functions
    double precision :: COORDS_elem,XI,R,dRdxi,detJ
    dimension COORDS_elem(3,MAXVAL(NNODE)),XI(3),R(MAXVAL(NNODE)),  &
     &     dRdxi(MAXVAL(NNODE),3)
    integer :: i,j,k,l,kk,NumPatch,num_elem,sctr,idim
    dimension sctr(MAXVAL(NNODE))

    !! Embedded entities
    double precision :: COORDSmap,Re,dRedxi,Rm,dRmdxi,ddRmddxi,BI,VI,   &
     &     dVI
    dimension COORDSmap(3,MAXVAL(NNODE)),BI(3,3),VI(3,3),   &
     &     Re(MAXVAL(NNODE)),dRedxi(MAXVAL(NNODE),3),   &
     &     Rm(MAXVAL(NNODE)),dRmdxi(MAXVAL(NNODE),3),   &
     &     ddRmddxi(MAXVAL(NNODE),6),dVI(3,6)
    integer          :: sctr_map,isave
    dimension sctr_map(MAXVAL(NNODE))
      
    !! grad
    double precision :: AIxAJ,AI,normV!, dA1dPe,dA2dPe,dA1dPm,dA2dPm,    &
    !  &     coef1,coef2
    dimension AIxAJ(3,3),AI(3,3)!!,     &
    !  &     dA1dPe(3,3,MAXVAL(NNODE)),dA2dPe(3,3,MAXVAL(NNODE)), &
    !  &     dA1dPm(3,3,MAXVAL(NNODE)),dA2dPm(3,3,MAXVAL(NNODE))
    double precision :: dAidPe, dAidPm
    double precision coefi
    dimension coefi(3)
    dimension dAidPe(3,3,3,MAXVAL(NNODE)), dAidPm(3,3,3,MAXVAL(NNODE))
      
    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(4,MAXVAL(NBINT))
      
    !! Compute Volume ....................................................
      
    isave = 0
    gradV(:,:) = zero
    !! Loop on patches
    do NumPatch = 1,nb_patch
        if (listpatch(numPatch) == 1) then
         
            call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv, &
                &        weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,  &
                &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
            
            if ((ELT_TYPE_patch == 'U30').or.(ELT_TYPE_patch == 'U10')) then
                i = int(PROPS_patch(2))
                call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
                    &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
            endif

            !! Get Gauss infos
            NbPtInt = int( NBINT(numPatch)**(1.0/float(dim_patch)) )
            if (NbPtInt**dim_patch<NBINT(numPatch)) NbPtInt = NbPtInt+1
            call Gauss(NbPtInt,dim_patch,   &
                &        GaussPdsCoord(:dim_patch+1,:NBINT(numPatch)),0)
            
            !! Loop on elements
            do num_elem = 1,nb_elem_patch(NumPatch)
                
                sctr(:nnode_patch) = IEN_patch(:,num_elem)
                do i = 1,nnode_patch
                COORDS_elem(:,i) = COORDS3D(:,sctr(i))
                enddo
                call extractNurbsElementInfos(num_elem)
                
                !! Loop on gauss points
                do n = 1,NBINT(numPatch)
                
                    !! Evaluate nurbs basis functions
                    XI(:) = zero
                    DetJ  = GaussPdsCoord(1,n)
                    do i = 1,dim_patch
                        XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))    &
                            & *GaussPdsCoord(1+i,n)     &
                            & +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
                        DetJ = DetJ * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
                    enddo
                
                    R(:) = zero
                    dRdxi(:,:) = zero
                    call evalnurbs(XI(:),R(:nnode_patch),       &
                        &   dRdxi(:nnode_patch,:))
                
                    !! Build tangent vectors

                    !! For embedded cases
                    if ((ELT_TYPE_patch == 'U30').or.(ELT_TYPE_patch == 'U10')) then
                        Re(:) = R(:)
                        dRedxi(:,:) = dRdxi(:,:)
                    
                        XI(:)   = zero
                        BI(:,:) = zero
                        do i = 1,nnode_patch
                            XI(:) = XI(:) + Re(i)*COORDS_elem(:,i)
                            do idim = 1, dim_patch
                                BI(:,idim) = BI(:,idim) + dRedxi(i,idim)*COORDS_elem(:,i)
                            enddo
                            !!BI(:,1) = BI(:,1) + dRedxi(i,1)*COORDS_elem(:,i)
                            !!BI(:,2) = BI(:,2) + dRedxi(i,2)*COORDS_elem(:,i)
                        enddo
                    
                        call updateMapElementNumber(XI(:))
                        call evalnurbs_mapping_w2ndDerv(XI(:),Rm(:nnode_map),   &
                            &   dRmdxi(:nnode_map,:),ddRmddxi(:nnode_map,:))
                    
                        ! extract coords mapping
                        if (isave /= current_map_elem) then
                            sctr_map(:nnode_map) = IEN_map(:,current_map_elem)
                        
                            do i = 1,nnode_map
                                COORDSmap(:,i) = COORDS3D(:,sctr_map(i))
                            enddo
                        
                            isave = current_map_elem
                        endif
                    
                        VI(:,:) = zero
                        do k = 1,dim_map
                            do i = 1,nnode_map
                                VI(:,k) = VI(:,k) + dRmdxi(i,k)*COORDSmap(:,i)
                            enddo
                        enddo
                    
                        dVi(:,:) = zero
                        do k = 1,6
                            do i = 1,nnode_map
                                dVI(:,k)=dVI(:,k) + ddRmddxi(i,k)*COORDSmap(:,i)
                            enddo
                        enddo
                    
                        AI(:,:) = zero
                        do k = 1,dim_patch
                            do i = 1,3
                                AI(:,k) = AI(:,k) + BI(i,k)*VI(:,i)
                            enddo
                        enddo

                        !! compute derivative versus the control points
                        !! - embedded CPs
                        ! dA1dPe(:,:,:) = zero
                        ! dA2dPe(:,:,:) = zero
                        dAidPe(:,:,:,:) = zero
                        do i = 1,nnode_patch
                            do l = 1,3
                                ! dA1dPe(:,l,i) = dRedxi(i,1)*VI(:,l)     &
                                !     &   + BI(l,1)*R(i)*dVI(:,l)
                                ! dA2dPe(:,l,i) = dRedxi(i,2)*VI(:,l)     &
                                !     &   + BI(l,2)*R(i)*dVI(:,l)
                                ! do k = 1,3
                                !     if (k /= l) then
                                !         kk = l+k+1
                                !         dA1dPe(:,l,i) = dA1dPe(:,l,i)       &
                                !             &   + BI(k,1)*R(i)*dVI(:,kk)
                                !         dA2dPe(:,l,i) = dA2dPe(:,l,i)       &
                                !             &   + BI(k,2)*R(i)*dVI(:,kk)
                                !     endif
                                ! enddo
                                do idim  = 1, dim_patch
                                    dAidPe(idim,:,l,i) = dRedxi(i,idim)*VI(:,l)     &
                                     &   + BI(l,idim)*R(i)*dVI(:,l)
                                    do k = 1, 3
                                        if (k /= l) then
                                            kk = l+k+1
                                            dAidPe(idim,:,l,i) = dAidPe(idim,:,l,i)       &
                                                &   + BI(k,idim)*R(i)*dVI(:,kk)
                                        endif
                                    enddo
                                enddo
                            enddo
                        enddo

                        !! - mapping CPs
                        ! dA1dPm(:,:,:) = zero
                        ! dA2dPm(:,:,:) = zero
                        dAidPm(:,:,:,:) = zero

                        do i = 1,nnode_map
                            ! coef1 = SUM( BI(:,1)*dRmdxi(i,:) )
                            ! coef2 = SUM( BI(:,2)*dRmdxi(i,:) )
                            do idim = 1, dim_patch
                                coefi(idim) = SUM( BI(:,idim)*dRmdxi(i,:) )
                            enddo
                            do k = 1,3
                                do idim = 1, dim_patch
                                    dAidPm(idim,k,k,i) = coefi(idim)
                                enddo
                                ! dA1dPm(k,k,i) = coef1
                                ! dA2dPm(k,k,i) = coef2
                            enddo
                        enddo

                    else
                        AI(:,:) = zero
                        do k = 1,dim_patch
                            do i = 1,nnode_patch
                                AI(:,k) = AI(:,k) + dRdxi(i,k)*COORDS_elem(:,i)
                            enddo
                        enddo
                    endif

                    !! Update volume
                    if     (dim_patch == 1) then
                        ! curve
                        ! --> to do
                    elseif (dim_patch == 2) then
                        ! surface
                        call cross(AI(:,1),AI(:,2), AI(:,3))
                        call norm(AI(:,3),3, normV)
                        AI(:,3) = AI(:,3)/normV
                        call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
                        call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
                    
                        if (ELT_TYPE_patch == 'U30') then
                            ! embedded surface
                            do i = 1,nnode_patch
                                k = sctr(i)
                                do l = 1,3
                                    gradV(l,k) = gradV(l,k) +   &
                                        &   ( SUM( AIxAJ(:,1) * dAidPe(1,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,2) * dAidPe(2,:,l,i) )     &
                                        &   ) * detJ
                                enddo
                            enddo

                            ! mapping
                            do i = 1,nnode_map
                                k = sctr_map(i)
                                do l = 1,3
                                    gradV(l,k) = gradV(l,k) +       &
                                        &   ( SUM( AIxAJ(:,1) * dAidPm(1,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,2) * dAidPm(2,:,l,i) )     &
                                        &   ) * detJ
                                enddo
                            enddo
                        else
                            do i = 1,nnode_patch
                                k = sctr(i)
                                gradV(:,k) = gradV(:,k) +       &
                                    &   ( AIxAJ(:,1)*dRdxi(i,1)     &
                                    &   + AIxAJ(:,2)*dRdxi(i,2)     &
                                    &   ) * detJ
                            enddo
                        endif


                    elseif (dim_patch == 3) then
                        ! volume
                        call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
                        call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
                        call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
                        if (ELT_TYPE_patch == 'U10') then
                            ! embedded volume
                            do i = 1,nnode_patch
                                k = sctr(i)
                                do l = 1,3
                                    gradV(l,k) = gradV(l,k) +   &
                                        &   ( SUM( AIxAJ(:,1) * dAidPe(1,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,2) * dAidPe(2,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,3) * dAidPe(3,:,l,i) )     &
                                        &   ) * detJ
                                enddo
                            enddo

                            ! mapping
                            do i = 1,nnode_map
                                k = sctr_map(i)
                                do l = 1,3
                                    gradV(l,k) = gradV(l,k) +       &
                                        &   ( SUM( AIxAJ(:,1) * dAidPm(1,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,2) * dAidPm(2,:,l,i) )     &
                                        &   + SUM( AIxAJ(:,3) * dAidPm(3,:,l,i) )     &
                                        &   ) * detJ
                                enddo
                            enddo
                        else
                            do i = 1,nnode_patch
                                k = sctr(i)
                                gradV(:,k) = gradV(:,k) +       &
                                    &   ( AIxAJ(:,3)*dRdxi(i,3)     &
                                    &   + AIxAJ(:,1)*dRdxi(i,1)     &
                                    &   + AIxAJ(:,2)*dRdxi(i,2)     &
                                    &   ) * detJ
                            enddo
                        endif
                    endif
                enddo
            enddo
        endif
    
        call finalizeNurbsPatch()
         
    enddo
      
end subroutine computeGradVolume
