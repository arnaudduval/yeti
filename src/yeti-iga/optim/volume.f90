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

!! Compute volume

      
subroutine computeVolume(VOL, listpatch,     &
            &  COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE, &
            &  PROPS,JPROPS,TENSOR,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp)
      
    use parameters
    use nurbspatch
    use embeddedMapping
      
    implicit None
      
      
    !! Variables declaration
      
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
    double precision, intent(in) :: PROPS
    integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,   &
        &     nb_elem_patch,IEN,JPROPS
    dimension PROPS(:), JPROPS(nb_patch), nb_elem_patch(nb_patch),  &
        &     IEN(:), NNODE(nb_patch), NBINT(nb_patch)
      
    !! Other infos
    integer, intent(in) :: listpatch
    dimension listpatch(nb_patch)
      
    !! Output variables : coefficients in diagonal of mass matrix
    !! ----------------
    double precision, intent(out) :: VOL
      
 
    !! Local Variables
    !! ---------------
            
    !! NURBS basis functions
    double precision :: COORDS_elem,XI,R,dRdxi,detJ,normV,AI,vectV
    dimension COORDS_elem(3,MAXVAL(NNODE)),XI(3),R(MAXVAL(NNODE)),  &
        &     dRdxi(MAXVAL(NNODE),3),AI(3,3),vectV(3)
    integer :: i,j,k,NumPatch,num_elem
      
    !! Embedded entities
    double precision :: COORDSmap,Re,dRedxi,BI
    dimension COORDSmap(3,MAXVAL(NNODE)),Re(MAXVAL(NNODE)), &
        &     dRedxi(MAXVAL(NNODE),3),BI(3,3)
    integer          :: sctr_map,isave
    dimension sctr_map(MAXVAL(NNODE))

    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(4,MAXVAL(NBINT))
      
      
    !! Compute Volume

    isave = 0
    VOL = zero
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
                COORDS_elem(:,:) = zero
                do i = 1,nnode_patch
                    COORDS_elem(:,i) = COORDS3D(:,IEN_patch(i,num_elem))
                enddo
                call extractNurbsElementInfos(num_elem)
            
                !! Loop on gauss points
                do n = 1,NBINT(numPatch)
               
                    !! Evaluate NURBS basis functions
                    XI(:) = zero
                    DetJ  = GaussPdsCoord(1,n)
                    do i = 1,dim_patch
                        XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))    &
                        &       *GaussPdsCoord(1+i,n)   &
                        &       +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0

                        DetJ = DetJ * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
                    enddo
               
                    call evalnurbs(XI(:),R(:nnode_patch),   &
                        &              dRdxi(:nnode_patch,:))
               
                    if ((ELT_TYPE_patch == 'U30').or.(ELT_TYPE_patch == 'U10')) then
                        !! For embedded cases
                        XI(:)   = zero
                        BI(:,:) = zero
                        do i = 1,nnode_patch
                            XI(:) = XI(:) + R(i)*COORDS_elem(:,i)
                            do j = 1, dim_patch
                                BI(:,j) = BI(:,j) + dRdxi(i,j)*COORDS_elem(:,i)
                            enddo
                        enddo
                  
                        dRedxi(:,:) = zero
                        call updateMapElementNumber(XI(:))
                        call evalnurbs_mapping(XI(:),Re(:nnode_map),    &
                            &                 dRedxi(:nnode_map,:))
                  
                        dRdxi(:,:) = zero
                        do i = 1,3
                            do j = 1, dim_patch
                                dRdxi(:,j) = dRdxi(:,j) + BI(i,j)*dRedxi(:,i)
                            enddo
                        enddo

                        !! extract COORDS
                        if (isave /= current_map_elem) then
                            sctr_map(:nnode_map) = IEN_map(:,current_map_elem)
                     
                            do i = 1,nnode_map
                                COORDSmap(:,i) = COORDS3D(:,sctr_map(i))
                            enddo
                     
                            isave = current_map_elem
                        endif
                  
                        AI(:,:) = zero
                        do k = 1,dim_patch
                            do i = 1,nnode_map
                                AI(:,k) = AI(:,k) + dRdxi(i,k)*COORDSmap(:,i)
                            enddo
                        enddo
                  
                    else
                        !! Classical patch
                        AI(:,:) = zero
                        do k = 1,dim_patch
                            do i = 1,nnode_patch
                                AI(:,k) = AI(:,k) + dRdxi(i,k)*COORDS_elem(:,i)
                            enddo
                        enddo

                    endif

                    !! Update volume
                    if (dim_patch == 1) then
                        ! curve
                        call norm(AI(:,1),3, normV)
                    elseif (dim_patch == 2) then
                        ! surface
                        call cross(AI(:,1),AI(:,2),vectV(:))
                        call norm(vectV(:),3, normV)
                    elseif (dim_patch == 3) then
                        ! volume
                        call cross(AI(:,1),AI(:,2),vectV(:))
                        call dot(  AI(:,3),vectV(:),normV)
                        call MatrixDet(AI, normV, 3)
                    endif
               
                    VOL = VOL + normV*detJ
                enddo
            enddo
        endif
         
        call finalizeNurbsPatch()

    enddo
      
end subroutine computeVolume
