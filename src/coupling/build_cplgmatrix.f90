!! Copyright 2016-2019 Thibaut Hirschler
!! Copyright 2020-2021 Arnaud Duval
!! Copyright 2021 Marie Guerder

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

!! Build sparse coupling matrix
!! Returns
!!      Cdata : vector containing values of elementary matrices
!!      Crow : vector containing row indices of non zero values
!!      Ccol : vector containing column indices of non zero values

subroutine cplg_matrix(Cdata,Crow,Ccol,nb_data, &
                &   COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,   &
                &   PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,ind_dof_free,   &
                &   nb_dof_free,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp,    &
                &   nb_dof_tot)
      
    use parameters
    use nurbspatch
    use embeddedMapping
      
    implicit none

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
    double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
    integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,IEN,   &
        &   nb_elem_patch, JPROPS
    dimension MATERIAL_PROPERTIES(2,nb_patch),  &
        &   PROPS(:),   &
        &   NNODE(nb_patch),    &
        &   IEN(:), &
        &   nb_elem_patch(nb_patch),    &
        &   JPROPS(nb_patch),   &
        &   NBINT(nb_patch)
      
      
    !! Degree Of Freedom
    integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
    dimension ind_dof_free(nb_dof_tot)
      
      
    !! Storage INFOS
    integer, intent(in) :: nb_data
      
      
    !! Output variables
    !! ----------------
    integer,          intent(out) :: Crow,Ccol
    double precision, intent(out) :: Cdata
    dimension Cdata(nb_data),Crow(nb_data),Ccol(nb_data)
      
      
    !! Local variables
    !! ---------------
      
    !! Coupling infos
    integer :: numDomain,numLgrge,ismaster,dimInterface,nb_gps, &
        &   dispORrot
      
    !! Allocatable quantities
    double precision, dimension(:,:),   allocatable :: saveXI,saveXIb,  &
        &   saveRl
    double precision, dimension(:,:,:), allocatable :: saveBI
    integer         , dimension(:),     allocatable :: saveEL
    integer         , dimension(:,:),   allocatable :: saveIEN
      
    !! Extract infos
    integer, allocatable          :: sctr(:),sctr_l(:)
    
    double precision, allocatable :: COORDS_elem(:,:)
    
      
    !! CPLG matrix assembly
    double precision, allocatable :: CMAT(:,:,:)
    
      
    integer :: i,n,dofi,dofj,cpi,cpj,kk,k,l,count,numPatch,num_elem
      
    !! Allocations
    allocate(sctr(MAXVAL(NNODE)),sctr_l(MAXVAL(NNODE)))
    allocate(CMAT(MCRD,MCRD,MAXVAL(NNODE)**2))
    allocate(COORDS_elem(MCRD,MAXVAL(NNODE)))

    !! Initialisation
    Cdata = zero
    Crow  = 0
    Ccol  = 0
  
    !! Start assembly..........
      
    count = 1
    do NumPatch = 1,nb_patch
         
        call extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,  &
                &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
        if (ELT_TYPE_patch .eq. 'U00') then
         
            CMAT(:,:,:) = zero
         
            call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
         
            numDomain = int(PROPS_patch(2))
            numLgrge  = int(PROPS_patch(3))
            ismaster  = int(PROPS_patch(4))
            dimInterface = dim_patch
         
            !! 1. Integration points through the immerged patch
            !!  - save gps parametric coords to further evaluate the lagrangian field
            !!  - compute gps position on parameter patch
            !!  - compute covariant vectors
         
            n = NBINT(numPatch)
            nb_gps = n*nb_elem_patch(NumPatch)
            if (allocated(saveXIb)) deallocate(saveXIb)
            if (allocated(saveXI )) deallocate(saveXI )
            if (allocated(saveBI )) deallocate(saveBI )
            allocate(saveXIb(4,          nb_gps))
            allocate(saveXI( 3,          nb_gps))
            allocate(saveBI( 3,dim_patch,nb_gps))
         
            kk = 0
            do num_elem = 1,nb_elem_patch(NumPatch)
            
                do i = 1,nnode_patch
                    COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
                enddo
                call extractNurbsElementInfos(num_elem)
            
                call getGPsOnParamSpace(    &
                        &   saveXIb( :,kk+1:kk+n),saveXI(:,kk+1:kk+n),  &
                        &   saveBI(:,:,kk+1:kk+n),dim_patch,MCRD,nnode_patch,   &
                        &   NBINT(numPatch),COORDS_elem(:,:nnode_patch)      )
                
                kk = kk + n
                
            enddo
            
            
            !! 2. Lagrangian field
            !!  - compute basis functions
            call extractNurbsPatchGeoInfos(numLgrge, Nkv,Jpqr,Nijk,Ukv, &
                    &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS,  &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
            dispORrot = int(PROPS_patch(2))

            if (allocated(saveRl))  deallocate(saveRl)
            if (allocated(saveEL))  deallocate(saveEL)
            if (allocated(saveIEN)) deallocate(saveIEN)
            allocate(saveRl(nnode_patch,nb_gps))
            allocate(saveEL(nb_gps))
            allocate(saveIEN(nnode_patch,nb_elem_patch(numLgrge)))
         
            saveIEN(:,:) = IEN_patch(:,:)
         
            do i = 1,nb_gps

                call updateElementNumber(saveXIb(:3,i))
                call evalLgrge(saveXIb(:3,i),saveRl(:,i))
                
                saveEL(i) = current_elem
                
            enddo

         
            !! 3. Domain to couple
            !!  - compute basis functions
            !!  - compute Jacobian for the integral
            !!  - build coupling matrix
            call extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv,    &
                    &   weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS, &
                    &   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
            if ((ELT_TYPE_patch .eq. 'U30') .or. (ELT_TYPE_patch .eq. 'U10')) then
                i = int(PROPS_patch(2))
                call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
                    &   weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
            endif
         
         
            do kk = 1,nb_gps
                
                CMAT(:,:,:) = zero
                
                call updateElementNumber(saveXI(:,kk))
                sctr(:nnode_patch) = IEN_patch(:,current_elem)
                do i = 1,nnode_patch
                    COORDS_elem(:,i) = COORDS3D(:MCRD,sctr(i))
                enddo
            
                n = nnode_patch*NNODE(numLgrge)
                if ((ELT_TYPE_patch .eq. 'U30') .or. (ELT_TYPE_patch .eq. 'U10')) then
                    if  (dispORrot .eq. 0) then
                        call CPLINGDISP_embedded(saveRl(:,kk),saveXIb(4,kk),    &
                            &   saveXI(:,kk),saveBI(:,:,kk),dim_patch,dimInterface, &
                            &   MCRD,nnode_patch,NNODE(numLgrge),nnode_map,nb_cp,   &
                            &   COORDS_elem(:,:nnode_patch),COORDS3D(:,:),  &
                            &   CMAT(:,:,:n))
                    else
                        call CPLINGROT_embedded( saveRl(:,kk),saveXIb(4,kk),    &
                            &   saveXI(:,kk),saveBI(:,:,kk),dim_patch,dimInterface, &
                            &   MCRD,nnode_patch,NNODE(numLgrge),nnode_map,nb_cp,   &
                            &   COORDS_elem(:,:nnode_patch),COORDS3D(:,:),  &
                            &   CMAT(:,:,:n))
                    endif
                else
                    if (dispORrot .eq. 0) then
                        call CPLINGDISP(saveRl(:,kk),saveXIb(4,kk),saveXI(:,kk),    &
                            &   saveBI(:,:,kk),dim_patch,dimInterface,MCRD, &
                            &   nnode_patch,NNODE(numLgrge),    &
                            &   COORDS_elem(:,:nnode_patch),CMAT(:,:,:n))
                    elseif (dispORrot .eq. 1) then
                        call CPLINGROT( saveRl(:,kk),saveXIb(4,kk),saveXI(:,kk),    &
                            &   saveBI(:,:,kk),dim_patch,dimInterface,MCRD, &
                            &   nnode_patch,NNODE(numLgrge),    &
                            &   COORDS_elem(:,:nnode_patch),CMAT(:,:,:n))
                    !!elseif (dispORrot .eq. 2) then
                    !!    call CPLINGDISPderv( saveRl(:,kk),saveXIb(4,kk),    &
                    !!        &   saveXI(:,kk),saveBI(:,:,kk),dim_patch,dimInterface, &
                    !!        &   MCRD,nnode_patch,NNODE(numLgrge),   &
                    !!        &   COORDS_elem(:,:nnode_patch),CMAT(:,:,:n))
                    endif
                endif
                
                if (ismaster .eq. 0) CMAT(:,:,:n) = -CMAT(:,:,:n)
                
                !! ASSEMBLING
                sctr(:nnode_patch)       = IEN_patch(:,current_elem)
                sctr_l(:NNODE(numLgrge)) = saveIEN(:,saveEL(kk))
            
                i = 0
                do cpj = 1,NNODE(numLgrge)
                    dofj = (sctr_l(cpj)-1)*MCRD
               
                    do cpi = 1,nnode_patch
                        dofi= (sctr(cpi)-1)*MCRD
                        i   = i + 1
                        do l = 1,MCRD
                            do k = 1,MCRD
                                Cdata(count) = CMAT(k,l,i)
                                Crow( count) = dofi + k - 1
                                Ccol( count) = dofj + l - 1
                                count = count + 1
                            enddo
                        enddo
                    enddo
                enddo
            
            enddo
        endif
         
        call deallocateMappingData()
        call finalizeNurbsPatch()
         
        if (allocated(saveXI ))  deallocate(saveXI )
        if (allocated(saveXIb))  deallocate(saveXIb)
        if (allocated(saveBI ))  deallocate(saveBI )
        if (allocated(saveRl ))  deallocate(saveRl )
        if (allocated(saveEL ))  deallocate(saveEL )
        if (allocated(saveIEN))  deallocate(saveIEN)
         
    enddo ! end loop on patch

    !! Deallocations
    deallocate(sctr,sctr_l)
    deallocate(CMAT)
    deallocate(COORDS_elem)
      
end subroutine cplg_matrix
