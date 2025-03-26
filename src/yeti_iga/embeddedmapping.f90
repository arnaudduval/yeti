!! Copyright 2016-2020 Thibaut Hirschler
!! Copyright 2020 Arnaud Duval

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

MODULE embeddedMapping
  
  use parameters
  use nurbspatch

  Implicit none
  
  ! Patch geo infos
  Integer :: current_map, nnode_map, dim_map, nbel_map
  Integer, dimension(3) :: Jpqr_map, Nkv_map
  Integer,          allocatable, dimension(:,:) :: Nijk_map
  Double precision, allocatable, dimension(:)   :: Ukv1_map,Ukv2_map,  &
       Ukv3_map
  Double precision, allocatable, dimension(:,:) :: weight_map
  
  ! Patch mech infos
  Integer :: JPROPS_map
  Integer,          allocatable, dimension(:,:) :: IEN_map
  Double precision, allocatable, dimension(:)   :: PROPS_map
  Character(:), allocatable  :: ELT_TYPE_map,TENSOR_map
  
  ! Element infos
  Integer :: current_map_elem
  Double precision, allocatable, dimension(:) :: weight_map_elem
  Double precision, dimension(2,3) :: Ukv_map_elem
  

CONTAINS
  
  Subroutine extractMappingInfos(num_patch4map,nb_elem_patch,Nkv,Jpqr, &
       Nijk,Ukv,weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)
    implicit none
    
    ! Declaration des variables ........................................
    
    ! Input variables :
    ! ---------------
    Integer, intent(in) :: num_patch4map

    ! - geo infos
    Integer,          dimension(:),   intent(in) :: nb_elem_patch
    Integer,          dimension(:,:), intent(in) :: Nkv, Jpqr, Nijk
    Double precision, dimension(:),   intent(in) :: Ukv, weight
    ! - mech infos
    Integer,          dimension(:),   intent(in) :: IEN,JPROPS,NNODE
    Double precision, dimension(:),   intent(in) :: PROPS
    Character(len=*),                 intent(in) :: ELT_TYPE,TENSOR
    

    ! Local variables :
    ! ---------------
    Integer :: i_save,n
    
    ! Fin declaration des variables ....................................
    !
    !
    ! Initialisation ...................................................
    
    i_save = current_patch
    CALL extractNurbsPatchGeoInfos(num_patch4map, Nkv,Jpqr,Nijk,Ukv,   &
         weight,nb_elem_patch)
    CALL extractNurbsPatchMechInfos(num_patch4map,IEN,PROPS,JPROPS,    &
         NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
    CALL deallocateMappingData()
    
    ! copy data from nurbspatch module
    current_map = num_patch4map
    dim_map     = dim_patch
    nnode_map   = nnode_patch
    nbel_map    = nbel_patch
        
    ! - geo
    Jpqr_map(:) = Jpqr_patch(:)
    Nkv_map(:)  = Nkv_patch(:)
    
    allocate(Ukv1_map(Nkv_map(1)))
    Ukv1_map(:) = Ukv1_patch(:)
    if (dim_map>1) then
       allocate(Ukv2_map(Nkv_map(2)))
       Ukv2_map(:) = Ukv2_patch(:)
    endif
    if (dim_map>2) then
       allocate(Ukv3_map(Nkv_map(3)))
       Ukv3_map(:) = Ukv3_patch(:)
    endif

    n = nb_elem_patch(num_patch4map)
    allocate( Nijk_map(3,n) )
    Nijk_map(:,:) = Nijk_patch(:,:)
    
    allocate(weight_map(nnode_map,n))
    weight_map(:,:) = weight_patch(:,:)
    allocate(weight_map_elem(nnode_map))
    
    
    ! - mech
    JPROPS_map = JPROPS_patch
    allocate(PROPS_map(JPROPS_map))
    PROPS_map(:) = PROPS_patch(:)

    allocate(IEN_map(nnode_map,n))
    IEN_map(:,:) = IEN_patch(:,:)
    
    ELT_TYPE_map = ELT_TYPE_patch
    TENSOR_map   = TENSOR_patch
        
    CALL extractNurbsPatchGeoInfos(i_save, Nkv,Jpqr,Nijk,Ukv,weight,   &
         nb_elem_patch)
    CALL extractNurbsPatchMechInfos(i_save,IEN,PROPS,JPROPS,NNODE,     &
         nb_elem_patch,ELT_TYPE,TENSOR)
    
  End Subroutine extractMappingInfos
  
  
    
  Subroutine deallocateMappingData()
    implicit none
    ! geo
    if (allocated(Ukv1_map))   deallocate(Ukv1_map)
    if (allocated(Ukv2_map))   deallocate(Ukv2_map)
    if (allocated(Ukv3_map))   deallocate(Ukv3_map)
    if (allocated(Nijk_map))   deallocate(Nijk_map)
    if (allocated(weight_map)) deallocate(weight_map)
    if (allocated(weight_map_elem)) deallocate(weight_map_elem)
    ! mech
    if (allocated(PROPS_map))  deallocate(PROPS_map)
    if (allocated(IEN_map))    deallocate(IEN_map)
    if (allocated(ELT_TYPE_map))    deallocate(ELT_TYPE_map)
    if (allocated(TENSOR_map)) deallocate(TENSOR_map)
    
  End Subroutine deallocateMappingData


  Subroutine extractMapElementInfos(num_elem)
    implicit none
    integer, intent(in)   :: num_elem
    integer, dimension(3) :: Ni
    integer :: i
    
    current_map_elem   = num_elem
    weight_map_elem(:) = weight_map(:,num_elem)
    
    Ni(:) = Nijk_map(:,num_elem)
    Ukv_map_elem(:,:) = zero
    Ukv_map_elem(:,1) = Ukv1_map(Ni(1):Ni(1)+1)
    if (dim_map>1) Ukv_map_elem(:,2) = Ukv2_map(Ni(2):Ni(2)+1)
    if (dim_map>2) Ukv_map_elem(:,3) = Ukv3_map(Ni(3):Ni(3)+1)
    
  End Subroutine extractMapElementInfos

  

  Subroutine updateMapElementNumber(Xi)
    implicit none
    Double precision, dimension(3), intent(in) :: Xi
    
    Integer :: i,j,k,p,n,nb_elemD,numel,nbfound
    Logical :: mask
    dimension nb_elemD(3),mask(nbel_map,3), numel(1)
    
    ! ...
    !Nijk_map(:,current_map_elem)
    mask(:,:) =.True.
    
    nb_elemD(:) = 1
    
    p = Jpqr_map(1)
    n =  Nkv_map(1)
    call BINARYSEARCH(Xi(1), ukv1_map(p+1:n-(p+1)), n-2*p, i)
    i = i + p
    !nb_elemD(1) = n-2*p - 1
    
    mask(:,1) = Nijk_map(1,:) == i
    
    j = 0 !1
    if (dim_map >1) then
       p = Jpqr_map(2)
       n =  Nkv_map(2)
       call BINARYSEARCH(Xi(2), ukv2_map(p+1:n-(p+1)), n-2*p, j)
       j = j + p
       
       mask(:,2) = Nijk_map(2,:) == j
       !nb_elemD(2) = n-2*p - 1
    end if
    
    k = 0 !1
    if (dim_map >2) then
       p = Jpqr_map(3)
       n =  Nkv_map(3)
       call BINARYSEARCH(Xi(3), ukv3_map(p+1:n-(p+1)), n-2*p, k)
       k = k + p
       
       mask(:,3) = Nijk_map(3,:) == k
       !nb_elemD(3) = n-2*p - 1
    end if
    
    nbfound = COUNT(ALL(MASK,2))
    if (nbfound /= 1) then
       print*,'Error in updating element number for XI:',XI
       print*,'-- variable current_elem set to 1'
       current_map_elem = 1
    else
       numel(:) = PACK((/(i, i=1,nbel_map,1)/),ALL(MASK,2))
       current_map_elem = numel(1)
    endif
    
    !current_map_elem = i
    !current_map_elem = current_map_elem + (j-1)*nb_elemD(1)
    !current_map_elem = current_map_elem + (k-1)*nb_elemD(1)*nb_elemD(2)
    !print*,'test mapping',current_map_elem
    weight_map_elem(:) = weight_map(:,current_map_elem)
    
  End subroutine updateMapElementNumber
  
END MODULE embeddedMapping
