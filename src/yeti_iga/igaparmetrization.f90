!! Copyright 2016-2020 Thibaut Hirschler

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

MODULE nurbspatch

  use parameters

  Implicit none

  ! Patch geo infos
  Integer :: current_patch, nnode_patch, dim_patch, nbel_patch
  !! TODO : Those arrays should allocatable to chosse their dimension
  Integer, dimension(3) :: Jpqr_patch, Nkv_patch
  Integer,          allocatable, dimension(:,:) :: Nijk_patch
  Double precision, allocatable, dimension(:)   :: Ukv1_patch,         &
       Ukv2_patch,Ukv3_patch
  Double precision, allocatable, dimension(:,:) :: weight_patch

  ! Patch mech infos
  Integer :: current_patch_mech, JPROPS_patch
  Integer,          allocatable, dimension(:,:) :: IEN_patch
  Double precision, allocatable, dimension(:)   :: PROPS_patch
  Character(:), allocatable  :: ELT_TYPE_patch,TENSOR_patch

  ! Element infos
  Integer :: current_elem
  Double precision, allocatable, dimension(:) :: weight_elem
  Double precision, dimension(2,3) :: Ukv_elem

CONTAINS

  Subroutine extractNurbsPatchGeoInfos(num_patch2extract,Nkv,Jpqr,Nijk,&
       Ukv, weight, nb_elem_patch)
    implicit none

    ! Declaration des variables ........................................

    ! Input variables :
    ! ---------------
    Integer, intent(in) :: num_patch2extract
    Integer,          dimension(:),   intent(in) :: nb_elem_patch
    Integer,          dimension(:,:), intent(in) :: Nkv, Jpqr, Nijk
    Double precision, dimension(:),   intent(in) :: Ukv, weight


    ! Local variables :
    ! ---------------
    Integer :: num_patch,num_elem,count,n,temp,i


    ! Fin declaration des variables ....................................
    !
    !
    ! Initalisation ....................................................

    ! Num. patch to extract
    current_patch = num_patch2extract

    ! Degree per direction
    Jpqr_patch(:) = Jpqr(:,num_patch2extract)

    ! Knot vector sizes
    Nkv_patch(:)  = Nkv(:,num_patch2extract)
    dim_patch = 1
    If (Nkv_patch(2)>0) dim_patch = dim_patch + 1
    If (Nkv_patch(3)>0) dim_patch = dim_patch + 1

    ! Knot vectors
    count = 1
    Do num_patch = 1,num_patch2extract-1
       count = count + SUM(Nkv(:,num_patch))
    End Do
    ! - xi
    n = Nkv_patch(1)
    if (allocated(Ukv1_patch)) deallocate(Ukv1_patch)
    allocate(Ukv1_patch(n))
    Ukv1_patch(:) = Ukv(count:count+n-1)
    count = count+n
    ! - eta
    if (dim_patch>1) then
       n = Nkv_patch(2)
       if (allocated(Ukv2_patch)) deallocate(Ukv2_patch)
       allocate(Ukv2_patch(n))
       Ukv2_patch(:) = Ukv(count:count+n-1)
       count = count+n
    end if
    ! - zeta
    if (dim_patch>2) then
       n = Nkv_patch(3)
       if (allocated(Ukv3_patch)) deallocate(Ukv3_patch)
       allocate(Ukv3_patch(n))
       Ukv3_patch(:) = Ukv(count:count+n-1)
    end if


    ! Element knot position
    n = nb_elem_patch(num_patch2extract)
    if (allocated(Nijk_patch)) deallocate(Nijk_patch)
    allocate( Nijk_patch(3,n) )

    nbel_patch = n

    count = 0
    Do num_patch = 1,num_patch2extract-1
       count = count + nb_elem_patch(num_patch)
    End Do
    Do num_elem = 1,n
       Nijk_patch(:,num_elem) = Nijk(:,count+num_elem)
    End Do

    ! Element weights
    nnode_patch = 1
    Do i = 1,3
       nnode_patch = nnode_patch*(Jpqr_patch(i)+1)
    End Do
    if (allocated(weight_patch)) deallocate(weight_patch)
    allocate(weight_patch(nnode_patch,n))
    if (allocated(weight_elem))  deallocate(weight_elem)
    allocate(weight_elem(nnode_patch))

    count = 0
    Do num_patch = 1,num_patch2extract-1
       temp = 1
       Do i = 1,3
          temp = temp*(Jpqr(i,num_patch)+1)
       End Do
       count = count + nb_elem_patch(num_patch)*temp
    End Do
    Do num_elem = 1,n
       weight_patch(:,num_elem) = weight(count+1:count+nnode_patch)
       count = count + nnode_patch
    End Do

  End Subroutine extractNurbsPatchGeoInfos






  Subroutine extractNurbsPatchMechInfos(num_patch2extract,IEN,PROPS,   &
    JPROPS,NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
    Implicit none

    ! Declaration des variables ........................................

    ! Input variables :
    ! ---------------
    Integer, intent(in) :: num_patch2extract
    Integer,          dimension(:), intent(in) :: IEN,JPROPS, NNODE, nb_elem_patch
    Double precision, dimension(:), intent(in) :: PROPS
    Character(len=*),               intent(in) :: ELT_TYPE, TENSOR

    ! Local variables :
    ! ---------------
    Integer :: count,n,num_patch,num_elem,i,j

    ! Fin declaration des variables ....................................
    !
    !
    ! Initalisation ....................................................

    ! Num. patch to extract
    current_patch_mech = num_patch2extract

    ! Properties
    JPROPS_patch = JPROPS(num_patch2extract)
    count = 0
    Do num_patch = 1,num_patch2extract-1
       count = count + JPROPS(num_patch)
    End Do
    if (allocated(PROPS_patch)) deallocate(PROPS_patch)
    allocate(PROPS_patch(JPROPS_patch))
    PROPS_patch(:) = PROPS(count+1:count+JPROPS_patch)

    ! IEN
    n = NNODE(num_patch2extract)
    if (allocated(IEN_patch)) deallocate(IEN_patch)
    allocate(IEN_patch(n,nb_elem_patch(num_patch2extract)))
    count = 0
    Do num_patch = 1,num_patch2extract-1
       count = count + NNODE(num_patch)*nb_elem_patch(num_patch)
    End Do
    Do num_elem = 1,nb_elem_patch(num_patch2extract)
       IEN_patch(:,num_elem) = IEN(count+1:count+n)
       count = count+n
    End Do

    ! ELT_TYPE
    count = 0
    Do num_patch = 1,num_patch2extract
       n = count
       i = INDEX(ELT_TYPE(count:),    'U')-1
       j = INDEX(ELT_TYPE(count+i+1:),'U')-1
       if (j<0) then
          j=LEN(ELT_TYPE)-count-i
       else
          count = count+i+j
       end if
    Enddo
    if (allocated(ELT_TYPE_patch)) deallocate(ELT_TYPE_patch)
    ELT_TYPE_patch = ELT_TYPE(n+i:n+i+j)

    ! TENSOR
    count = 0
    Do num_patch = 1,num_patch2extract
       n = count
       i = INDEX(TENSOR(count:),    '/')-1
       j = INDEX(TENSOR(count+i+1:),'/')-1
       if (j<0) then
          j=LEN(TENSOR)-count-i
       else
          count = count+i+j
       end if
    Enddo
    if (allocated(TENSOR_patch)) deallocate(TENSOR_patch)
    TENSOR_patch = TENSOR(n+i+1:n+i+j)


  End Subroutine extractNurbsPatchMechInfos



  Subroutine finalizeNurbsPatch()
    implicit none
    ! geo
    if (allocated(Ukv1_patch))   deallocate(Ukv1_patch)
    if (allocated(Ukv2_patch))   deallocate(Ukv2_patch)
    if (allocated(Ukv3_patch))   deallocate(Ukv3_patch)
    if (allocated(Nijk_patch))   deallocate(Nijk_patch)
    if (allocated(weight_patch)) deallocate(weight_patch)
    if (allocated(weight_elem))  deallocate(weight_elem)
    ! mech
    if (allocated(PROPS_patch))  deallocate(PROPS_patch)
    if (allocated(IEN_patch))    deallocate(IEN_patch)

  End Subroutine finalizeNurbsPatch



  ! ....................................................................
  !
  ! Element level


  Subroutine extractNurbsElementInfos(num_elem)
    implicit none
    integer, intent(in)   :: num_elem
    integer, dimension(3) :: Ni
    integer :: i

    current_elem   = num_elem
    weight_elem(:) = weight_patch(:,num_elem)

    Ni(:) = Nijk_patch(:,num_elem)
    Ukv_elem(:,:) = zero
    Ukv_elem(:,1) = Ukv1_patch(Ni(1):Ni(1)+1)
    if (dim_patch>1) Ukv_elem(:,2) = Ukv2_patch(Ni(2):Ni(2)+1)
    if (dim_patch>2) Ukv_elem(:,3) = Ukv3_patch(Ni(3):Ni(3)+1)

  End Subroutine extractNurbsElementInfos



  Subroutine INRANGE(val,left,right, bool)
    implicit none
    Double precision, intent(in)  :: val,left,right
    Logical,          intent(out) :: bool

    bool = .FALSE.
    IF ((val .GE. left) .AND. (val .LT. right)) bool = .TRUE.
  End Subroutine INRANGE

  Subroutine BINARYSEARCH(val,uknot,nkv, index)
    implicit none
    Integer,          intent(in) :: nkv
    Double precision, intent(in) :: val,uknot
    dimension uknot(nkv)
    Integer,          intent(out):: index
    Integer :: low,high,middle, found
    Logical :: test

    low=1;high=nkv
    found = 0
    middle= 1
    Do while (low.NE.high  .AND. found.EQ.0)
       middle = (low+high)/2

       CALL INRANGE(val,uknot(middle),uknot(middle+1), test)
       If (test) then
          found = 1
       Elseif (val < uknot(middle)) then
          high = middle
       Else
          low  = middle+1
       End If

    End Do
    index = middle

  End Subroutine BINARYSEARCH


  Subroutine updateElementNumber(Xi)
    implicit none
    Double precision, dimension(3), intent(in) :: Xi

    Integer :: i,j,k,p,n,numel,nbfound
    Logical :: mask
    dimension mask(nbel_patch,3), numel(1)

    ! ...

    mask(:,:) = .TRUE.

    p = Jpqr_patch(1)
    n =  Nkv_patch(1)
    call BINARYSEARCH(Xi(1), ukv1_patch(p+1:n-(p+1)), n-2*p, i)
    i = i + p

    mask(:,1) = Nijk_patch(1,:) == i

    j = 0
    if (dim_patch >1) then
       p = Jpqr_patch(2)
       n =  Nkv_patch(2)
       call BINARYSEARCH(Xi(2), ukv2_patch(p+1:n-(p+1)), n-2*p, j)
       j = j + p

       mask(:,2) = Nijk_patch(2,:) == j

    end if

    k = 0
    if (dim_patch >2) then
       p = Jpqr_patch(3)
       n =  Nkv_patch(3)
       call BINARYSEARCH(Xi(3), ukv3_patch(p+1:n-(p+1)), n-2*p, k)
       k = k + p

       mask(:,3) = Nijk_patch(3,:) == k

    end if

    nbfound = COUNT(ALL(MASK,2))
    if (nbfound /= 1) then
       print*,'Error in updating element number for XI:',XI
       print*,'-- variable current_elem set to 1'
       call extractNurbsElementInfos(1)
    else
       numel(:) = PACK((/(i, i=1,nbel_patch,1)/),ALL(MASK,2))
       call extractNurbsElementInfos(numel(1))
    endif

  End subroutine updateElementNumber






END MODULE nurbspatch
