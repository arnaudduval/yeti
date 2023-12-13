!! Copyright 2016-2020 Thibaut Hirschler
!! Copyright 2023 Arnaud Duval


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


module m_nurbspatch
    implicit none

    type :: nurbspatch
        !! Geometry
        integer :: current_patch
        integer :: nnode_patch
        integer :: dim_patch
        integer :: nbel_patch
        integer, dimension(3) :: Jpqr_patch
        integer, dimension(3) :: Nkv_patch
        integer, allocatable, dimension(:,:) :: Nijk_patch
        double precision, allocatable, dimension(:) :: Ukv1_patch
        double precision, allocatable, dimension(:) :: Ukv2_patch
        double precision, allocatable, dimension(:) :: Ukv3_patch
        double precision, allocatable, dimension(:,:) :: weight_patch

        !! Mechanics
        integer :: current_patch_mech
        integer :: jprops_patch
        integer, allocatable, dimension(:,:) :: ien_patch
        double precision, allocatable, dimension(:) :: props_patch
        character(:), allocatable :: elt_type_patch, tensor_patch

    contains
        procedure :: extractNurbsPatchGeoInfos
        procedure :: extractNurbsPatchMechInfos
        procedure :: finalizeNurbsPatch
    end type

    type :: nurbselement
        integer :: current_elem
        double precision, allocatable, dimension(:) :: weight_elem
        double precision, dimension(2,3) :: ukv_elem

    contains
        procedure :: extractNurbsElementInfos
        procedure :: print
    end type

contains

    subroutine extractNurbsPatchGeoInfos(self, num_patch2extract, Nkv, Jpqr, Nijk, &
        Ukv, weight, nb_elem_patch)

        class(nurbspatch), intent(inout) :: self
        integer, intent(in) :: num_patch2extract
        integer, dimension(:), intent(in) :: nb_elem_patch
        integer, dimension(:,:), intent(in) :: Nkv, Jpqr, Nijk
        double precision, dimension(:),   intent(in) :: Ukv, weight

        integer :: count, ipatch, ielem, n, i, temp

        !! Patch ID
        self%current_patch = num_patch2extract

        !! Degree per parametric direction
        self%Jpqr_patch = Jpqr(:, num_patch2extract)

        !! Knot vectors sizes
        self%Nkv_patch(:) = Nkv(:, num_patch2extract)
        self%dim_patch = 1
        if (self%Nkv_patch(2) > 0) self%dim_patch = self%dim_patch + 1
        if (self%Nkv_patch(3) > 0) self%dim_patch = self%dim_patch + 1

        !! Knot vectors
        count = 1
        do ipatch = 1, num_patch2extract-1
            count = count + sum(Nkv(:, ipatch))
        enddo

        !! xi
        n = self%Nkv_patch(1)
        if (allocated(self%Ukv1_patch)) deallocate(self%Ukv1_patch)
        allocate(self%Ukv1_patch(n))
        self%Ukv1_patch(:) = Ukv(count:count+n-1)
        count = count + n
        !! eta
        if (self%dim_patch > 1) then
            n = self%Nkv_patch(2)
            if (allocated(self%Ukv2_patch)) deallocate(self%Ukv2_patch)
            allocate(self%Ukv2_patch(n))
            self%Ukv2_patch(:) = Ukv(count:count+n-1)
            count = count + n
        endif
        !! zeta
        if (self%dim_patch > 2) then
            n = self%Nkv_patch(3)
            if (allocated(self%Ukv3_patch)) deallocate(self%Ukv3_patch)
            allocate(self%Ukv3_patch(n))
            self%Ukv3_patch(:) = Ukv(count:count+n-1)
            count = count + n
        endif

        !! Element knot positions
        n = nb_elem_patch(num_patch2extract)
        if (allocated(self%Nijk_patch)) deallocate(self%Nijk_patch)
        allocate(self%Nijk_patch(3,n))

        self%nbel_patch = n

        count = 0
        do ipatch = 1, num_patch2extract-1
            !! TODO : this loop could be replaced by sum() function
            count = count + nb_elem_patch(ipatch)
        enddo
        do ielem = 1, n
            self%Nijk_patch(:, ielem) = Nijk(:, count+ielem)
        enddo

        !! Element weights
        self%nnode_patch = 1
        do i=1, 3
            self%nnode_patch = self%nnode_patch*(self%Jpqr_patch(i)+1)
        enddo
        if (allocated(self%weight_patch)) deallocate(self%weight_patch)
        allocate(self%weight_patch(self%nnode_patch, n))

        count = 0
        do ipatch = 1, num_patch2extract-1
            temp = 1
            do i = 1, 3
                temp = temp*(Jpqr(i, ipatch)+1)
            enddo
            count = count + nb_elem_patch(ipatch)*temp
        enddo
        do ielem = 1, n
            self%weight_patch(:, ielem) = weight(count+1:count+self%nnode_patch)
            count = count + self%nnode_patch
        enddo

    end subroutine extractNurbsPatchGeoInfos

    subroutine extractNurbsPatchMechInfos(self, num_patch2extract, ien, props, &
        jprops, nnode, nb_elem_patch, elt_type, tensor)

        implicit none

        class(nurbspatch), intent(inout) :: self
        integer, intent(in) :: num_patch2extract
        integer, dimension(:), intent(in) :: ien, jprops, nnode, nb_elem_patch
        double precision, dimension(:), intent(in) :: props
        character(len=*), intent(in) :: elt_type, tensor

        integer :: count, ipatch, ielem, n, i, j

        !! ID of patch to extract
        self%current_patch_mech = num_patch2extract

        !! Properties
        self%jprops_patch = jprops(num_patch2extract)
        count = 0
        do ipatch = 1, num_patch2extract-1
            !! TODO : this loop could be replaced by sum() function
            count = count + jprops(ipatch)
        enddo
        if (allocated(self%props_patch)) deallocate(self%props_patch)
        allocate(self%props_patch(self%jprops_patch))
        self%props_patch(:) = props(count+1:count+self%jprops_patch)

        !! Connectivity table
        n = nnode(num_patch2extract)
        if (allocated(self%ien_patch)) deallocate(self%ien_patch)
        allocate(self%ien_patch(n, nb_elem_patch(num_patch2extract)))
        count = 0
        do ipatch = 1, num_patch2extract-1
            count = count + nnode(ipatch)*nb_elem_patch(ipatch)
        enddo
        do ielem = 1, nb_elem_patch(num_patch2extract)
            self%ien_patch(:, ielem) = ien(count+1:count+n)
            count = count + n
        enddo

        !! Element type
        count = 0
        do ipatch = 1, num_patch2extract
            n = count
            i = index(elt_type(count:), 'U')-1
            j = index(elt_type(count+i+1:), 'U')-1
            if (j < 0) then
                j = len(elt_type) - count - i
            else
                count = count + i + j
            endif
        enddo
        if (allocated(self%elt_type_patch)) deallocate(self%elt_type_patch)
        self%elt_type_patch = elt_type(n+i:n+i+j)

        !! Tensor type
        count = 0
        do ipatch = 1, num_patch2extract
            n = count
            i = index(tensor(count:), '/')-1
            j = index(tensor(count+i+1:), '/')-1
            if (j<0) then
                j=len(tensor) - count - i
            else
                count = count + i + j
            endif
        enddo
        if (allocated(self%tensor_patch)) deallocate(self%tensor_patch)
        self%tensor_patch = tensor(n+i+1:n+i+j)

    end subroutine extractNurbsPatchMechInfos


    subroutine finalizeNurbsPatch(self)

        implicit none

        class(nurbspatch), intent(inout) :: self

        if (allocated(self%Ukv1_patch))   deallocate(self%Ukv1_patch)
        if (allocated(self%Ukv2_patch))   deallocate(self%Ukv2_patch)
        if (allocated(self%Ukv3_patch))   deallocate(self%Ukv3_patch)
        if (allocated(self%Nijk_patch))   deallocate(self%Nijk_patch)
        if (allocated(self%weight_patch)) deallocate(self%weight_patch)
        if (allocated(self%props_patch))  deallocate(self%props_patch)
        if (allocated(self%ien_patch))    deallocate(self%ien_patch)

    end subroutine finalizeNurbsPatch

    subroutine finalizeNurbsElement(self)

        implicit none

        class(nurbselement), intent(inout) :: self

        if (allocated(self%weight_elem))  deallocate(self%weight_elem)

    end subroutine finalizeNurbsElement

    subroutine extractNurbsElementInfos(self, num_elem, patch)
        class(nurbselement), intent(inout) :: self
        integer, intent(in) :: num_elem
        type(nurbspatch), intent(in) :: patch

        integer, dimension(3) :: Ni
        integer :: i

        if (allocated(self%weight_elem))  deallocate(self%weight_elem)
        allocate(self%weight_elem(patch%nnode_patch))

        self%current_elem = num_elem
        self%weight_elem(:) = patch%weight_patch(:, num_elem)

        Ni(:) = patch%Nijk_patch(:, num_elem)
        self%ukv_elem(:,:) = 0.0
        self%ukv_elem(:, 1) = patch%Ukv1_patch(Ni(1):Ni(1)+1)
        if (patch%dim_patch>1) self%ukv_elem(:, 2) = patch%Ukv1_patch(Ni(2):Ni(2)+1)
        if (patch%dim_patch>2) self%ukv_elem(:, 3) = patch%Ukv1_patch(Ni(3):Ni(3)+1)

    end subroutine extractNurbsElementInfos

    subroutine print(self)
        class(nurbselement), intent(inout) :: self

        write(*,*) '--------'
        write(*,'(I7)') self%current_elem
        write(*,'(27F4.1)') self%weight_elem
        write(*,'(6F4.1)') self%ukv_elem
        write(*,*) '--------'
    end subroutine print
end module m_nurbspatch