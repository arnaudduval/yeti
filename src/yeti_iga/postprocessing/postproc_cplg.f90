!! Copyright 2021 Arnaud Duval

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

subroutine post_cplg(FILENAME, FieldOutput_flag, nb_ref, SOL, COORDS3D,&
                  &  IEN, nb_elem_patch, Nkv, Ukv, Nijk, weight, Jpqr, &
                  &  ELT_TYPE, PROPS, JPROPS, MATERIAL_PROPERTIES, TENSOR, &
                  &  NBPINT, NNODE, nb_patch, nb_elem, nb_cp, MCRD)

    use parameters
    use nurbspatch
    use embeddedMapping
    
    implicit none
    
    Character(len=*) , intent(in) :: FILENAME
    integer, intent(in) :: nb_ref
    !! TODO flag is not necessary ?
    logical, intent(in) :: FieldOutput_flag
    dimension FieldOutput_flag(3)
    
    !! Nurbs geometry
    integer, intent(in) :: nb_cp, nb_patch, nb_elem
    Double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3,nb_cp)
    
    double precision, intent(in) :: Ukv, weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem), &
           &  Ukv(:), weight(:)
    
    !! patches and elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
    integer, intent(in) :: MCRD, NBPINT, NNODE
    integer, intent(in) :: IEN, nb_elem_patch
    integer, intent(in) :: JPROPS
    
    dimension IEN(:), nb_elem_patch(nb_patch), NNODE(nb_patch)
    dimension PROPS(:), MATERIAL_PROPERTIES(2, nb_patch)
    dimension JPROPS(nb_patch)
    
    !! Analysis solution      
    Double precision, intent(in) :: SOL
    dimension SOL(MCRD,nb_cp)
    
    
    
    !! ----- Local variables -------
    integer :: numPatch, numDomain, numLgrge, ismaster, dispOrRot, diminterface
    integer :: nb_eval_pt, npt1D
    double precision, dimension(:,:), allocatable :: XI, XIbar
    double precision, dimension(:,:,:), allocatable :: BI
    double precision, dimension(:), allocatable :: evalpt1D
    
    double precision :: COORDS_elem
    dimension COORDS_elem(MCRD, maxval(NNODE))
    double precision :: SOL_elem
    dimension SOL_elem(MCRD, maxval(NNODE))
    
    integer :: sctr
    dimension sctr(MAXVAL(NNODE))
    
    double precision pos, disp, rot
    dimension pos(MCRD), disp(MCRD), rot(MCRD)
    
    integer :: i, j, kk, n, num_elem, ipt
    character(len=7) :: filename_ext
    
    do numPatch = 1, nb_patch
        
        call extractNurbsPatchMechInfos(numPatch,IEN,PROPS,JPROPS, &
                   &  NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
        if (ELT_TYPE_patch == 'U00') then
            write(*,*) numPatch
            call extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
                    
            numDomain = int(PROPS_patch(2))
            numLgrge  = int(PROPS_patch(3))
            ismaster  = int(PROPS_patch(4))
            dimInterface = dim_patch
            
            ! Points for evaluation
            npt1D = (2**nb_ref+1)
            n = npt1D**dimInterface
            nb_eval_pt = nb_elem_patch(numPatch)*n
            if (allocated(XIbar)) deallocate(XIbar)
            if (allocated(XI))  deallocate(XI)
            if (allocated(BI))  deallocate(BI)
            
            allocate(XIbar(3, nb_eval_pt))    !! Gauss pt in knotvector
            allocate(XI(3, nb_eval_pt))     !! Gauss pt in interface
            allocate(BI(3, dimInterface, nb_eval_pt))  !! covariant basis vector
            
            XIbar(:,:) = zero
            XI(:,:) = zero
            BI(:,:,:) = zero
            
            kk=0
            do num_elem = 1, nb_elem_patch(numPatch)
                do i = 1, nnode_patch
                    COORDS_elem(:,i) = COORDS3D(:MCRD, IEN_patch(i, num_elem))
                enddo

                call extractNurbsElementInfos(num_elem)

                call getEvalPtsOnParamSpace(XIbar(:,kk+1:kk+n), XI(:, kk+1:kk+n), &
                             & BI(:,:, kk+1:kk+n), dim_patch, MCRD, nnode_patch, &
                             & n, COORDS_elem(:,:nnode_patch))

                kk = kk + n
            enddo
            
            ! Lagrange field
            call extractNurbsPatchGeoInfos(numLgrge, Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS, &
                    &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
            dispORrot = int(PROPS_patch(2))
            
            ! Domain to couple
            call extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv, &
                    &        weight,nb_elem_patch)
            call extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS, &
                    &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
            if (ELT_TYPE_patch == 'U30') then
                i = int(PROPS_patch(2))
                call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
                    &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
            endif
            

            write (filename_ext,"(A4,I1)") FILENAME, NumPatch



            open(90, file='results/'//trim(filename_ext)//'.res', form='formatted')
            do kk = 1, nb_eval_pt
                !!write(*,*) XI(:,kk)
                call updateElementNumber(XI(:,kk))
                sctr(:nnode_patch) = IEN_patch(:, current_elem)
                do i = 1, nnode_patch
                    COORDS_elem(:,i) = COORDS3D(:MCRD, sctr(i))
                    SOL_elem(:,i) = SOL(:MCRD, sctr(i))
                enddo
                call evalRot(pos, disp, rot, XI(:,kk), dim_patch, MCRD, nnode_patch, COORDS_elem, SOL_elem)
            
                write(90,*) XI(:,kk), pos(:), disp(:), rot(:)
            enddo
            close(90)
        endif
    enddo

end subroutine post_cplg

subroutine getEvalPtsOnParamSpace(XIbar, XI, BI, dim, MCRD, NNODE, nb_pt, COORDS)

    use parameters
    use nurbspatch
    
    implicit none
    
    !! Input arguments
    integer, intent(in) :: dim, MCRD, NNODE, nb_pt
    double precision, intent(in) :: COORDS
    dimension COORDS(MCRD, NNODE)
    
    !! Output variables
    double precision, intent(out) :: XIbar, XI, BI
    dimension XIbar(3, nb_pt), XI(3, nb_pt), BI(3, dim, nb_pt)
    
    !! Local variables
    double precision :: R, dRdxi
    dimension R(NNODE), dRdxi(NNODE, 3)
    double precision :: pt_coord
    dimension pt_coord(dim, nb_pt)
    integer :: nb_pt1D
    integer :: i, j, ipt
    
    XIbar(:,:) = zero
    XI(:,:)    = zero
    BI(:,:,:)  = zero
    
    nb_pt1D = int( nb_pt**(1.0/float(dim)) )
    if (nb_pt1D**dim < nb_pt) nb_pt1D = nb_pt1D + 1
    
    !! create grid in both directions between 0 and 1
    if (dim .eq. 1) then
        do i = 1, nb_pt1D
            pt_coord(1, i) = (i-1.)/(nb_pt1D-1.)
        enddo
    elseif (dim .eq. 2) then
        ipt = 1
        do i = 1, nb_pt1D
            do j =1, nb_pt1D
                pt_coord(1, ipt) = (i-1.)/(nb_pt1D-1.)
                pt_coord(2, ipt) = (j-1.)/(nb_pt1D-1.)
                ipt = ipt+1
            enddo
        enddo
    endif

    do ipt = 1, nb_pt
        ! get XIbar
        do i = 1, dim
            XIbar(i,ipt) = Ukv_elem(1,i) + pt_coord(i,ipt)*(Ukv_elem(2,i) - Ukv_elem(1,i))
        enddo

        ! evaluate basis functions
        call evalnurbs(XIbar(:3,ipt), R(:), dRdxi(:,:))

        ! get position
        do i = 1, NNODE
            XI(:MCRD,ipt) = XI(:MCRD,ipt) + R(i)*COORDS(:,i)
        enddo

        ! get covariant basis vectors
        do j = 1, dim
            do i=1, NNODE
                BI(:MCRD, j, ipt) = BI(:MCRD, j, ipt) + dRdxi(i,j)*COORDS(:,i)
            enddo
        enddo
    enddo
end subroutine getEvalPtsOnParamSpace

subroutine evalRot(pos, disp, rot, XI, dim, MCRD, NNODE, COORDS, SOL)

    use parameters
    use nurbspatch
    
    implicit none
    
    !! Input arguments
    integer, intent(in) :: dim, MCRD, NNODE
    double precision, intent(in) :: XI, COORDS, SOL
    dimension COORDS(MCRD, NNODE), SOL(MCRD, NNODE), XI(MCRD)
    
    !! Output variables => TBD
    double precision, intent(out) :: pos, disp, rot
    dimension pos(MCRD), disp(MCRD), rot(MCRD)
    
    !! Local variables
    double precision :: R, dRdxi
    dimension R(NNODE), dRdxi(NNODE, 3)
    double precision :: BI, A1, A2, A3, norm, dudxi
    dimension BI(3, dim), A1(3), A2(3), A3(3), dudxi(3,dim)
    integer :: i, j, numCP
    
    pos(:) = zero
    disp(:) = zero
    rot(:) = zero
    BI(:,:) = zero
    dudxi(:,:) = zero
    
    
    call evalnurbs(XI(:3), R(:), dRdxi(:,:))
      
    !! get position
    do i = 1, NNODE
        pos(:MCRD) = pos(:MCRD) + R(i)*COORDS(:,i)
    enddo
    
    !! get displacement
    do i = 1, NNODE
        disp(:MCRD) = disp(:MCRD) + R(i)*SOL(:,i)
    enddo
    
!     !! get covariant basis vector
!     do j = 1, dim
!         do i=1, NNODE
!             BI(:MCRD, j) = BI(:MCRD, j) + dRdxi(i,j)*COORDS(:,i)
!         enddo
!     enddo
!     !! get rotation
!     do j = 1, dim
!         do i=1, NNODE
!             rot(:MCRD) = rot(:MCRD) + dRdxi(i,j)*SOL(:,i)*BI(3,j)
!         enddo
!     enddo

    A1(:) = zero
    A2(:) = zero
    do numCP = 1, NNODE
        A1(:) = A1(:) + dRdxi(numCP,1)*COORDS(:,numCP)
        A2(:) = A2(:) + dRdxi(numCP,2)*COORDS(:,numCP)
    enddo
    call cross(A1(:), A2(:), A3(:))
    norm = sqrt(A3(1)*A3(1) + A3(2)*A3(2) + A3(3)*A3(3))
    A3(:) = A3(:)/norm
    
    do j = 1, dim
        do i=1, NNODE
            dudxi(:MCRD, j) = dudxi(:MCRD, j) + dRdxi(i,j)*SOL(:,i)
        enddo
    enddo
    
    rot(1) = dudxi(1,1)*A3(1) + dudxi(2,1)*A3(2) + dudxi(3,1)*A3(3)
    rot(2) = dudxi(1,2)*A3(1) + dudxi(2,2)*A3(2) + dudxi(3,2)*A3(3)
    
end subroutine evalRot
