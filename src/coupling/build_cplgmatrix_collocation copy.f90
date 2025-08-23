subroutine cplg_matrix_collocation(xi_master, xi_slave, xi_inter, i_master, i_slave, &
            &  i_inter, d, ien, nb_elem_patch, elt_type, tensor, &
            &  props, jprops, nnode, nb_patch, nb_elem, &
            &  nkv, ukv, nijk, weight, jpqr, nb_cp, coords3d, mcrd, nb_data1, nb_data2, &
            &  Cdata1, Crow1, Ccol1, Cdata2, Crow2, Ccol2)

use parameters
use nurbspatch

implicit none

!! Input arguments
!! ---------------

!! NURBS geometry
integer, intent(in) :: nb_cp
double precision, intent(in) :: COORDS3D
dimension COORDS3D(3,nb_cp)

integer, intent(in) :: mcrd, i_master, i_slave, i_inter, nb_data1, nb_data2

double precision, intent(in), dimension(3) :: xi_master, xi_slave, xi_inter
double precision, intent(in) :: d

double precision, intent(in) :: Ukv, weight
integer, intent(in) :: Nkv, Jpqr, Nijk
dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),   &
    &     Ukv(:),weight(:)

!! Patches and Elements
character(len=*), intent(in) :: TENSOR, ELT_TYPE
double precision, intent(in) :: PROPS
integer, intent(in) :: NNODE,nb_elem,IEN,nb_patch,   &
    &   nb_elem_patch, JPROPS
dimension PROPS(:),   &
    &   NNODE(nb_patch),    &
    &   IEN(:), &
    &   nb_elem_patch(nb_patch),    &
    &   JPROPS(nb_patch)


!! Output variables
!! ----------------
integer,          intent(out) :: Crow1,Ccol1,Crow2,Ccol2
dimension Crow1(nb_data1),Ccol1(nb_data1),Crow2(nb_data2),Ccol2(nb_data2)

double precision, intent(out) :: Cdata1, Cdata2
dimension Cdata1(nb_data1), Cdata2(nb_data2)

!! Local variables
    !! ---------------
! double precision, dimension(mcrd,mcrd) :: dRdxi_master
! double precision, dimension(mcrd) :: R_master

integer :: icp, idim, nb_gps, idof, jdof
integer :: i, j, k, l, m, count
integer :: NNODE_1, NNODE_2, NNODE_l

!! Manage infos from projection algorithm
integer :: info

double precision, allocatable :: COORDS_elem(:,:)
integer, allocatable          :: sctr_1(:), sctr_2(:), sctr_l(:)


double precision, dimension(:,:), allocatable :: dRdxi_master, dRdxi_slave, dRdxi_lgrge
double precision, dimension(:,:), allocatable :: dxdxi_master, dxdxi_slave, dxdxi_lgrge

double precision, dimension(:,:), allocatable :: dxidx_master, dxidx_slave, dxidx_lgrge

double precision :: det_dxdxi_master, det_dxdxi_slave, det_dxdxi_lgrge 

double precision, dimension(:), allocatable :: R_master, R_slave, R_lgrge

double precision, dimension(nb_data1) :: CMAT1
double precision, dimension(nb_data2) :: CMAT2

double precision, dimension(3) :: xi_slave_inv, xi_inter_inv

allocate(COORDS_elem(MCRD,MAXVAL(NNODE)))
allocate(sctr_1(MAXVAL(NNODE)), sctr_2(MAXVAL(NNODE)), sctr_l(MAXVAL(NNODE)))





! Master
! ------------------------------

call extractNurbsPatchMechInfos(i_master,IEN,PROPS,JPROPS,  &
&   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

call extractNurbsPatchGeoInfos(i_master,Nkv,Jpqr,Nijk,Ukv, &
&        weight,nb_elem_patch)

if (allocated(R_master)) deallocate(R_master)
allocate(R_master(nnode_patch))
            !! Allocate array for shape function derivatives
if (allocated(dRdxi_master)) deallocate(dRdxi_master)
allocate(dRdxi_master(nnode_patch, 3))
if (allocated(dxdxi_master)) deallocate(dxdxi_master)
allocate(dxdxi_master(3, 3))

if (allocated(dxidx_master)) deallocate(dxidx_master)
allocate(dxidx_master(mcrd, mcrd))

NNODE_1 = nnode_patch

call updateElementNumber(xi_master)

write(*,*) 'Current_Elem_master'
write(*,*) current_elem

sctr_1(:NNODE_1) = IEN_patch(:,current_elem)

call evalnurbs(xi_master, R_master, dRdxi_master)

write(*,*) 'xi_master', xi_master

write(*,*) '**********************'

! write(*,*) 'R_master'

! write(*,*) R_master

! write(*,*) '**********************'

! write(*,*) dRdxi_master

! write(*,*) '**********************'

dxdxi_master(:,:) = 0.0
do icp = 1, nnode_patch
    COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
    do idim = 1, MCRD
        dxdxi_master(idim, :) = dxdxi_master(idim, :) + &
            &  dRdxi_master(icp, :)*COORDS_elem(idim, icp)
    enddo
enddo

! write(*,*) dxdxi_master

dxidx_master(:,:) = zero
call MatrixInv(dxidx_master(:mcrd,:mcrd), dxdxi_master(:mcrd,:mcrd), det_dxdxi_master, mcrd)

! write(*,*) '**********************'

! write(*,*) dxidx_master



! Slave
! ------------------------------

call extractNurbsPatchMechInfos(i_slave,IEN,PROPS,JPROPS,  &
&   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

call extractNurbsPatchGeoInfos(i_slave,Nkv,Jpqr,Nijk,Ukv, &
&        weight,nb_elem_patch)

if (allocated(R_slave)) deallocate(R_slave)
allocate(R_slave(nnode_patch))
            !! Allocate array for shape function derivatives
if (allocated(dRdxi_slave)) deallocate(dRdxi_slave)
allocate(dRdxi_slave(nnode_patch, 3))
if (allocated(dxdxi_slave)) deallocate(dxdxi_slave)
allocate(dxdxi_slave(3, 3))

if (allocated(dxidx_slave)) deallocate(dxidx_slave)
allocate(dxidx_slave(mcrd, mcrd))

call point_inversion_surface(xi_master, 5, COORDS3D, nb_cp,   &
    &   .false., xi_slave_inv, info)

NNODE_2 = nnode_patch

call updateElementNumber(xi_slave_inv)

sctr_2(:NNODE_2) = IEN_patch(:,current_elem)

call evalnurbs(xi_slave_inv, R_slave, dRdxi_slave)

! write(*,*) current_elem

write(*,*) 'xi_slave', xi_slave

write(*,*) '**********************'

write(*,*) 'xi_slave_inv', xi_slave_inv

write(*,*) '**********************'

! write(*,*) R_slave

! write(*,*) '**********************'

! write(*,*) dRdxi_slave

! write(*,*) '**********************'

dxdxi_slave(:,:) = 0.0
do icp = 1, nnode_patch
    COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
    do idim = 1, MCRD
        dxdxi_slave(idim, :) = dxdxi_slave(idim, :) + &
            &  dRdxi_slave(icp, :)*COORDS_elem(idim, icp)
    enddo
enddo

! write(*,*) dxdxi_slave


dxidx_slave(:,:) = zero
call MatrixInv(dxidx_slave(:mcrd,:mcrd), dxdxi_slave(:mcrd,:2), det_dxdxi_slave, mcrd)

! write(*,*) '**********************'
! write(*,*) dxidx_slave



! Interface
! ------------------------------

call extractNurbsPatchMechInfos(i_inter,IEN,PROPS,JPROPS,  &
&   NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

call extractNurbsPatchGeoInfos(i_inter,Nkv,Jpqr,Nijk,Ukv, &
&        weight,nb_elem_patch)

allocate(R_lgrge(nnode_patch))

            !! Allocate array for shape function derivatives
if (allocated(dRdxi_lgrge)) deallocate(dRdxi_lgrge)
allocate(dRdxi_lgrge(nnode_patch, 3))
if (allocated(dxdxi_lgrge)) deallocate(dxdxi_lgrge)
allocate(dxdxi_lgrge(3, 3))

if (allocated(dxidx_lgrge)) deallocate(dxidx_lgrge)
allocate(dxidx_lgrge(mcrd, mcrd))

NNODE_l = nnode_patch

call point_inversion_surface(xi_master, 5, COORDS3D, nb_cp,   &
    &   .false., xi_inter_inv, info)

call updateElementNumber(xi_inter_inv)
call evalnurbs(xi_inter_inv, R_lgrge, dRdxi_lgrge)

sctr_l(:NNODE_l) = IEN_patch(:,current_elem)

! write(*,*) current_elem

write(*,*) 'xi_inter', xi_inter

write(*,*) '**********************'

write(*,*) 'xi_inter_inv', xi_inter_inv

write(*,*) '**********************'

! write(*,*) R_lgrge

! write(*,*) '**********************'

! write(*,*) dRdxi_lgrge

! write(*,*) '**********************'

dxdxi_lgrge(:,:) = 0.0
do icp = 1, nnode_patch
    COORDS_elem(:, icp) = COORDS3D(:mcrd, IEN_patch(icp, current_elem))
    do idim = 1, MCRD
        dxdxi_lgrge(idim, :) = dxdxi_lgrge(idim, :) + &
            &  dRdxi_lgrge(icp, :)*COORDS_elem(idim, icp)
    enddo
enddo

! write(*,*) dxdxi_lgrge

dxidx_lgrge(:,:) = zero
call MatrixInv(dxidx_lgrge(:mcrd,:mcrd), dxdxi_lgrge(:mcrd,:mcrd), det_dxdxi_lgrge, mcrd)

! write(*,*) '**********************'
! write(*,*) dxidx_lgrge

! write(*,*) sctr_l
count = 1

do j = 1, NNODE_l
    jdof = (sctr_l(j)-1)*MCRD 
    do i = 1, NNODE_1
        idof = (sctr_1(i)-1)*MCRD
        CMAT1(count) = R_lgrge(j)*R_master(i)
        do m = 1, mcrd
            do l = 1, mcrd
                do k = 1,mcrd
                    CMAT1(count) = CMAT1(count) + &  
                    & ((d**2)/8)* (dRdxi_lgrge(j,m) * dxidx_lgrge(m,l) * dRdxi_master(i,k) * dxidx_master(k,l)) 
                enddo
            enddo
        enddo
        ! write(*,*) 'idof = ', idof
        ! write(*,*) 'jdof = ', jdof
        CMAT1(count+1) = CMAT1(count)
        Cdata1(count) = CMAT1(count)
        Cdata1(count+1) = CMAT1(count+1)
        Crow1(count) = idof 
        Ccol1(count) = jdof 
        Crow1(count+1) = idof + 1
        Ccol1(count+1) = jdof + 1  
        count = count + 2 
    enddo
enddo

! write(*,*) CMAT1

count = 1

! FACTORISER
do j = 1, NNODE_l
    jdof = (sctr_l(j)-1)*MCRD
    do i = 1, NNODE_2
        idof = (sctr_2(i)-1)*MCRD
        CMAT2(count) = -(R_lgrge(j)*R_slave(i))
        do m = 1,mcrd
            do l = 1, mcrd
                do k = 1, mcrd
                    CMAT2(count) = CMAT2(count) - &
                    & ((d**2)/8)* (dRdxi_lgrge(j,m) * dxidx_lgrge(m,l) * dRdxi_slave(i,k) * dxidx_slave(k,l)) 
                enddo
            enddo
        enddo
        CMAT2(count + 1) = CMAT2(count)
        Cdata2(count) = CMAT2(count)
        Cdata2(count+1) = CMAT2(count+1)
        Crow2(count) = idof 
        Ccol2(count) = jdof 
        Crow2(count+1) = idof + 1
        Ccol2(count+1) = jdof + 1  
        count = count + 2
    enddo
enddo




end subroutine cplg_matrix_collocation

