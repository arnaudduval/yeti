!!! Find the index of the face of hull object corresponding
!!! to the interface location (given on the embedded entity)

subroutine find_hull_face(iface, coords3D, nb_cp, ifaceHull)
    use parameters
    use nurbspatch
    use embeddedMapping
    
    implicit none

    !! Input variables
    integer, intent(in) :: iface
    double precision, intent(in) :: coords3D
    dimension coords3D(3,nb_cp)
    integer, intent(in) :: nb_cp

    !! Output variable
    integer, intent(out) :: ifaceHull

    !! Local variables
    double precision :: u, v  ! Parameters values to evaluate embedded entity
    double precision :: lst  ! List of parameters to evaluate in each direction
    dimension lst(3)
    double precision :: xi  ! Parameter to evaluate basisfuns
    dimension xi(3)
    double precision :: coords  ! Coords
    dimension coords(3,nnode_patch)
    double precision :: R  ! Basisfuns
    dimension R(nnode_patch)
    double precision :: pt  ! Evaluated point
    dimension pt(3)
    double precision :: S  ! List of evaluated surface points
    dimension S(9, 3)
    integer :: count, i_u, i_v, icp, i_dim  ! For loops
    double precision :: tol1 ! Tolerance for mask
    double precision :: tol2  ! Tolerance for mask
    logical :: mask  ! Mask for param. direction and side
    dimension mask(9)
    integer :: direction  ! Parametric direction
    integer :: side  ! Side of parametric direction
    integer :: face  ! Possible face ids.
    dimension :: face(3, 2)


    !! Create list of evaluated points
        
    lst = (/ 0.0, 0.5, 1.0 /) 
    pt(:) = zero
    count = 1
    S(:, :) = zero
    
    do i_v = 1, 3
        v = lst(i_v)  ! Value v-parameter
        do i_u = 1, 3
            u = lst(i_u)  ! Value u-parameter
            !! Generate 3D point to evaluate volume
            call point_on_solid_face(u,v,iface,xi)
            !! Get knot span and CP coordinates
            call updateElementNumber(xi)
            do icp = 1, nnode_patch
                coords(:,icp) = COORDS3D(:3,IEN_patch(icp,current_elem))
            enddo
            !! Evaluate basisfuns
            call evalnurbs_noder(xi, R)
            !! Compute surface value
            do icp = 1, nnode_patch
                pt(:) = pt(:) + R(icp)*coords(:, icp)
            enddo
            S(count, :) = pt(:)
            count = count + 1
        enddo
    enddo
    
    !! Build mask to find which parametric direction is involved and at which side
    
    face(1, :) = (/ 1, 2 /)
    face(2, :) = (/ 3, 4 /)
    face(3, :) = (/ 5, 6 /)
    
    tol1 = 1.D-6
    tol2 = 1 - tol1
    
    direction = 0
    
    do i_dim = 1, 3  ! Check if side = 0
        count = 1
        mask(:) = .false.
        do i_v = 1, 3 
            do i_u = 1, 3
                if (S(count, i_dim) .le. tol1) then
                    mask(count) = .true.
                endif
                count = count + 1
            enddo
        enddo
        if (all(mask)) then 
            direction = i_dim
            side = 0
        endif
    enddo
    
    if (direction .eq. 0) then  ! Check if side = 1
        do i_dim = 1, 3
            count = 1
            mask(:) = .false.
            do i_v = 1, 3 
                do i_u = 1, 3
                    if (S(count, i_dim) .ge. tol2) then
                        mask(count) = .true.
                    endif
                    count = count + 1
                enddo
            enddo
            if (all(mask)) then 
                direction = i_dim
                side = 1
            endif
        enddo
    endif

    if (direction .eq. 0) then
        write(*,*) "Warning : did not find face number for hull. Default face id is set to 6."
        direction = 3
        side = 1
    endif
    
    !! Determine face number  
    ifaceHull = face(direction, side + 1)
    
    end subroutine find_hull_face
    