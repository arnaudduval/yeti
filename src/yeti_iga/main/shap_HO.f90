! Compute shape functions and their 1st,2nd and 3rd derivative at given point
!! Given point has coordinates in isoparametric element [-1 1]
!! WARNING : Only B-Spline functions are used ( Weight=1)
!! WARNING : ONLY 2D analysis

subroutine shap_HO(R, dRdx, d2Rdx, d3Rdx, DetJac, COORDS, PtGauss, MCRD)

    use parameters
    use nurbspatch

    implicit none

    !! Input arguments
    !! ---------------
    integer :: MCRD    !! IN OUR ANALYSIS MCRD=2

    double precision :: COORDS, PtGauss
    dimension COORDS(MCRD, nnode_patch), PtGauss(MCRD)

    !! Output variables
    !! ----------------
    double precision :: R, dRdx, d2Rdx, d3Rdx, DetJac
    dimension R(nnode_patch)
    dimension dRdx(MCRD, nnode_patch)
    dimension d2Rdx(MCRD + 1, nnode_patch)
    dimension d3Rdx(2 * MCRD, nnode_patch)


    !! Local variables
    !! ---------------
    double precision :: FN, FM, dNdxi, dMdEta, dxdxi, Detdxdxi, dxidx, xi, dRdxi ,  &
        &           d2Ndxi, d2MdEta, d2Rdxi, d2xdxi, d2Radxi, M2mat, Inv_M2mat, DetM2mat,    &
        &           d3Ndxi, d3MdEta, d3xdxi, d3Rdxi, matB,  &
        &           matA, Inv_matA, DetmatA, &
        &           dXidtildexi, AJmat
    double precision :: d3Rbdxi

    dimension FN(Jpqr_patch(1) + 1)
    dimension dNdxi(Jpqr_patch(1) + 1)
    dimension d2Ndxi(Jpqr_patch(1) + 1)
    dimension d3Ndxi(Jpqr_patch(1) + 1)

    dimension FM(Jpqr_patch(2) + 1)
    dimension dMdEta(Jpqr_patch(2) + 1)
    dimension d2MdEta(Jpqr_patch(2) + 1)
    dimension d3MdEta(Jpqr_patch(2) + 1)

    dimension dxdxi(MCRD, MCRD), dxidx(MCRD, MCRD)
    dimension d2xdxi(MCRD, MCRD + 1)
    dimension d3xdxi(MCRD, 2 * MCRD)
    dimension xi(dim_patch)

    dimension dRdxi(nnode_patch, MCRD)
    dimension d2Rdxi(nnode_patch, MCRD + 1)
    dimension d3Rdxi(nnode_patch, 2 * MCRD)


    dimension d2Radxi(MCRD + 1, nnode_patch)

    dimension M2mat(MCRD + 1, MCRD + 1), Inv_M2mat(MCRD + 1, MCRD + 1),&
        &     matB(2 * MCRD, MCRD + 1)

    dimension d3Rbdxi(nnode_patch, 2 * MCRD)

    dimension matA(2 * MCRD, 2 * MCRD), Inv_matA(2 * MCRD, 2 * MCRD), &
        &     AJmat(MCRD, MCRD)

    dimension dXidtildexi(MCRD, MCRD)

    Integer ::  i, j, k, Ni, NumLoc, Na, Nb, Nc
    dimension Ni(dim_patch)

    !! Knot spans of current element
    do i = 1, dim_patch
        Ni(i) = Nijk_patch(i, current_elem)
    end do



    !! **************************** INITIALIZATION ***************************************

    d2Rdx(:,:)       = zero

    d3Rdx(:,:)       = zero

    dXidtildexi(:,:) = zero
    AJmat(:,:)       = zero

    !! ******************* PARAMETRIC COORDINATES ***************************************

    !! Compute parametric coordinates of the given point    !! Ukv_elem(2,i)
    do i = 1, dim_patch
        xi(i) = ((Ukv_elem(2, i) - Ukv_elem(1, i)) * PtGauss(i)  &
            &  +  (Ukv_elem(2, i) + Ukv_elem(1, i)) ) * 0.5d0
    end do

    !! *******************  CALL DERIVATIVES *******************************************

    !! Compute univariate B-Spline function derivative

    call dersbasisfuns3(Ni(1), Jpqr_patch(1), Nkv_patch(1), xi(1),  &
        &   Ukv1_patch(:), FN, dNdxi, d2Ndxi, d3Ndxi)

    call dersbasisfuns3(Ni(2), Jpqr_patch(2), Nkv_patch(2), xi(2),  &
        &   Ukv2_patch(:), FM, dMdEta, d2MdEta, d3MdEta)


    !! *******************BUILD DERIVATIVES *******************************************
    NumLoc   = 0

    do j = 0, Jpqr_patch(2)
        do i = 0, Jpqr_patch(1)

            NumLoc = NumLoc + 1

            R(NumLoc) = FN(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)

            !! first derivatives
            dRdxi(NumLoc, 1) = dNdxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)
            dRdxi(NumLoc, 2) = FN(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j)

            !! second derivatives
            d2Rdxi(NumLoc, 1) = d2Ndxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)
            d2Rdxi(NumLoc, 2) = FN(Jpqr_patch(1) + 1 - i) * d2MdEta(Jpqr_patch(2) + 1 - j)
            d2Rdxi(NumLoc, 3) = dNdxi(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j)

            !! third derivatives
            d3Rdxi(NumLoc, 1) = d3Ndxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)
            d3Rdxi(NumLoc, 2) = FN(Jpqr_patch(1) + 1 - i) * d3MdEta(Jpqr_patch(2) + 1 - j)
            d3Rdxi(NumLoc, 3) = d2Ndxi(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j)
            d3Rdxi(NumLoc, 4) = dNdxi(Jpqr_patch(1) + 1 - i) * d2MdEta(Jpqr_patch(2) + 1 - j)
        end do
    end do

    ! write(*,*) "R : ", R
    ! write(*,*) "dRdxi : ", dRdxi
    ! write(*,*) "d2Rdxi : "!!, d2Rdxi
    ! do i = 1, nnode_patch
    ! write(*,*) d2Rdxi(i, :)
    ! enddo
    ! write(*,*) "d3Rdxi : ", d3Rdxi



    !! *************************** FIRST DERIVATIVE (dRdx)*********************************
    !! Gradient of mapping from parameter space to physical space
    dxdxi(:,:)       = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dxdxi(Na, Nb) = dxdxi(Na, Nb) + COORDS(Na, NumLoc) * dRdxi(NumLoc, Nb)
            end do
        end do
    end do


    !! Compute inverse of gradient
    call MatrixInv(dxidx, dxdxi, Detdxdxi, MCRD)

    !! Compute first derivatives of basis functions with respect to physical coordinates
    dRdx(:,:)        = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dRdx(Na, NumLoc) = dRdx(Na, NumLoc) + dRdxi(NumLoc, Nb) * dxidx(Nb, Na)
            end do
        end do
    end do

    !! Gradient of mapping from parent element to parameter space
    do i = 1, dim_patch
        dXidtildexi(i, i) = 0.5d0 * ( Ukv_elem(2, i) - Ukv_elem(1, i) )
    end do

    do Na = 1, MCRD
        do Nb = 1, MCRD
            do Nc = 1, MCRD
                AJmat(Na, Nb) = AJmat(Na, Nb) + dxdxi(Na, Nc) * dXidtildexi(Nc, Nb)
            end do
        end do
    end do

    DetJac = AJmat(1, 1) * AJmat(2, 2) - AJmat(2, 1) * AJmat(1, 2)


    !! ***************************SECOND DERIVATIVE (d2Rdx)*************************************
    !! Hessian of mapping from parameter space to physical space
    d2xdxi(:,:)      = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD + 1
                d2xdxi(Na, Nb) = d2xdxi(Na, Nb) + COORDS(Na, NumLoc) * d2Rdxi(NumLoc, Nb)
            end do
        end do
    end do

    !! Intermediate vector
    d2Radxi(:,:) = zero

    do NumLoc = 1, nnode_patch
        do  Na = 1, MCRD
            do Nb = 1, MCRD + 1
                d2Radxi(Nb, NumLoc) = d2Rdxi(NumLoc, Nb) - d2xdxi(Na, Nb) * dRdx(Na, NumLoc)
            end do
        end do
    end do

    !!	Build  matrice which will be inverted

    M2mat(1, 1) = dxdxi(1, 1) ** 2.
    M2mat(1, 2) = dxdxi(2, 1) ** 2.
    M2mat(1, 3) = two * dxdxi(1, 1) * dxdxi(2, 1)
    M2mat(2, 1) = dxdxi(1, 2) ** 2.
    M2mat(2, 2) = dxdxi(2, 2) ** 2.
    M2mat(2, 3) = two * dxdxi(1, 2) * dxdxi(2, 2)
    M2mat(3, 1) = dxdxi(1, 1) * dxdxi(1, 2)
    M2mat(3, 2) = dxdxi(2, 1) * dxdxi(2, 2)
    M2mat(3, 3) = (dxdxi(1, 1) * dxdxi(2, 2)) + (dxdxi(1, 2) * dxdxi(2, 1))


    call MatrixInv(Inv_M2mat, M2mat, DetM2mat, 3)

    !! Compute second derivatives of basis functions with respect to physical coordinates

    do NumLoc = 1, nnode_patch
        call MulVect(Inv_M2mat, d2Radxi(:, NumLoc), d2Rdx(:, NumLoc), MCRD + 1, MCRD + 1)
    end do

    !! ***************************THIRD DERIVATIVE (d3Rdx)*************************************


    !! Hessian derivative of mapping from parameter space to physical space

    !! Loop per col (4 rows et 2 col)
    d3xdxi(:,:) = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD !! 2 columns
            do Nb = 1, 2 * MCRD !! 4 rows
                d3xdxi(Na, Nb) = d3xdxi(Na, Nb) + COORDS(Na, NumLoc) * d3Rdxi(NumLoc, Nb)
            end do
        end do
    end do

    !! second intermediate derivative(4 lignes x 3 colonnes)

    matB(1, 1) = three * d2xdxi(1, 1) * dxdxi(1, 1)
    matB(1, 2) = three * d2xdxi(2, 1) * dxdxi(2, 1)
    matB(1, 3) = three * (d2xdxi(1, 1) * dxdxi(2, 1) + d2xdxi(2, 1) * dxdxi(1, 1))

    matB(2, 1) = three * d2xdxi(1, 2) * dxdxi(1, 2)
    matB(2, 2) = three * d2xdxi(2, 2) * dxdxi(2, 2)
    matB(2, 3) = three * (d2xdxi(1, 2) * dxdxi(2, 2) + d2xdxi(2, 2) * dxdxi(1, 2))

    matB(3, 1) = d2xdxi(1, 1) * dxdxi(1, 2) + two * d2xdxi(1, 3) * dxdxi(1, 1)
    matB(3, 2) = d2xdxi(2, 1) * dxdxi(2, 2) + two * d2xdxi(2, 3) * dxdxi(2, 1)
    matB(3, 3) = d2xdxi(1, 1) * dxdxi(2, 2) + d2xdxi(2, 1) * dxdxi(1, 2)  &
        & + two * (d2xdxi(2, 3) * dxdxi(1, 1) + d2xdxi(1, 3) * dxdxi(2, 1))

    matB(4, 1) = d2xdxi(1, 2) * dxdxi(1, 1) + two * d2xdxi(1, 3) * dxdxi(1, 2)
    matB(4, 2) = d2xdxi(2, 2) * dxdxi(2, 1) + two * d2xdxi(2, 3) * dxdxi(2, 2)
    matB(4, 3) = d2xdxi(2, 2) * dxdxi(1, 1) + d2xdxi(1, 2) * dxdxi(2, 1) &
        & + two * (d2xdxi(1, 3) * dxdxi(2, 2) + d2xdxi(2, 3) * dxdxi(1, 2))


    d3Rbdxi(:,:) = zero
    do NumLoc = 1, nnode_patch
        do  Na = 1, 2 * MCRD
            d3Rbdxi(NumLoc, Na) = d3Rdxi(NumLoc, Na)
            do Nb = 1, MCRD
                d3Rbdxi(NumLoc, Na) = d3Rbdxi(NumLoc, Na) - d3xdxi(Nb, Na) * dRdx(Nb, NumLoc)
            end do
            do Nb = 1, MCRD + 1
                d3Rbdxi(NumLoc, Na) = d3Rbdxi(NumLoc, Na) - matB(Na, Nb) * d2Rdx(Nb, NumLoc)
            end do
        end do
    end do

    !! Build  matrice which will be inverted

    matA(1, 1) = dxdxi(1, 1) ** three
    matA(1, 2) = dxdxi(2, 1) ** three
    matA(1, 3) = three * dxdxi(2, 1) * (dxdxi(1, 1)) ** two
    matA(1, 4) = three * dxdxi(1, 1) * (dxdxi(2, 1)) ** two


    matA(2, 1) = dxdxi(1, 2) ** three
    matA(2, 2) = dxdxi(2, 2) ** three
    matA(2, 3) = three * dxdxi(2, 2) * (dxdxi(1, 2)) ** two
    matA(2, 4) = three * dxdxi(1, 2) * (dxdxi(2, 2)) ** two


    matA(3, 1) = dxdxi(1, 2) * (dxdxi(1, 1)) ** two
    matA(3, 2) = dxdxi(2, 2) * (dxdxi(2, 1)) ** two
    matA(3, 3) = two * dxdxi(1, 1) * dxdxi(1, 2) * dxdxi(2, 1) &
        &   + dxdxi(2, 2) * (dxdxi(1, 1)) ** two
    matA(3, 4) = two * dxdxi(1, 1) * dxdxi(2, 2) * dxdxi(2, 1) &
        &   + dxdxi(1, 2) * dxdxi(2, 1) ** two

    matA(4, 1) = dxdxi(1, 1) * (dxdxi(1, 2)) ** two
    matA(4, 2) = dxdxi(2, 1) * (dxdxi(2, 2)) ** two
    matA(4, 3) = two * dxdxi(1, 1) * dxdxi(1, 2) * dxdxi(2, 2) &
        &   + dxdxi(2, 1) * (dxdxi(1, 2)) ** two
    matA(4, 4) = two * dxdxi(1, 2) * dxdxi(2, 1) * dxdxi(2, 2)   &
        &   + dxdxi(1, 1) * (dxdxi(2, 2)) ** two

    call MatrixInv(Inv_matA, matA, DetmatA, 2 * MCRD)

    do NumLoc = 1, nnode_patch
        call MulVect(Inv_matA, d3Rbdxi(NumLoc,:), d3Rdx(:, NumLoc), 2 * MCRD, 2 * MCRD)
    end do




end subroutine shap_HO

!! Same as shap_HO with only 1st gradient of strain
subroutine shap_HO_1stG(R, dRdx, d2Rdx, DetJac, COORDS, PtGauss, MCRD)

    use parameters
    use nurbspatch

    implicit none

    !! Input arguments
    !! ---------------
    integer :: MCRD    !! IN OUR ANALYSIS MCRD=2

    double precision :: COORDS, PtGauss
    dimension COORDS(MCRD, nnode_patch), PtGauss(MCRD)

    !! Output variables
    !! ----------------
    double precision :: R, dRdx, d2Rdx, d3Rdx, DetJac
    dimension R(nnode_patch)
    dimension dRdx(MCRD, nnode_patch)
    dimension d2Rdx(MCRD + 1, nnode_patch)

    !! Local variables
    !! ---------------
    double precision :: FN, FM, dNdxi, dMdEta, dxdxi, Detdxdxi, dxidx, xi, dRdxi ,  &
        &           d2Ndxi, d2MdEta, d2Rdxi, d2xdxi, d2Radxi, M2mat, Inv_M2mat, DetM2mat,    &
        &           matB,  &
        &           matA, Inv_matA, DetmatA, &
        &           dXidtildexi, AJmat
    double precision :: d3Rbdxi

    dimension FN(Jpqr_patch(1) + 1)
    dimension dNdxi(Jpqr_patch(1) + 1)
    dimension d2Ndxi(Jpqr_patch(1) + 1)

    dimension FM(Jpqr_patch(2) + 1)
    dimension dMdEta(Jpqr_patch(2) + 1)
    dimension d2MdEta(Jpqr_patch(2) + 1)

    dimension dxdxi(MCRD, MCRD), dxidx(MCRD, MCRD)
    dimension d2xdxi(MCRD, MCRD + 1)
    dimension xi(dim_patch)

    dimension dRdxi(nnode_patch, MCRD)
    dimension d2Rdxi(nnode_patch, MCRD + 1)

    dimension d2Radxi(nnode_patch, MCRD + 1)

    dimension M2mat(MCRD + 1, MCRD + 1), Inv_M2mat(MCRD + 1, MCRD + 1),&
        &     matB(2 * MCRD, MCRD + 1)

    dimension matA(2 * MCRD, 2 * MCRD), Inv_matA(2 * MCRD, 2 * MCRD), &
        &     AJmat(MCRD, MCRD)

    dimension dXidtildexi(MCRD, MCRD)

    Integer ::  i, j, k, Ni, NumLoc, Na, Nb, Nc
    dimension Ni(dim_patch)

    !! Knot spans of current element
    do i = 1, dim_patch
        Ni(i) = Nijk_patch(i, current_elem)
    end do



    !! **************************** INITIALIZATION ***************************************

    d2Rdx(:,:)       = zero

    dXidtildexi(:,:) = zero
    AJmat(:,:)       = zero

    !! ******************* PARAMETRIC COORDINATES ***************************************

    !! Compute parametric coordinates of the given point    !! Ukv_elem(2,i)
    do i = 1, dim_patch
        xi(i) = ((Ukv_elem(2, i) - Ukv_elem(1, i)) * PtGauss(i)  &
            &  +  (Ukv_elem(2, i) + Ukv_elem(1, i)) ) * 0.5d0
    end do

    !! *******************  CALL DERIVATIVES *******************************************

    !! Compute univariate B-Spline function derivative

    call dersbasisfuns2(Ni(1), Jpqr_patch(1), Nkv_patch(1), xi(1),  &
        &   Ukv1_patch(:), FN, dNdxi, d2Ndxi)

    call dersbasisfuns2(Ni(2), Jpqr_patch(2), Nkv_patch(2), xi(2),  &
        &   Ukv2_patch(:), FM, dMdEta, d2MdEta)


    !! *******************BUILD DERIVATIVES *******************************************
    NumLoc   = 0

    do j = 0, Jpqr_patch(2)
        do i = 0, Jpqr_patch(1)

            NumLoc = NumLoc + 1

            R(NumLoc) = FN(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)

            !! first derivatives
            dRdxi(NumLoc, 1) = dNdxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)
            dRdxi(NumLoc, 2) = FN(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j)

            !! second derivatives
            d2Rdxi(NumLoc, 1) = d2Ndxi(Jpqr_patch(1) + 1 - i) * FM(Jpqr_patch(2) + 1 - j)
            d2Rdxi(NumLoc, 2) = FN(Jpqr_patch(1) + 1 - i) * d2MdEta(Jpqr_patch(2) + 1 - j)
            d2Rdxi(NumLoc, 3) = dNdxi(Jpqr_patch(1) + 1 - i) * dMdEta(Jpqr_patch(2) + 1 - j)
        end do
    end do

    ! write(*,*) "R : ", R
    ! write(*,*) "dRdxi : ", dRdxi
    ! write(*,*) "d2Rdxi : "!!, d2Rdxi
    ! do i = 1, nnode_patch
    ! write(*,*) d2Rdxi(i, :)
    ! enddo

    !! *************************** FIRST DERIVATIVE (dRdx)*********************************
    !! Gradient of mapping from parameter space to physical space
    dxdxi(:,:)       = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dxdxi(Na, Nb) = dxdxi(Na, Nb) + COORDS(Na, NumLoc) * dRdxi(NumLoc, Nb)
            end do
        end do
    end do


    !! Compute inverse of gradient
    call MatrixInv(dxidx, dxdxi, Detdxdxi, MCRD)

    !! Compute first derivatives of basis functions with respect to physical coordinates
    dRdx(:,:)        = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD
                dRdx(Na, NumLoc) = dRdx(Na, NumLoc) + dRdxi(NumLoc, Nb) * dxidx(Nb, Na)
            end do
        end do
    end do

    !! Gradient of mapping from parent element to parameter space
    do i = 1, dim_patch
        dXidtildexi(i, i) = 0.5d0 * ( Ukv_elem(2, i) - Ukv_elem(1, i) )
    end do

    do Na = 1, MCRD
        do Nb = 1, MCRD
            do Nc = 1, MCRD
                AJmat(Na, Nb) = AJmat(Na, Nb) + dxdxi(Na, Nc) * dXidtildexi(Nc, Nb)
            end do
        end do
    end do

    DetJac = AJmat(1, 1) * AJmat(2, 2) - AJmat(2, 1) * AJmat(1, 2)


    !! ***************************SECOND DERIVATIVE (d2Rdx)*************************************
    !! Hessian of mapping from parameter space to physical space
    d2xdxi(:,:)      = zero

    do NumLoc = 1, nnode_patch
        do Na = 1, MCRD
            do Nb = 1, MCRD + 1
                d2xdxi(Na, Nb) = d2xdxi(Na, Nb) + COORDS(Na, NumLoc) * d2Rdxi(NumLoc, Nb)
            end do
        end do
    end do

    !! Intermediate vector
    d2Radxi(:,:) = zero

    do NumLoc = 1, nnode_patch
        do  Na = 1, MCRD
            do Nb = 1, MCRD + 1
                d2Radxi(NumLoc, Nb) = d2Rdxi(NumLoc, Nb) - d2xdxi(Na, Nb) * dRdx(Na, NumLoc)
            end do
        end do
    end do

    !!	Build  matrice which will be inverted

    M2mat(1, 1) = dxdxi(1, 1) ** 2.
    M2mat(1, 2) = dxdxi(2, 1) ** 2.
    M2mat(1, 3) = two * dxdxi(1, 1) * dxdxi(2, 1)
    M2mat(2, 1) = dxdxi(1, 2) ** 2.
    M2mat(2, 2) = dxdxi(2, 2) ** 2.
    M2mat(2, 3) = two * dxdxi(1, 2) * dxdxi(2, 2)
    M2mat(3, 1) = dxdxi(1, 1) * dxdxi(1, 2)
    M2mat(3, 2) = dxdxi(2, 1) * dxdxi(2, 2)
    M2mat(3, 3) = (dxdxi(1, 1) * dxdxi(2, 2)) + (dxdxi(1, 2) * dxdxi(2, 1))


    call MatrixInv(Inv_M2mat, M2mat, DetM2mat, 3)

    !! Compute second derivatives of basis functions with respect to physical coordinates

    do NumLoc = 1, nnode_patch
        call MulVect(Inv_M2mat, d2Radxi(NumLoc,:), d2Rdx(:, NumLoc), MCRD + 1, MCRD + 1)
    end do

end subroutine shap_HO_1stG
