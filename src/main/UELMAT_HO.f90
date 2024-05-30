!! Compute elementary matrix and RHS vector
!! for solid 2D/3D elements with a higher order elastic material behaviour


subroutine UELMAT_HO_byCP(NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS,            &
        &   TENSOR, MATERIAL_PROPERTIES, n_mat_props, DENSITY, nb_load, indDLoad,    &
        &   load_target_nbelem, JDLType, ADLMAG, load_additionalInfos, &
        &   nb_load_additionalInfos, n_dist_elem, nb_n_dist, RHS, AMATRX)


    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, &
        &   DENSITY,    &
        &   n_dist_elem
    integer, intent(in) :: n_mat_props
    dimension COORDS(MCRD, NNODE), MATERIAL_PROPERTIES(n_mat_props),    &
        &     n_dist_elem(nb_n_dist, NNODE)

    integer, intent(in) :: indDLoad, load_target_nbelem, JDLType, &
        &     nb_load, nb_load_additionalInfos, nb_n_dist

    double precision, intent(in) :: ADLMAG, load_additionalInfos
    dimension ADLMAG(nb_load), indDLoad(SUM(load_target_nbelem)),    &
        &     load_target_nbelem(nb_load), JDLType(nb_load),         &
        &     load_additionalInfos(nb_load_additionalInfos)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD, MCRD, NNODE * (NNODE + 1) / 2)

    !! Local variables
    !! ---------------

    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(MCRD + 1, NBINT)

    !! Nurbs basis functions
    double precision :: R, dRdx, d2Rdx, d3Rdx, DetJac
    dimension R(NNODE), dRdx(MCRD, NNODE), d2Rdx(MCRD + 1, NNODE), d3Rdx(2 * MCRD, NNODE)

    !! Material behaviour
    double precision :: ddsdde
    dimension ddsdde(8 * MCRD + 1, 8 * MCRD + 1)

    !! Stiffness matrix
    integer :: k1, k2, ntens
    double precision :: stiff, dvol
    dimension stiff( MCRD, MCRD, NNODE * (NNODE + 1) / 2 )

    !! Load vector
    integer :: i, j, kk, KNumFace, KTypeDload, numCP, numI, k3, iField
    double precision :: FbL, VectNorm, y, f_mag
    dimension FbL(NDOFEL), VectNorm(MCRD)
    !! centrifugal load
    integer :: loadcount
    double precision :: pointGP, pointA, pointB, vectD, vectAG, vectR, scal
    dimension pointGP(MCRD), pointA(MCRD), pointB(MCRD), vectD(MCRD),  &
        &     vectAG(MCRD), vectR(MCRD)

    !! Initialization

    ntens   = 8 * MCRD + 1          !! Size of stiffness tensor
    NbPtInt = int( NBINT ** (1.0 / float(MCRD)) ) !! Nb of gauss pts per direction
    if (NbPtInt ** MCRD < NBINT) NbPtInt = NbPtInt + 1

    !! Compute Gauss points coordinates and weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)

    !! Stiffness matrix and load vector initialized to zero
    RHS(:)        = zero
    AMATRX(:,:,:) = zero

    !! Material behaviour

    call material_lib_HO(MATERIAL_PROPERTIES(:n_mat_props), n_mat_props, TENSOR, MCRD, ddsdde)
    ! write(*,*) "ddsdde : ", ddsdde

    !! Loop on integration points
    do n = 1, NBINT
        ! write(*,*) "=== IP : ", n, "==="
        !! Compute NURBS basis functions and derivatives
        call shap_HO(R, dRdx, d2Rdx, d3Rdx, DetJac, COORDS, GaussPdsCoord(2:, n), MCRD)
        ! write(*,*) "R", R
        ! write(*,*) "dRdx", dRdx
        ! write(*,*) "d2Rx", d2rdx
        ! write(*,*) "d3Rdx", d3Rdx
        !! Compute stiffness matrix
        call stiffmatrix_HO_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, d2Rdx, d3Rdx, stiff)
        ! write(*,*) "stiff : "
        ! write(*,*) stiff
        !! Assemble AMATRIX
        dvol = GaussPdsCoord(1, n) * detJac
        AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:) * dvol
        !! body load
        loadcount = 1
        do i = 1, nb_load
            if (JDLTYPE(i) == 101) then
                !! Centrifugal load
                !! Gauss point location
                pointGP(:) = zero
                do numCP = 1, NNODE
                    pointGP(:) = pointGP(:) + R(numCP) * COORDS(:, numCP)
                end do
                !! Distance to rotation axis
                pointA(:) = load_additionalInfos(loadcount:loadcount + MCRD)
                loadcount = loadcount + MCRD
                pointB(:) = load_additionalInfos(loadcount:loadcount + MCRD)
                loadcount = loadcount + MCRD

                vectD(:)  = pointB(:) - pointA(:)
                vectD(:)  = vectD(:) / SQRT(SUM(vectD(:) * vectD(:)))
                vectAG(:) = pointGP(:) - pointA(:)
                call dot(vectAG(:), vectD(:), scal)
                vectR(:)   = vectAG(:) - scal * vectD(:)
                !! Update load vector
                kk = 0
                do numCP = 1, NNODE
                    do j = 1, MCRD
                        kk = kk + 1
                        RHS(kk) = RHS(kk) + DENSITY * ADLMAG(i) ** two * vectR(j) * R(numCP) * dvol
                    end do
                end do
            end if
        end do
    end do   !! End of the loop on integration points
    ! write(*,*) "K"
    ! do i = 1, NNODE*(NNODE+1)/2
    ! write(*,*) i, AMATRX(:,:,i)
    ! enddo

    !! Loop for load : find boundary loads
    kk = 0
    do i = 1, nb_load
        if ((JDLTYPE(i) > 9 .AND. JDLTYPE(i) < 100) .AND.   &
            &   ANY(indDLoad(kk + 1:kk + load_target_nbelem(i)) == JELEM)) then
        !! Define Gauss points coordinates and weights on surf(3D)/edge(2D)
        call LectCle (JDLType(i), KNumFace, KTypeDload)
        if (KTypeDload == 4) then
            !! Get Idex of nodal distribution
            iField = int(load_additionalInfos(loadcount))
            loadcount = loadcount + 1
        end if
        call Gauss (NbPtInt, MCRD, GaussPdsCoord, KNumFace)

        FbL(:) = zero
        do n = 1, NbPtInt ** (MCRD - 1)

            call shapPress(R, VectNorm, DetJac, COORDS,        &
                &   GaussPdsCoord(2:, n), MCRD, KNumFace, KTypeDload,   &
                &   VectNorm)

            dvol = GaussPdsCoord(1, n) * DetJac

            !! Non-uniform pressure case
            if (KTypeDload == 4) then
                f_mag = 0
                do k3 = 1, NNODE
                    f_mag = f_mag + n_dist_elem(iField, k3) * R(k3)
                end do
                !! Uniform pressure case
            else
                f_mag = ADLMAG(i)
            end if

            do numCP = 1, NNODE
                numI = (numCP - 1) * MCRD
                do k2 = 1, MCRD
                    numI = numI + 1
                    FbL(numI) = FbL(numI) + R(numCP) * VectNorm(k2) * dvol * f_mag
                end do
            end do
        end do

        !! Assemble RHS
        do k1 = 1, NDOFEL
            RHS(k1) = RHS(k1) + FbL(k1)
        end do
    end if
    kk = kk + load_target_nbelem(i)
end do
end subroutine UELMAT_HO_byCP


!!!!!!!!!!!!!!!!!!!!!!!
!! Same ad UELMAT_HO_byCP with only 1st gradient of strain taken into account
subroutine UELMAT_HO_byCP_1stG(NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS,            &
        &   TENSOR, MATERIAL_PROPERTIES, n_mat_props, DENSITY, nb_load, indDLoad,    &
        &   load_target_nbelem, JDLType, ADLMAG, load_additionalInfos, &
        &   nb_load_additionalInfos, n_dist_elem, nb_n_dist, RHS, AMATRX)


    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, &
        &   DENSITY,    &
        &   n_dist_elem
    integer, intent(in) :: n_mat_props
    dimension COORDS(MCRD, NNODE), MATERIAL_PROPERTIES(n_mat_props),    &
        &     n_dist_elem(nb_n_dist, NNODE)

    integer, intent(in) :: indDLoad, load_target_nbelem, JDLType, &
        &     nb_load, nb_load_additionalInfos, nb_n_dist

    double precision, intent(in) :: ADLMAG, load_additionalInfos
    dimension ADLMAG(nb_load), indDLoad(SUM(load_target_nbelem)),    &
        &     load_target_nbelem(nb_load), JDLType(nb_load),         &
        &     load_additionalInfos(nb_load_additionalInfos)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD, MCRD, NNODE * (NNODE + 1) / 2)

    !! Local variables
    !! ---------------

    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(MCRD + 1, NBINT)

    !! Nurbs basis functions
    double precision :: R, dRdx, d2Rdx, DetJac
    dimension R(NNODE), dRdx(MCRD, NNODE), d2Rdx(MCRD + 1, NNODE)

    !! Material behaviour
    double precision :: ddsdde
    dimension ddsdde(9, 9)

    !! Stiffness matrix
    integer :: k1, k2, ntens
    double precision :: stiff, dvol
    dimension stiff( MCRD, MCRD, NNODE * (NNODE + 1) / 2 )

    !! Load vector
    integer :: i, j, kk, KNumFace, KTypeDload, numCP, numI, k3, iField
    double precision :: FbL, VectNorm, y, f_mag
    dimension FbL(NDOFEL), VectNorm(MCRD)
    !! centrifugal load
    integer :: loadcount
    double precision :: pointGP, pointA, pointB, vectD, vectAG, vectR, scal
    dimension pointGP(MCRD), pointA(MCRD), pointB(MCRD), vectD(MCRD),  &
        &     vectAG(MCRD), vectR(MCRD)

    !! Initialization

    ntens   = 9          !! Size of stiffness tensor
    NbPtInt = int( NBINT ** (1.0 / float(MCRD)) ) !! Nb of gauss pts per direction
    if (NbPtInt ** MCRD < NBINT) NbPtInt = NbPtInt + 1

    !! Compute Gauss points coordinates and weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)

    !! Stiffness matrix and load vector initialized to zero
    RHS(:)        = zero
    AMATRX(:,:,:) = zero

    !! Material behaviour

    call material_lib_HO_1stG(MATERIAL_PROPERTIES(:n_mat_props), n_mat_props, TENSOR, MCRD, ddsdde)
    ! write(*,*) "ddsdde : ", ddsdde

    !! Loop on integration points
    do n = 1, NBINT
        ! write(*,*) "=== IP : ", n, "==="
        !! Compute NURBS basis functions and derivatives
        call shap_HO_1stG(R, dRdx, d2Rdx, DetJac, COORDS, GaussPdsCoord(2:, n), MCRD)
        ! write(*,*) "R", R
        ! write(*,*) "dRdx", dRdx
        ! write(*,*) "d2Rx", d2rdx
        ! write(*,*) "d3Rdx", d3Rdx
        !! Compute stiffness matrix
        call stiffmatrix_HO_byCP_1stG(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, d2Rdx, stiff)
        ! write(*,*) "stiff : "
        ! write(*,*) stiff
        !! Assemble AMATRIX
        dvol = GaussPdsCoord(1, n) * detJac
        AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:) * dvol
        !! body load
        loadcount = 1
        do i = 1, nb_load
            if (JDLTYPE(i) == 101) then
                !! Centrifugal load
                !! Gauss point location
                pointGP(:) = zero
                do numCP = 1, NNODE
                    pointGP(:) = pointGP(:) + R(numCP) * COORDS(:, numCP)
                end do
                !! Distance to rotation axis
                pointA(:) = load_additionalInfos(loadcount:loadcount + MCRD)
                loadcount = loadcount + MCRD
                pointB(:) = load_additionalInfos(loadcount:loadcount + MCRD)
                loadcount = loadcount + MCRD

                vectD(:)  = pointB(:) - pointA(:)
                vectD(:)  = vectD(:) / SQRT(SUM(vectD(:) * vectD(:)))
                vectAG(:) = pointGP(:) - pointA(:)
                call dot(vectAG(:), vectD(:), scal)
                vectR(:)   = vectAG(:) - scal * vectD(:)
                !! Update load vector
                kk = 0
                do numCP = 1, NNODE
                    do j = 1, MCRD
                        kk = kk + 1
                        RHS(kk) = RHS(kk) + DENSITY * ADLMAG(i) ** two * vectR(j) * R(numCP) * dvol
                    end do
                end do
            end if
        end do
    end do   !! End of the loop on integration points
    ! write(*,*) "K"
    ! do i = 1, NNODE*(NNODE+1)/2
    ! write(*,*) i, AMATRX(:,:,i)
    ! enddo

    !! Loop for load : find boundary loads
    kk = 0
    do i = 1, nb_load
        if ((JDLTYPE(i) > 9 .AND. JDLTYPE(i) < 100) .AND.   &
            &   ANY(indDLoad(kk + 1:kk + load_target_nbelem(i)) == JELEM)) then
        !! Define Gauss points coordinates and weights on surf(3D)/edge(2D)
        call LectCle (JDLType(i), KNumFace, KTypeDload)
        if (KTypeDload == 4) then
            !! Get Idex of nodal distribution
            iField = int(load_additionalInfos(loadcount))
            loadcount = loadcount + 1
        end if
        call Gauss (NbPtInt, MCRD, GaussPdsCoord, KNumFace)

        FbL(:) = zero
        do n = 1, NbPtInt ** (MCRD - 1)

            call shapPress(R, VectNorm, DetJac, COORDS,        &
                &   GaussPdsCoord(2:, n), MCRD, KNumFace, KTypeDload,   &
                &   VectNorm)

            dvol = GaussPdsCoord(1, n) * DetJac

            !! Non-uniform pressure case
            if (KTypeDload == 4) then
                f_mag = 0
                do k3 = 1, NNODE
                    f_mag = f_mag + n_dist_elem(iField, k3) * R(k3)
                end do
                !! Uniform pressure case
            else
                f_mag = ADLMAG(i)
            end if

            do numCP = 1, NNODE
                numI = (numCP - 1) * MCRD
                do k2 = 1, MCRD
                    numI = numI + 1
                    FbL(numI) = FbL(numI) + R(numCP) * VectNorm(k2) * dvol * f_mag
                end do
            end do
        end do

        !! Assemble RHS
        do k1 = 1, NDOFEL
            RHS(k1) = RHS(k1) + FbL(k1)
        end do
    end if
    kk = kk + load_target_nbelem(i)
end do
end subroutine UELMAT_HO_byCP_1stG
