subroutine uelmat_byCP_omp(patch, element, ndofel, mcrd, jelem, nbint, coords,            &
    &   material_properties, density, nb_load, indDLoad,    &
    &   load_target_nbelem,JDLType,ADLMAG,load_additionalInfos, &
    &   len_load_additionalInfos, nb_load_additionalInfos, &
    &   n_dist_elem, nb_n_dist, RHS,AMATRX)

    use m_nurbspatch, only : nurbspatch, nurbselement

    implicit none

    !! Input arguments
    !! ---------------

    type(nurbspatch), intent(in) :: patch
    type(nurbselement), intent(in) :: element

    integer, intent(in) :: NDOFEL,MCRD,JELEM,NBINT
    double precision, intent(in) :: COORDS,MATERIAL_PROPERTIES, &
        &   DENSITY,    &
        &   n_dist_elem
    dimension COORDS(MCRD,patch%nnode_patch),MATERIAL_PROPERTIES(2),    &
        &     n_dist_elem(nb_n_dist,patch%nnode_patch)

    integer, intent(in) :: indDLoad,load_target_nbelem,JDLType, &
        &     nb_load,len_load_additionalInfos, nb_n_dist,  &
        &     nb_load_additionalInfos

    double precision, intent(in) :: ADLMAG,load_additionalInfos
    dimension ADLMAG(nb_load),indDLoad(SUM(load_target_nbelem)),    &
        &     load_target_nbelem(nb_load),JDLType(nb_load),         &
        &     load_additionalInfos(len_load_additionalInfos),       &
        &     nb_load_additionalInfos(nb_load)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD,MCRD,patch%nnode_patch*(patch%nnode_patch+1)/2)


    !! Local variables
    !! ---------------

    !! Gauss points
    integer :: nb_pt_int        !! Number of integration points per direction
    double precision, dimension(mcrd+1, NBINT) :: gauss_pds_coords    !! coordinates and weights of integration points
    integer :: igauss           !! Loop counter

    !! NURS basis functions
    double precision, dimension(patch%nnode_patch) :: R
    double precision, dimension(mcrd, patch%nnode_patch) :: dRdx
    double precision :: det_jac

    !! Material behaviour
    double precision, dimension(2*mcrd, 2*mcrd) :: ddsdde               !! Tangent modulus


    !! Stiffness matrix
    integer :: ntens
    double precision, dimension(mcrd, mcrd, patch%nnode_patch*(patch%nnode_patch+1)/2) :: stiff
    double precision :: dvol

    !! Load vector
    integer :: iload, kload, icp, ifield, numI
    integer :: k_numface, k_typedload
    double precision, dimension(NDOFEL) :: FbL
    double precision, dimension(mcrd) :: norm_vect
    double precision :: f_mag

    !! Loads requiring additional informations
    !! (centrifugal body force, distributed pressure, ...)
    integer :: load_addinfos_count      !! index for number of additional infos
    double precision, dimension(mcrd) :: pointGP, vectAG, vectR
    double precision, dimension(mcrd) :: pointA, pointB, vectD     !! Centrifugal force axis
    double precision :: scal

    !! counters
    integer :: kk, j, k3, k2, k1

    !! Initialization
    ntens = 2*mcrd      !! Stiffness tensor size
    nb_pt_int = int(NBINT**(1.0/float(mcrd)))       !! Number of Gauss point per direction
    !! TODO improve code to handle different degrees in different directions
    if (nb_pt_int**mcrd < NBINT) nb_pt_int = nb_pt_int + 1

    !! Compute Gauss points coordinates anb weights
    call Gauss(nb_pt_int, mcrd, gauss_pds_coords, 0)

    !! initialize matrix and load vector to zero
    RHS(:) = 0.0D0
    AMATRX(:,:,:) = 0.0d0

    !! Compute material beheviour
    !! TODO this must be moved inside Gauss points loop to handle plasticity
    call material_lib(MATERIAL_PROPERTIES, patch%tensor_patch, mcrd, ddsdde)

    !! Loop on integration points
    do igauss = 1, NBINT
        !! Compute NURBS basis functions and derivatives
        call shap_omp(patch, element, dRdx, R, det_jac, COORDS, gauss_pds_coords(2:, igauss), mcrd)

        !! Compute stiffness matrix
        call stiffmatrix_byCP(ntens, patch%nnode_patch, mcrd, NDOFEL, ddsdde,       &
                        &   dRdx, stiff)

        !! Assemble AMATRX
        dvol = gauss_pds_coords(1, igauss)*abs(det_jac)
        AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:)*dvol

        !! body load
        load_addinfos_count = 1
        kload = 0
        do iload = 1, nb_load
            if (JDLType(iload)==101) then
                if (any(indDLoad(kload+1:kload+load_target_nbelem(iload))==JELEM)) then
                    !! centrifugal force
                    ! Gauss point location in physical space
                    pointGP(:) = 0.0
                    do icp = 1, patch%nnode_patch
                        pointGP(:) = pointGP(:) + R(icp)*COORDS(:, icp)
                    enddo
                    ! Distance to rotation axis
                    pointA(:) = load_additionalInfos(load_addinfos_count:    &
                                &                    load_addinfos_count+MCRD)
                    pointB(:) = load_additionalInfos(load_addinfos_count+MCRD:   &
                                &                    load_addinfos_count+2*MCRD)
                    vectD(:)  = pointB(:) - pointA(:)
                    vectD(:)  = vectD(:)/SQRT(SUM(vectD(:)*vectD(:)))
                    vectAG(:) = pointGP(:) - pointA(:)
                    ! call dot(vectAG(:),vectD(:),scal)
                    scal = dot_product(vectAG(:), vectD(:))
                    vectR(:)   = vectAG(:) - scal*vectD(:)
                    !! Update load vector
                    kk = 0
                    do icp = 1, patch%nnode_patch
                        do j = 1,MCRD
                            kk = kk+1
                            RHS(kk) = RHS(kk) + DENSITY*ADLMAG(iload)**2.0*vectR(j)*R(icp)*dvol
                        enddo
                    enddo
                endif
            endif
            kload = kload + load_target_nbelem(iload)
            load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(iload)
        enddo
    enddo       !! End loop on integration points

    !! Loop for boundary load
    load_addinfos_count = 1
    kk = 0
    do iload = 1, nb_load
        if ((JDLTYPE(iload)>9 .AND. JDLTYPE(iload)<100) .AND.   &
                &   ANY(indDLoad(kk+1:kk+load_target_nbelem(iload))==JELEM)) then
            !! Define Gauss points coordinates and weights on surf(3D)/edge(2D)
            call LectCle (JDLType(iload),k_numface,k_typedload)

            if (k_typedload == 4) then
                !! Get Index of nodal distribution
                ifield = int(load_additionalInfos(load_addinfos_count))
            endif
            call Gauss (nb_pt_int,MCRD,gauss_pds_coords,k_numface)

            FbL(:) = 0.0
            do igauss = 1, nb_pt_int**(MCRD-1)

                call shapPress_omp(patch, element, R,norm_vect,det_jac,COORDS,        &
                    &   gauss_pds_coords(2:,igauss),MCRD,k_numface,k_typedload)

                dvol = gauss_pds_coords(1, igauss)*det_jac

                !! Non-uniform pressure case
                if (k_typedload==4) then
                    f_mag = 0
                    do k3 = 1, patch%nnode_patch
                        f_mag = f_mag + n_dist_elem(iField,k3) * R(k3)
                    enddo
                !! Uniform pressure case
                else
                    f_mag = ADLMAG(iload)
                endif

                do iCP = 1, patch%nnode_patch
                    numI = (iCP-1)*MCRD
                    do k2 = 1,MCRD
                        numI = numI + 1
                        FbL(numI)=FbL(numI) +R(iCP)*norm_vect(k2)*dvol*f_mag
                    enddo
                enddo
            enddo

            !! Assemble RHS
            do k1 = 1,NDOFEL
                RHS(k1) = RHS(k1) + FbL(k1)
            enddo
        endif
        kk = kk + load_target_nbelem(iload)
        load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(iload)
    enddo


end subroutine uelmat_byCP_omp