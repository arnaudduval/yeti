!! Copyright 2021 Marie Guerder
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
!! with Yeti. if not, see <https://www.gnu.org/licenses/>


subroutine UELMAT10(NDOFEL, MCRD, NNODE, NNODEmap, nb_cp, JELEM,            &
            &   NBINT, COORDS, COORDSall, TENSOR, MATERIAL_PROPERTIES,      &
            &   DENSITY, nb_load, indDLoad, load_target_nbelem, JDLType,    &
            &   ADLMAG, load_additionalInfos, len_load_additionalInfos,     &
            &   nb_load_additionalInfos,                                    &
            &   n_dist_elem, nb_n_dist, RHS, AMATRX)

    use parameters
    use embeddedMapping

    implicit none

    !! Input arguments
    !! ---------------

    !! Embedded solid
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    integer, intent(in) :: nb_n_dist
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES,    &
        &   DENSITY, n_dist_elem
      dimension COORDS(MCRD,NNODE), MATERIAL_PROPERTIES(2),         &
        &   n_dist_elem(nb_n_dist,NNODE)

    !! Hull object
    integer, intent(in) :: NNODEmap, nb_cp
    double precision, intent(in) :: COORDSall
    dimension COORDSall(3,nb_cp)

    integer, intent(in) :: indDLoad, load_target_nbelem, JDLType,   &
        &   nb_load, len_load_additionalInfos,                      &
        &   nb_load_additionalInfos

    double precision, intent(in) :: ADLMAG, load_additionalInfos
    dimension ADLMAG(nb_load),      &
        &   load_target_nbelem(nb_load), JDLType(nb_load),  &
        &   load_additionalInfos(len_load_additionalInfos), &
        &   nb_load_additionalInfos(nb_load)
    dimension indDLoad(SUM(load_target_nbelem))

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD, MCRD, NNODE*(NNODE+1)/2)

    !! Local variables
    !! ---------------

    !! For gauss points
    integer :: NbPtInt
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(MCRD+1, NBINT)

    !! Embedded solid
    !! - NURBS basis functions
    double precision :: R, dRdTheta, dRdxi
    dimension R(NNODE), dRdTheta(NNODE, 3), dRdxi(3, NNODE)
    !! - Parametric space
    double precision :: Theta
    dimension Theta(3)
    !! - Physical space
    double precision :: xi
    dimension xi(3)
    !! - Mapping: parametric space >> physical space
    double precision :: dxidTheta, dThetadxi, detdxidTheta
    dimension dxidTheta(3, 3), dThetadxi(3, 3)
    !! - Mapping: parent element >> parametric space
    double precision :: dtildexidTheta, dThetadtildexi,     &
        &   detdThetadtildexi
    dimension dtildexidTheta(3, 3), dThetadtildexi(3, 3)

    !! Hull object
    !! - NURBS basis functions
    double precision :: N, dNdxi
    dimension N(NNODEmap), dNdxi(NNODEmap, 3)
    !! - Mapping: parametric space >> physical space
    double precision :: dXdxi, dxidX, detdXdxi, dThetadX
    dimension dXdxi(3, 3), dxidX(3, 3), dThetadX(3, 3)
    !! - Mapping: parent element >> parametric space
    double precision :: dxidtildexi, detdxidtildexi
    dimension dxidtildexi(3, 3)
    !! - Element infos
    double precision :: COORDSmap
    dimension COORDSmap(MCRD, NNODEmap)
    integer          :: sctr_map
    dimension sctr_map(NNODEmap)

    !! Composition Hull object x Embedded solid
    double precision :: dRdX, AJmat, dRdXT
    dimension dRdX(3, NNODE), AJmat(3, 3), dRdXT(NNODE, 3)

    !! For material matrix
    double precision :: ddsdde
    dimension ddsdde(2*MCRD, 2*MCRD)

    !! For stiffness matrices and load vectors
    integer :: iField
    double precision :: stiff, dvol, FbL, VectNorm, y, f_mag,   &
        &   coef, DetJac
    dimension stiff(MCRD, MCRD,NNODE*(NNODE+1)/2), FbL(NDOFEL), &
        &   VectNorm(MCRD)

    !! - centrifugal load
    integer :: load_addinfos_count
    double precision :: pointGP, pointA, pointB, vectD, vectAG, vectR,  &
        &   scal
    dimension pointGP(MCRD), pointA(MCRD), pointB(MCRD), vectD(MCRD),   &
        &   vectAG(MCRD), vectR(MCRD)

    !! For loops
    integer ntens, isave
    integer nint, k1, k2, i, j, kk, numCP, Numi, dof, KTypeDload,       &
        &   KNumFace, k3, Na, Nb, Nc, kload


    !! Initialization
    !! --------------
    ntens = 2*MCRD                          ! size of stiffness tensor
    NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
    if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1

    !! Defining Gauss points coordinates and Gauss weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)

    !! Stiffness matrix and force vector are initialized to zero
    RHS(:) = zero
    AMATRX(:, :, :) = zero

    !! Material behavior
    call material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde)


    !! Computation
    !! -----------

    isave = 0

    !! Loop on integration points
    do nint = 1, NBINT

        !! 1. Embedded solid
        !! ..................

        !! - Compute parametric coordinates from parent element
        Theta(:) = zero
        do j = 1, 3
            coef = GaussPdsCoord(1+j, nint)
            Theta(j) = ((Ukv_elem(2, j) - Ukv_elem(1, j))*coef      &
                &   + (Ukv_elem(2, j) + Ukv_elem(1, j)))*0.5d0
        enddo

        !! - Compute NURBS basis functions and derivatives of the embedded solid
        call evalnurbs(Theta, R, dRdTheta)

        !! - Compute embedded solid physical position
        !!   NB: physical space (embedded) = parametric space (hull)
        xi(:) = zero
        do numCP = 1, NNODE
            xi(:) =  xi(:) + R(numCP)*COORDS(:, numCP)
        enddo

        !! - Gradient of mapping: parent element >> parameter space (embedded)
        dThetadtildexi(:, :) = zero
        do j = 1, dim_patch
            dThetadtildexi(j, j) = 0.5d0*(Ukv_elem(2, j) -      &
                &   Ukv_elem(1, j))
        enddo

        call MatrixDet(dThetadtildexi, detdThetadtildexi, 3)

        !! Gradient of mapping: parameter space (embedded) >> physical space (embedded)
        !!   NB: physical space (embedded) = parametric space (hull)
        dxidTheta(:, :) = zero
        do numCP = 1, NNODE
            dxidTheta(:, 1) = dxidTheta(:, 1) +     &
                &   dRdTheta(numCP, 1)*COORDS(:, numCP)
            dxidTheta(:, 2) = dxidTheta(:, 2) +     &
                &   dRdTheta(numCP, 2)*COORDS(:, numCP)
            dxidTheta(:, 3) = dxidTheta(:, 3) +     &
                &   dRdTheta(numCP, 3)*COORDS(:, numCP)
        enddo

        call MatrixInv(dThetadxi, dxidTheta, detdxidTheta, 3)

        !! 2. Hull object
        !! ..............

        !! - Get active element number
        call updateMapElementNumber(xi(:))

        !! - Evaluate NURBS basis functions and derivatives of the hull object
        call evalnurbs_mapping(xi(:), N(:), dNdxi(:, :))

        !! - Extract coordinates of the CPs of the hull object
        if (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:, current_map_elem)
            do numCP = 1, NNODEmap
                COORDSmap(:, numCP) =       &
                    &   COORDSall(:, sctr_map(numCP))
            enddo
            isave = current_map_elem
        endif

        !! - Gradient of mapping: parameter space (hull) >> physical space (hull)
        dXdxi(:, :) = zero
        do numCP = 1, NNODEmap
            dXdxi(:, 1) = dXdxi(:, 1) +     &
                &   dNdxi(numCP, 1)*COORDSmap(:, numCP)
            dXdxi(:, 2) = dXdxi(:, 2) +     &
                &   dNdxi(numCP, 2)*COORDSmap(:, numCP)
            dXdxi(:, 3) = dXdxi(:, 3) +     &
                &   dNdxi(numCP, 3)*COORDSmap(:, numCP)
        enddo

        call MatrixInv(dxidX, dXdxi, detdXdxi, 3)

        !! 3. Composition: hull object x embedded solid
        !! ............................................

        !! - Intermediate mapping determinant product
        call MulMat(dThetadxi, dxidX, dThetadX, 3, 3, 3)

        !! - Basis functions composition
        call Mulmat(dRdTheta, dThetadX, dRdXT, NNODE, 3, 3)

        !! - Transpose basis functions array for further usage
        dRdX(:, :) = zero
        do numCP = 1, NNODE
            dRdX(1, numCP) = dRdXT(numCP, 1)
            dRdX(2, numCP) = dRdXT(numCP, 2)
            dRdX(3, numCP) = dRdXT(numCP, 3)
        enddo

        !! - Compute product of all mappings gradients
        DetJac = detdXdxi*detdxidTheta*detdThetadtildexi

        !! 4. Stiffness matrix & assembly
        !! ..............................

        !! - Compute stiffness matrix
        call stiffmatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde,   &
            &   dRdX, stiff)

        !! - Assemble AMATRIX
        dvol = GaussPdsCoord(1, nint)*DetJac

        AMATRX(:, :, :) = AMATRX(:, :, :) + stiff(:, :, :)*dvol

        !! 5. Body loads
        !! .............
        load_addinfos_count = 1
        kload = 0
        do i = 1, nb_load
            !! Centrifugal load
            !! $f_{b} = \rho \, \omega^{2} \, r$
            if ((JDLTYPE(i)==101) .and.             &
                &   ANY(indDLoad(kload+1:           &
                &   kload+load_target_nbelem(i))==JELEM)) then
                !! Gauss point location
                pointGP(:) = zero
                do numCP = 1, NNODEmap
                    pointGP(:) = pointGP(:) +           &
                        &   N(numCP)*COORDSmap(:, numCP)
                enddo
                !! Compute distance to rotation axis
                !! - Start point
                pointA(:) =                 &
                    &   load_additionalInfos(load_addinfos_count:       &
                    &                        load_addinfos_count+MCRD)
                !! - End point
                pointB(:) =                 &
                    &   load_additionalInfos(load_addinfos_count+MCRD:  &
                    &                        load_addinfos_count+2*MCRD)
                !! - Direction vector
                vectD(:) = pointB(:) - pointA(:)
                vectD(:) = vectD(:) / SQRT(SUM(vectD(:)*vectD(:)))  ! Normalise
                !! - Distance from Gauss point to start point
                vectAG(:) = pointGP(:) - pointA(:)
                !! - Scalar product
                call dot(vectAG(:), vectD(:), scal)
                !! - Final vector
                vectR(:) = vectAG(:) - scal*vectD(:)
                !! Update load vector
                kk = 0
                do numCP = 1, NNODE
                    do j = 1, MCRD
                        kk = kk + 1
                        RHS(kk) = RHS(kk)                                   &
                            &   + DENSITY * ADLMAG(i)**two * vectR(j) *     &
                            &   R(numCP) * dvol
                    enddo
                enddo
            endif
            kload = kload + load_target_nbelem(i)
            load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(i)
        enddo

    enddo  !! End of the loop on integration points

    !! Loop for load : find boundary loads
    load_addinfos_count = 1
    kk = 0
    do i = 1, nb_load
        if ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND.       &
            &   ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then

            !! Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
            call LectCle(JDLType(i), KNumFace, KTypeDload)
            if (KTypeDload == 4) then
                !! Get Index of nodal distribution
                iField = int(load_additionalInfos(load_addinfos_count))
            endif
            call Gauss(NbPtInt, MCRD, GaussPdsCoord, KNumFace)

            FbL(:) = zero

            do nint = 1, NbPtInt**(MCRD-1)
                !! 1. Embedded solid
                !! .................

                !! - Compute parametric coordinates from parent element
                Theta(:) = zero
                do j = 1, 3
                    coef = GaussPdsCoord(1+j, nint)
                    Theta(j) = ((Ukv_elem(2, j) - Ukv_elem(1, j))*coef      &
                        &   + (Ukv_elem(2, j) + Ukv_elem(1, j)))*0.5d0
                enddo

                !! - Compute NURBS basis functions and derivatives of the embedded solid
                call evalnurbs(Theta, R, dRdTheta)

                !! - Compute embedded solid physical position
                !!   NB: physical space (embedded) = parametric space (hull)
                xi(:) = zero
                do numCP = 1, NNODE
                    xi(:) =  xi(:) + R(numCP)*COORDS(:, numCP)
                enddo

                !! Gradient of mapping: parameter space (embedded) >> physical space (embedded)
                dThetadtildexi(:, :) = zero
                do j = 1, dim_patch
                    dThetadtildexi(j, j) = 0.5d0*(Ukv_elem(2, j) -          &
                        &   Ukv_elem(1, j))
                enddo

                !! Gradient of mapping: parameter space (embedded) >> physical space (embedded)
                !!   NB: physical space (embedded) = parametric space (hull)
                dxidTheta(:, :) = zero
                do numCP = 1, NNODE
                    dxidTheta(:, 1) = dxidTheta(:, 1) +         &
                        &   dRdTheta(numCP, 1)*COORDS(:, numCP)
                    dxidTheta(:, 2) = dxidTheta(:, 2) +         &
                        &   dRdTheta(numCP, 2)*COORDS(:, numCP)
                    dxidTheta(:, 3) = dxidTheta(:, 3) +         &
                        &   dRdTheta(numCP, 3)*COORDS(:, numCP)
                enddo

                !! 2. Hull object
                !! ..............

                !! - Get active element number
                call updateMapElementNumber(xi(:))

                !! - Evaluate NURBS basis functions and derivatives of the hull object
                call evalnurbs_mapping(xi(:), N(:), dNdxi(:, :))

                !! - Extract coordinates of the CPs of the hull object
                if (isave /= current_map_elem) then
                    sctr_map(:) = IEN_map(:, current_map_elem)
                    do numCP = 1, NNODEmap
                        COORDSmap(:, numCP) =       &
                            &   COORDSall(:, sctr_map(numCP))
                    enddo
                    isave = current_map_elem
                endif

                !! - Gradient of mapping: parameter space (hull) >> physical space (hull)
                dXdxi(:, :) = zero
                do numCP = 1, NNODEmap
                    dXdxi(:, 1) = dXdxi(:, 1) +     &
                        &   dNdxi(numCP, 1)*COORDSmap(:, numCP)
                    dXdxi(:, 2) = dXdxi(:, 2) +     &
                        &   dNdxi(numCP, 2)*COORDSmap(:, numCP)
                    dXdxi(:, 3) = dXdxi(:, 3) +     &
                        &   dNdxi(numCP, 3)*COORDSmap(:, numCP)
                enddo

                !! 3. Composition: hull object x embedded solid
                !! ............................................

                !! - Intermediate product
                call MulMat(dxidTheta, dThetadtildexi, dxidtildexi,     &
                    &   3, 3, 3)

                !! - Final gradient matrix
                call MulMat(dXdxi, dxidtildexi, AJmat, 3, 3, 3)

                !! 4. Process loads
                !! ................

                !! - Compute element surface according to the face number
                select case (KNumFace+MCRD*10)
                    case(31, 32)               ! Faces 1 and 2
                        call SurfElem(AJmat(:, 2), AJmat(:,3), DetJac)
                    case(33, 34)               ! Faces 3 and 4
                        call SurfElem(AJmat(:, 3), AJmat(:,1), DetJac)
                    case(35, 36)               ! Faces 5 and 6
                        call SurfElem(AJmat(:, 1), AJmat(:,2), DetJac)
                end select

                !! - Compute normalized normal vector
                select case (KTypeDload+MCRD*10)
                    case(31)                  ! Constraint x direction
                        VectNorm(1) = one
                        VectNorm(2) = zero
                        VectNorm(3) = zero
                    case(32)                  ! Constraint y direction
                        VectNorm(1) = zero
                        VectNorm(2) = one
                        VectNorm(3) = zero
                    case(33)                  ! Constraint z direction
                        VectNorm(1) = zero
                        VectNorm(2) = zero
                        VectNorm(3) = one
                    case(30)                  ! Normal pressure in 3D
                        select case (KNumFace)
                            case(1)                ! Face 1
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 2), AJmat(:, 3), DetJac)
                            case(2)                ! Face 2
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 3), AJmat(:, 2), DetJac)
                            case(3)                ! Face 3
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 3), AJmat(:, 1), DetJac)
                            case(4)                ! Face 4
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 1), AJmat(:, 3), DetJac)
                            case(5)                ! Face 5
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 1),AJmat(:, 2), DetJac)
                            case(6)                ! Face 6
                                  call VectNormNorm(VectNorm,       &
                                    &   AJmat(:, 2), AJmat(:, 1), DetJac)
                        end select
                    case(34)                  ! Non-uniform normal pressure
                        select case (KNumFace)
                            case(1)                ! Face 1
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 2), AJmat(:, 3), DetJac)
                            case(2)                ! Face 2
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 3),AJmat(:, 2), DetJac)
                            case(3)                ! Face 3
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 3), AJmat(:, 1), DetJac)
                            case(4)                ! Face 4
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 1), AJmat(:, 3), DetJac)
                            case(5)                ! Face 5
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 1), AJmat(:, 2), DetJac)
                            case(6)                ! Face 6
                                call VectNormNorm(VectNorm,         &
                                    &   AJmat(:, 2), AJmat(:, 1), DetJac)
                        end select

                    end select

                    !! - Compute solid gradient
                    dvol = GaussPdsCoord(1, nint)*DetJac

                    !! - Compute force magnitude
                    !! Non-uniform pressure case
                    if (KTypeDload==4) then
                        f_mag = 0
                        do k3 = 1, NNODE
                            f_mag = f_mag + n_dist_elem(iField,k3) * R(k3)
                        enddo
                    !! Uniform pressure case
                    else
                        f_mag = ADLMAG(i)
                    endif

                    !! - Assembling force vector
                    do numCP = 1, NNODE
                        numI = (numCP-1)*MCRD
                        do k2 = 1, MCRD
                            numI = numI + 1
                            FbL(numI) = FbL(numI) +     &
                                &   R(numCP)*VectNorm(k2)*dvol*f_mag
                        enddo
                    enddo
                enddo

            !! - Assembling RHS
            do k1 = 1,NDOFEL
                RHS(k1) = RHS(k1) + FbL(k1)
            enddo

        endif

        kk = kk + load_target_nbelem(i)
        load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(i)

    enddo
end subroutine UELMAT10
