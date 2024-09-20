!! Copyright 2011 Florian Maurin
!! Copyright 2016-2020 Thibaut Hirschler
!! Copyright 2019-2023 Arnaud Duval

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

!! Compute elementary matrix and RHS vector
!! for solid 2D/3D elements


subroutine UELVEC_byCP(NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS,            &
        &   TENSOR, DENSITY, nb_load, indDLoad,    &
        &   load_target_nbelem, JDLType, ADLMAG, load_additionalInfos, &
        &   len_load_additionalInfos, nb_load_additionalInfos, &
        &   n_dist_elem, nb_n_dist, RHS)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, &
        &   DENSITY,    &
        &   n_dist_elem
    dimension COORDS(MCRD, NNODE),    &
        &     n_dist_elem(nb_n_dist, NNODE)

    integer, intent(in) :: indDLoad, load_target_nbelem, JDLType, &
        &     nb_load, len_load_additionalInfos, nb_n_dist,  &
        &     nb_load_additionalInfos

    double precision, intent(in) :: ADLMAG, load_additionalInfos
    dimension ADLMAG(nb_load), indDLoad(SUM(load_target_nbelem)),    &
        &     load_target_nbelem(nb_load), JDLType(nb_load),         &
        &     load_additionalInfos(len_load_additionalInfos),       &
        &     nb_load_additionalInfos(nb_load)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS
    dimension RHS(NDOFEL)


    !! Local variables
    !! ---------------

    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(MCRD + 1, NBINT)

    !! Nurbs basis functions
    double precision :: R, dRdx, DetJac
    dimension R(NNODE), dRdx(MCRD, NNODE)

    !! Material behaviour
    double precision :: ddsdde
    dimension ddsdde(2 * MCRD, 2 * MCRD)

    !! Stiffness matrix
    integer :: k1, k2, ntens
    double precision :: stiff, dvol
    dimension stiff( MCRD, MCRD, NNODE * (NNODE + 1) / 2 )

    !! Load vector
    integer :: i, j, kk, KNumFace, KTypeDload, numCP, numI, k3, iField
    integer :: kload, i_load
    double precision :: FbL, VectNorm, y, f_mag
    dimension FbL(NDOFEL), VectNorm(MCRD)

    !! loads requiring additional information
    !! (centrufugal body force, distributed pressure, ...)
    integer :: load_addinfos_count  ! Index for number of additional infos
    double precision :: pointGP, pointA, pointB, vectD, vectAG, vectR, scal
    dimension pointGP(MCRD), pointA(MCRD), pointB(MCRD), vectD(MCRD),  &
        &     vectAG(MCRD), vectR(MCRD)

    !! Initialization

    ntens   = 2 * MCRD          !! Size of stiffness tensor
    NbPtInt = int( NBINT ** (1.0 / float(MCRD)) ) !! Nb of gauss pts per direction
    if (NbPtInt ** MCRD < NBINT) NbPtInt = NbPtInt + 1

    !! Compute Gauss points coordinates and weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)

    !! Stiffness matrix and load vector initialized to zero
    RHS(:)        = zero

    !! Loop on integration points
    do n = 1, NBINT
        !! Compute NURBS basis functions and derivatives
        call shap(dRdx, R, DetJac, COORDS, GaussPdsCoord(2:, n), MCRD)

        !! body load
        load_addinfos_count = 1
        kload = 0
        do i_load = 1, nb_load
            if (JDLTYPE(i_load) == 101) then
                if (ANY(indDLoad(kload + 1:kload + load_target_nbelem(i_load)) == JELEM)) then
                    !! Centrifugal load
                    !! Gauss point location
                    pointGP(:) = zero
                    do numCP = 1, NNODE
                        pointGP(:) = pointGP(:) + R(numCP) * COORDS(:, numCP)
                    end do
                    !! Distance to rotation axis
                    pointA(:) = load_additionalInfos(load_addinfos_count:    &
                        &                    load_addinfos_count + MCRD)
                    pointB(:) = load_additionalInfos(load_addinfos_count + MCRD:   &
                        &                    load_addinfos_count + 2 * MCRD)

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
                            RHS(kk) = RHS(kk) + DENSITY * ADLMAG(i_load) ** two * vectR(j) * R(numCP) * dvol
                        end do
                    end do
                end if
            end if
            kload = kload + load_target_nbelem(i_load)
            load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(i_load)
        end do
    end do   !! End of the loop on integration points

    !! Loop for load : find boundary loads
    load_addinfos_count = 1
    kk = 0
    do i_load = 1, nb_load
        if ((JDLTYPE(i_load) > 9 .AND. JDLTYPE(i_load) < 100) .AND.   &
            &   ANY(indDLoad(kk + 1:kk + load_target_nbelem(i_load)) == JELEM)) then
        !! Define Gauss points coordinates and weights on surf(3D)/edge(2D)
        call LectCle (JDLType(i_load), KNumFace, KTypeDload)

        if (KTypeDload == 4) then
            !! Get Index of nodal distribution
            iField = int(load_additionalInfos(load_addinfos_count))
        end if
        call Gauss (NbPtInt, MCRD, GaussPdsCoord, KNumFace)

        FbL(:) = zero
        do n = 1, NbPtInt ** (MCRD - 1)

            call shapPress(R, VectNorm, DetJac, COORDS,        &
                &   GaussPdsCoord(2:, n), MCRD, KNumFace, KTypeDload)

            dvol = GaussPdsCoord(1, n) * DetJac

            !! Non-uniform pressure case
            if (KTypeDload == 4) then
                f_mag = 0
                do k3 = 1, NNODE
                    f_mag = f_mag + n_dist_elem(iField, k3) * R(k3)
                end do
                !! Uniform pressure case
            else
                f_mag = ADLMAG(i_load)
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
    kk = kk + load_target_nbelem(i_load)
    load_addinfos_count = load_addinfos_count + nb_load_additionalInfos(i_load)
end do

end subroutine UELVEC_byCP

