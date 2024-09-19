!! Copyright 2022 Arnaud Duval

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


!! Compute gradient for the embedded solid element U10
!! Implemented for 3D case
!! TODO : should be adapted to work with 2D case

!! Details can be found in :
!! Marie Guerder. IsoGeometric analysis and shape optimisation of aircraft compressor 
!! blades. Mechanics [physics.med-ph]. Université de Lyon, 2022.
!! https://theses.hal.science/tel-03859679

subroutine gradUELMAT10adj(activeElementMap, nb_elemMap,    &
                    &      sctr,        &
                    &      Uelem, UAelem,               &
                    &      nadj, mcrd, nnode, nnodemap, nb_cp, jelem, nbint, &
                    &      coords, coordsall, &
                    &      tensor, material_properties,     &
                    &      density, nb_load, indDload, load_target_nbelem, jdltype,        &
                    &      adlmag, load_additionalInfos, nb_load_additionalInfos, &
                    &      computeWint, computeWext,    &
                    &      gradWint, gradWext)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! Input arguments
    integer, intent(in) :: nadj, mcrd, nnode, nnodemap, nb_cp, nbint, jelem
    double precision, intent(in) :: coords
    dimension coords(mcrd, nnode)
    double precision, intent(in) :: coordsall
    dimension coordsall(3, nb_cp)       !! nb_cp = total number of control points
    character(len=*), intent(in) :: tensor
    double precision, intent(in) :: material_properties, density
    dimension material_properties(2)

    integer, intent(in) :: sctr
    dimension sctr(nnode)

    double precision, intent(in) :: Uelem, UAelem
    dimension Uelem(3, nnode), UAelem(3, nnode, nadj)

    integer, intent(in) :: activeElementMap, nb_elemMap
    dimension activeElementMap(nb_elemMap)

    logical, intent(in) :: computeWint, computeWext

    integer, intent(in) :: nb_load, indDload, load_target_nbelem, jdltype
    integer, intent(in) :: nb_load_additionalInfos
    dimension load_target_nbelem(nb_load), jdltype(nb_load)
    dimension indDLoad(sum(load_target_nbelem))

    double precision, intent(in) :: adlmag, load_additionalInfos
    dimension adlmag(nb_load), load_additionalInfos(nb_load_additionalInfos)

    !! Outputs
    double precision, intent(out) :: gradWint, gradWext
    dimension gradWint(nadj, 3, nb_cp)
    dimension gradWext(nadj, 3, nb_cp)

    !! local variables
    integer :: ntens
    integer :: isave    !! index of last extracted hull element

    !! Quadrature
    integer :: nbPtInt
    double precision :: GaussPdsCoord, PtGauss
    dimension GaussPdsCoord(mcrd+1, nbint), PtGauss(mcrd+1)
    integer :: igp

    double precision :: detjac, dvol

    !! Materiel behaviour
    double precision :: E, nu, lambda, mu

    !! Embedded solid
    double precision :: R, dRdtheta, ddRddtheta
    dimension R(nnode), dRdtheta(nnode, 3), ddRddtheta(nnode, 6)
    double precision :: theta
    dimension theta(3)
    double precision :: dRdx
    dimension dRdx(nnode, 3)

    !! Hull
    double precision :: N, dNdxi, ddNddxi
    dimension N(nnodemap), dNdxi(nnodemap, 3), ddNddxi(nnodemap, 6)

    double precision :: xi
    dimension xi(3)

    !! Mapping parent element -> embedded parameter space
    double precision :: dthetadtildexi, det_dthetadtildexi
    dimension dthetadtildexi(3, 3)

    !! Mapping embedded parametric space -> hull parametric space
    double precision :: dxidtheta, dthetadxi, det_dxidtheta
    dimension dxidtheta(3, 3), dthetadxi(3, 3)

    !! Elements infos
    double precision :: coordsmap
    dimension coordsmap(mcrd, nnodemap)
    integer :: sctr_map
    dimension sctr_map(nnodemap)

    !! Mapping hull parametric space -> physical space
    double precision :: dxdxi, dxidx, det_dxdxi
    dimension dxdxi(3, 3), dxidx(3, 3)

    !! Mapping embedded parametric space -> physical space
    double precision :: dthetadx
    dimension dthetadx(3, 3)
    double precision :: dxdtheta
    dimension dxdtheta(3, 3)

    !! disp/strain/stress fields
    double precision :: ddsdde
    dimension ddsdde(2*mcrd, 2*mcrd)
    double precision :: dUdtheta, dUAdtheta
    dimension dUdtheta(3,3), dUAdtheta(3, 3, nadj)
    double precision :: strain, stress
    dimension strain(2*mcrd), stress(2*mcrd)
    double precision :: strainAdj, stressAdj
    dimension strainAdj(2*mcrd, nadj), stressAdj(2*mcrd, nadj)
    double precision :: coef1, coef2
    double precision :: work
    dimension work(nadj)
    double precision :: UA
    dimension UA(3, nadj)
    double precision :: tr_strain

    !! Derivatives w.r.t embbeded control points P
    double precision :: DdxdxiDP, DdxidthetaDP      !! mappings
    dimension DdxdxiDP(3, 3, 3), DdxidthetaDP(3, 3, 3)
    double precision :: DdxidxDP, DdthetadxiDP      !! inverse mappings
    dimension DdxidxDP(3, 3, 3), DdthetadxiDP(3, 3, 3)
    double precision :: dJdP        !! jacobian determinant
    dimension dJdP(3, nnode)
    double precision :: DdRdxDP
    dimension DdRdxDP(nnode, 3, 3)

    double precision :: DdUAdxDP, DdUdxDP
    dimension DdUAdxDP(3, 3, 3, nadj), DdUdxDP(3, 3, 3)

    double precision :: dEAdP, dEdP, dSdP
    dimension dEAdP(2*mcrd, 3, nadj), dEdP(2*mcrd, 3), dSdP(2*mcrd, 3)

    double precision :: dEAdP_S, dSdP_EA
    dimension dEAdP_S(3, nadj), dSdP_EA(3, nadj)

    double precision dkronkdP, DdxdthetaDP
    dimension dkronkdP(3, 3, 3), DdxdthetaDP(3, 3, 3)

    !! Derivatives w.r.t hull control points Q
    double precision :: DdxdxiDQ, DdxidxDQ
    dimension DdxdxiDQ(3, 3, 3), DdxidxDQ(3, 3, 3)

    double precision :: dJdQ        !! jacobian determinant
    dimension dJdQ(3, nnodemap)
    double precision :: DdRdxDQ
    dimension DdRdxDQ(nnode, 3, 3)

    double precision :: DdUAdxDQ, DdUdxDQ
    dimension DdUAdxDQ(3, 3, 3, nadj), DdUdxDQ(3, 3, 3)

    double precision dkronkdQ
    dimension dkronkdQ(3, 3, 3)

    double precision :: dEAdQ, dEdQ, dSdQ
    dimension dEAdQ(2*mcrd, 3, nadj), dEdQ(2*mcrd, 3), dSdQ(2*mcrd, 3)

    double precision :: dEAdQ_S, dSdQ_EA
    dimension dEAdQ_S(3, nadj), dSdQ_EA(3, nadj)

    !! loading
    integer :: loadcount, kload
    double precision :: pointGP, vectR, vectAG, vectD, pointA, pointB, scal
    dimension pointGP(mcrd), vectR(mcrd), vectAG(mcrd), vectD(mcrd), pointA(mcrd), pointB(mcrd)
    double precision :: dvectRdP, dFdP, dxdP, dxdP_x_D
    dimension dvectRdP(mcrd, mcrd), dFdP(mcrd, mcrd), dxdP(mcrd, mcrd), dxdP_x_D(mcrd)
    double precision :: dxdQ, dxdQ_x_D, dvectRdQ, dFdQ
    dimension dxdQ(mcrd, mcrd), dxdQ_x_D(mcrd), dvectRdQ(mcrd, mcrd), dFdQ(mcrd, mcrd)



    !! Voigt convention
    integer :: voigt
    dimension voigt(6,2)

    !! Temporary storage
    double precision :: temp, temp1, temp2, tempa, tempb, tempc, tempd
    dimension temp(3, 3), temp1(nnode, 3), temp2(nnode, 3)
    dimension tempa(3, 3), tempb(3, 3), tempc(3, 3), tempd(3, 3)

    !! Various loop variables
    integer :: i, j, k, ij, icp, inodemap, iA, b, idim, jdim, kdim, iload
    integer :: imcp

    !! Initialization

    voigt(:,1) = (/ 1,2,3,1,1,2 /)
    voigt(:,2) = (/ 1,2,3,2,3,3 /)

    ntens = 2*mcrd      ! Size of stiffness matrix
    nbPtInt = int(nbint**(1.0/float(mcrd)))     ! Nb of quadrature pts per direction
    if (nbPtInt**mcrd < nbint) nbPtInt = nbPtInt + 1

    !! Compute Gauss pts coordinates and weights
    call Gauss(nbPtInt, mcrd, GaussPdsCoord, 0)

    !! Gradients
    ! gradWint(:,:,:) = zero
    ! gradWext(:,:,:) = zero

    !! Material behaviour
    E = material_properties(1)
    nu = material_properties(2)
    lambda = E*nu/(one+nu)/(one-two*nu)
    mu     = E/two/(one+nu)
    !! Will be necessary for further 2D implementations
    if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)

    !! Material behavior
    call material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde)


    !! Computation
    isave = 0

    !! Loop on integration points
    do igp = 1, nbint
        !! Embedded solid
        !! ==============

        ptGauss(:) = GaussPdsCoord(2:, igp)
        !! DetJac = GaussPdsCoord(1, igp) ! a gerer de façon plus lisible

        !! Compute parametric coordinates in embedded element from parent element
        !!theta(:) = zero
        do i = 1, mcrd
            theta(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i)) * ptGauss(i)  &
                &     + (Ukv_elem(2,i) + Ukv_elem(1,i))) * 0.5
        enddo

        !! Compute NURBS basis function and derivative
        call evalnurbs_w2ndDerv(theta, R, dRdtheta, ddRddtheta)

        !! Compute coordinates in hull parametric space
        xi(:) = zero
        do icp = 1, nnode
            xi(:) = xi(:) + R(icp)*coords(:, icp)
        enddo

        !! Gradient of mapping : parent element -> embedded parametric space
        !! tildexi : parametric coord in parent element [-1,1]^d
        dthetadtildexi(:,:) = zero
        do i = 1, dim_patch
            dthetadtildexi(i,i) = 0.5d0 * (Ukv_elem(2, i) - Ukv_elem(1, i))
        enddo

        call MatrixDet(dthetadtildexi, det_dthetadtildexi, mcrd)

        !! Gradient of mapping : embedded parametric space -> hull parametric space
        dxidtheta(:,:) = zero
        do icp = 1, nnode
            do i = 1, dim_patch
                dxidtheta(:, i) = dxidtheta(:, i) + dRdtheta(icp, i) * coords(:, icp)
            enddo
        enddo

        call MatrixInv(dthetadxi, dxidtheta, det_dxidtheta, mcrd)

        !! Hull
        !! ====

        !! Get active element number
        call updateMapElementNumber(xi(:))

        !! Evaluate NURBS basis functions and derivative of hull
        call evalnurbs_mapping_w2ndDerv(xi(:), N(:), dNdxi(:,:), ddNddxi(:,:))

        !! Extract CP coordinates of the hull
        if (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:, current_map_elem)
            do icp = 1, nnodemap
                coordsmap(:, icp) = coordsall(:, sctr_map(icp))
            enddo
            isave = current_map_elem
        endif

        !! Gradient of mapping : hull paraletric space -> physical space
        dxdxi(:, :) = zero
        do icp = 1, nnodemap
            do i = 1, dim_patch
                dxdxi(:, i) = dxdxi(:, i) + dNdxi(icp, i) * coordsmap(:, icp)
            enddo
        enddo

        call MatrixInv(dxidx, dxdxi, det_dxdxi, mcrd)

        !! Composition : hull x embedded solid
        !! ===================================

        !! Mapping
        call MulMat(dthetadxi, dxidx, dthetadx, mcrd, mcrd, mcrd)

        !! Basis functions composition
        call Mulmat(dRdtheta, dthetadx, dRdx, nnode, mcrd, mcrd)

        !! Compute product of all mapping determinants
        detjac = det_dxdxi * det_dxidtheta * det_dthetadtildexi

        call MulMat(dxdxi(:,:), dxidtheta(:,:), dxdtheta(:,:), mcrd, mcrd, mcrd)

        !! Compute disp and adjoint derivatives
        dUdtheta(:,:) = zero
        do i = 1, mcrd
            do icp = 1, nnode
                dUdtheta(:, i) = dUdtheta(:, i) + dRdtheta(icp, i) * Uelem(:, icp)
            enddo
        enddo

        dUAdtheta(:,:,:) = zero
        do iA = 1, nadj
            do i = 1, mcrd
                do icp = 1, nnode
                    dUAdtheta(:, i, iA) = dUAdtheta(:, i, iA) + dRdtheta(icp, i) * UAelem(:, icp, iA)
                enddo
            enddo
        enddo

        !! Compute state strain and stress
        strain(:) = zero
        stress(:) = zero
        do ij = 1, ntens
            i = voigt(ij, 1); j = voigt(ij, 2)
            if (i==j) then
                call dot(dUdtheta(i, :), dthetadx(:, i), strain(ij))
            else
                call dot(dUdtheta(i, :), dthetadx(:, j), coef1)
                call dot(dUdtheta(j, :), dthetadx(:, i), coef2)
                strain(ij) = coef1 + coef2
            endif
        enddo
        tr_strain = strain(1) + strain(2) + strain(3)
        call MulVect(ddsdde, strain, stress, ntens, ntens)

        !! Compute adjoint strain and stress
        !! A PRIORI PAS BESOIN DE STRESSADJ
        strainAdj(:,:) = zero
        stressAdj(:,:) = zero
        do iA = 1, nadj
            do ij = 1, ntens
                i = voigt(ij, 1); j = voigt(ij, 2)
                if (i==j) then
                    call dot(dUAdtheta(i, :, iA), dthetadx(:, i), strainAdj(ij, iA))
                else
                    call dot(dUAdtheta(i, :, iA), dthetadx(:, j), coef1)
                    call dot(dUAdtheta(j, :, iA), dthetadx(:, i), coef2)
                    strainAdj(ij, iA) = coef1 + coef2
                endif
            enddo
            call MulVect(ddsdde, strainAdj(:, iA), stressAdj(:, iA), ntens, ntens)
        enddo

        !! Compute local work
        work(:) = zero
        do ij = 1, ntens
            work(:) = work(:) + strainAdj(ij, :)*stress(ij)
        enddo

        !! Derivatives w.r.t embedded CP
        !! =============================

        do icp = 1, nnode
            !! Compute mapping derivatives
            DdxidthetaDP(:,:,:) = zero
            do i = 1, mcrd
                do j = 1, mcrd
                    DdxidthetaDP(i,j,i) = dRdtheta(icp, j)
                enddo
            enddo

            DdxdxiDP(:,:,:) = zero
            do inodemap = 1, nnodemap
                do i = 1, mcrd
                    do j = 1, mcrd
                        do k = 1, mcrd
                            if (j==k) then
                                DdxdxiDP(i,j,k) = DdxdxiDP(i,j,k) + ddNddxi(inodemap, j) * coordsmap(i, inodemap)
                            else
                                DdxdxiDP(i,j,k) = DdxdxiDP(i,j,k) + ddNddxi(inodemap, j+k+1) * coordsmap(i, inodemap)
                            endif
                        enddo
                    enddo
                enddo
            enddo
            DdxdxiDP(:,:,:) = DdxdxiDP(:,:,:) * R(icp)

            !! Compute inverse mapping derivatives
            DdxidxDP(:,:,:) = zero
            DdthetadxiDP(:,:,:) = zero

            do k = 1, mcrd    !! Loop on CP coordinates
                call MulMat(dxidx(:,:), DdxdxiDP(:, :, k), temp(:,:), mcrd, mcrd, mcrd)
                call mulmat(temp(:,:), dxidx(:,:), DdxidxDP(:, :, k), mcrd, mcrd, mcrd)
                call MulMat(dthetadxi(:,:), DdxidthetaDP(:, :, k), temp(:,:), mcrd, mcrd, mcrd)
                call MulMat(temp(:,:), dthetadxi(:,:), DdthetadxiDP(:, :, k), mcrd, mcrd, mcrd)
            enddo
            DdxidxDP(:, :, :) = -1.D0 * DdxidxDP(:, :, :)
            DdthetadxiDP(:, :, :) = -1.D0 * DdthetadxiDP(:, :, :)

            !! Compute derivative of jacobian determinant
            !! Value of dJdP is stored for each control point for later use in gradWext computation
            dJdP(:, icp) = zero
            do i = 1, mcrd
                do j = 1, mcrd
                    do k = 1, mcrd
                        dJdP(k, icp) = dJdP(k, icp) + dxidx(i, j)*DdxdxiDP(j, i, k) + dthetadxi(i, j)*DdxidthetaDP(j, i, k)
                    enddo
                enddo
            enddo
            dJdP(:, icp) = dJdP(:, icp) * detjac

            do iA = 1, nadj
                gradWint(iA, : , sctr(icp)) &
                 & = gradWint(iA, : , sctr(icp)) - work(iA)*dJdP(:, icp)*GaussPdsCoord(1, igp)
            enddo

            !! Compute derivative of embbeded element shape function gradient
            DdRdxDP(:, :, :) = zero
            do k = 1, mcrd     !! Loop on coordinates of current CP
                call MulMat(dRdtheta(:,:), DdthetadxiDP(:,:,k), temp1(:,:), nnode, mcrd, mcrd)
                call MulMat(temp1(:,:), dxidx(:,:), DdRdxDP(:,:,k), nnode, mcrd, mcrd)

                !! TO DO : this line does not depend on k and can be computed outside of the loop
                call MulMat(dRdtheta(:,:), dthetadxi(:,:), temp1(:,:), nnode, mcrd, mcrd)
                call MulMat(temp1(:,:), DdxidxDP(:,:,k), temp2(:,:), nnode, mcrd, mcrd)
                DdRdxDP(:,:,k) = DdRdxDP(:,:,k) + temp2(:,:)
            enddo

            !! Compute derivative of adjoint displacement gradient
            DdUAdxDP(:, :, :, :) = zero
            do iA = 1, nadj
                do k = 1, mcrd     !! Loop on coordinates of current CP
                    do b = 1, nnode     !! Loop on CP where adjoint disp is supported
                        do i = 1, mcrd
                            do j = 1, mcrd
                                DdUAdxDP(i, j, k, iA) = DdUAdxDP(i, j, k, iA) + DdRdxDP(b, j, k) * UAelem(i, b, iA)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            !! Compute derivative of displacement gradient
            DdUdxDP(:, :, :) = zero
            do k = 1, mcrd     !! Loop on coordinates of current CP
                do b = 1, nnode     !! Loop on CP where disp is supported
                    do i = 1, mcrd
                        do j = 1, mcrd
                            DdUdxDP(i, j, k) = DdUdxDP(i, j, k) + DdRdxDP(b, j, k) * Uelem(i, b)
                        enddo
                    enddo
                enddo
            enddo

            !! Compute strain and adjoint strain derivatives
            do ij = 1, ntens
                i = voigt(ij, 1)
                j = voigt(ij, 2)
                do k = 1, mcrd     !! Loop on cordinates of current CP
                    if (i == j) then
                        dEdP(ij, k) = DdUdxDP(i, i, k)
                        do iA = 1, nadj
                            dEAdP(ij, k, iA) =  DdUAdxDP(i, i, k, iA)
                        enddo
                    else
                        dEdP(ij, k) = DdUdxDP(i, j, k) + DdUdxDP(j, i, k)
                        do iA = 1, nadj
                            dEAdP(ij, k, iA) =  DdUAdxDP(i, j, k, iA) + DdUAdxDP(j, i, k, iA)
                        enddo
                    endif
                enddo
            enddo

            !! Compute derivative of kronecker operator for material law

            DdxdthetadP(:,:,:) = zero
            dkronkdP(:,:,:) = zero

            do k = 1, mcrd
                call MulMat(DdxdxiDP(:,:,k), dxidtheta(:,:), tempa(:,:), mcrd, mcrd, mcrd)
                call MulMat(dxdxi, DdxidthetaDP(:,:,k), tempb(:,:), mcrd, mcrd, mcrd)
                call MulMat(dthetadx(:,:), tempa(:,:) + tempb(:,:), dkronkdP(:,:,k), mcrd, mcrd, mcrd)

                call MulMat(dthetadxi(:,:), DdxidthetaDP(:,:,k), tempa(:,:), mcrd, mcrd, mcrd)
                call MulMat(tempa(:,:), dthetadx(:,:), tempb(:,:), mcrd, mcrd, mcrd)

                call MulMat(dthetadx(:,:), DdxdxiDP(:,:,k), tempa(:,:), mcrd, mcrd, mcrd)
                call MulMat(tempa(:,:), dxidx(:,:), tempc(:,:), mcrd, mcrd, mcrd)
                call MulMat(-1.0*(tempb(:,:)+tempc(:,:)), dxdtheta(:,:), tempd(:,:), mcrd, mcrd, mcrd)

                dkronkdP(:,:,k) = dkronkdP(:,:,k) + tempd(:,:)

            enddo


            !! Compute stress and its derivative
            do k = 1, mcrd
                call MulVect(ddsdde(:,:), dEdP(:, k), dSdP(:, k), ntens, ntens)
            enddo

            dEAdP_S(:,:) = zero
            dSdP_EA(:,:) = zero
            do iA = 1, nadj
                do k = 1, mcrd
                    do ij = 1, ntens
                        dEAdP_S(k, iA) = dEAdP_S(k, iA) + dEAdP(ij, k, iA)*stress(ij)
                        dSdP_EA(k, iA) = dSdP_EA(k, iA) + strainAdj(ij, iA)*dSdP(ij, k)
                    enddo
                enddo
            enddo

            do iA = 1, nadj
                gradWint(iA, : , sctr(icp)) = gradWint(iA, : , sctr(icp)) &
                 & - dEAdP_S(:, iA) * detJac * GaussPdsCoord(1, igp) &
                 & - dSdP_EA(:, iA) * detJac * GaussPdsCoord(1, igp)
            enddo

            !! Add contribution of the derivative of material law
            !! re-use variable dSdP_EA
            !! TODO : verify if dkronk is symmetric (if it's the case, voigt notation should be used)
            dSdP_EA(:,:) = zero
            do iA = 1, nadj
                do k = 1, mcrd
                    do i = 1, mcrd
                        do j =1, mcrd
                            if (i == j) then
                                dSdP_EA(k, iA) = dSdP_EA(k, iA) + strainAdj(i, iA)*lambda*tr_strain*dkronkdP(i,j,k)
                            else
                                dSdP_EA(k, iA) = dSdP_EA(k, iA) + 0.5*strainAdj(i+j+1, iA)*lambda*tr_strain*dkronkdP(i,j,k)
                            endif
                        enddo
                    enddo
                enddo
            enddo

            do iA = 1, nadj
                gradWint(iA, : , sctr(icp)) = gradWint(iA, : , sctr(icp)) &
                    & - dSdP_EA(:, iA) * detJac * GaussPdsCoord(1, igp)
            enddo
        enddo   !! End loop on embedded control points icp

        !! Derivatives w.r.t hull CP
        !! =========================

        do imcp = 1, nnode_map
            !! Compute mapping derivatives
            DdxdxiDQ(:,:,:) = zero
            do i = 1, mcrd
                do j = 1, mcrd
                    DdxdxiDQ(i,j,i) = dNdxi(imcp, j)
                enddo
            enddo

            !! Compute inverse mapping derivatives
            DdxidxDQ(:,:,:) = zero
            do k = 1, mcrd     !! Loop on CP coordinates
                call MulMat(dxidx(:,:), DdxdxiDQ(:, :, k), temp(:,:), mcrd, mcrd, mcrd)
                call mulmat(temp(:,:), dxidx(:,:), DdxidxDQ(:, :, k), mcrd, mcrd, mcrd)
            enddo
            DdxidxDQ(:, :, :) = -1.D0 * DdxidxDQ(:, :, :)

            !! Compute derivative of jacobian determinant
            !! Value of dJdQ is storef for each control point for later use in gradWext computation
            dJdQ(:, imcp) = zero
            do i =1, mcrd
                do j =1, mcrd
                    do k = 1, mcrd
                        dJdQ(k, imcp) = dJdQ(k, imcp) + dxidx(i, j) * DdxdxiDQ(j, i, k)
                    enddo
                enddo
            enddo
            dJdQ(:, imcp) = dJdQ(:, imcp) * detjac

            do iA = 1, nadj
                gradWint(iA, : , sctr_map(imcp)) &
                 & = gradWint(iA, : , sctr_map(imcp)) - work(iA)*dJdQ(:, imcp)*GaussPdsCoord(1, igp)
            enddo

            !! Compute derivative of embedded element shape function gradient
            DdRdxDQ(:, :, :) = zero
            call MulMat(dRdtheta(:,:), dthetadxi(:,:), temp1, nnode, mcrd, mcrd)
            do k = 1, mcrd !! Loop on coordinates of current hull CP
                call MulMat(temp1(:,:), DdxidxDQ(:,:,k), DdRdxDQ(:,:,k), nnode, mcrd, mcrd)
            enddo

            !! Compute derivative of adjoint displacement gradient
            DdUAdxDQ(:, :, :, :) = zero
            do iA = 1, nadj
                do k = 1, mcrd     !! Loop on coordinates of current hull CP
                    do b = 1, nnode     !! Loop on CP where adjoint disp is supported
                        do i = 1, mcrd
                            do j = 1, mcrd
                                DdUAdxDQ(i, j, k, iA) = DdUAdxDQ(i, j, k, iA) + DdRdxDQ(b, j, k) * UAelem(i, b, iA)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            !! Compute derivative of displacement gradient
            DdUdxDQ(:,:,:) = zero
            do k = 1, mcrd     !! Loop on coordinates of current CP
                do b = 1, nnode     !! Loop on CP where disp is supported
                    do i = 1, mcrd
                        do j = 1, mcrd
                            DdUdxDQ(i, j, k) = DdUdxDQ(i, j, k) + DdRdxDQ(b, j, k) * Uelem(i, b)
                        enddo
                    enddo
                enddo
            enddo

            !! Compute strain and adjoint strain derivatives
            do ij = 1, ntens
                i = voigt(ij, 1)
                j = voigt(ij, 2)
                do k = 1, mcrd     !! Loop on cordinates of current CP
                    if (i == j) then
                        dEdQ(ij, k) = DdUdxDQ(i, i, k)
                        do iA = 1, nadj
                            dEAdQ(ij, k, iA) =  DdUAdxDQ(i, i, k, iA)
                        enddo
                    else
                        dEdQ(ij, k) = DdUdxDQ(i, j, k) + DdUdxDQ(j, i, k)
                        do iA = 1, nadj
                            dEAdQ(ij, k, iA) =  DdUAdxDQ(i, j, k, iA) + DdUAdxDQ(j, i, k, iA)
                        enddo
                    endif
                enddo
            enddo

            !! Compute derivative of kronecker operator for material law
            dkronkdQ(:,:,:) = zero

            do k = 1, mcrd
                call MulMat(dthetadxi(:,:), DdxidxDQ(:,:,k), tempa(:,:), mcrd, mcrd, mcrd)
                call MulMat(tempa(:,:), dxdtheta(:,:), dkronkdQ(:,:,k), mcrd, mcrd, mcrd)

                call MulMat(dthetadx(:,:), DdxdxiDQ(:,:,k), tempb, mcrd, mcrd, mcrd)
                call MulMat(tempb(:,:), dxidtheta(:,:), tempc(:,:), mcrd, mcrd, mcrd)

                dkronkdQ(:,:,k) = dkronkdQ(:,:,k) + tempc(:,:)
            enddo

            !! Compute stress and its derivative
            do k = 1, mcrd
                call MulVect(ddsdde(:,:), dEdQ(:, k), dSdQ(:, k), ntens, ntens)
            enddo

            dEAdQ_S(:,:) = zero
            dSdQ_EA(:,:) = zero
            do iA = 1, nadj
                do k = 1, mcrd
                    do ij = 1, ntens
                        dEAdQ_S(k, iA) = dEAdQ_S(k, iA) + dEAdQ(ij, k, iA)*stress(ij)
                        dSdQ_EA(k, iA) = dSdQ_EA(k, iA) + strainAdj(ij, iA)*dSdQ(ij, k)
                    enddo
                enddo
            enddo

            do iA = 1, nadj
                gradWint(iA, : , sctr_map(imcp)) = gradWint(iA, : , sctr_map(imcp)) &
                 & - dEAdQ_S(:, iA) * detJac * GaussPdsCoord(1, igp) &
                 & - dSdQ_EA(:, iA) * detJac * GaussPdsCoord(1, igp)
            enddo

            !! Add contribution of the derivative of material law
            !! re-use variable dSdQ_EA
            !! TODO : verify if dkronk is symmetric (if it's the case, voigt notation should be used)
            dSdQ_EA(:,:) = zero
            do iA = 1, nadj
                do k = 1, mcrd
                    do i = 1, mcrd
                        do j =1, mcrd
                            if (i == j) then
                                dSdQ_EA(k, iA) = dSdQ_EA(k, iA) + strainAdj(i, iA)*lambda*tr_strain*dkronkdQ(i,j,k)
                            else
                                dSdQ_EA(k, iA) = dSdQ_EA(k, iA) + 0.5*strainAdj(i+j+1, iA)*lambda*tr_strain*dkronkdQ(i,j,k)
                            endif
                        enddo
                    enddo
                enddo
            enddo

            do iA = 1, nadj
                gradWint(iA, : , sctr_map(imcp)) = gradWint(iA, : , sctr_map(imcp)) &
                    & - dSdQ_EA(:, iA) * detJac * GaussPdsCoord(1, igp)
            enddo

        enddo


        !! Body loads

        if(computeWext) then
            !! Compute adjoint solution
            UA(:,:) = zero
            do iA = 1, nadj
                do icp = 1, nnode
                    UA(:, iA) = UA(:, iA) + R(icp)*UAelem(:, icp, iA)
                enddo
            enddo

            loadcount = 1
            kload = 0
            do iload = 1, nb_load
                if((JDLTYPE(iload) == 101) .and.        &
                    &   any(indDLoad(kload+1: kload+load_target_nbelem(iload)) == JELEM)) then
                    !! Centrifugal body force
                    !! Gauss point location in physical space
                    pointGP(:) = zero
                    do icp = 1, nnodemap
                        pointGP(:) = pointGP(:) + N(icp)*COORDSmap(:, icp)
                    enddo
                    !! Distance to rotation axis
                    pointA(:) = zero
                    pointA(:mcrd) = load_additionalInfos(loadcount:loadcount+mcrd)
                    loadcount = loadcount + mcrd
                    pointB(:) = zero
                    pointB(:mcrd) = load_additionalInfos(loadcount:loadcount+mcrd)
                    loadcount = loadcount + mcrd

                    vectD(:) = pointB(:) - pointA(:)
                    vectD(:) = vectD(:)/sqrt(sum(vectD(:)*vectD(:)))
                    vectAG(:) = pointGP(:) - pointA(:)
                    call dot(vectAG(:), vectD(:), scal)
                    vectR(:) = vectAG(:) - scal*vectD(:)


                    !! Derivatives / control points of embedded entity
                    !! -----------------------------------------------
                    do icp = 1, nnode

                        dxdP(:,:) = zero
                        do idim = 1, mcrd
                            do jdim = 1, mcrd
                                dxdP(idim, jdim) = dxdP(idim, jdim) + R(icp)*dxdxi(idim, jdim)
                            enddo
                        enddo

                        dxdP_x_D(:) = zero
                        do idim = 1, dim_patch
                            dxdP_x_D(:) = dxdP_x_D(:) + dxdP(idim,:)*vectD(idim)
                        enddo

                        do idim = 1, dim_patch
                            dvectRdP(idim,:) = dxdP(idim,:) - dxdP_x_D(:)*vectD(idim)
                        enddo


                        do idim = 1, dim_patch
                            do jdim = 1, dim_patch
                                dFdP(idim,jdim) = density*(adlmag(iload)**two)*       &
                                    &   (dvectRdP(idim,jdim)*detjac + vectR(idim)*dJdP(jdim, icp)) &
                                    &       * GaussPdsCoord(1,igp)
                            enddo
                        enddo

                        do iA = 1, nadj
                            do idim = 1, mcrd
                                gradWext(iA,:,sctr(icp)) = gradWext(iA,:,sctr(icp)) + &
                                    &   UA(idim, iA)*dFdP(idim, :)
                            enddo
                        enddo

                    enddo

                    !! Derivative / control points of hull
                    !! -----------------------------------
                    do imcp = 1, nnode_map
                        dxdQ(:,:) = zero
                        do idim = 1, mcrd
                            dxdQ(idim, idim) = N(imcp)
                        enddo

                        dxdQ_x_D(:) = zero
                        do idim = 1, dim_map
                            dxdQ_x_D(:) = dxdQ_x_D(:) + dxdQ(idim,:)*vectD(idim)
                        enddo

                        do idim = 1, dim_map
                            dvectRdQ(idim,:) = dxdQ(idim,:) - dxdQ_x_D(:)*vectD(idim)
                        enddo

                        do idim = 1, dim_patch
                            do jdim = 1, dim_map
                                dFdQ(idim,jdim) = density*(adlmag(iload)**two)*       &
                                    &   (dvectRdQ(idim,jdim)*detjac + vectR(idim)*dJdQ(jdim, imcp)) &
                                    &       * GaussPdsCoord(1,igp)
                            enddo
                        enddo

                        do iA = 1, nadj
                            do idim = 1, mcrd
                                gradWext(iA,:,sctr_map(imcp)) = gradWext(iA,:,sctr_map(imcp)) + &
                                    &   UA(idim, iA)*dFdQ(idim, :)
                            enddo
                        enddo

                    enddo

                endif
                kload = kload + load_target_nbelem(iload)
            enddo


        endif

    enddo       !! End loop on Gauss point igp
end subroutine gradUELMAT10adj


