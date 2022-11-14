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

subroutine gradUELMAT10adj(Uelem, UAelem,               &
                    &      nadj, mcrd, nnode, nnodemap, nb_cp, nbint, &
                    &      coords, coordsall, &
                    &      tensor, material_properties,     &
                    &      gradWint_elem, gradWext_elem)

    use parameters
    use nurbspatch
    use embeddedMapping

    implicit none

    !! TODO : DOCUMENTER À QUOI CORRESPOND nb_cp -> EST-CE POUR UN ELEMENT 
    !! DE L'ENVELOPPE OU POUR L'ENVELOPPE ENTIÈRE ?

    !! Input arguments
    integer, intent(in) :: nadj, mcrd, nnode, nnodemap, nb_cp, nbint
    double precision, intent(in) :: coords
    dimension coords(mcrd, nnode)
    double precision, intent(in) :: coordsall
    dimension coordsall(3, nb_cp)
    character(len=*), intent(in) :: tensor
    double precision, intent(in) :: material_properties
    dimension material_properties(2)

    double precision, intent(in) :: Uelem, UAelem
    dimension Uelem(3, nnode), UAelem(3, nnode, nadj)

    !! Outputs
    double precision, intent(out) :: gradWint_elem, gradWext_elem
    dimension gradWint_elem(nadj, 3, nnode)
    dimension gradWext_elem(nadj, 3, nnode)

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

    !! Derivatives w.r.t embbeded control points P
    double precision :: DdxdxiDP, DdxidthetaDP      !! mappings
    dimension DdxdxiDP(3, 3, 3), DdxidthetaDP(3, 3, 3)
    double precision :: DdxidxDP, DdthetadxiDP      !! inverse mappings
    dimension DdxidxDP(3, 3, 3), DdthetadxiDP(3, 3, 3)
    double precision :: dJdP        !! jacobia determinant
    dimension dJdP(3)

    double precision :: DdUAdxDP
    dimension DdUAdxDP(3, 3, 3, nadj)
    
    double precision :: dEAdP
    dimension dEAdP(2*mcrd, 3, nadj)

    !! Voigt convention
    integer :: voigt
    dimension voigt(6,2)

    !! Temporary storage
    double precision :: temp, temp2
    dimension temp(3, 3), temp2(3, 3)

    !! Various loop variables
    integer :: i, j, k, ij, icp, inodemap, iA

    !! Initialization

    voigt(:,1) = (/ 1,2,3,1,1,2 /)
    voigt(:,2) = (/ 1,2,3,2,3,3 /)

    ntens = 2*mcrd      ! Size of stiffness matrix
    nbPtInt = int(nbint**(1.0/float(mcrd)))     ! Nb of quadrature pts per direction
    if (nbPtInt**mcrd < nbint) nbPtInt = nbPtInt + 1

    !! Compute Gauss pts coordinates and weights
    call Gauss(nbPtInt, mcrd, GaussPdsCoord, 0)

    !! Gradients
    gradWint_elem(:,:,:) = zero
    gradWext_elem(:,:,:) = zero

    !! Material behaviour
    !! TODO : use material_lib subroutine
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
            theta(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i)) * PtGauss(i)  &
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

        call MatrixDet(dthetadtildexi, det_dthetadtildexi, 3)

        !! Gradient of mapping : embedded parametric space -> hull parametric space
        dxidtheta(:,:) = zero
        do icp = 1, nnode
            do i = 1, dim_patch
                dxidtheta(:, i) = dxidtheta(:, i) + dRdtheta(icp, i) * coords(:, icp)
            enddo
        enddo

        call MatrixInv(dthetadxi, dxidtheta, det_dxidtheta, 3)

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

        call MatrixInv(dxidx, dxdxi, det_dxdxi, 3)

        !! Composition : hull x embedded solid
        !! ===================================

        !! Mapping
        call MulMat(dthetadxi, dxidx, dthetadx, 3, 3, 3)

        !! Basis functions composition
        call Mulmat(dRdtheta, dthetadx, dRdx, nnode, 3, 3)

        !! Compute product of all mapping determinants
        detjac = det_dxdxi * det_dxidtheta*det_dthetadtildexi

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
                call dot(dUdtheta(i, :), dthetadx(:, i), coef1)
            else
                call dot(dUdtheta(i, :), dthetadx(:, j), coef1)
                call dot(dUdtheta(j, :), dthetadx(:, i), coef2)
                strain(ij) = coef1 + coef2
            endif
        enddo
        call MulVect(ddsdde, strain, stress, ntens, ntens)

        !! Compute adjoint strain and stress
        strainAdj(:,:) = zero
        stressAdj(:,:) = zero
        do iA = 1, nadj
            do ij = 1, ntens
                i = voigt(ij, 1); j = voigt(ij, 2)
                if (i==j) then
                    call dot(dUAdtheta(i, :, iA), dthetadx(:, i), coef1)
                else
                    call dot(dUAdtheta(i, :, iA), dthetadx(:, j), coef1)
                    call dot(dUAdtheta(j, :, iA), dthetadx(:, i), coef2)
                    strain(ij) = coef1 + coef2
                endif
            enddo
            call MulVect(ddsdde, strainAdj(:, iA), stressAdj(:, iA), ntens, ntens)
        enddo

        !! Compute local work
        !! TODO : verify if is needed
        work(:) = zero
        do ij = 1, ntens
            work(:) = work(:) + strainAdj(ij, :)*stress(ij)
        enddo





        !! Derivatives
        !! ===========

        do icp = 1, nnode
            !! Compute mapping derivatives
            DdxidthetaDP(:,:,:) = zero
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        if (i == k) then
                            DdxidthetaDP(i,j,k) = dRdtheta(icp, j)
                        endif
                    enddo
                enddo
            enddo

            DdxdxiDP(:,:,:) = zero
            do inodemap = 1, nnodemap
                do i = 1, 3
                    do j = 1, 3
                        do k = 1, 3
                            if (i==j) then
                                DdxdxiDP(i,j,k) = DdxdxiDP(i,j,k) + ddNddxi(inodemap, i) * coordsmap(k, inodemap)
                            else
                                DdxdxiDP(i,j,k) = DdxdxiDP(i,j,k) + ddNddxi(inodemap, i+j+1) * coordsmap(k, inodemap)
                            endif
                        enddo
                    enddo
                enddo
            enddo
            DdxdxiDP(:,:,:) = DdxdxiDP(:,:,:) * R(icp)

            !! Compute inverse mapping derivatives
            DdxidxDP(:,:,:) = zero
            DdthetadxiDP(:,:,:) = zero

            do i = 1, 3
                call MulMat(dxidx(:,:), DdxdxiDP(i, :, :), temp(:,:), 3, 3, 3)
                call mulmat(temp(:,:), dxidx(:,:), DdxidxDP(i, :, :), 3, 3, 3)
                call MulMat(dthetadxi(:,:), DdxidthetaDP(i, :, :), temp(:,:), 3, 3, 3)
                call MulMat(temp(:,:), dthetadxi(:,:), DdthetadxiDP(i, :, :), 3, 3, 3)
            enddo

            DdxidxDP(:, :, :) = -1.D0 * DdxdxiDP(:, :, :)
            DdthetadxiDP(:, :, :) = -1.D0 * DdthetadxiDP(:, :, :)

            
            !! Compute derivative of jacobian determinant
            dJdP(:) = zero
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3
                        dJdP(i) = dJdP(i) + dxidx(j, k)*DdxdxiDP(i, k, j) + dthetadxi(j, k)*DdxidthetaDP(i, k, j)
                    enddo
                enddo
            enddo
            dJdP(:) = dJdP * detjac

            !! Compute derivative of adjoint displacement
            DdUAdxDP(:, :, :, :) = zero
            do iA = 1, nadj
                do i = 1, 3
                    !! ATTENTION : ERREUR PROBABLE AVEC INVERSION DU 2EME ET
                    !! 3EME INDICE SUR LES DERIVÉES PAR RAPPORT À P (TABLEAUX À 4 INDICES)
                    
                    call MulMat(DdxidthetaDP(i, :, :), dthetadx(:,:), temp(:,:), 3, 3, 3)
                    DdUAdxDP(i, :, :, iA) = DdUAdxDP(i, :, :, iA) + temp(:,:)

                    call MulMat(dxidtheta(:,:), DdthetadxiDP(i, :, :), temp(:,:), 3, 3, 3)
                    call MulMat(temp(:,:), dxidx(:,:), temp2(:,:), 3, 3, 3)
                    DdUAdxDP(i, :, :, iA) = DdUAdxDP(i, :, :, iA) + temp2(:,:)

                    DdUAdxDP(i, :, :, iA) = DdUAdxDP(i, :, :, iA) + DdxidxDP(i, :, :)

                                ! DdxidthetaDP x dthetadxi x dxidx
                                ! + dxidtheta x DdthetadxiDP x dxidx
                                ! + dxidtheta x dthetadxi x DdxidxDP
                enddo
            enddo
            
            !! Compute derivative of adjoint strain
            dEAdP(:,:,:) = zero
            



        enddo

        




        
        dvol = GaussPdsCoord(i, igp)*detjac



    enddo

    call evalnurbs_mapping_w2ndDerv


end subroutine gradUELMAT10adj


