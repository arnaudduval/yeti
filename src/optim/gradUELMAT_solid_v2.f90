!! Copyright 2019-2020 Thibaut Hirschler
!! Copyright 2019 Arnaud Duval

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


!! NOTE : this routine is never called


subroutine gradUELMAT1adjv2(activeDir, Uelem, UAelem, &
    &     nadj, NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS, &
    &     TENSOR, MATERIAL_PROPERTIES, DENSITY, PROPS, JPROPS, nb_load, &
    &     indDLoad, load_target_nbelem, JDLType, ADLMAG, &
    &     load_additionalInfos, nb_load_additionalInfos, &
    &     computeWint, computeWext, &
    &     gradWint_elem, gradWext_elem)
    
    use parameters
    use nurbspatch
    
    implicit none
    
    !! Input arguments
    !! ---------------
    integer, intent(in) :: nadj, NDOFEL, MCRD, NNODE, JELEM, NBINT, JPROPS
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, DENSITY,PROPS
    dimension COORDS(3, NNODE), MATERIAL_PROPERTIES(2), PROPS(JPROPS)
    
    integer, intent(in) :: indDLoad, load_target_nbelem, JDLType, &
        &   nb_load, nb_load_additionalInfos
    double precision, intent(in) :: ADLMAG, load_additionalInfos
    dimension ADLMAG(nb_load), &
        &   load_target_nbelem(nb_load), JDLType(nb_load), &
        &   load_additionalInfos(nb_load_additionalInfos)
    dimension indDLoad(SUM(load_target_nbelem))
      
    double precision, intent(in) :: Uelem, UAelem
    integer, intent(in) :: activeDir
    dimension Uelem(3, NNODE), UAelem(3, NNODE, nadj), activeDir(3)
    
    logical, intent(in) :: computeWint,computeWext
      
    !! Output variables
    !! ----------------
    double precision, intent(out) :: gradWint_elem, gradWext_elem
    dimension gradWint_elem(nadj, 3, NNODE), gradWext_elem(nadj, 3, NNODE)
      
      
    !! Local variables
    !! ---------------
      
    !! Gauss points
    integer :: NbPtInt, nint
    double precision :: GaussPdsCoord, PtGauss, WtGauss
    dimension GaussPdsCoord(MCRD+1, NBINT), PtGauss(MCRD+1)
    
    !! NURBS basis functions
    double precision :: Xi, R, dRdXi
    dimension Xi(3), R(NNODE), dRdXi(NNODE, 3)
    
    !! Mapping: parent element >> parametric space
    double precision :: dXidtildeXi, dtildeXidXi, detdXidtildeXi
    dimension dXidtildeXi(3, 3), dtildeXidXi(3, 3)
    
    !! Mapping: parametric space >> physical space
    double precision :: dXdXi, dXidX, detdXdXi
    dimension dXdXi(3, 3), dXidX(3, 3)
    
    !! Integrals transformation
    double precision :: dVol
    
    !! For material matrix
    integer :: voigt
    dimension voigt(6, 2)
    double precision :: E, nu, lambda, mu, coef
    
    !! State fields
    !! - Displacement
    double precision :: dUdXi
    dimension dUdXi(3, 3)
    !! - Strain
    double precision :: strain, traceStrain, kronk
    dimension strain(2*MCRD), kronk(MCRD, MCRD)
    !! - Stress
    double precision :: stress
    dimension stress(2*MCRD)
    
    !! Adjoint fields
    !! - Displacement
    double precision :: dUAdxi
    dimension dUAdxi(3, 3, nadj)
    !! - Strain
    double precision :: strainAdj, traceStrainAdj
    dimension strainAdj(2*MCRD, nadj)
    !! - Stress
    double precision :: stressAdj
    dimension stressAdj(2*MCRD, nadj)
    !! - Work
    double precision :: work
    dimension work(nadj)
    !! - Solution
    double precision :: UA
    dimension UA(3, nadj)    
    
    !! Derivatives 
    !! - Jacobian determinant
    double precision :: dJdP, tracedXidX
    dimension dJdP(3)
    !! - Strain
    double precision :: dEdP, dEAdP_S, tracedEdP
    dimension dEdP(3, 2*MCRD), dEAdP_S(3, nadj), tracedEdP(3)    
    !! - Stress
    double precision :: dSdP, dKronkdP, dKronkInvdP, dSdP_EA
    dimension dSdP(3, nadj), dKronkdP(3, 3, 3), &
        &   dKronkInvdP(3, 3, 3), dSdP_EA(3, nadj)
    
    !! Misc. for computing strain partial terms
    double precision :: coef1, coef2

    !! Loads
    integer :: nb_load_bnd, nb_load_srf, ind_load_loc, numLoad
    dimension ind_load_loc(nb_load)
    double precision :: VectNorm_U, normV
    !! - Centrifugal load
    integer :: loadcount, nl
    double precision :: pointGP, pointA, pointB, vectD, vectAG, &
        &   vectR, scal, loadF, scalFUA, vectDDUA, dfdP_UA
    dimension pointGP(3), pointA(3), pointB(3), vectD(3), vectAG(3), &
        &   vectR(3), loadF(3), scalFUA(nadj), vectDDUA(3, nadj), &
        &   dfdP_UA(3, nadj)
       
    !! For loops
    integer ntens
    integer k1, k2, i, j, k, l, kk, ll, ij, kl, numCP, &
        &   KTypeDload, KNumFace, idim, iadj
    double precision :: temp, temp1, temp2

    
    !! Initialization
    !! --------------
    
    ntens = 2*MCRD  ! size of stiffness tensor
    NbPtInt = int(NBINT**(1.0 / float(MCRD)))  ! nb gauss pts per dir.
    if (NbPtInt**MCRD < NBINT) NbPtInt = NbPtInt + 1
    
    !! Defining Gauss points coordinates and weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)
    
    !! grad
    gradWint_elem(:, :, :) = zero
    gradWext_elem(:, :, :) = zero
    
    !! Material behavior
    E = MATERIAL_PROPERTIES(1)
    nu = MATERIAL_PROPERTIES(2)
    lambda = E * nu / (one + nu) / (one - two * nu)
    mu = E / two / (one + nu)

    !! Computation
    !! -----------
    
    !! Loop on integration points
    do nint = 1, NBINT
       
        !! Preliminary quantities
        !! ----------------------
        
        !! - Gauss point coordinates and weight
        WtGauss = GaussPdsCoord(1, nint)
        PtGauss(:) = GaussPdsCoord(2:, nint)
        
        !! - NURBS basis functions and derivatives
        R(:) = zero
        dRdXi(:, :) = zero
        Xi(:) = zero
        do idim = 1, MCRD
            Xi(idim) = ((Ukv_elem(2, idim) - Ukv_elem(1, idim)) * PtGauss(idim) + &
                &   (Ukv_elem(2, idim) + Ukv_elem(1, idim))) * 0.5d0
        enddo
        
        call evalnurbs(Xi, R, dRdXi)
        
        !! - Gradient of mapping: parent element >> parametric space
        dXidtildeXi(:, :) = zero
        do idim = 1, MCRD
            dXidtildeXi(idim, idim) = 0.5d0 * (Ukv_elem(2, idim) - Ukv_elem(1, idim))
        enddo
        
        call MatrixDet(dXidtildeXi, detdXidtildeXi, 3)
        
        !! - Gradient of mapping: parametric space >> physical space
        dXdXi(:, :) = zero
        do idim = 1, MCRD
            do numCP = 1, NNODE
                dXdXi(:, idim) = dXdXi(:, idim) + dRdXi(numCP, idim)*COORDS(:, numCP)
            enddo
        enddo
               
        call MatrixInv(dXidX, dXdXi, detdXdXi, 3)
        
        tracedXidX = dXidX(1, 1) + dXidX(2, 2) + dXidX(3, 3)
        
        !! - Transformation of integrals
        dVol = WtGauss * detdXdXi * detdXidtildeXi
        
        if (computeWint) then 
            
            !! Computing disp. and adjoint derivative
            dUdXi(:, :) = zero
            do idim = 1, MCRD
                do numCP = 1, NNODE
                    dUdXi(:, idim) = dUdXi(:, idim) + dRdXi(numCP, idim)*Uelem(:, numCP)
                enddo
            enddo
            
            dUAdxi(:, :, :) = zero
            do iadj = 1, nadj
                do idim = 1, MCRD
                    do numCP = 1, NNODE
                        dUAdxi(:, idim, iadj) = &
                            &   dUAdxi(:, idim, iadj) + &
                            &   dRdXi(numCP, idim) * UAelem(:, numCP, iadj)
                    enddo
                enddo
            enddo
            
            !! Computing state strain
            strain(:) = zero
            do ij = 1, ntens
                i = voigt(ij, 1)
                j = voigt(ij, 2)
                if (i==j) then
                    call dot(dXdXi(:, i), dUdXi(:, i), strain(ij))
                else
                    call dot(dXdXi(:, i), dUdXi(:, j), coef1)
                    call dot(dXdXi(:, j), dUdXi(:, i), coef2)
                    strain(ij) = coef1 + coef2
                endif
            enddo
            
            !! Computing Kronecker delta
            kronk(:, :) = zero
            do j = 1, MCRD
                do i = 1, MCRD
                    call dot(dXidX(i, :), dXidX(j, :), coef)
                    kronk(i, j) = coef
                enddo
            enddo
            
            !! Computing state stress & Kronecker delta
            stress(:) = zero
            traceStrain = strain(1) + strain(2) + strain(3)
            do ij = 1, ntens
                i = voigt(ij, 1)
                j = voigt(ij, 2)                
                stress(ij) = lambda * traceStrain * kronk(i, j) + &
                    &   2 * mu * strain(ij)
            enddo
            
            !! Computing adjoint quantities
            strainAdj(:, :) = zero
            stressAdj(:, :) = zero
            do iadj = 1, nadj
                ! Adjoint strain
                do ij = 1, ntens
                    i = voigt(ij, 1)
                    j = voigt(ij, 2)
                    if (i==j) then
                        call dot(dXdXi(:, i), dUAdxi(:, i, iadj), &
                            &   strainAdj(ij, iadj))
                    else
                        call dot(dXdXi(:, i), dUAdxi(:, j, iadj), coef1)
                        call dot(dXdXi(:, j), dUAdxi(:, i, iadj), coef2)
                        strainAdj(ij, iadj) = coef1 + coef2
                    endif
                enddo 
                ! Adjoint stress
                traceStrainAdj = strainAdj(1, iadj) + strainAdj(2, iadj) + &
                    &   strainAdj(3, iadj)
                do ij = 1, ntens
                    i = voigt(ij, 1)
                    j = voigt(ij, 2)
                    stressAdj(ij, iadj) = lambda * traceStrainAdj * & 
                    &   kronk(i, j) + 2 * mu * strainAdj(ij, iadj)
                enddo
            enddo
            
            !! Computing local work
            work(:) = zero
            do ij = 1, ntens
                work(:) = work(:) + strainAdj(ij, :) * stress(ij)
            enddo
            
        endif ! test computeWint is True
        
        !! Derivatives
        if (computeWint) then
            
            do numCP = 1, NNODE
                
                ! 1. Derivatives of the Jacobian determinant
                dJdP(:) = zero
                do idim = 1, MCRD
                    dJdP(idim) = detdXdXi * tracedXidX * dRdXi(numCP, idim)
                enddo
                
                ! Contribution to gradWint
                do iadj = 1, nadj
                    gradWint_elem(iadj, :, numCP) = &
                        &    gradWint_elem(iadj, :, numCP) - &
                        &    work(iadj) * dJdP(:) * WtGauss * detdXidtildeXi
                enddo

                ! 2. Derivatives of the adjoint strain 
                !    (with dble prod. by stresses)
                dEAdP_S(:, :) = zero
                do iadj = 1, nadj
                    do ij = 1, ntens
                        i = voigt(ij, 1)
                        j = voigt(ij, 2)
                        if (i==j) then
                            dEAdP_S(:, iadj) = dEAdP_S(:, iadj) + &
                                & stress(ij) * dUAdxi(:, i, iadj) * dRdXi(numCP, i)
                        else
                            dEAdP_S(:, iadj) = dEAdP_S(:, iadj) + &
                                &   stress(ij) * dUAdxi(:, i, iadj) * dRdXi(numCP, j) + &
                                &   stress(ij) * dUAdxi(:, j, iadj) * dRdXi(numCP, i)
                        endif
                    enddo
                
                ! Contribution to gradWint
                gradWint_elem(iadj, :, numCP) = gradWint_elem(iadj,:,numCP) - &
                    &   dEAdP_S(:, iadj) * dVol
                enddo
                
                ! 3. Derivatives of the stress 
                !    (with dble prod. by adjoint strains)
                ! - Derivative of the inverse of Kronecker delta
                dKronkInvdP(:, :, :) = zero
                do j = 1, MCRD
                    do i = 1, MCRD
                        dKronkInvdP(:, i, j) = dRdXi(numCP, i) * dXdXi(:, j) + &
                            &   dXdXi(:, i) * dRdXi(numCP, j)
                    enddo
                enddo
                ! - Derivative of Kronecker delta
                dKronkdP(:, :, :) = zero
                do j = 1, MCRD
                    do i = 1, MCRD
                        do l = 1, MCRD
                            do k = 1, MCRD
                                dKronkdP(:, i, j) = dKronkdP(:, i, j) - &
                                    &   kronk(i, k) * dKronkInvdP(:, k, l) * kronk(l, j)
                            enddo
                        enddo
                    enddo
                enddo
                ! - Derivative of the strain
                dEdP(:, :) = zero
                do ij = 1, ntens
                    i = voigt(ij, 1)
                    j = voigt(ij, 2)
                    if (i==j) then
                        dEdP(:, ij) = dEdP(:, ij) + &
                            &   dUdXi(:, i) * dRdXi(numCP, i)
                    else
                        dEdP(:, ij) = dEdP(:, ij) + & 
                            &   dUdXi(:, i) * dRdXi(numCP, j) + &
                            &   dUdXi(:, j) * dRdXi(numCP, i)
                    endif
                enddo
                ! - Compute trace 
                tracedEdP(:) = zero
                do i = 1, MCRD
                    tracedEdP(i) = dEdP(i, 1) + dEdP(i, 2) + dEdP(i, 3)
                enddo
                ! - Derivatives of the stress
                dSdP(:, :) = zero
                do ij = 1, ntens
                    i = voigt(ij, 1)
                    j = voigt(ij, 2)
                    do idim = 1, MRCD
                        dSdP(idim, ij) = &
                            &   lambda * tracedEdP(idim) * kronk(i, j) + &
                            &   lambda * traceStrain * dKronkdP(idim, i, j) + &
                            &   2 * mu * dEdP(idim, ij)
                    enddo
                enddo
                ! - Double product with adjoint strains
                dSdP_EA(:, :) = zero
                do kl = 1, ntens
                    k = voigt(kl, 1)
                    l = voigt(kl, 2)
                    do ij = 1, ntens
                        i = voigt(ij, 1)
                        j = voigt(ij, 2)
                        do iadj = 1, nadj
                            dSdP_EA(:, iadj) = dSdP_EA(:, iadj) + &
                                &   dSdP(:, ij) * strainAdj(kl, iadj)
                        enddo
                    enddo
                enddo
                
                ! Contribution to gradWint
                do iadj = 1, nadj
                    gradWint_elem(iadj, :, numCP) = gradWint_elem(iadj, :, numCP) - & 
                        &   dSdP_EA(:, iadj) * dVol
                enddo
                
            enddo  ! End loop on control points
            
        endif  ! test computeWint is True

        !! Body loads
        if (computeWext) then
            !! - Computing adjoint solution
            UA(:, :) = zero
            do iadj = 1, nadj
                do numCP = 1, NNODE
                    !!! --- INDEX I NOT DEFINED --- !!!
                    UA(:, iadj) = UA(:, iadj) + dRdXi(numCP, i) * UAelem(:, numCP, iadj)
                enddo
            enddo
            
            loadcount = 1
            do nl = 1, nb_load
                if (JDLTYPE(i)==101) then
                    !! - Centrifugal load
                    ! Gauss point location
                    pointGP(:) = zero
                    do numCP = 1, NNODE
                        pointGP(:) = pointGP(:) + R(numCP) * COORDS(:, numCP)
                    enddo
                    ! Distance to rotation axis
                    pointA(:) = zero
                    pointA(:MCRD) = &
                        &   load_additionalInfos(loadcount: loadcount+MCRD)
                    loadcount = loadcount + MCRD
                    pointB(:) = zero
                    pointB(:MCRD) = &
                        &   load_additionalInfos(loadcount: loadcount+MCRD)
                    loadcount = loadcount + MCRD
                    vectD(:) = pointB(:) - pointA(:)
                    vectD(:) = vectD(:) / SQRT(SUM(vectD(:) * vectD(:)))
                    vectAG(:) = pointGP(:) - pointA(:)
                    call dot(vectAG(:), vectD(:), scal)
                    vectR(:) = vectAG(:) - scal * vectD(:)
                    ! Save data
                    loadF(:) = DENSITY * ADLMAG(nl)**two * vectR(:)
                    scalFUA(:) = zero
                    vectDDUA(:, :) = zero
                    do iadj = 1, nadj
                        call dot(loadF(:), UA(:, iadj), scalFUA(iadj))
                        call dot(vectD(:), UA(:, iadj), coef1)
                        vectDDUA(:, iadj) = vectD(:) * coef1
                    enddo
                                       
                    do numCP = 1, NNODE
                        ! 1. Derivatives of the Jacobian determinant
                        dJdP(:) = zero
                        do idim = 1, MCRD
                            dJdP(idim) = detdXdXi * tracedXidX * dRdXi(numCP, idim)
                        enddo
                        do iadj = 1, nadj
                            gradWint_elem(iadj, :, numCP) = &
                                &   gradWint_elem(iadj, :, numCP) + &
                                &   scalFUA(iadj) * dJdP(:) * WtGauss * detdXidtildeXi
                        enddo
                        ! 2. Derivatives of the body force
                        dfdP_UA(:, :) = zero
                        do iadj = 1, nadj
                            dfdP_UA(:, iadj) = DENSITY * ADLMAG(nl)**two * & 
                                &   R(numCP) * (UA(:, iadj) - vectDDUA(:, iadj))
                        enddo
                        do iadj = 1, nadj
                            gradWext_elem(iadj, :, numCP) = &
                                &   gradWext_elem(iadj, :, numCP) + &
                                &   dfdP_UA(:, iadj) * dVol
                        enddo
                    enddo
                endif
            
            enddo  ! End loop on loads
            
        endif  ! test computeWext is True
         
    enddo  ! End of the loop on integration points
    
    !! Boundary loads
    if (computeWext) then
        kk = 0
        do i = 1, nb_load
            if ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND. &
                &   ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then
                !! Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
                call LectCle (JDLType(i), KNumFace, KTypeDload)
                call Gauss (NbPtInt, MCRD, GaussPdsCoord, KNumFace)
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!                           !!
                !!           TO DO           !!
                !!                           !!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
            endif
            kk = kk + load_target_nbelem(i)
        enddo
        
    endif ! test computeWext is True
      
end subroutine gradUELMAT1adjv2
