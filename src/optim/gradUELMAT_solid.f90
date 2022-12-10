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
!! with Yeti. If not, see <https://www.gnu.org/licenses/>

      
subroutine gradUELMAT1adj(activeDir,Uelem,UAelem,       &
                &   NADJ,NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,      &
                &   TENSOR,MATERIAL_PROPERTIES,DENSITY,PROPS,JPROPS,nb_load,    &
                &   indDLoad,load_target_nbelem,JDLType,ADLMAG,     &
                &   load_additionalInfos,nb_load_additionalInfos,   &
                &   computeWint,computeWext,gradWint_elem,gradWext_elem)
      
    use parameters
    use nurbspatch
      
    implicit none
      
    !! Input arguments :
    !! ---------------
    integer, intent(in) :: NADJ,NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES,    &
        &   DENSITY,PROPS
    dimension COORDS(3,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
    integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,     &
        &   nb_load,nb_load_additionalInfos
    double precision, intent(in) :: ADLMAG,load_additionalInfos
    dimension ADLMAG(nb_load),      &
        &   load_target_nbelem(nb_load),JDLType(nb_load),       &
        &   load_additionalInfos(nb_load_additionalInfos)
    dimension indDLoad(SUM(load_target_nbelem))
      
    double precision, intent(in) :: Uelem,UAelem
    integer, intent(in)          :: activeDir
    dimension Uelem(3,NNODE),UAelem(3,NNODE,NADJ),activeDir(3)
      
    logical, intent(in) :: computeWint,computeWext
      
    !! Output variables :
    !! ----------------
    double precision, intent(out) :: gradWint_elem, gradWext_elem
    dimension gradWint_elem(NADJ,3,NNODE),gradWext_elem(NADJ,3,NNODE)
      
      
    !! Local variables :
    !! ---------------
      
    !! For gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord,PtGauss
    dimension GaussPdsCoord(MCRD+1,NBINT),PtGauss(MCRD+1)
      
    !! For nurbs basis functions
    double precision :: XI, R, dRdxi, DetJac
    dimension R(NNODE), dRdxi(NNODE,3), XI(3)
      
    !! For curvilinear coordinate objects
    double precision :: AI,AIxAJ,AAI,AAE,det,Area
    dimension AI(3,3),AIxAJ(3,3),AAI(3,3),AAE(3,3)
      
    !! For material matrix
    integer :: voigt
    dimension voigt(6,2)
    double precision :: E,nu,lambda,mu, matH,coef
    dimension matH(2*MCRD,2*MCRD)
      
    !! For disp/strain/stress fields (state and adjoint)
    double precision :: dUdxi,dUAdxi,strain,stress,strainAdj,   &
        &   stressAdj,work,coef1,coef2
    dimension dUdxi(3,3),dUAdxi(3,3,NADJ),strain(2*MCRD),       &
        &   stress(2*MCRD),strainAdj(2*MCRD,nadj),stressAdj(2*MCRD,nadj),   &
        &   work(nadj)
      
    !! For derivatives
    double precision :: dAAIdP,dAAEdP,dJdP,dEAdP_S,dEdP_SA,dCdP,        &
        &   dSdP_EA
    dimension dAAIdP(3,3,3),dAAEdP(3,3,3),dJdP(3),      &
        &   dEAdP_S(3,nadj),dEdP_SA(3,nadj),dCdP(3),dSdP_EA(3,nadj)
      
    !! For loads
    integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
    dimension ind_load_loc(nb_load)
    double precision :: VectNorm_U,normV
    !! centrifugal load
    integer :: loadcount,nl, kload
    double precision :: pointGP, pointA,pointB,vectD,vectAG,vectR,  &
        &   scal,loadF,scalFUA,vectDDUA,dfdP_UA
    dimension pointGP(3),pointA(3),pointB(3),vectD(3),vectAG(3),    &
        &   vectR(3),loadF(3),scalFUA(nadj),vectDDUA(3,nadj),       &
        &   dfdP_UA(3,nadj)
    !! adjoint solution
    double precision :: UA
    dimension UA(3,NADJ)
      
    !! For loops
    integer ntens
    integer k1,k2,i,j,k,l,iA,kk,ll,ij,kl,cp, KTypeDload, KNumFace
    double precision :: temp,temp1,temp2
    integer :: idim, jdim
      
      
    !! Initialization
    !! --------------
    ntens   = 2*MCRD          ! size of stiffness tensor
    NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
    if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1
    
    !! Defining Gauss points coordinates and weights
    call Gauss(NbPtInt,MCRD,GaussPdsCoord,0)
      
    !! grad
    gradWint_elem(:,:,:)  = zero
    gradWext_elem(:,:,:)  = zero
      
    !! Material behaviour
    E  = MATERIAL_PROPERTIES(1)
    nu = MATERIAL_PROPERTIES(2)
    lambda = E*nu/(one+nu)/(one-two*nu)
    mu     = E/two/(one+nu)
    if (TENSOR == 'PSTRESS') lambda = two*lambda*mu/(lambda+two*mu)
      
    !! Computation
    !! -----------
      
    !! Loop on integration points
    do n = 1,NBINT
        
        !! PRELIMINARY QUANTITES
        !!  Computing NURBS basis functions and derivatives
        R(:)       = zero
        dRdxi(:,:) = zero
         
        XI(:)      = zero
        PtGauss(:) = GaussPdsCoord(2:,n)
        DetJac     = GaussPdsCoord( 1,n)
        do i = 1,MCRD
            XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i) &
                &   +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
            DetJac  = DetJac * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
        enddo
        call evalnurbs(XI,R,dRdxi)
         
        !! Computing Covariant basis vectors
        AI(:,:) = zero
        do i = 1,MCRD
            do cp = 1,NNODE
                AI(:,i) = AI(:,i) + dRdxi(cp,i)*COORDS(:,cp)
            enddo
        enddo
         
        call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
        if (MCRD==2) then
            call norm( AIxAJ(:,3),3, Area)
            AI(:,3) = AIxAJ(:,3)/Area
        else
            call dot(AIxAJ(:,3),AI(:,3),Area)
        endif
        call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
        call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
         
        if (computeWint) then   
            !! Computing metrics
            AAI(:,:) = zero
            AAE(:,:) = zero
            AAI(3,3) = one
            AAE(3,3) = one
            do i = 1,MCRD
                do j = 1,MCRD
                    call dot(AI(:,i), AI(:,j), AAI(i,j))
                enddo
            enddo
            call MatrixInv(AAE(:MCRD,:MCRD), AAI(:MCRD,:MCRD), det, MCRD)
         
         
            !! Computing material matrix
            voigt(:,1) = (/ 1,2,3,1,1,2 /)
            voigt(:,2) = (/ 1,2,3,2,3,3 /)
         
            matH(:,:) = zero
            do kl = 1,ntens
                k=voigt(kl,1); l=voigt(kl,2)
                do ij = 1,kl
                    i=voigt(ij,1); j=voigt(ij,2)

                    matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)      &
                        &   + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
               
                enddo
            enddo
            do kl = 2,ntens
                matH(kl:,kl-1) = matH(kl-1,kl:)
            enddo
         
            !! Computing disp and adjoint derivative
            dUdxi(:,:) = zero
            do i = 1,MCRD
                do cp = 1,NNODE
                    dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
                enddo
            enddo
         
            dUAdxi(:,:,:) = zero
            do iA = 1,nadj
                do i = 1,MCRD
                    do cp = 1,NNODE
                        dUAdxi(:,i,iA) =        &
                            &   dUAdxi(:,i,iA) + dRdxi(cp,i)*UAelem(:,cp,iA)
                    enddo
                enddo
            enddo
         
            !! Computing state strain and stress
            strain(:) = zero
            stress(:) = zero
            do ij = 1,ntens
                i=voigt(ij,1); j=voigt(ij,2)
                if (i==j) then
                    call dot(AI(:,i),dUdxi(:,i), strain(ij))
                else
                    call dot(AI(:,i),dUdxi(:,j), coef1)
                    call dot(AI(:,j),dUdxi(:,i), coef2)
                    strain(ij) = coef1 + coef2
                endif
            enddo
            call MulVect(matH,strain,stress,ntens,ntens)
         
            !! Computing adjoint strain and stress
            strainAdj(:,:) = zero
            stressAdj(:,:) = zero
            do iA = 1,nadj
                do ij = 1,ntens
                    i=voigt(ij,1); j=voigt(ij,2)
                    if (i==j) then
                        call dot(AI(:,i),dUAdxi(:,i,iA), strainAdj(ij,iA))
                    else
                        call dot(AI(:,i),dUAdxi(:,j,iA), coef1)
                        call dot(AI(:,j),dUAdxi(:,i,iA), coef2)
                        strainAdj(ij,iA) = coef1 + coef2
                    endif
                enddo 
                call MulVect(matH, strainAdj(:,iA), stressAdj(:,iA),        &
                    &   ntens, ntens)
            enddo
         
            !! Computing local work
            work(:)  = zero
            do ij = 1,ntens
                work(:) = work(:) + strainAdj(ij,:)*stress(ij)
            enddo
        endif ! test computeWint is True
         
        !! --
        !! Derivatives
        if (computeWint) then
            do cp = 1,NNODE
                !! 1. derivatives of the jacobian
                dJdP(:) = zero
                do i = 1,MCRD
                    dJdP(:) = dJdP(:) + AIxAJ(:,i)*dRdxi(cp,i)
                enddo

                do iA = 1,nadj
                    gradWint_elem(iA,:,cp)      &
                    &   = gradWint_elem(iA,:,cp) - work(iA)*dJdP(:)*detJac
                enddo
            
                !! 2. derivatives of the adjoint strain 
                !!    (with dble prod. by stress)
                dEAdP_S(:,:) = zero
                do iA = 1,nadj
                    do ij = 1,ntens
                        i=voigt(ij,1); j=voigt(ij,2)
                        if (i==j) then
                            dEAdP_S(:,iA) = dEAdP_S(:,iA)       &
                                &   + stress(ij)*dUAdxi(:,i,iA)*dRdxi(cp,i)
                        else
                            dEAdP_S(:,iA) = dEAdP_S(:,iA)       &
                                &   + stress(ij)*dUAdxi(:,i,iA)*dRdxi(cp,j) &
                                &   + stress(ij)*dUAdxi(:,j,iA)*dRdxi(cp,i)
                        endif
                    enddo
                    gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp)         &
                        &   - dEAdP_S(:,iA)*Area*detJac
                enddo

                !! 3. derivatives of the stress 
                !!    (with dble prod. by adjoint strain)
                dSdP_EA(:,:) = zero
                !! - derivatives of covariant metrics
                dAAIdP(:,:,:) = zero
                do j = 1,MCRD
                    do i = 1,MCRD
                        dAAIdP(:,i,j) = dRdxi(cp,i)*AI(:,j)     &
                            &   + AI(:,i)*dRdxi(cp,j)
                    enddo
                enddo
                !! - derivatives of contravariant metrics
                dAAEdP(:,:,:) = zero
                do j = 1,MCRD
                    do i = 1,MCRD
                        do l = 1,MCRD
                            do k = 1,MCRD
                                dAAEdP(:,i,j) = dAAEdP(:,i,j)       &
                                    &   - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                            enddo
                        enddo
                    enddo
                enddo
                !! - first subterm (material dot derivative strain)
                dEdP_SA(:,:) = zero
                do iA = 1,nadj
                    do ij = 1,ntens
                        i=voigt(ij,1); j=voigt(ij,2)
                        if (i==j) then
                            dEdP_SA(:,iA) = dEdP_SA(:,iA)       &
                                &   + stressAdj(ij,iA)*dUdxi(:,i)*dRdxi(cp,i)
                        else
                            dEdP_SA(:,iA) = dEdP_SA(:,iA)       &
                                &   + stressAdj(ij,iA)*dUdxi(:,i)*dRdxi(cp,j)       &
                                &   + stressAdj(ij,iA)*dUdxi(:,j)*dRdxi(cp,i)
                        endif
                    enddo
                enddo
                dSdP_EA(:,:) = dEdP_SA(:,:)
            
                !! - second subterm (derivative material tensor)
                do kl = 1,ntens
                    k=voigt(kl,1); l=voigt(kl,2)
               
                    do ij = 1,ntens
                        i=voigt(ij,1); j=voigt(ij,2)
                  
                        dCdP(:) =       &
                        &   lambda*dAAEdP(:,i,j)*AAE(k,l)       &
                        &   + lambda*AAE(i,j)*dAAEdP(:,k,l)     &
                        &   + mu*dAAEdP(:,i,k)*AAE(j,l)         &
                        &   + mu*AAE(i,k)*dAAEdP(:,j,l)         &
                        &   + mu*dAAEdP(:,i,l)*AAE(j,k)         &
                        &   + mu*AAE(i,l)*dAAEdP(:,j,k)
                  
                        do iA = 1,nadj
                            dSdP_EA(:,iA) = dSdP_EA(:,iA)       &
                            &   + dCdP(:)*strain(ij)*strainAdj(kl,iA)
                        enddo
                    enddo
                enddo

                do iA = 1,nadj
                    gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp)     &
                        &   - dSdP_EA(:,iA)*Area*detJac
                enddo
            enddo
        endif  ! test computeWint is True

        !! --
        !! Body loads
        if (computeWext) then
            !! - computing adjoint solution
            UA(:,:) = zero
            do iA = 1,nadj
                do cp = 1,NNODE
                    UA(:,iA)  = UA(:,iA) + R(cp)*UAelem(:,cp,iA)
                enddo
            enddo  
         
            loadcount = 1
            kload = 0
            do nl = 1,nb_load
                if ((JDLTYPE(nl)==101) .and.        &
                    &   ANY(indDLoad(kload+1:kload+load_target_nbelem(nl))==JELEM))     &
                    &   then

                    !! - centrifugal load
                    !! Gauss point location
                    pointGP(:) = zero
                    do cp = 1,NNODE
                        pointGP(:) = pointGP(:) + R(cp)*COORDS(:,cp)
                    enddo
                    !! Distance to rotation axis
                    pointA(:) = zero
                    pointA(:MCRD) =         &
                        &   load_additionalInfos(loadcount:loadcount+MCRD)
                    loadcount = loadcount+MCRD
                    pointB(:) = zero
                    pointB(:MCRD) =         &
                        &   load_additionalInfos(loadcount:loadcount+MCRD)
                    loadcount = loadcount+MCRD
               
                    vectD(:)  = pointB(:) - pointA(:)
                    vectD(:)  = vectD(:)/SQRT(SUM(vectD(:)*vectD(:)))
                    vectAG(:) = pointGP(:) - pointA(:)
                    call dot(vectAG(:),vectD(:),scal)
                    vectR(:)   = vectAG(:) - scal*vectD(:)
                    !! Save data
                    loadF(:) = DENSITY*ADLMAG(nl)**two*vectR(:)
                    scalFUA(:) = zero
                    vectDDUA(:,:) = zero
                    do iA = 1,nadj
                        call dot(loadF(:),UA(:,iA), scalFUA(iA))
                        call dot(vectD(:),UA(:,iA), coef1)
                        vectDDUA(:,iA) = vectD(:)*coef1
                    enddo
               
                    do cp = 1,NNODE
                        !! 1. derivatives of the jacobian
                        dJdP(:) = zero
                        do i = 1,MCRD
                            dJdP(:) = dJdP(:) + AIxAJ(:,i)*dRdxi(cp,i)
                        enddo

                        do iA = 1,nadj
                            gradWext_elem(iA,:,cp) = gradWext_elem(iA,:,cp)     &
                                &   + scalFUA(iA)*dJdP(:)*detJac
                        enddo

                        !! 2. derivatives of the body force
                        dfdP_UA(:,:) = zero
                        do iA = 1, nadj
                            do idim = 1, mcrd
                                do jdim = 1, mcrd
                                    if (idim == jdim) then
                                        dfdP_UA(jdim,iA) = dfdP_UA(jdim,iA) +   &
                                            &   UA(idim, iA)                    &
                                            &   *(one - vectD(idim)*vectD(jdim))
                                    else
                                        dfdP_UA(jdim,iA) = dfdP_UA(jdim,iA) -   &
                                            &   UA(idim, iA)                    &
                                            &   *(vectD(idim)*vectD(jdim))  
                                    endif
                                enddo
                            enddo
                        enddo
                        dfdP_UA(:,:) = dfdP_UA(:,:)     &
                            &   *DENSITY*ADLMAG(nl)**two*R(cp)

                        do iA = 1,nadj
                            gradWext_elem(iA,:,cp) = gradWext_elem(iA,:,cp)         &
                                &   + dfdP_UA(:,iA)*Area*detJac
                        enddo
                    enddo
                endif
                kload = kload + load_target_nbelem(nl)
            enddo
        endif  ! test computeWext is True
         
    enddo
    !! End of the loop on integration points

      
    if (computeWext) then
        kk = 0
        do i = 1,nb_load
            if ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND.       &
            &   ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then
                !! Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
                call LectCle (JDLType(i),KNumFace,KTypeDload)
                call Gauss (NbPtInt,MCRD,GaussPdsCoord,KNumFace)
            
                !! to do
            
            endif
            kk = kk + load_target_nbelem(i)
        enddo
    endif   !! test computeWext is True
      
end subroutine gradUELMAT1adj


