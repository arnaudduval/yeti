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


subroutine UELMAT_byCP(NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS,            &
        &   TENSOR, MATERIAL_PROPERTIES, DENSITY, nb_load, indDLoad,    &
        &   load_target_nbelem, JDLType, ADLMAG, load_additionalInfos, &
        &   len_load_additionalInfos, nb_load_additionalInfos, &
        &   n_dist_elem, nb_n_dist, RHS, AMATRX)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, &
        &   DENSITY,    &
        &   n_dist_elem
    dimension COORDS(MCRD, NNODE), MATERIAL_PROPERTIES(2),    &
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
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD, MCRD, NNODE * (NNODE + 1) / 2)


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
    AMATRX(:,:,:) = zero

    !! Material behaviour
    call material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde)

    !! Loop on integration points
    do n = 1, NBINT
        !! Compute NURBS basis functions and derivatives
        call shap(dRdx, R, DetJac, COORDS, GaussPdsCoord(2:, n), MCRD)

        !! Compute stiffness matrix
        call stiffmatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx,  &
            &        stiff)

        !! Assemble AMATRIX
        dvol = GaussPdsCoord(1, n) * ABS(detJac)
        AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:) * dvol

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
            else if (JDLTYPE(i_load) == 102) then
                if (ANY(indDLoad(kload + 1:kload + load_target_nbelem(i_load)) == JELEM)) then
                    !! Volumic force
                    vectR(:) = load_additionalInfos(load_addinfos_count:    &
                        &                    load_addinfos_count + MCRD)
                    !! Update load vector
                    kk = 0
                    do numCP = 1, NNODE
                        do j = 1, MCRD
                            kk = kk + 1
                            RHS(kk) = RHS(kk) + ADLMAG(i_load) * vectR(j) * R(numCP) * dvol
                        end do
                    end do
                endif
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

end subroutine UELMAT_byCP



!! ***************************************************************************
!! *       THE FOLLOWING ROUTINES ARE OBSOLETE AND NOT USED ANYMORE          *
!! ***************************************************************************

!
!
!
!
! SUBROUTINE UELMAT(NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,TENSOR,
! 1     MATERIAL_PROPERTIES,nb_load,indDLoad,load_target_nbelem,
! 2     JDLType,ADLMAG,RHS,AMATRX)
!
! use parameters
!
! Implicit None
!
! c     Input arguments :
! c     ---------------
! Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT
! Character(len=*), intent(in) :: TENSOR
! Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES
! dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2)
!
! Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
!
! Double precision, intent(in) :: ADLMAG
! dimension ADLMAG(nb_load),
! &     load_target_nbelem(nb_load),JDLType(nb_load)
! dimension indDLoad(SUM(load_target_nbelem))
!
! c     Output variables :
! c     ----------------
! Double precision, intent(out) :: RHS, AMATRX
! dimension RHS(NDOFEL), AMATRX(NDOFEL,NDOFEL)
!
!
! c     Local variables :
! c     ---------------
!
! !     For gauss points
! Integer :: NbPtInt, n
! Double precision :: GaussPdsCoord
! dimension GaussPdsCoord(MCRD+1,NBINT)
!
! !     For nurbs basis functions
! Double precision :: R, dRdx, DetJac
! dimension R(NNODE),dRdx(MCRD,NNODE)
!
! !     For material behaviour
! Double precision :: ddsdde
! dimension ddsdde(2*MCRD,2*MCRD)
!
! !     For stiffness matrix
! Integer :: k1,k2,ntens
! Double precision :: stiff,dvol
! dimension stiff(NDOFEL,NDOFEL)
!
! !     For load vector
! Integer :: i,kk,KNumFace,KTypeDload,numPC,numI
! Double precision :: FbL,VectNorm, y
! dimension FbL(NDOFEL),VectNorm(MCRD)
!
!
!
! C     ------------------------------------------------------------------
!
! c     Initialization :
! c     --------------
! ntens   = 2*MCRD          ! size of stiffness tensor
! NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
! if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1
!
! c     Defining Gauss points coordinates and weights
! call Gauss(NbPtInt,MCRD,GaussPdsCoord,0)
!
! c     Stiffness matrix and load vector initialized to zero
! RHS(:) = zero
! AMATRX(:,:) = zero
!
! c     Material behaviour
! call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde)
! c
! c     ..................................................................
! c
! C     Computation :
! c     -----------
!
!
! c     Loop on integration points
! Do n = 1,NBINT
! c     Computing NURBS basis functions and derivatives
! call shap (dRdx,R,DetJac,COORDS,GaussPdsCoord(2:,n),MCRD)
!
! c     Computing stiffness matrix
! call stiffmatrix(ntens,NNODE,MCRD,NDOFEL,ddsdde,dRdx,stiff)
!
! c     Assembling AMATRIX
! dvol = GaussPdsCoord(1,n)*detJac
! Do k2 = 1,NDOFEL
! Do k1 = 1,k2 !1,NDOFEL
! AMATRX(k1,k2) = AMATRX(k1,k2) + stiff(k1,k2)*dvol
! Enddo
! Enddo
! Enddo
!
! c     Symmetry
! Do k2 = 1,NDOFEL-1
! Do k1 = k2+1,NDOFEL
! AMATRX(k1,k2) = AMATRX(k2,k1)
! Enddo
! Enddo
!
! c     End of the loop on integration points on main surf
! c
! c     ..................................................................
! c
! c     Loop for load : find boundary loads
! kk = 0
! Do i = 1,nb_load
! If ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND.
! &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then
! c     Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
! call LectCle (JDLType(i),KNumFace,KTypeDload)
! call Gauss (NbPtInt,MCRD,GaussPdsCoord,KNumFace)
!
! FbL(:) = zero
! Do n = 1,NbPtInt**(MCRD-1)
!
! If (KTypeDload==7) then
!
! call shapPress(R,VectNorm,DetJac,COORDS,
! &              GaussPdsCoord(2:,n),MCRD,KNumFace,2,
! &              VectNorm)
! y = 0.d0
! Do numPC = 1,NNODE
! y = y + R(numPC)*COORDS(2,numPC)
! Enddo
! dvol = -3.d0/40.d0*ADLMAG(i)*(1.d0-y*y/100.d0)
! &              *GaussPdsCoord(1,n)*DetJac
! Elseif (KTypeDload==6) then
! ! snow load
! call shapPress(R,VectNorm,DetJac,COORDS,
! &                 GaussPdsCoord(2:,n),MCRD,KNumFace,0,
! &                 VectNorm)
! VectNorm(1:2) = zero
! dvol = ADLMAG(i)*GaussPdsCoord(1,n)*DetJac
!
! Else
!
! call shapPress(R,VectNorm,DetJac,COORDS,
! &                 GaussPdsCoord(2:,n),MCRD,KNumFace,KTypeDload,
! &                 VectNorm)
! dvol = ADLMAG(i)*GaussPdsCoord(1,n)*DetJac
! Endif
!
!
!
! Do numPC = 1,NNODE
! numI = (NumPC-1)*MCRD
! Do k2 = 1,MCRD
! numI = numI + 1
! FbL(numI)=FbL(numI) + R(numPC)*VectNorm(k2)*dvol
! !FbL(numI)=FbL(numI) + VectNorm(k2)*dvol
! Enddo
! Enddo
! Enddo
!
! c     Assembling RHS
! Do k1 = 1,NDOFEL
! RHS(k1) = RHS(k1) + FbL(k1)
! Enddo
! Endif
! kk = kk + load_target_nbelem(i)
! Enddo
!
! End SUBROUTINE UELMAT
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
! SUBROUTINE UELMAT_wCentri(NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,
! 1     TENSOR,MATERIAL_PROPERTIES,nb_load,indDLoad,
! 2     load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX)
!
! use parameters
!
! Implicit None
!
! c     Input arguments :
! c     ---------------
! Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT
! Character(len=*), intent(in) :: TENSOR
! Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES
! dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2)
!
! Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
!
! Double precision, intent(in) :: ADLMAG
! dimension ADLMAG(nb_load),
! &     load_target_nbelem(nb_load),JDLType(nb_load)
! dimension indDLoad(SUM(load_target_nbelem))
!
! c     Output variables :
! c     ----------------
! Double precision, intent(out) :: RHS, AMATRX
! dimension RHS(NDOFEL), AMATRX(MCRD,MCRD,NNODE*(NNODE+1)/2)
!
!
! c     Local variables :
! c     ---------------
!
! !     For gauss points
! Integer :: NbPtInt, n
! Double precision :: GaussPdsCoord
! dimension GaussPdsCoord(MCRD+1,NBINT)
!
! !     For nurbs basis functions
! Double precision :: R, dRdx, DetJac
! dimension R(NNODE),dRdx(MCRD,NNODE)
!
! !     For material behaviour
! Double precision :: ddsdde
! dimension ddsdde(2*MCRD,2*MCRD)
!
! !     For stiffness matrix
! Integer :: k1,k2,ntens
! Double precision :: stiff,dvol
! dimension stiff( MCRD,MCRD,NNODE*(NNODE+1)/2 )
!
! !     For load vector
! Integer :: i,kk,KNumFace,KTypeDload,numPC,numI
! Double precision :: FbL,VectNorm, y
! dimension FbL(NDOFEL),VectNorm(MCRD)
!
! !      For centrifugal load
! Integer :: cp
! Double precision :: pointA,pointB,pointGP,vectU,normU,vectV,normV,
! &     vectW,vectD,distR,omegaRot,density
! dimension pointA(3),pointB(3),pointGP(3),vectU(3),vectV(3),
! &     vectW(3),vectD(3)
!
!
!
! C     ------------------------------------------------------------------
!
! c     Initialization :
! c     --------------
! ntens   = 2*MCRD          ! size of stiffness tensor
! NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
! if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1
!
! pointA(:) = (/ zero,zero,zero /)
! pointB(:) = (/  one,zero,zero /)
! omegaRot  = 1000.d0
! density   = 2750.d0
! vectU(:)  = pointB(:) - pointA(:)
! call norm(vectU(:),3,normU)
! vectU(:) = vectU(:)/normU
!
!
! c     Defining Gauss points coordinates and weights
! call Gauss(NbPtInt,MCRD,GaussPdsCoord,0)
!
! c     Stiffness matrix and load vector initialized to zero
! RHS(:)        = zero
! AMATRX(:,:,:) = zero
!
! c     Material behaviour
! call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde)
! c
! c     ..................................................................
! c
! C     Computation :
! c     -----------
!
! c     Loop on integration points
! Do n = 1,NBINT
! c     Computing NURBS basis functions and derivatives
! call shap (dRdx,R,DetJac,COORDS,GaussPdsCoord(2:,n),MCRD)
!
! c     Computing stiffness matrix
! call stiffmatrix_byCP(ntens,NNODE,MCRD,NDOFEL,ddsdde,dRdx,
! &        stiff)
!
! c     Assembling AMATRIX
! dvol = GaussPdsCoord(1,n)*detJac
! AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:)*dvol
!
! c     Centrifugal load
! ! Gauss point location
! pointGP(:) = zero
! Do cp = 1,NNODE
! pointGP(:) = pointGP(:) + R(cp)*COORDS(:,cp)
! Enddo
!
! ! Distance to rotation axis
! vectD(:) = zero
! call CROSS(pointGP(:)-pointA(:),vectU(:),vectD(:))
! call norm(vectD(:),3,distR)
! vectW(:) = vectD(:)/distR
! call CROSS(vectW(:),vectU(:),vectV(:))
!
! ! update load vector
! kk = 0
! Do cp = 1,NNODE
! Do i = 1,MCRD
! kk = kk+1
! RHS(kk) = RHS(kk)
! &              - density*omegaRot**two*distR*vectV(i)*R(cp)*dvol
! Enddo
! Enddo
!
!
! Enddo
!
! c     End of the loop on integration points on main surf
! c
! c     ..................................................................
! c
! c     Loop for load : find boundary loads
! kk = 0
! Do i = 1,nb_load
! If ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND.
! &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then
! c     Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
! call LectCle (JDLType(i),KNumFace,KTypeDload)
! call Gauss (NbPtInt,MCRD,GaussPdsCoord,KNumFace)
!
! FbL(:) = zero
! Do n = 1,NbPtInt**(MCRD-1)
!
! call shapPress(R,VectNorm,DetJac,COORDS,
! &              GaussPdsCoord(2:,n),MCRD,KNumFace,KTypeDload,
! &              VectNorm)
! dvol = ADLMAG(i)*GaussPdsCoord(1,n)*DetJac
!
!
! Do numPC = 1,NNODE
! numI = (NumPC-1)*MCRD
! Do k2 = 1,MCRD
! numI = numI + 1
! FbL(numI)=FbL(numI) + R(numPC)*VectNorm(k2)*dvol
! !FbL(numI)=FbL(numI) + VectNorm(k2)*dvol
! Enddo
! Enddo
! Enddo
!
! c     Assembling RHS
! Do k1 = 1,NDOFEL
! RHS(k1) = RHS(k1) + FbL(k1)
! Enddo
! Endif
! kk = kk + load_target_nbelem(i)
! Enddo
!
! End SUBROUTINE UELMAT_wCentri
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
! c     ------------------------------------------------------------------
! c
! c     Formulation solide en utilisant les coordonnees curvilignes
! c
!
! SUBROUTINE UELMAT_curvilinear(NDOFEL,MCRD,NNODE,JELEM,NBINT,
! 1     COORDS,TENSOR,MATERIAL_PROPERTIES,nb_load,indDLoad,
! 2     load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX)
!
! use parameters
! use nurbspatch
!
! Implicit None
!
! c     Input arguments :
! c     ---------------
! Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT
! Character(len=*), intent(in) :: TENSOR
! Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES
! dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2)
!
! Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
!
! Double precision, intent(in) :: ADLMAG
! dimension ADLMAG(nb_load),
! &     load_target_nbelem(nb_load),JDLType(nb_load)
! dimension indDLoad(SUM(load_target_nbelem))
!
! c     Output variables :
! c     ----------------
! Double precision, intent(out) :: RHS, AMATRX
! dimension RHS(NDOFEL), AMATRX(MCRD,MCRD,NNODE*(NNODE+1)/2)
!
!
! c     Local variables :
! c     ---------------
!
! !     For gauss points
! Integer :: NbPtInt, n
! Double precision :: GaussPdsCoord,PtGauss
! dimension GaussPdsCoord(MCRD+1,NBINT),PtGauss(MCRD+1)
!
! !     For nurbs basis functions
! Double precision :: XI, R, dRdxi, DetJac
! dimension R(NNODE),dRdxi(NNODE,3),XI(3)
!
! !     For curvilinear quantities
! Double precision :: AI,AAI,AAE,det,scal
! dimension AI(3,MCRD),AAI(MCRD,MCRD),AAE(MCRD,MCRD)
!
! !     For material behaviour
! Double precision :: ddsdde,coef,coef1,coef2,E,nu
! dimension ddsdde(2*MCRD,2*MCRD)
!
! !     For stiffness matrix
! Integer :: i,j, k1,k2,ntens,cp
! Double precision :: stiff,dvol,vectV,normV
! dimension stiff( MCRD,MCRD,NNODE*(NNODE+1)/2 ),vectV(3)
!
! !     For load vector
! Integer :: kk,KNumFace,KTypeDload,numPC,numI
! Double precision :: FbL,VectNorm, y
! dimension FbL(NDOFEL),VectNorm(MCRD)
!
!
!
! C     ------------------------------------------------------------------
!
! c     Initialization :
! c     --------------
! ntens   = 2*MCRD          ! size of stiffness tensor
! NbPtInt = int( NBINT**(1.0/float(MCRD)) ) ! nb gauss pts per dir.
! if (NbPtInt**MCRD<NBINT) NbPtInt = NbPtInt + 1
!
! c     Defining Gauss points coordinates and weights
! call Gauss(NbPtInt,MCRD,GaussPdsCoord,0)
!
! c     Stiffness matrix and load vector initialized to zero
! RHS(:)        = zero
! AMATRX(:,:,:) = zero
!
! c     Material behaviour
! !call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde)
! E  = MATERIAL_PROPERTIES(1)
! nu = MATERIAL_PROPERTIES(2)
! coef  = E/(one-nu*nu)
! coef1 = E/two/(one+nu)
! coef2 = two*nu/(one-two*nu)
!
! c
! c     ..................................................................
! c
! C     Computation :
! c     -----------
!
! c     Loop on integration points
! Do n = 1,NBINT
! c     Computing NURBS basis functions and derivatives
! R(:)       = zero
! dRdxi(:,:) = zero
!
! XI(:)      = zero
! PtGauss(:) = GaussPdsCoord(2:,n)
! DetJac     = GaussPdsCoord( 1,n)
! Do i = 1,MCRD
! XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
! &            +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
! DetJac  = DetJac * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
! End do
! call evalnurbs(XI,R,dRdxi)
!
! c     Computing Covariant basis vectors
! AI(:,:) = zero
! Do i = 1,MCRD
! Do cp = 1,NNODE
! AI(:MCRD,i) = AI(:MCRD,i) + dRdxi(cp,i)*COORDS(:,cp)
! Enddo
! Enddo
!
! c     Computing material stiffness tensor
! AAI(:,:) = zero
! Do i = 1,MCRD
! Do j = 1,MCRD
! call dot(AI(:,i), AI(:,j), scal)
! AAI(i,j) = scal
! Enddo
! Enddo
! call MatrixInv(AAE, AAI, det, MCRD)
!
!
!
! ddsdde(:,:) = zero
! if (TENSOR == 'PSTRESS') then
!
! ddsdde(1,1) = AAE(1,1)*AAE(1,1)
! ddsdde(2,2) = AAE(2,2)*AAE(2,2)
! ddsdde(4,4) = 0.5d0*
! &        ((one-nu)*AAE(1,1)*AAE(2,2) + (one+nu)*AAE(1,2)*AAE(1,2))
! ddsdde(1,2) = nu*AAE(1,1)*AAE(2,2) + (one-nu)*AAE(1,2)*AAE(1,2)
! ddsdde(1,4) = AAE(1,1)*AAE(1,2)
! ddsdde(2,4) = AAE(2,2)*AAE(1,2)
! ddsdde(2,1) = ddsdde(1,2)
! ddsdde(4,1) = ddsdde(1,4)
! ddsdde(4,2) = ddsdde(2,4)
! ddsdde(:,:) = coef*ddsdde(:,:)
!
! else
!
! Do i = 1,MCRD
! ddsdde(i,i) = two*AAE(i,i)*AAE(i,i) + coef2*AAE(i,i)*AAE(i,i)
! Do j = i+1,MCRD
! ddsdde(i,j) = two*AAE(i,j)*AAE(i,j) + coef2*AAE(i,i)*AAE(j,j)
! Enddo
!
! ddsdde(i,4) = two*AAE(i,1)*AAE(i,2) + coef2*AAE(i,i)*AAE(1,2)
! if (MCRD==3) then
! ddsdde(i,5) = two*AAE(i,1)*AAE(i,3) + coef2*AAE(i,i)*AAE(1,3)
! ddsdde(i,6) = two*AAE(i,2)*AAE(i,3) + coef2*AAE(i,i)*AAE(2,3)
! Endif
! Enddo
!
! ddsdde(4,4) = AAE(1,1)*AAE(2,2) + AAE(1,2)*AAE(2,1)
! &               + coef2*AAE(1,2)*AAE(1,2)
! if (MCRD==3) then
! ddsdde(5,5) = AAE(1,1)*AAE(3,3) + AAE(1,3)*AAE(3,1)
! &               + coef2*AAE(1,3)*AAE(1,3)
! ddsdde(6,6) = AAE(2,2)*AAE(3,3) + AAE(2,3)*AAE(3,2)
! &               + coef2*AAE(2,3)*AAE(2,3)
!
! ddsdde(4,5) = AAE(1,1)*AAE(2,3) + AAE(1,3)*AAE(2,1)
! &               + coef2*AAE(1,2)*AAE(1,3)
! ddsdde(4,6) = AAE(1,2)*AAE(2,3) + AAE(1,3)*AAE(2,2)
! &               + coef2*AAE(1,2)*AAE(2,3)
! ddsdde(5,6) = AAE(1,2)*AAE(3,3) + AAE(1,3)*AAE(3,2)
! &               + coef2*AAE(1,3)*AAE(2,3)
! Endif
!
! ! - symmetry
! Do i = 2,ntens
! Do j = 1,i-1
! ddsdde(i,j) = ddsdde(j,i)
! Enddo
! Enddo
! ddsdde(:,:) = coef1*ddsdde(:,:)
!
! Endif
!
!
! c     Computing stiffness matrix
! call stiffmatrix_curv(ntens,NNODE,MCRD,NDOFEL,ddsdde,AI,
! &        dRdxi(:,:MCRD),stiff)
!
! c     Assembling AMATRIX
! call cross(AI(:,1),AI(:,2),vectV(:))
! if (MCRD==2) then
! call norm(vectV(:),3, normV)
! else
! call dot(  AI(:,3),vectV(:),normV)
! endif
! dvol   = normV*detJac
! AMATRX(:,:,:) = AMATRX(:,:,:) + stiff(:,:,:)*dvol
!
! Enddo
!
! c     End of the loop on integration points on main surf
! c
! c     ..................................................................
! c
! c     Loop for load : find boundary loads
! kk = 0
! Do i = 1,nb_load
! If ((JDLTYPE(i)>9 .AND. JDLTYPE(i)<100) .AND.
! &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i))==JELEM)) then
! c     Defining Gauss points coordinates and weights on surf(3D)/edge(2D)
! call LectCle (JDLType(i),KNumFace,KTypeDload)
! call Gauss (NbPtInt,MCRD,GaussPdsCoord,KNumFace)
!
! FbL(:) = zero
! Do n = 1,NbPtInt**(MCRD-1)
!
! call shapPress(R,VectNorm,DetJac,COORDS,
! &              GaussPdsCoord(2:,n),MCRD,KNumFace,KTypeDload,
! &              VectNorm)
! dvol = ADLMAG(i)*GaussPdsCoord(1,n)*DetJac
!
!
! Do numPC = 1,NNODE
! numI = (NumPC-1)*MCRD
! Do k2 = 1,MCRD
! numI = numI + 1
! FbL(numI)=FbL(numI) + R(numPC)*VectNorm(k2)*dvol
! !FbL(numI)=FbL(numI) + VectNorm(k2)*dvol
! Enddo
! Enddo
! Enddo
!
! c     Assembling RHS
! Do k1 = 1,NDOFEL
! RHS(k1) = RHS(k1) + FbL(k1)
! Enddo
! Endif
! kk = kk + load_target_nbelem(i)
! Enddo
!
! !print*,''
! !Do i = 1,ntens
! !   print*,ddsdde(i,:ntens)
! !Enddo
!
! End SUBROUTINE UELMAT_curvilinear
!
!


