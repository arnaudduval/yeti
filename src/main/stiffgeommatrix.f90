!! Compute elementary matrix at Gauss point, for 2D/3D solids
!! ---
!! Only upper triangular part is built. Lower part is kpt at zero
!! symmetry is taken into account during sum at Gauss points, at the
!! end of elementary built (in UELTANMAT.f90)


!! Storage without building the matrix : return 3x3 matrices
!! for each couple of control points (upper triangle part only)
subroutine stiffgeommatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, Uelem,       &
     &     stiffgeom)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL
    double precision, intent(in) :: ddsdde, dRdx
    dimension ddsdde(ntens, ntens), dRdx(MCRD, NNODE)
    double precision, intent(in) ::  Uelem
    dimension Uelem(MCRD,NNODE)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: stiffgeom
    dimension stiffgeom(MCRD, MCRD, NNODE*(NNODE+1)/2)

    !! Local variables
    !! ---------------
    double precision :: BoJ, dNjdx, dNidx, stiffLoc, SBoJ, Bl, Bnl, Ui, Uj, PiJ
    dimension BoJ(ntens, MCRD), dNjdx(MCRD), Bl(ntens, MCRD), Bnl(ntens, MCRD)
    dimension dNidx(MCRD), stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD), Ui(MCRD), Uj(MCRD), piJ(ntens)

    integer :: nodj, nodi, jdim, idim, i, j, k, count
    double precision :: res, voigt

    stiffgeom(:,:,:) = zero
    BoJ(:,:)     = zero
    PiJ(:)       = zero
    Bl(:,:)      = zero
    voigt        = sqrt(2.0D0)/2


     ! Loop on control points
    do nodj = 1, NNODE
        dNjdx(:) = dRdx(:, nodj)
        Uj(:) = Uelem(:,nodj)

        ! Compute BlJ
        do i = 1,MCRD
            Bl(i, i) = dNjdx(i)
        enddo
        Bl(4,1) = voigt*dNjdx(2)
        Bl(4,2) = voigt*dNjdx(1)
        Bl(5,1) = voigt*dNjdx(3)
        Bl(5,3) = voigt*dNjdx(1)
        Bl(6,2) = voigt*dNjdx(3)
        Bl(6,3) = voigt*dNjdx(2)

        BoJ(:,:)  = Bl(:,:)
        
        do nodi=1,NNODE
            dNidx(:) = dRdx(:, nodi)
            Ui(:) = Uelem(:,nodi)

            !! Compute Bnl
            Bnl(1,1) = half*Ui(1)*dNidx(1)*dNjdx(1)
            Bnl(1,2) = half*Ui(2)*dNidx(1)*dNjdx(1)
            Bnl(1,3) = half*Ui(3)*dNidx(1)*dNjdx(1)

            Bnl(2,1) = half*Ui(1)*dNidx(2)*dNjdx(2)
            Bnl(2,2) = half*Ui(2)*dNidx(2)*dNjdx(2)
            Bnl(2,3) = half*Ui(3)*dNidx(2)*dNjdx(2)

            Bnl(3,1) = half*Ui(1)*dNidx(3)*dNjdx(3)
            Bnl(3,2) = half*Ui(2)*dNidx(3)*dNjdx(3)
            Bnl(3,3) = half*Ui(3)*dNidx(3)*dNjdx(3)

            Bnl(4,1) = voigt*half*Ui(1)*(dNidx(1)*dNjdx(2)+dNidx(2)*dNjdx(1))
            Bnl(4,2) = voigt*half*Ui(2)*(dNidx(1)*dNjdx(2)+dNidx(2)*dNjdx(1))
            Bnl(4,3) = voigt*half*Ui(3)*(dNidx(1)*dNjdx(2)+dNidx(2)*dNjdx(1))

            Bnl(5,1) = voigt*half*Ui(1)*(dNidx(1)*dNjdx(3)+dNidx(3)*dNjdx(1))
            Bnl(5,2) = voigt*half*Ui(2)*(dNidx(1)*dNjdx(3)+dNidx(3)*dNjdx(1))
            Bnl(5,3) = voigt*half*Ui(3)*(dNidx(1)*dNjdx(3)+dNidx(3)*dNjdx(1))

            Bnl(6,1) = voigt*half*Ui(1)*(dNidx(2)*dNjdx(3)+dNidx(3)*dNjdx(2))
            Bnl(6,2) = voigt*half*Ui(2)*(dNidx(2)*dNjdx(3)+dNidx(3)*dNjdx(2))
            Bnl(6,3) = voigt*half*Ui(3)*(dNidx(2)*dNjdx(3)+dNidx(3)*dNjdx(2))

            BoJ(:,:)  = BoJ(:,:) + Bnl(:,:)
        end do

        !! Compute product ddsdde*BoJ
        SBoJ(:, :) = zero
        do j=1,MCRD
            do i=1,MCRD
                res = 0
                do k=1,MCRD
                    res = res + ddsdde(i,k)*BoJ(k,j)
                end do
                sBoJ(i,j) = res
            enddo
        enddo
        SBoJ(4, 1) = ddsdde(4, 4)*BoJ(4, 1)
        SBoJ(4, 2) = ddsdde(4, 4)*BoJ(4, 2)
        SBoJ(4, 3) = ddsdde(4, 4)*BoJ(4, 3)
        SBoJ(5, 1) = ddsdde(5, 5)*BoJ(5, 1)
        SBoJ(5, 2) = ddsdde(5, 5)*BoJ(5, 2)
        SBoJ(5, 3) = ddsdde(5, 5)*BoJ(5, 3)
        SBoJ(6, 1) = ddsdde(6, 6)*BoJ(6, 1)
        SBoJ(6, 2) = ddsdde(6, 6)*BoJ(6, 2)
        SBoJ(6, 3) = ddsdde(6, 6)*BoJ(6, 3)

        PiJ(1) = PiJ(1) + SBoJ(1,1)*Uj(1) + SBoJ(1,2)*Uj(2) + SBoJ(1,3)*Uj(3)
        PiJ(2) = PiJ(2) + SBoJ(2,1)*Uj(1) + SBoJ(2,2)*Uj(2) + SBoJ(2,3)*Uj(3)
        PiJ(3) = PiJ(3) + SBoJ(3,1)*Uj(1) + SBoJ(3,2)*Uj(2) + SBoJ(3,3)*Uj(3)
        PiJ(4) = PiJ(4) + SBoJ(4,1)*Uj(1) + SBoJ(4,2)*Uj(2) + SBoJ(4,3)*Uj(3)
        PiJ(5) = PiJ(5) + SBoJ(5,1)*Uj(1) + SBoJ(5,2)*Uj(2) + SBoJ(5,3)*Uj(3)
        PiJ(6) = PiJ(6) + SBoJ(6,1)*Uj(1) + SBoJ(6,2)*Uj(2) + SBoJ(6,3)*Uj(3)
    enddo

    count = 1
    !! Loop on control points
    do nodj = 1, NNODE
        dNjdx(:) = dRdx(:, nodj)

        !! Second loop on control points
        do nodi = 1, nodj
            dNidx(:) = dRdx(:, nodi)
        
            !! Computation and Assembly of GEOMATRIX
            do k =1,MCRD
                stiffLoc(k,k) = PiJ(1)*dNidx(1)*dNjdx(1) + PiJ(2)*dNidx(2)*dNjdx(2) + PiJ(3)*dNidx(3)*dNjdx(3) &
                              + PiJ(4)*voigt*(dNidx(1)*dNjdx(2)+dNidx(2)*dNjdx(1)) &
                              + PiJ(5)*voigt*(dNidx(1)*dNjdx(3)+dNidx(3)*dNjdx(1)) &
                              + PiJ(6)*voigt*(dNidx(2)*dNjdx(3)+dNidx(3)*dNjdx(2))
            enddo

            !! Assembly
            stiffgeom(:, :, count) = stiffLoc(:, :)

            count = count + 1
        enddo
    enddo

end subroutine stiffgeommatrix_byCP



