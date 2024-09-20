!! Compute elementary matrix at Gauss point, for 2D/3D solids
!! ---
!! Only upper triangular part is built. Lower part is kpt at zero
!! symmetry is taken into account during sum at Gauss points, at the
!! end of elementary built (in UELTANMAT.f90)


!! Storage without building the matrix : return 3x3 matrices
!! for each couple of control points (upper triangle part only)
subroutine geomvector_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, Uelem,       &
     &     Fgeom)

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
    double precision, intent(out) :: Fgeom
    dimension Fgeom(MCRD, NNODE)

    !! Local variables
    !! ---------------
    double precision :: BoJ, BoIt, dNjdx, dNidx, dNkdx, SBoJ, BlJ, BnlJ, BlIt, BnlIt, Uk, Uj, SoJ
    double precision :: stiffLoc
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), BlIt(MCRD, ntens), BnlIt(MCRD, ntens), BlJ(ntens, MCRD), BnlJ(ntens, MCRD)
    dimension stiffLoc(MCRD, MCRD)
    dimension dNidx(MCRD), dNjdx(MCRD), dNkdx(MCRD), SBoJ(ntens, MCRD), Uk(MCRD), Uj(MCRD), SoJ(ntens)
    integer :: nodj, nodi, nodk, jdim, idim, i, j, k, count
    double precision :: res, voigt


    Fgeom(:,:) = zero
    BlJ(:,:)   = zero
    BlIt(:,:)  = zero
    SoJ(:)  = zero
    voigt  = sqrt(2.0D0)/2

    ! Loop on control points
    do nodj = 1, NNODE
        dNjdx(:) = dRdx(:, nodj)
        Uj(:) = Uelem(:,nodj)

        ! Compute BlJ
        do i = 1,MCRD
            BlJ(i, i) = dNjdx(i)
        enddo
        BlJ(4,1) = voigt*dNjdx(2)
        BlJ(4,2) = voigt*dNjdx(1)
        BlJ(5,1) = voigt*dNjdx(3)
        BlJ(5,3) = voigt*dNjdx(1)
        BlJ(6,2) = voigt*dNjdx(3)
        BlJ(6,3) = voigt*dNjdx(2)

        BoJ(:,:)  = BlJ(:,:)
        
        do nodk=1,NNODE
            dNkdx(:) = dRdx(:, nodk)
            Uk(:) = Uelem(:,nodk)

            ! Compute BnlJ
            BnlJ(1,1) = half*Uk(1)*dNkdx(1)*dNjdx(1)
            BnlJ(1,2) = half*Uk(2)*dNkdx(1)*dNjdx(1)
            BnlJ(1,3) = half*Uk(3)*dNkdx(1)*dNjdx(1)

            BnlJ(2,1) = half*Uk(1)*dNkdx(2)*dNjdx(2)
            BnlJ(2,2) = half*Uk(2)*dNkdx(2)*dNjdx(2)
            BnlJ(2,3) = half*Uk(3)*dNkdx(2)*dNjdx(2)

            BnlJ(3,1) = half*Uk(1)*dNkdx(3)*dNjdx(3)
            BnlJ(3,2) = half*Uk(2)*dNkdx(3)*dNjdx(3)
            BnlJ(3,3) = half*Uk(3)*dNkdx(3)*dNjdx(3)

            BnlJ(4,1) = voigt*half*Uk(1)*(dNkdx(1)*dNjdx(2)+dNkdx(2)*dNjdx(1))
            BnlJ(4,2) = voigt*half*Uk(2)*(dNkdx(1)*dNjdx(2)+dNkdx(2)*dNjdx(1))
            BnlJ(4,3) = voigt*half*Uk(3)*(dNkdx(1)*dNjdx(2)+dNkdx(2)*dNjdx(1))

            BnlJ(5,1) = voigt*half*Uk(1)*(dNkdx(1)*dNjdx(3)+dNkdx(3)*dNjdx(1))
            BnlJ(5,2) = voigt*half*Uk(2)*(dNkdx(1)*dNjdx(3)+dNkdx(3)*dNjdx(1))
            BnlJ(5,3) = voigt*half*Uk(3)*(dNkdx(1)*dNjdx(3)+dNkdx(3)*dNjdx(1))

            BnlJ(6,1) = voigt*half*Uk(1)*(dNkdx(2)*dNjdx(3)+dNkdx(3)*dNjdx(2))
            BnlJ(6,2) = voigt*half*Uk(2)*(dNkdx(2)*dNjdx(3)+dNkdx(3)*dNjdx(2))
            BnlJ(6,3) = voigt*half*Uk(3)*(dNkdx(2)*dNjdx(3)+dNkdx(3)*dNjdx(2))

            BoJ(:,:)  = BoJ(:,:) + BnlJ(:,:)
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

        SoJ(1) = SoJ(1) + SBoJ(1,1)*Uj(1) + SBoJ(1,2)*Uj(2) + SBoJ(1,3)*Uj(3)
        SoJ(2) = SoJ(2) + SBoJ(2,1)*Uj(1) + SBoJ(2,2)*Uj(2) + SBoJ(2,3)*Uj(3)
        SoJ(3) = SoJ(3) + SBoJ(3,1)*Uj(1) + SBoJ(3,2)*Uj(2) + SBoJ(3,3)*Uj(3)
        SoJ(4) = SoJ(4) + SBoJ(4,1)*Uj(1) + SBoJ(4,2)*Uj(2) + SBoJ(4,3)*Uj(3)
        SoJ(5) = SoJ(5) + SBoJ(5,1)*Uj(1) + SBoJ(5,2)*Uj(2) + SBoJ(5,3)*Uj(3)
        SoJ(6) = SoJ(6) + SBoJ(6,1)*Uj(1) + SBoJ(6,2)*Uj(2) + SBoJ(6,3)*Uj(3)
    enddo

     !! Loop on control points
    do nodi = 1, NNODE
        dNidx(:) = dRdx(:, nodi)

        !! Compute BlIt
        do i = 1,MCRD
            BlIt(i, i) = dNidx(i)
        enddo
        BlIt(1,4) = voigt*dNidx(2)
        BlIt(2,4) = voigt*dNidx(1)
        BlIt(1,5) = voigt*dNidx(3)
        BlIt(3,5) = voigt*dNidx(1)
        BlIt(2,6) = voigt*dNidx(3)
        BlIt(3,6) = voigt*dNidx(2)

        BoIt(:,:)  = BlIt(:,:)

        !! Loop on control points
        do nodj = 1, NNODE
            dNjdx(:) = dRdx(:, nodj)
            Uj(:) = Uelem(:,nodj)

            !! Compute BnlJt
            BnlIt(1,1) = Uj(1)*dNjdx(1)*dNidx(1)
            BnlIt(2,1) = Uj(2)*dNjdx(1)*dNidx(1)
            BnlIt(3,1) = Uj(3)*dNjdx(1)*dNidx(1)

            BnlIt(1,2) = Uj(1)*dNjdx(2)*dNidx(2)
            BnlIt(2,2) = Uj(2)*dNjdx(2)*dNidx(2)
            BnlIt(3,2) = Uj(3)*dNjdx(2)*dNidx(2)

            BnlIt(1,3) = Uj(1)*dNjdx(3)*dNidx(3)
            BnlIt(2,3) = Uj(2)*dNjdx(3)*dNidx(3)
            BnlIt(3,3) = Uj(3)*dNjdx(3)*dNidx(3)

            BnlIt(1,4) = voigt*Uj(1)*(dNjdx(1)*dNidx(2)+dNjdx(2)*dNidx(1))
            BnlIt(2,4) = voigt*Uj(2)*(dNjdx(1)*dNidx(2)+dNjdx(2)*dNidx(1))
            BnlIt(3,4) = voigt*Uj(3)*(dNjdx(1)*dNidx(2)+dNjdx(2)*dNidx(1))

            BnlIt(1,5) = voigt*Uj(1)*(dNjdx(1)*dNidx(3)+dNjdx(3)*dNidx(1))
            BnlIt(2,5) = voigt*Uj(2)*(dNjdx(1)*dNidx(3)+dNjdx(3)*dNidx(1))
            BnlIt(3,5) = voigt*Uj(3)*(dNjdx(1)*dNidx(3)+dNjdx(3)*dNidx(1))

            BnlIt(1,6) = voigt*Uj(1)*(dNjdx(2)*dNidx(3)+dNjdx(3)*dNidx(2))
            BnlIt(2,6) = voigt*Uj(2)*(dNjdx(2)*dNidx(3)+dNjdx(3)*dNidx(2))
            BnlIt(3,6) = voigt*Uj(3)*(dNjdx(2)*dNidx(3)+dNjdx(3)*dNidx(2))

            BoIt(:,:)  = BoIt(:,:) + BnlIt(:,:)
        end do 

        !! Assembly
        Fgeom(1, nodi) =  BoIt(1,1)*SoJ(1) + BoIt(1,2)*SoJ(2) + BoIt(1,3)*SoJ(3) + BoIt(1,4)*SoJ(4) &
        &               + BoIt(1,5)*SoJ(5) + BoIt(1,6)*SoJ(6) 
        Fgeom(2, nodi) =  BoIt(2,1)*SoJ(1) + BoIt(2,2)*SoJ(2) + BoIt(2,3)*SoJ(3) + BoIt(2,4)*SoJ(4) &
        &               + BoIt(2,5)*SoJ(5) + BoIt(2,6)*SoJ(6) 
        Fgeom(3, nodi) =  BoIt(3,1)*SoJ(1) + BoIt(3,2)*SoJ(2) + BoIt(3,3)*SoJ(3) + BoIt(3,4)*SoJ(4) &
        &               + BoIt(3,5)*SoJ(5) + BoIt(3,6)*SoJ(6) 

    enddo
    
end subroutine geomvector_byCP


