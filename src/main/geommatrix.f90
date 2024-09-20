!! Compute elementary matrix at Gauss point, for 2D/3D solids
!! ---
!! Only upper triangular part is built. Lower part is kpt at zero
!! symmetry is taken into account during sum at Gauss points, at the
!! end of elementary built (in UELTANMAT.f90)


!! Storage without building the matrix : return 3x3 matrices
!! for each couple of control points (upper triangle part only)
subroutine geommatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, Uelem,       &
     &     geom)

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
    double precision, intent(out) :: geom
    dimension geom(MCRD, MCRD, NNODE*(NNODE+1)/2)

    !! Local variables
    !! ---------------
    double precision :: BoJ, BoI, dNjdx, dNidx, dNkdx, geomLoc, SBoI, BlI, BnlI, BlJt, BnlJt, Uk
    dimension BoI(ntens, MCRD), BoJ(MCRD, ntens), BlJt(MCRD, ntens), BnlJt(MCRD, ntens), BlI(ntens, MCRD), BnlI(ntens, MCRD)
    dimension dNidx(MCRD), dNjdx(MCRD), dNkdx(MCRD), geomLoc(MCRD, MCRD), SBoI(ntens, MCRD), Uk(MCRD)
    integer :: nodj, nodi, nodk, jdim, idim, i, j, k, count
    double precision :: res, voigt

    geom(:,:,:) = zero
    BlI(:,:)      = zero
    BlJt(:,:)      = zero
    voigt = sqrt(2.0D0)/2


    count = 1
    !! Loop on control points
    do nodi = 1, NNODE
        dNidx(:) = dRdx(:, nodi)

        ! Compute BlJ
        do i = 1,MCRD
            BlI(i, i) = dNidx(i)
        enddo
        BlI(4,1) = voigt*dNidx(2)
        BlI(4,2) = voigt*dNidx(1)
        BlI(5,1) = voigt*dNidx(3)
        BlI(5,3) = voigt*dNidx(1)
        BlI(6,2) = voigt*dNidx(3)
        BlI(6,3) = voigt*dNidx(2)

        BoI(:,:)  = BlI(:,:)

        !! Second loop on control points
        do nodk = 1, NNODE
            dNkdx(:) = dRdx(:, nodk)
            Uk(:) = Uelem(:,nodk)

        !! Compute BnlJ
            BnlI(1,1) = Uk(1)*dNkdx(1)*dNidx(1)
            BnlI(1,2) = Uk(2)*dNkdx(1)*dNidx(1)
            BnlI(1,3) = Uk(3)*dNkdx(1)*dNidx(1)

            BnlI(2,1) = Uk(1)*dNkdx(2)*dNidx(2)
            BnlI(2,2) = Uk(2)*dNkdx(2)*dNidx(2)
            BnlI(2,3) = Uk(3)*dNkdx(2)*dNidx(2)

            BnlI(3,1) = Uk(1)*dNkdx(3)*dNidx(3)
            BnlI(3,2) = Uk(2)*dNkdx(3)*dNidx(3)
            BnlI(3,3) = Uk(3)*dNkdx(3)*dNidx(3)

            BnlI(4,1) = voigt*Uk(1)*(dNkdx(1)*dNidx(2)+dNkdx(2)*dNidx(1))
            BnlI(4,2) = voigt*Uk(2)*(dNkdx(1)*dNidx(2)+dNkdx(2)*dNidx(1))
            BnlI(4,3) = voigt*Uk(3)*(dNkdx(1)*dNidx(2)+dNkdx(2)*dNidx(1))

            BnlI(5,1) = voigt*Uk(1)*(dNkdx(1)*dNidx(3)+dNkdx(3)*dNidx(1))
            BnlI(5,2) = voigt*Uk(2)*(dNkdx(1)*dNidx(3)+dNkdx(3)*dNidx(1))
            BnlI(5,3) = voigt*Uk(3)*(dNkdx(1)*dNidx(3)+dNkdx(3)*dNidx(1))

            BnlI(6,1) = voigt*Uk(1)*(dNkdx(2)*dNidx(3)+dNkdx(3)*dNidx(2))
            BnlI(6,2) = voigt*Uk(2)*(dNkdx(2)*dNidx(3)+dNkdx(3)*dNidx(2))
            BnlI(6,3) = voigt*Uk(3)*(dNkdx(2)*dNidx(3)+dNkdx(3)*dNidx(2))

            BoI(:,:)  = BoI(:,:) + BnlI(:,:)
        enddo

         !! Compute product ddsdde*BoJ
        SBoI(:, :) = zero
        do j=1,MCRD
            do i=1,MCRD
                res = 0
                do k=1,MCRD
                    res = res + ddsdde(i,k)*BoI(k,j)
                end do
                sBoI(i,j) = res
            enddo
        enddo
        SBoI(4, 1) = ddsdde(4, 4)*BoI(4, 1)
        SBoI(4, 2) = ddsdde(4, 4)*BoI(4, 2)
        SBoI(4, 3) = ddsdde(4, 4)*BoI(4, 3)
        SBoI(5, 1) = ddsdde(5, 5)*BoI(5, 1)
        SBoI(5, 2) = ddsdde(5, 5)*BoI(5, 2)
        SBoI(5, 3) = ddsdde(5, 5)*BoI(5, 3)
        SBoI(6, 1) = ddsdde(6, 6)*BoI(6, 1)
        SBoI(6, 2) = ddsdde(6, 6)*BoI(6, 2)
        SBoI(6, 3) = ddsdde(6, 6)*BoI(6, 3)


        !! Second loop on control points
        do nodj = 1, nodi
            dNjdx(:) = dRdx(:, nodj)


            !! Compute BlJt
            do i = 1,MCRD
                BlJt(i, i) = dNjdx(i)
            enddo
            BlJt(1,4) = voigt*dNjdx(2)
            BlJt(2,4) = voigt*dNjdx(1)
            BlJt(1,5) = voigt*dNjdx(3)
            BlJt(3,5) = voigt*dNjdx(1)
            BlJt(2,6) = voigt*dNjdx(3)
            BlJt(3,6) = voigt*dNjdx(2)

            BoJ(:,:)  = BlJt(:,:)

            !! Third loop on control points
            do nodk = 1, NNODE
                dNkdx(:) = dRdx(:, nodk)
                Uk(:) = Uelem(:,nodk)

                !! Compute BnlJ
                BnlJt(1,1) = Uk(1)*dNjdx(1)*dNkdx(1)
                BnlJt(2,1) = Uk(2)*dNjdx(1)*dNkdx(1)
                BnlJt(3,1) = Uk(3)*dNjdx(1)*dNkdx(1)

                BnlJt(1,2) = Uk(1)*dNjdx(2)*dNkdx(2)
                BnlJt(2,2) = Uk(2)*dNjdx(2)*dNkdx(2)
                BnlJt(3,2) = Uk(3)*dNjdx(2)*dNkdx(2)

                BnlJt(1,3) = Uk(1)*dNjdx(3)*dNkdx(3)
                BnlJt(2,3) = Uk(2)*dNjdx(3)*dNkdx(3)
                BnlJt(3,3) = Uk(3)*dNjdx(3)*dNkdx(3)

                BnlJt(1,4) = voigt*Uk(1)*(dNjdx(1)*dNkdx(2)+dNjdx(2)*dNkdx(1))
                BnlJt(2,4) = voigt*Uk(2)*(dNjdx(1)*dNkdx(2)+dNjdx(2)*dNkdx(1))
                BnlJt(3,4) = voigt*Uk(3)*(dNjdx(1)*dNkdx(2)+dNjdx(2)*dNkdx(1))

                BnlJt(1,5) = voigt*Uk(1)*(dNjdx(1)*dNkdx(3)+dNjdx(3)*dNkdx(1))
                BnlJt(2,5) = voigt*Uk(2)*(dNjdx(1)*dNkdx(3)+dNjdx(3)*dNkdx(1))
                BnlJt(3,5) = voigt*Uk(3)*(dNjdx(1)*dNkdx(3)+dNjdx(3)*dNkdx(1))

                BnlJt(1,6) = voigt*Uk(1)*(dNjdx(2)*dNkdx(3)+dNjdx(3)*dNkdx(2))
                BnlJt(2,6) = voigt*Uk(2)*(dNjdx(2)*dNkdx(3)+dNjdx(3)*dNkdx(2))
                BnlJt(3,6) = voigt*Uk(3)*(dNjdx(2)*dNkdx(3)+dNjdx(3)*dNkdx(2))

                BoJ(:,:)  = BoJ(:,:) + BnlJt(:,:)
            enddo

        call MulMat(BoJ, SBoI, geomLoc, MCRD, MCRD, ntens)

        !! Assembly
        geom(:, :, count) = geomLoc(:, :)
        count = count + 1
        enddo
    enddo
end subroutine geommatrix_byCP



