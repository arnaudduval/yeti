!! Calcul au point de gauss de la matrice elementaire pour le cas des
!! elements classiques solide 2D
!! -
!!
!! Stockage sans construire la matrice : renvoie les matrices 3x3 de
!! chaque couple de points de controle (partie triangulaire
!! supperieure uniquement)

subroutine stiffmatrix_HO_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, d2Rdx, d3Rdx, stiff) !! not use NDOFEL

    use parameters

    implicit None

    !! Input arguments :
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL   !! NDOFEL=MCRD*NNODE
    double precision, intent(in) :: ddsdde, dRdx, d2Rdx, d3Rdx
    dimension ddsdde(ntens, ntens)
    dimension dRdx(MCRD, NNODE)
    dimension d2Rdx(MCRD + 1, NNODE)
    dimension d3Rdx(2 * MCRD, NNODE)

    !! Output variables :
    !! ----------------
    double precision, intent(out) :: stiff
    dimension stiff(MCRD, MCRD, NNODE * (NNODE + 1) / 2)


    !! Local variables :
    !! ---------------
    double precision :: BoJ, BoIt, dNjdx, d2Njdx, d3Njdx, stiffLoc, SBoJ, dNidx, d2Nidx, d3Nidx
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), dNjdx(MCRD), d2Njdx(MCRD + 1), &
        &      d3Njdx(2 * MCRD), stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD), &
        &      dNidx(MCRD), d2Nidx(MCRD + 1), d3Nidx(2 * MCRD)

    integer :: nodi, nodj, idim, jdim, i, j, k, count

    !! Fin declaration des variables ....................................

    !! Initialisation ...................................................

    !! Initialisation des matrices stiff, BoI, BoJ

    stiff(:,:,:) = zero


    !! Boucle points de controle
    count = 1
    ! d3Njdx(:) = zero
    ! d3Nidx(:) = zero

    do nodj = 1, NNODE

        do jdim = 1, MCRD
            dNjdx(jdim)  = dRdx(jdim, nodj)
        end do
        do jdim = 1, MCRD + 1
            d2Njdx(jdim) = d2Rdx(jdim, nodj)
        end do
        do jdim = 1, 2 * MCRD
            d3Njdx(jdim) = d3Rdx(jdim, nodj)
        end do

        !! Calcul de la matrice BoJ
        BoJ(:,:) = zero

        BoJ(1, 1)  = dNjdx(1)
        BoJ(2, 2)  = dNjdx(2)
        BoJ(3, 1)  = (dNjdx(2)) / Ractwo
        BoJ(3, 2)  = (dNjdx(1)) / Ractwo

        BoJ(4, 1)  = d2Njdx(1)
        BoJ(5, 1)  = d2Njdx(2)
        BoJ(6, 1)  = (d2Njdx(3)) * Ractwo
        BoJ(7, 2)  = d2Njdx(1)
        BoJ(8, 2)  = d2Njdx(2)
        BoJ(9, 2)  = (d2Njdx(3)) * Ractwo

        BoJ(10, 1) = d3Njdx(1)
        BoJ(11, 1) = d3Njdx(2)
        BoJ(12, 1) = (d3Njdx(3)) * Racthree
        BoJ(13, 1) = (d3Njdx(4)) * Racthree
        BoJ(14, 2) = d3Njdx(1)
        BoJ(15, 2) = d3Njdx(2)
        BoJ(16, 2) = (d3Njdx(3)) * Racthree
        BoJ(17, 2) = (d3Njdx(4)) * Racthree

        !! Calcul du produit ddsdde*BoJ

        SBoJ(:,:) = zero
        do i = 1, ntens
            do j = 1, MCRD
                do k = 1, ntens
                    SBoJ(i, j) = SBoJ(i, j) + ddsdde(i, k) * BoJ(k, j)
                end do
            end do
        end do


        !! Calcul de la matrice BoIt
        do nodi = 1, nodj
            do idim = 1, MCRD
                dNidx(idim) = dRdx(idim, nodi)
            end do
            do idim = 1, MCRD + 1
                d2Nidx(idim) = d2Rdx(idim, nodi)
            end do
            do idim = 1, 2 * MCRD
                d3Nidx(idim)  = d3Rdx(idim, nodi)
            end do

            BoIt(:,:) = zero

            BoIt(1, 1)  = dNidx(1)
            BoIt(2, 2)  = dNidx(2)
            BoIt(1, 3)  = (dNidx(2)) / Ractwo
            BoIt(2, 3)  = (dNidx(1)) / Ractwo

            BoIt(1, 4)  = d2Nidx(1)
            BoIt(1, 5)  = d2Nidx(2)
            BoIt(1, 6)  = (d2Nidx(3)) * Ractwo
            BoIt(2, 7)  = d2Nidx(1)
            BoIt(2, 8)  = d2Nidx(2)
            BoIt(2, 9)  = (d2Nidx(3)) * Ractwo

            BoIt(1, 10) = d3Nidx(1)
            BoIt(1, 11) = d3Nidx(2)
            BoIt(1, 12) = (d3Nidx(3)) * Racthree
            BoIt(1, 13) = (d3Nidx(4)) * Racthree
            BoIt(2, 14) = d3Nidx(1)
            BoIt(2, 15) = d3Nidx(2)
            BoIt(2, 16) = (d3Nidx(3)) * Racthree
            BoIt(2, 17) = (d3Nidx(4)) * Racthree

            call MulMat(BoIt, SBoJ, stiffLoc, MCRD, MCRD, ntens)


            !! Assemblage
            stiff(:,:, count) = stiffLoc(:,:)
            count = count + 1

        end do
    end do
end Subroutine stiffmatrix_HO_byCP

!! Same as stiffmatrix_HO_byCP with only 1st gradient of strain taken into account
subroutine stiffmatrix_HO_byCP_1stG(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, d2Rdx, stiff) !! not use NDOFEL

    use parameters

    implicit None

    !! Input arguments :
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL   !! NDOFEL=MCRD*NNODE
    double precision, intent(in) :: ddsdde, dRdx, d2Rdx
    dimension ddsdde(ntens, ntens)
    dimension dRdx(MCRD, NNODE)
    dimension d2Rdx(MCRD + 1, NNODE)

    !! Output variables :
    !! ----------------
    double precision, intent(out) :: stiff
    dimension stiff(MCRD, MCRD, NNODE * (NNODE + 1) / 2)


    !! Local variables :
    !! ---------------
    double precision :: BoJ, BoIt, dNjdx, d2Njdx, stiffLoc, SBoJ, dNidx, d2Nidx
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), dNjdx(MCRD), d2Njdx(MCRD + 1), &
        &      stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD), &
        &      dNidx(MCRD), d2Nidx(MCRD + 1)

    integer :: nodi, nodj, idim, jdim, i, j, k, count

    !! Fin declaration des variables ....................................

    !! Initialisation ...................................................

    !! Initialisation des matrices stiff, BoI, BoJ

    stiff(:,:,:) = zero


    !! Boucle points de controle
    count = 1
    ! d3Njdx(:) = zero
    ! d3Nidx(:) = zero

    do nodj = 1, NNODE

        do jdim = 1, MCRD
            dNjdx(jdim)  = dRdx(jdim, nodj)
        end do
        do jdim = 1, MCRD + 1
            d2Njdx(jdim) = d2Rdx(jdim, nodj)
        end do

        !! Calcul de la matrice BoJ
        BoJ(:,:) = zero

        BoJ(1, 1)  = dNjdx(1)
        BoJ(2, 2)  = dNjdx(2)
        BoJ(3, 1)  = (dNjdx(2)) * inv_ractwo
        BoJ(3, 2)  = (dNjdx(1)) * inv_ractwo

        BoJ(4, 1)  = d2Njdx(1)
        BoJ(5, 1)  = d2Njdx(2)
        BoJ(6, 1)  = (d2Njdx(3)) * Ractwo
        BoJ(7, 2)  = d2Njdx(1)
        BoJ(8, 2)  = d2Njdx(2)
        BoJ(9, 2)  = (d2Njdx(3)) * Ractwo

        !! Calcul du produit ddsdde*BoJ

        SBoJ(:,:) = zero
        do i = 1, ntens
            do j = 1, MCRD
                do k = 1, ntens
                    SBoJ(i, j) = SBoJ(i, j) + ddsdde(i, k) * BoJ(k, j)
                end do
            end do
        end do


        !! Calcul de la matrice BoIt
        do nodi = 1, nodj
            do idim = 1, MCRD
                dNidx(idim) = dRdx(idim, nodi)
            end do
            do idim = 1, MCRD + 1
                d2Nidx(idim) = d2Rdx(idim, nodi)
            end do

            BoIt(:,:) = zero

            BoIt(1, 1)  = dNidx(1)
            BoIt(2, 2)  = dNidx(2)
            BoIt(1, 3)  = (dNidx(2)) * inv_ractwo
            BoIt(2, 3)  = (dNidx(1)) * inv_ractwo

            BoIt(1, 4)  = d2Nidx(1)
            BoIt(1, 5)  = d2Nidx(2)
            BoIt(1, 6)  = (d2Nidx(3)) * Ractwo
            BoIt(2, 7)  = d2Nidx(1)
            BoIt(2, 8)  = d2Nidx(2)
            BoIt(2, 9)  = (d2Nidx(3)) * Ractwo

            call MulMat(BoIt, SBoJ, stiffLoc, MCRD, MCRD, ntens)


            !! Assemblage
            stiff(:,:, count) = stiffLoc(:,:)
            count = count + 1

        end do
    end do
end Subroutine stiffmatrix_HO_byCP_1stG

