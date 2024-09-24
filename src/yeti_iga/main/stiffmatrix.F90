!! Copyright 2011 Florian Maurin
!! Copyright 2016-2019 Thibaut Hirschler

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


!! Compute elementary matrix at Gauss point, for 2D/3D solids
!! ---
!! Only upper triangular part is built. Lower part is kpt at zero
!! symmetry is taken into account during sum at Gauss points, at the
!! end of elementary built (in UELMAT.f90)

subroutine stiffmatrix(ntens,NNODE,MCRD,NDOFEL,ddsdde,dRdx,stiff)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL
    double precision, intent(in) :: ddsdde, dRdx
    dimension ddsdde(ntens, ntens), dRdx(MCRD, NNODE)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: stiff
    dimension stiff(NDOFEL, NDOFEL)

    !! Local variables
    !! ---------------
    double precision :: BoJ, BoIt, dNjdx, dNidx, stiffLoc, SBoJ
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), dNjdx(MCRD)
    dimension dNidx(MCRD), stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD)

    integer :: nodj, nodi, jdim, idim, i, j, jcol, irow

    stiff(:,:) = zero
    BoJ(:,:) = zero
    BoIt(:,:) = zero

    !! Loop on control points
    jcol = 0
    do nodj = 1, NNODE
        do jdim = 1, MCRD
            dNjdx(jdim) = dRdx(jdim, nodj)
        enddo

        !! Comput BoJ matrix
        do i = 1, MCRD
            BoJ(i, i) = dNjdx(i)
        enddo
        BoJ(4, 1) = dNjdx(2)
        BoJ(4, 2) = dNjdx(1)
        if (MCRD == 3) then
            BoJ(5, 1) = dNjdx(3)
            BoJ(5, 3) = dNjdx(1)
            BoJ(6, 2) = dNjdx(3)
            BoJ(6, 3) = dNjdx(2)
        endif

        !! Compute product ddsdde*BoJ
        SBoJ(:,:) = zero
        do j = 1,MCRD
            do i = 1,MCRD
                SBoJ(i, j) = ddsdde(i, j)*BoJ(j, j)
            enddo
            SBoJ(4, 1) = ddsdde(4, 4)*BoJ(4, 1)
            SBoJ(4, 2) = ddsdde(4, 4)*BoJ(4, 2)
            if (MCRD == 3) then
                SBoJ(5, 1) = ddsdde(5, 5)*BoJ(5, 1)
                SBoJ(5, 3) = ddsdde(5, 5)*BoJ(5, 3)
                SBoJ(6, 2) = ddsdde(4, 4)*BoJ(6, 2)
                SBoJ(6, 3) = ddsdde(4, 4)*BoJ(6, 3)
            endif
        enddo
        !SBoJ = MATMUL(ddsdde,BoJ)

        !! Second loop on control points
        irow = 0
        do nodi = 1, nodj !,NNODE
            do idim = 1, MCRD
                dNidx(idim) = dRdx(idim, nodi)
            enddo

            !! Compute BoIt
            do i = 1,MCRD
               BoIt(i,i) = dNidx(i)
            Enddo
            BoIt(1, 4) = dNidx(2)
            BoIt(2, 4) = dNidx(1)
            if (MCRD == 3) then
                BoIt(1, 5) = dNidx(3)
                BoIt(3, 5) = dNidx(1)
                BoIt(2, 6) = dNidx(3)
                BoIt(3, 6) = dNidx(2)
            endif


            !! Compute stiffLoc
            call MulMat(BoIt, SBoJ, stiffLoc, MCRD, MCRD, ntens)

            !! Assembly
            stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)            &
     &           = stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)     &
     &           + stiffLoc(:,:)

            irow = irow + MCRD
        enddo
        jcol = jcol + MCRD
    enddo

end subroutine stiffmatrix



!! Storage without building the matrix : return 3x3 matrices
!! for each couple of control points (upper triangle part only)
subroutine stiffmatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx,       &
     &     stiff)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL
    double precision, intent(in) :: ddsdde, dRdx
    dimension ddsdde(ntens, ntens), dRdx(MCRD, NNODE)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: stiff
    dimension stiff(MCRD, MCRD, NNODE*(NNODE+1)/2)

    !! Local variables
    !! ---------------
    double precision :: BoJ, BoIt, dNjdx, dNidx, stiffLoc, SBoJ
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), dNjdx(MCRD)
    dimension dNidx(MCRD), stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD)

    integer :: nodj, nodi, jdim, idim, i, j, count

    stiff(:,:,:) = zero
    BoJ(:,:)     = zero
    BoIt(:,:)    = zero

    !! Loop on control points
    count = 1
    do nodj = 1, NNODE
        do jdim = 1, MCRD
            dNjdx(jdim) = dRdx(jdim, nodj)
        enddo

        !! Compute BoJ matrix
        do i = 1,MCRD
            BoJ(i, i) = dNjdx(i)
        enddo
        BoJ(4,1) = dNjdx(2)
        BoJ(4,2) = dNjdx(1)
        if (MCRD == 3) then
            BoJ(5, 1) = dNjdx(3)
            BoJ(5, 3) = dNjdx(1)
            BoJ(6, 2) = dNjdx(3)
            BoJ(6, 3) = dNjdx(2)
        endif

        !! Compute product ddsdde*BoJ
        SBoJ(:, :) = zero
        do j = 1, MCRD
            do i = 1, MCRD
                SBoJ(i,j) = ddsdde(i, j)*BoJ(j, j)
            enddo
            SBoJ(4, 1) = ddsdde(4, 4)*BoJ(4, 1)
            SBoJ(4, 2) = ddsdde(4, 4)*BoJ(4, 2)
            if (MCRD == 3) then
                SBoJ(5, 1) = ddsdde(5, 5)*BoJ(5, 1)
                SBoJ(5, 3) = ddsdde(5, 5)*BoJ(5, 3)
                SBoJ(6, 2) = ddsdde(4, 4)*BoJ(6, 2)
                SBoJ(6, 3) = ddsdde(4, 4)*BoJ(6, 3)
            endif
        enddo

        !! Second loop on control points
        do nodi = 1, nodj
            do idim = 1, MCRD
                dNidx(idim) = dRdx(idim, nodi)
            enddo

            !! Compute BoIt
            do i = 1, MCRD
                BoIt(i, i) = dNidx(i)
            enddo
            BoIt(1, 4) = dNidx(2)
            BoIt(2, 4) = dNidx(1)
            if (MCRD == 3) then
                BoIt(1,5) = dNidx(3)
                BoIt(3,5) = dNidx(1)
                BoIt(2,6) = dNidx(3)
                BoIt(3,6) = dNidx(2)
            endif

            !! Compute stiffLoc
            call MulMat(BoIt, SBoJ, stiffLoc, MCRD, MCRD, ntens)

            !! Assembly
            stiff(:, :, count) = stiffLoc(:, :)

            count = count + 1
        enddo
    enddo

end subroutine stiffmatrix_byCP


!! Case os curvlinear cordinates formulation
subroutine stiffmatrix_curv(ntens, NNODE, MCRD, NDOFEL, ddsdde, AI,     &
     &     dRdxi, stiff)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: ntens, NNODE, MCRD, NDOFEL
    double precision, intent(in) :: ddsdde, dRdxi, AI
    dimension ddsdde(ntens, ntens), dRdxi(NNODE, MCRD), AI(3, MCRD)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: stiff
    dimension stiff(MCRD, MCRD, NNODE*(NNODE+1)/2)


    double precision :: BoJ, BoIt, dNjdxi, dNidxi, stiffLoc, SBoJ
    dimension BoJ(ntens, MCRD), BoIt(MCRD, ntens), dNjdxi(MCRD)
    dimension dNidxi(MCRD), stiffLoc(MCRD, MCRD), SBoJ(ntens, MCRD)

    integer :: nodj, nodi, jdim, idim, i, j, count

    stiff(:,:,:) = zero
    BoJ(:, :) = zero
    BoIt(:, :) = zero

    !! Loop on control points
    count = 1
    do nodj = 1, NNODE
        do jdim = 1, MCRD
            dNjdxi(jdim) = dRdxi(nodj, jdim)
        enddo

        !! Compute BoJ matrix
        do i = 1,MCRD
            BoJ(i,:) = dNjdxi(i)*AI(:MCRD, i)
        enddo
        BoJ(4,:) = dNjdxi(2)*AI(:MCRD, 1) + dNjdxi(1)*AI(:MCRD, 2)
        if (MCRD == 3) then
            BoJ(5, :) = dNjdxi(3)*AI(:MCRD, 1) + dNjdxi(1)*AI(:MCRD, 3)
            BoJ(6, :) = dNjdxi(3)*AI(:MCRD, 2) + dNjdxi(2)*AI(:MCRD, 3)
        endif

        !! Compute product ddsdde*BoJ
        SBoJ(:, :) = zero
        call MulMat(ddsdde, BoJ, SBoJ, ntens, MCRD, ntens)

        !! Second loop on control points
        do nodi = 1,nodj
            do idim = 1,MCRD
                dNidxi(idim) = dRdxi(nodi, idim)
            enddo

            !! COmpute BoIt
            do i = 1, MCRD
                BoIt(:, i) = dNidxi(i)*AI(:MCRD, i)
            enddo
            BoIt(:, 4) = dNidxi(2)*AI(:MCRD, 1) + dNidxi(1)*AI(:MCRD, 2)
            if (MCRD == 3) then
                BoIt(:, 5) = dNidxi(3)*AI(:MCRD, 1) + dNidxi(1)*AI(:MCRD, 3)
                BoIt(:, 6) = dNidxi(3)*AI(:MCRD, 2) + dNidxi(2)*AI(:MCRD, 3)
            endif

            !! Compute stiff loc
            call MulMat(BoIt, SBoJ, stiffLoc, MCRD, MCRD, ntens)

            !! Assembly
            stiff(:, :, count) = stiffLoc(:, :)

            count = count + 1
        enddo
    enddo

end subroutine stiffmatrix_curv

