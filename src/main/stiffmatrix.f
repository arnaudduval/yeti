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

c     --
c     Calcul au point de gauss de la matrice elementaire pour le cas des
c     elements classiques solide 2D/3D
c     -
c     On construit uniquement la partie triangulaire superieure. La
c     partie inferieure est laissee nulle. La symmetrie est prise en
c     compte lors de la sommation sur les points de gauss, a la fin de
c     la construiction elementaire (dans UELMAT.f).
c     --

      Subroutine stiffmatrix(ntens,NNODE,MCRD,NDOFEL,ddsdde,dRdx,stiff)

      Implicit None

c     Input arguments :
c     ---------------
      Integer, intent(in) :: ntens,NNODE,MCRD,NDOFEL
      Double precision, intent(in) :: ddsdde,dRdx
      dimension ddsdde(ntens,ntens),dRdx(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(NDOFEL,NDOFEL)

c     Local variables :
c     ---------------
      Double precision :: BoJ,BoIt,dNjdx,dNidx,stiffLoc,SBoJ
      dimension BoJ(ntens,MCRD),BoIt(MCRD,ntens),dNjdx(MCRD),
     &     dNidx(MCRD),stiffLoc(MCRD,MCRD),SBoJ(ntens,MCRD)

      Integer :: nodj,nodi, jdim,idim, i,j, jcol,irow

C     Fin declaration des variables ....................................
c
c
c
c
c     Initialisation ...................................................

c     Initialisation des matrices stiff, BoI, BoJ
      stiff(:,:) = 0.0
      BoJ(:,:) = 0.0
      !BoI(:,:) = zero
      BoIt(:,:) = 0.0

c     Boucle points de controle
      jcol = 0
      Do nodj = 1,NNODE
         Do jdim = 1,MCRD
            dNjdx(jdim) = dRdx(jdim,nodj)
         Enddo

c     Calcul de la matrice BoJ
         Do i = 1,MCRD
            BoJ(i,i) = dNjdx(i)
         Enddo
         BoJ(4,1) = dNjdx(2)
         BoJ(4,2) = dNjdx(1)
         If (MCRD==3) then
            BoJ(5,1) = dNjdx(3)
            BoJ(5,3) = dNjdx(1)
            BoJ(6,2) = dNjdx(3)
            BoJ(6,3) = dNjdx(2)
         Endif

c     Calcul du produit ddsdde*BoJ
         SBoJ(:,:) = 0.0
         Do j = 1,MCRD
            Do i = 1,MCRD
               SBoJ(i,j) = ddsdde(i,j)*BoJ(j,j)
            Enddo
            SBoJ(4,1) = ddsdde(4,4)*BoJ(4,1)
            SBoJ(4,2) = ddsdde(4,4)*BoJ(4,2)
            If (MCRD==3) then
               SBoJ(5,1) = ddsdde(5,5)*BoJ(5,1)
               SBoJ(5,3) = ddsdde(5,5)*BoJ(5,3)
               SBoJ(6,2) = ddsdde(4,4)*BoJ(6,2)
               SBoJ(6,3) = ddsdde(4,4)*BoJ(6,3)
            Endif
         Enddo
         !SBoJ = MATMUL(ddsdde,BoJ)

c     Deuxieme boucle points de controle
         irow = 0
         Do nodi = 1,nodj !,NNODE
            Do idim = 1,MCRD
               dNidx(idim) = dRdx(idim,nodi)
            Enddo

c     Calcul BoI
            Do i = 1,MCRD
               !BoI(i,i) = dNidx(i)
               BoIt(i,i) = dNidx(i)
            Enddo
            !BoI(4,1) = dNidx(2)
            !BoI(4,2) = dNidx(1)
            BoIt(1,4) = dNidx(2)
            BoIt(2,4) = dNidx(1)
            If (MCRD==3) then
               !BoI(5,1) = dNidx(3)
               !BoI(5,3) = dNidx(1)
               BoIt(1,5) = dNidx(3)
               BoIt(3,5) = dNidx(1)
               !BoI(6,2) = dNidx(3)
               !BoI(6,3) = dNidx(2)
               BoIt(2,6) = dNidx(3)
               BoIt(3,6) = dNidx(2)
            Endif


c     Calcul stiffLoc
            call MulMat(BoIt,SBoJ, stiffLoc, MCRD,MCRD,ntens)

c     Assemblage
            stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)
     &           = stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)
     &           + stiffLoc(:,:)

            irow = irow + MCRD
         Enddo
         jcol = jcol + MCRD
      Enddo

c     Symmetry : lower part
c      Do j = 1,NDOFEL-1
c         Do i = j+1,NDOFEL
c            stiff(i,j) = stiff(j,i)
c         Enddo
c      Enddo

      End SUBROUTINE stiffmatrix











c     --
c     Stockage sans construire la matrice : renvoie les matrices 3x3 de
c     chaque couple de points de controle (partie triangulaire
c     supperieure uniquement)
c     --

      Subroutine stiffmatrix_byCP(ntens,NNODE,MCRD,NDOFEL,ddsdde,dRdx,
     &     stiff)

      use parameters

      Implicit None

c     Input arguments :
c     ---------------
      Integer, intent(in) :: ntens,NNODE,MCRD,NDOFEL
      Double precision, intent(in) :: ddsdde,dRdx
      dimension ddsdde(ntens,ntens),dRdx(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(MCRD,MCRD,NNODE*(NNODE+1)/2)

c     Local variables :
c     ---------------
      Double precision :: BoJ,BoIt,dNjdx,dNidx,stiffLoc,SBoJ
      dimension BoJ(ntens,MCRD),BoIt(MCRD,ntens),dNjdx(MCRD),
     &     dNidx(MCRD),stiffLoc(MCRD,MCRD),SBoJ(ntens,MCRD)

      Integer :: nodj,nodi, jdim,idim, i,j, count

C     Fin declaration des variables ....................................
c
c
c
c
c     Initialisation ...................................................

c     Initialisation des matrices stiff, BoI, BoJ
      stiff(:,:,:) = zero
      BoJ(:,:)     = zero
      BoIt(:,:)    = zero

c     Boucle points de controle
      count = 1
      Do nodj = 1,NNODE
         Do jdim = 1,MCRD
            dNjdx(jdim) = dRdx(jdim,nodj)
         Enddo

c     Calcul de la matrice BoJ
         Do i = 1,MCRD
            BoJ(i,i) = dNjdx(i)
         Enddo
         BoJ(4,1) = dNjdx(2)
         BoJ(4,2) = dNjdx(1)
         If (MCRD==3) then
            BoJ(5,1) = dNjdx(3)
            BoJ(5,3) = dNjdx(1)
            BoJ(6,2) = dNjdx(3)
            BoJ(6,3) = dNjdx(2)
         Endif

c     Calcul du produit ddsdde*BoJ
         SBoJ(:,:) = zero
         Do j = 1,MCRD
            Do i = 1,MCRD
               SBoJ(i,j) = ddsdde(i,j)*BoJ(j,j)
            Enddo
            SBoJ(4,1) = ddsdde(4,4)*BoJ(4,1)
            SBoJ(4,2) = ddsdde(4,4)*BoJ(4,2)
            If (MCRD==3) then
               SBoJ(5,1) = ddsdde(5,5)*BoJ(5,1)
               SBoJ(5,3) = ddsdde(5,5)*BoJ(5,3)
               SBoJ(6,2) = ddsdde(4,4)*BoJ(6,2)
               SBoJ(6,3) = ddsdde(4,4)*BoJ(6,3)
            Endif
         Enddo

c     Deuxieme boucle points de controle
         Do nodi = 1,nodj
            Do idim = 1,MCRD
               dNidx(idim) = dRdx(idim,nodi)
            Enddo

c     Calcul BoI
            Do i = 1,MCRD
               BoIt(i,i) = dNidx(i)
            Enddo
            BoIt(1,4) = dNidx(2)
            BoIt(2,4) = dNidx(1)
            If (MCRD==3) then
               BoIt(1,5) = dNidx(3)
               BoIt(3,5) = dNidx(1)
               BoIt(2,6) = dNidx(3)
               BoIt(3,6) = dNidx(2)
            Endif


c     Calcul stiffLoc
            call MulMat(BoIt,SBoJ, stiffLoc, MCRD,MCRD,ntens)

c     Assemblage
            stiff(:,:,count) = stiffLoc(:,:)

            count = count + 1
         Enddo
      Enddo

      End SUBROUTINE stiffmatrix_byCP














c     --
c     Cas de la formulation en coordonnees curvilignes
c     --

      Subroutine stiffmatrix_curv(ntens,NNODE,MCRD,NDOFEL,ddsdde,AI,
     &     dRdxi,stiff)

      use parameters

      Implicit None

c     Input arguments :
c     ---------------
      Integer, intent(in) :: ntens,NNODE,MCRD,NDOFEL
      Double precision, intent(in) :: ddsdde,dRdxi,AI
      dimension ddsdde(ntens,ntens),dRdxi(NNODE,MCRD),AI(3,MCRD)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(MCRD,MCRD,NNODE*(NNODE+1)/2)

c     Local variables :
c     ---------------
      Double precision :: BoJ,BoIt,dNjdxi,dNidxi,stiffLoc,SBoJ
      dimension BoJ(ntens,MCRD),BoIt(MCRD,ntens),dNjdxi(MCRD),
     &     dNidxi(MCRD),stiffLoc(MCRD,MCRD),SBoJ(ntens,MCRD)

      Integer :: nodj,nodi, jdim,idim, i,j, count

C     Fin declaration des variables ....................................
c
c
c
c
c     Initialisation ...................................................

c     Initialisation des matrices stiff, BoI, BoJ
      stiff(:,:,:) = zero
      BoJ(:,:)     = zero
      BoIt(:,:)    = zero

c     Boucle points de controle
      count = 1
      Do nodj = 1,NNODE
         Do jdim = 1,MCRD
            dNjdxi(jdim) = dRdxi(nodj,jdim)
         Enddo

c     Calcul de la matrice BoJ
         Do i = 1,MCRD
            BoJ(i,:) = dNjdxi(i)*AI(:MCRD,i)
         Enddo
         BoJ(4,:) = dNjdxi(2)*AI(:MCRD,1) + dNjdxi(1)*AI(:MCRD,2)
         If (MCRD==3) then
            BoJ(5,:) = dNjdxi(3)*AI(:MCRD,1) + dNjdxi(1)*AI(:MCRD,3)
            BoJ(6,:) = dNjdxi(3)*AI(:MCRD,2) + dNjdxi(2)*AI(:MCRD,3)
         Endif

c     Calcul du produit ddsdde*BoJ
         SBoJ(:,:) = zero
         call MulMat(ddsdde,BoJ, SBoJ, ntens,MCRD,ntens)

c     Deuxieme boucle points de controle
         Do nodi = 1,nodj
            Do idim = 1,MCRD
               dNidxi(idim) = dRdxi(nodi,idim)
            Enddo

c     Calcul BoI
            Do i = 1,MCRD
               BoIt(:,i) = dNidxi(i)*AI(:MCRD,i)
            Enddo
            BoIt(:,4) = dNidxi(2)*AI(:MCRD,1) + dNidxi(1)*AI(:MCRD,2)
            If (MCRD==3) then
               BoIt(:,5) = dNidxi(3)*AI(:MCRD,1) + dNidxi(1)*AI(:MCRD,3)
               BoIt(:,6) = dNidxi(3)*AI(:MCRD,2) + dNidxi(2)*AI(:MCRD,3)
            Endif

c     Calcul stiffLoc
            call MulMat(BoIt,SBoJ, stiffLoc, MCRD,MCRD,ntens)

c     Assemblage
            stiff(:,:,count) = stiffLoc(:,:)

            count = count + 1
         Enddo
      Enddo


      End SUBROUTINE stiffmatrix_curv

