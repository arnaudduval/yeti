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

c     Evaluation des fonctions Nurbs et des derivees non nulles en un
c     point xi de l'espace parametrique

      subroutine evalLgrge(xi,R)

      !! TODO : this routine may be a doublon of an existing routine computing NURBS function -> to verify

      use parameters
      use nurbspatch

      Implicit none

c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)

c     Output variables :
c     ----------------
      Double precision :: R
      dimension R(nnode_patch)


c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, SumTot,
     &     SumTot_inv, SumXi
      dimension FN(Jpqr_patch(1)+1),dNdxi(Jpqr_patch(1)+1),
     &     FM(Jpqr_patch(2)+1),dMdEta(Jpqr_patch(2)+1),
     &     FL(Jpqr_patch(3)+1),dLdZeta(Jpqr_patch(3)+1)

      Integer ::  i,j,k, Ni, NumLoc, Na,Nb,Nc
      dimension Ni(dim_patch)

c     Fin declaration des variables
c
c     ..................................................................
c
c     Initialisation:

c     Pour cas 1D et 2D
      FM(1)      = one
      dMdEta(1)  = zero
      FL(1)      = one
      dLdZeta(1) = zero

c     Intervalles de l'element courant
      Do i = 1,dim_patch
         Ni(i) = Nijk_patch(i,current_elem)
      End do

c     Fin initialisation
c
c     ..................................................................
c
c     Calcul:

c     - Calculate univariate B-spline function
      CALL dersbasisfuns(Ni(1),Jpqr_patch(1),Nkv_patch(1),xi(1),
     &     Ukv1_patch(:),FN,dNdxi)
      If (dim_patch>1) then
         CALL dersbasisfuns(Ni(2),Jpqr_patch(2),Nkv_patch(2),xi(2),
     &        Ukv2_patch(:),FM,dMdeta)
      Endif
      If (dim_patch>2) then
         CALL dersbasisfuns(Ni(3),Jpqr_patch(3),Nkv_patch(3),xi(3),
     &        Ukv3_patch(:),FL,dLdZeta)
      Endif

c     - Build numerators and denominators
      NumLoc     = 0
      SumTot     = zero
      Do k = 0,Jpqr_patch(3)
         Do j = 0,Jpqr_patch(2)
            Do i = 0,Jpqr_patch(1)
               NumLoc = NumLoc+1

               R(NumLoc) = FN(Jpqr_patch(1)+1-i)*FM(Jpqr_patch(2)+1-j)
     &              *FL(Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)

               SumTot = SumTot + R(NumLoc)

            Enddo
         Enddo
      Enddo

c     - Divide by denominator to complete definition of fct and deriv.
      SumTot_inv = one/SumTot
      R(:) = R(:)*SumTot_inv

      End subroutine evalLgrge
