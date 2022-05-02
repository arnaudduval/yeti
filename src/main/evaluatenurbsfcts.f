!! Copyright 2016-2018 Thibaut Hirschler
!! Copyright 2021 Arnaud Duval

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


C     Evaluation of non zero NURBS functions, at a given point of parametric space
C     Derivative is not computed
      subroutine evalnurbs_noder(xi, R)
        use parameters
        use nurbspatch

        implicit none

C       input arguments
        double precision :: xi
        dimension xi(dim_patch)

C       output variables
        double precision :: R
        dimension R(nnode_patch)

C       Local variables
        double precision :: FN, FM, FL
        double precision :: SumTot, SumTot_inv, SumXi
        dimension FN(1, Jpqr_patch(1)+1), FM(1, Jpqr_patch(2)+1),
     &            FL(1, Jpqr_patch(3)+1)
        dimension SumXi(3)
        integer :: Ni
        dimension Ni(dim_patch)
        integer :: i, j, k, NumLoc

C       Initialization
        FM(1,1) = one
        FL(1,1) = one

C       Get knot span
        do i = 1, dim_patch
            Ni(i) = Nijk_patch(i, current_elem)
        enddo

        call gen_dersbasisfuns(Ni(1), Jpqr_patch(1), Nkv_patch(1), 
     &      xi(1), Ukv1_patch(:), 0, FN)
        if (dim_patch .gt. 1) then
            call gen_dersbasisfuns(Ni(2), Jpqr_patch(2), Nkv_patch(2), 
     &      xi(2), Ukv2_patch(:), 0, FM)
        endif 
        if (dim_patch .gt. 2) then
            call gen_dersbasisfuns(Ni(3), Jpqr_patch(3), Nkv_patch(3), 
     &      xi(3), Ukv3_patch(:), 0, FL)
        endif

C       Build numerators and denominators
        NumLoc = 0
        SumTot = zero
        SumXi(:) = zero
        do k = 0, Jpqr_patch(3)
            do j = 0, Jpqr_patch(2)
                do i = 0, Jpqr_patch(1)
                    NumLoc = Numloc+1

                    R(NumLoc) = FN(1,Jpqr_patch(1)+1-i)
     &                  *FM(1,Jpqr_patch(2)+1-j)
     &                  *FL(1,Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)

                    SumTot = SumTot + R(NumLoc)
                enddo
            enddo
        enddo

C       Divide by dnominator to complete definition of function
        SumTot_inv = one/SumTot
        do NumLoc = 1, nnode_patch
            R(NumLoc) = R(NumLoc)*SumTot_inv
        enddo

      end subroutine evalnurbs_noder


c     Evaluation des fonctions Nurbs et des derivees non nulles en un 
c     point xi de l'espace parametrique
      
      subroutine evalnurbs(xi,R,dRdxi)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)
      
c     Output variables :
c     ----------------
      Double precision :: R,dRdxi
      dimension R(nnode_patch),dRdxi(nnode_patch,3)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, SumTot,
     &     SumTot_inv, SumXi
      dimension FN(Jpqr_patch(1)+1),dNdxi(Jpqr_patch(1)+1),
     &     FM(Jpqr_patch(2)+1),dMdEta(Jpqr_patch(2)+1),
     &     FL(Jpqr_patch(3)+1),dLdZeta(Jpqr_patch(3)+1),SumXi(3)
      
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
      
c     - Compute univariate B-spline function
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
      SumXi(:)   = zero
      dRdxi(:,:) = zero
      Do k = 0,Jpqr_patch(3)
         Do j = 0,Jpqr_patch(2)
            Do i = 0,Jpqr_patch(1)
               NumLoc = NumLoc+1
               
               R(NumLoc) = FN(Jpqr_patch(1)+1-i)*FM(Jpqr_patch(2)+1-j)
     &              *FL(Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)
               
               SumTot = SumTot + R(NumLoc)
c     
               dRdxi(NumLoc,1) = 
     &              dNdxi(Jpqr_patch(1)+1-i)*FM(Jpqr_patch(2)+1-j)
     &              *FL(Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)
               SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)
c     
               dRdxi(NumLoc,2)=
     &              FN(Jpqr_patch(1)+1-i)*dMdEta(Jpqr_patch(2)+1-j)
     &              *FL(Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)
               SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)
c     
               dRdxi(NumLoc,3) = 
     &              FN(Jpqr_patch(1)+1-i)*FM(Jpqr_patch(2)+1-j)
     &              *dLdZeta(Jpqr_patch(3)+1-k)*Weight_elem(NumLoc)
               SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
            Enddo
         Enddo
      Enddo
      
c     - Divide by denominator to complete definition of fct and deriv.
      SumTot_inv = one/SumTot
      Do NumLoc = 1,nnode_patch
         R(NumLoc) = R(NumLoc)*SumTot_inv
         Do i = 1,dim_patch
            dRdxi(NumLoc,i)
     &           = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
         Enddo
      Enddo
            
      End subroutine evalnurbs




      
c     Evaluation des fonctions Nurbs et des derivees non nulles en un 
c     point xi de l'espace parametrique
      
      subroutine evalnurbs_mapping(xi,R,dRdxi)
      
      use parameters
      use embeddedMapping
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)
      
c     Output variables :
c     ----------------
      Double precision :: R,dRdxi
      dimension R(nnode_map),dRdxi(nnode_map,3)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, SumTot,
     &     SumTot_inv, SumXi
      dimension FN(Jpqr_map(1)+1),dNdxi(Jpqr_map(1)+1),
     &     FM(Jpqr_map(2)+1),dMdEta(Jpqr_map(2)+1),
     &     FL(Jpqr_map(3)+1),dLdZeta(Jpqr_map(3)+1),SumXi(3)
      
      Integer ::  i,j,k, Ni, NumLoc, Na,Nb,Nc
      dimension Ni(dim_map)
      
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
      Do i = 1,dim_map
         Ni(i) = Nijk_map(i,current_map_elem)
      End do
      
c     Fin initialisation
c     
c     ..................................................................
c     
c     Calcul:
      
c     - Calculate univariate B-spline function
      CALL dersbasisfuns(Ni(1),Jpqr_map(1),Nkv_map(1),xi(1),
     &     Ukv1_map(:),FN,dNdxi)
      If (dim_map>1) then
         CALL dersbasisfuns(Ni(2),Jpqr_map(2),Nkv_map(2),xi(2),
     &        Ukv2_map(:),FM,dMdeta)
      Endif
      If (dim_map>2) then
         CALL dersbasisfuns(Ni(3),Jpqr_map(3),Nkv_map(3),xi(3),
     &        Ukv3_map(:),FL,dLdZeta)
      Endif
      
c     - Build numerators and denominators
      NumLoc     = 0
      SumTot     = zero
      SumXi(:)   = zero
      dRdxi(:,:) = zero
      Do k = 0,Jpqr_map(3)
         Do j = 0,Jpqr_map(2)
            Do i = 0,Jpqr_map(1)
               NumLoc = NumLoc+1
               
               R(NumLoc) = FN(Jpqr_map(1)+1-i)*FM(Jpqr_map(2)+1-j)
     &              *FL(Jpqr_map(3)+1-k)*Weight_map_elem(NumLoc)
               
               SumTot = SumTot + R(NumLoc)
c     
               dRdxi(NumLoc,1) = 
     &              dNdxi(Jpqr_map(1)+1-i)*FM(Jpqr_map(2)+1-j)
     &              *FL(Jpqr_map(3)+1-k)*Weight_map_elem(NumLoc)
               SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)
c     
               dRdxi(NumLoc,2)=
     &              FN(Jpqr_map(1)+1-i)*dMdEta(Jpqr_map(2)+1-j)
     &              *FL(Jpqr_map(3)+1-k)*Weight_map_elem(NumLoc)
               SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)
c     
               dRdxi(NumLoc,3) = 
     &              FN(Jpqr_map(1)+1-i)*FM(Jpqr_map(2)+1-j)
     &              *dLdZeta(Jpqr_map(3)+1-k)*Weight_map_elem(NumLoc)
               SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
            Enddo
         Enddo
      Enddo
      
c     - Divide by denominator to complete definition of fct and deriv.
      SumTot_inv = one/SumTot
      Do NumLoc = 1,nnode_map
         R(NumLoc) = R(NumLoc)*SumTot_inv
         Do i = 1,dim_map
            dRdxi(NumLoc,i)
     &           = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
         Enddo
      Enddo
            
      End subroutine evalnurbs_mapping

      
      


      
c     Evaluation des fonctions Nurbs et des derivees non nulles en un 
c     point xi de l'espace parametrique pour element immerge
      
      subroutine evalnurbs_mapping_noder(xi, R)
      
      use parameters
      use embeddedMapping
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)
      
c     Output variables :
c     ----------------
      Double precision :: R
      dimension R(nnode_map)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN, FM, FL, SumTot, SumTot_inv
      dimension FN(Jpqr_map(1)+1), FM(Jpqr_map(2)+1), FL(Jpqr_map(3)+1)
      
      Integer ::  i, j, k, Ni, NumLoc, Na, Nb, Nc
      dimension Ni(dim_map)
      
c     Fin declaration des variables
c
c     ..................................................................
c
c     Initialisation:
      
c     Pour cas 1D et 2D
      FM(1) = one
      FL(1) = one
      
c     Intervalles de l'element courant
      do i = 1, dim_map
         Ni(i) = Nijk_map(i, current_map_elem)
      enddo
      
c     Fin initialisation
c     
c     ..................................................................
c     
c     Calcul:
      
c     - Calculate univariate B-spline function
      call gen_dersbasisfuns(Ni(1), Jpqr_map(1), Nkv_map(1), xi(1),
     &     Ukv1_map(:), 0, FN)
      if (dim_map .gt. 1) then
          call gen_dersbasisfuns(Ni(2), Jpqr_map(2), Nkv_map(2), xi(2),
     &    Ukv2_map(:), 0, FM)
      endif
      if (dim_map .gt. 2) then
          call gen_dersbasisfuns(Ni(3), Jpqr_map(3), Nkv_map(3), xi(3),
     &    Ukv3_map(:), 0, FL)
      endif
      
c     - Build numerators and denominators
      NumLoc = 0
      SumTot = zero

      Do k = 0, Jpqr_map(3)
         Do j = 0, Jpqr_map(2)
            Do i = 0, Jpqr_map(1)
               NumLoc = NumLoc + 1
               
               R(NumLoc) = FN(Jpqr_map(1)+1-i)
     &             *FM(Jpqr_map(2)+1-j)
     &             *FL(Jpqr_map(3)+1-k)*Weight_map_elem(NumLoc)
               
               SumTot = SumTot + R(NumLoc)

            Enddo
         Enddo
      Enddo
      
c     - Divide by denominator to complete definition of fct
      SumTot_inv = one/SumTot
      do NumLoc = 1, nnode_map
         R(NumLoc) = R(NumLoc)*SumTot_inv
      enddo
            
      end subroutine evalnurbs_mapping_noder
      
      
      


      
c     Evaluation des fonctions Nurbs et des derivees premieres et 
c     secondes en un point xi pour l'element courrant (definit dans le 
c     module nurbspatch)
      
      subroutine evalnurbs_w2ndDerv(xi,R,dRdxi,ddRddxi)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)
      
c     Output variables :
c     ----------------
      Double precision :: R,dRdxi,ddRddxi
      dimension R(nnode_patch),dRdxi(nnode_patch,3),
     &     ddRddxi(nnode_patch,6)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, SumTot,
     &     SumTot_inv, SumXi, ddNddxi,ddMddEta,ddLddZeta,SumXixi
      dimension FN(Jpqr_patch(1)+1),dNdxi(Jpqr_patch(1)+1),
     &     FM(Jpqr_patch(2)+1),dMdEta(Jpqr_patch(2)+1),
     &     FL(Jpqr_patch(3)+1),dLdZeta(Jpqr_patch(3)+1),SumXi(3),
     &     ddNddxi(Jpqr_patch(1)+1),ddMddEta(Jpqr_patch(2)+1),
     &     ddLddZeta(Jpqr_patch(3)+1),SumXixi(6)
      
      Integer ::  i,j,k, Ni, NumLoc, Na,Nb,Nc
      dimension Ni(dim_patch)
      
c     Fin declaration des variables
c
c     ..................................................................
c
c     Initialisation:
      
c     Pour cas 1D et 2D
      FM(1)       = one
      dMdEta(1)   = zero
      ddMddEta(1) = zero
      FL(1)       = one
      dLdZeta(1)  = zero
      ddLddZeta(1)= zero
      
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
      CALL    dersbasisfuns2(Ni(1),Jpqr_patch(1),Nkv_patch(1),xi(1),
     &        Ukv1_patch(:),FN,  dNdxi,  ddNddxi)
      If (dim_patch>1) then
         CALL dersbasisfuns2(Ni(2),Jpqr_patch(2),Nkv_patch(2),xi(2),
     &        Ukv2_patch(:),FM, dMdeta, ddMddEta)
      Endif
      If (dim_patch>2) then
         CALL dersbasisfuns2(Ni(3),Jpqr_patch(3),Nkv_patch(3),xi(3),
     &        Ukv3_patch(:),FL,dLdZeta,ddLddZeta)
      Endif
      
c     - Build numerators and denominators
      NumLoc     = 0
      SumTot     = zero
      SumXi(:)   = zero
      dRdxi(:,:) = zero
      SumXixi(:) = zero
      ddRddxi(:,:) = zero
      Do k = Jpqr_patch(3)+1,1,-1
         Do j = Jpqr_patch(2)+1,1,-1
            Do i = Jpqr_patch(1)+1,1,-1
               NumLoc = NumLoc+1
!     Nurbs functions
               R(NumLoc) = FN(i)*FM(j)*FL(k)*Weight_elem(NumLoc)
               
               SumTot = SumTot + R(NumLoc)

c     Nurbs first derivatives
               dRdxi(NumLoc,1) = 
     &              dNdxi(i)*FM(j)*FL(k)*Weight_elem(NumLoc)
               SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)
c     
               dRdxi(NumLoc,2)=
     &              FN(i)*dMdEta(j)*FL(k)*Weight_elem(NumLoc)
               SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)
c     
               dRdxi(NumLoc,3) = 
     &              FN(i)*FM(j)*dLdZeta(k)*Weight_elem(NumLoc)
               SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
               
c     Nurbs second derivatives
               ddRddxi(NumLoc,1) = 
     &              ddNddxi(i)*FM(j)*FL(k)*Weight_elem(NumLoc)
               SumXixi(1) = SumXixi(1) + ddRddxi(NumLoc,1)
c     
               ddRddxi(NumLoc,2) = 
     &              FN(i)*ddMddEta(j)*FL(k)*Weight_elem(NumLoc)
               SumXixi(2) = SumXixi(2) + ddRddxi(NumLoc,2)
c     
               ddRddxi(NumLoc,3) = 
     &              FN(i)*FM(j)*ddLddZeta(k)*Weight_elem(NumLoc)
               SumXixi(3) = SumXixi(3) + ddRddxi(NumLoc,3)
c     
               ddRddxi(NumLoc,4) = 
     &              dNdxi(i)*dMdEta(j)*FL(k)*Weight_elem(NumLoc)
               SumXixi(4) = SumXixi(4) + ddRddxi(NumLoc,4)
c     
               ddRddxi(NumLoc,5) = 
     &              dNdxi(i)*FM(j)*dLdZeta(k)*Weight_elem(NumLoc)
               SumXixi(5) = SumXixi(5) + ddRddxi(NumLoc,5)
c     
               ddRddxi(NumLoc,6) = 
     &              FN(i)*dMdEta(j)*dLdZeta(k)*Weight_elem(NumLoc)
               SumXixi(6) = SumXixi(6) + ddRddxi(NumLoc,6)
               
            Enddo
         Enddo
      Enddo
      
c     - Divide by denominator to complete definition of fct and deriv.
      SumTot_inv = one/SumTot
      R(:) = R(:)*SumTot_inv
      Do i = 1,dim_patch
         Do NumLoc = 1,nnode_patch
            dRdxi(NumLoc,i)
     &           = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
         Enddo
      Enddo

c     - Complete definition of 2nd derivatives
      Do i = 1,dim_patch
         Do NumLoc = 1,nnode_patch
            ddRddxi(NumLoc,i) = 
     &           ( ddRddxi(NumLoc,i) - R(NumLoc)*SumXixi(i)
     &             - two*dRdxi(NumLoc,i)*SumXi(i) 
     &           )*SumTot_inv
         Enddo
      Enddo
      If (dim_patch>1) then
         Do NumLoc = 1,nnode_patch
            ddRddxi(NumLoc,4) = 
     &           ( ddRddxi(NumLoc,4) - R(NumLoc)*SumXixi(4)
     &           - dRdxi(NumLoc,1)*SumXi(2) - dRdxi(NumLoc,2)*SumXi(1)
     &           )*SumTot_inv
         Enddo
      Endif
      If (dim_patch>2) then
         Do NumLoc = 1,nnode_patch
            ddRddxi(NumLoc,5) = 
     &           ( ddRddxi(NumLoc,5) - R(NumLoc)*SumXixi(5)
     &           - dRdxi(NumLoc,1)*SumXi(3) -dRdxi(NumLoc,3)*SumXi(1)
     &           )*SumTot_inv
c     
            ddRddxi(NumLoc,6) = 
     &           ( ddRddxi(NumLoc,6) - R(NumLoc)*SumXixi(6)
     &           - dRdxi(NumLoc,2)*SumXi(3) -dRdxi(NumLoc,3)*SumXi(2)
     &           )*SumTot_inv
         Enddo
      Endif

      End subroutine evalnurbs_w2ndDerv





c     Evaluation des fonctions Nurbs et des derivees premieres et 
c     secondes en un point xi pour l'element courrant (definit dans le 
c     module nurbspatch)
      
      subroutine evalnurbs_mapping_w2ndDerv(xi,R,dRdxi,ddRddxi)
      
      use parameters
      use embeddedMapping
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Double precision :: xi
      dimension xi(3)
      
c     Output variables :
c     ----------------
      Double precision :: R,dRdxi,ddRddxi
      dimension R(nnode_map),dRdxi(nnode_map,3),
     &     ddRddxi(nnode_map,6)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, SumTot,
     &     SumTot_inv, SumXi, ddNddxi,ddMddEta,ddLddZeta,SumXixi
      dimension FN(Jpqr_map(1)+1),dNdxi(Jpqr_map(1)+1),
     &     FM(Jpqr_map(2)+1),dMdEta(Jpqr_map(2)+1),
     &     FL(Jpqr_map(3)+1),dLdZeta(Jpqr_map(3)+1),SumXi(3),
     &     ddNddxi(Jpqr_map(1)+1),ddMddEta(Jpqr_map(2)+1),
     &     ddLddZeta(Jpqr_map(3)+1),SumXixi(6)
      
      Integer ::  i,j,k, Ni, NumLoc, Na,Nb,Nc
      dimension Ni(dim_map)
      
c     Fin declaration des variables
c
c     ..................................................................
c
c     Initialisation:
      
c     Pour cas 1D et 2D
      FM(1)       = one
      dMdEta(1)   = zero
      ddMddEta(1) = zero
      FL(1)       = one
      dLdZeta(1)  = zero
      ddLddZeta(1)= zero
      
c     Intervalles de l'element courant
      Do i = 1,dim_map
         Ni(i) = Nijk_map(i,current_map_elem)
      End do
      
c     Fin initialisation
c     
c     ..................................................................
c     
c     Calcul:
      
c     - Calculate univariate B-spline function
      CALL    dersbasisfuns2(Ni(1),Jpqr_map(1),Nkv_map(1),xi(1),
     &        Ukv1_map(:),FN,  dNdxi,  ddNddxi)
      If (dim_map>1) then
         CALL dersbasisfuns2(Ni(2),Jpqr_map(2),Nkv_map(2),xi(2),
     &        Ukv2_map(:),FM, dMdeta, ddMddEta)
      Endif
      If (dim_map>2) then
         CALL dersbasisfuns2(Ni(3),Jpqr_map(3),Nkv_map(3),xi(3),
     &        Ukv3_map(:),FL,dLdZeta,ddLddZeta)
      Endif
      
c     - Build numerators and denominators
      NumLoc     = 0
      SumTot     = zero
      SumXi(:)   = zero
      dRdxi(:,:) = zero
      SumXixi(:) = zero
      ddRddxi(:,:) = zero
      Do k = Jpqr_map(3)+1,1,-1
         Do j = Jpqr_map(2)+1,1,-1
            Do i = Jpqr_map(1)+1,1,-1
               NumLoc = NumLoc+1
!     Nurbs functions
               R(NumLoc) = FN(i)*FM(j)*FL(k)*Weight_map_elem(NumLoc)
               
               SumTot = SumTot + R(NumLoc)

c     Nurbs first derivatives
               dRdxi(NumLoc,1) = 
     &              dNdxi(i)*FM(j)*FL(k)*Weight_map_elem(NumLoc)
               SumXi(1) = SumXi(1) + dRdxi(NumLoc,1)
c     
               dRdxi(NumLoc,2)=
     &              FN(i)*dMdEta(j)*FL(k)*Weight_map_elem(NumLoc)
               SumXi(2) = SumXi(2) + dRdxi(NumLoc,2)
c     
               dRdxi(NumLoc,3) = 
     &              FN(i)*FM(j)*dLdZeta(k)*Weight_map_elem(NumLoc)
               SumXi(3) = SumXi(3) + dRdxi(NumLoc,3)
               
c     Nurbs second derivatives
               ddRddxi(NumLoc,1) = 
     &              ddNddxi(i)*FM(j)*FL(k)*Weight_map_elem(NumLoc)
               SumXixi(1) = SumXixi(1) + ddRddxi(NumLoc,1)
c     
               ddRddxi(NumLoc,2) = 
     &              FN(i)*ddMddEta(j)*FL(k)*Weight_map_elem(NumLoc)
               SumXixi(2) = SumXixi(2) + ddRddxi(NumLoc,2)
c     
               ddRddxi(NumLoc,3) = 
     &              FN(i)*FM(j)*ddLddZeta(k)*Weight_map_elem(NumLoc)
               SumXixi(3) = SumXixi(3) + ddRddxi(NumLoc,3)
c     
               ddRddxi(NumLoc,4) = 
     &              dNdxi(i)*dMdEta(j)*FL(k)*Weight_map_elem(NumLoc)
               SumXixi(4) = SumXixi(4) + ddRddxi(NumLoc,4)
c     
               ddRddxi(NumLoc,5) = 
     &              dNdxi(i)*FM(j)*dLdZeta(k)*Weight_map_elem(NumLoc)
               SumXixi(5) = SumXixi(5) + ddRddxi(NumLoc,5)
c     
               ddRddxi(NumLoc,6) = 
     &              FN(i)*dMdEta(j)*dLdZeta(k)*Weight_map_elem(NumLoc)
               SumXixi(6) = SumXixi(6) + ddRddxi(NumLoc,6)
               
            Enddo
         Enddo
      Enddo
      
c     - Divide by denominator to complete definition of fct and deriv.
      SumTot_inv = one/SumTot
      R(:) = R(:)*SumTot_inv
      Do i = 1,dim_map
         Do NumLoc = 1,nnode_map
            dRdxi(NumLoc,i)
     &           = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
         Enddo
      Enddo
      
c     - Complete definition of 2nd derivatives
      Do i = 1,dim_map
         Do NumLoc = 1,nnode_map
            ddRddxi(NumLoc,i) = 
     &           ( ddRddxi(NumLoc,i) - R(NumLoc)*SumXixi(i)
     &             - two*dRdxi(NumLoc,i)*SumXi(i) 
     &           )*SumTot_inv
         Enddo
      Enddo
      If (dim_map>1) then
         Do NumLoc = 1,nnode_map
            ddRddxi(NumLoc,4) = 
     &           ( ddRddxi(NumLoc,4) - R(NumLoc)*SumXixi(4)
     &           - dRdxi(NumLoc,1)*SumXi(2) - dRdxi(NumLoc,2)*SumXi(1)
     &           )*SumTot_inv
         Enddo
      Endif
      If (dim_map>2) then
         Do NumLoc = 1,nnode_map
            ddRddxi(NumLoc,5) = 
     &           ( ddRddxi(NumLoc,5) - R(NumLoc)*SumXixi(5)
     &           - dRdxi(NumLoc,1)*SumXi(3) -dRdxi(NumLoc,3)*SumXi(1)
     &           )*SumTot_inv
c     
            ddRddxi(NumLoc,6) = 
     &           ( ddRddxi(NumLoc,6) - R(NumLoc)*SumXixi(6)
     &           - dRdxi(NumLoc,2)*SumXi(3) -dRdxi(NumLoc,3)*SumXi(2)
     &           )*SumTot_inv
         Enddo
      Endif

      End subroutine evalnurbs_mapping_w2ndDerv
