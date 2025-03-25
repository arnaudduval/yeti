!! Copyright 2011 Florian Maurin
!! Copyright 2016-2018 Thibaut Hirschler
!! Copyright 2020 Arnaud Duval

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

      subroutine shapPress (R,Vect,DetJac, COORDS,PtGauss,MCRD,KNumFace,
     1     KTypeDload)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
c     Input arguments :
c     ---------------
      Integer :: MCRD,KNumFace,KTypeDload
      
      Double precision :: COORDS,PtGauss
      dimension COORDS(MCRD,nnode_patch),PtGauss(MCRD)
      
c     Output variables :
c     ----------------
      Double precision :: R,Vect,DetJac
      dimension R(nnode_patch),Vect(MCRD)
      
      
c     Local variables :
c     ---------------
      Double precision :: FN,FM,FL, dNdxi,dMdEta,dLdZeta, dxdxi,
     &     dXidtildexi, AJmat, xi, dRdxi, Detdxdxi, SumTot,
     &     SumTot_inv, SumXi
      dimension FN(Jpqr_patch(1)+1),dNdxi(Jpqr_patch(1)+1),
     &     FM(Jpqr_patch(2)+1),dMdEta(Jpqr_patch(2)+1),
     &     FL(Jpqr_patch(3)+1),dLdZeta(Jpqr_patch(3)+1),
     &     dxdxi(MCRD,MCRD),dXidtildexi(MCRD,MCRD),AJmat(MCRD,MCRD),
     &     xi(dim_patch),dRdxi(nnode_patch,3),SumXi(3)
      
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
      
c     Initialisation des differentes matrices
      dxdxi(:,:) = zero
      dXidtildexi(:,:) = zero
      AJmat(:,:) = zero
      
c     Fin initialisation
c     
c     ..................................................................
c     
c     Calcul:

c     - Calculate paraetrique coordinates from parents element
      Do i = 1,dim_patch
         xi(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &        +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
      End do
      
c     - Calculate univariate B-spline function
      CALL dersbasisfuns(Ni(1),Jpqr_patch(1),Nkv_patch(1),xi(1),
     &     Ukv1_patch(:),FN,dNdxi)
      CALL dersbasisfuns(Ni(2),Jpqr_patch(2),Nkv_patch(2),xi(2),
     &     Ukv2_patch(:),FM,dMdeta)
      If (dim_patch==3) then
         CALL dersbasisfuns(Ni(3),Jpqr_patch(3),Nkv_patch(3),xi(3),
     &        Ukv3_patch(:),FL,dLdZeta)
      Endif
      
c     - Build numerators and denominators
      NumLoc   = 0
      SumTot   = zero
      SumXi(:) = zero
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
         Do i = 1,MCRD
            dRdxi(NumLoc,i)
     &           = (dRdxi(NumLoc,i)-R(NumLoc)*SumXi(i))*SumTot_inv
         Enddo
      Enddo
      
      
c     - Gradient of mapping from parameter space to physical space
      Do NumLoc= 1,nnode_patch
         Do Na = 1,MCRD
            Do Nb = 1,MCRD
               dxdxi(Na,Nb) = dxdxi(Na,Nb)
     &              + COORDS(Na,NumLoc)*dRdxi(NumLoc,Nb)
            Enddo
         Enddo
      Enddo
      
c     - Gradient of mapping from parent element to parameter space
      Do i = 1,dim_patch
         dXidtildexi(i,i) = 0.5d0*( Ukv_elem(2,i) - Ukv_elem(1,i) )
      Enddo
      
      Do Na = 1,MCRD
         Do Nb = 1,MCRD
            Do Nc = 1,MCRD
               AJmat(Na,Nb)=AJmat(Na,Nb)+dxdxi(Na,Nc)*dXidtildexi(Nc,Nb)
            Enddo
         Enddo
      Enddo
      
      
c     Determination de la surface �l�mentaire en fct du num de la face
      SELECT CASE (KNumFace+MCRD*10)
      case(21,22)               !Cot� 1 et 2
         DetJac = sqrt(AJmat(1,2)**two+AJmat(2,2)**two)
      case(23,24)               !Cot� 3 et 4
         DetJac = sqrt(AJmat(1,1)**two+AJmat(2,1)**two)
      case(31,32)               !Face 1 et 2
         call SurfElem(AJmat(:,2),AJmat(:,3),DetJac)
      case(33,34)               !Face 3 et 4
         call SurfElem(AJmat(:,3),AJmat(:,1),DetJac)
      case(35,36)               !Face 5 et 6
         call SurfElem(AJmat(:,1),AJmat(:,2),DetJac)
      END SELECT

      

c     Determination du vecteur normalise
      SELECT CASE (KTypeDload+MCRD*10)
      case(21)                  !Cas contrainte direction x
         Vect(1)=One
         Vect(2)=zero
      case(22)                  !Cas contrainte direction y
         Vect(1)=zero
         Vect(2)=one
      case(31)                  !Cas contrainte direction x
         Vect(1)=one
         Vect(2)=zero         
         Vect(3)=zero
      case(32)                  !Cas contrainte direction y
         Vect(1)=zero
         Vect(2)=one        
         Vect(3)=zero 
      case(33)                  !Cas contrainte direction z
         Vect(1)=zero
         Vect(2)=zero         
         Vect(3)=one    
      case(20)                  !Cas des pressions normales en 2D
         SELECT CASE (KNumFace)
         case(1)                ! Side 1     
            Vect(1)= AJmat(2,2)/DetJac
            Vect(2)=-AJmat(1,2)/DetJac 
         case(2)                ! Side 2
            Vect(1)=-AJmat(2,2)/DetJac
            Vect(2)= AJmat(1,2)/DetJac
         case(3)                ! Side 3
            Vect(1)=-AJmat(2,1)/DetJac
            Vect(2)= AJmat(1,1)/DetJac
         case(4)                ! Side 4
            Vect(1)= AJmat(2,1)/DetJac
            Vect(2)=-AJmat(1,1)/DetJac
         END SELECT
      case(29)                  !Cas des pressions tengencielles en 2D
         SELECT CASE (KNumFace)
         case(1)                !Cote 1
            Vect(1)= AJmat(2,2)/DetJac
            Vect(2)= AJmat(1,2)/DetJac 
         case(2)                !Cote 2
            Vect(1)=-AJmat(2,2)/DetJac
            Vect(2)=-AJmat(1,2)/DetJac
         case(3)                !Cote 3
            Vect(1)= AJmat(2,1)/DetJac
            Vect(2)= AJmat(1,1)/DetJac
         case(4)                !Cote 4
            Vect(1)=-AJmat(2,1)/DetJac
            Vect(2)=-AJmat(1,1)/DetJac
         END SELECT
      case(30)                  !Cas des pressions normales en 3D
         SELECT CASE (KNumFace)
         case(1)                !Face 1
            call VectNormNorm(Vect,AJmat(:,2),AJmat(:,3),DetJac)
         case(2)                !Face 2
            call VectNormNorm(Vect,AJmat(:,3),AJmat(:,2),DetJac)
         case(3)                !Face 3
            call VectNormNorm(Vect,AJmat(:,3),AJmat(:,1),DetJac)
         case(4)                !Face 4
            call VectNormNorm(Vect,AJmat(:,1),AJmat(:,3),DetJac)
         case(5)                !Face 5
            call VectNormNorm(Vect,AJmat(:,1),AJmat(:,2),DetJac)
         case(6)                !Face 6
            call VectNormNorm(Vect,AJmat(:,2),AJmat(:,1),DetJac)
         END SELECT     
      case(34)                  !Pression non uniforme normale
         SELECT CASE (KNumFace)
         case(1)                !Face 1
            call VectNormNorm(Vect,AJmat(:,2),AJmat(:,3),DetJac)
         case(2)                !Face 2
            call VectNormNorm(Vect,AJmat(:,3),AJmat(:,2),DetJac)
         case(3)                !Face 3
            call VectNormNorm(Vect,AJmat(:,3),AJmat(:,1),DetJac)
         case(4)                !Face 4
            call VectNormNorm(Vect,AJmat(:,1),AJmat(:,3),DetJac)
         case(5)                !Face 5
            call VectNormNorm(Vect,AJmat(:,1),AJmat(:,2),DetJac)
         case(6)                !Face 6
            call VectNormNorm(Vect,AJmat(:,2),AJmat(:,1),DetJac)
         END SELECT     
		 
      END SELECT
      
      end SUBROUTINE shapPress
      
      
      
