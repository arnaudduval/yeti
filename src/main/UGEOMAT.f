!! Copyright 2017-2018 Thibaut Hirschler

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

c     Construction elementaire de la matrice de raideur geometrique pour
c     pour les elements classiques solides
c     
C     ******************************************************************
      
      Subroutine UGEOMAT(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,COORDS,
     1     MATERIAL_PROPERTIES,TENSOR,Uelem,GEOMATRX)
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NumPatch,NBINT
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, Uelem
      dimension COORDS(MCRD,NNODE), MATERIAL_PROPERTIES(2),
     &     Uelem(MCRD,NNODE)
      Character(len=*), intent(in) :: TENSOR
      
c     Output variables : matrice de masse elementaire
c     ----------------
      Double precision, intent(out) :: GEOMATRX
      dimension GEOMATRX(NDOFEL,NDOFEL)
      
      
      
c     Local variables :
c     ---------------
      
!     indices
      Integer :: n,k,nodi,nodj,dof_i,dof_j
      
!     Parameters
      Double precision :: zero, one, two
      Parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      
!     for material_lib.f
      Integer :: ntens
      Double precision :: ddsdde
      dimension ddsdde(2*MCRD,2*MCRD)
      
!     for gauss.f
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(MCRD+1,NBINT)
      
!     for shap.f
      Double precision :: R,dRdx,dRidx,dRjdx, DetJac, PtGauss
      dimension R(NNODE),dRdx(MCRD,NNODE),dRidx(MCRD),dRjdx(MCRD),
     &     PtGauss(MCRD)
      
!     internal procedure
      Double precision :: strain,stress,Ui,dvol,temp1,temp2,temp3,
     &     temp41,temp42,temp51,temp53,temp62,temp63,res
      dimension strain(MCRD*MCRD),stress(MCRD*MCRD),Ui(MCRD)

      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Initialize GEOMATRX to zero
      GEOMATRX(:,:) = zero

!     Number of integration point in each direction
      if (MCRD==2) then
         NbPtInt = int(NBINT**(1.0/2.0))
      else
         NbPtInt = int(NBINT**(1.0/3.0))
      endif
      
!     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)
      
!     Material behaviour
      ntens = 2*MCRD
      call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde)
      
c     Fin initialisaton ................................................
c     
c     
c     
c     Fin calcul .......................................................
      
!     Loop on integration points
      Do n = 1,NBINT
!     Computing Nurbs basis functions at integration points
         dvol = zero
         
!     NURBS basis functions and derivatives
         call shap(COORDS,dRdx,NNODE,R,GaussPdsCoord(2:,n),DetJac,
     &        NDOFEL,MCRD,JELEM,NumPatch)
         dvol = GaussPdsCoord(1,n)*DetJac
         
!     Compute strain field at current gauss point
         strain(:) = zero
         Do nodi = 1,NNODE
            dRidx(:) = dRdx(:,nodi)
            Ui(:) = Uelem(:,nodi)
            
            strain(1) = strain(1) + dRidx(1)*Ui(1)
            strain(2) = strain(2) + dRidx(2)*Ui(2)
            strain(4) = strain(4) + dRidx(2)*Ui(1)+dRidx(1)*Ui(2)
            If (MCRD==3) then
               strain(3) = strain(3) + dRidx(3)*Ui(3)
               strain(5) = strain(5) + dRidx(3)*Ui(1)+dRidx(1)*Ui(3)
               strain(6) = strain(6) + dRidx(3)*Ui(2)+dRidx(2)*Ui(3)
            Endif
         Enddo
         
!     Get stress field at current gauss point
         call MulVect(ddsdde(:MCRD,:MCRD), strain(:MCRD), stress(:MCRD), 
     &        MCRD,MCRD)
         stress(4) = ddsdde(4,4)*strain(4)
         stress(5) = ddsdde(5,5)*strain(5)
         stress(6) = ddsdde(6,6)*strain(6)
         
!     Computation and Assembly of CGEOMATRIX
         stress(:) = stress(:)*dvol
         Do nodj = 1,NNODE
            dRjdx(:) = dRdx(:,nodj)
            
            temp1 = stress(1)*dRjdx(1)
            temp2 = stress(2)*dRjdx(2)
            temp3 = stress(3)*dRjdx(3)
            temp41= stress(4)*dRjdx(1); temp42= stress(4)*dRjdx(2)
            temp51= stress(5)*dRjdx(1); temp53= stress(5)*dRjdx(3)
            temp62= stress(6)*dRjdx(2); temp63= stress(6)*dRjdx(3)
            
            dof_j = (nodj-1)*MCRD
            Do nodi = nodj,NNODE
               dRidx(:) = dRdx(:,nodi)
               
               res = temp1*dRidx(1) + temp2*dRidx(2) + temp3*dRidx(3)
     &              + temp41*dRidx(2) + temp42*dRidx(1)
     &              + temp51*dRidx(3) + temp53*dRidx(1)
     &              + temp62*dRidx(3) + temp63*dRidx(2)
               
               dof_i = (nodi-1)*MCRD
               Do k=1,MCRD
                  GEOMATRX(dof_i+k,dof_j+k)
     &                 = GEOMATRX(dof_i+k,dof_j+k) + res
               Enddo
            Enddo
         Enddo
      Enddo
!     End loop on integration points
      
c     Symmetry
      Do dof_j = 2,NDOFEL
         Do dof_i = 1,dof_j-1
            GEOMATRX(dof_i,dof_j) = GEOMATRX(dof_j,dof_i)
         Enddo
      Enddo
      
      
c     Fin calcul .......................................................
      
      End Subroutine UGEOMAT
