!! Copyright 2017-2019 Thibaut Hirschler

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

c     Construction elementaire de la matrice de masse 
C     --> "consistent mass matrix"
      
C     ******************************************************************
      
      Subroutine UMASSMAT(NDOFEL,MCRD,NNODE,NBINT,COORDS,ELT_TYPE,
     1     DENSITY,PROPS,JPROPS,CMASSMATRX)
      
      use parameters
!     use gauss
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: NDOFEL, MCRD, NNODE, NBINT, JPROPS
      Double precision, intent(in) :: COORDS, DENSITY, PROPS
      dimension COORDS(MCRD,NNODE), PROPS(JPROPS)
      
      Character(len=*), intent(in) :: ELT_TYPE
      
      
c     Output variables : matrice de masse elementaire
c     ----------------
      Double precision, intent(out) :: CMASSMATRX
      dimension CMASSMATRX(NDOFEL,NDOFEL)
      
      
      
c     Local variables :
c     ---------------
      
      Integer :: i,j,n,k,numLoc,nodi,nodj,dof_i,dof_j
      Double precision :: dvol,temp1,temp2
      
!     for gauss.f
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(MCRD+1,NBINT)
      
!     for shap.f
      Double precision :: R, dRdx, DetJac
      dimension R(NNODE), dRdx(MCRD,NNODE)
      
!     for nurbsbasis.f
      Double precision :: A1,A2,A3, A00, dRdxi, ddRddxi, Area,Area2, h,
     &     h3_12, m00,mij, A3tA3, Ui,Uj, A1xA3,A2xA3
      dimension A1(MCRD), A2(MCRD), A3(MCRD), dRdxi(NNODE,2), 
     &     ddRddxi(NNODE,3), A3tA3(MCRD,MCRD), Ui(MCRD), Uj(MCRD),
     &     A1xA3(MCRD), A2xA3(MCRD)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Initialize CMASSMATRX to zero
      CMASSMATRX(:,:) = zero
      h = zero; h3_12 = zero
      
!     Number of integration point in each direction
      if (MCRD==2 .OR. ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         NbPtInt = int(NBINT**(1.0/2.0))
         if (NbPtInt*NbPtInt<NBINT) NbPtInt = NbPtInt + 1
      else if (MCRD==3) then
         NbPtInt = int(NBINT**(1.0/3.0))
         if (NbPtInt**3<NBINT)      NbPtInt = NbPtInt + 1
      endif
      
      
!     Defining Gauss points coordinates and Gauss weights
      If (ELT_TYPE=='U1') then
         call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)
      elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         call Gauss(NbPtInt, 2, GaussPdsCoord(:MCRD,:), 0)
         h = PROPS(2)
         h3_12 = h*h*h/12.d0
      Endif
      
      
c     Fin initialisaton ................................................
c     
c     
c     
c     Fin calcul .......................................................
      
!     Loop on integration points
      Do n = 1,NBINT
!     Computing Nurbs basis functions at integration points
         dvol = zero
         
c     --
c     SOLID element
         If (ELT_TYPE=='U1') then
            CALL shap(dRdx,R,DetJac,COORDS,GaussPdsCoord(2:,n),MCRD)
            dvol = GaussPdsCoord(1,n)*DetJac
         
!     Assembling CMASSMATRIX
            dvol = DENSITY*dvol
            Do nodj = 1,NNODE
               temp1 = dvol*R(nodj)
               dof_j = (nodj-1)*MCRD
               Do nodi = nodj,NNODE
                  temp2 = temp1*R(nodi)
                  dof_i = (nodi-1)*MCRD
                  Do k=1,MCRD
                     CMASSMATRX(dof_i+k,dof_j+k)
     &                    = CMASSMATRX(dof_i+k,dof_j+k) + temp2
                  Enddo
               Enddo
            Enddo
            
            
c     --
c     Kirchhoff-Love PLATE and SHELL
         elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
            call nurbsbasis(R,dRdxi,ddRddxi,DetJac,
     &           GaussPdsCoord(2:MCRD,n))
            
c     Computing Curvilinear Coordinate objects
            A1(:) = zero; A2(:) = zero
            Do numLoc = 1,NNODE
               A1(:) = A1(:) + dRdxi(numLoc,1)*COORDS(:,numLoc)
               A2(:) = A2(:) + dRdxi(numLoc,2)*COORDS(:,numLoc)
            Enddo
            call cross(A1,A2, A3)
            Area2= A3(1)*A3(1) + A3(2)*A3(2) + A3(3)*A3(3)
            Area = sqrt( Area2 )
            A3(:) = A3(:)/Area
            
            call cross(A1,A3, A1xA3)
            call cross(A2,A3, A2xA3)
                        
            A3tA3(:,:) = zero
            Do j = 1,MCRD
               Do i = 1,MCRD
                  A3tA3(i,j) = A3(i)*A3(j)
               Enddo
            Enddo
            
c     Computation of redundant quantities
            dvol = DENSITY * GaussPdsCoord(1,n)*DetJac*Area
            A00 = h * dvol
            dvol= h3_12/Area2 * dvol
            
c     Assembling
            Do nodj = 1,NNODE
               Uj(:)= dRdxi(nodj,2)*A1xA3(:) - dRdxi(nodj,1)*A2xA3(:)
               Do nodi = nodj,NNODE
                  m00  = A00 * R(nodi)*R(nodj)
                  Ui(:)= dRdxi(nodi,2)*A1xA3(:) - dRdxi(nodi,1)*A2xA3(:)
                  
                  call dot(Ui,Uj, mij)
                  mij = mij * dvol
                  Do j = 1,MCRD
                     dof_j = (nodj-1)*MCRD + j
                     dof_i = (nodi-1)*MCRD + j
                     CMASSMATRX(dof_i,dof_j)= CMASSMATRX(dof_i,dof_j)
     &                    + m00
                     Do i = 1,MCRD
                        dof_i = (nodi-1)*MCRD + i
                        CMASSMATRX(dof_i,dof_j)= CMASSMATRX(dof_i,dof_j)
     &                       + mij * A3tA3(i,j)
                     Enddo
                  Enddo
               Enddo
            Enddo
         Endif
         
                  
      Enddo
!     End loop on integration points
      
c     Symmetry
      Do dof_j = 2,NDOFEL
         Do dof_i = 1,dof_j-1
            CMASSMATRX(dof_i,dof_j) = CMASSMATRX(dof_j,dof_i)
         Enddo
      Enddo
      
      
c     Fin calcul .......................................................
      
      End Subroutine UMASSMAT



















c     Construction elementaire de la matrice de masse 
C     --> "consistent mass matrix"
c     Stockage sans construire la matrice : renvoie les matrices 3x3 de
c     chaque couple de points de controle (partie triangulaire 
c     supperieure uniquement)
      
C     ******************************************************************
      
      Subroutine UMASSMAT_byCP(NDOFEL,MCRD,NNODE,NBINT,COORDS,ELT_TYPE,
     1     DENSITY,PROPS,JPROPS,CMASSMATRX)
      
      use parameters
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: NDOFEL, MCRD, NNODE, NBINT, JPROPS
      Double precision, intent(in) :: COORDS, DENSITY, PROPS
      dimension COORDS(MCRD,NNODE), PROPS(JPROPS)
      
      Character(len=*), intent(in) :: ELT_TYPE
      
      
c     Output variables : matrice de masse elementaire
c     ----------------
      Double precision, intent(out) :: CMASSMATRX
      dimension CMASSMATRX(MCRD,MCRD,NNODE*(NNODE+1)/2)
      
      
      
c     Local variables :
c     ---------------
      
      Integer :: i,j,n,k,numLoc,nodi,nodj,dof_i,dof_j,count
      Double precision :: dvol,temp1,temp2
      
!     for gauss.f
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(MCRD+1,NBINT)
      
!     for shap.f
      Double precision :: R, dRdx, DetJac
      dimension R(NNODE), dRdx(MCRD,NNODE)
      
!     for nurbsbasis.f
      Double precision :: A1,A2,A3, A00, dRdxi, ddRddxi, Area,Area2, h,
     &     h3_12, m00,mij, A3tA3, Ui,Uj, A1xA3,A2xA3
      dimension A1(MCRD), A2(MCRD), A3(MCRD), dRdxi(NNODE,2), 
     &     ddRddxi(NNODE,3), A3tA3(MCRD,MCRD), Ui(MCRD), Uj(MCRD),
     &     A1xA3(MCRD), A2xA3(MCRD)
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Initialize CMASSMATRX to zero
      CMASSMATRX(:,:,:) = zero
      h = zero; h3_12 = zero
      
!     Number of integration point in each direction
      if (MCRD==2 .OR. ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         NbPtInt = int(NBINT**(1.0/2.0))
         if (NbPtInt*NbPtInt<NBINT) NbPtInt = NbPtInt + 1
      else if (MCRD==3) then
         NbPtInt = int(NBINT**(1.0/3.0))
         if (NbPtInt**3<NBINT)      NbPtInt = NbPtInt + 1
      endif
      
      
!     Defining Gauss points coordinates and Gauss weights
      If (ELT_TYPE=='U1') then
         call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)
      elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
         call Gauss(NbPtInt, 2, GaussPdsCoord(:MCRD,:), 0)
         h = PROPS(2)
         h3_12 = h*h*h/12.d0
      Endif
      
      
c     Fin initialisaton ................................................
c     
c     
c     
c     Fin calcul .......................................................
      
!     Loop on integration points
      Do n = 1,NBINT
!     Computing Nurbs basis functions at integration points
         dvol = zero
         
c     --
c     SOLID element
         If (ELT_TYPE=='U1') then
            CALL shap(dRdx,R,DetJac,COORDS,GaussPdsCoord(2:,n),MCRD)
            dvol = GaussPdsCoord(1,n)*DetJac
         
!     Assembling CMASSMATRIX
            dvol = DENSITY*dvol
            count= 1
            Do nodj = 1,NNODE
               temp1 = dvol*R(nodj)
               Do nodi = 1,nodj
                  temp2 = temp1*R(nodi)
                  Do k=1,MCRD
                     CMASSMATRX(k,k,count)
     &                    = CMASSMATRX(k,k,count) + temp2
                  Enddo
                  count = count+1
               Enddo
            Enddo
            
            
c     --
c     Kirchhoff-Love PLATE and SHELL
         elseif (ELT_TYPE=='U2' .OR. ELT_TYPE=='U3') then
            call nurbsbasis(R,dRdxi,ddRddxi,DetJac,
     &           GaussPdsCoord(2:MCRD,n))
            
c     Computing Curvilinear Coordinate objects
            A1(:) = zero; A2(:) = zero
            Do numLoc = 1,NNODE
               A1(:) = A1(:) + dRdxi(numLoc,1)*COORDS(:,numLoc)
               A2(:) = A2(:) + dRdxi(numLoc,2)*COORDS(:,numLoc)
            Enddo
            call cross(A1,A2, A3)
            Area2= A3(1)*A3(1) + A3(2)*A3(2) + A3(3)*A3(3)
            Area = sqrt( Area2 )
            A3(:) = A3(:)/Area
            
            call cross(A1,A3, A1xA3)
            call cross(A2,A3, A2xA3)
                        
            A3tA3(:,:) = zero
            Do j = 1,MCRD
               Do i = 1,MCRD
                  A3tA3(i,j) = A3(i)*A3(j)
               Enddo
            Enddo
            
c     Computation of redundant quantities
            dvol = DENSITY * GaussPdsCoord(1,n)*DetJac*Area
            A00 = h * dvol
            dvol= h3_12/Area2 * dvol
            
c     Assembling
            count = 1
            Do nodj = 1,NNODE
               Uj(:)= dRdxi(nodj,2)*A1xA3(:) - dRdxi(nodj,1)*A2xA3(:)
               Do nodi = 1,nodj
                  m00  = A00 * R(nodi)*R(nodj)
                  Ui(:)= dRdxi(nodi,2)*A1xA3(:) - dRdxi(nodi,1)*A2xA3(:)
                  
                  call dot(Ui,Uj, mij)
                  mij = mij * dvol
                  Do j = 1,MCRD
                     CMASSMATRX(j,j,count) = CMASSMATRX(j,j,count) + m00
                     Do i = 1,MCRD
                        CMASSMATRX(i,j,count) = CMASSMATRX(i,j,count)
     &                       + mij * A3tA3(i,j)
                     Enddo
                  Enddo
                  count = count+1
               Enddo
            Enddo
         Endif
         
                  
      Enddo
!     End loop on integration points
           
      
c     Fin calcul .......................................................
      
      End Subroutine UMASSMAT_byCP

