!! Copyright 2016-2018 Thibaut Hirschler

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

c     Include iga subroutines
      include "./shap.f"
      include "./Gauss.f"
      include "./operateurs.f"
      include "./dersbasisfuns.f"
      include "./DLMMAT.f"
      include "./shell/nurbsbasisfuns.f"
      include "./plate/dersbasisfuns_4KL.f"
      
      
      
C     ******************************************************************
      
      
      Subroutine build_DLMMatrix(A, COORDS3D,IEN,nb_elem_patch, Nkv_e,
     1     Ukv_e,Nijk_e,weight_e,Jpqr_e,ELT_TYPE,MCRD,NBINT,nb_patch,
     2     nb_elem,nnode,nb_cp)
      
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(nb_cp,3)
      
      Integer Lknot, Lpatch, Lelement, Lnode
      Parameter (Lelement=10000, Lnode=10000, Lknot=10000, Lpatch=100)
      Double precision, intent(in) :: Ukv_e, weight_e
      Integer, intent(in) :: Nkv_e, Jpqr_e, Nijk_e
      dimension Ukv_e(Lknot,3,Lpatch),
     &     Nkv_e(3,Lpatch),
     &     weight_e(Lelement,Lnode),
     &     Jpqr_e(3),
     &     Nijk_e(Lelement,3)
      
      
!     Parches and Elements
      Integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,
     &     nb_elem_patch, IEN
      dimension nb_elem_patch(nb_patch), IEN(nb_elem,NNODE)
      Character(len=*), intent(in) :: ELT_TYPE
            
      
      
c     Output variables : coefficient diag matrice de masse
c     ----------------
      Double precision, intent(out) :: A
      dimension A(nb_cp)
      
      
      
      
c     Local Variables :
c     ---------------
      
!     Parameters and COMMON variables
      Double precision, parameter :: zero=0.0D0, one=1.0D0, two=2.0D0

      Common /NurbsParameter/ Ukv,weight, Nkv, Jpqr,Nijk      
      Double precision :: Ukv, weight
      Integer :: Nkv, Jpqr, Nijk
      dimension Ukv(Lknot,3,Lpatch), Nkv(3,Lpatch),
     &     weight(Lelement,Lnode), Jpqr(3), Nijk(Lelement,3)
      
      ! Diagonally-Lumped Mapping Matrix
      Integer :: num_elem,sctr
      dimension sctr(NNODE)
      
      ! Variables only used in subroutines
      Double precision :: COORDS_elem
      dimension COORDS_elem(MCRD, NNODE)
      integer i,j,JELEM,NumPatch
      
      ! for DLMMAT.f
      Double precision :: DLMMATRX
      dimension DLMMATRX(NNODE)
            
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Assign common variables
      Ukv  = Ukv_e
      Nkv  = Nkv_e
      Jpqr = Jpqr_e
      Nijk = Nijk_e
      weight = weight_e  
      
!     Initialize A to zero
      A(:) = zero
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      JELEM = 0
!     Loop on patches
      Do NumPatch = 1,nb_patch
!     Loop on elements
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
            JELEM = JELEM + 1
            
            Do i = 1,NNODE
               Do j = 1,MCRD 
                  COORDS_elem(j,i) = COORDS3D(IEN(JELEM,i),j)
               Enddo
            Enddo
            
            
!     Build elementary Matrix
            call DLMMAT(MCRD,NNODE,JELEM,NBINT,NumPatch,COORDS_elem,
     &           ELT_TYPE,DLMMATRX)
            
            
!     Assemble DLMMATRX to global matrix A    
            sctr = IEN(JELEM,:)
            Do i = 1,NNODE
               A(sctr(i)) = A(sctr(i)) + DLMMATRX(i)
            Enddo
            
         Enddo
      Enddo
      
      
c     Fin Assemblage ...................................................
      
      End subroutine build_DLMMatrix
