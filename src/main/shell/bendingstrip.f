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

      
C     ******************************************************************
      
      subroutine add_bendingstrips(K,COORDS,nb_cp,MCRD,nb_bc,bc_target,
     1     bc_target_nbelem,bc_values,MATERIAL_PROPERTIES,PROPS,JPROPS)
      
      Implicit None
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
      Integer, intent(in) :: nb_cp, MCRD, nb_bc, bc_target,
     &     bc_target_nbelem
      dimension bc_target(1000,1000), bc_target_nbelem(1000)
      
      Double precision, intent(in) :: COORDS, bc_values
      dimension COORDS(nb_cp,3), bc_values(1000,2)
      
      Integer, intent(in) :: JPROPS
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      dimension MATERIAL_PROPERTIES(2), PROPS(10)
      
      
c     Output variables :
c     ----------------
      Double precision, intent(inout) :: K
      dimension K(nb_cp*MCRD,nb_cp*MCRD)
      
      
c     Local variables :
c     ---------------

!     Parameters and COMMON variables
      Double precision :: zero, one, two
      parameter (zero=0.0D0, one=1.0D0, two=2.0D0)
      Integer, parameter :: Lelement=10000, LNODE=10000, Lknot=10000,
     &     Lpatch=100
      
      Common /BndStripParameter/ Ukv_bs, Nkv_bs, Jpqr_bs, Nijk_bs
      Double precision :: Ukv_bs
      Integer :: Nkv_bs, Jpqr_bs, Nijk_bs
      dimension Ukv_bs(Lknot,3,Lpatch), Nkv_bs(3,Lpatch), Jpqr_bs(3),
     &     Nijk_bs(Lelement,3)
      
!     UELMAT_bndstrip.f
      Integer :: NBINT
      Double precision :: COORDS_elem,AMATRX
      dimension COORDS_elem(MCRD,6), AMATRX(MCRD*6,MCRD*6)
      
!     Other
      Integer :: num_bc, num_bndstrip, nb_bndstrip, nb_cpD_eta, JELEM,
     &     num_elem, nb_elem_bndstrip, IEN_bndstrip, i, j, kk, NNODE,
     &     sctr, sctrB, num_knot
      dimension nb_elem_bndstrip(Lpatch), IEN_bndstrip(Lelement,6),
     &     sctr(6), sctrB(MCRD*6)
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation variables communes ................................
      
      Ukv_bs(:,:,:) = 0
      Nkv_bs(:,:)   = 0
      Jpqr_bs(:)    = (/ 2, 1, 0 /)
      Nijk_bs(:,:)  = 0
      
      num_bndstrip = 0
      JELEM = 0
      Do num_bc = 1,nb_bc
         If (bc_values(num_bc,1) == 11.) then
            num_bndstrip = num_bndstrip + 1
            
            nb_cpD_eta = bc_target_nbelem(num_bc)
            
!     --
!     Vecteurs nodaux
!     Direction xi
            Nkv_bs(1,num_bndstrip) = 6
            Ukv_bs(1:Nkv_bs(1,num_bndstrip),1,num_bndstrip)
     &           = (/ zero, zero, zero, one, one, one /)
            
!     Direction eta
            Nkv_bs(2,num_bndstrip) = nb_cpD_eta + 2
            Ukv_bs(1,2,num_bndstrip) = zero
            Ukv_bs(2,2,num_bndstrip) = zero
            Do num_knot = 1,nb_cpD_eta-2
               Ukv_bs(num_knot+2,2,num_bndstrip) = 
     &              float(num_knot)/float((nb_cpD_eta-1))
            Enddo
            Ukv_bs(nb_cpD_eta+1,2,num_bndstrip) = one
            Ukv_bs(nb_cpD_eta+2,2,num_bndstrip) = one
            
!     --
!     Definition des elements
            nb_elem_bndstrip(num_bndstrip) = nb_cpD_eta-1
            Do num_elem = 1,nb_elem_bndstrip(num_bndstrip)
               JELEM = JELEM + 1
               
!     Delimitation des elements
               Nijk_bs(JELEM,1:2) = (/ 3, num_elem+1 /)
               
!     Points de controle associes aux elements
               IEN_bndstrip(JELEM,1) = bc_target(num_bc+600, num_elem+1)
               IEN_bndstrip(JELEM,2) = bc_target(num_bc, num_elem+1)
               IEN_bndstrip(JELEM,3) = bc_target(num_bc+500, num_elem+1)
               IEN_bndstrip(JELEM,4) = bc_target(num_bc+600, num_elem)
               IEN_bndstrip(JELEM,5) = bc_target(num_bc, num_elem)
               IEN_bndstrip(JELEM,6) = bc_target(num_bc+500, num_elem)
               
            Enddo
         Endif
      Enddo
      
      nb_bndstrip = num_bndstrip
      
c     Fin variables communes  ..........................................
c     
c     Ajout des bending strip a la matrice K ...........................
            
      NBINT = 9
      NNODE = 6
      JELEM = 0
      Do num_bndstrip = 1,nb_bndstrip
         Do num_elem = 1,nb_elem_bndstrip(num_bndstrip)
            JELEM = JELEM + 1
            
!     Get element control point coordinates
            Do i = 1,NNODE
               Do j = 1,MCRD
                  COORDS_elem(j,i) = COORDS(IEN_bndstrip(JELEM,i),j)
               Enddo
            Enddo
            
!     Compute elementary stiffness matrix
            AMATRX(:,:) = zero
            call UELMATbndstrip(MCRD,NNODE,JELEM,NBINT,num_bndstrip,
     &           COORDS_elem, MATERIAL_PROPERTIES,PROPS,JPROPS,AMATRX)
            
!     Compute ddl infos
            sctr(:) = IEN_bndstrip(JELEM,:)
            Do kk = 1,MCRD
               j = 1
               Do i = kk,MCRD*NNODE,MCRD
                  sctrB(i) = MCRD*(sctr(j)-1) + kk
                  j = j+1
               Enddo
            Enddo
            
!     Update global stiffness matrix
            Do i = 1,MCRD*NNODE
               Do j = 1,MCRD*NNODE
                  K(sctrB(i),sctrB(j)) = K(sctrB(i),sctrB(j)) 
     &                 + AMATRX(i,j)
               Enddo
            Enddo
         Enddo
      Enddo
      
      End subroutine add_bendingstrips
