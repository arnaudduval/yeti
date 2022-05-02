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

C     --
C     Assemblage de la matrice de rigidite geometrique pour le calcul
c     de flambage des structures minces
c     -
c     Ref. Bouclier2014 - PhD Thesis
c          Tiso2006  - PdH Thesis
c          Legay2003 - http://doi.wiley.com/10.1002/nme.728
c     --
      
c     Include iga subroutines
      include "./shap.f"
      include "./Gauss.f"
      include "./operateurs.f"
      include "./dersbasisfuns.f"
      include "./material_lib.f"
      include "./UGEOMAT.f"
      
      include "./shell/UGEOMAT_shell.f"
      include "./shell/curvilinearCoordinates.f"
      include "./shell/nurbsbasisfuns.f"
      include "./plate/dersbasisfuns_4KL.f"
      include "./shell/ComputeMemStrain_shell.f"
      include "./shell/ComputeBndStrain_shell.f"
      
      
      
C     ******************************************************************
      
      
      Subroutine build_GEOStiffMatrix(Ks,COORDS3D,SOL,IEN,nb_elem_patch,
     1     ind_dof_free,Nkv_e,Ukv_e,Nijk_e,weight_e,Jpqr_e,ELT_TYPE,
     2     PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,nb_dof_free,NBINT,
     3     MCRD,nb_patch,nb_elem,nnode,nb_cp,nb_dof_tot)
      
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D, SOL
      dimension COORDS3D(nb_cp,3), SOL(nb_cp,MCRD)
      
      Integer Lknot, Lpatch, Lelement, Lnode
      Parameter (Lelement=10000, Lnode=10000, Lknot=10000, Lpatch=100)
      Double precision, intent(in) :: Ukv_e, weight_e
      Integer, intent(in) :: Nkv_e, Jpqr_e, Nijk_e
      dimension Ukv_e(Lknot,3,Lpatch),
     &     Nkv_e(3,Lpatch),
     &     weight_e(Lelement,Lnode),
     &     Jpqr_e(3),
     &     Nijk_e(Lelement,3)
      
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,
     &     nb_elem_patch,IEN,JPROPS
      dimension MATERIAL_PROPERTIES(nb_patch,2), PROPS(nb_patch,10),
     &     JPROPS(nb_patch), nb_elem_patch(nb_patch), IEN(nb_elem,NNODE)
      
      
!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)
            
      
      
c     Output variables : coefficient diag matrice de masse
c     ----------------
      Double precision, intent(out) :: Ks
      dimension Ks(nb_dof_free,nb_dof_free)
      
      
      
      
c     Local Variables :
c     ---------------
      
      ! Parameters and COMMON variables
      Double precision :: zero, one, two
      Parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      
      Common /NurbsParameter/ Ukv,weight, Nkv, Jpqr,Nijk
      Double precision :: Ukv, weight
      Integer :: Nkv, Jpqr, Nijk
      dimension Ukv(Lknot,3,Lpatch), Nkv(3,Lpatch),
     &     weight(Lelement,Lnode), Jpqr(3), Nijk(Lelement,3)
      
      ! for assembly
      Integer :: num_elem,num_cp,sctr,k,l,kk,ll,k_dof_free,l_dof_free,
     &     tab,k_dof_loc,l_dof_loc
      dimension sctr(NNODE),tab(nb_dof_tot)
      
      ! Variables only used in subroutines
      Double precision :: COORDS_elem, Uelem
      dimension COORDS_elem(MCRD,NNODE), Uelem(MCRD,NNODE)
      Integer :: i,j,JELEM,NumPatch
      
      ! for UGEOMAT.f
      Integer :: JPROPS_patch,NDOFEL
      Double precision :: GEOMATRX,PROPS_patch,MAT_patch
      dimension GEOMATRX(NNODE*MCRD,NNODE*MCRD),PROPS_patch(10),
     &     MAT_patch(2)
      
      
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
      
!     Initialize mass matrix
      NDOFEL = NNODE*MCRD
      Ks(:,:) = zero
      
!     Create tab to link dof_tot <--> dof_free
!     tab(i) = | 0  if dof is blocked
!              | ii if dof is free (ii: index for the reduced matrix)
      tab(:) = 0
      Do i = 1,nb_dof_free
         tab(ind_dof_free(i)) = i
      Enddo
      
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      JELEM = 0
c     Loop on patches
      Do NumPatch = 1,nb_patch
         MAT_patch(:) = MATERIAL_PROPERTIES(NumPatch,:)
         JPROPS_patch = JPROPS(NumPatch)
         PROPS_patch(:) = PROPS(NumPatch,:)
         
c     Loop on elements
         Do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            Do j = 1,MCRD
               Do i = 1,NNODE
                  num_cp = IEN(JELEM,i)
                  COORDS_elem(j,i) = COORDS3D(num_cp,j)
                  Uelem(j,i) = SOL(num_cp,j)
               Enddo
            Enddo
            
            if (ELT_TYPE == 'U1') then
               ! 'Element classique solide'
               call UGEOMAT(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,
     &              COORDS_elem,MAT_patch,TENSOR,Uelem,GEOMATRX)
            elseif (ELT_TYPE == 'U3') then
               ! 'Element coque Kirchhoff-Love'
               call UGEOMAT3(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,
     &              COORDS_elem,TENSOR,MAT_patch,PROPS_patch,
     &              JPROPS_patch,Uelem,GEOMATRX)
            else
               print*,"Element type not available yet (return zero)"
               GEOMATRX(:,:) = zero
            endif
            
c     Assemble DLMMATRX to global matrix A
            k_dof_loc = 0
            l_dof_loc = 0
            sctr = IEN(JELEM,:)
            Do l = 1,NNODE
               ll = MCRD*(sctr(l)-1)
               Do k = 1,NNODE
                  kk = MCRD*(sctr(k)-1)
                  Do j = 1,MCRD
                     l_dof_free = tab(ll+j)
                     l_dof_loc  = MCRD*(l-1)+j
                     Do i = 1,MCRD
                        k_dof_free = tab(kk+i)
                        k_dof_loc  = MCRD*(k-1)+i
                        If (k_dof_free>0 .AND. l_dof_free>0) then
                           Ks(k_dof_free,l_dof_free)
     &                          = Ks(k_dof_free,l_dof_free)
     &                          + GEOMATRX(k_dof_loc,l_dof_loc)
                        Endif
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
      Enddo
      
      
c     Fin Assemblage ...................................................
      
      End subroutine build_GEOStiffMatrix
