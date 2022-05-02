!! Copyright 2018 Thibaut Hirschler

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
      
      
      subroutine bandstorageInfo(nb_diag,NUMDiag2Dof,NUMDof2Diag,
     1     IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,PROPS,
     2     JPROPS,TENSOR,ind_dof_free,nb_dof_free,MCRD,NBINT,nb_patch,
     3     nb_elem,nnode,nb_dof_tot)
      
      use parameters
      use nurbspatch
      
      Implicit none
            
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      
      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: PROPS
      Integer,          intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN,nb_elem_patch, JPROPS
      dimension PROPS(:),NNODE(nb_patch),IEN(:),nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),NBINT(nb_patch)
      
!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)
            
      
c     Output variables :
c     ----------------
      Integer, intent(out):: nb_diag,NUMDiag2Dof,NUMDof2Diag
      dimension NUMDiag2Dof(nb_dof_tot),NUMDof2Diag(nb_dof_tot)
      
      
      
c     Local variables :
c     ---------------
            
!     Loops
      Integer :: NumPatch,num_elem,sctr,cpi,i,j,kk
      dimension sctr(MAXVAL(NNODE))
      
!     For test
      Logical :: MASK
      dimension MASK(nb_dof_tot)
      Integer :: diag,range,num_diag
      dimension diag(nb_dof_tot), range(nb_dof_tot)
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Start ............................................................
      
c     Donnee pour le stockage en bandes
      
      diag(:) = 0
      Do NumPatch = 1,nb_patch
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         Do num_elem = 1,nb_elem_patch(NumPatch)
            sctr(:nnode_patch) = IEN_patch(:,num_elem)
            do i = 1,nnode_patch
               cpi  = sctr(i)
               do j = i+1,nnode_patch
                  num_diag = abs(sctr(j) - cpi)*MCRD
                  
                  do kk = num_diag-MCRD+2,num_diag+MCRD
                     diag(kk) = kk
                  enddo
               enddo
            enddo
         Enddo
      Enddo
      Do kk = 1,MCRD
         diag(kk) = kk
      Enddo
      
      MASK    = diag.NE.0
      nb_diag = COUNT(MASK)
      
      range(:nb_diag) = (/(i, i=1,nb_diag, 1)/)
      NUMDiag2Dof(:nb_diag) = PACK(diag, MASK)
      NUMDof2Diag(:) = 0
      NUMDof2Diag(:) = UNPACK(range(:nb_diag), MASK, NUMDof2Diag(:))
      
c     ..................................................................
      
      end subroutine bandstorageInfo
