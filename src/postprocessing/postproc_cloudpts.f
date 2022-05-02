!! Copyright 2019 Thibaut Hirschler

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

c     POST-PROCESSING : evaluate NURBS at multiple points
c      - output.txt file
      
      subroutine generateCloudPtsTXT(FILENAME,ActivePatch,
     1     nb_pts, COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,MCRD,
     3     NBINT,nnode,nb_patch,nb_elem,nb_cp)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer,          intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer,          intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN,nb_elem_patch,JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),
     &     PROPS(:),
     &     NNODE(nb_patch),
     &     IEN(:),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),
     &     NBINT(nb_patch)
      
      
!     Output INFOS
      Character(len=*), intent(in) :: FILENAME
      Integer,          intent(in) :: nb_pts,ActivePatch
      dimension ActivePatch(nb_patch)
      
      
c     Local variables :
c     ---------------
      
!     Extract infos
      Integer          :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: COORDS_elem
      dimension COORDS_elem(MCRD,MAXVAL(NNODE))
      
!     Compute quantities
      Integer          :: nb_xi,nb_eta,nb_zeta,i_xi,i_eta,i_zeta
      Double precision :: h,XIb,XI,R,dRdxi,ddRddxi,pts
      dimension XIb(3),XI(3),R(MAXVAL(NNODE)), dRdxi(MAXVAL(NNODE),3),
     &     pts(MCRD)
      
      Integer :: i,numPatch,numel
      

      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Start  ...........................................................
      
      
c     Retour ecran
      write(*,*)'Post processing point cloud ...'
      
c     File
      Open(90,file='results/'// FILENAME //'.txt',form='formatted')
      
      Do NumPatch = 1,nb_patch
         
         If (ActivePatch(NumPatch) == 1) then
            
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         Do numel = 1,nb_elem_patch(NumPatch)
            
c     Extract element solution
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,numel))
            Enddo
            CALL extractNurbsElementInfos(numel)
            
c     Initialization
            XIb(:) = zero
            h = two/dble(MAX(2,nb_pts)-1)
            
            nb_xi  = MAX(1,nb_pts)
            nb_eta = 1
            nb_zeta= 1
            XIb(1) = -one
            if (dim_patch>1) then
               nb_eta = MAX(1,nb_pts)
               XIb(2) = -one
            endif
            if (dim_patch>2) then
               nb_zeta= MAX(1,nb_pts)
               XIb(3) = -one
            Endif
c     Calculation
            Do i_zeta= 1,nb_zeta
            Do i_eta = 1,nb_eta
            Do i_xi  = 1,nb_xi
               XI(:) = zero
               Do i = 1,dim_patch
                  XI(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*XIb(i)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
               End do
               
               CALL evalnurbs(XI(:),R(:nnode_patch),
     &              dRdxi(:nnode_patch,:))
               
               pts(:) = zero
               Do i = 1,nnode_patch
                  pts(:) = pts(:) + R(i)*COORDS_elem(:,i)
               Enddo
               write(90,*) pts(:)
               
            XIb(1) = XIb(1) + h
            Enddo
            XIb(2) = XIb(2) + h
            Enddo
            XIb(3) = XIb(3) + h
            Enddo
            
         Enddo                  ! loop elem
         
         call finalizeNurbsPatch()
         
      Endif                     ! active patch
      
      Enddo                     ! end loop on patch
      
      close(90)
      print*,' generated file: ',FILENAME,'.txt'
      
c     ..................................................................
            
      end subroutine generateCloudPtsTXT
