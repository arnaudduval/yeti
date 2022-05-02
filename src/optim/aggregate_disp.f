!! Copyright 2020 Thibaut Hirschler

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

      
      subroutine pnormDispMagn(dispPN,
     1     Pnorm,PTSEVAL,nb_pts_patch,nb_pts,SOL,
     2     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     3     ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,
     4     nb_patch,nb_elem,nb_cp,MCRD)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(:),
     &     JPROPS(nb_patch),NNODE(nb_patch),IEN(:),
     &     nb_elem_patch(nb_patch)
      
!     Analysis solution      
      Double precision, intent(in) :: SOL
      dimension SOL(MCRD,nb_cp)
            
!     Type of norm
      Integer, intent(in) :: Pnorm

!     Discrete pts where to evaluate the VM stress
      Integer,          intent(in) :: nb_pts_patch,nb_pts
      Double precision, intent(in) :: PTSEVAL
      dimension PTSEVAL(3,nb_pts),nb_pts_patch(nb_patch)
      

C     Output variables :
c     ----------------
      Double precision, intent(out) :: dispPN
      
      
C     Local variables :
c     ---------------
      ! nurbs basis
      Integer :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: XI,R,dRdxi
      dimension XI(3), R(MAXVAL(NNODE)), dRdxi(MAXVAL(NNODE),3)
      
      ! displacement field
      Double precision :: Uelem,Disp,DispNorm
      dimension Uelem(MCRD,MAXVAL(NNODE)),Disp(3)
      
      ! other
      Integer :: i,cp,NumPatch,count,pts
      
      
C     ------------------------------------------------------------------
      
      
      
      dispPN  = zero
      count = 0
      Do NumPatch = 1,nb_patch
         
         ! Extract patch infos
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         
         Do pts = 1,nb_pts_patch(numPatch)
            
            count = count + 1
            XI(:) = PTSEVAL(:,count)
            
            ! Extract current element infos
            CALL updateElementNumber(XI)
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Uelem(:,:)      = zero
            Do cp = 1,nnode_patch
               Uelem(:,cp)      = sol(:,sctr(cp))
            Enddo

            ! Compute disp
            call evalnurbs(Xi,R,dRdxi)
            Disp(:) = zero
            Do cp = 1,nnode_patch
               Disp(:MCRD) = Disp(:MCRD) + R(cp)*Uelem(:,cp)
            Enddo
            DispNorm = (SUM(Disp(:)**two))**0.5d0
            
            dispPN = dispPN + abs(dispNorm)**dble(Pnorm)
            

         Enddo
               
         call finalizeNurbsPatch()
         
      Enddo ! end loop on patch
      
      
      dispPN = dispPN**(one/dble(Pnorm))

C     ------------------------------------------------------------------
      
      end subroutine pnormDispMagn
