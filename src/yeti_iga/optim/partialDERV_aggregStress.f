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

c     Partial derivatives of the stress aggregation with respect to the
c     control points 
      
      Subroutine partialDervCPpnormStress(partialGradAS,
     1     Pnorm,PTSEVAL,nb_pts_patch,nb_pts,SOL,
     1     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,
     3     nb_patch,nb_elem,nb_cp,MCRD)
      
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
      Double precision, intent(out) :: partialGradAS
      dimension partialGradAS(6,3,nb_cp)
      
      
C     Local variables :
c     ---------------
      ! nurbs basis
      Integer :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: XI
      dimension XI(3)
      
      ! displacement field
      Double precision :: Uelem
      dimension Uelem(MCRD,MAXVAL(NNODE))
      
      ! curvilinear quantities
      Double precision :: COORDSelem
      dimension COORDSelem(3,MAXVAL(NNODE))
      
      ! material behavior
      Double precision :: MAT_patch
      dimension MAT_patch(2)

      ! strain/stress
      Integer          :: ntens
      Double precision :: stress,gradQstress,stressPN,gradElem,coef
      dimension stress(6),gradQstress(6,3,MAXVAL(NNODE)),stressPN(6),
     &     gradElem(6,3,MAXVAL(NNODE))
      
      ! other
      Integer :: i,k,cp,pts,NumPatch,num_elem,count
      
      
C     ------------------------------------------------------------------
      
      
      partialGradAS(:,:,:) = zero
      stressPN(:) = zero
      count = 0
      Do NumPatch = 1,nb_patch
         
         ! Extract infos
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         
         MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
         
         
         Do pts = 1,nb_pts_patch(numPatch)
            
            count = count + 1
            XI(:) = PTSEVAL(:,count)
            
            ! Extract current element infos
            CALL updateElementNumber(XI)
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            COORDSelem(:,:) = zero
            Uelem(:,:)      = zero
            Do cp = 1,nnode_patch
               COORDSelem(:,cp) = COORDS3D(:,sctr(cp))
               Uelem(:,cp)      = sol(:,sctr(cp))
            Enddo
            
            ! Compute derivative of stress w.r.t. the control points
            gradQstress(:,:,:)  = zero
            If (ELT_TYPE_patch == 'U1') then
               ! 'Element solide'
               call partialDervStress1(XI,COORDSelem(:,:nnode_patch),
     &              Uelem(:,:nnode_patch),TENSOR_patch,MAT_patch,
     &              PROPS_patch,JPROPS_patch,nnode_patch,MCRD,
     &              stress(:),gradQstress(:,:,:nnode_patch))
               ElseIf (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               call partialDervStress3(XI,COORDSelem(:,:nnode_patch),
     &              Uelem(:,:nnode_patch),TENSOR_patch,MAT_patch,
     &              PROPS_patch,JPROPS_patch,nnode_patch,
     &              stress(:),gradQstress(:,:3,:nnode_patch))
               
            Endif
            
            
            stressPN(:) = stressPN(:) + abs(stress(:))**dble(Pnorm)
            Do i = 1,6
               coef = one
               If (MOD(Pnorm,2) == 1) coef = SIGN(one,stress(i))
               gradElem(i,:,:nnode_patch) = 
     &              gradQstress(i,:,:nnode_patch)
     &              *coef*stress(i)**dble(Pnorm-1)
            Enddo
            
            Do cp = 1,nnode_patch
               k = sctr(cp)
               partialGradAS(:,:,k) = partialGradAS(:,:,k) 
     &              + gradElem(:,:,cp)
            Enddo
            
         Enddo ! end loop on element
         
         call finalizeNurbsPatch()
         
      Enddo ! end loop on patch
      
      Do i = 1,6
         partialGradAS(i,:,:) = 
     &        partialGradAS(i,:,:)*stressPN(i)**(one/dble(Pnorm)-one)
      Enddo
      
C     ------------------------------------------------------------------
      
      End Subroutine partialDervCPpnormStress




























c     Partial derivatives of the Von Mises stress aggregation with 
c     respect to the control points 
      
      Subroutine partialDervCPpnormVM(partialGradAVM,
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
      Double precision, intent(out) :: partialGradAVM
      dimension partialGradAVM(3,nb_cp)
      
      
C     Local variables :
c     ---------------
      ! nurbs basis
      Integer :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: XI
      dimension XI(3)
      
      ! displacement field
      Double precision :: Uelem
      dimension Uelem(MCRD,MAXVAL(NNODE))
      
      ! curvilinear quantities
      Double precision :: COORDSelem
      dimension COORDSelem(3,MAXVAL(NNODE))
      
      ! material behavior
      Double precision :: MAT_patch
      dimension MAT_patch(2)

      ! strain/stress
      Integer          :: ntens
      Double precision :: stress,stressVM,gradQstress,vonmisesPN,
     &     gradElem,coef,gradQvonmises
      dimension stress(6),gradQstress(6,3,MAXVAL(NNODE)),
     &     gradElem(3,MAXVAL(NNODE)),gradQvonmises(3,MAXVAL(NNODE))
      
      ! other
      Integer :: i,k,cp,pts,NumPatch,num_elem,count
      
      
C     ------------------------------------------------------------------
      
      
      partialGradAVM(:,:) = zero
      vonmisesPN  = zero
      count = 0
      Do NumPatch = 1,nb_patch
         
         ! Extract infos
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         
         MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
         
         
         Do pts = 1,nb_pts_patch(numPatch)
            
            count = count + 1
            XI(:) = PTSEVAL(:,count)
            
            ! Extract current element infos
            CALL updateElementNumber(XI)
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            COORDSelem(:,:) = zero
            Uelem(:,:)      = zero
            Do cp = 1,nnode_patch
               COORDSelem(:,cp) = COORDS3D(:,sctr(cp))
               Uelem(:,cp)      = sol(:,sctr(cp))
            Enddo
            
            ! Compute derivative of stress w.r.t. the control points
            gradElem(:,:)  = zero
            If (ELT_TYPE_patch == 'U1') then
               ! 'Element solide'
               call partialDervStress1(XI,COORDSelem(:,:nnode_patch),
     &              Uelem(:,:nnode_patch),TENSOR_patch,MAT_patch,
     &              PROPS_patch,JPROPS_patch,nnode_patch,MCRD,
     &              stress(:),gradQstress(:,:,:nnode_patch))
               ElseIf (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               call partialDervStress3(XI,COORDSelem(:,:nnode_patch),
     &              Uelem(:,:nnode_patch),TENSOR_patch,MAT_patch,
     &              PROPS_patch,JPROPS_patch,nnode_patch,
     &              stress(:),gradQstress(:,:3,:nnode_patch))
               
            Endif
            
            if (dim_patch==2) then ! plane stress
               if (MCRD==3) stress(4) = stress(3)
               stress(3)  = zero
               stress(5:) = zero
               
               if (MCRD==3) gradQstress(4,:,:) = gradQstress(3,:,:)
               gradQstress( 3,:,:) = zero
               gradQstress(5:,:,:) = zero
            endif
            stressVM = (stress(1)-stress(2))**2
     &           + (stress(1)-stress(3))**2
     &           + (stress(2)-stress(3))**2
     &           + 6.d0*(stress(4)**2 + stress(5)**2 + stress(6)**2)
            stressVM = SQRT(0.5d0*stressVM)
            
            vonmisesPN = vonmisesPN
     &           + abs(stressVM)**dble(Pnorm)
            
            gradQvonmises(:,:) = one/two/stressVM*(  
     &           (gradQstress(1,:,:)-gradQstress(2,:,:))
     &           *(stress(1)-stress(2))
     &           + (gradQstress(2,:,:)-gradQstress(3,:,:))
     &           *(stress(2)-stress(3))
     &           + (gradQstress(3,:,:)-gradQstress(1,:,:))
     &           *(stress(3)-stress(1))
     &           + 6.0d0*(gradQstress(4,:,:)*stress(4)
     &           +gradQstress(5,:,:)*stress(5)
     &           +gradQstress(6,:,:)*stress(6))
     &           )
            
            gradElem(:,:) =gradQvonmises(:,:)*stressVM**dble(Pnorm-1)
            
            Do cp = 1,nnode_patch
               k = sctr(cp)
               partialGradAVM(:,k)=partialGradAVM(:,k) + gradElem(:,cp)
            Enddo
            
         Enddo ! end loop on element
         
         call finalizeNurbsPatch()
         
      Enddo ! end loop on patch
      
      partialGradAVM(:,:) = 
     &     partialGradAVM(:,:)*vonmisesPN**(one/dble(Pnorm)-one)
      
C     ------------------------------------------------------------------
      
      End Subroutine partialDervCPpnormVM
