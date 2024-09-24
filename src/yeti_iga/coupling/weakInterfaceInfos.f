!! Copyright 2016-2019 Thibaut Hirschler

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
c     
c     Generation du tableau donnant des infos sur les interfaces 
c     associes a des multiplicateurs de lagrange
c     
      
      subroutine tabWEAKinfos(tabWEAK,nb_pts,nb_data,nb_interface,
     1     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,
     2     PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,MCRD,NBINT,nb_patch,
     3     nb_elem,nnode,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
      Implicit none
      
c     Declaration des Variables ........................................
      
c     Input arguments :
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
      Integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(:),
     &     NNODE(nb_patch),IEN(:),nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),NBINT(nb_patch)
            
!     Output size
      Integer, intent(in) :: nb_pts,nb_data,nb_interface
      
      
c     Output variables : system lineaire a resoudre
c     ----------------
      Double precision, intent(out) :: tabWEAK
      dimension tabWEAK(nb_data,nb_interface)
      
c     Local variables :
c     ---------------
      
!     Coupling infos
      Integer :: numDomain,numLgrge,ismaster,dimInterface,dispORrot
      
!     NURBS evalutation
      Double precision :: XI,saveXI,R,dRdxi,umin,umax,step,step2
      dimension XI(3),saveXI(3,nb_pts**2),R(MAXVAL(NNODE)),
     &     dRdxi(MAXVAL(NNODE),3),umin(2),umax(2)
      
!     Extract infos
      Integer          :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: COORDS_elem
      dimension COORDS_elem(3,MAXVAL(NNODE))
      
      Integer :: i,kl,k,l,count,numPatch,num_elem
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     Recuperation des infos ...........................................
      
      count = 0
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         If (ELT_TYPE_patch == 'U00') then
         
         count = count+1
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         numDomain = int(PROPS_patch(2))
         numLgrge  = int(PROPS_patch(3))
         ismaster  = int(PROPS_patch(4))
         dimInterface = dim_patch
         
         
         
c     1. Compute pts through the immerse curve/surface
         ! Bounds
         umin(1) = Ukv1_patch(1)
         umax(1) = Ukv1_patch(Nkv_patch(1))
         If (dimInterface==2) then
            umin(2) = Ukv2_patch(1)
            umax(2) = Ukv2_patch(Nkv_patch(2))
         Endif
         
         ! interface 1D
         IF (dimInterface==1) THEN
         XI(:) = zero
         XI(1) = umin(1)
         step  = (umax(1)-umin(1))/dble(nb_pts-1)
         Do k = 1,nb_pts
            call updateElementNumber(XI)
            call evalnurbs(XI,R(:nnode_patch),dRdxi(:nnode_patch,:))
            
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
            Enddo
            
            saveXI(:,k) = zero
            Do i = 1,nnode_patch
               saveXI(:,k) = saveXI(:,k) + R(i)*COORDS_elem(:,i)
            Enddo
            XI(1) = XI(1) + step
         Enddo

         
         ! interface 2D
         ELSEIF (dimInterface==2) THEN
         XI(:) = zero
         XI(2) = umin(2)
         step  = (umax(1)-umin(1))/dble(nb_pts-1)
         step2 = (umax(2)-umin(2))/dble(nb_pts-1)
         kl = 1
         Do k = 1,nb_pts
            XI(1) = umin(1)
            Do l = 1,nb_pts
               call updateElementNumber(XI)
               call evalnurbs(XI,R(:nnode_patch),dRdxi(:nnode_patch,:))
               
               sctr(:nnode_patch) = IEN_patch(:,current_elem)
               Do i = 1,nnode_patch
                  COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               Enddo
               
               saveXI(:,kl) = zero
               Do i = 1,nnode_patch
                  saveXI(:,kl) = saveXI(:,kl) + R(i)*COORDS_elem(:,i)
               Enddo
               
               kl = kl+1
               XI(1) = XI(1) + step
            Enddo
            XI(2) = XI(2) + step2
         Enddo
         ENDIF
         
         
c     2. Lagrangian field: type of constraint (disp/rot)
         CALL extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         dispORrot = int(PROPS_patch(2))
         
         
c     3. Domain to couple
         CALL extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         If (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            CALL extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif
         
         Do kl = 1,nb_pts**dimInterface
            call updateElementNumber(saveXI(:,kl))
            call evalnurbs(saveXI(:,kl),R(:nnode_patch),
     &           dRdxi(:nnode_patch,:))
            
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
            Enddo
            
            saveXI(:,kl) = zero
            Do i = 1,nnode_patch
               saveXI(:,kl) = saveXI(:,kl) + R(i)*COORDS_elem(:,i)
            Enddo
         Enddo
         
         ! Cas domain immerge - transformation supplementaire
         IF (ELT_TYPE_patch == 'U30') then
         Do kl = 1,nb_pts**dimInterface
            call updateMapElementNumber(saveXI(:,kl))
            call evalnurbs_mapping(saveXI(:,kl),R(:nnode_map),
     &           dRdxi(:nnode_map,:))
            
            sctr(:nnode_map) = IEN_map(:,current_map_elem)
            Do i = 1,nnode_map
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
            Enddo
            
            saveXI(:,kl) = zero
            Do i = 1,nnode_map
               saveXI(:,kl) = saveXI(:,kl) + R(i)*COORDS_elem(:,i)
            Enddo
         Enddo
         ENDIF
         
         
c     4. Merge infos into the output table
         Do kl = 1,nb_pts**dimInterface
            tabWEAK((kl-1)*3+1:kl*3,count) = saveXI(:,kl)
         Enddo
         tabWEAK(nb_data-3,count) = dble(ismaster )
         tabWEAK(nb_data-2,count) = dble(dispORrot)
         tabWEAK(nb_data-1,count) = dble(numLgrge )
         tabWEAK(nb_data  ,count) = dble(numDomain)
         
         Endif ! -- test if elt_type is 'U00'
         
         call deallocateMappingData()
         call finalizeNurbsPatch()
                  
      Enddo ! end loop on patch
      
c     Fin Assemblage ...................................................
      
      end subroutine tabWEAKinfos
