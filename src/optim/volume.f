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

C     Calcul du volume
      
C     ******************************************************************
      
      
      Subroutine computeVolume(VOL, listpatch,
     1     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,
     2     PROPS,JPROPS,TENSOR,MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer, intent(in)          :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in)          :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,
     &     nb_elem_patch,IEN,JPROPS
      dimension PROPS(:), JPROPS(nb_patch), nb_elem_patch(nb_patch),
     &     IEN(:), NNODE(nb_patch), NBINT(nb_patch)
      
!     Other infos
      Integer, intent(in) :: listpatch
      dimension listpatch(nb_patch)
      
c     Output variables : coefficient diag matrice de masse
c     ----------------
      Double precision, intent(out) :: VOL
      
      
      
      
c     Local Variables :
c     ---------------
            
!     For Nurbs basis functions
      Double precision :: COORDS_elem,XI,R,dRdxi,detJ,normV,AI,vectV
      dimension COORDS_elem(3,MAXVAL(NNODE)),XI(3),R(MAXVAL(NNODE)),
     &     dRdxi(MAXVAL(NNODE),3),AI(3,3),vectV(3)
      Integer :: i,k,NumPatch,num_elem
      
!     For embedded entities
      Double precision :: COORDSmap,Re,dRedxi,BI
      dimension COORDSmap(3,MAXVAL(NNODE)),Re(MAXVAL(NNODE)),
     &     dRedxi(MAXVAL(NNODE),3),BI(3,3)
      Integer          :: sctr_map,isave
      dimension sctr_map(MAXVAL(NNODE))

!     For gauss points
      Integer :: NbPtInt, n
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(4,MAXVAL(NBINT))
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Calcul Volume ....................................................
      
      isave = 0
      VOL = zero
c     Loop on patches
      Do NumPatch = 1,nb_patch
         
         IF (listpatch(numPatch) == 1) then

         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         if (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif

c     Get gauss infos
         NbPtInt = int( NBINT(numPatch)**(1.0/float(dim_patch)) )
         if (NbPtInt**dim_patch<NBINT(numPatch)) NbPtInt = NbPtInt+1
         call Gauss(NbPtInt,dim_patch,
     &        GaussPdsCoord(:dim_patch+1,:NBINT(numPatch)),0)
         
         
c     Loop on elements
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
            COORDS_elem(:,:) = zero
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            
c     Loop on gauss points
            Do n = 1,NBINT(numPatch)
               
c     Evaluate nurbs basis functions
               XI(:) = zero
               DetJ  = GaussPdsCoord(1,n)
               Do i = 1,dim_patch
                  XI(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))
     &                 *GaussPdsCoord(1+i,n)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
                  DetJ = DetJ * 0.5d0*(Ukv_elem(2,i) - Ukv_elem(1,i))
               Enddo
               
               call evalnurbs(XI(:),R(:nnode_patch),
     &              dRdxi(:nnode_patch,:))
               
c     For embedded cases
               IF (ELT_TYPE_patch == 'U30') then
                  
                  XI(:)   = zero
                  BI(:,:) = zero
                  Do i = 1,nnode_patch
                     XI(:) = XI(:) + R(i)*COORDS_elem(:,i)
                     BI(:,1) = BI(:,1) + dRdxi(i,1)*COORDS_elem(:,i)
                     BI(:,2) = BI(:,2) + dRdxi(i,2)*COORDS_elem(:,i)
                  Enddo
                  
                  dRedxi(:,:) = zero
                  call updateMapElementNumber(XI(:))
                  call evalnurbs_mapping(XI(:),Re(:nnode_map),
     &                 dRedxi(:nnode_map,:))
                  
                  dRdxi(:,:) = zero
                  Do i = 1,3
                     dRdxi(:,1) = dRdxi(:,1) + BI(i,1)*dRedxi(:,i)
                     dRdxi(:,2) = dRdxi(:,2) + BI(i,2)*dRedxi(:,i)
                  Enddo

!                 extract COORDS
                  If (isave /= current_map_elem) then
                     sctr_map(:nnode_map) = IEN_map(:,current_map_elem)
                     
                     Do i = 1,nnode_map
                        COORDSmap(:,i) = COORDS3D(:,sctr_map(i))
                     Enddo
                     
                     isave = current_map_elem

                  Endif
                  
                  AI(:,:) = zero
                  Do k = 1,dim_patch
                     Do i = 1,nnode_map
                        AI(:,k) = AI(:,k) + dRdxi(i,k)*COORDSmap(:,i)
                     Enddo
                  Enddo
                  
               Else
c     Classical patch
                  
                  AI(:,:) = zero
                  Do k = 1,dim_patch
                     Do i = 1,nnode_patch
                        AI(:,k) = AI(:,k) + dRdxi(i,k)*COORDS_elem(:,i)
                     Enddo
                  Enddo

               Endif


               
c     Update volume
               If     (dim_patch == 1) then
                  ! courbe
                  call norm(AI(:,1),3, normV)
               Elseif (dim_patch == 2) then
                  ! surface
                  call cross(AI(:,1),AI(:,2),vectV(:))
                  call norm(vectV(:),3, normV)
               Elseif (dim_patch == 3) then
                  ! volume
                  call cross(AI(:,1),AI(:,2),vectV(:))
                  call dot(  AI(:,3),vectV(:),normV)
               Endif
               
               VOL = VOL + normV*detJ
               
            Enddo
         Enddo
         
         Endif
         
         CALL finalizeNurbsPatch()

      Enddo
      
      
c     Fin Assemblage ...................................................
      
      End subroutine computeVolume
