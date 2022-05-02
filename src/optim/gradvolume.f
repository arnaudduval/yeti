!! Copyright 2018-2020 Thibaut Hirschler

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
      
      
      Subroutine computeGradVolume(gradV, listpatch,
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
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
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
      Double precision, intent(out) :: gradV
      dimension gradV(3,nb_cp)
      
      
      
c     Local Variables :
c     ---------------
            
!     For Nurbs basis functions
      Double precision :: COORDS_elem,XI,R,dRdxi,detJ
      dimension COORDS_elem(3,MAXVAL(NNODE)),XI(3),R(MAXVAL(NNODE)),
     &     dRdxi(MAXVAL(NNODE),3)
      Integer :: i,k,l,kk,NumPatch,num_elem,sctr
      dimension sctr(MAXVAL(NNODE))
      
!     For embedded entities
      Double precision :: COORDSmap,Re,dRedxi,Rm,dRmdxi,ddRmddxi,BI,VI,
     &     dVI
      dimension COORDSmap(3,MAXVAL(NNODE)),BI(3,3),VI(3,3),
     &     Re(MAXVAL(NNODE)),dRedxi(MAXVAL(NNODE),3),
     &     Rm(MAXVAL(NNODE)),dRmdxi(MAXVAL(NNODE),3),
     &     ddRmddxi(MAXVAL(NNODE),6),dVI(3,6)
      Integer          :: sctr_map,isave
      dimension sctr_map(MAXVAL(NNODE))
      
!     For grad
      Double precision :: AIxAJ,AI,normV, dA1dPe,dA2dPe,dA1dPm,dA2dPm,
     &     coef1,coef2
      dimension AIxAJ(3,3),AI(3,3),  
     &     dA1dPe(3,3,MAXVAL(NNODE)),dA2dPe(3,3,MAXVAL(NNODE)),
     &     dA1dPm(3,3,MAXVAL(NNODE)),dA2dPm(3,3,MAXVAL(NNODE))
      
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
      gradV(:,:) = zero
c     Loop on patches
      Do NumPatch = 1,nb_patch
         
         IF (listpatch(numPatch) == 1) then
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         If (ELT_TYPE_patch == 'U30') then
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
            
            sctr(:nnode_patch) = IEN_patch(:,num_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
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
               
               R(:) = zero
               dRdxi(:,:) = zero
               call evalnurbs(XI(:),R(:nnode_patch),
     &              dRdxi(:nnode_patch,:))
               
c     Build tangent vectors

c     For embedded cases
               IF (ELT_TYPE_patch == 'U30') then
                  
                  Re(:) = R(:)
                  dRedxi(:,:) = dRdxi(:,:)
                  
                  XI(:)   = zero
                  BI(:,:) = zero
                  Do i = 1,nnode_patch
                     XI(:) = XI(:) + Re(i)*COORDS_elem(:,i)
                     BI(:,1) = BI(:,1) + dRedxi(i,1)*COORDS_elem(:,i)
                     BI(:,2) = BI(:,2) + dRedxi(i,2)*COORDS_elem(:,i)
                  Enddo
                  
                  call updateMapElementNumber(XI(:))
                  call evalnurbs_mapping_w2ndDerv(XI(:),Rm(:nnode_map),
     &                 dRmdxi(:nnode_map,:),ddRmddxi(:nnode_map,:))
                  
                  ! extract coords mapping
                  If (isave /= current_map_elem) then
                     sctr_map(:nnode_map) = IEN_map(:,current_map_elem)
                     
                     Do i = 1,nnode_map
                        COORDSmap(:,i) = COORDS3D(:,sctr_map(i))
                     Enddo
                     
                     isave = current_map_elem
                  Endif
                  
                  VI(:,:) = zero
                  Do k = 1,dim_map
                     Do i = 1,nnode_map
                        VI(:,k) = VI(:,k) + dRmdxi(i,k)*COORDSmap(:,i)
                     Enddo
                  Enddo
                  
                  dVi(:,:) = zero
                  Do k = 1,6
                     Do i = 1,nnode_map
                        dVI(:,k)=dVI(:,k) + ddRmddxi(i,k)*COORDSmap(:,i)
                     Enddo
                  Enddo
                  
                  AI(:,:) = zero
                  Do k = 1,dim_patch
                     Do i = 1,3
                        AI(:,k) = AI(:,k) + BI(i,k)*VI(:,i)
                     Enddo
                  Enddo

                  ! compute derivative versus the control points
                  ! - embedded CPs
                  dA1dPe(:,:,:) = zero
                  dA2dPe(:,:,:) = zero
                  Do i = 1,nnode_patch
                     Do l = 1,3
                        dA1dPe(:,l,i) = dRedxi(i,1)*VI(:,l) 
     &                       + BI(l,1)*R(i)*dVI(:,l)
                        dA2dPe(:,l,i) = dRedxi(i,2)*VI(:,l)
     &                       + BI(l,2)*R(i)*dVI(:,l)
                        Do k = 1,3
                        If (k /= l) then
                           kk = l+k+1
                           dA1dPe(:,l,i) = dA1dPe(:,l,i)
     &                          + BI(k,1)*R(i)*dVI(:,kk)
                           dA2dPe(:,l,i) = dA2dPe(:,l,i)
     &                          + BI(k,2)*R(i)*dVI(:,kk)
                        Endif
                        Enddo
                     Enddo
                  Enddo
                                    
                  ! - mapping CPs
                  dA1dPm(:,:,:) = zero
                  dA2dPm(:,:,:) = zero

                  Do i = 1,nnode_map
                     coef1 = SUM( BI(:,1)*dRmdxi(i,:) )
                     coef2 = SUM( BI(:,2)*dRmdxi(i,:) )
                     Do k = 1,3
                        dA1dPm(k,k,i) = coef1
                        dA2dPm(k,k,i) = coef2
                     Enddo
                  Enddo

                  
               ELSE
                  
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
                  ! --> to do
               Elseif (dim_patch == 2) then
                  ! surface
                  call cross(AI(:,1),AI(:,2), AI(:,3))
                  call norm(AI(:,3),3, normV)
                  AI(:,3) = AI(:,3)/normV
                  call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
                  call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
                  
                  If (ELT_TYPE_patch == 'U30') then

                  ! embedded surface
                  Do i = 1,nnode_patch
                     k = sctr(i)
                     Do l = 1,3
                        gradV(l,k) = gradV(l,k) + 
     &                       ( SUM( AIxAJ(:,1) * dA1dPe(:,l,i) )
     &                       + SUM( AIxAJ(:,2) * dA2dPe(:,l,i) )
     &                       ) * detJ
                     Enddo
                  Enddo
                  
                  ! mapping
                  Do i = 1,nnode_map
                     k = sctr_map(i)
                     Do l = 1,3
                        gradV(l,k) = gradV(l,k) + 
     &                       ( SUM( AIxAJ(:,1) * dA1dPm(:,l,i) )
     &                       + SUM( AIxAJ(:,2) * dA2dPm(:,l,i) )
     &                       ) * detJ
                     Enddo
                  Enddo
                  
                  Else
                  
                  Do i = 1,nnode_patch
                     k = sctr(i)
                     gradV(:,k) = gradV(:,k) + 
     &                    ( AIxAJ(:,1)*dRdxi(i,1)
     &                    + AIxAJ(:,2)*dRdxi(i,2)
     &                    ) * detJ
                  Enddo
                     
                  Endif


               Elseif (dim_patch == 3) then
                  ! volume
                  call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
                  call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
                  call cross(AI(:,3),AI(:,1), AIxAJ(:,2))

                  Do i = 1,nnode_patch
                     k = sctr(i)
                     gradV(:,k) = gradV(:,k) + 
     &                    ( AIxAJ(:,3)*dRdxi(i,3)
     &                    + AIxAJ(:,1)*dRdxi(i,1)
     &                    + AIxAJ(:,2)*dRdxi(i,2)
     &                    ) * detJ
                  Enddo
               Endif
               
            Enddo
         Enddo
         
         Endif

         CALL finalizeNurbsPatch()
         
      Enddo
      
      
c     Fin Assemblage ...................................................
      
      End subroutine computeGradVolume
