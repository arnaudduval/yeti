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

c     ...
c     Premiere etape pour l'analyse des structures immergees:
c      - transformation des points de gauss vers l'espace parametrique
c        du mapping a l'aide de la description de l'entite immergee
c      - calcul des vecteurs covariants
c      - table pour le lien numero point de gauss/numero de l'element
c        ou il se situe.
c     ...
      
      subroutine embeddedEntities2Param(XI,BI,dBI,GPsEltLOCATION,
     1     nb_points,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,
     3     NBINT,nb_patch,nb_elem,nb_cp,MCRD,
     4     ien_size,props_size,ukv_size,weight_size)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------

      Integer, intent(in) :: ien_size,props_size,ukv_size,weight_size
!     Geometry NURBS
      Integer,          intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(ukv_size),weight(weight_size)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer,          intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN, nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(props_size),
     &     JPROPS(nb_patch),NNODE(nb_patch),NBINT(nb_patch),
     &     IEN(ien_size),nb_elem_patch(nb_patch)
      
!     Output infos
      Integer,          intent(in) :: nb_points

!     Output arguments :
!     ----------------
      Double precision, intent(out):: XI,BI,dBI
      dimension XI(3,nb_points),BI(3,3,nb_points),dBI(3,6,nb_points)
      Integer,          intent(out):: GPsEltLOCATION
      dimension GPsEltLOCATION(nb_points)
      
      
C     Local variables :
c     ---------------
      
!     For gauss points
      Integer :: NbPtInt, n
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(4,MAXVAL(NBINT))

!     For nurbs basis functions
      Double precision :: COORDS_elem, R,dRdxi,ddRddxi, xi_bar
      dimension COORDS_elem(MCRD,MAXVAL(NNODE)), R(MAXVAL(NNODE)),
     &     dRdxi(MAXVAL(NNODE),3), ddRddxi(MAXVAL(NNODE),6), xi_bar(3)
      
!     For loops
      Integer :: numel,NumPatch, count,i,j, JPOINT
      
      
C     ------------------------------------------------------------------
      
C     Initialization :
c     --------------
      
      XI(:,:)    = zero
      BI(:,:,:)  = zero
      dBI(:,:,:) = zero
      
      
C     ------------------------------------------------------------------
      
C     Computing output data :
c     ---------------------
      
      count  = 0
      JPOINT = 0
      Do NumPatch = 1,nb_patch
         
c     Test if patch defines a stiffener or not ...
         n = count
         i = INDEX(ELT_TYPE(count:),    'U')-1
         j = INDEX(ELT_TYPE(count+i+1:),'U')-1
         if (j<0) then
            j=LEN(ELT_TYPE)-count-i
         else
            count = count+i+j
         end if
         If (ELT_TYPE(n+i:n+i+j) == 'U31') then
c     ... Test is true

            CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &           weight,nb_elem_patch)
            CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &           NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

            i = int(PROPS_patch(2))
            call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)


c     Get gauss infos
            NbPtInt = int( NBINT(numPatch)**(1.0/float(dim_patch)) )
            if (NbPtInt**dim_patch<NBINT(numPatch)) NbPtInt = NbPtInt+1
            call Gauss(NbPtInt,dim_patch,
     &           GaussPdsCoord(:dim_patch+1,:NBINT(numPatch)),0)
            
            
c     Loop on element
            Do numel = 1,nb_elem_patch(NumPatch)
               
c     Extract element infos
               Do i = 1,nnode_patch
                  COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,numel))
               Enddo
               CALL extractNurbsElementInfos(numel)
               
c     Loop on Gauss pts
               Do n = 1,NBINT(numPatch)
                  JPOINT = JPOINT + 1
                  
c     Map GP to embedded entities param. space
                  xi_bar(:) = zero
                  Do i = 1,dim_patch
                     xi_bar(i) = ((Ukv_elem(2,i) - Ukv_elem(1,i))
     &                    *GaussPdsCoord(1+i,n)
     &                    +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
                  Enddo
                  
c     Compute nurbs basis fcts and derivatives
                  call evalnurbs_w2ndDerv(xi_bar(:),R(:nnode_patch),
     &                 dRdxi(:nnode_patch,:),ddRddxi(:nnode_patch,:))
                  
c     Nurbs transformation
                  Do i = 1,nnode_patch
                     XI(:,JPOINT) = XI(:,JPOINT) + R(i)*COORDS_elem(:,i)
                  Enddo
c     Covariant vectors
                  Do j = 1,dim_patch
                     Do i = 1,nnode_patch
                        BI(:,j,JPOINT) = BI(:,j,JPOINT)
     &                       +   dRdxi(i,j)*COORDS_elem(:,i)
                     Enddo
                  Enddo
c     Derivative of covariant vectors
                  Do j = 1,dim_patch**2
                     Do i = 1,nnode_patch
                        dBI(:,j,JPOINT)= dBI(:,j,JPOINT)
     &                       + ddRddxi(i,j)*COORDS_elem(:,i)
                     Enddo
                  Enddo
                  
                  
c     Find location
                  call updateMapElementNumber(XI(:,JPOINT))
                  GPsEltLOCATION(JPOINT) = current_map_elem
                  
               Enddo
            Enddo
            
            CALL finalizeNurbsPatch()
            CALL deallocateMappingData()
            
         Endif
      Enddo
      
      
C     ------------------------------------------------------------------
      
      end subroutine embeddedEntities2Param
