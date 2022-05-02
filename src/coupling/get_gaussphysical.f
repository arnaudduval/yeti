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


C     CALCUL DE LA POSITION DANS L'ESPACE PHYSIQUE DES POINTS DE GAUSS
C     D'UNE PARAMETRISATION IGA
            
      
C     ******************************************************************
      
!      include 'evaluatenurbsfcts.f'
      
      Subroutine gaussPhysicalPosition(GaussCoords3D,COORDS3D,IEN,
     1     nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,PROPS,JPROPS,
     2     TENSOR,KNumFace,NBINT,nb_gauss_tot,MCRD,nb_patch,nb_elem,
     3     nb_cp,NNODE)
      
      use parameters
      use nurbspatch

      Implicit None
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer,          intent(in) :: MCRD,nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: PROPS
      Integer, intent(in) :: NNODE,nb_patch,nb_elem,nb_elem_patch,IEN,
     &     JPROPS
      dimension NNODE(nb_patch),nb_elem_patch(nb_patch),IEN(:),PROPS(:),
     &     JPROPS(nb_patch)
      
!     Gaussian points
      Integer,          intent(in) :: nb_gauss_tot,NBINT,KNumFace
      dimension NBINT(nb_patch)
      
c     Output variables : operateur de couplage
c     ----------------
      Double precision, intent(out) :: GaussCoords3D
      dimension GaussCoords3D(3,nb_gauss_tot)
      
      
      
      
c     Local Variables :
c     ---------------
      
      ! extract element infos
      Double precision :: COORDS_elem,R,dRdxi
      Integer          :: JELEM
      dimension COORDS_elem(3,MAXVAL(NNODE)),R(MAXVAL(NNODE)),
     &     dRdxi(MAXVAL(NNODE),3)
      
      ! gauss points
      Integer          :: num_gauss,NbPtInt,NBINT_face
      Double precision :: PtGauss,xi,GaussPdsCoord
      dimension PtGauss(3), xi(3), GaussPdsCoord(4,MAXVAL(NBINT))
            
      ! for local loops
      Integer          :: i,j,n,numPatch,num_elem,numCP
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Programme principal ..............................................
      
      GaussCoords3D(:,:) = zero

      num_gauss = 0
      JELEM     = 0
c     Loop on patches
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)


C     Gauss info
         NbPtInt = int(float(NBINT(NumPatch))**(one/float(dim_patch)))
         if (NbPtInt**dim_patch<NBINT(NumPatch)) NbPtInt = NbPtInt+1
         if (KNumFace==0) then
            j=0
         else
            j=1
         endif
         NBINT_face = NbPtInt**(dim_patch-j)

         GaussPdsCoord(:,:) = zero
         call gauss(NbPtInt,dim_patch,
     &        GaussPdsCoord(:dim_patch+1,:NBINT_face),KNumFace)
         
c     Loop on elements
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
c     Extract element infos of the geometry
            JELEM = JELEM + 1
            
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            
c     Loop on gauss points .............................................
            
            Do n = 1,NBINT_face
               num_gauss = num_gauss+1
               
!     Get gauss coordinates
               PtGauss(:) = GaussPdsCoord(2:,n)
               xi(:) = zero
               Do i  = 1,dim_patch
                  xi(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
               Enddo
               
!     Compute nurbs basis functions
               call evalnurbs(xi,R(:nnode_patch),dRdxi(:nnode_patch,:))
               
!     Get gauss physical coords
               Do numCP = 1,nnode_patch
                  GaussCoords3D(:,num_gauss)
     &                 = GaussCoords3D(:,num_gauss)
     &                 + R(numCP)*COORDS_elem(:,numCP)
               Enddo
               
            Enddo
!     End loop gauss points ............................................
         Enddo
!     End loop on Lagrange mult. elements

         CALL finalizeNurbsPatch()

      Enddo
!     End loop on Lagrange mult. patches
      
      
      End subroutine gaussPhysicalPosition
