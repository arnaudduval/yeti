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

      
      subroutine generateCPLGInterfaceTXT(FILENAME,FieldOutput_flag,
     1     nb_pts,SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,MCRD,
     3     NBINT,nnode,nb_patch,nb_elem,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
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
      
!     Analysis solution      
      Double precision, intent(in) :: SOL
      dimension SOL(3,nb_cp)
      
!     Output INFOS
      Character(len=*), intent(in) :: FILENAME
      Integer,          intent(in) :: nb_pts
      Logical,          intent(in) :: FieldOutput_flag
      dimension FieldOutput_flag(3)
      
      
c     Local variables :
c     ---------------
      
!     Coupling infos
      Integer :: numDomain,numLgrge,ismaster,dimInterface,nb_tot,
     &     dispORrot

!     Allocatable quantities
      Double precision, dimension(:,:),   allocatable :: saveXI,saveXIb
      Double precision, dimension(:,:,:), allocatable :: saveBI
      
!     Extract infos
      Integer          :: sctr,sctr_l
      dimension sctr(MAXVAL(NNODE)),sctr_l(MAXVAL(NNODE))
      Double precision :: COORDS_elem
      dimension COORDS_elem(MCRD,MAXVAL(NNODE))

!     Compute quantities
      Double precision :: R,dRdxi,U,AI,T,normA,invA,normT,coef1,coef2
      dimension R(MAXVAL(NNODE)), dRdxi(MAXVAL(NNODE),3),U(3),AI(3,3),
     &     T(3)
      
      Integer :: i,j,n,kk,count,numPatch,num_elem
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
c     Retour ecran
      write(*,*)'Post processsing coupling quantities ...'
      
c     File
      Open(90,file='results/'// FILENAME //'.txt',form='formatted')
            
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      count = 1
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         If ((ELT_TYPE_patch == 'U00') .AND. (dim_patch == 1)) then
         
         numDomain = int(PROPS_patch(2))
         numLgrge  = int(PROPS_patch(3))
         ismaster  = int(PROPS_patch(4))
         dimInterface = dim_patch
         
c     1. Integration points through the immerged patch
c        - save gps parametric coords to further evaluate the lagrangian field
c        - compute gps position on parameter patch
c        - compute covariant vectors
         
         nb_tot = nb_pts*nb_elem_patch(NumPatch)
         if (allocated(saveXIb)) deallocate(saveXIb)
         if (allocated(saveXI )) deallocate(saveXI )
         if (allocated(saveBI )) deallocate(saveBI )
         allocate(saveXIb(4,          nb_tot))
         allocate(saveXI( 3,          nb_tot))
         allocate(saveBI( 3,dim_patch,nb_tot))
         
         saveXIb(:,:)  = zero
         saveXI(:,:)   = zero
         saveBI(:,:,:) = zero
         
         kk = 0
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            
            Do n = 1,nb_pts
               kk = kk + 1
               saveXIb(1,kk) = 
     &              (Ukv_elem(2,1) - Ukv_elem(1,1))/dble(nb_pts-1)
     &              * dble(n-1) + Ukv_elem(1,1)
               
               
!     - evaluate basis functions
               call evalnurbs(saveXIb(:3,kk),R(:nnode_patch),
     &              dRdxi(:nnode_patch,:))

!     - get position
               Do i = 1,nnode_patch
                  saveXI(:MCRD,kk) = saveXI(:MCRD,kk)
     &                 + R(i)*COORDS_elem(:,i)
               Enddo
               
!     - get covariant basis vectors
               Do j = 1,dimInterface
                  Do i = 1,nnode_patch
                     saveBI(:MCRD,j,kk) = saveBI(:MCRD,j,kk) 
     &                    + dRdxi(i,j)*COORDS_elem(:,i)
                  Enddo
               Enddo
            Enddo
         Enddo
         

c     2. Lagrangian field
c        - compute basis functions
         
         CALL extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         dispORrot = int(PROPS_patch(2))
         
         IF (ismaster==1) then
         CALL extractNurbsPatchGeoInfos(numLgrge, Nkv,Jpqr,Nijk,Ukv,
     &           weight,nb_elem_patch)
         Do kk = 1,nb_tot

            CALL updateElementNumber(saveXIb(:3,kk))
            CALL evalLgrge(saveXIb(:3,kk),R(:nnode_patch))
            
            U(:) = zero
            Do i = 1,nnode_patch
               U(:) = U(:) + R(i)*SOL(:,IEN_patch(i,current_elem))
            Enddo
            
            write(90,*) saveXIb(1,kk),',',U(1),',',U(2),',',U(3)
            
         Enddo
         
         ENDIF
         
         
c     3. Domain to couple
c        - compute basis functions
c        - compute Jacobian for the integral
c        - build coupling matrix
         CALL extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         IF (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            CALL extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif
         
         
         Do kk = 1,nb_tot
            
            CALL updateElementNumber(saveXI(:,kk))
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,sctr(i))
            Enddo
            
            call evalnurbs(saveXI(:3,kk),R(:nnode_patch),
     &           dRdxi(:nnode_patch,:))
            
            
            IF (dispORrot==0) then
!           - champ de deplacement
               U(:) = zero
               Do i = 1,nnode_patch
                  U(:) = U(:) + R(i)*SOL(:,IEN_patch(i,current_elem))
               Enddo
               
               write(90,*) saveXIb(1,kk),',',U(1),',',U(2),',',U(3)
            
               
            Elseif (dispORrot==1) then
!           - rotation
               AI(:,:) = zero
               Do j = 1,dim_patch
                  Do i = 1,nnode_patch
                     AI(:MCRD,j) = AI(:MCRD,j) 
     &                    + dRdxi(i,j)*COORDS_elem(:,i)
                  Enddo
               Enddo
               call cross(AI(:,1),AI(:,2),AI(:,3))
               call norm(AI(:,3),3, normA)
               AI(:,3) = AI(:,3)/normA
               invA = one/normA

               T(:) = zero
               Do i = 1,dim_patch
                  T(:) = T(:) + saveBI(i,1,kk)*AI(:,i)
               Enddo
               call norm(T(:),3, normT)
               T(:) = T(:)/normT
               
               call dot(AI(:,1),T(:), coef1)
               call dot(AI(:,2),T(:), coef2)
               
               U(:) = zero
               Do i = 1,nnode_patch
                  U(1) = U(1) + invA
     &                 *SUM( AI(:,3)*SOL(:,IEN_patch(i,current_elem)) )
     &                 *( dRdxi(i,2)*coef1 - dRdxi(i,1)*coef2 )
               Enddo
               
               write(90,*) saveXIb(1,kk),',',U(1),',0.,0.'
               
            Endif
            
         Enddo
         Endif
         
         call deallocateMappingData()
         call finalizeNurbsPatch()
         
         if (allocated(saveXI ))  deallocate(saveXI )
         if (allocated(saveXIb))  deallocate(saveXIb)
         if (allocated(saveBI ))  deallocate(saveBI )
         
      Enddo ! end loop on patch
      
c     ..................................................................
      
c     Writing End of file
      
      close(90)
      
      write(*,*)' --> File ' // FILENAME // '.txt has been created.'
      
      end subroutine generateCPLGInterfaceTXT






































      subroutine displayAngleDeformConfig(FILENAME,FieldOutput_flag,
     1     nb_pts,SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,MCRD,
     3     NBINT,nnode,nb_patch,nb_elem,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
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
      
!     Analysis solution      
      Double precision, intent(in) :: SOL
      dimension SOL(3,nb_cp)
      
!     Output INFOS
      Character(len=*), intent(in) :: FILENAME
      Integer,          intent(in) :: nb_pts
      Logical,          intent(in) :: FieldOutput_flag
      dimension FieldOutput_flag(3)
      
      
c     Local variables :
c     ---------------
      
!     Coupling infos
      Integer :: numDomain,numLgrge,ismaster,dimInterface,nb_tot,
     &     dispORrot

!     Allocatable quantities
      Double precision, dimension(:,:),   allocatable :: saveXI,saveXIb
      Double precision, dimension(:,:,:), allocatable :: saveBI
      
!     Extract infos
      Integer          :: sctr,sctr_l
      dimension sctr(MAXVAL(NNODE)),sctr_l(MAXVAL(NNODE))
      Double precision :: COORDS_elem
      dimension COORDS_elem(MCRD,MAXVAL(NNODE))

!     Compute quantities
      Double precision :: R,dRdxi,U,AI,T,vN,normA,invA,normT,coef1,coef2
      dimension R(MAXVAL(NNODE)), dRdxi(MAXVAL(NNODE),3),U(3),AI(3,3),
     &     T(3),vN(3)
      
      Integer :: i,j,n,kk,count,numPatch,num_elem
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
c     Retour ecran
      write(*,*)'Post processsing angle after deformation ...'
      
c     File
      Open(90,file='results/'// FILENAME //'.txt',form='formatted')
            
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      count = 1
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         If ((ELT_TYPE_patch == 'U00') .AND. (dim_patch == 1)) then
         
         numDomain = int(PROPS_patch(2))
         numLgrge  = int(PROPS_patch(3))
         ismaster  = int(PROPS_patch(4))
         dimInterface = dim_patch
         
c     1. Integration points through the immerged patch
c        - save gps parametric coords to further evaluate the lagrangian field
c        - compute gps position on parameter patch
c        - compute covariant vectors
         
         nb_tot = nb_pts*nb_elem_patch(NumPatch)
         if (allocated(saveXIb)) deallocate(saveXIb)
         if (allocated(saveXI )) deallocate(saveXI )
         if (allocated(saveBI )) deallocate(saveBI )
         allocate(saveXIb(4,          nb_tot))
         allocate(saveXI( 3,          nb_tot))
         allocate(saveBI( 3,dim_patch,nb_tot))
         
         saveXIb(:,:)  = zero
         saveXI(:,:)   = zero
         saveBI(:,:,:) = zero
         
         kk = 0
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            
            Do n = 1,nb_pts
               kk = kk + 1
               saveXIb(1,kk) = 
     &              (Ukv_elem(2,1) - Ukv_elem(1,1))/dble(nb_pts-1)
     &              * dble(n-1) + Ukv_elem(1,1)
               
               
!     - evaluate basis functions
               call evalnurbs(saveXIb(:3,kk),R(:nnode_patch),
     &              dRdxi(:nnode_patch,:))

!     - get position
               Do i = 1,nnode_patch
                  saveXI(:MCRD,kk) = saveXI(:MCRD,kk)
     &                 + R(i)*COORDS_elem(:,i)
               Enddo
               
!     - get covariant basis vectors
               Do j = 1,dimInterface
                  Do i = 1,nnode_patch
                     saveBI(:MCRD,j,kk) = saveBI(:MCRD,j,kk) 
     &                    + dRdxi(i,j)*COORDS_elem(:,i)
                  Enddo
               Enddo
            Enddo
         Enddo
         

c     2. Lagrangian field
c        - compute basis functions
         
         CALL extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         dispORrot = int(PROPS_patch(2))
         
         IF (ismaster==1) then
         CALL extractNurbsPatchGeoInfos(numLgrge, Nkv,Jpqr,Nijk,Ukv,
     &           weight,nb_elem_patch)
         Do kk = 1,nb_tot

            CALL updateElementNumber(saveXIb(:3,kk))
            CALL evalLgrge(saveXIb(:3,kk),R(:nnode_patch))
            
            U(:) = zero
            Do i = 1,nnode_patch
               U(:) = U(:) + R(i)*SOL(:,IEN_patch(i,current_elem))
            Enddo
            
            !write(90,*) saveXIb(1,kk),',',U(1),',',U(2),',',U(3)
            
         Enddo
         
         ENDIF
         
         
c     3. Domain to couple
c        - compute basis functions
c        - compute Jacobian for the integral
c        - build coupling matrix
         CALL extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         IF (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            CALL extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif
         
         
         Do kk = 1,nb_tot
            
            CALL updateElementNumber(saveXI(:,kk))
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,sctr(i)) 
     &              + SOL(:MCRD,sctr(i))
            Enddo
            
            call evalnurbs(saveXI(:3,kk),R(:nnode_patch),
     &           dRdxi(:nnode_patch,:))
            
            
            If (dispORrot==1) then
!           - rotation
               AI(:,:) = zero
               Do j = 1,dim_patch
                  Do i = 1,nnode_patch
                     AI(:MCRD,j) = AI(:MCRD,j) 
     &                    + dRdxi(i,j)*COORDS_elem(:,i)
                  Enddo
               Enddo
               call cross(AI(:,1),AI(:,2),AI(:,3))
               call norm(AI(:,3),3, normA)
               AI(:,3) = AI(:,3)/normA
               invA = one/normA
               
               T(:) = zero
               Do i = 1,dim_patch
                  T(:) = T(:) + saveBI(i,1,kk)*AI(:,i)
               Enddo
               call norm(T(:),3, normT)
               T(:) = T(:)/normT
               
               call cross(AI(:,3),T(:),vN(:))
               call norm(vN(:),3, normT)
               vN(:) = vN(:)/normT
               
               !vN(:) = AI(:,3)

               write(90,*) saveXIb(1,kk),',',vN(1),',',vN(2),',',vN(3)
               
            Endif
            
         Enddo
         Endif
         
         call deallocateMappingData()
         call finalizeNurbsPatch()
         
         if (allocated(saveXI ))  deallocate(saveXI )
         if (allocated(saveXIb))  deallocate(saveXIb)
         if (allocated(saveBI ))  deallocate(saveBI )
         
      Enddo ! end loop on patch
      
c     ..................................................................
      
c     Writing End of file
      
      close(90)
      
      write(*,*)' --> File ' // FILENAME // '.txt has been created.'
      
      end subroutine displayAngleDeformConfig
