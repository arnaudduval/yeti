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

C     Gradient analytique de l'operateur de couplage
      
C     ******************************************************************
            
      Subroutine gradCplg_AN(gradC, listpatch,
     1     SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,PROPS,JPROPS,TENSOR,MCRD,NBINT,nb_patch,nb_elem,
     3     nnode,nb_cp)
      
      use parameters
      use nurbspatch
      use embeddedMapping
      
      Implicit None
            
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: SOL,COORDS3D
      dimension SOL(3,nb_cp),COORDS3D(3,nb_cp)
      
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
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: gradC
      dimension gradC(3,nb_cp)
      
      
      
c     Local Variables :
c     ---------------
      
!     Coupling infos
      Integer :: numDomain,numLgrge,ismaster,dimInterface,nb_gps,
     &     dispORrot

!     Interface
      Integer :: nnode_interface
      Double precision,dimension(:,:),  allocatable :: saveXI,saveXIb,
     &     saveRI
      Double precision,dimension(:,:,:),allocatable :: saveBI,save_dRdxi
      Integer         ,dimension(:),    allocatable :: saveEL
      Integer         ,dimension(:,:),  allocatable :: saveIEN
      
!     Lagrange mult
      Double precision,dimension(:,:),  allocatable :: lmbda
      Double precision :: R
      dimension R(MAXVAL(NNODE))
      
!     Domain patch
      Double precision :: U_elem,COORDS_elem
      dimension U_elem(3,MAXVAL(NNODE)),COORDS_elem(3,MAXVAL(NNODE))
      
!     Common variables
      Integer          :: i,k,n,kk,NumPatch,num_elem,sctr
      dimension sctr(MAXVAL(NNODE))
      
!     For grad
      Double precision :: coef,dCdP,dCdPI,dCdPm
      dimension dCdP(3,MAXVAL(NNODE)),dCdPI(3,MAXVAL(NNODE)),
     &     dCdPm(3,MAXVAL(NNODE))
      
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Calcul Volume ....................................................
            
      gradC(:,:) = zero
c     Loop on patches
      Do NumPatch = 1,nb_patch

         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         If (ELT_TYPE_patch == 'U00') then
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         numDomain = int(PROPS_patch(2))
         numLgrge  = int(PROPS_patch(3))
         ismaster  = int(PROPS_patch(4))
         dimInterface = dim_patch

         If ((listpatch(numPatch)==1).OR.(listpatch(numDomain)==1)) then
         
c     1. Interface curve/surface
c        - save gps parametric coords to further evaluate the lagrangian field
c        - compute gps position on parameter patch
c        - compute covariant vectors
c        - compute derivative of covariant vectors
         
         n = NBINT(numPatch)
         nb_gps = n*nb_elem_patch(NumPatch)
         nnode_interface = nnode_patch
         if (allocated(saveXIb)) deallocate(saveXIb)
         if (allocated(saveXI )) deallocate(saveXI )
         if (allocated(saveBI )) deallocate(saveBI )
         if (allocated(saveRI )) deallocate(saveRI )         
         if (allocated(save_dRdxi))  deallocate(save_dRdxi)
         if (allocated(saveEL))  deallocate(saveEL)
         if (allocated(saveIEN)) deallocate(saveIEN)
         allocate(saveXIb(4,          nb_gps))
         allocate(saveXI( 3,          nb_gps))
         allocate(saveBI( 3,dim_patch,nb_gps))
         allocate(saveRI(nnode_patch, nb_gps))
         allocate(save_dRdxi(nnode_patch,dim_patch,nb_gps))
         allocate(saveEL(nb_gps))
         allocate(saveIEN(nnode_patch,nb_elem_patch(NumPatch)))
         
         saveIEN(:,:) = IEN_patch(:,:)
         
         kk = 0
         Do num_elem = 1,nb_elem_patch(NumPatch)
            
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            saveEL(kk+1:kk+n) = current_elem
            
            CALL getGPsOnParamSpace_wBasis(
     &           saveXIb( :,kk+1:kk+n),saveXI(:,kk+1:kk+n),
     &           saveBI(:,:,kk+1:kk+n),saveRI(:,kk+1:kk+n), 
     &           save_dRdxi(:,:,kk+1:kk+n),
     &           dim_patch,MCRD,nnode_patch,NBINT(numPatch),
     &           COORDS_elem(:,:nnode_patch)                     )
                        
            kk = kk + n
            
         Enddo
         
         
c     2. Lagrangian field
c        - evaluate lagrangian field at GPs
         CALL extractNurbsPatchGeoInfos(numLgrge, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(numLgrge,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         if (allocated(lmbda)) deallocate(lmbda)
         allocate(lmbda(3,nb_gps))
         
         dispORrot = int(PROPS_patch(2))
         lmbda(:,:)= zero
         Do kk = 1,nb_gps
            
            CALL updateElementNumber(saveXIb(:3,kk))
            CALL evalLgrge(saveXIb(:,kk),R(:nnode_patch))
            
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               U_elem(:,i) = SOL(:,sctr(i))
               lmbda(:,kk) = lmbda(:,kk) + R(i)*U_elem(:MCRD,i)
            Enddo
            
         Enddo
         
         
c     3. Subdomain
c        - evaluate disp field at GPs
c        - derivative of tangent vector to the interface
         CALL extractNurbsPatchGeoInfos(numDomain, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(numDomain,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         IF (ELT_TYPE_patch == 'U30') then
            i = int(PROPS_patch(2))
            CALL extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)         
         Endif
         
         
         Do kk = 1,nb_gps
            
            CALL updateElementNumber(saveXI(:,kk))
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
               U_elem(:,i)      = SOL(:,sctr(i))
            Enddo
            
            
            If (ELT_TYPE_patch=='U30') then
               IF  (dispORrot==0) then
                  call gradCPLGdispEMBD(saveXI(:,kk),saveRI(:,kk),
     &                 save_dRdxi(:,:,kk),saveBI(:,:,kk),dim_map,
     &                 dim_patch,dimInterface,MCRD,nnode_map,
     &                 nnode_patch,nnode_interface,nb_cp,
     &                 COORDS_elem(:,:nnode_patch),COORDS3D(:,:),
     &                 U_elem(:,:nnode_patch),lmbda(:,kk),
     &                 dCdPI(:,:nnode_interface),dCdP(:,:nnode_patch),
     &                 dCdPm(:,:nnode_map) )
               Else
                  call gradCPLGrotEMBD(saveXI(:,kk),saveRI(:,kk),
     &                 save_dRdxi(:,:,kk),saveBI(:,:,kk),dim_map,
     &                 dim_patch,dimInterface,MCRD,nnode_map,
     &                 nnode_patch,nnode_interface,nb_cp,
     &                 COORDS_elem(:,:nnode_patch),COORDS3D(:,:),
     &                 U_elem(:,:nnode_patch),lmbda(:,kk),
     &                 dCdPI(:,:nnode_interface),dCdP(:,:nnode_patch),
     &                 dCdPm(:,:nnode_map) )
               Endif
            Else
               IF     (dispORrot==0) then
                  call gradCPLGdisp(saveXI(:,kk),saveRI(:,kk),
     &                 save_dRdxi(:,:,kk),saveBI(:,:,kk),dim_patch,
     &                 dimInterface,MCRD,nnode_patch,nnode_interface,
     &                 COORDS_elem(:,:nnode_patch),
     &                 U_elem(:,:nnode_patch),lmbda(:,kk),
     &                 dCdPI(:,:nnode_interface),dCdP(:,:nnode_patch) )
               Else
                  call gradCPLGrot(saveXI(:,kk),saveRI(:,kk),
     &                 save_dRdxi(:,:,kk),saveBI(:,:,kk),dim_patch,
     &                 dimInterface,MCRD,nnode_patch,nnode_interface,
     &                 COORDS_elem(:,:nnode_patch),
     &                 U_elem(:,:nnode_patch),lmbda(:,kk),
     &                 dCdPI(:,:nnode_interface),dCdP(:,:nnode_patch) )
               Endif
            Endif
            
            
c     Assembly
!     Update gradient
            ! - interface CPs
            coef = saveXIb(4,kk)
            If (ismaster == 0) coef = -coef
            Do i = 1,nnode_interface
               k = saveIEN(i,saveEL(kk))
               gradC(:,k) = gradC(:,k) + coef*dCdPI(:,i)
            Enddo
            
            ! - patch CPs
            Do i = 1,nnode_patch
               k = sctr(i)
               gradC(:,k) = gradC(:,k) + coef*dCdP(:,i)
            Enddo
            
            ! - mapping CPs (for embedded entities)
            If (ELT_TYPE_patch=='U30') then
            sctr(:nnode_map) = IEN_map(:,current_map_elem)
            Do i = 1,nnode_map
               k = sctr(i)
               gradC(:,k) = gradC(:,k) + coef*dCdPm(:,i)
            Enddo
            Endif

         Enddo
            
         ENDIF ! test if active patches
         ENDIF ! interface
         
         call deallocateMappingData()
         CALL finalizeNurbsPatch()
         
         if (allocated(saveXIb)) deallocate(saveXIb)
         if (allocated(saveXI )) deallocate(saveXI )
         if (allocated(saveBI )) deallocate(saveBI )
         if (allocated(saveRI )) deallocate(saveRI )
         if (allocated(save_dRdxi))  deallocate(save_dRdxi)
         if (allocated(saveEL))  deallocate(saveEL)
         if (allocated(saveIEN)) deallocate(saveIEN)
         if (allocated(lmbda))   deallocate(lmbda)
         
      Enddo
            
c     Fin Assemblage ...................................................
      
      End subroutine gradCplg_AN
