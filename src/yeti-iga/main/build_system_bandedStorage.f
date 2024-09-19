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

      
      
      
C     ******************************************************************
      
      
      subroutine sys_linmat_lindef_static(K, F,
     1     nb_diag,NUMDof2Diag,activeElement,
     2     COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     3     Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,
     4     indDLoad,JDLType,ADLMAG,load_target_nbelem,
     5     bc_values,nb_bc,bc_target,bc_target_nbelem,ind_dof_free,
     6     nb_dof_free,MCRD,NBINT,nb_load,
     7     nb_patch,nb_elem,nnode,nb_cp,nb_dof_tot)
      
      use parameters
      use nurbspatch
      
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
      dimension MATERIAL_PROPERTIES(2,nb_patch),
     &     PROPS(:),
     &     NNODE(nb_patch),
     &     IEN(:),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),
     &     NBINT(nb_patch)
      
      
!     Loads
      Double precision, intent(in) :: ADLMAG
      Integer, intent(in) :: nb_load,indDLoad,JDLType,load_target_nbelem
      dimension ADLMAG(nb_load),
     &     indDLoad(:),
     &     JDLType(nb_load),
     &     load_target_nbelem(nb_load)
      
      
!     Boundary Conditions
      Double precision, intent(in) :: bc_values
      Integer, intent(in) :: nb_bc,bc_target,bc_target_nbelem
      dimension bc_values(2,nb_bc),
     &     bc_target(:),
     &     bc_target_nbelem(nb_bc)
      
      
!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)
!      Logical, intent(in) :: COUPLG_flag, BNDSTRIP_flag
      
      
!     Banded storage INFOS
      Integer, intent(in) :: nb_diag,NUMDof2Diag, activeElement
      dimension NUMDof2Diag(nb_dof_tot),activeElement(nb_elem)
      
      
c     Output variables : system lineaire a resoudre
c     ----------------
      Double precision, intent(out) :: K, F
      dimension K(nb_dof_tot,nb_diag), F(nb_dof_tot)
            
      
c     Local variables :
c     ---------------
      
!     UELMAT.f
      Integer :: NDOFEL
      Double precision :: COORDS_elem, RHS, AMATRX, MAT_patch
      dimension COORDS_elem(MCRD,MAXVAL(NNODE)),RHS(MCRD*MAXVAL(NNODE)),
     &     AMATRX(MCRD*MAXVAL(NNODE),MCRD*MAXVAL(NNODE)),MAT_patch(2)
      
!     Global stiffness matrix and force vector
      Integer :: num_elem,i,j,JELEM,Numpatch, sctr,num_load,
     &     num_cp, ddl
      dimension sctr(MAXVAL(NNODE))
      
c     Test banded storage
      Integer :: n,m,dofi,dofj,ddli,ddlj,cpi,cpj, kk,ll
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
!     Initialize K to zero and F to concentrated loads
      K = zero
      F = zero
      kk = 0
      Do num_load = 1,nb_load
         i = JDLType(num_load)
         If (i/10 < 1) then
            Do num_cp = 1,load_target_nbelem(num_load)
               ddl = (indDLoad(kk+num_cp)-1)*MCRD + i
               F(ddl) = ADLMAG(num_load)
            Enddo
         Endif
         kk = kk + load_target_nbelem(num_load)
      Enddo
      
c     Fin Initialisation ...............................................
c     
c     
c     
c      
c     Debut Assemblage .................................................
      
      JELEM = 0
      Do NumPatch = 1,nb_patch
         
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         
         NDOFEL = nnode_patch*MCRD
c     Loop on element
         Do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1
            
            IF (activeElement(JELEM)==1) then

            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
            Enddo
            CALL extractNurbsElementInfos(num_elem)
            
            
c     Compute elementary matrix and load vector
            RHS    = zero
            AMATRX = zero
            MAT_patch(:) = MATERIAL_PROPERTIES(:,NumPatch)
            if (ELT_TYPE_patch == 'U1') then
               ! 'Element classique solide'
               call UELMAT(NDOFEL,MCRD,nnode_patch,JELEM,
     1              NBINT(NumPatch),COORDS_elem(:,:nnode_patch),
     2              TENSOR_patch,MAT_patch,nb_load,indDLoad,
     3              load_target_nbelem,JDLType,ADLMAG,RHS(:NDOFEL),
     4              AMATRX(:NDOFEL,:NDOFEL))
!               
!            elseif (ELT_TYPE(NumPatch) == 'U2') then
!               ! 'Element plaque'
!               call UELMAT2(NDOFEL,MCRD,nnode_patch,JELEM,NBINT,
!     1              NumPatch,COORDS_elem,TENSOR,MAT_patch,PROPS_patch,
!     2              JPROPS_patch,nb_dload,indDLoad,load_target_nbelem,
!     3              JDLType,ADLMAG,RHS,AMATRX(:NDOFEL,:NDOFEL))
!            
            elseif (ELT_TYPE_patch == 'U3') then
               ! 'Element coque'
               call UELMAT3(NDOFEL,MCRD,nnode_patch,JELEM,
     1              NBINT(NumPatch),COORDS_elem,TENSOR_patch,MAT_patch,
     2              PROPS_patch,JPROPS_patch,nb_load,indDLoad,
     3              load_target_nbelem,JDLType,ADLMAG,RHS(:NDOFEL),
     4              AMATRX(:NDOFEL,:NDOFEL))
            else
!               !print*, 'Element'//ELT_TYPE_patch//' not available.'
!               
            endif
            
            
c     Assemble AMATRX to global stiffness matrix K    
            sctr(:nnode_patch) = IEN_patch(:,num_elem)
            
            dofj = 0
            Do j = 1,nnode_patch
               cpj  = sctr(j)
               ddlj = (cpj-1)*MCRD
               
               Do kk = 1,MCRD
                  Do ll = kk,MCRD
                     K(ddlj+ll,kk) = K(ddlj+ll,kk) 
     &                    + AMATRX(dofj+ll-kk+1,dofj+ll)
                  Enddo
               Enddo
               
               
               dofi   = dofj + mcrd
               Do i   = j+1,nnode_patch
                  cpi = sctr(i)
                  
                  !n = abs(cpi - cpj)
                  n = cpj - cpi
                  m = NUMDof2Diag(abs(n)*MCRD + 1)
                  
                  if (n<0) then
                     n   = -n
                     ddl = (cpi-1)*MCRD
                     
                     Do kk = 2,MCRD
                        Do ll = kk,MCRD
                           K(ddl+ll,m+kk-1) 
     &                          = K(ddl+ll,m+kk-1)
     &                          + AMATRX(dofi+ll,dofj+1-kk+ll)
                           K(ddl+1-kk+ll, m-kk+1) 
     &                          = K(ddl+1-kk+ll,m-kk+1)
     &                          + AMATRX(dofi+1-kk+ll,dofj+ll)
                        Enddo
                     Enddo
                  else
                     ddl = ddlj
                     
                     Do kk = 2,MCRD
                        Do ll = kk,MCRD
                           K(ddl+ll,m+kk-1) 
     &                          = K(ddl+ll,m+kk-1)
     &                          + AMATRX(dofi+1-kk+ll,dofj+ll)
                           K(ddl+1-kk+ll, m-kk+1) 
     &                          = K(ddl+1-kk+ll,m-kk+1)
     &                          + AMATRX(dofi+ll,dofj+1-kk+ll)
                           
                        Enddo
                     Enddo
                  endif
                  
                  Do kk = 1,MCRD
                     K(ddl+kk,m) = K(ddl+kk,m)
     &                    + AMATRX(dofi+kk,dofj+kk)
                  Enddo
                  
                  dofi = dofi + MCRD
               Enddo
               
               dofj = dofj + MCRD
            Enddo
            
c     ..
c     Update Load Vector
            dofi = 0
            Do i = 1,nnode_patch
               ddli  = (sctr(i)-1)*MCRD
               Do kk = 1,MCRD
                  F(ddli+kk) = F(ddli+kk) + RHS(dofi+kk)
               Enddo
               dofi = dofi + MCRD
            Enddo
            
            Endif
            
         Enddo ! end loop on element
      Enddo ! end loop on patch

c     Fin Assemblage ...................................................
c     
c     Debut Modifications ..............................................
      
c     Build system to solve
!      call apply_dispBC(K,F,nb_cp,MCRD,nb_bc,bc_target,bc_target_nbelem,
!     &     bc_values,nb_dof_tot,U)
      
!      F_inv = zero
!      K_inv = zero
!
!      do i = 1,0 ! nb_dof_free
!         F_inv(i) = F(ind_dof_free(i))
!         do j = 1,nb_dof_free
!            K_inv(i,j) = K(ind_dof_free(i),ind_dof_free(j))
!         enddo
!      enddo
      
c     Fin Modifications ................................................
      
      end subroutine sys_linmat_lindef_static
