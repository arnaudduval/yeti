!! Copyright 2016-2018 Thibaut Hirschler

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

c     Include iga subroutines
      include "./shap.f"
      include "./Gauss.f"
      include "./operateurs.f"
      include "./dersbasisfuns.f"
      include "./stiffmatrix.f"
      include "./UELMAT.f"
      include "./material_lib.f"
      include "./applyBC.f"
      include "./applyDispBC.f"
      include "./resolution.f"
      include "./shapPress.f" 
      include "./strong_coupling.f"
      include "./indDOF.f"
      include "./reconstruction.f"
      
      include "./plate/shap_4KL.f"
      include "./plate/dersbasisfuns_4KL.f"
      include "./plate/UELMAT_4KL.f"
      include "./plate/shapPress_4KL.f" 
      include "./plate/USFMEM.f"
      include "./plate/USFBND.f"
      
      include "./shell/curvilinearCoordinates.f"
      include "./shell/nurbsbasisfuns.f"
      include "./shell/UELMAT_shell.f"
      include "./shell/USFBND_shell.f" 
      include "./shell/USFMEM_shell.f"
      include "./shell/bendingstrip.f"
      include "./shell/bsplinebasisfuns.f"
      include "./shell/UELMAT_bndstrip_v2.f"
      
      

c     Include lapack subroutines
      include "../ext/lapack/lapack_routine/ieeeck.f"
      include "../ext/lapack/lapack_routine/ilaenv.f"
      include "../ext/lapack/lapack_routine/iparmq.f"
      include "../ext/lapack/lapack_routine/lsame.f"
      include "../ext/lapack/lapack_routine/xerbla.f"
      include "../ext/lapack/double/disnan.f"
      include "../ext/lapack/double/dlaisnan.f"
      include "../ext/lapack/double/dposv.f"
      include "../ext/lapack/double/dpotf2.f"
      include "../ext/lapack/double/dpotrf.f"
      include "../ext/lapack/double/dpotrs.f"  



C     ******************************************************************


      subroutine iga_linmat_lindef_static(COORDS3D,IEN,nb_elem_patch ,
     1     Nkv_e,Ukv_e,Nijk_e,weight_e,Jpqr_e,ELT_TYPE,PROPS,JPROPS  ,
     2     MATERIAL_PROPERTIES,TENSOR, bc_target,sol,indDLoad,JDLType,
     3     ADLMAG,bc_target_nbelem,load_target_nbelem,bc_values,MCRD ,
     4     NBINT,nb_bc,nb_load,nb_cload,nb_patch,nb_elem,nnode,nb_cp )
      
       
      Implicit none
      
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------

!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(nb_cp,3)
      
      Integer Lknot, Lpatch, Lelement, Lnode
      Parameter (Lelement=10000, Lnode=10000, Lknot=10000, Lpatch=100)
      Double precision, intent(in) :: Ukv_e, weight_e
      Integer, intent(in) :: Nkv_e, Jpqr_e, Nijk_e
      dimension Ukv_e(Lknot,3,Lpatch),
     &     Nkv_e(3,Lpatch),
     &     weight_e(Lelement,Lnode),
     &     Jpqr_e(3),
     &     Nijk_e(Lelement,3)
      
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(nb_patch,2),
     &     PROPS(nb_patch,10),
     &     IEN(nb_elem,NNODE),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch)

      
      
!     Loads
      Double precision, intent(in) :: ADLMAG
      Integer, intent(in) :: nb_load,nb_cload,indDLoad,JDLType,
     &     load_target_nbelem
      dimension ADLMAG(1000),
     &     indDLoad(1000,1000),
     &     JDLType(1000),
     &     load_target_nbelem(1000)
      
      
!     Boundary Conditions
      Double precision, intent(in) :: bc_values
      Integer, intent(in) :: nb_bc,bc_target,bc_target_nbelem
      dimension bc_values(1000,2),
     &     bc_target(1000,1000),
     &     bc_target_nbelem(1000)
      
      
      
      
c     Output variables : champ des deplacements
c     ----------------
      Double precision, intent(out) :: sol
      dimension sol(nb_cp,MCRD)
      
      
      
      
c     Local variables :
c     ---------------
      
!     Parameters and COMMON variables
      Double precision, parameter :: zero=0.0D0, one=1.0D0, two=2.0D0

      Common /NurbsParameter/ Ukv,weight, Nkv, Jpqr,Nijk      
      Double precision :: Ukv, weight
      Integer :: Nkv, Jpqr, Nijk
      dimension Ukv(Lknot,3,Lpatch), Nkv(3,Lpatch),
     &     weight(Lelement,Lnode), Jpqr(3), Nijk(Lelement,3)
      
!     indDOF.f and applyBC.f
      Logical :: COUPLG_flag, BNDSTRIP_flag
      Integer nb_dof_free, nb_dof_bloq, nb_dof_tot
      Integer, dimension(nb_cp*MCRD) :: ind_dof_bloq, ind_dof_free
      
!     UELMAT.f
      Integer :: NDOFEL, JPROPS_patch, nb_dload
      Double precision :: COORDS_elem, RHS, AMATRX, dvol_save,MAT_patch,
     &     PROPS_patch, MAT_bndstrip,PROPS_bndstrip
      dimension COORDS_elem(MCRD, NNODE), RHS(MCRD*NNODE),
     &     AMATRX(MCRD*NNODE,MCRD*NNODE), dvol_save(NBINT),MAT_patch(2),
     &     PROPS_patch(10), MAT_bndstrip(2),PROPS_bndstrip(10)
      
!     Global stiffness matrix and force vector
      Double precision, dimension(:,:), allocatable :: K_inv
      Double precision, dimension(:), allocatable :: F_inv, U_inv
      Double precision :: K, F, U, u_elem
      dimension K(nb_cp*MCRD,nb_cp*MCRD), F(nb_cp*MCRD), U(nb_cp*MCRD),
     &     u_elem(NNODE*MCRD)
      Integer :: num_elem,i,j,k1,k2,JELEM,Numpatch, sctr,sctrB,num_load,
     &     num_cp, ddl, kk
      dimension sctr(NNODE), sctrB(NNODE*MCRD)
      
!     for CPU TIME
      Double precision t2,t1,TOTAL_TIME_2, TOTAL_TIME_1
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
      call CPU_TIME(TOTAL_TIME_1)

!     Assign common variables
      Ukv = Ukv_e
      Nkv = Nkv_e
      Jpqr = Jpqr_e
      Nijk = Nijk_e
      weight = weight_e  
      
!     Initialize K and F to cload
      K = zero
      F = zero
      nb_dload = nb_load-nb_cload
      do num_load = nb_dload+1, nb_load
         do num_cp = 1,load_target_nbelem(num_load)
            ddl = (indDLoad(num_load,num_cp)-1)*MCRD + JDLType(num_load)
            F(ddl) = ADLMAG(num_load)
         enddo
      enddo
      
      
!     Compute DOF information
      NDOFEL = MCRD*NNODE
      nb_dof_tot = nb_cp*MCRD
      call getindDOF(MCRD,nb_bc,bc_target,bc_target_nbelem, bc_values,
     &     nb_dof_tot,nb_dof_bloq,nb_dof_free,ind_dof_bloq,ind_dof_free,
     &     COUPLG_flag,BNDSTRIP_flag)
      
      
      
!     Print informations
 1    format(X,A30,I6)
      write(*,*) ''
      write(*,'(A)') 
     &     '...........................................................'
      write(*,*) 'ISOGEOMETRIC ANALYSIS'
      write(*,*) ''
      write(*,*) '.....................................'
      write(*,*) ' Summary of principal data'
      write(*,*) '-------------------------------------'
      write(*,*) ' Element Type : ', TENSOR,', ',ELT_TYPE
      write(*,1) 'Dimension                    ', MCRD
      write(*,1) 'Number of patches            ', nb_patch
      write(*,1) 'Number of elements           ', nb_elem
      write(*,1) 'Number of CP by element      ', NNODE
      write(*,1) 'Number of DOF by element     ', NDOFEL
      write(*,1) 'Number of IP by element      ', NBINT
      write(*,*) '--'
      write(*,1) 'Total number of CP           ', nb_cp
      write(*,1) 'Total number of DOF          ', nb_dof_tot
      write(*,1) 'Number of Boundary conditions', nb_bc
      write(*,1) 'Number of fixed DOF          ', nb_dof_bloq
      write(*,1) 'Number of free DOF           ', nb_dof_free
      write(*,1) 'Number of Loads              ', nb_load
      write(*,*) '-------------------------------------'
      write(*,*) ''
      
c     Fin Initialisation ...............................................
c     
c     
c     
c     
c     Debut Analyse ....................................................
c      
c     Debut Assemblage .................................................
      
      write(*,'(A)',advance='no') ' Begin assembly...'
      call CPU_TIME(t1)
      
      JELEM = 0
      Do NumPatch = 1,nb_patch  ! boucle patch
         Do num_elem = 1,nb_elem_patch(NumPatch) ! boucle element
            
            JELEM = JELEM + 1
            do i = 1,NNODE
               do j = 1,MCRD
                  COORDS_elem(j,i) = COORDS3D(IEN(JELEM,i),j)
               enddo
            enddo
            
c     Compute elementary matrix and load vector
            RHS = zero
            AMATRX = zero
            dvol_save = zero
            MAT_patch(:) = MATERIAL_PROPERTIES(NumPatch,:)
            JPROPS_patch = JPROPS(NumPatch)
            PROPS_patch(:) = PROPS(NumPatch,:)
            if (ELT_TYPE == 'U1') then
               ! 'Element classique solide'
               call UELMAT(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch ,
     1              COORDS_elem,TENSOR,MAT_patch,nb_dload,indDLoad,
     2              load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX  ,
     3              dvol_save)
            
            elseif (ELT_TYPE == 'U2') then
               ! 'Element plaque'
               call UELMAT2(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,
     1              COORDS_elem,TENSOR,MAT_patch,PROPS_patch,
     2              JPROPS_patch,nb_dload,indDLoad,load_target_nbelem,
     3              JDLType,ADLMAG,RHS,AMATRX,dvol_save)
            
            elseif (ELT_TYPE == 'U3') then
               ! 'Element coque'
               call UELMAT3(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,
     1              COORDS_elem,TENSOR,MAT_patch,PROPS_patch,
     2              JPROPS_patch,nb_dload,indDLoad,load_target_nbelem,
     3              JDLType,ADLMAG,RHS,AMATRX,dvol_save)               
            else
               print*, 'Element' // ELT_TYPE // ' not availble.'
               
            endif
            
c     Assemble AMATRX to global stiffness matrix K    
            sctr = IEN(JELEM,:)
            Do kk = 1,MCRD
               j = 1
               Do i = kk,MCRD*NNODE,MCRD
                  sctrB(i) = MCRD*(sctr(j)-1) + kk
                  j = j+1
               Enddo
            Enddo
            
            Do i = 1,NDOFEL
               F(sctrB(i)) = F(sctrB(i)) + RHS(i)
               Do j = 1,NDOFEL
                  K(sctrB(i),sctrB(j)) = K(sctrB(i),sctrB(j))
     1                 + AMATRX(i,j)
               Enddo
            Enddo
         Enddo
      Enddo
      call CPU_TIME(t2)
      write(*,*)'Assembly DONE.'
      write(*,'(3X,A,X,F8.2,X,A)') 'CPU Time for assembly : ',t2-t1,'s'
      write(*,*) ''
      
c     Fin Assemblage ...................................................
c     
c     Debut Resolution .................................................
      
c     Bending strip
      If (BNDSTRIP_flag) then
         write(*,'(A)',advance='no') ' Begin add bndstrip...'
         call CPU_TIME(t1)
         
         MAT_bndstrip(:)   = zero
         PROPS_bndstrip(:) = zero
         Do numPatch = 1,nb_patch
            MAT_bndstrip(1) =
     &           max(MAT_bndstrip(1), MATERIAL_PROPERTIES(NumPatch,1))
            PROPS_bndstrip(2) =
     &           max(PROPS_bndstrip(2), PROPS(numPatch,2))
         Enddo
         MAT_bndstrip(1)   = MAT_bndstrip(1)  / float(nb_patch)
         PROPS_bndstrip(2) = PROPS_bndstrip(2)/ float(nb_patch)
         PROPS_bndstrip(1) = float(nb_patch+1)
         
         call add_bendingstrips(K,COORDS3D,nb_cp,MCRD,nb_bc,bc_target,
     &        bc_target_nbelem,bc_values,MAT_bndstrip,PROPS_bndstrip,
     &        JPROPS(1))
         call CPU_TIME(t2)
         write(*,*)'BndStrip DONE.'
         write(*,'(3X,A,X,F8.2,X,A)') 'CPU Time bndstrip : ',t2-t1,'s'
         write(*,*) ''
      Endif
      
      
c     Couplage de DDL si necessaire
      If (COUPLG_flag) then
         write(*,'(A)',advance='no') ' Begin coupling...'
         call CPU_TIME(t1)
         call coupling(K,F,nb_cp,MCRD,nb_bc,bc_target,bc_target_nbelem,
     &        bc_values,COORDS3D)
         call CPU_TIME(t2)
         write(*,*)'Coupling DONE.'
         write(*,'(3X,A,X,F8.2,X,A)') 'CPU Time coupling : ',t2-t1,'s'
         write(*,*) ''
      Endif
      
      
c     Build system to solve
      allocate(K_inv(nb_dof_free,nb_dof_free))
      allocate(F_inv(nb_dof_free))
      allocate(U_inv(nb_dof_free))
      write(*,*)'Begin resolution...'
      write(*,'(A)',advance='no')'  Apply BCs...'
      
c      call apply_BC(K,F,nb_cp,MCRD,nb_bc,bc_target,bc_target_nbelem,
c     1     bc_values,nb_dof_tot,nb_dof_bloq,nb_dof_free,ind_dof_bloq,
c     2     ind_dof_free,U)
      call apply_dispBC(K,F,nb_cp,MCRD,nb_bc,bc_target,bc_target_nbelem,
     &     bc_values,nb_dof_tot,U)
      write(*,*) 'DONE.'
      
      do i = 1,nb_dof_free
         F_inv(i) = F(ind_dof_free(i))
         do j = 1,nb_dof_free
            K_inv(i,j) = K(ind_dof_free(i),ind_dof_free(j))
         enddo
      enddo
      
      
c     Resolution
      call resolution(K_inv,F_inv,nb_dof_free,U_inv)
      write(*,*)'Resolution DONE.'
      write(*,*)''
      
c     Reconstruction of the solution
      call reconstruction(sol,U,U_inv,COORDS3D,ind_dof_free,nb_bc,
     &     bc_target,bc_target_nbelem,bc_values,COUPLG_flag,MCRD,nb_cp,
     &     nb_dof_free)
      
      
c     Fin Analyse ......................................................
c     
c     Ecriture resultats ...............................................
      
      if (.False.) then
         open(90, file="output/raideur.post", form="formatted")
         open(91, file="output/secondmembre.post", form="formatted")
         open(92, file="output/solution.post", form="formatted")
         open(93, file="output/raideurinv.post", form="formatted")
         open(94, file="output/secondmembreinv.post", form="formatted")
         open(95, file="output/numerotation.post", form="formatted")
         
         do i = 1,nb_dof_tot
            write(90,*)K(i,:)
            write(91,*)F(i)
            write(92,*)U(i)
         enddo
         do i = 1,nb_dof_free
            write(93,*)K_inv(i,:)
            write(94,*)F_inv(i)
         enddo
         
!     Numerotation globale de K_inv
         do i = 1,nb_dof_free
            k1 = floor(real(ind_dof_free(i)+1)/MCRD)
            k2 = mod(ind_dof_free(i),MCRD)
            if (k2==0) then
               write(95,*)k1,MCRD
            else 
               write(95,*)k1,k2
            endif
         enddo
         
         close(90)
         close(91)
         close(92)
         close(93)
         close(94)
         close(95)
      endif
      
c     Fin ecriture .....................................................
      
      deallocate(K_inv, F_inv, U_inv)
      
      call CPU_TIME(TOTAL_TIME_2)
      write(*,'(X,A,F6.2,A)')
     1     'Total time of analysis (CPU TIME) :',
     2     TOTAL_TIME_2 - TOTAL_TIME_1,' s'
      write(*,'(A)') 
     &     '...........................................................'
      write(*,*) ''
      
      
      
      end subroutine iga_linmat_lindef_static
