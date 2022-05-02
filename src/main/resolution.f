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

!>    Routine resolvant le systeme lineaire KU=F avec DPOSV de Lapack
!>    \param K_inv : matrice a inverser
!>    \param F_inv : second membre
!>    \param U_inv : solution
!>    \param nb_dof_free : taille de la matrice K_inv







      subroutine resolution(K_inv,F_inv,nb_dof_free,U_inv)
      implicit none

      INTEGER   INFO, LDA, LDB, N, NRHS

      integer, intent(in) :: nb_dof_free

        
      

      double precision, intent(in) :: K_inv, F_inv
      dimension K_inv(nb_dof_free,nb_dof_free)
      dimension F_inv(nb_dof_free)                                                                                                
      
      double precision, intent(inout) :: U_inv
      dimension U_inv(nb_dof_free)

      integer, dimension(:), allocatable :: IPIV
      double precision, dimension(:,:),allocatable :: K_invc
      double precision, dimension(:),allocatable :: F_invc 
    
      double precision t1,t2

      CHARACTER  :: UPLO
      
       
      N = nb_dof_free
      NRHS = 1
      LDA = N
      LDB = N

      allocate(K_invc(LDA,N))
      allocate(F_invc(LDB)) 
           
      K_invc=K_inv
      F_invc=F_inv
      
     

      UPLO='U'

      call CPU_TIME(t1)
      call DPOSV(UPLO,N,NRHS,K_invc,LDA,F_invc,LDB,INFO)
      call CPU_TIME(t2)


      if (INFO==0) then
         write(*,*)'  The linear system has been resolved successfully'
         write(*,'(3X,A,X,F8.2,X,A)')
     &        'CPU TIME for solving : ', t2-t1,' s'
      else
         write(*,*)'  Failure during resolution of linear 
     1        system with dposv'
         
      endif

      
      U_inv = F_invc
      
      
      deallocate(K_invc, F_invc)
      
      
      end subroutine resolution
