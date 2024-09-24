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

      SUBROUTINE getGPsOnParamSpace(XIbar,XI,BI,dim,MCRD,NNODE,NBINT,
     &     COORDS)
      
      use parameters
      use nurbspatch

      Implicit None
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: dim,MCRD,NNODE,NBINT
      Double precision, intent(in) :: COORDS
      dimension COORDS(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: XIbar,XI,BI
      dimension XIbar(4,NBINT),XI(3,NBINT),BI(3,dim,NBINT)
      
c     Local variables :
c     ---------------
      
!     Nurbs basis fcts
      Double precision :: R,dRdxi
      dimension R(NNODE),dRdxi(NNODE,3)
      
!     For gauss points
      Integer          :: NbPtInt
      Double precision :: GaussPdsCoord,PtGauss
      dimension GaussPdsCoord(dim+1,NBINT),PtGauss(dim)
      
!     for loops
      Integer :: n,i,j

C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      
      XIbar(:,:) = zero
      XI(:,:)    = zero
      BI(:,:,:)  = zero
      
!     Defining Gauss points coordinates and weights
      NbPtInt = int( NBINT**(1.0/float(dim)) )
      if (NbPtInt**dim<NBINT) NbPtInt = NbPtInt + 1
      call Gauss(NbPtInt,dim,GaussPdsCoord,0)
      
c     
c     ..................................................................
c
C     Computation :
c     -----------
      
      Do n = 1,NBINT
         
!     - get XIbar
         PtGauss(:) = GaussPdsCoord(2:,n)
         XIbar(4,n) = GaussPdsCoord(1 ,n)
         Do i = 1,dim
            XIbar(i,n) = ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
            XIbar(4,n) = XIbar(4,n)
     &                 * ( Ukv_elem(2,i) - Ukv_elem(1,i))*0.5d0
         End do
         
!     - evaluate basis functions
         call evalnurbs(XIbar(:3,n),R(:),dRdxi(:,:))
         
!     - get position
         Do i = 1,NNODE
            XI(:MCRD,n) = XI(:MCRD,n) + R(i)*COORDS(:,i)
         Enddo

!     - get covariant basis vectors
         Do j = 1,dim
            Do i = 1,NNODE
               BI(:MCRD,j,n) = BI(:MCRD,j,n) + dRdxi(i,j)*COORDS(:,i)
            Enddo
         Enddo
         
      Enddo
      
c     
c     ..................................................................
c
      
      END SUBROUTINE getGPsOnParamSpace












      



      
      SUBROUTINE getGPsOnParamSpace_wBasis(XIbar,XI,BI,RI,dRdxi,dim,
     &     MCRD,NNODE,NBINT,COORDS)
      
      use parameters
      use nurbspatch

      Implicit None
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: dim,MCRD,NNODE,NBINT
      Double precision, intent(in) :: COORDS
      dimension COORDS(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: XIbar,XI,BI,RI,dRdxi
      dimension XIbar(4,NBINT),XI(3,NBINT),BI(3,dim,NBINT),
     &     RI(NNODE,NBINT),dRdxi(NNODE,dim,NBINT)
      
c     Local variables :
c     ---------------
      
!     Nurbs basis fcts
      Double precision :: R,tempdRdxi
      dimension R(NNODE),tempdRdxi(NNODE,3)
      
!     For gauss points
      Integer          :: NbPtInt
      Double precision :: GaussPdsCoord,PtGauss
      dimension GaussPdsCoord(dim+1,NBINT),PtGauss(dim)
      
!     for loops
      Integer :: n,i,j

C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      
      XIbar(:,:)   = zero
      XI(:,:)      = zero
      BI(:,:,:)    = zero
      RI(:,:)      = zero
      dRdxi(:,:,:) = zero
      
!     Defining Gauss points coordinates and weights
      NbPtInt = int( NBINT**(1.0/float(dim)) )
      if (NbPtInt**dim<NBINT) NbPtInt = NbPtInt + 1
      call Gauss(NbPtInt,dim,GaussPdsCoord,0)
      
c     
c     ..................................................................
c
C     Computation :
c     -----------
      
      Do n = 1,NBINT
         
!     - get XIbar
         PtGauss(:) = GaussPdsCoord(2:,n)
         XIbar(4,n) = GaussPdsCoord(1 ,n)
         Do i = 1,dim
            XIbar(i,n) = ((Ukv_elem(2,i) - Ukv_elem(1,i))*PtGauss(i)
     &                 +  (Ukv_elem(2,i) + Ukv_elem(1,i)) ) * 0.5d0
            XIbar(4,n) = XIbar(4,n)
     &                 * ( Ukv_elem(2,i) - Ukv_elem(1,i))*0.5d0
         End do
         
!     - evaluate basis functions
         call evalnurbs(XIbar(:3,n),R(:),tempdRdxi(:,:))
         dRdxi(:,:,n) = tempdRdxi(:,:dim)
         RI(:,n) = R(:)

!     - get position
         Do i = 1,NNODE
            XI(:MCRD,n) = XI(:MCRD,n) + R(i)*COORDS(:,i)
         Enddo

!     - get covariant basis vectors
         Do j = 1,dim
            Do i = 1,NNODE
               BI(:MCRD,j,n) = BI(:MCRD,j,n) +tempdRdxi(i,j)*COORDS(:,i)
            Enddo
         Enddo
         
      Enddo
      
c     
c     ..................................................................
c
      
      END SUBROUTINE getGPsOnParamSpace_wBasis
