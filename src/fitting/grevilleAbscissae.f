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

      Subroutine grevilleAbscissae(Uknot,p,n,nknots,greAbsc)
      
      Implicit none
      
      ! Inputs
      Integer,          intent(in) :: p,n,nknots
      Double precision, intent(in) :: Uknot
      dimension Uknot(nknots)

      ! Output
      Double precision, intent(out):: greAbsc
      dimension greAbsc(n)
      
      ! Local var
      Integer :: i,k
      
!     --
      greAbsc(:) = 0.0d0
      Do i = 1,n
         Do k = 1,p
            greAbsc(i) = greAbsc(i) + Uknot(i+k)
         Enddo
      Enddo
      greAbsc(:) = greAbsc(:)/dble(p)
      
      End subroutine grevilleAbscissae
      
