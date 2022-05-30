!! Copyright 2022 Arnaud Duval

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


!! Routine for Python binding defining number of integration points explicitely
!! Parameters:
subroutine GaussPts(norder,dim,GaussWtCoord,iFace,n_gps)
    implicit none

    integer, intent(in) :: norder, dim, iFace,n_gps
    double precision, intent(out) :: GaussWtCoord
    dimension GaussWtCoord(dim+1,n_gps)

    call Gauss(norder, dim, GaussWtCoord, iFace)

end subroutine GaussPts
