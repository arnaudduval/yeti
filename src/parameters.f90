!! Copyright 2016-2020 Thibaut Hirschler

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

module parameters
  implicit none
  double precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
  double precision, parameter :: three = 3.d0, four = 4.d0
  double precision, parameter :: half=0.5d0, third=1.d0/3.d0
  double precision, parameter :: ractwo = sqrt(2.d0)
  double precision, parameter :: racthree = sqrt(3.d0)
  double precision, parameter :: inv_ractwo = (1.0d0)/sqrt(2.0d0)
  double precision, parameter :: inv_racthree = (1.0d0)/sqrt(3.0d0)
  double precision, parameter :: inv_racsix = (1.0d0)/sqrt(3.0d0)

end module parameters
