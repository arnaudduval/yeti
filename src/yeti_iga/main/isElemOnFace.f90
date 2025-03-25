!! Copyright 2021 Arnaud Duval

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

!! Returns true if given element Nijk is on given face
function IsElemOnFace(face, Nijk, Jpqr, Nkv, dim)
    implicit none
    integer :: face, Nijk, Jpqr, Nkv
    integer :: dim
    !! TODO adapt array dimension to problem dimension
    ! dimension Nijk(dim), Jpqr(dim), Nkv(dim)
    dimension Nijk(3), Jpqr(3), Nkv(3)
    logical :: IsElemOnFace

    IsElemOnFace = .false.
    if ((face == 1) .and. (Nijk(1) == (Jpqr(1) + 1))) IsElemOnFace = .true.
    if ((face == 2) .and. (Nijk(1) == (Nkv(1) - Jpqr(1) - 1))) IsElemOnFace = .true.
    if ((face == 3) .and. (Nijk(2) == (Jpqr(2) + 1))) IsElemOnFace = .true.
    if ((face == 4) .and. (Nijk(2) == (Nkv(2) - Jpqr(2) - 1))) IsElemOnFace = .true.
    if (dim > 2) then
        if ((face == 5) .and. (Nijk(3) == (Jpqr(3) + 1))) IsElemOnFace = .true.
        if ((face == 6) .and. (Nijk(3) == (Nkv(3) - Jpqr(3) - 1))) IsElemOnFace = .true.
    end if
end function IsElemOnFace

!! Return the number of elements on the given face of a patch
!! WARNING must be tested with repeated knots
function NbElemOnFace(face, Jpqr, Nkv)
    implicit none
    integer :: face, Jpqr, Nkv
    dimension Jpqr(3), Nkv(3)
    integer :: NbElemOnFace

    select case (face)
        case (1, 2)
        NbElemOnFace = (Nkv(2) - 2 * Jpqr(2) - 1) * (Nkv(3) - 2 * Jpqr(3) - 1)
        case (3, 4)
        NbElemOnFace = (Nkv(1) - 2 * Jpqr(1) - 1) * (Nkv(3) - 2 * Jpqr(3) - 1)
        case (5, 6)
        NbElemOnFace = (Nkv(1) - 2 * Jpqr(1) - 1) * (Nkv(2) - 2 * Jpqr(2) - 1)
    end select

end function NbElemOnFace
