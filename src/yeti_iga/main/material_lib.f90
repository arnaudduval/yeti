!! Copyright 2016-2018 Thibaut Hirschler
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

!! Compute material behaviour tangent matrix ddsdde
!! WARNING : only isotropic elastic behaviour is taken into account
!! Conventions :
!! 3D
!! stress : {sigma_11 sigma_22 sigma_33 sigma_12 sigma_13 sigma_33}
!! strain : {eps_11   eps_22   eps_33   2*eps_12 2*eps_13 2*eps_33}
!! 2D
!! stress : {sigma_11 sigma_22 sigma_33 sigma_12}
!! strain : {eps_11   eps_22   eps_33   2*eps_12}


subroutine material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde)

    implicit none

    double precision, intent(in) :: MATERIAL_PROPERTIES
    character(len=*), intent(in) :: TENSOR
    integer, intent(in) :: MCRD

    double precision, intent(out) :: ddsdde

    dimension MATERIAL_PROPERTIES(2)
    dimension ddsdde(2 * MCRD, 2 * MCRD)

    double precision zero, one, two
    parameter(zero = 0.0D0, one = 1.0D0, two = 2.0D0)

    double precision :: E, nu, lambda, mu, coef
    integer i, j

    !! Initialization
    ddsdde(:,:) = zero

    E  = MATERIAL_PROPERTIES(1)
    nu = MATERIAL_PROPERTIES(2)
    lambda = (E * nu) / ((one + nu) * (one - two * nu))
    mu = E / (two * (one + nu))

    if (TENSOR == 'PSTRAIN') then
        ddsdde(1, 1) = lambda + two * mu
        ddsdde(2, 2) = lambda + two * mu
        ddsdde(3, 3) = lambda + two * mu
        ddsdde(1, 2) = lambda
        ddsdde(2, 1) = lambda
        ddsdde(1, 3) = lambda
        ddsdde(3, 1) = lambda
        ddsdde(2, 3) = lambda
        ddsdde(3, 2) = lambda
        ddsdde(4, 4) = mu

    else if (TENSOR == 'PSTRESS') then
        coef = E / (one - nu * nu)
        ddsdde(1, 1) = coef
        ddsdde(2, 2) = coef
        ddsdde(1, 2) = nu * coef
        ddsdde(2, 1) = nu * coef
        ddsdde(4, 4) = coef * (one - nu) * 0.5d0

    else if (TENSOR == 'THREED') then
        ddsdde(1, 1) = lambda + two * mu
        ddsdde(2, 2) = lambda + two * mu
        ddsdde(3, 3) = lambda + two * mu
        ddsdde(1, 2) = lambda
        ddsdde(2, 1) = lambda
        ddsdde(1, 3) = lambda
        ddsdde(3, 1) = lambda
        ddsdde(2, 3) = lambda
        ddsdde(3, 2) = lambda
        ddsdde(4, 4) = mu
        ddsdde(5, 5) = mu
        ddsdde(6, 6) = mu

    else if (TENSOR == 'THREED_bnd') then
        ddsdde(:,:) = zero
        ddsdde(5, 5) = mu * 10.d0 ** 10
        ddsdde(6, 6) = mu * 10.d0 ** 10

    else
        write(*, *) 'TENSOR ', TENSOR, ' NOT AVAILABLE'
    end if

end subroutine material_lib
