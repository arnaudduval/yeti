!! Compute material behaviour tangent matrix ddsdde
!! WARNING : only isotropic elastic behaviour is taken into account and 2D analysis
!! Conventions :
!! 2D
!! stress : {sigma_11, sigma_22, √2*sigma_12,

!! sigma_111, sigma_221, √2*sigma_121, √2sigma_112, sigma_222, √2*sigma_122

!! sigma_1111, sigma_2221, √3*sigma_1121, √3*sigma_1221, sigma_1112
!! sigma_2222,  √3*sigma_1122, √3*sigma_1222}


!! strain : {eps_11, eps_22, √2*eps_12,

!! eps_111, eps_221, √2*eps_121, √2eps_112, eps_222, √2*eps_122

!! eps_1111, eps_2221, √3*eps_1121, √3*eps_1221, eps_1112, eps_2222, √3*eps_1122, √3*eps_1222}


!! Material properties are given in the following order :
!! Definition from
!! R.D. Mindlin
!! Second gradient of strain and suirface tension in linear elasticity
!! Int. J. Solids Strucrures, 1965, Vol. 1, pp. 417 to 438.
!! [1] : lambda
!! [2] : mu
!! [3->7] : a1 -> a5
!! [8->15] : b0 -> b7
!! [16->18] : c1 -> c3

subroutine material_lib_HO(MATERIAL_PROPERTIES, n_mat_props, TENSOR, MCRD, ddsdde)
    use parameters
    implicit none

    double precision, intent(in) :: MATERIAL_PROPERTIES
    integer, intent(in) :: n_mat_props
    character(len=*), intent(in) :: TENSOR
    integer, intent(in) :: MCRD

    double precision, intent(out) :: ddsdde

    dimension :: MATERIAL_PROPERTIES(n_mat_props)

    double precision :: ddsdde1, ddsdde2, ddsdde3, ddsdde4, ddsdde5

    dimension ddsdde1(3, 3)
    dimension ddsdde2(6, 6)
    dimension ddsdde3(8, 8)
    dimension ddsdde4(3, 8)
    dimension ddsdde5(8, 3)
    dimension ddsdde(17, 17)

    double precision :: lambda, mu, a1, a2, a3, a4, a5, c1, c2, c3, b1, b2, b3, b4, b5, b6, b7, b0
    double precision :: Sa, Ta
    double precision :: Nb, Qb, Rb, Sb, Tb, Ub, Vb, Wb, Xb
    double precision :: Sc
    integer i, j

    !! Initialization

    lambda  = MATERIAL_PROPERTIES(1)
    mu = MATERIAL_PROPERTIES(2)

    a1 = MATERIAL_PROPERTIES(3)
    a2 = MATERIAL_PROPERTIES(4)
    a3 = MATERIAL_PROPERTIES(5)
    a4 = MATERIAL_PROPERTIES(6)
    a5 = MATERIAL_PROPERTIES(7)

    b0 = MATERIAL_PROPERTIES(8)
    b1 = MATERIAL_PROPERTIES(9)
    b2 = MATERIAL_PROPERTIES(10)
    b3 = MATERIAL_PROPERTIES(11)
    b4 = MATERIAL_PROPERTIES(12)
    b5 = MATERIAL_PROPERTIES(13)
    b6 = MATERIAL_PROPERTIES(14)
    b7 = MATERIAL_PROPERTIES(15)

    c1 = MATERIAL_PROPERTIES(16)
    c2 = MATERIAL_PROPERTIES(17)
    c3 = MATERIAL_PROPERTIES(18)

    !! *****************KO**********************************

    ddsdde1(:,:) = zero

    ddsdde1(1, 1) = lambda + two * mu
    ddsdde1(2, 2) = lambda + two * mu
    ddsdde1(1, 2) = lambda
    ddsdde1(2, 1) = lambda
    ddsdde1(3, 3) = two * mu

    !! *****************   Ka    **********************************

    Sa  = a5 + a4 + a3 + a2 + a1
    Ta  = a5 + two * a4 + a1

    ddsdde2(:,:) = zero
    !! Define triangular sup
    ddsdde2(1, 1) = two * Sa
    ddsdde2(1, 2) = two * a3 + a2
    ddsdde2(1, 6) = (a2 + two * a1) * Inv_Ractwo
    ddsdde2(2, 2) = two * (a4 + a3)
    ddsdde2(2, 6) = (two * a5 + a2) * Inv_Ractwo
    ddsdde2(3, 3) = Ta
    ddsdde2(3, 4) = (two * a5 + a2) * Inv_Ractwo
    ddsdde2(3, 5) = (a2 + two * a1) * Inv_Ractwo
    ddsdde2(4, 4) = two * (a4 + a3)
    ddsdde2(4, 5) = two * a3 + a2
    ddsdde2(5, 5) = two * Sa
    ddsdde2(6, 6) = Ta

    !! WARNING : pas performant : copie beaucoup de zeros
    do j = 1, 6
        do i = j + 1, 6
            ddsdde2(i, j) = ddsdde2(j, i)
        end do
    end do

    !! *****************  Kb   **********************************

    Nb  =  two * b4 + b3 + four * b2
    Qb  =  two * b7 + b3 + b1
    Rb  =  two * (b6 + b5)
    Sb  =  b7 + b6 + b5 + b4 + b3 + b2 + b1
    Tb  =  b7 + three * b6 + b5 + b4 + b2 + b1
    Ub  =  two * b7 + three * b6 + b5 + b3 + two * b2
    Vb  =  two * b5 + two * b4 + b3 + two * b1
    Wb  =  b3 + two * b2 + two * b1
    Xb  =  two * b7 + two * b4 + b3

    ddsdde3(:,:) = zero

    ddsdde3(1, 1) = two * Sb
    ddsdde3(1, 4) = Vb * Inv_Racthree
    ddsdde3(1, 6) = two * b1
    ddsdde3(1, 7) = Wb * Inv_Racthree
    ddsdde3(2, 2) = Rb
    ddsdde3(2, 3) = (two * b5 + b3) * Inv_Racthree
    ddsdde3(2, 5) = two * b4
    ddsdde3(2, 8) = Xb * Inv_Racthree
    ddsdde3(3, 3) = two * Ub * third
    ddsdde3(3, 5) = Xb * Inv_Racthree
    ddsdde3(3, 8) = Nb * third
    ddsdde3(4, 4) = two * Tb * third
    ddsdde3(4, 6) = Wb * Inv_Racthree
    ddsdde3(4, 7) = two * Qb * third
    ddsdde3(5, 5) = Rb
    ddsdde3(5, 8) = (two * b5 + b3) * Inv_Racthree
    ddsdde3(6, 6) = two * Sb
    ddsdde3(6, 7) = Vb * Inv_Racthree
    ddsdde3(7, 7) = two * Tb * third
    ddsdde3(8, 8) = two * Ub * third


    do j = 1, 8
        do i = j + 1, 8
            ddsdde3(i, j) = ddsdde3(j, i)
        end do
    end do

    !! *****************  Kc & Kc transpose  **********************************

    Sc  = c3 + c2 + c1

    ddsdde4(:,:) = zero

    ddsdde4(1, 1) = Sc
    ddsdde4(1, 4) = (c3 + c1) * Inv_Racthree
    ddsdde4(1, 6) = c1
    ddsdde4(1, 7) = (c2 + c1) * Inv_Racthree
    ddsdde4(2, 1) = c1
    ddsdde4(2, 4) = (c2 + c1) * Inv_Racthree
    ddsdde4(2, 6) = Sc
    ddsdde4(2, 7) = (c3 + c1) * Inv_Racthree
    ddsdde4(3, 2) = c3 * Inv_Ractwo
    ddsdde4(3, 3) = (c3 + two * c2) * Inv_RacSix
    ddsdde4(3, 5) = c3 * Inv_Ractwo
    ddsdde4(3, 8) = (c3 + two * c2) * Inv_RacSix

    ddsdde5(:,:) = zero
    ddsdde5(:,:) = TRANSPOSE(ddsdde4(:,:))

    !! *****************       global C       ****************************************

    ddsdde(:,:)  = zero

    ddsdde(1:3, 1:3)     = ddsdde1(:,:)
    ddsdde(1:3, 10:17)   = ddsdde4(:,:)
    ddsdde(4:9, 4:9)     = ddsdde2(:,:)
    ddsdde(10:17, 1:3)   = ddsdde5(:,:)
    ddsdde(10:17, 10:17) = ddsdde3(:,:)

    ! write(*,*) "C1"
    ! do i = 1, 3
    ! write(*,*) ddsdde1(i, :)
    ! enddo
    ! write(*,*) "C2"
    ! do i = 1, 6
    ! write(*,*) ddsdde2(i, :)
    ! enddo
    ! write(*,*) "C3"
    ! do i = 1, 8
    ! write(*,*) ddsdde3(i, :)
    ! enddo
    ! write(*,*) "C4"
    ! do i = 1, 3
    ! write(*,*) ddsdde4(i, :)
    ! enddo
    ! write(*,*) "C5"
    ! do i = 1, 8
    ! write(*,*) ddsdde5(i, :)
    ! enddo
    ! write(*,*) "ddsdde"
    ! do i = 1, 17
    ! write(*,*) ddsdde(i, :)
    ! enddo

end subroutine material_lib_HO


!! Same as material_lib_HO but with only 1st gradient of strain
subroutine material_lib_HO_1stG(MATERIAL_PROPERTIES, n_mat_props, TENSOR, MCRD, ddsdde)

    use parameters

    implicit none

    double precision, intent(in) :: MATERIAL_PROPERTIES
    integer, intent(in) :: n_mat_props
    character(len=*), intent(in) :: TENSOR
    integer, intent(in) :: MCRD

    double precision, intent(out) :: ddsdde

    dimension :: MATERIAL_PROPERTIES(n_mat_props)

    double precision :: ddsdde1, ddsdde2

    dimension ddsdde1(3, 3)
    dimension ddsdde2(6, 6)
    dimension ddsdde(9, 9)

    double precision :: lambda, mu, a1, a2, a3, a4, a5, c1, c2, c3, b1, b2, b3, b4, b5, b6, b7, b0
    double precision :: Sa, Ta
    double precision :: Nb, Qb, Rb, Sb, Tb, Ub, Vb, Wb, Xb
    double precision :: Sc
    integer i, j

    !! Initialization

    lambda  = MATERIAL_PROPERTIES(1)
    mu = MATERIAL_PROPERTIES(2)

    a1 = MATERIAL_PROPERTIES(3)
    a2 = MATERIAL_PROPERTIES(4)
    a3 = MATERIAL_PROPERTIES(5)
    a4 = MATERIAL_PROPERTIES(6)
    a5 = MATERIAL_PROPERTIES(7)

    b0 = MATERIAL_PROPERTIES(8)
    b1 = MATERIAL_PROPERTIES(9)
    b2 = MATERIAL_PROPERTIES(10)
    b3 = MATERIAL_PROPERTIES(11)
    b4 = MATERIAL_PROPERTIES(12)
    b5 = MATERIAL_PROPERTIES(13)
    b6 = MATERIAL_PROPERTIES(14)
    b7 = MATERIAL_PROPERTIES(15)

    c1 = MATERIAL_PROPERTIES(16)
    c2 = MATERIAL_PROPERTIES(17)
    c3 = MATERIAL_PROPERTIES(18)

    !! *****************KO**********************************

    ddsdde1(:,:) = zero

    ddsdde1(1, 1) = lambda + two * mu
    ddsdde1(2, 2) = lambda + two * mu
    ddsdde1(1, 2) = lambda
    ddsdde1(2, 1) = lambda
    ddsdde1(3, 3) = two * mu

    !! *****************   Ka    **********************************

    Sa  = a5 + a4 + a3 + a2 + a1
    Ta  = a5 + two * a4 + a1

    ddsdde2(:,:) = zero
    !! Define triangular sup
    ddsdde2(1, 1) = two * Sa
    ddsdde2(1, 2) = two * a3 + a2
    ddsdde2(1, 6) = (a2 + two * a1) * Inv_Ractwo
    ddsdde2(2, 2) = two * (a4 + a3)
    ddsdde2(2, 6) = (two * a5 + a2) * Inv_Ractwo
    ddsdde2(3, 3) = Ta
    ddsdde2(3, 4) = (two * a5 + a2) * Inv_Ractwo
    ddsdde2(3, 5) = (a2 + two * a1) * Inv_Ractwo
    ddsdde2(4, 4) = two * (a4 + a3)
    ddsdde2(4, 5) = two * a3 + a2
    ddsdde2(5, 5) = two * Sa
    ddsdde2(6, 6) = Ta

    !! WARNING : pas performant : copie beaucoup de zeros
    do j = 1, 6
        do i = j + 1, 6
            ddsdde2(i, j) = ddsdde2(j, i)
        end do
    end do

    !! *****************       global C       ****************************************

    ddsdde(:,:)  = zero

    ddsdde(1:3, 1:3)     = ddsdde1(:,:)
    ddsdde(4:9, 4:9)     = ddsdde2(:,:)

end subroutine material_lib_HO_1stG








