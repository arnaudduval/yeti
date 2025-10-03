module linalg3x3
contains
    function solve3x3(A, b) result(x)
        implicit none
        real(8), intent(in) :: A(3,3)
        real(8), intent(in) :: b(3)
        real(8) :: x(3)
        real(8) :: detA, det

        ! Déterminant de A
        detA = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
            A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
            A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

        if (detA == 0.0d0) stop "Matrice singulière (pas de solution unique)."

        ! x(1) : remplacer la première colonne par b
        det = b(1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
            A(1,2)*(b(2)*A(3,3) - A(2,3)*b(3)) + &
            A(1,3)*(b(2)*A(3,2) - A(2,2)*b(3))
        x(1) = det / detA

        ! x(2) : remplacer la deuxième colonne par b
        det = A(1,1)*(b(2)*A(3,3) - A(2,3)*b(3)) - &
            b(1)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
            A(1,3)*(A(2,1)*b(3) - b(2)*A(3,1))
        x(2) = det / detA

        ! x(3) : remplacer la troisième colonne par b
        det = A(1,1)*(A(2,2)*b(3) - b(2)*A(3,2)) - &
            A(1,2)*(A(2,1)*b(3) - b(2)*A(3,1)) + &
            b(1)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
        x(3) = det / detA

    end function solve3x3

end module linalg3x3