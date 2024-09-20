subroutine UELTANMAT_byCP(NDOFEL, MCRD, NNODE, JELEM, NBINT, COORDS,            &
        &   TENSOR, MATERIAL_PROPERTIES, Uelem, AMATRX)

    use parameters

    implicit None

    !! Input arguments
    !! ---------------
    integer, intent(in) :: NDOFEL, MCRD, NNODE, JELEM, NBINT
    character(len=*), intent(in) :: TENSOR
    double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, Uelem
    dimension COORDS(MCRD, NNODE), MATERIAL_PROPERTIES(2),  Uelem(MCRD,NNODE)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: AMATRX
    dimension AMATRX(MCRD, MCRD, NNODE * (NNODE + 1) / 2)


    !! Local variables
    !! ---------------

    !! Gauss points
    integer :: NbPtInt, n
    double precision :: GaussPdsCoord
    dimension GaussPdsCoord(MCRD + 1, NBINT)

    !! Nurbs basis functions
    double precision :: R, dRdx, DetJac
    dimension R(NNODE), dRdx(MCRD, NNODE)

    !! Material behaviour
    double precision :: ddsdde
    dimension ddsdde(2 * MCRD, 2 * MCRD)

    !! Stiffness matrix
    integer :: k1, k2, ntens
    double precision :: stiffgeom, geom, dvol
    dimension stiffgeom( MCRD, MCRD, NNODE * (NNODE + 1) / 2 ), geom( MCRD, MCRD, NNODE * (NNODE + 1) / 2 )

    !! Initialization
    ntens   = 2 * MCRD          !! Size of stiffness tensor
    NbPtInt = int( NBINT ** (1.0 / float(MCRD)) ) !! Nb of gauss pts per direction
    if (NbPtInt ** MCRD < NBINT) NbPtInt = NbPtInt + 1

    !! Compute Gauss points coordinates and weights
    call Gauss(NbPtInt, MCRD, GaussPdsCoord, 0)

    AMATRX(:,:,:) = zero
    !! Material behaviour
    call material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde)

    !! Loop on integration points
    do n = 1, NBINT
        !! Compute NURBS basis functions and derivatives
        call shap(dRdx, R, DetJac, COORDS, GaussPdsCoord(2:, n), MCRD)

        !! Compute AMATRIX
        call stiffgeommatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, Uelem, stiffgeom)
        call geommatrix_byCP(ntens, NNODE, MCRD, NDOFEL, ddsdde, dRdx, Uelem, geom)

        ! Assemble AMATRIX
        dvol = GaussPdsCoord(1, n) * ABS(detJac)
        AMATRX(:,:,:) = AMATRX(:,:,:) + (stiffgeom(:,:,:) + geom(:,:,:)) * dvol
end do

end subroutine UELTANMAT_byCP

