!!     Calcul des deformations, contraintes et deplacements au niveau des
!!     frontiere delimitant un element nurbs. Cela permettra de creer le
!!     fichier VTU pour la visualisation.
subroutine compute_svars_Q1_HO_1stG(COORDS,sol,svars,nsvint,Output_FLAG, &
                &   nb_vertice,nb_REF,MCRD,NNODE,MATERIAL_PROPERTIES,TENSOR)
    use parameters

    implicit none

    !!  Input arguments :
    !!  ---------------
    double precision, intent(in) :: COORDS,sol,MATERIAL_PROPERTIES
    dimension COORDS(MCRD,NNODE)
    dimension sol(MCRD,NNODE)
    dimension MATERIAL_PROPERTIES(2)
    character(len=*), intent(in) :: TENSOR

    integer,intent(in) :: nsvint,MCRD,NNODE,nb_vertice,nb_REF
    dimension nb_REF(3)
      
    logical, intent(in) :: Output_FLAG
    dimension Output_FLAG(3)
      
    !! Output variables :
    !! ----------------
    double precision, intent(inout) :: svars
    dimension svars(nsvint*nb_vertice)
      
    !! Local variables :
    !! ---------------
    double precision :: vertice
    double precision :: R, dRdx, d2Rdx, DetJac
    dimension R(NNODE)
    dimension dRdx(MCRD,NNODE)
    dimension d2Rdx(3*(MCRD-1), NNODE)
    double precision :: stran, stress,    &
        &     svarsip,coords_ip,u_ip,ddsdde
    double precision :: dstrandx
    !! WARNING : only 2D is implemented for dstrandx
    dimension :: dstrandx(6)
    double precision :: dNidx,dNidy,dNidz
    double precision :: d2Nidx
    dimension d2Nidx(3*(MCRD-1))
    dimension vertice(MCRD,nb_vertice)

    dimension stran(2*MCRD),stress(2*MCRD),svarsip(nsvint),coords_ip(3),    &
        &     u_ip(3),ddsdde(2*MCRD,2*MCRD)

    integer :: NDOFEL,n,k1,ntens,nb_xi,nb_eta,nb_zeta,i_xi,     &
        &     i_eta,i_zeta,offset,i


    !! Initialization :
    !! --------------

    NDOFEL = NNODE*MCRD  

    !! Get material behaviour tensor
    !! TODO : ajuster pour avoir le comportement materiau d'un element HO
    call material_lib(MATERIAL_PROPERTIES,TENSOR,MCRD,ddsdde) 

    !! Defining element bounds : coords in parent space
    if (MCRD==2) then
        vertice(:,1) = (/-one, -one/)
        vertice(:,2) = (/ one, -one/)
        vertice(:,3) = (/ one,  one/)
        vertice(:,4) = (/-one,  one/)
    else if (MCRD==3) then
        vertice(:,1) = (/-one, -one, -one/)
        vertice(:,2) = (/-one,  one, -one/)
        vertice(:,3) = (/-one, -one,  one/)
        vertice(:,4) = (/-one,  one,  one/)
        vertice(:,5) = (/ one, -one, -one/)
        vertice(:,6) = (/ one,  one, -one/)
        vertice(:,7) = (/ one, -one,  one/)
        vertice(:,8) = (/ one,  one,  one/)
    endif

    !!  Compute disp., stress, strain :
    !! -----------------------------
    R(:) = zero
    dRdx(:,:) = zero
    d2Rdx(:,:) = zero
    svars(:) = zero

    nb_xi  = 2**max(nb_REF(1)-1,0)+1
    nb_eta = 2**max(nb_REF(2)-1,0)+1
    nb_zeta= 2**max(nb_REF(3)-1,0)+1
    if (MCRD==2) nb_zeta= 1
    do i_zeta= 1,nb_zeta
    do i_eta = 1,nb_eta
    do i_xi  = 1,nb_xi
        n = (i_zeta-1)*nb_eta*nb_xi + (i_eta-1)*nb_xi + i_xi
        vertice(1,n) = two/dble(nb_xi -1)*dble(i_xi -1) - one
        vertice(2,n) = two/dble(nb_eta-1)*dble(i_eta-1) - one
        if (MCRD==3) then
            vertice(3,n) = two/dble(nb_zeta-1)*dble(i_zeta-1) - one
        endif
        call shap_HO_1stG(R,dRdx,d2Rdx,DetJac,COORDS,vertice(1:,n),MCRD)

        !! Get intergration points coordinates and displacements
        coords_ip = zero
        u_ip = zero
        do k1 = 1,NNODE
            coords_ip(:MCRD) = coords_ip(:MCRD)+R(k1)*coords(:,k1)
            if (Output_FLAG(1)) then
                u_ip(:MCRD) = u_ip(:MCRD) + R(k1)*sol(:,k1)
            endif
        enddo

        !! Get strain
        ntens = 2*MCRD
        stran(:) = zero
        dstrandx(:) = zero
        if (Output_FLAG(2) .OR. Output_FLAG(3)) then
            do k1 = 1,NNODE
                dNidx = dRdx(1,k1)
                dNidy = dRdx(2,k1)
                d2Nidx(:) = d2Rdx(:,k1)
                stran(1) = stran(1) + dNidx*sol(1,k1)
                stran(2) = stran(2) + dNidy*sol(2,k1)
                stran(4) = stran(4) + dNidy*sol(1,k1)+dNidx*sol(2,k1)
                if (MCRD==3) then
                    dNidz = dRdx(3,k1)
                    stran(3) = stran(3) + dNidz*sol(3,k1)
                    stran(5) = stran(5) + dNidz*sol(1,k1)+dNidx*sol(3,k1)
                    stran(6) = stran(6) + dNidy*sol(3,k1)+dNidz*sol(2,k1)
                endif
                !! WARNING : only 2D is implemented for dstrandx
                dstrandx(1) = dstrandx(1) + d2Nidx(1)*sol(1,k1) 
                dstrandx(2) = dstrandx(2) + d2Nidx(2)*sol(1,k1)
                dstrandx(3) = dstrandx(3) + d2Nidx(3)*sol(1,k1)
                dstrandx(4) = dstrandx(4) + d2Nidx(1)*sol(2,k1)
                dstrandx(5) = dstrandx(5) + d2Nidx(2)*sol(2,k1)
                dstrandx(6) = dstrandx(6) + d2Nidx(3)*sol(2,k1)
            enddo
            !! Get stress
            if (Output_FLAG(2)) then
                call MulVect(ddsdde, stran, stress, ntens, ntens)
            endif
        endif

        !! Sum up all variables into svarsip
        svarsip = zero
        offset = 1

        svarsip(offset:offset+2) = coords_ip(:)
        offset = offset + 3

        if (Output_FLAG(1)) then
            svarsip(offset:offset+2) = u_ip(:)
            offset = offset + 3
        endif

        if (Output_FLAG(2)) then
            svarsip(offset:offset+ntens-1) = stress(:)
            offset = offset + ntens
        endif

        if (Output_FLAG(3)) then
            svarsip(offset:offset+ntens-1) = stran(:)
            offset = offset + ntens
            svarsip(offset:offset+6-1) = dstrandx(:)
        endif

        !! Update global variable : all variables at each intergration point
        !!write(*,*) n, " : ", svarsip(:)
        do i = 1,nsvint
            svars(nsvint*(n-1)+i) = svarsip(i)
        enddo
    enddo
    enddo
    enddo
end subroutine compute_svars_Q1_HO_1stG
