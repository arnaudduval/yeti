! subroutine uelmat_bycp_omp(plop)
!     implicit none

!     integer :: plop

!     write(*,*) plop

! end subroutine uelmat_bycp_omp


subroutine uelmat_byCP_omp(patch, element, ndofel, mcrd, jelem, nbint, coords,            &
    &   material_properties, density, nb_load, indDLoad,    &
    &   load_target_nbelem,JDLType,ADLMAG,load_additionalInfos, &
    &   len_load_additionalInfos, nb_load_additionalInfos, &
    &   n_dist_elem, nb_n_dist, RHS,AMATRX)

    use m_nurbspatch, only : nurbspatch, nurbselement

    implicit none

    !! Input arguments
    !! ---------------

    type(nurbspatch), intent(in) :: patch
    type(nurbselement), intent(in) :: element

    integer, intent(in) :: NDOFEL,MCRD,JELEM,NBINT
    double precision, intent(in) :: COORDS,MATERIAL_PROPERTIES, &
        &   DENSITY,    &
        &   n_dist_elem
    dimension COORDS(MCRD,patch%nnode_patch),MATERIAL_PROPERTIES(2),    &
        &     n_dist_elem(nb_n_dist,patch%nnode_patch)

    integer, intent(in) :: indDLoad,load_target_nbelem,JDLType, &
        &     nb_load,len_load_additionalInfos, nb_n_dist,  &
        &     nb_load_additionalInfos

    double precision, intent(in) :: ADLMAG,load_additionalInfos
    dimension ADLMAG(nb_load),indDLoad(SUM(load_target_nbelem)),    &
        &     load_target_nbelem(nb_load),JDLType(nb_load),         &
        &     load_additionalInfos(len_load_additionalInfos),       &
        &     nb_load_additionalInfos(nb_load)

    !! Output variables
    !! ----------------
    double precision, intent(out) :: RHS, AMATRX
    dimension RHS(NDOFEL), AMATRX(MCRD,MCRD,patch%nnode_patch*(patch%nnode_patch+1)/2)


end subroutine uelmat_byCP_omp