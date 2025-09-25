!! Copyright 2018-2019 Thibaut Hirschler
!! Copyright 2020-2023 Arnaud Duval
!! Copyright 2021 Marie Guerder

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

!! Build global sparse stiffness matrix for scipy.sparse
!! Inputs are handled by Python object IGAparametrization
!! Outputs : Kdata : vector containing element matrices
!!           Krow : matrix line indices
!!           Kcol : matrix column indices

!! New data model to avoid global variables that can cause problems with openMP

subroutine sys_linmat_lindef_static_omp(Kdata,Krow,Kcol,F,  &
    &   activeElement,nb_data, &
    &   COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight, &
    &   Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,n_mat_props,RHO,TENSOR,  &
    &   indDLoad,JDLType,ADLMAG,load_target_nbelem, &
    &   load_additionalInfos, nb_load_additionalInfos, bc_values,   &
    &   nb_bc,bc_target, &
    &   bc_target_nbelem,ind_dof_free,nb_dof_free,MCRD,NBINT,nb_load,   &
    &   nb_patch,nb_elem,nnode,nb_cp,nb_dof_tot,nodal_dist, &
    &   nb_n_dist, nb_cp_n_dist, n_mat_props_max,           &
    &   len_Ukv, len_weight, len_PROPS, len_IEN, len_indDLoad, &
    &   len_load_additionalInfos, len_nb_load_additionalInfos, &
    &   len_bc_target)

    use m_nurbspatch, only : nurbspatch, nurbselement
    use omp_lib


    implicit none

    !! Input arguments
    !! ---------------

    !! NURBS geometry
    integer, intent(in) :: nb_cp
    double precision, intent(in) :: COORDS3D
    dimension COORDS3D(3,nb_cp)

    integer, intent(in) :: len_Ukv, len_weight
    double precision, dimension(len_Ukv), intent(in) :: Ukv
    double precision, dimension(len_weight), intent(in) :: weight
    integer, intent(in) :: Nkv, Jpqr, Nijk
    dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem)


    !! Patches and Elements
    character(len=*), intent(in) :: TENSOR, ELT_TYPE
    double precision, intent(in) :: MATERIAL_PROPERTIES,RHO
    integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,   &
        &     nb_elem_patch, JPROPS
    integer, intent(in) :: n_mat_props      !! Number of material properties per patch
    dimension n_mat_props(nb_patch)
    integer, intent(in) :: n_mat_props_max  !! Maximum value of n_mat_props
    integer, intent(in) :: len_PROPS, len_IEN
    double precision, dimension(len_PROPS), intent(in) :: PROPS
    integer, dimension(len_IEN), intent(in) :: IEN
    dimension MATERIAL_PROPERTIES(n_mat_props_max,nb_patch),  &
        &     RHO(nb_patch),   &
        &     NNODE(nb_patch), &
        &     nb_elem_patch(nb_patch), &
        &     JPROPS(nb_patch),    &
        &     NBINT(nb_patch)


    !! Loads
    double precision, intent(in) :: ADLMAG, &
        &     nodal_dist
    integer, intent(in) :: nb_load,JDLType,load_target_nbelem
    integer, intent(in) :: nb_n_dist, nb_cp_n_dist
    integer, intent(in) :: len_indDLoad, len_nb_load_additionalInfos, len_load_additionalInfos
    integer, dimension(len_nb_load_additionalInfos), intent(in) :: nb_load_additionalInfos
    integer, dimension(len_indDLoad), intent(in) :: indDLoad
    double precision, dimension(len_load_additionalInfos), intent(in) :: load_additionalInfos
    dimension ADLMAG(nb_load),  &
        &     nodal_dist(nb_n_dist, nb_cp_n_dist),  &
        &     JDLType(nb_load), &
        &     load_target_nbelem(nb_load)


    !! Boundary Conditions
    double precision, intent(in) :: bc_values
    integer, intent(in) :: nb_bc,bc_target_nbelem
    dimension bc_values(2,nb_bc),   &
        &     bc_target_nbelem(nb_bc)

    integer, intent(in) :: len_bc_target
    integer, dimension(len_bc_target), intent(in) :: bc_target


    !! Degrees Of Freedom
    integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
    dimension ind_dof_free(nb_dof_tot)


    !! Storage infos
    integer(kind=8), intent(in) :: nb_data
    integer, dimension(nb_elem), intent(in) :: activeElement


    !! Output variables
    !! ----------------

    !! linear system to solve
    integer, dimension(nb_data), intent(out) :: Krow,Kcol
    double precision, dimension(nb_data), intent(out) :: Kdata
    double precision, dimension(nb_dof_tot), intent(out) :: F

    !! Local variables
    !! ---------------

    integer :: kk, ll, iload, icp, jcp, idof, jdof, i, ipatch, ielem
    type(nurbspatch) :: patch
    type(nurbselement) :: element

    double precision, dimension(maxval(n_mat_props)) :: mat_patch           !! material properties of a given patch
    double precision, dimension(mcrd,mcrd, maxval(nnode)*(maxval(nnode)+1)/2) :: amatrx         !! elementary stiffness matrix
    double precision, dimension(mcrd*maxval(nnode)) :: rhs                  !! elementary RHS vector
    double precision, dimension(mcrd, maxval(nnode)) :: coords_elem         !! CP coordinates of a given element
    double precision, dimension(nb_n_dist, maxval(nnode)) :: n_dist_elem    !! nodal distribution for a given element

    integer, dimension(maxval(nnode)) :: sctr

    integer :: ndofel, nnodeSum, nb_data_elem, count, jelem
    integer(kind=8) :: loccount, idx


    Kdata(:) = 0.0
    Krow(:) = 0.0
    Kcol(:) = 0.0
    F(:) = 0.0

    kk = 0
    do iload = 1,nb_load
        i = JDLType(iload)
        if (i/10 < 1) then
            do icp = 1,load_target_nbelem(iload)
                idof = (indDLoad(kk+icp)-1)*MCRD + i
                F(idof) = ADLMAG(iload)
            enddo
        endif
        kk = kk + load_target_nbelem(iload)
    enddo

    count = 1
    jelem = 0

    do ipatch = 1, nb_patch
        call patch%extractNurbsPatchGeoInfos(ipatch, Nkv, Jpqr, Nijk, Ukv, weight, nb_elem_patch)
        call patch%extractNurbsPatchMechInfos(ipatch, ien, props, jprops, nnode, nb_elem_patch, elt_type, tensor)

        if ((patch%elt_type_patch .eq. 'U30') .or. (patch%elt_type_patch .eq. 'U10')) then
            i = int(patch%props_patch(2))
            write(*,*) 'ERROR : THIS CASE IS NOT IMPLEMENTED YET'
            call exit(1)
            !! TODO : implement exreact mapping function for nurbspatch type
            ! call extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv, &
            ! &           weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)
        endif

        !! TODO : material should be defined as a derived type
        mat_patch(:n_mat_props(ipatch)) = material_properties(:n_mat_props(ipatch), ipatch)

        ndofel = patch%nnode_patch*mcrd
        nnodeSum = patch%nnode_patch*(patch%nnode_patch + 1)/2
        nb_data_elem = mcrd*mcrd*(patch%nnode_patch*(patch%nnode_patch-1))/2      &
            &           + patch%nnode_patch * (mcrd * (mcrd+1))/2
        write(*,*) 'nnode_patch : ', patch%nnode_patch

        !! Loop on elements
        !$omp parallel
        !$omp do private(element, coords_elem, n_dist_elem, rhs, amatrx), &
        !$omp& private(sctr, loccount, i, idx, icp, jcp, idof, jdof, ll, kk), &
        !$omp& reduction(+ : F)
        do ielem = 1, nb_elem_patch(ipatch)

            if(activeElement(jelem+ielem)==1) then
                call element%extractNurbsElementInfos(ielem, patch)

                do i = 1, patch%nnode_patch
                    coords_elem(:,i) = COORDS3D(:mcrd, patch%ien_patch(i, ielem))
                    n_dist_elem(:,i) = nodal_dist(:, patch%ien_patch(i, ielem))
                enddo

                rhs(:) = 0.0
                amatrx(:,:,:) = 0.0

                if (patch%elt_type_patch == 'U0') then
                    !! Void element --> do nothing
                elseif (patch%elt_type_patch == 'U1') then
                    !! Solid element
                    call uelmat_byCP_omp(patch, element, ndofel, mcrd, jelem+ielem, nbint(ipatch),   &
                                    &   coords_elem(:, :patch%nnode_patch),                &
                                    &   mat_patch(:2), rho(ipatch), nb_load,               &
                                    &   indDLoad, load_target_nbelem, JDLType, ADLMAG,     &
                                    &   load_additionalInfos, size(load_additionalInfos),  &
                                    &   nb_load_additionalInfos,                           &
                                    &   n_dist_elem, nb_n_dist,                            &
                                    &   rhs(:ndofel), amatrx(:,:,:nnodeSum))
                else
                    write(*,*) 'Element '//patch%elt_type_patch//' not available'
                    call exit(-1)
                endif

                !! Assemble amatrx to globall stiffness matrix K
                sctr(:patch%nnode_patch) = patch%ien_patch(:,ielem)
                loccount = 0
                i = 0
                !! start index for current element
                idx = count + (ielem-1) * nb_data_elem
                do jcp = 1, patch%nnode_patch
                    jdof= (sctr(jcp)-1)*MCRD
                    !! case cpi < cpj
                    do icp = 1,jcp-1
                        idof= (sctr(icp)-1)*MCRD
                        i   = i + 1
                        do ll = 1,MCRD
                            do kk = 1,MCRD
                                ! write(*,*) idx+loccount
                                Kdata(idx + loccount) = AMATRX(kk,ll,i)
                                Krow(idx + loccount) = idof + kk - 1
                                Kcol(idx + loccount) = jdof + ll - 1
                                loccount = loccount + 1
                            enddo
                        enddo
                    enddo

                    !! case cpi == cpj
                    i = i + 1
                    do ll = 1,MCRD
                        AMATRX(ll,ll,i) = AMATRX(ll,ll,i)*0.5d0
                        do kk = 1,ll
                            ! write(*,*) idx+loccount
                            Kdata(idx + loccount) = AMATRX(kk,ll,i)
                            Krow(idx + loccount) = jdof + kk - 1
                            Kcol(idx + loccount) = jdof + ll - 1
                            loccount = loccount + 1
                        enddo
                    enddo

                enddo

                !! Update Load Vector
                !! ********* ATTENTION, FB non traité *********
                idof = 0
                do i = 1, patch%nnode_patch
                    jdof  = (sctr(i)-1)*MCRD
                    do kk = 1,MCRD
                        F(jdof+kk) = F(jdof+kk) + RHS(idof+kk)
                    enddo
                    idof = idof + MCRD
                enddo
                !! ********* ATTENTION, FB non traité *********


            endif
        enddo
        !$omp end do
        !$omp end parallel

        jelem = jelem + nb_elem_patch(ipatch)
        count = count + nb_elem_patch(ipatch) * nb_data_elem

        !! TODO : reimplement these function for better memory management
        ! call deallocateMappingData()
        call patch%finalizeNurbsPatch()

    enddo

end subroutine sys_linmat_lindef_static_omp