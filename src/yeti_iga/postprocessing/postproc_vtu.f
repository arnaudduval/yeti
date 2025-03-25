!! Copyright 2018 Thibaut Hirschler
!! Copyright 2021 Arnaud Duval
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
!! with Yeti. if not, see <https://www.gnu.org/licenses/>

      !!include './compute_svars_solid.f'
      !!include './compute_svars_shell.f'
      !!include './compute_svars_embdedshell.f'

      subroutine generateVTU(FILENAME, output_path, FieldOutput_flag,
     1     nb_refinement,
     2     SOL, COORDS3D, IEN, nb_elem_patch, Nkv, Ukv, Nijk, weight,
     3     Jpqr, ELT_TYPE, MATERIAL_PROPERTIES, TENSOR, PROPS, JPROPS,
     4     NNODE, nb_patch, nb_elem, nb_cp, MCRD)

      use parameters
      use nurbspatch
      use embeddedMapping

      implicit none

      !! Input arguments
      !! ---------------

      !! NURBS Geometry
      integer, intent(in) :: nb_cp
      double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3, nb_cp)

      double precision, intent(in) :: Ukv, weight
      integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem),
     &     Ukv(:), weight(:)

      !! Patches and Elements
      character(len=*), intent(in) :: TENSOR, ELT_TYPE
      double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2, nb_patch), PROPS(:),
     &     JPROPS(nb_patch), NNODE(nb_patch), IEN(:),
     &     nb_elem_patch(nb_patch)

      !! Analysis solution
      double precision, intent(in) :: SOL
      dimension SOL(MCRD, nb_cp)

      !! Output infos
      character(len=*), intent(in) :: FILENAME
      character(len=*), intent(in) :: output_path
      logical, intent(in) :: FieldOutput_flag
      dimension FieldOutput_flag(3)
      integer, intent(in) :: nb_refinement
      dimension nb_refinement(3)

      !! Local variables
      !! ---------------
      double precision, dimension(:, :), allocatable :: svars

      double precision COORDS_elem, Sol_elem
      dimension COORDS_elem(MCRD, MAXVAL(NNODE))
      dimension sol_elem(MCRD, MAXVAL(NNODE))

      integer :: numel, ntens, nsvint, i_patch, JELEM, i, nb_vertice,
     &     nb_node, nb_elemVTU, nb_xi, nb_eta, nb_zeta, n

      ! test embedded mapping
      integer :: k, l, saveActiveElem, count, iter
      dimension saveActiveElem(nb_elem)
      double precision xi, R, dRdxi, COORDS_map, SOL_map
      dimension xi(3), R(MAXVAL(NNODE)), dRdxi(MAXVAL(NNODE), 3),
     &     COORDS_map(MCRD, MAXVAL(NNODE)), SOL_map(MCRD, MAXVAL(NNODE))

      !! Initialization
      !! --------------

      !! Test embedded stiff
      saveActiveElem(:) = 0

      !! Useful data when writing .vtu file
      nb_xi = 2**(max(nb_refinement(1)-1, 0)) + 1
      nb_eta = 2**(max(nb_refinement(2)-1, 0)) + 1
      nb_zeta = 2**(max(nb_refinement(3)-1, 0)) + 1
      if (MCRD==2) nb_zeta= 1

      !! Printing
    !   write(*, *)'Post processing ...'
    !   write(*, *)' Map data on FEM mesh ...'
    !   write(*, *)' Writing .vtu file ...'

      !! File
      Open(90, file=output_path // '/'// FILENAME //'.vtu',
     1    form='formatted')

      !! Write header
      write(90, *)'<VTKFile type="UnstructuredGrid"  version="0.1"   >'
      write(90, *)'<UnstructuredGrid>'

      !! Computation
      !! -----------

      JELEM = 0
      do i_patch = 1, nb_patch

          call extractNurbsPatchGeoInfos(i_patch, Nkv, Jpqr, Nijk, Ukv,
     &        weight, nb_elem_patch)
          call extractNurbsPatchMechInfos(i_patch, IEN, PROPS, JPROPS,
     &        NNODE, nb_elem_patch, ELT_TYPE, TENSOR)

          if (ELT_TYPE_patch == 'U30' .or. ELT_TYPE_patch == 'U10') then
              i = int(PROPS_patch(2))
              call extractMappingInfos(i, nb_elem_patch, Nkv, Jpqr,
     &            Nijk, Ukv, weight, IEN, PROPS, JPROPS, NNODE,
     &            ELT_TYPE, TENSOR)
          endif

          if (ELT_TYPE_patch == 'U1'
     &        .or. ELT_TYPE_patch == 'U2'
     &        .or. ELT_TYPE_patch == 'U3'
     &        .or. ELT_TYPE_patch == 'U30'
     &        .or. ELT_TYPE_patch == 'U10'
     &        ) then

              !! Constant defining the size of tensors
              ntens = 2*MCRD         ! =4 in 2D and =6 in 3D

              !! Number of variable of interest per integration point
              !! coord + U + ntens stress + ntens strain
              nsvint = 3
              if (FieldOutput_flag(1)) nsvint = nsvint + 3 ! DISP
              if (FieldOutput_flag(2)) nsvint = nsvint + ntens ! STRESS
              if (FieldOutput_flag(3)) nsvint = nsvint + ntens ! STRAIN

              !! Svars contains all values computed on an element
              nb_elemVTU = 1
              nb_vertice = 1
              do i = 1, dim_patch
                  nb_elemVTU = nb_elemVTU *
     &                2**(max(nb_refinement(i)-1, 0))
                  nb_vertice = nb_vertice *
     &                (2**(max(nb_refinement(i)-1, 0))+1)
              enddo
              nb_node = nb_elem_patch(i_patch) * nb_vertice
              allocate(svars(nsvint*nb_vertice, nb_elem_patch(i_patch)))

              !! Loop on element
              svars = 0.0D0

              count = 1
              do iter = 1, count

                  do numel = 1, nb_elem_patch(i_patch)
                      JELEM = JELEM +1

                      !! Extract element solution
                      do i = 1,nnode_patch
                          COORDS_elem(:, i) =
     &                        COORDS3D(:MCRD, IEN_patch(i, numel))
                          sol_elem(:, i) =
     &                        sol(:MCRD, IEN_patch(i, numel))
                      enddo
                      call extractNurbsElementInfos(numel)

c                     Calcul des variables d'interet aux frontieres de l'element
                      if (ELT_TYPE_patch == 'U1') then
                          call compute_svars_Q1(COORDS_elem, sol_elem,
     1                        svars(:, numel), nsvint, FieldOutput_flag,
     2                        nb_vertice, nb_refinement, MCRD,
     3                        NNODE(i_patch),
     4                        MATERIAL_PROPERTIES(:, i_patch),
     5                        TENSOR_patch)

                      elseif (ELT_TYPE_patch == 'U10') then
                          call compute_svars_Q1_embdedsolid(COORDS_elem,
     1                        COORDS3D, sol_elem, svars(:, numel),
     2                        nsvint, FieldOutput_flag, nb_vertice,
     3                        nb_refinement, MCRD, NNODE(i_patch),
     4                        nnode_map, nb_cp,
     5                        MATERIAL_PROPERTIES(:, i_patch),
     6                        TENSOR_patch)
c                     elseif (ELT_TYPE == 'U2') then
c                         call compute_svars_Q1_plate(COORDS_elem,sol_elem,u_elem,
c     1                       svars,nsvint,nb_vertice,MCRD,NNODE,JELEM,
c     2                       i_patch,MATERIAL_PROPERTIES(i_patch,:),TENSOR)

                      elseif (ELT_TYPE_patch=='U2' .OR.
     &                        ELT_TYPE_patch=='U3') then
                          call compute_svars_Q1_shell(COORDS_elem,
     1                        sol_elem, svars(:, numel), nsvint,
     2                        FieldOutput_flag, nb_vertice,
     3                        nb_refinement, MCRD, NNODE(i_patch),
     4                        MATERIAL_PROPERTIES(:, i_patch),
     5                        TENSOR_patch, PROPS_patch, JPROPS_patch)

                      elseif (ELT_TYPE_patch=='U30') then
                          call compute_svars_Q1_embdedshell(COORDS_elem,
     1                        COORDS3D, sol_elem, svars(:, numel),
     2                        nsvint, FieldOutput_flag, nb_vertice,
     3                        nb_refinement, MCRD, NNODE(i_patch),
     4                        nnode_map, nb_cp,
     5                        MATERIAL_PROPERTIES(:, i_patch),
     6                        TENSOR_patch, PROPS_patch, JPROPS_patch)

                      elseif (ELT_TYPE_patch == 'U31') then
                          ! Embedded surface
                          n = NNODE(i_patch)
                          call compute_svars_Q1_shell(
     1                        COORDS_elem(:, :n), sol_elem(:, :n),
     2                        svars(:, numel), nsvint,
     3                        (/.False., .False., .False./),
     4                        nb_vertice, nb_refinement, MCRD, n,
     5                        MATERIAL_PROPERTIES(:, i_patch),
     6                        TENSOR_patch, PROPS_patch, JPROPS_patch)
                          ! Mapping
                          if (iter == 1) then
                              l = 0
                              do n = 1, nb_vertice
                                  xi(:) = svars(l+1: l+3, numel)
                                  call updateMapElementNumber(xi)
                                  saveActiveElem(current_map_elem) = 1
                                  !print*,'current_map_elem',current_map_elem
                                  do i = 1, nnode_map
                                      k = IEN_map(i, current_map_elem)
                                      COORDS_map(:, i) =
     &                                    COORDS3D(:MCRD, k)
                                      SOL_map(:, i) = SOL(:MCRD, k)
                                  enddo

                                  call evalnurbs_mapping(xi,
     &                                R(:nnode_map),
     &                                dRdxi(:nnode_map, :))

                                  svars(l+1: l+3, numel) = zero
                                  do i = 1, nnode_map
                                      svars(l+1: l+3, numel) =
     &                                    svars(l+1: l+3, numel)
     &                                    + R(i)*COORDS_map(:, i)
                                      svars(l+4: l+6, numel) =
     &                                    svars(l+4: l+6, numel)
     &                                    + R(i)*SOL_map(:, i)
                                      svars(l+7: l+nsvint, numel) = zero
                                  enddo
                                  !print*,'Rm',MINVAL(R),'RM',MAXVAL(R)
                                  l = l + nsvint
                              enddo
                          endif

                      else
                          print*, ' /!\ Element type'// ELT_TYPE_patch
     &                         // 'not avalaible for postprocessing.'
                      endif
                  enddo

                  ! Writting patch results
                  call writepatch(90, FieldOutput_flag, svars, nb_node,
     &                nb_elem_patch(i_patch), dim_patch, MCRD, nsvint,
     &                ntens, nb_elemVTU, nb_xi, nb_eta, nb_zeta,
     &                nb_vertice)

              enddo

              deallocate(svars)

          endif

          call finalizeNurbsPatch()
          call deallocateMappingData()

      enddo

c     test embedded stiff
      !print*,'Active Elements'
      !do i = 1,nb_elem_patch(current_map)
      !   if (saveActiveElem(i)>0) print*,i-1,','
      !enddo

      !! Writing End of file
      write(90, *)'</UnstructuredGrid>'
      write(90, *)'</VTKFile>'
      close(90)

    !   write(*, *)' --> File ' // output_path // '/' // FILENAME //
    !  1     '.vtu has been created.'

      end subroutine generateVTU
