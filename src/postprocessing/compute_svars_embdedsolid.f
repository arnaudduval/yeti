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
!!
!! Compute displacement, strain and stress at the boundaries of a NURBS
!! embedded solid element. Used to further fill in a .vtu file for 
!! visualisation.
!!
      subroutine compute_svars_Q1_embdedsolid(COORDS, COORDSall, sol,
     &     svars, nsvint, Output_FLAG, nb_vertice, nb_REF, MCRD, NNODE,
     &     NNODEmap, nb_cp, MATERIAL_PROPERTIES, TENSOR)
      
      use parameters
      use embeddedMapping
      
      implicit none
           
      !! Input arguments
      !! ---------------
      double precision, intent(in) :: COORDS, COORDSall, sol, 
     &     MATERIAL_PROPERTIES
      dimension COORDS(MCRD, NNODE)
      dimension COORDSall(3, nb_cp)
      dimension sol(MCRD, NNODE)
      dimension MATERIAL_PROPERTIES(2)
      character(len=*), intent(in) :: TENSOR
      
      integer,intent(in) :: nsvint, MCRD, NNODE, NNODEmap, nb_cp, 
     &     nb_vertice, nb_REF
      dimension nb_REF(3)
      
      logical, intent(in) :: Output_FLAG
      dimension Output_FLAG(3)
      
      !! Output variables
      !! ----------------
      double precision, intent(inout) :: svars
      dimension svars(nsvint*nb_vertice)
      
      !! Local variables
      !! ---------------
      !! For gauss points
      double precision :: vertice
      dimension vertice(MCRD, nb_vertice)
      
      !! Embedded Volume
      ! - NURBS basis functions
      double precision :: R, dRdTheta, dRdxi
      dimension R(NNODE), dRdTheta(NNODE, 3), dRdxi(3, NNODE)
      ! - Parametric space
      double precision :: Theta
      dimension Theta(3)
      ! - Physical space
      double precision :: xi
      dimension xi(3)
      ! - From parametric to physical space
      double precision :: dxidTheta, dThetadxi, detdxidTheta
      dimension dxidTheta(3, 3), dThetadxi(3, 3)
      ! - From parent element to pametric space
      double precision :: dtildexidTheta, dThetadtildexi 
      dimension dtildexidTheta(3, 3), dThetadtildexi(3, 3)
      
      !! Mapping
      ! - NURBS basis functions
      double precision :: N, dNdxi
      dimension N(NNODEmap), dNdxi(NNODEmap, 3)
      ! - From parametric to physical space
      double precision :: dXdxi, dxidX, detdXdxi, dThetadX
      dimension dXdxi(3, 3), dxidX(3, 3), dThetadX(3, 3)
      ! - From parent element to pametric space
      double precision :: dxidtildexi
      dimension dxidtildexi(3, 3)
      ! - Element infos
      double precision :: COORDSmap
      dimension COORDSmap(MCRD, NNODEmap)
      integer :: sctr_map
      dimension sctr_map(NNODEmap)

      !! Composition Mapping+Volume
      double precision :: dRdX, dRdXT
      dimension dRdX(3, NNODE), dRdXT(NNODE, 3)
      
      !! For material matrix
      double precision :: ddsdde
      dimension ddsdde(2*MCRD, 2*MCRD)
      
      !! Other
      double precision :: stran, stress, coef,
     &     svarsip, coords_ip, u_ip, dNidx, dNidy, dNidz
      dimension stran(2*MCRD), stress(2*MCRD), svarsip(nsvint), 
     &     coords_ip(3), u_ip(3)
      
      integer :: NDOFEL, n_p, k1, ntens, nb_xi, nb_eta, nb_zeta, i_xi,
     &     i_eta, i_zeta, offset, i, isave, numCP, j
      
      
      !! Initialization 
      !! --------------
      
      NDOFEL = NNODE*MCRD  
      ntens = 2*MCRD
      
      !! Get material behaviour tensor
      call material_lib(MATERIAL_PROPERTIES, TENSOR, MCRD, ddsdde) 
      
      !! Defining element bounds : coords in parent space
      if (MCRD==2) then
          vertice(:, 1) = (/-one, -one/)
          vertice(:, 2) = (/ one, -one/)
          vertice(:, 3) = (/ one,  one/)
          vertice(:, 4) = (/-one,  one/)
      else if (MCRD==3) then
          vertice(:, 1) = (/-one, -one, -one/)
          vertice(:, 2) = (/-one,  one, -one/)
          vertice(:, 3) = (/-one, -one,  one/)
          vertice(:, 4) = (/-one,  one,  one/)
          vertice(:, 5) = (/ one, -one, -one/)
          vertice(:, 6) = (/ one,  one, -one/)
          vertice(:, 7) = (/ one, -one,  one/)
          vertice(:, 8) = (/ one,  one,  one/)
      endif
      
      !! Compute disp., stress, strain
      !!-----------------------------
      
      isave = 0
      
      svars = zero
      nb_xi  = 2**max(nb_REF(1)-1, 0) + 1
      nb_eta = 2**max(nb_REF(2)-1, 0) + 1
      nb_zeta= 2**max(nb_REF(3)-1, 0) + 1
      if (MCRD==2) nb_zeta= 1
      
      do i_zeta = 1, nb_zeta
          do i_eta = 1, nb_eta
              do i_xi = 1, nb_xi
                  n_p = (i_zeta-1)*nb_eta*nb_xi + (i_eta-1)*nb_xi + i_xi
                  vertice(1, n_p) = two/dble(nb_xi -1)*dble(i_xi -1) -
     &                one
                  vertice(2, n_p) = two/dble(nb_eta-1)*dble(i_eta-1) - 
     &                one
                  if (MCRD==3) then
                      vertice(3, n_p) = 
     &                    two/dble(nb_zeta-1)*dble(i_zeta-1) - one
                  endif

                  !! Compute parametric coordinates from parent element
                  Theta(:) = zero
                  do j = 1, 3
                      coef = vertice(j, n_p)
                      Theta(j) = ((Ukv_elem(2, j) - Ukv_elem(1, j))*coef
     &                    + (Ukv_elem(2, j) + Ukv_elem(1, j)))*0.5d0
                  enddo
                  
                  !! Compute NURBS basis functions and derivatives of the
                  !! embedded entity (in parametric space)                 
                  call evalnurbs(Theta, R, dRdTheta) 
                  
                  !! Compute mapping parametric position
                  xi(:) = zero
                  do numCP = 1, NNODE
                      xi(:) =  xi(:) + R(numCP)*COORDS(:, numCP)
                  enddo                  
                  
                  !! Computing NURBS basis functions and derivatives of the
                  !! mapping
                  
                  ! Get active element number
                  call updateMapElementNumber(xi(:))       
                  
                  ! Evaluate functions and derivatives
                  call evalnurbs_mapping(xi(:), N(:), dNdxi(:, :))
                  
                  !! Get mapping coordinates                  
                  
                  ! Extract coordinates of the CPs of the mapping
                  if (isave /= current_map_elem) then
                      sctr_map(:) = IEN_map(:, current_map_elem)
                      do numCP = 1, NNODEmap 
                          COORDSmap(:, numCP) = 
     &                             COORDSall(:, sctr_map(numCP))
                      enddo
                      isave = current_map_elem
                  endif
                  
                  !! Get integration points coordinates and displacements
                  coords_ip = zero
                  do k1 = 1, NNODEmap
                      coords_ip(:MCRD) = coords_ip(:MCRD) + 
     &                    N(k1)*COORDSmap(:, k1)
                  enddo                  
                  
                  u_ip = zero
                  if (Output_FLAG(1)) then
                      do k1 = 1, NNODE   
                          u_ip(:MCRD) = u_ip(:MCRD) + R(k1)*sol(:, k1)
                      enddo
                  endif

                  !! Stress and strain processing
                  if (Output_FLAG(2) .OR. Output_FLAG(3)) then
                      
                      !! Composition of the basis funs. for stress and strain                      
                      
                      ! Gradient of mapping from parent element to parameter space
                      dThetadtildexi(:, :) = zero
                      do j = 1, dim_patch
                          dThetadtildexi(j, j) = 0.5d0*(Ukv_elem(2, j) -
     &                        Ukv_elem(1, j))
                      enddo
                  
                      !! Gradient of mapping from parameter space to physical space
                            
                      ! Gradient from embedded volume param. space to mapping param. space
                      dxidTheta(:, :) = zero
                      do numCP = 1, NNODE
                          dxidTheta(:, 1) = dxidTheta(:, 1) + 
     &                        dRdTheta(numCP, 1)*COORDS(:, numCP)
                          dxidTheta(:, 2) = dxidTheta(:, 2) + 
     &                        dRdTheta(numCP, 2)*COORDS(:, numCP)
                          dxidTheta(:, 3) = dxidTheta(:, 3) + 
     &                        dRdTheta(numCP, 3)*COORDS(:, numCP)
                      enddo
          
                      call MatrixInv(dThetadxi, dxidTheta, 
     &                    detdxidTheta, 3)
          
                      ! Gradient from mapping param. space to mapping phys. space                 
                      dXdxi(:, :) = zero
                      do numCP = 1, NNODEmap
                          dXdxi(:, 1) = dXdxi(:, 1) + 
     &                        dNdxi(numCP, 1)*COORDSmap(:, numCP)
                          dXdxi(:, 2) = dXdxi(:, 2) + 
     &                        dNdxi(numCP, 2)*COORDSmap(:, numCP)
                          dXdxi(:, 3) = dXdxi(:, 3) + 
     &                        dNdxi(numCP, 3)*COORDSmap(:, numCP)
                      enddo
          
                      call MatrixInv(dxidX, dXdxi, detdXdxi, 3)
          
                      ! Intermediate product
                      call MulMat(dThetadxi, dxidX, dThetadX, 3, 3, 3)
          
                      ! Basis functions composition
                      call Mulmat(dRdTheta, dThetadX, dRdXT, 
     &                    NNODE, 3, 3)
          
                      dRdX(:, :) = zero
                      do numCP = 1, NNODE
                          dRdX(1, numCP) = dRdXT(numCP, 1)
                          dRdX(2, numCP) = dRdXT(numCP, 2)
                          dRdX(3, numCP) = dRdXT(numCP, 3)
                      enddo

                      !! Compute strain
                      stran = zero
                      do k1 = 1,NNODE
                          dNidx = dRdX(1, k1)
                          dNidy = dRdX(2, k1)
                          stran(1) = stran(1) + dNidx*sol(1, k1)
                          stran(2) = stran(2) + dNidy*sol(2, k1)
                          stran(4) = stran(4) + dNidy*sol(1, k1) + 
     &                        dNidx*sol(2, k1)
                          if (MCRD==3) then
                              dNidz = dRdX(3,k1)
                              stran(3) = stran(3) + dNidz*sol(3, k1)
                              stran(5) = stran(5) + dNidz*sol(1, k1) + 
     &                            dNidx*sol(3, k1)
                              stran(6) = stran(6) + dNidy*sol(3, k1) + 
     &                            dNidz*sol(2, k1)
                          endif
                      enddo
                      
                      !! Compute stress
                      if (Output_FLAG(2)) then
                          call MulVect(ddsdde, stran, stress, 
     &                        ntens, ntens)
                      endif
                  endif
                  
                  !! Sum up all variables into svarsip
                  svarsip = zero
                  offset = 1
                  !! - Coordinates
                  svarsip(offset:offset+2) = coords_ip(:)
                  offset = offset + 3
                  !! - Displacement
                  if (Output_FLAG(1)) then
                      svarsip(offset:offset+2) = u_ip(:)
                      offset = offset + 3
                  endif
                  !! - Stress
                  if (Output_FLAG(2)) then
                      svarsip(offset:offset+ntens-1) = stress(:)
                      offset = offset + ntens
                  endif
                  !! - Strains
                  if (Output_FLAG(3)) then
                      svarsip(offset:offset+ntens-1) = stran(:)
                  endif
                  
                  !! Update global variable : all variables at each 
                  !! integration point
                  do i = 1, nsvint
                      svars(nsvint*(n_p-1)+i) = svarsip(i)
                  enddo
                  
              enddo
          enddo
      enddo
      
C     ------------------------------------------------------------------
      
      end subroutine compute_svars_Q1_embdedsolid