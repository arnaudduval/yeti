!! Copyright 2019 Thibaut Hirschler

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


      !! Evaluate displacement at give parametric coordinates from solution given
      !! et control points

      subroutine evalDisp(disp,NumPatch,Xi,
     1     SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,
     3     nb_patch,nb_elem,nb_cp,MCRD)

      use parameters
      use nurbspatch

      Implicit none

C     ------------------------------------------------------------------

C     Input arguments :
c     ---------------
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)

      Double precision, intent(in) :: Ukv, weight
      Integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)

!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),PROPS(:),
     &     JPROPS(nb_patch),NNODE(nb_patch),IEN(:),
     &     nb_elem_patch(nb_patch)

!     Analysis solution
      Double precision, intent(in) :: SOL
      dimension SOL(MCRD,nb_cp)

!     Infos
      Integer,          intent(in) :: NumPatch
      Double precision, intent(in) :: Xi
      dimension Xi(3)


C     Output variables :
c     ----------------
      Double precision, intent(out) :: disp
      dimension disp(3)


C     Local variables :
c     ---------------
      Integer :: sctr,i
      dimension sctr(NNODE(numPatch))
      Double precision :: R,dRdxi
      dimension R(NNODE(numPatch)), dRdxi(NNODE(numPatch),3)


C     ------------------------------------------------------------------


c     Retour ecran
      !write(*,*)'Evaluate displacement '

      CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &     weight,nb_elem_patch)
      CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &     NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

      CALL updateElementNumber(Xi)



      call evalnurbs(Xi,R,dRdxi)
      sctr(:) = IEN_patch(:,current_elem)
      disp(:) = zero
      Do i = 1,NNODE(numPatch)
         disp(:MCRD) = disp(:MCRD) + R(i)*sol(:,sctr(i))
      Enddo

      CALL finalizeNurbsPatch()

C     ------------------------------------------------------------------

      end subroutine evalDisp

      subroutine evalDispMulti(disp,NumPatch,n_xi, n_eta, n_zeta,
     1     SOL,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,MATERIAL_PROPERTIES,TENSOR,PROPS,JPROPS,NNODE,
     3     nb_patch,nb_elem,nb_cp,MCRD)

      use parameters
      use nurbspatch

      implicit none

C     Input arguments
c     ---------------

!     NURBS geometry
      integer, intent(in) :: nb_cp
      double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3, nb_cp)

      double precision, intent(in) :: Ukv, weight
      integer, intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3, nb_patch), Jpqr(3, nb_patch), Nijk(3, nb_elem),
     &      Ukv(:), weight(:)

!     Patches and elements
      character(len=*), intent(in) :: TENSOR, ELT_TYPE
      double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      integer, intent(in) :: MCRD, NNODE, nb_patch, nb_elem, IEN,
     &      nb_elem_patch, JPROPS
      dimension MATERIAL_PROPERTIES(2, nb_patch), PROPS(:),
     &      JPROPS(nb_patch), NNODE(nb_patch), IEN(:),
     &      nb_elem_patch(nb_patch)

!     Analysis solution
      double precision, intent(in) :: SOL
      dimension SOL(MCRD, nb_cp)

!     Infos
      integer, intent(in) :: NumPatch
      integer, intent(in) :: n_xi, n_eta, n_zeta


!     Output variables
!     ----------------
      double precision, intent(out) :: disp
      dimension disp(n_xi, n_eta, n_zeta, 3)

!     Local variables
!     ---------------
      integer :: sctr, i, ixi, ieta, izeta
      dimension sctr(NNODE(numPatch))
      double precision :: R, dRdxi, Xi
      dimension R(NNODE(numPatch)), dRdxi(NNODE(numPatch),3)
      dimension Xi(3)


      CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &     weight,nb_elem_patch)
      CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &     NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

      do ixi = 1, n_xi
        do ieta = 1, n_eta
          do izeta = 1, n_zeta
!!          Assertion : knot vector are between 0 and 1
            Xi(1) = (ixi-1.D0)*1.D0/(n_xi-1.D0)
            Xi(2) = (ieta-1.D0)*1.D0/(n_eta-1.D0)
            Xi(3) = (izeta-1.D0)*1.D0/(n_zeta-1.D0)

            CALL updateElementNumber(Xi)

            call evalnurbs(Xi,R,dRdxi)
            sctr(:) = IEN_patch(:,current_elem)
            disp(ixi, ieta, izeta, :) = zero
            Do i = 1,NNODE(numPatch)
              disp(ixi, ieta, izeta, :MCRD) =
     &             disp(ixi, ieta, izeta, :MCRD)
     &             + R(i)*sol(:,sctr(i))
            Enddo
          enddo
        enddo
      enddo
      CALL finalizeNurbsPatch()
      end subroutine evalDispMulti

