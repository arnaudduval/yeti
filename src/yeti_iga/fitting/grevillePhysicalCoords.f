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

c     ...
c     Construction des  matrices pour l'interpolation de B-Spline avec
c     comme coordonnees parametriques les abscisses de greville.
c      --> Renvoie les matrices au format creux CSR
c     ...
      Subroutine grevilleAbscissae(Uknot,p,n,nknots,greAbsc)

        Implicit none

        ! Inputs
        Integer,          intent(in) :: p,n,nknots
        Double precision, intent(in) :: Uknot
        dimension Uknot(nknots)

        ! Output
        Double precision, intent(out):: greAbsc
        dimension greAbsc(n)

        ! Local var
        Integer :: i,k

    !     --
        greAbsc(:) = 0.0d0
        Do i = 1,n
        Do k = 1,p
            greAbsc(i) = greAbsc(i) + Uknot(i+k)
        Enddo
        Enddo
        greAbsc(:) = greAbsc(:)/dble(p)

      End subroutine grevilleAbscissae






      subroutine getGrevAbscPhysicalCoords(COORDSgrev,nb_grev,
     1     activePatch,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,TENSOR,PROPS,JPROPS,NNODE,nb_patch,nb_elem,
     3     nb_cp)

      use parameters
      use nurbspatch
      use embeddedMapping

      Implicit none

C     ------------------------------------------------------------------

C     Input arguments :
c     ---------------
!     Geometry NURBS
      Integer, intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)

      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)

!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: PROPS
      Integer,          intent(in) :: NNODE,nb_patch,nb_elem,IEN,
     &     nb_elem_patch, JPROPS
      dimension PROPS(:),
     &     JPROPS(nb_patch),NNODE(nb_patch),IEN(:),
     &     nb_elem_patch(nb_patch)

!     Additional infos
      Integer,intent(in) :: nb_grev,activePatch
      dimension activePatch(nb_patch)

!     Output arguments :
!     ----------------
      Double precision, intent(out):: COORDSgrev
      dimension COORDSgrev(3,nb_grev)


C     Local variables :
c     ---------------

!     Nurbs basis functions
      Integer          :: n,p,ng,mg,ngtot,direction,IND,cpbyDir
      dimension IND(nb_cp),cpbyDir(3)
      Double precision :: xi_bar,Ukv_dir,XI,R,dRdxi
      dimension xi_bar(MAXVAL(Nkv)),Ukv_dir(MAXVAL(Nkv)),XI(3),
     &     R(MAXVAL(NNODE)),dRdxi(MAXVAL(NNODE),3)

      Integer          :: sctr,sctr_map
      dimension sctr(MAXVAL(NNODE)),sctr_map(MAXVAL(NNODE))
      Double precision :: COORDS_elem,COORDSmap,VECT
      dimension COORDS_elem(3,MAXVAL(NNODE)),COORDSmap(3,nb_cp),VECT(3)

      Integer :: i,j,k,numPatch,count


C     ------------------------------------------------------------------


      COORDSgrev(:,:) = zero
      count = 0

      DO numPatch = 1,nb_patch

         IF (activePatch(numPatch) == 1) then

         ! extract infos
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

         If (ELT_TYPE_patch == 'U30') then
         i = int(PROPS_patch(2))
         CALL extractMappingInfos(i,nb_elem_patch,Nkv,Jpqr,Nijk,Ukv,
     &        weight,IEN,PROPS,JPROPS,NNODE,ELT_TYPE,TENSOR)
         Endif

!     1. Get Greville abscissae

         cpbyDir(:) = 1
         cpbyDir(:dim_patch) =  Nkv_patch(:dim_patch)
     &                       -(Jpqr_patch(:dim_patch)+1)
         ngtot = cpbyDir(1)*cpbyDir(2)*cpbyDir(3)
         Do direction = 1,dim_patch
            p = Jpqr_patch(direction)
            n =  Nkv_patch(direction)
            ng= cpbyDir(direction)
            If     (direction == 1) then
               Ukv_dir(:n) = Ukv1_patch(:)

               mg = cpbyDir(2)*cpbyDir(3)
               IND(:ngtot) = (/( (/(i,i=1,ng)/), j=1,mg )/)

            elseif (direction == 2) then
               Ukv_dir(:n) = Ukv2_patch(:)

               k = 0
               Do j = 1,cpbyDir(3)
                  Do i = 1,cpbyDir(2)
                     IND(k+1:k+cpbyDir(1)) = i
                     k = k + cpbyDir(1)
                  Enddo
               Enddo

            else
               Ukv_dir(:n) = Ukv3_patch(:)

               mg = cpbyDir(1)*cpbyDir(2)
               IND(:ngtot) = (/( (/(0,i=1,mg)/)+j, j=1,cpbyDir(3) )/)

            Endif
            xi_bar(:) = zero
            call grevilleAbscissae(Ukv_dir(:n),p,ng,n,xi_bar(:ng))

            Do k = 1,ngtot
               COORDSgrev(direction,count+k) = xi_bar(IND(k))
            Enddo
         Enddo ! loop on parametric direction


!     2. Compute Physical position
         Do k = 1,ngtot
            XI(:) = COORDSgrev(:,count+k)

            call updateElementNumber(XI)
            sctr(:nnode_patch) = IEN_patch(:,current_elem)
            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:,sctr(i))
            Enddo

            call evalnurbs(XI(:),R(:nnode_patch),dRdxi(:nnode_patch,:))
            VECT(:) = zero
            Do i = 1,nnode_patch
               VECT(:) = VECT(:) + R(i)*COORDS_elem(:,i)
            Enddo

            If (ELT_TYPE_patch == 'U30') then
               call updateMapElementNumber(VECT(:))
               sctr_map(:nnode_map) = IEN_map(:,current_map_elem)
               Do i = 1,nnode_map
                  COORDSmap(:,i) = COORDS3D(:,sctr_map(i))
               Enddo

               call evalnurbs_mapping(VECT(:),R(:nnode_map),
     &              dRdxi(:nnode_map,:))
               VECT(:) = zero
               Do i = 1,nnode_map
                  VECT(:) = VECT(:) + R(i)*COORDSmap(:,i)
               Enddo
            Endif

            COORDSgrev(:,count+k) = VECT(:)

         Enddo

         count = count + ngtot

         CALL finalizeNurbsPatch()
         CALL deallocateMappingData()

         ENDIF ! active patch
      ENDDO ! loop on patch


C     ------------------------------------------------------------------

      End subroutine getGrevAbscPhysicalCoords

