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


    !   include "./grevilleAbscissae.f"

      subroutine buildGrevInterpolMat(Mdata,Mindices,Mindptr,
     1     nb_data,nb_row, activePatch,
     2     IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,TENSOR,
     3     PROPS,JPROPS,NNODE,nb_patch,nb_elem)

      use parameters
      use nurbspatch

      Implicit none

C     ------------------------------------------------------------------

C     Input arguments :
c     ---------------
!     Geometry NURBS
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
      Integer,intent(in) :: nb_data,nb_row,activePatch
      dimension activePatch(nb_patch)

!     Output arguments :
!     ----------------
      Double precision, intent(out):: Mdata
      Integer,          intent(out):: Mindices,Mindptr
      dimension Mdata(nb_data),Mindices(nb_data),Mindptr(nb_row)


C     Local variables :
c     ---------------

!     Nurbs basis functions
      Integer          :: n,p,ng,direction
      Double precision :: xi_bar,Ukv_dir,FN,dNdxi
      dimension xi_bar(MAXVAL(Nkv)),Ukv_dir(MAXVAL(Nkv)),
     &     FN(MAXVAL(Jpqr)+1),dNdxi(MAXVAL(Jpqr)+1)

      Integer :: i,j,k,numPatch,rcount,dcount


C     ------------------------------------------------------------------


      Mdata(:)    = zero
      Mindices(:) = 0
      Mindptr(:)  = 0

      rcount = 2
      dcount = 1

      DO numPatch = 1,nb_patch

         IF (activePatch(numPatch) == 1) then

         ! extract infos
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

         Do direction = 1,dim_patch

!     Get Greville abscissae
            p = Jpqr_patch(direction)
            n =  Nkv_patch(direction)
            ng= n - (p+1)
            If     (direction == 1) then
               Ukv_dir(:n) = Ukv1_patch(:)
            elseif (direction == 2) then
               Ukv_dir(:n) = Ukv2_patch(:)
            else
               Ukv_dir(:n) = Ukv3_patch(:)
            Endif
            xi_bar(:) = zero
            call grevilleAbscissae(Ukv_dir(:n),p,n,ng,xi_bar(:ng))


!     Evaluate bspline basis functions
            Mdata(dcount)    = one
            Mindices(dcount) = 0
            Mindptr(rcount)  = Mindptr(rcount-1) + 1
            dcount = dcount + 1
            rcount = rcount + 1

            i = p+1
            Do k = 2,ng-1

               Do while (
     &              (xi_bar(k) >= Ukv_dir(i+1)).AND.(i<ng))
                  i = i+1 ! find knot interval
               Enddo
               CALL dersbasisfuns(i,p,n,xi_bar(k),Ukv_dir(:n),FN(:p+1),
     &              dNdxi(:p+1))

               Mdata(dcount:dcount+p)    = FN(:p+1)
               Mindices(dcount:dcount+p) = (/(j,j=i,i+p)/) - (p+1)
               Mindptr(rcount) = Mindptr(rcount-1) + p+1
               dcount = dcount + p+1
               rcount = rcount + 1
            Enddo

            Mdata(dcount)    = one
            Mindices(dcount) = ng-1
            Mindptr(rcount)  = Mindptr(rcount-1) + 1
            dcount = dcount + 1
            rcount = rcount + 1

         Enddo ! loop on parametric direction

         call finalizeNurbsPatch()

         ENDIF ! active patch
      ENDDO ! loop on patch


C     ------------------------------------------------------------------

      End subroutine buildGrevInterpolMat
















c     ...
c     Construction des  matrices pour l'interpolation/approximation de
c     B-Spline avec comme coordonnees parametriques des donnees difinies
c     par l'utilisateur
c      --> Renvoie les matrices au format creux CSR
c     ...

      subroutine buildUserInterpolMat(Mdata,Mindices,Mindptr,
     1     nb_data,nb_row, activePatch, xi_bar,nb_rowByDir,nbactiveDir,
     2     IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,ELT_TYPE,TENSOR,
     3     PROPS,JPROPS,NNODE,nb_patch,nb_elem)

      use parameters
      use nurbspatch

      Implicit none

C     ------------------------------------------------------------------

C     Input arguments :
c     ---------------
!     Geometry NURBS
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
      Integer,intent(in) :: nb_data,nb_row,activePatch
      dimension activePatch(nb_patch)

      Integer,          intent(in) :: nb_rowByDir,nbactiveDir
      Double precision, intent(in) :: xi_bar
      dimension xi_bar(nb_row),nb_rowByDir(nbactiveDir)

!     Output arguments :
!     ----------------
      Double precision, intent(out):: Mdata
      Integer,          intent(out):: Mindices,Mindptr
      dimension Mdata(nb_data),Mindices(nb_data),Mindptr(nb_row+1)


C     Local variables :
c     ---------------

!     Nurbs basis functions
      Integer          :: n,p,ng,direction
      Double precision :: Ukv_dir,FN,dNdxi
      dimension Ukv_dir(MAXVAL(Nkv)),FN(MAXVAL(Jpqr)+1),
     &     dNdxi(MAXVAL(Jpqr)+1)

      Integer :: i,j,k,numPatch,rcount,dcount,countGlo


C     ------------------------------------------------------------------


      Mdata(:)    = zero
      Mindices(:) = 0
      Mindptr(:)  = 0

      rcount   = 1
      dcount   = 1
      countGlo = 1
      DO numPatch = 1,nb_patch

         IF (activePatch(numPatch) == 1) then

         ! extract infos
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

         Do direction = 1,dim_patch

!     Get Greville abscissae
            p = Jpqr_patch(direction)
            n =  Nkv_patch(direction)
            ng= n - (p+1)
            If     (direction == 1) then
               Ukv_dir(:n) = Ukv1_patch(:)
            elseif (direction == 2) then
               Ukv_dir(:n) = Ukv2_patch(:)
            else
               Ukv_dir(:n) = Ukv3_patch(:)
            Endif

!     Evaluate bspline basis functions
            i = p+1
            Do k = 1,nb_rowByDir(countGlo)

               Do while (
     &              (xi_bar(rcount) >= Ukv_dir(i+1)).AND.(i<ng))
                  i = i+1 ! find knot interval
               Enddo
               CALL dersbasisfuns(i,p,n,xi_bar(rcount),Ukv_dir(:n),
     &              FN(:p+1),dNdxi(:p+1))

               Mdata(dcount:dcount+p)    = FN(:p+1)
               Mindices(dcount:dcount+p) = (/(j,j=i,i+p)/) - (p+1)
               Mindptr(rcount+1) = Mindptr(rcount) + p+1
               dcount = dcount + p+1
               rcount = rcount + 1


            Enddo
            countGlo = countGlo + 1

         Enddo ! loop on parametric direction

         call finalizeNurbsPatch()

         ENDIF ! active patch
      ENDDO ! loop on patch


C     ------------------------------------------------------------------

      End subroutine buildUserInterpolMat







