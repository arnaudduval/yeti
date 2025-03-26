!! Copyright 2017-2019 Thibaut Hirschler

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

C     ASSEMBLAGE DE LA MATRICE DE MASSE

c     Include iga subroutines
!     include "./shap.f"
!     include "./Gauss.f"
!     include "./operateurs.f"
!     include "./dersbasisfuns.f"
!     include "./UMASSMAT.f"



C     ******************************************************************


      Subroutine build_CMASSMatrix(Mdata,Mrow,Mcol,activeElement,
     1     nb_data,COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,Jpqr,
     2     ELT_TYPE,PROPS,JPROPS,RHO,TENSOR,ind_dof_free,nb_dof_free,
     3     MCRD,NBINT,nb_patch,nb_elem,nnode,nb_cp,nb_dof_tot)

      use parameters
      use nurbspatch

      Implicit None


c     Declaration des Variables ........................................

c     Input arguments :
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
      Character(len=*), intent(in) :: TENSOR,ELT_TYPE
      Double precision, intent(in) :: RHO,PROPS
      Integer, intent(in) :: MCRD,NNODE,nb_patch,NBINT,nb_elem,
     &     nb_elem_patch,IEN,JPROPS
      dimension RHO(nb_patch), PROPS(:), JPROPS(nb_patch),
     &     nb_elem_patch(nb_patch),IEN(:), NNODE(nb_patch),
     &     NBINT(nb_patch)


!     Degree Of Freedom
      Integer, intent(in) :: nb_dof_tot, nb_dof_free, ind_dof_free
      dimension ind_dof_free(nb_dof_tot)

!     Storage INFOS
      Integer, intent(in) :: nb_data,activeElement
      dimension activeElement(nb_elem)



c     Output variables : coefficient diag matrice de masse
c     ----------------
      Double precision, intent(out) :: Mdata
      integer, intent(out) :: Mrow, Mcol
      dimension Mdata(nb_data),Mrow(nb_data),Mcol(nb_data)




c     Local Variables :
c     ---------------

      ! for UMASSMAT.f
      Integer :: NDOFEL
      Double precision :: COORDS_elem,CMASSMATRX
      dimension CMASSMATRX(MCRD,MCRD,MAXVAL(NNODE)*(MAXVAL(NNODE)+1)/2),
     &     COORDS_elem(MCRD,MAXVAL(NNODE))

!     Assembly
      Integer :: num_elem,i,j,JELEM,NumPatch, sctr
      dimension sctr(MAXVAL(NNODE))
      Integer :: dofi,dofj,cpi,cpj, kk,ll, nnodeSum,count

C     Fin declaration des variables ....................................
c
c
c
c
c     Initialisation ...................................................

!     Initialize mass matrix
      Mdata = zero
      Mrow  = 0
      Mcol  = 0


c     Fin Initialisation ...............................................
c
c
c
c
c     Debut Assemblage .................................................

      count = 1
      JELEM = 0
c     Loop on patches
      Do NumPatch = 1,nb_patch

         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)

         NDOFEL   = nnode_patch*MCRD
         nnodeSum = nnode_patch*(nnode_patch+1)/2

c     Loop on elements
         Do num_elem = 1,nb_elem_patch(NumPatch)
            JELEM = JELEM + 1

            IF (activeElement(JELEM)==1) then

            Do i = 1,nnode_patch
               COORDS_elem(:,i) = COORDS3D(:MCRD,IEN_patch(i,num_elem))
            Enddo

            CALL extractNurbsElementInfos(num_elem)

c     Build elementary Matrix
            CMASSMATRX = zero
            call UMASSMAT_byCP(NDOFEL,MCRD,nnode_patch,NBINT,
     &           COORDS_elem(:,:nnode_patch),ELT_TYPE_patch,
     &           RHO(NumPatch),PROPS_patch,JPROPS_patch,
     &           CMASSMATRX(:,:,:nnodeSum))


c     Assemble DLMMATRX to global matrix M
            sctr(:nnode_patch) = IEN_patch(:,num_elem)

            i = 0
            Do cpj = 1,nnode_patch

               dofj= (sctr(cpj)-1)*MCRD

!     cas cpi < cpj
               Do cpi = 1,cpj-1
                  dofi= (sctr(cpi)-1)*MCRD
                  i   = i + 1
                  Do ll = 1,MCRD
                     Do kk = 1,MCRD
                        Mdata(count) = CMASSMATRX(kk,ll,i)
                        Mrow( count) = dofi + kk - 1
                        Mcol( count) = dofj + ll - 1
                        count = count + 1
                     Enddo
                  Enddo
               Enddo

!     cas cpi == cpj
               i = i + 1
               Do ll = 1,MCRD
                  CMASSMATRX(ll,ll,i) = CMASSMATRX(ll,ll,i)*0.5d0
                  Do kk = 1,ll
                     Mdata(count) = CMASSMATRX(kk,ll,i)
                     Mrow( count) = dofj + kk - 1
                     Mcol( count) = dofj + ll - 1
                     count = count + 1
                  Enddo
               Enddo

            Enddo

            Endif
         Enddo

         CALL finalizeNurbsPatch()

      Enddo


c     Fin Assemblage ...................................................

      End subroutine build_CMASSMatrix
