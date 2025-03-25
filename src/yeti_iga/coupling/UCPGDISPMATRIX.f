!! Copyright 2016-2019 Thibaut Hirschler

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

      SUBROUTINE CPLINGDISP(Rl,detJ,XI,BI,dim,dimInterface,MCRD,NNODE,
     &     NNODE_l,COORDS,CMAT)

      use parameters

      Implicit None

c     Input arguments :
c     ---------------
      Integer,         intent(in) :: dim,dimInterface,MCRD,NNODE,NNODE_l
      Double precision,intent(in) :: Rl,detJ,BI,COORDS
      double precision, intent(inout) :: XI
      dimension Rl(NNODE_l),XI(3),BI(3,dimInterface),COORDS(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: CMAT
      dimension CMAT(MCRD,MCRD,NNODE*NNODE_l)

c     Local variables :
c     ---------------

!     Nurbs basis fcts
      Double precision :: R,dRdxi
      dimension R(NNODE),dRdxi(NNODE,3)

!     Curivilinear basis vectos
      Double precision :: AI,TI,vectV,normV
      dimension AI(3,3),TI(3,3),vectV(3)

!     Assembly
      Double precision :: temp,temp2
      Integer          :: i,j,k,count

C     ------------------------------------------------------------------

c     Initialization :
c     --------------




c
c     ..................................................................
c
C     Computation :
c     -----------

c     Force XI values in interval [0, 1]
      XI = max(0.0, min(1.0, XI))

      call evalnurbs(XI(:),R(:),dRdxi(:,:))

      AI(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            AI(:MCRD,j) = AI(:MCRD,j) + dRdxi(i,j)*COORDS(:,i)
         Enddo
      Enddo


      TI(:,:) = zero
      Do j = 1,dimInterface
         Do i = 1,dim
            TI(:MCRD,j) = TI(:MCRD,j) + BI(i,j)*AI(:MCRD,i)
         Enddo
      Enddo

      If     (dimInterface == 1) then
         call norm(TI(:,1),3, normV)
      Elseif (dimInterface == 2) then
         call cross(TI(:,1),TI(:,2),vectV(:))
         call norm(vectV(:),3, normV)
      Else
         call cross(TI(:,1),TI(:,2),vectV(:))
         call dot(  TI(:,3),vectV(:),normV)
      Endif


      ! Assembling
      count = 1
      Do j  = 1,NNODE_l
         temp = Rl(j)*normV*detJ
         Do i = 1,NNODE
            temp2 = temp*R(i)
            Do k = 1,MCRD
               CMAT(k,k,count) = temp2
               if (temp2<0.) then
                  write(*,*) "[WARNING] Negative contribution to ",
     &                    "coupling matrix detected in cplingdisp"
                  write(*,*) "(computed at xi = ", XI(:), " )"
               endif
            Enddo
            count = count + 1
         Enddo
      Enddo



c
c     ..................................................................
c

      END SUBROUTINE CPLINGDISP





























      SUBROUTINE CPLINGDISP_embedded(Rl,detJ,XI,BI,dim,dimInterface,
     &     MCRD,NNODE,NNODE_l,NNODEmap,nb_cp,COORDS,COORDSall,CMAT)

      use parameters
      use embeddedMapping

      Implicit None

c     Input arguments :
c     ---------------
      Integer,         intent(in) :: dim,dimInterface,MCRD,NNODE,
     &     NNODE_l,NNODEmap,nb_cp
      Double precision,intent(in) :: Rl,detJ,XI,BI,COORDS,COORDSall
      dimension Rl(NNODE_l),XI(3),BI(3,dimInterface),COORDS(MCRD,NNODE),
     &     COORDSall(3,nb_cp)


c     Output variables :
c     ----------------
      Double precision, intent(out) :: CMAT
      dimension CMAT(MCRD,MCRD,NNODE*NNODE_l)

c     Local variables :
c     ---------------

!     Nurbs basis fcts
      Double precision :: R,dRdxi
      dimension R(NNODE),dRdxi(NNODE,3)

!     Curivilinear basis vectors
      Double precision :: AI,TI,VI,vectV,normV
      dimension AI(3,3),TI(3,3),VI(3,3),vectV(3)

!     Mapping
      Integer          :: isave,sctr_map
      dimension sctr_map(NNODEmap)
      Double precision :: Re,dRedxi,COORDSmap
      dimension Re(NNODEmap),dRedxi(NNODEmap,3),COORDSmap(MCRD,NNODEmap)

!     Assembly
      Double precision :: temp,temp2
      Integer          :: i,j,k,count

C     ------------------------------------------------------------------

c     Initialization :
c     --------------




c
c     ..................................................................
c
C     Computation :
c     -----------

      call evalnurbs(XI(:),R(:),dRdxi(:,:))

      vectV(:) = zero
      Do i = 1,NNODE
         vectV(:MCRD) = vectV(:MCRD) + R(i)*COORDS(:,i)
      Enddo

      AI(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            AI(:,j) = AI(:,j) + dRdxi(i,j)*COORDS(:,i)
         Enddo
      Enddo

c     Computing NURBS basis functions and derivatives of the mapping
!     get active element number
      call updateMapElementNumber(VectV(:))
      call evalnurbs_mapping(VectV(:),Re(:),dRedxi(:,:))

!     extract COORDS
      sctr_map(:) = IEN_map(:,current_map_elem)
      Do i = 1,NNODEmap
         COORDSmap(:,i) = COORDSall(:,sctr_map(i))
      Enddo

      VI(:,:) = zero
      Do j = 1,dim_map
         Do i = 1,NNODEmap
            VI(:,j) = VI(:,j) + dRedxi(i,j)*COORDSmap(:,i)
         Enddo
      Enddo


      TI(:,:) = zero
      Do j = 1,dimInterface

         VectV(:) = zero
         Do i = 1,dim
            VectV(:) = VectV(:) + BI(i,j)*AI(:,i)
         Enddo

         Do k = 1,dim_map
            TI(:,j) = TI(:,j) + VectV(k)*VI(:,k)
         Enddo

      Enddo

      If     (dimInterface == 1) then
         call norm(TI(:,1),3, normV)
      Elseif (dimInterface == 2) then
         call cross(TI(:,1),TI(:,2),vectV(:))
         call norm(vectV(:),3, normV)
      Else
         call cross(TI(:,1),TI(:,2),vectV(:))
         call dot(  TI(:,3),vectV(:),normV)
      Endif


      ! Assembling
      count = 1
      Do j  = 1,NNODE_l
         temp = Rl(j)*normV*detJ
         Do i = 1,NNODE
            temp2 = temp*R(i)
            Do k = 1,MCRD
               CMAT(k,k,count) = temp2
            Enddo
            count = count + 1
         Enddo
      Enddo



c
c     ..................................................................
c

      END SUBROUTINE CPLINGDISP_embedded
























      SUBROUTINE CPLINGDISPderv(Rl,detJ,XI,BI,dim,dimInterface,MCRD,
     &     NNODE,NNODE_l,COORDS,CMAT)

      use parameters

      Implicit None

c     Input arguments :
c     ---------------
      Integer,         intent(in) :: dim,dimInterface,MCRD,NNODE,NNODE_l
      Double precision,intent(in) :: Rl,detJ,XI,BI,COORDS
      dimension Rl(NNODE_l),XI(3),BI(3,dimInterface),COORDS(MCRD,NNODE)

c     Output variables :
c     ----------------
      Double precision, intent(out) :: CMAT
      dimension CMAT(MCRD,MCRD,NNODE*NNODE_l)

c     Local variables :
c     ---------------

!     Nurbs basis fcts
      Double precision :: R,dRdxi
      dimension R(NNODE),dRdxi(NNODE,3)

!     Curivilinear basis vectos
      Double precision :: AI,TI,vectV,normV, normB,NI
      dimension AI(3,3),TI(3,3),vectV(3),NI(3)

!     Assembly
      Double precision :: temp,temp2
      Integer          :: i,j,k,count

C     ------------------------------------------------------------------

c     Initialization :
c     --------------




c
c     ..................................................................
c
C     Computation :
c     -----------

      call evalnurbs(XI(:),R(:),dRdxi(:,:))

      AI(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            AI(:MCRD,j) = AI(:MCRD,j) + dRdxi(i,j)*COORDS(:,i)
         Enddo
      Enddo


      TI(:,:) = zero
      Do j = 1,dimInterface
         Do i = 1,dim
            TI(:MCRD,j) = TI(:MCRD,j) + BI(i,j)*AI(:MCRD,i)
         Enddo
      Enddo

      If     (dimInterface == 1) then
         call norm(TI(:,1),3, normV)
      Elseif (dimInterface == 2) then
         call cross(TI(:,1),TI(:,2),vectV(:))
         call norm(vectV(:),3, normV)
      Else
         call cross(TI(:,1),TI(:,2),vectV(:))
         call dot(  TI(:,3),vectV(:),normV)
      Endif

      call norm(BI(:,1),3, normB)
      NI(:) = BI(:,1)/normB

      ! Assembling
      count = 1
      Do j  = 1,NNODE_l
         temp = Rl(j)*normV*detJ
         Do i = 1,NNODE
            temp2 = temp*(NI(1)*dRdxi(i,1) + NI(2)*dRdxi(i,2))
            Do k = 1,MCRD
               CMAT(k,k,count) = temp2
            Enddo
            count = count + 1
         Enddo
      Enddo



c
c     ..................................................................
c

      END SUBROUTINE CPLINGDISPderv
