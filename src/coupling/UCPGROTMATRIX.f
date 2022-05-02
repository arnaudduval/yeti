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

      SUBROUTINE CPLINGROT(Rl,detJ,XI,BI,dim,dimInterface,MCRD,NNODE,
     &     NNODE_l,COORDS,CMAT)
      
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
      Double precision :: AI,TI,normA
      dimension AI(3,3),TI(3,3)
      
!     Rotation field
      Double precision :: vectT,normT,coef1,coef2,invA
      dimension vectT(3)

!     Assembly
      Double precision :: temp
      Integer          :: i,j,k,count

C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      
      
      CMAT(:,:,:) = zero

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
      call cross(AI(:,1),AI(:,2),AI(:,3))

      call norm(AI(:,3),3, normA)
      AI(:,3) = AI(:,3)/normA
      invA = one/normA
      
      TI(:,:) = zero
      Do i = 1,dim
         TI(:MCRD,1) = TI(:MCRD,1) + BI(i,1)*AI(:MCRD,i)
      Enddo
      call norm(TI(:,1),3, normT)
      vectT(:) = TI(:,1)/normT
            
      call dot(AI(:,1),vectT(:), coef1)
      call dot(AI(:,2),vectT(:), coef2)
      
      
      ! Assembling
      count = 1
      Do j  = 1,NNODE_l
         temp = Rl(j)*normT*detJ
         Do i = 1,NNODE
            CMAT(:,1,count) = invA*AI(:,3)*
     &           ( dRdxi(i,2)*coef1 - dRdxi(i,1)*coef2 ) * temp
            
            count = count + 1
         Enddo
      Enddo
      
      
      
c     
c     ..................................................................
c
      
      END SUBROUTINE CPLINGROT





























      SUBROUTINE CPLINGROT_embedded(Rl,detJ,XI,BI,dim,dimInterface,
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

!     Rotation field
      Double precision :: vectT,normT,coef1,coef2,invA,A1,A2,A3,normA
      dimension vectT(3),A1(3),A2(3),A3(3)
      
!     Assembly
      Double precision :: temp
      Integer          :: i,j,k,count

C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      
      CMAT(:,:,:) = zero
      

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
      
      A1(:) = zero
      A2(:) = zero
      Do j = 1,dim_map
         A1(:) = A1(:) + AI(j,1)*VI(:,j)
         A2(:) = A2(:) + AI(j,2)*VI(:,j)
      Enddo
      call cross(A1(:),A2(:),A3(:))
      call norm(A3(:),3, normA)
      A3(:) = A3(:)/normA
      invA = one/normA
      
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
      call norm(TI(:,1),3, normT)
      vectT(:) = TI(:,1)/normT
      
      call dot(A1(:),vectT(:), coef1)
      call dot(A2(:),vectT(:), coef2)
      
      
      ! Assembling
      count = 1
      Do j  = 1,NNODE_l
         temp = Rl(j)*normT*detJ
         Do i = 1,NNODE
            CMAT(:,1,count) = invA*A3(:)*
     &           ( dRdxi(i,2)*coef1 - dRdxi(i,1)*coef2 ) * temp
            
            count = count + 1
         Enddo
      Enddo
      
      
      
c     
c     ..................................................................
c
      
      END SUBROUTINE CPLINGROT_embedded
