!! Copyright 2017-2018 Thibaut Hirschler

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

      
C     ******************************************************************
      
      Subroutine usfbnd_shell(NNODE,MCRD,NDOFEL,matH,AI,dAI1dxi,dAI2dxi,
     1     dRdxi,ddRddxi,stiff)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD, NDOFEL
      Double precision, intent(in) :: matH, AI, dAI1dxi, dAI2dxi, dRdxi,
     &     ddRddxi
      dimension matH(3,3), AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3),
     &     dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(NDOFEL, NDOFEL)
      
c     Local variables :
c     ---------------
      
      Double precision :: zero
      parameter (zero=0.0D0)
      
      Double precision :: stiffLoc, BoJ, SBoJ, BoI, dNjdxi, ddNjddxi,
     &     dNidxi, ddNiddxi, B1,B2,B3, Area, dA1d1_A2,dA2d2_A2,dA1d2_A2,
     &     A1_dA1d1, A1_dA2d2, A1_dA1d2, A2_A3, A3_A1, A3dA1d1, A3dA2d2,
     &     A3dA1d2
      dimension stiffLoc(MCRD,MCRD), BoJ(MCRD,MCRD), SBoJ(MCRD,MCRD),
     &     BoI(MCRD,MCRD), dNjdxi(2),ddNjddxi(3), dNidxi(2),ddNiddxi(3),
     &     B1(3), B2(3), B3(3), dA1d1_A2(3), dA2d2_A2(3), dA1d2_A2(3),
     &     A1_dA1d1(3), A1_dA2d2(3), A1_dA1d2(3), A2_A3(3), A3_A1(3)
      
      Integer :: nodj, incr_col, jdim, jdof, jcol, nodi, incr_row,
     &     idim, icol, idof, irow 
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
c     Initialisation des matrices stiff, BOI, BOJ
      stiff(:,:)= zero
      BoJ(:,:)  = zero
      BoI(:,:)  = zero 

c     Calculs preliminaire
      call SurfElem(AI(1,:), AI(2,:), Area)
      call cross(dAI1dxi(1,:), AI(2,:), dA1d1_A2(:))
      call cross(dAI2dxi(2,:), AI(2,:), dA2d2_A2(:))
      call cross(dAI1dxi(2,:), AI(2,:), dA1d2_A2(:))
      call cross(AI(1,:), dAI1dxi(1,:), A1_dA1d1(:))
      call cross(AI(1,:), dAI2dxi(2,:), A1_dA2d2(:))
      call cross(AI(1,:), dAI1dxi(2,:), A1_dA1d2(:))
      call cross(AI(2,:), AI(3,:), A2_A3(:))
      call cross(AI(3,:), AI(1,:), A3_A1(:))
      call dot(AI(3,:), dAI1dxi(1,:), A3dA1d1)
      call dot(AI(3,:), dAI2dxi(2,:), A3dA2d2)
      call dot(AI(3,:), dAI1dxi(2,:), A3dA1d2)

c     
c     Calcul rigidite flexion elementaire ..............................
c     
c     Boucle element
      jcol = 0
      do nodj = 1,NNODE
!     incr_col = (nodj-1)*MCRD
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
            ddNjddxi(jdim) = ddRddxi(nodj,jdim)
         enddo
         ddNjddxi(3) = ddRddxi(nodj,3)
         
c     Calcul matrice BoJ
         B1(:) = -ddNjddxi(1)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d1_A2(:) + dNjdxi(2)*A1_dA1d1(:)
     &        + A3dA1d1*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B2(:) = -ddNjddxi(2)*AI(3,:) +
     &        ( dNjdxi(1)*dA2d2_A2(:) + dNjdxi(2)*A1_dA2d2(:)
     &        + A3dA2d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B3(:) = -ddNjddxi(3)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d2_A2(:) + dNjdxi(2)*A1_dA1d2(:)
     &        + A3dA1d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         BoJ(1,:) = B1(:)
         BoJ(2,:) = B2(:)
         BoJ(3,:) = 2.0D0*B3(:)
         
c     Produit matH*BoJ
         call MulMat(matH,BoJ, SBoJ, MCRD,MCRD,MCRD)

         
c     Deuxieme boucle element
         irow = 0
         do nodi = 1,nodj !NNODE
!     incr_row = (nodi-1)*MCRD
            do idim = 1,2
               dNidxi(idim) = dRdxi(nodi,idim)
               ddNiddxi(idim) = ddRddxi(nodi,idim)
            enddo
            ddNiddxi(3) = ddRddxi(nodi,3)
            
c     Calcul matrice BoI
            B1(:) = -ddNiddxi(1)*AI(3,:) +
     &           ( dNidxi(1)*dA1d1_A2(:) + dNidxi(2)*A1_dA1d1(:)
     &           + A3dA1d1*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            B2(:) = -ddNiddxi(2)*AI(3,:) +
     &           ( dNidxi(1)*dA2d2_A2(:) + dNidxi(2)*A1_dA2d2(:)
     &           + A3dA2d2*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            B3(:) = -ddNiddxi(3)*AI(3,:) +
     &           ( dNidxi(1)*dA1d2_A2(:) + dNidxi(2)*A1_dA1d2(:)
     &           + A3dA1d2*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            BoI(:,1) = B1(:)
            BoI(:,2) = B2(:)
            BoI(:,3) = 2.0D0*B3(:)
            
c     Calcul stiffLoc
            !call MulATBA(stiffLoc,BoI,matH,BoJ,MCRD,MCRD)
            call MulMat(BoI,SBoJ, stiffLoc, MCRD,MCRD,MCRD)
            
c     Actualisation matrice elementaire
            stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)
     &           = stiff(irow+1:irow+MCRD,jcol+1:jcol+MCRD)
     &           + stiffLoc(:,:)
            
            irow = irow + MCRD

!            do jdof = 1,MCRD
!               icol = jdof + incr_col
!               do idof = 1,MCRD
!                  irow = idof + incr_row
!                  stiff(irow,icol) = stiff(irow,icol)
!     &                 + stiffLoc(idof,jdof)
!               enddo
!            enddo
         enddo
         jcol = jcol + MCRD
      enddo
      
      return
      End subroutine usfbnd_shell






























C     ******************************************************************
      
      Subroutine usfbnd_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dAI1dxi,
     1     dAI2dxi,dRdxi,ddRddxi,stiff)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD, NDOFEL
      Double precision, intent(in) :: matH, AI, dAI1dxi, dAI2dxi, dRdxi,
     &     ddRddxi
      dimension matH(3,3), AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3),
     &     dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(MCRD,MCRD,NNODE*(NNODE+1)/2)
      
c     Local variables :
c     ---------------
      
      Double precision :: zero
      parameter (zero=0.0D0)
      
      Double precision :: stiffLoc, BoJ, SBoJ, BoI, dNjdxi, ddNjddxi,
     &     dNidxi, ddNiddxi, B1,B2,B3, Area, dA1d1_A2,dA2d2_A2,dA1d2_A2,
     &     A1_dA1d1, A1_dA2d2, A1_dA1d2, A2_A3, A3_A1, A3dA1d1, A3dA2d2,
     &     A3dA1d2
      dimension stiffLoc(MCRD,MCRD), BoJ(MCRD,MCRD), SBoJ(MCRD,MCRD),
     &     BoI(MCRD,MCRD), dNjdxi(2),ddNjddxi(3), dNidxi(2),ddNiddxi(3),
     &     B1(3), B2(3), B3(3), dA1d1_A2(3), dA2d2_A2(3), dA1d2_A2(3),
     &     A1_dA1d1(3), A1_dA2d2(3), A1_dA1d2(3), A2_A3(3), A3_A1(3)
      
      Integer :: nodj,jdim,nodi,idim,count
      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Initialisation ...................................................
      
c     Initialisation des matrices stiff, BOI, BOJ
      stiff(:,:,:)= zero
      BoJ(:,:)    = zero
      BoI(:,:)    = zero 

c     Calculs preliminaire
      call SurfElem(AI(1,:), AI(2,:), Area)
      call cross(dAI1dxi(1,:), AI(2,:), dA1d1_A2(:))
      call cross(dAI2dxi(2,:), AI(2,:), dA2d2_A2(:))
      call cross(dAI1dxi(2,:), AI(2,:), dA1d2_A2(:))
      call cross(AI(1,:), dAI1dxi(1,:), A1_dA1d1(:))
      call cross(AI(1,:), dAI2dxi(2,:), A1_dA2d2(:))
      call cross(AI(1,:), dAI1dxi(2,:), A1_dA1d2(:))
      call cross(AI(2,:), AI(3,:), A2_A3(:))
      call cross(AI(3,:), AI(1,:), A3_A1(:))
      call dot(AI(3,:), dAI1dxi(1,:), A3dA1d1)
      call dot(AI(3,:), dAI2dxi(2,:), A3dA2d2)
      call dot(AI(3,:), dAI1dxi(2,:), A3dA1d2)

c     
c     Calcul rigidite flexion elementaire ..............................
c     
c     Boucle cps
      count = 1
      do nodj = 1,NNODE
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
            ddNjddxi(jdim) = ddRddxi(nodj,jdim)
         enddo
         ddNjddxi(3) = ddRddxi(nodj,3)
         
c     Calcul matrice BoJ
         B1(:) = -ddNjddxi(1)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d1_A2(:) + dNjdxi(2)*A1_dA1d1(:)
     &        + A3dA1d1*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B2(:) = -ddNjddxi(2)*AI(3,:) +
     &        ( dNjdxi(1)*dA2d2_A2(:) + dNjdxi(2)*A1_dA2d2(:)
     &        + A3dA2d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         B3(:) = -ddNjddxi(3)*AI(3,:) +
     &        ( dNjdxi(1)*dA1d2_A2(:) + dNjdxi(2)*A1_dA1d2(:)
     &        + A3dA1d2*( dNjdxi(1)*A2_A3(:) + dNjdxi(2)*A3_A1(:) )
     &        )/Area
         BoJ(1,:) = B1(:)
         BoJ(2,:) = B2(:)
         BoJ(3,:) = 2.0D0*B3(:)
         
c     Produit matH*BoJ
         call MulMat(matH,BoJ, SBoJ, MCRD,MCRD,MCRD)

         
c     Deuxieme boucle cps
         do nodi = 1,nodj
            do idim = 1,2
               dNidxi(idim) = dRdxi(nodi,idim)
               ddNiddxi(idim) = ddRddxi(nodi,idim)
            enddo
            ddNiddxi(3) = ddRddxi(nodi,3)
            
c     Calcul matrice BoI
            B1(:) = -ddNiddxi(1)*AI(3,:) +
     &           ( dNidxi(1)*dA1d1_A2(:) + dNidxi(2)*A1_dA1d1(:)
     &           + A3dA1d1*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            B2(:) = -ddNiddxi(2)*AI(3,:) +
     &           ( dNidxi(1)*dA2d2_A2(:) + dNidxi(2)*A1_dA2d2(:)
     &           + A3dA2d2*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            B3(:) = -ddNiddxi(3)*AI(3,:) +
     &           ( dNidxi(1)*dA1d2_A2(:) + dNidxi(2)*A1_dA1d2(:)
     &           + A3dA1d2*( dNidxi(1)*A2_A3(:) + dNidxi(2)*A3_A1(:))
     &           )/Area
            BoI(:,1) = B1(:)
            BoI(:,2) = B2(:)
            BoI(:,3) = 2.0D0*B3(:)
            
c     Calcul stiffLoc
            call MulMat(BoI,SBoJ, stiffLoc, MCRD,MCRD,MCRD)
            
c     Actualisation matrice elementaire
            stiff(:,:,count) = stiffLoc(:,:)
            
            count = count + 1
         enddo
      enddo
      
      return
      End subroutine usfbnd_shell_byCP
