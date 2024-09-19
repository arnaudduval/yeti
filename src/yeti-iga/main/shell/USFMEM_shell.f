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
      
      Subroutine usfmem_shell(NNODE,MCRD,NDOFEL,matH,AI,dRdxi,stiff)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD, NDOFEL
      Double precision, intent(in) :: matH, AI, dRdxi
      dimension matH(3,3), AI(3,3), dRdxi(NNODE,2)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(NDOFEL, NDOFEL)
      
c     Local variables :
c     ---------------
      
      Double precision :: zero
      parameter (zero=0.0D0)
      
      Double precision :: BoJ,SBoJ,BoI, dNidxi, dNjdxi, stiffLoc
      dimension BoJ(MCRD,MCRD),SBoJ(MCRD,MCRD),BoI(MCRD,MCRD),dNidxi(2),
     &      dNjdxi(2),stiffLoc(MCRD,MCRD)
      
      Integer :: i, j, nodj, incr_col, jdim, jdof, jcol, nodi, incr_row,
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
c     
c     Calcul rigidite flexion elementaire ..............................
c     
c     Boucle element
      jcol = 0
      do nodj = 1,NNODE
!     incr_col = (nodj-1)*MCRD
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
         enddo  
         
c     Calcul matrice BoJ
         BoJ(1,:) = dNjdxi(1)*AI(1,:)
         BoJ(2,:) = dNjdxi(2)*AI(2,:)
         BoJ(3,:) = dNjdxi(2)*AI(1,:) + dNjdxi(1)*AI(2,:)
         
c     Produit matH*BoJ
         call MulMat(matH,BoJ, SBoJ, MCRD,MCRD,MCRD)
         
c     Deuxieme boucle element
         irow = 0
         do nodi = 1,nodj !NNODE
!     incr_row = (nodi-1)*MCRD
            do idim = 1,2
               dNidxi(idim) = dRdxi(nodi,idim)
            enddo
            
c     Calcul matrice BoI
            BoI(:,1) = dNidxi(1)*AI(1,:)
            BoI(:,2) = dNidxi(2)*AI(2,:)
            BoI(:,3) = dNidxi(2)*AI(1,:) + dNidxi(1)*AI(2,:)
            
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
      End subroutine usfmem_shell
      






















      
      
C     ******************************************************************
      
      Subroutine usfmem_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dRdxi,
     &     stiff)
      
      Implicit None

c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NNODE, MCRD, NDOFEL
      Double precision, intent(in) :: matH, AI, dRdxi
      dimension matH(3,3), AI(3,3), dRdxi(NNODE,2)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: stiff
      dimension stiff(MCRD,MCRD,NNODE*(NNODE+1)/2)
      
c     Local variables :
c     ---------------
      
      Double precision :: zero
      parameter (zero=0.0D0)
      
      Double precision :: BoJ,SBoJ,BoI, dNidxi, dNjdxi, stiffLoc
      dimension BoJ(MCRD,MCRD),SBoJ(MCRD,MCRD),BoI(MCRD,MCRD),dNidxi(2),
     &      dNjdxi(2),stiffLoc(MCRD,MCRD)
      
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
c     
c     Calcul rigidite flexion elementaire ..............................
c     
c     Boucle cps
      count = 1
      do nodj = 1,NNODE
         do jdim = 1,2
            dNjdxi(jdim) = dRdxi(nodj,jdim)
         enddo  
         
c     Calcul matrice BoJ
         BoJ(1,:) = dNjdxi(1)*AI(1,:)
         BoJ(2,:) = dNjdxi(2)*AI(2,:)
         BoJ(3,:) = dNjdxi(2)*AI(1,:) + dNjdxi(1)*AI(2,:)
         
c     Produit matH*BoJ
         call MulMat(matH,BoJ, SBoJ, MCRD,MCRD,MCRD)
         
c     Deuxieme boucle cps
         do nodi = 1,nodj
            do idim = 1,2
               dNidxi(idim) = dRdxi(nodi,idim)
            enddo
            
c     Calcul matrice BoI
            BoI(:,1) = dNidxi(1)*AI(1,:)
            BoI(:,2) = dNidxi(2)*AI(2,:)
            BoI(:,3) = dNidxi(2)*AI(1,:) + dNidxi(1)*AI(2,:)
            
c     Calcul stiffLoc
            call MulMat(BoI,SBoJ, stiffLoc, MCRD,MCRD,MCRD)
            
c     Actualisation matrice elementaire
            stiff(:,:,count) = stiffLoc(:,:)
            
            count = count + 1
         enddo
      enddo
      
      return
      End subroutine usfmem_shell_byCP
      

