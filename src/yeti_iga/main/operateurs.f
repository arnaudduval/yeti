!! Copyright 2011 Thomas Elguedj
!! Copyright 2011 Florian Maurin
!! Copyright 2016-2018 Thibaut Hirschler
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
!! with Yeti. If not, see <https://www.gnu.org/licenses/>

cccccccccccccccccccc   Multiplication de matrice ccccccccccccccccc
c     
      SUBROUTINE MulMat(AA, BB, CC, mm, nn, ll)
      
      
      double precision AA, BB, CC
      Dimension AA(mm,ll)
      Dimension BB(ll,nn)
      Dimension CC(mm,nn)  
      
      integer ii, jj, kk, mm, nn, ll
      
      do ii=1,mm
            do jj=1,nn
                CC(ii,jj)=0.d0
                do kk=1,ll
                CC(ii,jj)=CC(ii,jj)+AA(ii,kk)*BB(kk,jj)                
                enddo     
            enddo        
      enddo
      RETURN
      END SUBROUTINE MulMat
c      
cccccccccccccccccccc   Multiplication de vecteur ccccccccccccccccc
c     
      SUBROUTINE MulVect(AA, BB,CC, mm, nn)
      
      
      double precision AA, BB, CC
      Dimension AA(mm,nn)
      Dimension BB(nn)
      Dimension CC(mm)
            
      integer ii, jj, mm, nn
      
      CC(:) = 0.d0
      do ii=1,mm
         CC(ii)=0.d0
         do jj=1,nn
            CC(ii) = CC(ii)+AA(ii,jj)*BB(jj)                
         enddo
      enddo
      RETURN
      END SUBROUTINE MulVect
c      
ccccccccccccc   Multiplication de vecteur type v=tA*u cccccccccccccccccc
c     
      SUBROUTINE MulATVect(AA, BB,CC, mm, nn)
            
      double precision AA, BB, CC
      Dimension AA(mm,nn)
      Dimension BB(mm)
      Dimension CC(nn)
            
      integer ii, jj, mm, nn
      
      CC(:) = 0.d0
      do jj=1,mm
         do ii=1,nn
            CC(ii) = CC(ii) + AA(jj,ii)*BB(jj)
         enddo
      enddo
      RETURN
      END SUBROUTINE MulATVect
c     
cccccccccccccccccccc   Initialisation matrice ccccccccccccccccc
c
     
      SUBROUTINE ZerMat(AA, mm, nn)
      
      
      double precision AA
      Dimension AA(mm,nn)
      
      integer ii, jj, mm, nn
      
      do ii=1,mm
        do jj=1,nn
        AA(ii,jj)=0.d0
        enddo
      enddo
      RETURN
      END SUBROUTINE ZerMat
c

cccccccccccccccccc  Initialize a matrix to identity   cccccccccc

      subroutine SetMatrixToIdentity(AA, nn)
      
      integer i, j, nn
      
      double precision AA
      dimension AA(nn,nn)
      
      AA(:,:) = 0.d0
      do i=1,nn
         AA(i,i) = 1.d0
      enddo
      
      end subroutine SetMatrixToIdentity


cccccccccccccccccccc   Surface Elementaire ccccccccccccccccc
c     
      SUBROUTINE SurfElem(BB, CC, DDetJ)
      
      Double precision :: DDetJ, AA, BB, CC
      Dimension AA(3), BB(3), CC(3)  
      
      AA(1) = BB(2)*CC(3) - BB(3)*CC(2)
      AA(2) = BB(3)*CC(1) - BB(1)*CC(3)
      AA(3) = BB(1)*CC(2) - BB(2)*CC(1)
      
!     DDetJ=sqrt((AA(1)**2.d0)+(AA(2)**2.d0)+(AA(3)**2.d0))
      DDetJ = sqrt(AA(1)*AA(1) + AA(2)*AA(2) + AA(3)*AA(3))
      
      RETURN
      END SUBROUTINE SurfElem
cccccccccccccccccccc   Vecteur tengent normalis� ccccccccccccccccc
c     
      SUBROUTINE VectTengNorm(AA, BB, CC,DDetJ)
      
      
      Double precision DDetJ, AA, BB, CC
      Dimension AA(3)
      Dimension BB(3)
      Dimension CC(3)  
      
      AA(1)=BB(2)*CC(3)-BB(3)*CC(2)
      AA(2)=BB(3)*CC(1)-BB(1)*CC(3)
      AA(3)=BB(1)*CC(2)-BB(2)*CC(1)
      DDetJ=sqrt((AA(1)**2.d0)+(AA(2)**2.d0)+(AA(3)**2.d0))
      AA(1)=BB(1)/DDetJ
      AA(2)=BB(2)/DDetJ
      AA(3)=BB(3)/DDetJ
      
      RETURN
      END SUBROUTINE VectTengNorm
c
cccccccccccccccccccc   Vecteur normal normalis� ccccccccccccccccc
c     
      SUBROUTINE VectNormNorm(AA, BB, CC, DDetJ)
      
      
      Double precision DDetJ, AA, BB, CC
      Dimension AA(3)
      Dimension BB(3)
      Dimension CC(3)  
      AA(1)=(BB(2)*CC(3)-BB(3)*CC(2))/DDetJ
      AA(2)=(BB(3)*CC(1)-BB(1)*CC(3))/DDetJ
      AA(3)=(BB(1)*CC(2)-BB(2)*CC(1))/DDetJ
      
      RETURN
      END SUBROUTINE VectNormNorm
c
cccccccccccccccccccc   Initialisation vecteur ccccccccccccccccc
c     
      SUBROUTINE ZerVec(AA, mm)
      
      
      integer mm, ii
      double precision AA
      Dimension AA(mm)
      do ii=1,mm
        AA(ii)=0.d0
      enddo
      RETURN
      END SUBROUTINE ZerVec
      
cccccccccccccccccccccccccc    SubRoutine d'afichage   cccccccccccccccccc 
      SUBROUTINE AFFICH (AA, mm, nn)
      
      
      integer mm, nn
      double precision AA
      Dimension AA(mm,nn)
      
      integer ii, jj
      
      DO ii=1, mm
      write(*,*) (AA(ii, jj), jj=1, nn)
1000  format(10F13.3)
      enddo
      RETURN
      END SUBROUTINE AFFICH
      

      
cccccccccccccccccccccccccc SubRoutine d'afichage dans un fichier cccccccccc

      SUBROUTINE Text (AA, mm, nn)
      
      
      integer mm, nn
      double precision AA
      Dimension AA(mm,nn)
      
      integer ii, jj, ios
      
c      open( unit=10,  file = "C:\Temp\UEL\new.txt",  iostat = ios)
      open( unit=10,  file = "temp/mat.txt",  iostat = ios)
      DO ii=1, mm
         DO jj=1, nn

            if (jj.ne.nn) then
               if (AA(ii, jj) == 0.0D0) then
                  write(10,'(2X)',ADVANCE="NO")
               else
                  write(10,1001,ADVANCE="NO") AA(ii, jj)
               end if
            else
               if (AA(ii, jj) == 0.0D0) then
                  write(10,'(2X)',ADVANCE="YES")
               else
                  write(10,1001,ADVANCE="YES") AA(ii, jj)
               end if
            endif 
1000        format(10F25.10)
1001        format(X,ES10.3) !EN1.0)
         enddo
      enddo
      close(10)
      RETURN
      END SUBROUTINE Text

cccccccccccccccccccc   Multiplication de type A=Bt*C*D ccccccccccccccccc
c     
      SUBROUTINE MulATBA(AA, BB, CC, DD, mm, nn)
      
      
      integer nn, mm, ii, jj, kk
      double precision AA, BB, CC, DD, CD
      
      Dimension AA(nn,nn)
      Dimension BB(mm,nn)
      Dimension CC(mm,mm)
      Dimension DD(mm,nn)
      Dimension CD(mm,nn)
      
      do ii = 1,mm
         do jj = 1,nn
            CD(ii,jj) = 0.d0
            do kk = 1,mm
               CD(ii,jj) = CD(ii,jj) + CC(ii,kk)*DD(kk,jj)
            enddo
         enddo
      enddo
      
c      AA(:,:) = 0.d0
c      do kk = 1,mm
      do ii = 1,nn
         do jj = 1,nn
            AA(ii,jj) = 0.d0
            do kk = 1,mm
               AA(ii,jj) = AA(ii,jj) + BB(kk,ii)*CD(kk,jj)
            enddo
         enddo
      enddo
      RETURN
      END SUBROUTINE MulATBA



      Subroutine MulATBA_intrinsic(AA, BB, CC, DD, mm, nn)
      
      implicit none
      
      integer :: mm,nn
      double precision, intent(in) :: BB, CC, DD
      double precision, intent(out):: AA
      Dimension AA(nn,nn),BB(mm,nn),CC(mm,mm),DD(mm,nn)
      
      AA = MATMUL(TRANSPOSE(BB),MATMUL(CC,DD))
      
      End subroutine MulATBA_intrinsic
      


      
ccccccccccccccccccc   calcul de raideur g�om�trique locale  ccccccccccccccc
      subroutine CalcStiffLoc(stiff, dNidx, sigma, dNjdx, nn)      
      
      
      double precision stiff, dNidx, sigma, dNjdx
      integer nn
      
      dimension stiff(nn,nn)
      dimension dNidx(nn)
      dimension sigma(nn,nn)
      dimension dNjdx(nn)
      
      double precision a
      integer kk, ll, ii, jj
      
      
      a = 0.d0
      do kk = 1,nn
        do ll=1, nn
          a = a + dNidx(kk)*sigma(kk,ll)*dNjdx(ll)
        enddo
      enddo

      
      do ii=1,nn
        do jj=1,nn
          if(ii.EQ.jj) then
            stiff(ii,jj) = a
          else
            stiff(ii,jj) = 0.d0
          endif
        enddo
      enddo

      return      
      end subroutine CalcStiffLoc
      
cccccccccc   Multiplication d'une matrice par sa transpos�e A = B.Bt ccccc
      subroutine mulBBt(AA, BB, nn)
        
        
        double precision AA, BB
        integer nn
        
        dimension AA(nn,nn)
        dimension BB(nn,nn)
        
        integer ii, jj, kk
        
        do ii=1, nn
          do jj=1, nn
            AA(ii,jj) = 0.d0
            do kk=1, nn
              AA(ii,jj) = AA(ii,jj)+BB(ii,kk)*BB(jj,kk)
            enddo
          enddo
        enddo

        return
      end subroutine mulBBt
      
cccccccccc   Multiplication d'une matrice par sa transpos�e A = Bt.B ccccc
      subroutine mulBtB(AA, BB, nn)
        
        
        double precision AA, BB
        integer nn
        
        dimension AA(nn,nn)
        dimension BB(nn,nn)
        
        integer ii, jj, kk
        
        do ii=1, nn
          do jj=1, nn
            AA(ii,jj) = 0.d0
            do kk=1, nn
              AA(ii,jj) = AA(ii,jj)+BB(kk,ii)*BB(kk,jj)
            enddo
          enddo
        enddo

        return
      end subroutine mulBtB      
      
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      SUBROUTINE MatrixInv(AA, BB,DetBB, MCRD)
      
      implicit none
      
      integer, intent(in) :: MCRD
      double precision, intent(out) :: AA
      double precision, intent(in) :: BB
      double precision, intent(out) :: DetBB
            
      dimension AA(MCRD,MCRD)
      dimension BB(MCRD,MCRD)
      
      

      if (MCRD==2) then
        DetBB=(BB(1,1)*BB(2,2))-(BB(2,1)*BB(1,2))
        AA(1,1)=BB(2,2)/DetBB
        AA(2,2)=BB(1,1)/DetBB
        AA(1,2)=-BB(1,2)/DetBB
        AA(2,1)=-BB(2,1)/DetBB
      else if (MCRD==3) then
       DetBB=BB(1,2)*BB(2,3)*BB(3,1)
     1           -BB(1,3)*BB(2,2)*BB(3,1)
     2           +BB(1,3)*BB(2,1)*BB(3,2)
     3           -BB(1,1)*BB(2,3)*BB(3,2)
     4           +BB(1,1)*BB(2,2)*BB(3,3)
     5           -BB(1,2)*BB(2,1)*BB(3,3)
c        if (abs(DetB).lt.1.d-40) then
c          write(*,'("Erreur: JAAcoBBien singulier")')
c        end if
       AA(1,1)=(BB(2,2)*BB(3,3)-BB(2,3)*BB(3,2))/DetBB
       AA(1,2)=(BB(1,3)*BB(3,2)-BB(1,2)*BB(3,3))/DetBB
       AA(1,3)=(BB(1,2)*BB(2,3)-BB(1,3)*BB(2,2))/DetBB
       AA(2,1)=(BB(2,3)*BB(3,1)-BB(2,1)*BB(3,3))/DetBB
       AA(2,2)=(BB(1,1)*BB(3,3)-BB(1,3)*BB(3,1))/DetBB
       AA(2,3)=(BB(1,3)*BB(2,1)-BB(1,1)*BB(2,3))/DetBB
       AA(3,1)=(BB(2,1)*BB(3,2)-BB(2,2)*BB(3,1))/DetBB
       AA(3,2)=(BB(1,2)*BB(3,1)-BB(1,1)*BB(3,2))/DetBB
       AA(3,3)=(BB(1,1)*BB(2,2)-BB(1,2)*BB(2,1))/DetBB
      else if (MCRD==4) then
       DetBB=BB(1,1)*BB(2,2)*BB(3,3)*BB(4,4)
     *         +BB(1,1)*BB(2,3)*BB(3,4)*BB(4,2)
     *         + BB(1,1)*BB(2,4)*BB(3,2)*BB(4,3)
     *         + BB(1,2)*BB(2,1)*BB(3,4)*BB(4,3)
     *         + BB(1,2)*BB(2,3)*BB(3,1)*BB(4,4)
     *         + BB(1,2)*BB(2,4)*BB(3,3)*BB(4,1)
     *         + BB(1,3)*BB(2,1)*BB(3,2)*BB(4,4)
     *         + BB(1,3)*BB(2,2)*BB(3,4)*BB(4,1)
     *         + BB(1,3)*BB(2,4)*BB(3,1)*BB(4,2)
     *         + BB(1,4)*BB(2,1)*BB(3,3)*BB(4,1)
     *        + BB(1,4)*BB(2,2)*BB(3,1)*BB(4,3)
     *        + BB(1,4)*BB(2,3)*BB(3,2)*BB(4,1)
     *        - BB(1,1)*BB(2,2)*BB(3,4)*BB(4,3)
     *        - BB(1,1)*BB(2,3)*BB(3,2)*BB(4,4)
     *        - BB(1,1)*BB(2,4)*BB(3,3)*BB(4,2)
     *        - BB(1,2)*BB(2,1)*BB(3,3)*BB(4,4)
     *        - BB(1,2)*BB(2,3)*BB(3,4)*BB(4,1)
     *        - BB(1,2)*BB(2,4)*BB(3,1)*BB(4,3)
     *        - BB(1,3)*BB(2,1)*BB(3,4)*BB(4,2)
     *        - BB(1,3)*BB(2,2)*BB(3,1)*BB(4,4)
     *        - BB(1,3)*BB(2,4)*BB(3,2)*BB(4,1)
     *        - BB(1,4)*BB(2,1)*BB(3,2)*BB(4,3)
     *        - BB(1,4)*BB(2,2)*BB(3,3)*BB(4,1)
     *        - BB(1,4)*BB(2,3)*BB(3,1)*BB(4,2)
c        if (abs(DetB).lt.1.d-40) then
c          write(*,'("Erreur: JAAcoBBien singulier")')
c        end if
      AA(1,1) = ( BB(2,2)*BB(3,3)*BB(4,4) + BB(2,3)*BB(3,4)*BB(4,2) 
     * + BB(2,4)*BB(3,2)*BB(4,3) - BB(2,2)*BB(3,4)*BB(4,3) 
     * - BB(2,3)*BB(3,2)*BB(4,4) - BB(2,4)*BB(3,3)*BB(4,2) ) * (1/DetBB)
      AA(1,2) = (BB(1,2)*BB(3,4)*BB(4,3) + BB(1,3)*BB(3,2)*BB(4,4) 
     * + BB(1,4)*BB(3,3)*BB(4,2) - BB(1,2)*BB(3,3)*BB(4,4) 
     * - BB(1,3)*BB(3,4)*BB(4,2) - BB(1,4)*BB(3,2)*BB(4,3) ) * (1/DetBB)
      AA(1,3) = (BB(1,2)*BB(2,3)*BB(4,4) + BB(1,3)*BB(2,4)*BB(4,2) 
     * + BB(1,4)*BB(2,2)*BB(4,3) - BB(1,2)*BB(2,4)*BB(4,3) 
     * - BB(1,3)*BB(2,2)*BB(4,4) - BB(1,4)*BB(2,3)*BB(4,2) ) * (1/DetBB)    
      AA(1,4) = (BB(1,2)*BB(2,4)*BB(3,3) + BB(1,3)*BB(2,2)*BB(3,4) 
     * + BB(1,4)*BB(2,3)*BB(3,2) - BB(1,2)*BB(2,3)*BB(3,4) 
     * - BB(1,3)*BB(2,4)*BB(3,2) - BB(1,4)*BB(2,2)*BB(3,3) ) * (1/DetBB) 
      AA(2,1) = (BB(2,1)*BB(3,4)*BB(4,3) + BB(2,3)*BB(3,1)*BB(4,4) 
     * + BB(2,4)*BB(3,3)*BB(4,1) - BB(2,1)*BB(3,3)*BB(4,4) 
     * - BB(2,3)*BB(3,4)*BB(4,1) - BB(2,4)*BB(3,1)*BB(4,3) ) * (1/DetBB)
      AA(2,2) = (BB(1,1)*BB(3,3)*BB(4,4) + BB(1,3)*BB(3,4)*BB(4,1) 
     * + BB(1,4)*BB(3,1)*BB(4,3) - BB(1,1)*BB(3,4)*BB(4,3) 
     * - BB(1,3)*BB(3,1)*BB(4,4) - BB(1,4)*BB(3,3)*BB(4,1) ) * (1/DetBB)
      AA(2,3) = (BB(1,1)*BB(2,4)*BB(4,3) + BB(1,3)*BB(2,1)*BB(4,4) 
     * + BB(1,4)*BB(2,3)*BB(4,1) - BB(1,1)*BB(2,3)*BB(4,4) 
     * - BB(1,3)*BB(2,4)*BB(4,1) - BB(1,4)*BB(2,1)*BB(4,3) ) * (1/DetBB)
      AA(2,4) = (BB(1,1)*BB(2,3)*BB(3,4) + BB(1,3)*BB(2,4)*BB(3,1) 
     * + BB(1,4)*BB(2,1)*BB(3,3) - BB(1,1)*BB(2,4)*BB(3,3) 
     * - BB(1,3)*BB(2,1)*BB(3,4) - BB(1,4)*BB(2,3)*BB(3,1) ) * (1/DetBB)
      AA(3,1) = (BB(2,1)*BB(3,2)*BB(4,4) + BB(2,2)*BB(3,4)*BB(4,1) 
     * + BB(2,4)*BB(3,1)*BB(4,2) - BB(2,1)*BB(3,4)*BB(4,2) 
     * - BB(2,1)*BB(3,3)*BB(4,4) - BB(2,4)*BB(3,2)*BB(4,1) ) * (1/DetBB)
      AA(3,2) = (BB(1,1)*BB(3,4)*BB(4,2) + BB(1,2)*BB(3,1)*BB(4,4) 
     * + BB(1,4)*BB(3,2)*BB(4,1) - BB(1,1)*BB(3,2)*BB(4,4) 
     * - BB(1,2)*BB(3,2)*BB(4,4) - BB(1,2)*BB(3,4)*BB(4,1) ) * (1/DetBB)
      AA(3,3) = (BB(1,1)*BB(2,2)*BB(4,4) + BB(1,2)*BB(2,4)*BB(4,1) 
     * + BB(1,4)*BB(2,1)*BB(4,2) - BB(1,1)*BB(2,4)*BB(4,2) 
     * - BB(1,2)*BB(2,1)*BB(4,4) - BB(1,4)*BB(2,2)*BB(4,1) ) * (1/DetBB)
      AA(3,4) = (BB(1,1)*BB(2,4)*BB(3,2) + BB(1,2)*BB(2,1)*BB(3,4) 
     * + BB(1,4)*BB(2,2)*BB(3,1) - BB(1,1)*BB(2,2)*BB(3,4) 
     * - BB(1,2)*BB(2,4)*BB(3,1) - BB(1,4)*BB(2,1)*BB(3,2) ) * (1/DetBB)
      AA(4,1) = (BB(2,1)*BB(3,3)*BB(4,2) + BB(2,2)*BB(3,1)*BB(4,3) 
     * + BB(2,3)*BB(3,2)*BB(4,1) - BB(2,1)*BB(3,2)*BB(4,3) 
     * - BB(2,2)*BB(3,3)*BB(4,1) - BB(2,3)*BB(3,1)*BB(4,2) ) * (1/DetBB)
       AA(4,2) = (BB(1,1)*BB(3,2)*BB(4,3) + BB(1,2)*BB(3,3)*BB(4,1) 
     * + BB(1,3)*BB(3,1)*BB(4,2) - BB(1,1)*BB(3,3)*BB(4,2) 
     * - BB(1,2)*BB(3,1)*BB(4,3) - BB(1,3)*BB(3,2)*BB(4,1) ) * (1/DetBB)
      AA(4,3) = (BB(1,1)*BB(2,3)*BB(4,2) + BB(1,2)*BB(2,1)*BB(4,3) 
     * + BB(1,3)*BB(2,2)*BB(4,3) - BB(1,1)*BB(2,2)*BB(4,3) 
     * - BB(1,2)*BB(2,3)*BB(4,1) - BB(1,3)*BB(2,1)*BB(4,2) ) * (1/DetBB)
      AA(4,4) = (BB(1,1)*BB(2,2)*BB(3,3) + BB(1,2)*BB(2,3)*BB(3,1) 
     * + BB(1,3)*BB(2,1)*BB(3,2) - BB(1,1)*BB(2,3)*BB(3,2) 
     * - BB(1,2)*BB(2,1)*BB(3,3) - BB(1,3)*BB(2,2)*BB(3,1) ) * (1/DetBB)





      endif
      
       
      RETURN
      END SUBROUTINE MatrixInv

Ccccccccccccccccccccccccc Matrix determinant ccccccccccccccccccccccccccc
      
      SUBROUTINE MatrixDet(AA, Det, MCRD)
      
      
      integer MCRD
      double precision AA, Det
            
      dimension AA(MCRD, MCRD)
      
      if (MCRD==2) then
          Det=(AA(1,1)*AA(2,2))-(AA(2,1)*AA(1,2))
      else if (MCRD==3) then
          Det=AA(1,2)*AA(2,3)*AA(3,1)
     1        -AA(1,3)*AA(2,2)*AA(3,1)
     2        +AA(1,3)*AA(2,1)*AA(3,2)
     3        -AA(1,1)*AA(2,3)*AA(3,2)
     4        +AA(1,1)*AA(2,2)*AA(3,3)
     5        -AA(1,2)*AA(2,1)*AA(3,3)
      else if (MCRD==4) then
          Det=AA(1,1)*AA(2,2)*AA(3,3)*AA(4,4)
     *        +AA(1,1)*AA(2,3)*AA(3,4)*AA(4,2)
     *        + AA(1,1)*AA(2,4)*AA(3,2)*AA(4,3)
     *        + AA(1,2)*AA(2,1)*AA(3,4)*AA(4,3)
     *        + AA(1,2)*AA(2,3)*AA(3,1)*AA(4,4)
     *        + AA(1,2)*AA(2,4)*AA(3,3)*AA(4,1)
     *        + AA(1,3)*AA(2,1)*AA(3,2)*AA(4,4)
     *        + AA(1,3)*AA(2,2)*AA(3,4)*AA(4,1)
     *        + AA(1,3)*AA(2,4)*AA(3,1)*AA(4,2)
     *        + AA(1,4)*AA(2,1)*AA(3,3)*AA(4,1)
     *        + AA(1,4)*AA(2,2)*AA(3,1)*AA(4,3)
     *        + AA(1,4)*AA(2,3)*AA(3,2)*AA(4,1)
     *        - AA(1,1)*AA(2,2)*AA(3,4)*AA(4,3)
     *        - AA(1,1)*AA(2,3)*AA(3,2)*AA(4,4)
     *        - AA(1,1)*AA(2,4)*AA(3,3)*AA(4,2)
     *        - AA(1,2)*AA(2,1)*AA(3,3)*AA(4,4)
     *        - AA(1,2)*AA(2,3)*AA(3,4)*AA(4,1)
     *        - AA(1,2)*AA(2,4)*AA(3,1)*AA(4,3)
     *        - AA(1,3)*AA(2,1)*AA(3,4)*AA(4,2)
     *        - AA(1,3)*AA(2,2)*AA(3,1)*AA(4,4)
     *        - AA(1,3)*AA(2,4)*AA(3,2)*AA(4,1)
     *        - AA(1,4)*AA(2,1)*AA(3,2)*AA(4,3)
     *        - AA(1,4)*AA(2,2)*AA(3,3)*AA(4,1)
     *        - AA(1,4)*AA(2,3)*AA(3,1)*AA(4,2)
      endif
      
      
      RETURN
      END SUBROUTINE MatrixDet      
      
      
cccccccccccccccccccccccccc Lecteur de cl� cccccccccc
c
c     KCle: cl� � lire (Ex:1)
c     KNumFace: Premier indice, num�ro de la face
c     KTypeDlaod:Deuxi�me indice, type de chargement surfacique

      SUBROUTINE LectCle (KCle,KNumFace,KTypeDload)
      character*2 s
      integer KCle, KNumFace, KTypeDload
      
      write (s,FMT='(I2)') KCle
      read(s(1:1),FMT='(I1)')   KNumFace
      read(s(2:2),FMT='(I1)')   KTypeDload
      RETURN
      END SUBROUTINE LectCle
      
      
      SUBROUTINE GetLoadDescription(JDLType,KNumFace,KTypeDload)
      Implicit None
      Integer, intent(in) :: JDLType
      Integer, intent(out):: KNumFace,KTypeDload

      If (JDLTYPE<10) then
         ! concentrated load
         KNumFace = -1
         KTypeDload = JDLType
      elseif (JDLTYPE>9 .AND. JDLTYPE<100) then
         ! surface load
         call LectCle(JDLType,KNumFace,KTypeDload)
      else
         ! body load
         KNumFace = 0
         KTypeDload = JDLType
      Endif
      
      END SUBROUTINE GetLoadDescription

Cccccccccccc Contruction Matrice de Passage Globlal 2 Local cccccccccccc
      
      Subroutine transposeMat(M,Mt,n)
      Implicit None
      Integer, intent(in) :: n
      Double precision, dimension(n,n), intent(in) :: M
      Double precision, dimension(n,n), intent(out):: Mt
      
      Integer :: i,j
      
      Mt(:,:) = 0.0D0
      do i = 1,n
         do j = 1,n
            Mt(j,i) = M(i,j)
         enddo
      enddo
      
      End subroutine transposeMat



      Subroutine dot(u,v,w)
      Implicit None
      double precision, dimension(3), intent(in) :: u,v
      double precision, intent(out):: w
      
      w = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      
      End subroutine dot
      
      Function DOTPROD(u,v) result(w)
      Implicit None
      double precision, dimension(3), intent(in) :: u,v
      double precision :: w
      w = SUM(u(:)*v(:))
      End function DOTPROD
      
      
      Subroutine norm(u,n, res)
      Implicit None
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: u
      double precision, intent(out):: res
      
      integer :: i

      res = 0.0d0
      do i = 1,n
         res = res + u(i)*u(i)
      enddo
      res = sqrt(res)
      
      End subroutine norm
      
      
      
      
      Subroutine cross(u,v,w)
      Implicit None
      double precision, dimension(3), intent(in) :: u,v
      double precision, dimension(3), intent(out):: w
      
      w(1) = u(2)*v(3) - u(3)*v(2) 
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      
      End subroutine cross
      
      Function CROSSPROD(u,v) result(w)
      Implicit None
      double precision, dimension(3), intent(in) :: u,v
      double precision, dimension(3) :: w
      w(1) = u(2)*v(3) - u(3)*v(2) 
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      End function CROSSPROD
      
      Function SCATRIPLEPROD(u,v,w) result(scalar)
      Implicit None
      double precision, dimension(3), intent(in) :: u,v,w
      double precision :: scalar
      double precision, dimension(3) :: tmp
      CALL CROSS(u,v,tmp)
      scalar = SUM(tmp(:)*w(:))
      End function SCATRIPLEPROD
      
      
      Subroutine getTransformationMatrix(P1,P2,P3,Mat)
      Implicit None
      
!     Inputs
      double precision, dimension(3), intent(in) :: P1,P2,P3
      
!     Outputs
      double precision, dimension(3,3), intent(out) :: Mat
      
!     Local variables
      double precision, dimension(3) :: x,y,z
      double precision :: norm
      
      x(:) = P2(:) - P1(:)
      call dot(x,x,norm)
      x(:) = x(:)/sqrt(norm)
      
      call cross(x, P3-P1, z)
      call dot(z,z,norm)
      z(:) = z(:)/sqrt(norm)
      
      call cross(z,x, y)
      
      Mat(1,:) = x(:)
      Mat(2,:) = y(:)
      Mat(3,:) = z(:)
      
      End subroutine getTransformationMatrix
      
      
      
      Subroutine matA_tilde(u, A)
      implicit none
      double precision, intent(in) :: u
      double precision, intent(out):: A
      dimension u(3), A(3,3)
      
      A(:,1) = (/  0.d0, -u(3),  u(2) /)
      A(:,2) = (/  u(3),  0.d0, -u(1) /)
      A(:,3) = (/ -u(2),  u(1),  0.d0 /)
      End subroutine matA_tilde
      
      Subroutine matA_tildeT(u, At)
      implicit none
      double precision, intent(in) :: u
      double precision, intent(out):: At
      dimension u(3), At(3,3)
      
      At(:,1) = (/  0.d0,  u(3), -u(2) /)
      At(:,2) = (/ -u(3),  0.d0,  u(1) /)
      At(:,3) = (/  u(2), -u(1),  0.d0 /)
      End subroutine matA_tildeT
      
      
      
      Subroutine matA_bar(u,v, A)
      implicit none
      double precision, intent(in) :: u,v
      double precision, intent(out):: A
      dimension u(3), v(3), A(3,3)
      integer j
      
      Do j = 1,3
         A(:,j) = u(:)*v(j)
      Enddo
      End subroutine matA_bar
      
      Subroutine matA_barT(u,v, At)
      implicit none
      double precision, intent(in) :: u,v
      double precision, intent(out):: At
      dimension u(3), v(3), At(3,3)
      integer j
      
      Do j = 1,3
         At(:,j) = v(:)*u(j)
      Enddo
      End subroutine matA_barT
      
      
      
      Subroutine matA_dbletidle(u,v, A)
      implicit none
      double precision, intent(in) :: u,v
      double precision, intent(out):: A
      dimension u(3), v(3), A(3,3)
      
      A(:,1) = (/ u(2)*v(2)+u(3)*v(3), -u(1)*v(2), -u(1)*v(3) /)
      A(:,2) = (/ -u(2)*v(1), u(1)*v(1)+u(3)*v(3), -u(2)*v(3) /)
      A(:,3) = (/ -u(3)*v(1), -u(3)*v(2), u(1)*v(1)+u(2)*v(2) /)
      End subroutine matA_dbletidle
      
      Subroutine matA_dbletidleT(u,v, At)
      implicit none
      double precision, intent(in) :: u,v
      double precision, intent(out):: At
      dimension u(3), v(3), At(3,3)
      
      At(:,1) = (/ u(2)*v(2)+u(3)*v(3), -u(2)*v(1), -u(3)*v(1) /)
      At(:,2) = (/ -u(1)*v(2), u(1)*v(1)+u(3)*v(3), -u(3)*v(2) /)
      At(:,3) = (/ -u(1)*v(3), -u(2)*v(3), u(1)*v(1)+u(2)*v(2) /)
      End subroutine matA_dbletidleT
      
