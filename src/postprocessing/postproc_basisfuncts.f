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

c     POST-PROCESSING : display NURBS basis functions and derivatives
c      - multiple .txt files (per patch, per order of derivation)
            
      subroutine generateNurbsBasisTXT(FILENAME,ActivePatch,DervOrder,
     1     nb_pts, COORDS3D,IEN,nb_elem_patch,Nkv,Ukv,Nijk,weight,
     2     Jpqr,ELT_TYPE,PROPS,JPROPS,MATERIAL_PROPERTIES,TENSOR,MCRD,
     3     NBINT,nnode,nb_patch,nb_elem,nb_cp)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
c     Declaration des Variables ........................................
      
c     Input arguments :
c     ---------------
      
!     Geometry NURBS
      Integer,          intent(in) :: nb_cp
      Double precision, intent(in) :: COORDS3D
      dimension COORDS3D(3,nb_cp)
      
      Double precision, intent(in) :: Ukv, weight
      Integer,          intent(in) :: Nkv, Jpqr, Nijk
      dimension Nkv(3,nb_patch), Jpqr(3,nb_patch), Nijk(3,nb_elem),
     &     Ukv(:),weight(:)
      
      
!     Patches and Elements
      Character(len=*), intent(in) :: TENSOR, ELT_TYPE
      Double precision, intent(in) :: MATERIAL_PROPERTIES, PROPS
      Integer,          intent(in) :: MCRD,NNODE,nb_patch,nb_elem,NBINT,
     &     IEN,nb_elem_patch,JPROPS
      dimension MATERIAL_PROPERTIES(2,nb_patch),
     &     PROPS(:),
     &     NNODE(nb_patch),
     &     IEN(:),
     &     nb_elem_patch(nb_patch),
     &     JPROPS(nb_patch),
     &     NBINT(nb_patch)
      
      
!     Output INFOS
      Character(len=*), intent(in) :: FILENAME
      Integer,          intent(in) :: nb_pts,ActivePatch,DervOrder
      dimension ActivePatch(nb_patch)
      
      
c     Local variables :
c     ---------------
      
!     Extract infos
      Integer          :: sctr
      dimension sctr(MAXVAL(NNODE))
      Double precision :: COORDS_elem
      dimension COORDS_elem(MCRD,MAXVAL(NNODE))
      
!     Compute quantities
      Integer          :: Ni,cp_dir,dir
      Double precision :: h,XI,R,dRdxi,ddRddxi,Ukv_dir
      dimension R(MAXVAL(Jpqr)+1), dRdxi(MAXVAL(Jpqr)+1),
     &     ddRddxi(MAXVAL(Jpqr)+1), Ukv_dir(MAXVAL(Nkv))
      
      Integer :: i,j,kk,count,numPatch,num_elem,nb_elem_dir
      
      Character(len=50) :: formatdata,filenamePatch,filenamePatchDerv,
     &     filenamePatchDDerv

      
C     Fin declaration des variables ....................................
c     
c     
c     
c     
c     Start  ...........................................................
      
      
c     Retour ecran
      write(*,*)'Post processsing 1D B-Spline basis functions ...'
      
      count = 1
      Do NumPatch = 1,nb_patch
         
         If (ActivePatch(NumPatch) == 1) then

         CALL extractNurbsPatchMechInfos(NumPatch,IEN,PROPS,JPROPS,
     &        NNODE,nb_elem_patch,ELT_TYPE,TENSOR)
         CALL extractNurbsPatchGeoInfos(NumPatch, Nkv,Jpqr,Nijk,Ukv,
     &        weight,nb_elem_patch)
         
         Do dir = 1,dim_patch
            
            cp_dir = Jpqr_patch(dir)+1
            write(filenamePatch,'(A,A,A,I0,A,I0,A)')
     &           'results/',FILENAME,'_patch',numPatch,'_dir',dir,'.txt'
            Open(90,file=filenamePatch,form='formatted')

            If (DervOrder>0) then
               write(filenamePatchDerv,'(A,A,A,I0,A,I0,A)')
     &              'results/',FILENAME,'Der1_patch',numPatch,'_dir',
     &              dir,'.txt'
               Open(91,file=filenamePatchDerv,form='formatted')
            Endif
            
            If (DervOrder>1) then
               write(filenamePatchDDerv,'(A,A,A,I0,A,I0,A)')
     &              'results/',FILENAME,'Der2_patch',numPatch,'_dir',
     &              dir,'.txt'
               Open(92,file=filenamePatchDDerv,form='formatted')
            Endif
            
            write(formatdata,'("(1E14.6,"I0,"(2X,1E14.6))")') cp_dir
            
            
            
            nb_elem_dir = Nkv_patch(dir) - 2*Jpqr_patch(dir) - 1
            Ni = Jpqr_patch(dir)+1
            
            If (dir==1) Ukv_dir(:Nkv_patch(dir))=Ukv1_patch(:)
            If (dir==2) Ukv_dir(:Nkv_patch(dir))=Ukv2_patch(:)
            If (dir==3) Ukv_dir(:Nkv_patch(dir))=Ukv3_patch(:)
            
            kk = 0
            j  = 1
            Do num_elem = 1,nb_elem_dir
               
               h  = (Ukv_dir(Ni+1) - Ukv_dir(Ni))/dble(nb_pts-1)
               If (h > zero) then
               
               XI = Ukv_dir(Ni)
               Do i = 1,nb_pts
!     - evaluate basis functions
                  
                  If (DervOrder>1) then
                     CALL dersbasisfuns2(Ni,Jpqr_patch(dir),
     &                    Nkv_patch(dir),XI,Ukv_dir(:Nkv_patch(dir)),
     &                    R(:cp_dir),dRdxi(:cp_dir),ddRddxi(:cp_dir))
                  else
                     CALL dersbasisfuns(Ni,Jpqr_patch(dir),
     &                    Nkv_patch(dir),XI,Ukv_dir(:Nkv_patch(dir)),
     &                    R(:cp_dir),dRdxi(:cp_dir))
                  Endif
               
!     - write TXT files
                  write(90,fmt=formatdata) XI,R(cp_dir+1-j:cp_dir),
     &                 R(1:cp_dir-j)
                  If (DervOrder>0) then
                  write(91,fmt=formatdata) XI,dRdxi(cp_dir+1-j:cp_dir),
     &                 dRdxi(1:cp_dir-j)
                  Endif
                  If (DervOrder>1) then
                  write(92,fmt=formatdata)XI,ddRddxi(cp_dir+1-j:cp_dir),
     &                 ddRddxi(1:cp_dir-j)
                  Endif
                  
                  XI = XI + h
               Enddo
               Endif

               Ni = Ni+1
               j  = j+1
               if (j==cp_dir+1) j=1

            Enddo ! loop elem
            
            close(90)
            print*,' generated file: ',filenamePatch
            If (DervOrder>0) then
               close(91)
               print*,' generated file: ',filenamePatchDerv
            Endif
            If (DervOrder>1) then
               close(92)
               print*,' generated file: ',filenamePatchDDerv
            Endif
            
         Enddo ! loop dir
         
         call finalizeNurbsPatch()
         
         Endif ! active patch
         
      Enddo ! end loop on patch
      
c     ..................................................................
            
      end subroutine generateNurbsBasisTXT
