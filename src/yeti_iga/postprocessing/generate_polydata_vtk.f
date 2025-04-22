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

C     Cette routine permet de generer un fichier vtk pour visualiser le
c     polygone de controle d'une geometrie NURBS avec Paraview


      Subroutine generate_vtk(case_name, output_path, COORDS3D,
     &        nb_cp_dir, nb_cp)

C     ----------
c     Input Data
      Character(len=*), intent(in) :: case_name
      character(len=*), intent(in) :: output_path
      Integer, dimension(3),  intent(in) :: nb_cp_dir
      Integer, intent(in) :: nb_cp

      Double precision, intent(in) :: COORDS3D
      Dimension COORDS3D(3,nb_cp)



C     ---------------
c     Local Variables
      Integer :: num_cp,num_cp_xi,num_cp_eta,num_cp_zeta,i,ind_cp0
      Integer :: nb_lines, nb_data4lines
c      Integer, dimension(3) :: nb_cp_dir

      Character(len=500) :: ind_cpthisline, fmt


C     Format pour ecrire les coordonnes des noeuds
 10   Format(3(X,E12.4))




C     --
c     Ecriture du fichier
      open(90, file= output_path// '/' //case_name // '.vtk',
     &      form='formatted')

      write(90,'(A26)') '# vtk DataFile Version 3.0'
      write(90,'(A10)') 'vtk output'
      write(90,'(A5)')  'ASCII'
      write(90,'(A16)') 'DATASET POLYDATA'
      write(90,'(A6,X,I8,X,A6)') 'POINTS',nb_cp,'double'

C     Coordonnes des points de controle
      Do num_cp = 1,nb_cp
         write(90,10) COORDS3D(1,num_cp), COORDS3D(2,num_cp),
     &        COORDS3D(3,num_cp)
      Enddo


C     Definition des lignes reliant les points de controle
      nb_lines = nb_cp_dir(2)*nb_cp_dir(3)
     &     + nb_cp_dir(1)*nb_cp_dir(3)
     &     + nb_cp_dir(1)*nb_cp_dir(2)
      nb_data4lines = nb_lines + 3*nb_cp

      write(90,'(A5,X,I5,X,I6)') 'LINES', nb_lines, nb_data4lines


c     Ligne xi
      Do num_cp_zeta = 1,nb_cp_dir(3)
         Do num_cp_eta = 1, nb_cp_dir(2)
            ind_cp0 = (num_cp_zeta-1)*nb_cp_dir(1)*nb_cp_dir(2)
     &           + (num_cp_eta-1)*nb_cp_dir(1)

            write(ind_cpthisline,'(I5)') nb_cp_dir(1)
            Do i = 1,nb_cp_dir(1)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline, ind_cp0 + i-1
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(1)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


c     Ligne eta
      Do num_cp_zeta = 1,nb_cp_dir(3)
         Do num_cp_xi = 1,nb_cp_dir(1)
            ind_cp0 = (num_cp_zeta-1)*nb_cp_dir(1)*nb_cp_dir(2)
     &           + num_cp_xi - 1

            write(ind_cpthisline,'(I5)') nb_cp_dir(2)
            Do i = 1,nb_cp_dir(2)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline,
     &              ind_cp0+(i-1)*nb_cp_dir(1)
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(2)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


c     Ligne zeta
      Do num_cp_eta = 1,nb_cp_dir(2)
         Do num_cp_xi = 1,nb_cp_dir(1)
            ind_cp0 = (num_cp_eta-1)*nb_cp_dir(1)
     &           + num_cp_xi-1

            write(ind_cpthisline,'(I5)') nb_cp_dir(3)
            Do i = 1,nb_cp_dir(3)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline,
     &              ind_cp0+(i-1)*nb_cp_dir(1)*nb_cp_dir(2)
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(3)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


      close(90)

C     --

      End subroutine generate_vtk















C     Cette routine permet de generer un fichier vtk pour visualiser le
c     polygone de controle d'une geometrie NURBS avec Paraview
c     --> cas avec solution


      Subroutine generate_vtk_wSOL(case_name, output_path, COORDS3D,
     &     SOL3D,nb_cp_dir,nb_cp)

C     ----------
c     Input Data
      Character(len=*), intent(in) :: case_name
      character(len=*), intent(in) :: output_path
      Integer, dimension(3),  intent(in) :: nb_cp_dir
      Integer, intent(in) :: nb_cp

      Double precision, intent(in) :: COORDS3D, SOL3D
      Dimension COORDS3D(3,nb_cp), SOL3D(3,nb_cp)



C     ---------------
c     Local Variables
      Integer :: num_cp,num_cp_xi,num_cp_eta,num_cp_zeta,i,ind_cp0
      Integer :: nb_lines, nb_data4lines
c      Integer, dimension(3) :: nb_cp_dir

      Character(len=500) :: ind_cpthisline, fmt


C     Format pour ecrire les coordonnes des noeuds
 10   Format(3(X,E12.4))




C     --
c     Ecriture du fichier
      open(90, file=output_path// '/' // case_name // '.vtk',
     &     form='formatted')

      write(90,'(A26)') '# vtk DataFile Version 3.0'
      write(90,'(A10)') 'vtk output'
      write(90,'(A5)')  'ASCII'
      write(90,'(A16)') 'DATASET POLYDATA'
      write(90,'(A6,X,I8,X,A6)') 'POINTS',nb_cp,'double'

C     Coordonnes des points de controle
      Do num_cp = 1,nb_cp
         write(90,10) COORDS3D(1,num_cp), COORDS3D(2,num_cp),
     &        COORDS3D(3,num_cp)
      Enddo


C     Definition des lignes reliant les points de controle
      nb_lines = nb_cp_dir(2)*nb_cp_dir(3)
     &     + nb_cp_dir(1)*nb_cp_dir(3)
     &     + nb_cp_dir(1)*nb_cp_dir(2)
      nb_data4lines = nb_lines + 3*nb_cp

      write(90,'(A5,X,I5,X,I6)') 'LINES', nb_lines, nb_data4lines


c     Ligne xi
      Do num_cp_zeta = 1,nb_cp_dir(3)
         Do num_cp_eta = 1, nb_cp_dir(2)
            ind_cp0 = (num_cp_zeta-1)*nb_cp_dir(1)*nb_cp_dir(2)
     &           + (num_cp_eta-1)*nb_cp_dir(1)

            write(ind_cpthisline,'(I5)') nb_cp_dir(1)
            Do i = 1,nb_cp_dir(1)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline, ind_cp0 + i-1
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(1)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


c     Ligne eta
      Do num_cp_zeta = 1,nb_cp_dir(3)
         Do num_cp_xi = 1,nb_cp_dir(1)
            ind_cp0 = (num_cp_zeta-1)*nb_cp_dir(1)*nb_cp_dir(2)
     &           + num_cp_xi - 1

            write(ind_cpthisline,'(I5)') nb_cp_dir(2)
            Do i = 1,nb_cp_dir(2)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline,
     &              ind_cp0+(i-1)*nb_cp_dir(1)
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(2)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


c     Ligne zeta
      Do num_cp_eta = 1,nb_cp_dir(2)
         Do num_cp_xi = 1,nb_cp_dir(1)
            ind_cp0 = (num_cp_eta-1)*nb_cp_dir(1)
     &           + num_cp_xi-1

            write(ind_cpthisline,'(I5)') nb_cp_dir(3)
            Do i = 1,nb_cp_dir(3)
               write(fmt,'(A2,I5,A6)') '(A',i*6,',X,I5)'
               write(ind_cpthisline, fmt) ind_cpthisline,
     &              ind_cp0+(i-1)*nb_cp_dir(1)*nb_cp_dir(2)
            Enddo
            write(fmt,'(A2,I5,A2)') '(A',(nb_cp_dir(3)+1)*6,')'
            write(90,fmt) ind_cpthisline
         Enddo
      Enddo


c     Solution (champ vecteur)
      write(90,'(A10,X,I8)') 'POINT_DATA',nb_cp
      write(90,'(A22)') 'VECTORS solution float'
      Do num_cp = 1,nb_cp
         write(90,10) SOL3D(1,num_cp), SOL3D(2,num_cp), SOL3D(3,num_cp)
      Enddo



      close(90)

C     --

      End subroutine generate_vtk_wSOL

