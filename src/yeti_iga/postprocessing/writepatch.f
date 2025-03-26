!! Copyright 2018-2019 Thibaut Hirschler

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


      Subroutine writePatch(file,FieldOutput_flag,svars,nb_node,nb_elem,
     1     dim,MCRD,nsvint,ntens,nb_elemVTU,nb_xi,nb_eta,nb_zeta,
     2     nb_vertice)
      
      use parameters
      use nurbspatch
      
      Implicit none
      
C     ------------------------------------------------------------------
      
C     Input arguments :
c     ---------------
!     Data to write
      Integer,          intent(in) :: file
      Logical,          intent(in) :: FieldOutput_flag
      Double precision, intent(in) :: svars
      dimension FieldOutput_flag(3)
      dimension svars(nsvint*nb_vertice,nb_elem)
      
!     Size infos
      Integer, intent(in) :: nb_node,nb_elem,dim,MCRD,nsvint,ntens
      Integer, intent(in) :: nb_elemVTU,nb_xi,nb_eta,nb_zeta,nb_vertice

      
C     Local variables :
c     ---------------
      Double precision :: VM,stress,strain
      Integer          :: i,j,i_xi,i_eta,i_zeta,comp,compt,offset,numel
      dimension stress(6),strain(6)
      
C     ------------------------------------------------------------------
c     Extract data
      




c     ------------------------------------------------------------------
c     Write

c     Start piece
      write(file,*)'<Piece  NumberOfPoints="  ', nb_node,
     1     '"  NumberOfCells=" ', nb_elem*nb_elemVTU, '">'  
      
      
c     Writing DATA points
      write(file,*)'<Points>'
      write(file,*)'<DataArray  type="Float64"' 
     1     //' NumberOfComponents="3"  format="ascii" >'
      Do numel = 1,nb_elem
         offset = 1
         Do i  = 1,nb_vertice
            write(file,*)svars(offset:offset+2,numel)
            offset = offset+nsvint
         Enddo
      Enddo
      write(file,*)'</DataArray>'
      write(file,*)'</Points>'
      
      
C     ---
c     Writing cell connectivity
      
      write(file,*)'<Cells>'
      write(file,*)'<DataArray  type="Int32"  Name="connectivity"'
     1     //'  format="ascii">'
      
      comp = 0
      if (dim==2) then
         do i = 1,nb_elem
            do i_eta = 1,nb_eta-1
               do i_xi = 1,nb_xi-1
                  comp = (i-1)*nb_vertice
                  comp = comp + (i_eta-1)*nb_xi + i_xi - 1
                  write(file,*)comp,comp+1,comp+nb_xi+1,comp+nb_xi
               enddo
            enddo
         enddo
      else if (dim==3) then
         do i = 1,nb_elem
            do i_zeta = 1,nb_zeta-1
               do i_eta = 1,nb_eta-1
                  do i_xi = 1,nb_xi-1
                     comp = (i-1)*nb_vertice
                     comp = comp + (i_zeta-1)*nb_xi*nb_eta
     &                    + (i_eta-1)*nb_xi +  i_xi - 1 
                     compt= comp + nb_xi*nb_eta
                     write(file,*)comp,comp+1,comp+nb_xi+1,comp+nb_xi,
     &                    compt,compt+1,compt+nb_xi+1,compt+nb_xi
                  enddo
               enddo
            enddo
         enddo
      endif
      write(file,*)'</DataArray>'
      
      
c     Writing cell offsets
      write(file,*)'<DataArray  type="Int32"  Name="offsets"'
     1     // '  format="ascii"> '
      offset = 0
      do i = 1,nb_elem*nb_elemVTU
         offset = offset + 2**dim
         write(file,*)offset
      enddo
      write(file,*)'</DataArray>'
      
      
c     writing cell types
      write(file,*)'<DataArray  type="UInt8"  Name="types"'
     1     // '  format="ascii">'
      do i = 1,nb_elem*nb_elemVTU
         if (dim==2) then
           write(file,*)'9'
        else if (dim==3) then
           write(file,*)'12'
        endif 
      enddo
      write(file,*) '</DataArray>'
      write(file,*) '</Cells> '
      
c     End cell informations
c     ---
      
      
      
c     ---
c     Ecriture variables d'interets
      write(file,*)'<PointData>'
      
      
c     Displacement field
      If (FieldOutput_flag(1)) then
      write(file,*)'<DataArray  type="Float64"'   
     1     // ' Name="disp" NumberOfComponents="3" format="ascii">'
      Do numel = 1,nb_elem
         offset= 4
         Do i  = 1,nb_vertice
            write(file,*)svars(offset:offset+2,numel)
            offset = offset+nsvint
         Enddo
      Enddo
      write(file,*)'</DataArray>'
      Endif
      
c     Stress
      If (FieldOutput_flag(2)) then
      compt = 4
      if (FieldOutput_flag(1)) compt = compt+3
      if (MCRD==2) then
         write(file,*)'<DataArray  type="Float64"'   
     1        // ' Name="stress" NumberOfComponents="3" format="ascii">'
         Do numel = 1,nb_elem
            offset= compt
            Do i  = 1,nb_vertice
               stress(:ntens) = svars(offset:offset+ntens-1,numel)
               write(file,*)stress(1),stress(2),stress(4)
               offset = offset+nsvint
            Enddo
         Enddo
      elseif (MCRD==3) then
         write(file,*)'<DataArray  type="Float64"'   
     1        // ' Name="stress" NumberOfComponents="6" format="ascii">'
         Do numel = 1,nb_elem
            offset= compt
            Do i  = 1,nb_vertice
               stress(:) = svars(offset:offset+ntens-1,numel)
               write(file,*)stress(1),stress(2),stress(3),stress(4),
     &              stress(5),stress(6)
               offset = offset+nsvint
            Enddo
         Enddo
      endif
      write(file,*)'</DataArray>'
      
c     Von Mises
      write(file,*)'<DataArray  type="Float64"'   
     1     // ' Name="VM" NumberOfComponents="1" format="ascii">'
      Do numel = 1,nb_elem
         offset= compt
         Do i  = 1,nb_vertice
            stress(:ntens) = svars(offset:offset+ntens-1,numel)
            if (dim==2) then    ! plane stress
               if (MCRD==3) stress(4) = stress(3)
               stress(3)  = zero
               stress(5:) = zero
            endif
            VM = (stress(1)-stress(2))**2
     &         + (stress(1)-stress(3))**2
     &         + (stress(2)-stress(3))**2
     &         + 6.d0*(stress(4)**2 + stress(5)**2 + stress(6)**2)
            write(file,*)SQRT(0.5d0*VM)
            offset = offset+nsvint
         Enddo
      Enddo
      write(file,*)'</DataArray>'
      Endif      
      
c     Strain
      If (FieldOutput_flag(3)) then
      compt = 4
      if (FieldOutput_flag(1)) compt = compt+3
      if (FieldOutput_flag(2)) compt = compt+ntens
      if (MCRD==2) then
         write(file,*)'<DataArray  type="Float64"'   
     1        // ' Name="strain" NumberOfComponents="3" format="ascii">'
         Do numel = 1,nb_elem
            offset= compt
            Do i  = 1,nb_vertice
               strain(:ntens) = svars(offset:offset+ntens-1,numel)
               write(file,*)strain(1),strain(2),strain(4)
               offset = offset+nsvint
            Enddo
         Enddo
         write(file,*)'</DataArray>'
         if (ELT_TYPE_patch == "U98") then
            !! Write strain gradient
            !! WARNING : only 2D case is implemented
            write(file,*)'<DataArray  type="Float64"'   
     1   // ' Name="dstraindx" NumberOfComponents="6" format="ascii">'
            do numel = 1, nb_elem
                offset = offset + ntens
                do i = 1, nb_vertice
                    !! Pas de variable temporaire
                    !!dstraindx(:6) = svars(offset:offset+6-1,numel)
                    !!write(file,*) (svars(offset+j, numel), j=0, 6-1)
                    !! **************************************
                    !! FAIRE L'ECRITURE DES VALEURS DE SVARS
                    !! *************************************
                    write(file,*) (0., j=0, 6-1)
                    !! write(*,*) nsvint
                    offset = offset+nsvint
                enddo
            enddo
            write(file,*)'</DataArray>'
         endif
      elseif (MCRD==3) then
         write(file,*)'<DataArray  type="Float64"'   
     1        // ' Name="strain" NumberOfComponents="6" format="ascii">'
         Do numel = 1,nb_elem
            offset= compt
            Do i  = 1,nb_vertice
               strain(:) = svars(offset:offset+ntens-1,numel)
               write(file,*)strain(1),strain(2),strain(3),strain(4),
     &              strain(5),strain(6)
               offset = offset+nsvint
            Enddo
         Enddo
         write(file,*)'</DataArray>'
      endif
      Endif
      
      write(file,*)'</PointData>'
c     Fin ecriture donnees
c     ---
      
      
c     Writing End of file
      write(file,*)'</Piece>'
      
      
      
      End subroutine writePatch
