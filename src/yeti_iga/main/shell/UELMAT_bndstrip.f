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

c     --
c     Contribution a la matrice de rigite globale pour l'ajout d'element
c     de flexion.
c     --
      
      SUBROUTINE UELMATbndstrip(MCRD,NNODE,JELEM,NBINT,NumPatch,COORDS,
     1     MATERIAL_PROPERTIES,PROPS,JPROPS,AMATRX)
      
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: MCRD,NNODE,JELEM,NumPatch,NBINT,JPROPS
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(MCRD,NNODE)
      dimension MATERIAL_PROPERTIES(2)
      dimension PROPS(10)
      
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: AMATRX
      dimension AMATRX(MCRD*NNODE,MCRD*NNODE)
            
      
c     Local variables :
c     ---------------
      Double precision zero, one, two
      parameter (zero=0.d0, one=1.d0, two=2.d0)
      
!     For gauss points
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(3,NBINT)
      
!     For nurbs basis functions
      Double precision :: R, dRdxi, ddRddxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI, dAI1dxi, dAI2dxi, AAE
      dimension AI(3,3), dAI1dxi(2,3), dAI2dxi(2,3), AAE(2,2)
      
!     For material matrix
      Double precision :: E, nu, h, matH, coef
      dimension matH(3,3)
      
!     For stiffness matrices and load vectors
      Double precision :: bendingStiff,dvol,Area
      dimension bendingStiff(MCRD*NNODE,MCRD*NNODE)
      
!     Other
      Integer ntens
      Integer n, k1, k2, numPC, numi
      
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
      
c     Stiffness matrix and force vector are initialized to zero
      AMATRX(:,:) = zero
      
c     Material behaviour
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1) * 10.d0**4.d0
      nu= zero
      coef = E/(one-nu**two)
      matH = zero
c
c     ..................................................................
c
C     Computation :
c     -----------
      
c     Loop on integration points on main surface
      do n = 1,NBINT
         
c     Computing NURBS basis functions and derivatives
         R = zero
         dRdxi = zero
         ddRddxi = zero
         DetJac = zero
         call bsplinebasis(R,dRdxi,ddRddxi,DetJac,NNODE,
     &        GaussPdsCoord(2:,n),JELEM,NumPatch)
         
         
c     Computing Curvilinear Coordinate objects
         call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,R,dRdxi,ddRddxi,MCRD,
     &        NNODE,COORDS)
         
         
c     Computing material matrix
         matH(:,:) = zero
         matH(1,1) = coef*AAE(1,1)*AAE(1,1)
         
         
c     Computing stiffness matrix
         bendingStiff = zero
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         call usfbnd_shell(NNODE,MCRD,NNODE*MCRD,matH,AI,dAI1dxi,
     &        dAI2dxi,dRdxi,ddRddxi,bendingStiff)
         
c     
c     Assembling RHS and AMATRIX
         call SurfElem(AI(1,:), AI(2,:), Area)
         dvol = GaussPdsCoord(1,n)*DetJac*Area
         do k1 = 1,MCRD*NNODE
            do k2 = 1,MCRD*NNODE
               amatrx(k1,k2) = amatrx(k1,k2) + bendingStiff(k1,k2)*dvol
            enddo
         enddo

      enddo
      
c     End of the loop on integration points on main surf
c
c     ..................................................................
      
      
      End SUBROUTINE UELMATbndstrip
