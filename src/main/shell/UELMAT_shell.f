!! Copyright 2017-2018 Thibaut Hirschler
!! Copyright 2019 Arnaud Duval

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
c     Construction de la matrice de rigite et du second membre 
c     elementaire pour la formulation coque type Kirchhoff-Love
c     --
      
      SUBROUTINE UELMAT3(NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,TENSOR,
     1     MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,indDLoad,
     2     load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX)
      
      use parameters
    
      implicit none      
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
      Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
      Double precision, intent(in) :: ADLMAG
      dimension ADLMAG(nb_load),
     &     load_target_nbelem(nb_load),JDLType(nb_load)
      dimension indDLoad(SUM(load_target_nbelem))
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: RHS, AMATRX
      dimension RHS(NDOFEL), AMATRX(NDOFEL,NDOFEL)
      
      
c     Local variables :
c     ---------------
      
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
      Double precision :: membraneStiff,bendingStiff,dvol,Area,normV,
     &     VectNorm
      dimension VectNorm(MCRD), membraneStiff(NDOFEL,NDOFEL),
     &     bendingStiff(NDOFEL,NDOFEL)
      
!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,kk,numPC,Numi,dof, KTypeDload, KNumFace
      
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
      
c     Stiffness matrix and force vector are initialized to zero
      RHS(:) = zero
      AMATRX(:,:) = zero
      
c     Material behaviour
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      coef = E/(one-nu*nu)
      matH = zero

c     Loads
      kk = 0
      nb_load_bnd = 0
      nb_load_srf = 0
      ind_load_loc(:) = 0
      Do i = 1,nb_load
         If (JDLTYPE(i)/10>0 .AND. 
     &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i)) == JELEM))then
            If (JDLType(i)/10 < 5) then
               nb_load_bnd = nb_load_bnd + 1
               ind_load_loc(nb_load_bnd) = i
            Else
               nb_load_srf = nb_load_srf + 1
               ind_load_loc(nb_load+1-nb_load_srf) = i
            Endif
         Endif
         kk = kk + load_target_nbelem(i)
      Enddo
      
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
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         
         
c     Computing Curvilinear Coordinate objects
         call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,R,dRdxi,ddRddxi,MCRD,
     &        NNODE,COORDS)
         
         
c     Computing material matrix
         matH(:,:) = zero
         matH(1,1) = AAE(1,1)*AAE(1,1)
         matH(2,2) = AAE(2,2)*AAE(2,2)
         matH(3,3) = 0.5d0*
     &        ((one-nu)*AAE(1,1)*AAE(2,2) + (one+nu)*AAE(1,2)*AAE(1,2))
         matH(1,2) = nu*AAE(1,1)*AAE(2,2) + (one-nu)*AAE(1,2)*AAE(1,2)
         matH(1,3) = AAE(1,1)*AAE(1,2)
         matH(2,3) = AAE(2,2)*AAE(1,2)
         matH(2,1) = matH(1,2)
         matH(3,1) = matH(1,3)
         matH(3,2) = matH(2,3)
         matH(:,:) = coef*matH(:,:)
         
         
c     Computing stiffness matrix
         membraneStiff = zero
         matH(:,:) = h*matH(:,:)
         call usfmem_shell(NNODE,MCRD,NDOFEL,matH,AI,dRdxi,
     &        membraneStiff)
         
         bendingStiff = zero
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         call usfbnd_shell(NNODE,MCRD,NDOFEL,matH,AI,dAI1dxi,dAI2dxi,
     &        dRdxi,ddRddxi,bendingStiff)
         
c     
c     Assembling RHS and AMATRIX
         call SurfElem(AI(1,:), AI(2,:), Area)
         dvol = GaussPdsCoord(1,n)*DetJac*Area
         do k2 = 1,NDOFEL
            do k1 = 1,k2 !NDOFEL
               amatrx(k1,k2) = amatrx(k1,k2) 
     1              + (membraneStiff(k1,k2)+bendingStiff(k1,k2)) * dvol
            enddo
         enddo
         
c     Surface loads
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            if (KTypeDload == 0) then
               if (KNumFace == 5) then
                  VectNorm(:) = AI(3,:)
               else
                  VectNorm(:) =-AI(3,:)
               endif
            elseif (KTypeDload == 1) then
               VectNorm(:) = (/ one, zero, zero/)
            elseif (KTypeDload == 2) then
               VectNorm(:) = (/ zero, one, zero/)
            elseif (KTypeDload == 3) then
               VectNorm(:) = (/ zero, zero, one/)
            elseif (KTypeDload == 4) then
               call norm(AI(1,:), MCRD, normV)
               VectNorm(:) = AI(1,:)/normV
            elseif (KTypeDload == 5) then
               call norm(AI(2,:), MCRD, normV)
               VectNorm(:) = AI(2,:)/normV
            elseif (KTypeDload == 6) then
               !VectNorm(:) = AI(3,:)
               VectNorm(:) = (/ zero, zero, AI(3,3) /) ! snow load
            endif
            
            do numPC = 1,NNODE
               do k2 = 1,MCRD
                  Numi = (numPC-1)*MCRD+k2
                  RHS(Numi) = RHS(Numi) + 
     1                 ADLMAG(numLoad)*R(numPC)*VectNorm(k2)*dvol
               enddo
            enddo
            
         Enddo
c     End load traitement
         
      enddo
c     End of the loop on integration points on main surf
      
c     Symmetry
      Do k2 = 1,NDOFEL-1
         Do k1 = k2+1,NDOFEL
            AMATRX(k1,k2) = AMATRX(k2,k1)
         Enddo
      Enddo
      
      
      
c     Build of element stiffness matrix done.
c
c     ..................................................................
c
c     Loop on load : find boundary loads
      Do i = 1,nb_load_bnd
         numLoad = ind_load_loc(i)
         call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
         
         call Gauss(NbPtInt,2,GaussPdsCoord,KNumFace)
         Do n = 1,NbPtInt
            call shapPress3(VectNorm,R,DetJac, KTypeDload,KNumFace,
     &           COORDS,GaussPdsCoord(2:,n),MCRD)
            
            dvol = h*ADLMAG(numLoad)*GaussPdsCoord(1,n)*DetJac
            R(:) = R(:)*dvol
            Numi = 0
            Do numPC = 1,NNODE
               Do k2 = 1,MCRD
                  dof = Numi+k2
                  RHS(dof) = RHS(dof) + R(numPC)*VectNorm(k2)
               Enddo
               Numi = Numi + MCRD
            Enddo
         Enddo
      Enddo
      
      
      End SUBROUTINE UELMAT3






































      SUBROUTINE UELMAT3_byCP(NDOFEL,MCRD,NNODE,JELEM,NBINT,COORDS,
     1     TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,indDLoad,
     2     load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX)
      
      use parameters
    
      implicit none
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
      Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
      Double precision, intent(in) :: ADLMAG
      dimension ADLMAG(nb_load),
     &     load_target_nbelem(nb_load),JDLType(nb_load)
      dimension indDLoad(SUM(load_target_nbelem))
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: RHS, AMATRX
      dimension RHS(NDOFEL), AMATRX(MCRD,MCRD,NNODE*(NNODE+1)/2)
      
      
c     Local variables :
c     ---------------
      
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
      Double precision :: membraneStiff,bendingStiff,dvol,Area,normV,
     &     VectNorm
      dimension VectNorm(MCRD),
     &     membraneStiff(MCRD,MCRD,NNODE*(NNODE+1)/2 ),
     &     bendingStiff( MCRD,MCRD,NNODE*(NNODE+1)/2 )
      
!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,kk,numPC,Numi,dof, KTypeDload, KNumFace
      
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
      
c     Stiffness matrix and force vector are initialized to zero
      RHS(:)        = zero
      AMATRX(:,:,:) = zero
      
c     Material behaviour
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      coef = E/(one-nu*nu)
      matH = zero

c     Loads
      kk = 0
      nb_load_bnd = 0
      nb_load_srf = 0
      ind_load_loc(:) = 0
      Do i = 1,nb_load
         If (JDLTYPE(i)/10>0 .AND. 
     &        ANY(indDLoad(kk+1:kk+load_target_nbelem(i)) == JELEM))then
            If (JDLType(i)/10 < 5) then
               nb_load_bnd = nb_load_bnd + 1
               ind_load_loc(nb_load_bnd) = i
            Else
               nb_load_srf = nb_load_srf + 1
               ind_load_loc(nb_load+1-nb_load_srf) = i
            Endif
         Endif
         kk = kk + load_target_nbelem(i)
      Enddo
      
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
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         
         
c     Computing Curvilinear Coordinate objects
         call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,R,dRdxi,ddRddxi,MCRD,
     &        NNODE,COORDS)
         
         
c     Computing material matrix
         matH(:,:) = zero
         matH(1,1) = AAE(1,1)*AAE(1,1)
         matH(2,2) = AAE(2,2)*AAE(2,2)
         matH(3,3) = 0.5d0*
     &        ((one-nu)*AAE(1,1)*AAE(2,2) + (one+nu)*AAE(1,2)*AAE(1,2))
         matH(1,2) = nu*AAE(1,1)*AAE(2,2) + (one-nu)*AAE(1,2)*AAE(1,2)
         matH(1,3) = AAE(1,1)*AAE(1,2)
         matH(2,3) = AAE(2,2)*AAE(1,2)
         matH(2,1) = matH(1,2)
         matH(3,1) = matH(1,3)
         matH(3,2) = matH(2,3)
         matH(:,:) = coef*matH(:,:)
         
         
c     Computing stiffness matrix
         membraneStiff = zero
         matH(:,:) = h*matH(:,:)
         call usfmem_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dRdxi,
     &        membraneStiff)
         
         bendingStiff = zero
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         call usfbnd_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dAI1dxi,
     &        dAI2dxi,dRdxi,ddRddxi,bendingStiff)
         
c     
c     Assembling RHS and AMATRIX
         call SurfElem(AI(1,:), AI(2,:), Area)
         dvol = GaussPdsCoord(1,n)*DetJac*Area
         
         AMATRX(:,:,:) = AMATRX(:,:,:)
     &        + ( membraneStiff(:,:,:) + bendingStiff(:,:,:) ) * dvol
         
c     Surface loads
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            if (KTypeDload == 0) then
               if (KNumFace == 5) then
                  VectNorm(:) = AI(3,:)
               else
                  VectNorm(:) =-AI(3,:)
               endif
            elseif (KTypeDload == 1) then
               VectNorm(:) = (/ one, zero, zero/)
            elseif (KTypeDload == 2) then
               VectNorm(:) = (/ zero, one, zero/)
            elseif (KTypeDload == 3) then
               VectNorm(:) = (/ zero, zero, one/)
            elseif (KTypeDload == 4) then
               call norm(AI(1,:), MCRD, normV)
               VectNorm(:) = AI(1,:)/normV
            elseif (KTypeDload == 5) then
               call norm(AI(2,:), MCRD, normV)
               VectNorm(:) = AI(2,:)/normV
            elseif (KTypeDload == 6) then
               !VectNorm(:) = AI(3,:)
               VectNorm(:) = (/ zero, zero, AI(3,3) /) ! snow load
            endif
            
            do numPC = 1,NNODE
               do k2 = 1,MCRD
                  Numi = (numPC-1)*MCRD+k2
                  RHS(Numi) = RHS(Numi) + 
     1                 ADLMAG(numLoad)*R(numPC)*VectNorm(k2)*dvol
               enddo
            enddo
            
         Enddo
c     End load traitement
         
      enddo
c     End of the loop on integration points on main surf
      
      
      
c     Build of element stiffness matrix done.
c
c     ..................................................................
c
c     Loop on load : find boundary loads
      Do i = 1,nb_load_bnd
         numLoad = ind_load_loc(i)
         call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
         
         call Gauss(NbPtInt,2,GaussPdsCoord,KNumFace)
         Do n = 1,NbPtInt
            call shapPress3(VectNorm,R,DetJac, KTypeDload,KNumFace,
     &           COORDS,GaussPdsCoord(2:,n),MCRD)
            
            dvol = h*ADLMAG(numLoad)*GaussPdsCoord(1,n)*DetJac
            R(:) = R(:)*dvol
            Numi = 0
            Do numPC = 1,NNODE
               Do k2 = 1,MCRD
                  dof = Numi+k2
                  RHS(dof) = RHS(dof) + R(numPC)*VectNorm(k2)
               Enddo
               Numi = Numi + MCRD
            Enddo
         Enddo
      Enddo
      
      
      End SUBROUTINE UELMAT3_byCP


