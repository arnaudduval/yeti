!! Copyright 2018 Thibaut Hirschler

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
c     Cas ou le champ des deplacements est approche par la param. de
c     la surface immergee
c     --
      
      SUBROUTINE UELMAT30(NDOFEL,MCRD,NNODE,NNODEmap,nb_cp,JELEM,NBINT,
     1     COORDS,COORDSall,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,
     2     nb_load,indDLoad,load_target_nbelem,JDLType,ADLMAG,
     3     RHS,AMATRX)
      
      use parameters
      use embeddedMapping
      
      implicit none      
      
c     Input arguments :
c     ---------------
!     Surface
      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(MCRD,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
!     Mapping
      Integer, intent(in) :: NNODEmap, nb_cp
      Double precision, intent(in) :: COORDSall
      dimension COORDSall(3,nb_cp)
      
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
      
!     Embedded Surface
      ! - nurbs basis functions
      Double precision :: Rs, dRdTheta, ddRddTheta, DetJac
      dimension Rs(NNODE), dRdTheta(NNODE,3), ddRddTheta(NNODE,6)
      ! - curvilinear coordinate objects
      Double precision :: XI,BI,dBI
      dimension XI(3),BI(3,2),dBI(3,3)
      
!     Mapping
      ! - nurbs basis functions
      Double precision :: R, dRdxi, ddRddxi
      dimension R(NNODEmap), dRdxi(NNODEmap,3), ddRddxi(NNODEmap,6)
      ! - element infos
      Double precision :: COORDSmap
      dimension COORDSmap(MCRD,NNODEmap)
      Integer          :: sctr_map
      dimension sctr_map(NNODEmap)
            
!     Composition Mapping+Surface
      Double precision :: dRRdTheta, ddRRddTheta, AI, dAI1dTheta,
     &     dAI2dTheta, AAE
      dimension dRRdTheta(NNODEmap,2),ddRRddTheta(NNODEmap,3),AI(3,3),
     &     dAI1dTheta(2,3), dAI2dTheta(2,3), AAE(2,2)
      
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
      Integer ntens,isave
      Integer n,k1,k2,i,j,kk,numCP,Numi,dof, KTypeDload, KNumFace
      
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
      h = PROPS(3)
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
         If (ANY(indDLoad(kk+1:kk+load_target_nbelem(i)) == JELEM)) then
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
      
      isave = 0

c     Loop on integration points on main surface
      do n = 1,NBINT
         
c     Computing NURBS basis functions and derivatives
         DetJac = zero
         Rs(:)  = zero
         dRdTheta(:,:)   = zero
         ddRddTheta(:,:) = zero
         call nurbsbasis(Rs,dRdTheta(:,:2),ddRddTheta(:,:3),DetJac,
     &        GaussPdsCoord(2:,n))
         
c     Find mapping parametric position
         XI(:)    = zero
         BI(:,:)  = zero
         dBI(:,:) = zero
         Do numCP = 1,NNODE
            XI(:)    =  XI(:)   +           Rs(numCP)*COORDS(:,numCP)
            BI(:,1)  =  BI(:,1) +   dRdTheta(numCP,1)*COORDS(:,numCP)
            BI(:,2)  =  BI(:,2) +   dRdTheta(numCP,2)*COORDS(:,numCP)
            dBI(:,1) = dBI(:,1) + ddRddTheta(numCP,1)*COORDS(:,numCP)
            dBI(:,2) = dBI(:,2) + ddRddTheta(numCP,2)*COORDS(:,numCP)
            dBI(:,3) = dBI(:,3) + ddRddTheta(numCP,3)*COORDS(:,numCP)
         Enddo

c     Computing NURBS basis functions and derivatives of the mapping
         ! get active element number
         call updateMapElementNumber(XI(:))

         call evalnurbs_mapping_w2ndDerv(XI(:),R(:),dRdxi(:,:),
     &        ddRddxi(:,:))
         
c     Composition of the basis functions
         dRRdTheta(:,:) = zero
         Do i = 1,3
            dRRdTheta(:,1) = dRRdTheta(:,1) + BI(i,1)*dRdxi(:,i)
            dRRdTheta(:,2) = dRRdTheta(:,2) + BI(i,2)*dRdxi(:,i)
         Enddo
         
         ddRRddTheta(:,:) = zero
         Do i = 1,3
            ddRRddTheta(:,1) = ddRRddTheta(:,1)
     &           + dBI(i,1)*dRdxi(:,i) + BI(i,1)*BI(i,1)*ddRddxi(:,i)
            ddRRddTheta(:,2) = ddRRddTheta(:,2) 
     &           + dBI(i,2)*dRdxi(:,i) + BI(i,2)*BI(i,2)*ddRddxi(:,i)
            ddRRddTheta(:,3) = ddRRddTheta(:,3) 
     &           + dBI(i,3)*dRdxi(:,i) + BI(i,1)*BI(i,2)*ddRddxi(:,i)
            Do j = i+1,3
               kk = i+j+1
               ddRRddTheta(:,1) = ddRRddTheta(:,1)
     &              + two*BI(i,1)*BI(j,1)*ddRddxi(:,kk)
               ddRRddTheta(:,2) = ddRRddTheta(:,2)
     &              + two*BI(i,2)*BI(j,2)*ddRddxi(:,kk)
               ddRRddTheta(:,3) = ddRRddTheta(:,3)
     &              +(BI(i,1)*BI(j,2) + BI(j,1)*BI(i,2))*ddRddxi(:,kk)
            Enddo
         Enddo

                  
c     Computing Curvilinear Coordinate objects
         ! extract COORDS
         If (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:,current_map_elem)
            
            Do numCP = 1,NNODEmap
               COORDSmap(:,numCP) = COORDSall(:,sctr_map(numCP))
            Enddo
            
            isave = current_map_elem
         Endif
         
         call curvilinear(AI,dAI1dTheta,dAI2dTheta,AAE,R,dRRdTheta,
     &        ddRRddTheta,MCRD,NNODEmap,COORDSmap)
         
         
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
         call usfmem_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dRdTheta,
     &        membraneStiff)
         
         bendingStiff = zero
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         call usfbnd_shell_byCP(NNODE,MCRD,NDOFEL,matH,AI,dAI1dTheta,
     &        dAI2dTheta,dRdTheta,ddRddTheta,bendingStiff)
         
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
            
            do numCP = 1,NNODE
               do k2 = 1,MCRD
                  Numi = (numCP-1)*MCRD+k2
                  RHS(Numi) = RHS(Numi) + 
     1                 ADLMAG(numLoad)*Rs(numCP)*VectNorm(k2)*dvol
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
      isave = 0
      Do j = 1,nb_load_bnd
         numLoad = ind_load_loc(j)
         call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
 
         !print*,'Warning: load case',JDLType(numLoad),'is not available'
        
         call Gauss(NbPtInt,2,GaussPdsCoord,KNumFace)
         Do n = 1,NbPtInt
            
            XI(:) = zero
            Do i = 1,2
               coef = GaussPdsCoord(1+i,n)
               XI(i)= ((Ukv_elem(2,i) - Ukv_elem(1,i))*coef
     &                +(Ukv_elem(2,i) + Ukv_elem(1,i))     )*0.5d0
            Enddo
            
            call evalnurbs(XI,Rs,dRdTheta)
            
            XI(:)    = zero
            BI(:,:)  = zero
            Do numCP = 1,NNODE
               XI(:)    =  XI(:)   +         Rs(numCP)*COORDS(:,numCP)
               BI(:,1)  =  BI(:,1) + dRdTheta(numCP,1)*COORDS(:,numCP)
               BI(:,2)  =  BI(:,2) + dRdTheta(numCP,2)*COORDS(:,numCP)
            Enddo
            
            !print*,'XI',XI

c     Computing NURBS basis functions and derivatives of the mapping
            ! get active element number
            call updateMapElementNumber(XI(:))
            
            call evalnurbs_mapping(XI(:),R(:),dRdxi(:,:))
            
c     Composition of the basis functions
            dRRdTheta(:,:) = zero
            Do i = 1,3
               dRRdTheta(:,1) = dRRdTheta(:,1) + BI(i,1)*dRdxi(:,i)
               dRRdTheta(:,2) = dRRdTheta(:,2) + BI(i,2)*dRdxi(:,i)
            Enddo
            
c     Computing Curvilinear Coordinate objects
            ! extract COORDS
            If (isave /= current_map_elem) then
               sctr_map(:) = IEN_map(:,current_map_elem)
               
               Do numCP = 1,NNODEmap
                  COORDSmap(:,numCP) = COORDSall(:,sctr_map(numCP))
               Enddo
               
               isave = current_map_elem
            Endif
            
            AI(:,:) = zero
            Do numCP = 1,NNODEmap
               AI(:,1) = AI(:,1) + dRRdTheta(numCP,1)*COORDSmap(:,numCP)
               AI(:,2) = AI(:,2) + dRRdTheta(numCP,2)*COORDSmap(:,numCP)
            Enddo
            call cross(AI(:,1),AI(:,2),AI(:,3))
            call norm(AI(:,3), MCRD, normV)
            AI(:,3) = AI(:,3)/normV
            
c     Compute jacobian from parent element to physical space
            SELECT CASE (KNumFace)
            case(1,2)
               call norm(AI(:,2),3,normV)
               DetJac = 0.5d0*( Ukv_elem(2,2) - Ukv_elem(1,2) )*normV
            case(3,4)
               call norm(AI(:,1),3,normV)
               DetJac = 0.5d0*( Ukv_elem(2,1) - Ukv_elem(1,1) )*normV
            END SELECT
            
c     Compute load direction
            SELECT CASE (KTypeDload)
!     -
!     Pression normale a la surface chargee
            case(0)
               SELECT CASE (KNumFace)
               case(1)
                  call cross(AI(:,2), AI(:,3), VectNorm(:))
               case(2)
                  call cross(AI(:,3), AI(:,2), VectNorm(:))
               case(3)
                  call cross(AI(:,3), AI(:,1), VectNorm(:))
               case(4)
                  call cross(AI(:,1), AI(:,3), VectNorm(:))
               END SELECT
!     -
!     Pression tangentielle dans le plan median
            case(9)
               SELECT CASE (KNumFace)
               case(1,2)
                  VectNorm(:) = AI(:,2)
               case(3,4)
                  VectNorm(:) = AI(:,1)
               END SELECT
!     -
!     Pression tangentielle hors plan
            case(8)
               VectNorm(:) = AI(:,3)
!     -
!     Pression dans direction x, y, ou z
            case(1:3)
               VectNorm(:) = zero
               VectNorm(KTypeDload) = one
            END SELECT
            
            call norm(VectNorm(:),3,normV)
            VectNorm(:) = VectNorm(:)/normV
            
c     Update
            dvol = h*ADLMAG(numLoad)*GaussPdsCoord(1,n)*DetJac
            Rs(:) = Rs(:)*dvol
            Numi = 0
            Do numCP = 1,NNODE
               Do k2 = 1,MCRD
                  dof = Numi+k2
                  RHS(dof) = RHS(dof) + Rs(numCP)*VectNorm(k2)
               Enddo
               Numi = Numi + MCRD
            Enddo
         Enddo
      Enddo
      
      End SUBROUTINE UELMAT30

 




















      




c     --
c     Cas ou le champ des deplacements est approche par la param. du
c     mapping
c     --
      
      SUBROUTINE UELMAT31(XI,BI,dBI,DetJac,NDOFEL,MCRD,NNODE,JELEM,
     1     COORDS,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,
     2     indDLoad,load_target_nbelem,JDLType,ADLMAG,RHS,AMATRX)
      
      use parameters
    
      implicit none      
      
c     Input arguments :
c     ---------------
      Double precision, intent(in) :: XI,BI,dBI,DetJac
      dimension XI(3),BI(3,3),dBI(3,6)

      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,JPROPS
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
      
!     For nurbs basis functions
      Double precision :: R,dRdxi,ddRddxi,dRdTheta,ddRddTheta
      dimension R(NNODE), dRdxi(NNODE,3), ddRddxi(NNODE,6),
     &     dRdTheta(NNODE,2),ddRddTheta(NNODE,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI, dAI1dTheta, dAI2dTheta, AAE,AAI,scal,det
      dimension AI(3,3),dAI1dTheta(2,3),dAI2dTheta(2,3),AAE(2,2),
     &     AAI(2,2)
      
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
      Integer n,k1,k2,i,j,k,kk,numPC,Numi,dof, KTypeDload, KNumFace
      
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      
      
c     Stiffness matrix and force vector are initialized to zero
      RHS(:) = zero
      AMATRX(:,:) = zero
      
c     Material behaviour
      h = PROPS(3)
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
         If (ANY(indDLoad(kk+1:kk+load_target_nbelem(i)) == JELEM)) then
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
      
c     Computing NURBS basis functions and derivatives
      R = zero
      dRdxi = zero
      ddRddxi = zero
      call evalnurbs_mapping_w2ndDerv(xi,R(:),dRdxi(:,:),ddRddxi(:,:))
      
      dRdTheta(:,:)   = zero
      Do i = 1,3
         dRdTheta(:,1) = dRdTheta(:,1) + BI(i,1)*dRdxi(:,i)
         dRdTheta(:,2) = dRdTheta(:,2) + BI(i,2)*dRdxi(:,i)
      Enddo
      
      ddRddTheta(:,:) = zero
      Do i = 1,3
         ddRddTheta(:,1) = ddRddTheta(:,1)
     &        + dBI(i,1)*dRdxi(:,i) + BI(i,1)*BI(i,1)*ddRddxi(:,i)
         ddRddTheta(:,2) = ddRddTheta(:,2) 
     &        + dBI(i,2)*dRdxi(:,i) + BI(i,2)*BI(i,2)*ddRddxi(:,i)
         ddRddTheta(:,3) = ddRddTheta(:,3) 
     &        + dBI(i,4)*dRdxi(:,i) + BI(i,1)*BI(i,2)*ddRddxi(:,i)
         Do j = i+1,3
            k = i+j+1
            ddRddTheta(:,1) = ddRddTheta(:,1)
     &           + two*BI(i,1)*BI(j,1)*ddRddxi(:,k)
            ddRddTheta(:,2) = ddRddTheta(:,2)
     &           + two*BI(i,2)*BI(j,2)*ddRddxi(:,k)
            ddRddTheta(:,3) = ddRddTheta(:,3)
     &           + two*BI(i,1)*BI(j,2)*ddRddxi(:,k)
         Enddo
      Enddo

c     Curvilinear Basis Vectors
      AI(:,:) = zero
      Do i = 1,NNODE
         AI(1,:) = AI(1,:) + dRdTheta(i,1)*COORDS(:,i)
         AI(2,:) = AI(2,:) + dRdTheta(i,2)*COORDS(:,i)
      Enddo
      call cross(AI(1,:), AI(2,:), AI(3,:))
      area    =sqrt(AI(3,1)*AI(3,1) + AI(3,2)*AI(3,2) + AI(3,3)*AI(3,3))
      AI(3,:) = AI(3,:)/area
      
c     Curvilinear Derivatives
      dAI1dTheta(:,:) = zero
      dAI2dTheta(:,:) = zero
      Do i = 1,NNODE
         dAI1dTheta(1,:) = dAI1dTheta(1,:) + ddRddTheta(i,1)*COORDS(:,i)
         dAI1dTheta(2,:) = dAI1dTheta(2,:) + ddRddTheta(i,3)*COORDS(:,i)
         
         dAI2dTheta(1,:) = dAI2dTheta(1,:) + ddRddTheta(i,3)*COORDS(:,i)
         dAI2dTheta(2,:) = dAI2dTheta(2,:) + ddRddTheta(i,2)*COORDS(:,i)
      Enddo
      
c     Mertic coefficients
      AAI = zero
      Do i = 1,2
         Do j = 1,2
            call dot(AI(i,:), AI(j,:), scal)
            AAI(i,j) = scal
         Enddo
      Enddo
      call MatrixInv(AAE, AAI, det, 2)
      
!     call curvilinear(AI,dAI1dTheta,dAI2dTheta,AAE,R,dRdxi,ddRddxi,MCRD,
!     &     NNODE,COORDS)
      
      
c     Computing material matrix
      matH(:,:) = zero
      matH(1,1) = AAE(1,1)*AAE(1,1)
      matH(2,2) = AAE(2,2)*AAE(2,2)
      matH(3,3) = 0.5d0*
     &     ((one-nu)*AAE(1,1)*AAE(2,2) + (one+nu)*AAE(1,2)*AAE(1,2))
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
      call usfmem_shell(NNODE,MCRD,NDOFEL,matH,AI,dRdTheta,
     &     membraneStiff)
      
      bendingStiff = zero
      matH(:,:) = h**two/12.0D0 * matH(:,:)
      call usfbnd_shell(NNODE,MCRD,NDOFEL,matH,AI,dAI1dTheta,dAI2dTheta,
     &     dRdTheta,ddRddTheta,bendingStiff)
      
c     
c     Assembling RHS and AMATRIX
      dvol = DetJac*Area
      do k2 = 1,NDOFEL
         do k1 = 1,k2           !NDOFEL
            amatrx(k1,k2) = amatrx(k1,k2) 
     1           + (membraneStiff(k1,k2)+bendingStiff(k1,k2)) * dvol
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
!     VectNorm(:) = AI(3,:)
            VectNorm(:) = (/ zero, zero, AI(3,3) /) ! snow load
         endif
         
         do numPC = 1,NNODE
            do k2 = 1,MCRD
               Numi = (numPC-1)*MCRD+k2
               RHS(Numi) = RHS(Numi) + 
     1              ADLMAG(numLoad)*R(numPC)*VectNorm(k2)*dvol
            enddo
         enddo
         
      Enddo
c     End load traitement
      
c     Symmetry
      Do k2 = 1,NDOFEL-1
         Do k1 = k2+1,NDOFEL
            AMATRX(k1,k2) = AMATRX(k2,k1)
         Enddo
      Enddo
      
c     Build of element stiffness matrix done.
c     
      End SUBROUTINE UELMAT31
