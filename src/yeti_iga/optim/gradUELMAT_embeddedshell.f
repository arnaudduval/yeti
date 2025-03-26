!! Copyright 2018-2020 Thibaut Hirschler

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
c     Construction de la sentitivite de  la matrice de rigite et du
c     second membre elementaire pour la formulation coque type 
c     Kirchhoff-Love immergeee
c      - sensibilite par rapport aux PCs de la surface immergee
c     --
      
      SUBROUTINE gradUELMAT30(activeDir,Uelem,NDOFEL,MCRD,NNODE,
     1     NNODEmap,nb_cp,JELEM,NBINT,COORDS,COORDSall,TENSOR,
     2     MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,indDLoad,
     3     load_target_nbelem,JDLType,ADLMAG,gradC_elem)
      
      use parameters
      use embeddedMapping
            
      Implicit None
      
c     Input arguments :
c     ---------------
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
 
      Double precision, intent(in) :: Uelem
      Integer, intent(in)          :: activeDir
      dimension Uelem(3,NNODE),activeDir(3)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: gradC_elem
      dimension gradC_elem(MCRD,NNODE)
      
      
c     Local variables :
c     ---------------
      
!     For gauss points
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(3,NBINT)
      
!     For nurbs basis functions
!     - embedded surface
      Double precision :: XI, R, dRdxi, ddRddxi, DetJac
      dimension XI(3), R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
!     Mapping
      ! - nurbs basis functions
      Double precision :: Rm, dRmdxi, ddRmddxi
      dimension Rm(NNODEmap), dRmdxi(NNODEmap,3), ddRmddxi(NNODEmap,6)
      ! - element infos
      Double precision :: COORDSmap
      dimension COORDSmap(MCRD,NNODEmap)
      Integer          :: sctr_map,isave
      dimension sctr_map(NNODEmap)
      
!     Composition Mapping+Surface
      Double precision :: dRRdxi, ddRRddxi
      dimension dRRdxi(NNODEmap,2),ddRRddxi(NNODEmap,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI,dAI1dxi,dAI2dxi,AAE,AIxAJ,BI,dBI,VI,dVI
      dimension AI(3,3),dAI1dxi(2,3),dAI2dxi(2,3),AAE(2,2),AIxAJ(3,3),
     &     BI(3,2),dBI(3,3),VI(3,3),dVI(3,6)

!     For material matrix
      Double precision :: E, nu, h, matH, coef
      dimension matH(3,3)
      
      
!     For Membrane term
      Double precision :: Ua,Va,Xa,Ub,Vb,Xb,mBU,mHBU,dAAIdPe,dAAEdPe,
     &     dHdPe,dHdPe_mBU,mBU_dHdPe_mBU,dJdPe,sKE,dvol,Area,dA1dPe,
     &     dA2dPe
      dimension Ua(3),Va(3),Xa(3),Ub(3),Vb(3),Xb(3),mBU(3,NNODE),
     &     mHBU(3,NNODE),dAAIdPe(3,2,2),dAAEdPe(3,2,2),
     &     dHdPe(3,NNODE,3,3),dHdPe_mBU(3,NNODE,3),
     &     mBU_dHdPe_mBU(3,NNODE),dJdPe(3,NNODE),dA1dPe(3,3,NNODE),
     &     dA2dPe(3,3,NNODE)
      
!     For Bending term
      Double precision :: B1,B2,B3,dA1d1_A2,dA2d2_A2,dA1d2_A2,A1_dA1d1,
     &     A1_dA2d2,A1_dA1d2,A3dA1d1,A3dA2d2,A3dA1d2,fBU,fHBU,dHdPe_fBU,
     &     fBU_dHdPe_fBU,Wa,Ya,Wb,Yb,UxA2,A1xU,A3xU,dA3dPe_U,V,VxA2,
     &     A1xV,dfBdPeU,UxdA1d1,UxdA2d2,UxdA1d2,save1,save2,save3,save4,
     &     UA2xA2,A1xUA2,A1UxA2,A1xA1U,d_A3dA1d1_dPe,d_A3dA2d2_dPe,
     &     d_A3dA1d2_dPe,ddA1d1Pe,ddA2d2Pe,ddA1d2Pe
      dimension B1(3),B2(3),B3(3),dA1d1_A2(3),dA2d2_A2(3),dA1d2_A2(3),
     &     A1_dA1d1(3),A1_dA2d2(3),A1_dA1d2(3),fBU(3,NNODE),
     &     fHBU(3,NNODE),dHdPe_fBU(3,NNODE,3),
     &     fBU_dHdPe_fBU(3,NNODE),Wa(3),Ya(3),Wb(3),Yb(3),UxA2(3),
     &     A1xU(3),A3xU(3),dA3dPe_U(3),V(3),VxA2(3),A1xV(3),
     &     dfBdPeU(3,NNODE,3,NNODE),UxdA1d1(3),UxdA2d2(3),UxdA1d2(3),
     &     UA2xA2(3),A1xUA2(3),A1UxA2(3),A1xA1U(3),
     &     d_A3dA1d1_dPe(3,NNODE),d_A3dA2d2_dPe(3,NNODE),
     &     d_A3dA1d2_dPe(3,NNODE),ddA1d1Pe(3,3,NNODE),
     &     ddA2d2Pe(3,3,NNODE),ddA1d2Pe(3,3,NNODE)


!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      Double precision :: VectNorm_U,normV
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,j,k,kk,ll,cp,cpa,cpb, KTypeDload, KNumFace
      Double precision :: temp,temp1,temp2
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
            
c     Stiffness matrix and force vector are initialized to zero
      gradC_elem  = zero
      
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
         
c     PRELIMINARY QUANTITES
c     Computing NURBS basis functions and derivatives
         R      = zero
         dRdxi  = zero
         ddRddxi= zero
         DetJac = zero
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         DetJac = DetJac * GaussPdsCoord(1,n)
         
c     Find mapping parametric position
         XI(:)    = zero
         BI(:,:)  = zero
         dBI(:,:) = zero
         Do cp = 1,NNODE
            XI(:)    =  XI(:)   +         R(cp)*COORDS(:,cp)
            BI(:,1)  =  BI(:,1) +   dRdxi(cp,1)*COORDS(:,cp)
            BI(:,2)  =  BI(:,2) +   dRdxi(cp,2)*COORDS(:,cp)
            dBI(:,1) = dBI(:,1) + ddRddxi(cp,1)*COORDS(:,cp)
            dBI(:,2) = dBI(:,2) + ddRddxi(cp,2)*COORDS(:,cp)
            dBI(:,3) = dBI(:,3) + ddRddxi(cp,3)*COORDS(:,cp)
         Enddo
         
c     Computing NURBS basis functions and derivatives of the mapping
         ! get active element number
         call updateMapElementNumber(XI(:))
         
         call evalnurbs_mapping_w2ndDerv(XI(:),Rm(:),dRmdxi(:,:),
     &        ddRmddxi(:,:))

c     Composition of the basis functions
         dRRdxi(:,:) = zero
         Do i = 1,3
            dRRdxi(:,1) = dRRdxi(:,1) + BI(i,1)*dRmdxi(:,i)
            dRRdxi(:,2) = dRRdxi(:,2) + BI(i,2)*dRmdxi(:,i)
         Enddo
         
         ddRRddxi(:,:) = zero
         Do i = 1,3
            ddRRddxi(:,1) = ddRRddxi(:,1)
     &           + dBI(i,1)*dRmdxi(:,i) + BI(i,1)*BI(i,1)*ddRmddxi(:,i)
            ddRRddxi(:,2) = ddRRddxi(:,2) 
     &           + dBI(i,2)*dRmdxi(:,i) + BI(i,2)*BI(i,2)*ddRmddxi(:,i)
            ddRRddxi(:,3) = ddRRddxi(:,3) 
     &           + dBI(i,3)*dRmdxi(:,i) + BI(i,1)*BI(i,2)*ddRmddxi(:,i)
            Do j = i+1,3
               kk = i+j+1
               ddRRddxi(:,1) = ddRRddxi(:,1)
     &              + two*BI(i,1)*BI(j,1)*ddRmddxi(:,kk)
               ddRRddxi(:,2) = ddRRddxi(:,2)
     &              + two*BI(i,2)*BI(j,2)*ddRmddxi(:,kk)
               ddRRddxi(:,3) = ddRRddxi(:,3)
     &              +(BI(i,1)*BI(j,2) + BI(j,1)*BI(i,2))*ddRmddxi(:,kk)
            Enddo
         Enddo
         
c     Computing Curvilinear Coordinate objects
         ! extract COORDS
         If (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:,current_map_elem)
            
            Do cp = 1,NNODEmap
               COORDSmap(:,cp) = COORDSall(:,sctr_map(cp))
            Enddo
            
            isave = current_map_elem
         Endif
         
         call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,Rm,dRRdxi,ddRRddxi,
     &        MCRD,NNODEmap,COORDSmap)
         
         VI(:,:) = zero
         Do k = 1,3
            Do cp = 1,NNODEmap
               VI(:,k) = VI(:,k) + dRmdxi(cp,k)*COORDSmap(:,cp)
            Enddo
         Enddo
         
         dVI(:,:) = zero
         Do k = 1,6
            Do cp = 1,NNODEmap
               dVI(:,k) = dVI(:,k) + ddRmddxi(cp,k)*COORDSmap(:,cp)
            Enddo
         Enddo
         
         call cross(AI(2,:),AI(3,:), AIxAJ(:,1))
         call cross(AI(3,:),AI(1,:), AIxAJ(:,2))
         call SurfElem(AI(1,:), AI(2,:), Area)
         
         call cross(dAI1dxi(1,:), AI(2,:), dA1d1_A2(:))
         call cross(dAI2dxi(2,:), AI(2,:), dA2d2_A2(:))
         call cross(dAI1dxi(2,:), AI(2,:), dA1d2_A2(:))
         call cross(AI(1,:), dAI1dxi(1,:), A1_dA1d1(:))
         call cross(AI(1,:), dAI2dxi(2,:), A1_dA2d2(:))
         call cross(AI(1,:), dAI1dxi(2,:), A1_dA1d2(:))
         call   dot(AI(3,:), dAI1dxi(1,:), A3dA1d1)
         call   dot(AI(3,:), dAI2dxi(2,:), A3dA2d2)
         call   dot(AI(3,:), dAI1dxi(2,:), A3dA1d2)
         
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
         matH(:,:) = (h*coef)*matH(:,:)
         
C     --
c     Preliminary derivatives
         ! - covariant basis vectors
         dA1dPe(:,:,:) = zero
         dA2dPe(:,:,:) = zero
         Do cp = 1,NNODE
            Do i = 1,3
               dA1dPe(:,i,cp) = dRdxi(cp,1)*VI(:,i) 
     &              + BI(i,1)*R(cp)*dVI(:,i)
               dA2dPe(:,i,cp) = dRdxi(cp,2)*VI(:,i)
     &              + BI(i,2)*R(cp)*dVI(:,i)
               Do j = 1,3
               If (j /= i) then
                  kk = i+j+1
                  dA1dPe(:,i,cp) = dA1dPe(:,i,cp)
     &                 + BI(j,1)*R(cp)*dVI(:,kk)
                  dA2dPe(:,i,cp) = dA2dPe(:,i,cp)
     &                 + BI(j,2)*R(cp)*dVI(:,kk)
               Endif
               Enddo
            Enddo
         Enddo
         
         ddA1d1Pe(:,:,:) = zero
         ddA1d2Pe(:,:,:) = zero
         ddA2d2Pe(:,:,:) = zero
         Do cp = 1,NNODE
            Do i = 1,3
               ddA1d1Pe(:,i,cp) = ddRddxi(cp,1)*VI(:,i)
     &              + dBI(i,1)*R(cp)*dVI(:,i)
     &              + two*dRdxi(cp,1)*BI(i,1)*dVI(:,i)
               ddA2d2Pe(:,i,cp) = ddRddxi(cp,2)*VI(:,i)
     &              + dBI(i,2)*R(cp)*dVI(:,i)
     &              + two*dRdxi(cp,2)*BI(i,2)*dVI(:,i)
               ddA1d2Pe(:,i,cp) = ddRddxi(cp,3)*VI(:,i)
     &              + dBI(i,3)*R(cp)*dVI(:,i)
     &              +(dRdxi(cp,1)*BI(i,2)+dRdxi(cp,2)*BI(i,1))*dVI(:,i)
               
               Do j = 1,3
               If (j /= i) then
                  kk = i+j+1
                  ddA1d1Pe(:,i,cp) = ddA1d1Pe(:,i,cp)
     &                 + dBI(j,1)*R(cp)*dVI(:,kk)
     &                 + two*dRdxi(cp,1)*BI(j,1)*dVI(:,kk)
                  ddA2d2Pe(:,i,cp) = ddA2d2Pe(:,i,cp)
     &                 + dBI(j,2)*R(cp)*dVI(:,kk)
     &                 + two*dRdxi(cp,2)*BI(j,2)*dVI(:,kk)
                  ddA1d2Pe(:,i,cp) = ddA1d2Pe(:,i,cp)
     &                 + dBI(j,3)*R(cp)*dVI(:,kk)
     &                 +( dRdxi(cp,1)*BI(j,2) 
     &                   + BI(j,1)*dRdxi(cp,2) )*dVI(:,kk)
               Endif
               Enddo
            Enddo
         Enddo
         
         ! - Jacobian
         dJdPe(:,:) = zero
         Do cp = 1,NNODE
            Do i = 1,3
               dJdPe(i,cp) = SUM( AIxAJ(:,1)*dA1dPe(:,i,cp) )
     &                     + SUM( AIxAJ(:,2)*dA2dPe(:,i,cp) ) 
            Enddo
         Enddo
         
         ! - for bending term (a_{3}*a_{i,j})
         d_A3dA1d1_dPe(:,:) = zero
         d_A3dA2d2_dPe(:,:) = zero
         d_A3dA1d2_dPe(:,:) = zero
         Do cp = 1,NNODE
            Do i = 1,3
               temp = dJdPe(i,cp)*A3dA1d1
     &              + SUM( dA1d1_A2(:)*dA1dPe(:,i,cp) )
     &              + SUM( A1_dA1d1(:)*dA2dPe(:,i,cp) )
               
               d_A3dA1d1_dPe(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA1d1Pe(:,i,cp) )

               temp = dJdPe(i,cp)*A3dA2d2
     &              + SUM( dA2d2_A2(:)*dA1dPe(:,i,cp) )
     &              + SUM( A1_dA2d2(:)*dA2dPe(:,i,cp) )
               
               d_A3dA2d2_dPe(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA2d2Pe(:,i,cp) )
               
               temp = dJdPe(i,cp)*A3dA1d2
     &              + SUM( dA1d2_A2(:)*dA1dPe(:,i,cp) )
     &              + SUM( A1_dA1d2(:)*dA2dPe(:,i,cp) )
               
               d_A3dA1d2_dPe(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA1d2Pe(:,i,cp) )
            Enddo
         Enddo
         
         

c     Matrices and vectors
         mBU(:,:)  = zero
         mHBU(:,:) = zero
         fBU(:,:)  = zero
         fHBU(:,:) = zero
         dHdPe(:,:,:,:) = zero
         Do cp = 1,NNODE
            
            ! membrane strain
            temp1 = SUM( AI(1,:)*Uelem(:,cp) )
            temp2 = SUM( AI(2,:)*Uelem(:,cp) )
            mBU(1,cp) = dRdxi(cp,1) * temp1
            mBU(2,cp) = dRdxi(cp,2) * temp2
            mBU(3,cp) = dRdxi(cp,1) * temp2 + dRdxi(cp,2) * temp1
            
            ! membrane stress
            mHBU(1,cp) = SUM( matH(:,1) * mBU(:,cp) )
            mHBU(2,cp) = SUM( matH(:,2) * mBU(:,cp) )
            mHBU(3,cp) = SUM( matH(:,3) * mBU(:,cp) )
            
            ! bending strain
            B1(:) = -ddRddxi(cp,1)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA1d1_A2(:) + dRdxi(cp,2)*A1_dA1d1(:)
     &           + A3dA1d1*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            B2(:) = -ddRddxi(cp,2)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA2d2_A2(:) + dRdxi(cp,2)*A1_dA2d2(:)
     &           + A3dA2d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            B3(:) = -ddRddxi(cp,3)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA1d2_A2(:) + dRdxi(cp,2)*A1_dA1d2(:)
     &           + A3dA1d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            fBU(1,cp) = SUM( B1(:)*Uelem(:,cp) )
            fBU(2,cp) = SUM( B2(:)*Uelem(:,cp) )
            fBU(3,cp) = SUM( B3(:)*Uelem(:,cp) ) * 2.0d0
            
            ! bending stress
            fHBU(1,cp) = h*h/12.0D0 * SUM( matH(:,1) * fBU(:,cp) )
            fHBU(2,cp) = h*h/12.0D0 * SUM( matH(:,2) * fBU(:,cp) )
            fHBU(3,cp) = h*h/12.0D0 * SUM( matH(:,3) * fBU(:,cp) )
            
            
            
            ! derivations of bending strain
            call cross(Uelem(:,cp),AI(2,:), UxA2(:))
            call cross(AI(1,:),Uelem(:,cp), A1xU(:))
            call cross(Uelem(:,cp),dAI1dxi(1,:), UxdA1d1(:))
            call cross(Uelem(:,cp),dAI2dxi(2,:), UxdA2d2(:))
            call cross(Uelem(:,cp),dAI1dxi(2,:), UxdA1d2(:))
            call cross(UxA2(:),AI(2,:), UA2xA2(:) )
            call cross(AI(1,:),UxA2(:), A1xUA2(:) )
            call cross(A1xU(:),AI(2,:), A1UxA2(:) )
            call cross(AI(1,:),A1xU(:), A1xA1U(:) )
            
            
            save1 = SUM( Uelem(:,cp)*(
     &           dRdxi(cp,1)*dA1d1_A2(:) + dRdxi(cp,2)*A1_dA1d1(:)
     &           + A3dA1d1*( dRdxi(cp,1)*AIxAJ(:,1)
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ) )/Area/Area
            save2 = SUM( Uelem(:,cp)*( 
     &           dRdxi(cp,1)*dA2d2_A2(:) + dRdxi(cp,2)*A1_dA2d2(:)
     &           + A3dA2d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ))/Area/Area
            save3 = SUM( Uelem(:,cp)*( 
     &           dRdxi(cp,1)*dA1d2_A2(:) + dRdxi(cp,2)*A1_dA1d2(:)
     &           + A3dA1d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ))/Area/Area * two
            
            save4 = one/Area * SUM( Uelem(:,cp)*(
     &           dRdxi(cp,1)*AIxAJ(:,1) + dRdxi(cp,2)*AIXAJ(:,2) ) )
            
            
            Do kk = 1,NNODE
               ! 1st term
               Do k = 1,3
               dA3dPe_U(k) = dJdPe(k,kk)*SUM( AI(3,:)*Uelem(:,cp) )
     &              + SUM( UxA2(:)*dA1dPe(:,k,kk) )
     &              + SUM( A1xU(:)*dA2dPe(:,k,kk) )
               Enddo
               dA3dPe_U(:) = -dA3dPe_U(:)/Area
               
               dfBdPeU(:,kk,1,cp) = -ddRddxi(cp,1)*dA3dPe_U(:)
               dfBdPeU(:,kk,2,cp) = -ddRddxi(cp,2)*dA3dPe_U(:)
               dfBdPeU(:,kk,3,cp) = -ddRddxi(cp,3)*dA3dPe_U(:)*two
               
               ! 2nd term dJ**(-1)dP
               dfBdPeU(:,kk,1,cp)=dfBdPeU(:,kk,1,cp) - dJdPe(:,kk)*save1
               dfBdPeU(:,kk,2,cp)=dfBdPeU(:,kk,2,cp) - dJdPe(:,kk)*save2
               dfBdPeU(:,kk,3,cp)=dfBdPeU(:,kk,3,cp) - dJdPe(:,kk)*save3
               
               ! 3rd term
               Do i = 1,3
               dfBdPeU(i,kk,1,cp) = dfBdPeU(i,kk,1,cp) + one/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA1d1Pe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA1d1Pe(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA1d1(:) * dA2dPe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA1d1(:) * dA1dPe(:,i,kk) ))
               dfBdPeU(i,kk,2,cp) = dfBdPeU(i,kk,2,cp) + one/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA2d2Pe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA2d2Pe(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA2d2(:) * dA2dPe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA2d2(:) * dA1dPe(:,i,kk) ))
               dfBdPeU(i,kk,3,cp) = dfBdPeU(i,kk,3,cp) + two/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA1d2Pe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA1d2Pe(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA1d2(:) * dA2dPe(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA1d2(:) * dA1dPe(:,i,kk) ))
               Enddo
               
               
               ! 4th term
               ! - scalar product
               dfBdPeU(:,kk,1,cp) = dfBdPeU(:,kk,1,cp) 
     &                            + save4*d_A3dA1d1_dPe(:,kk)
               dfBdPeU(:,kk,2,cp) = dfBdPeU(:,kk,2,cp) 
     &                            + save4*d_A3dA2d2_dPe(:,kk)
               dfBdPeU(:,kk,3,cp) = dfBdPeU(:,kk,3,cp) 
     &                            + save4*d_A3dA1d2_dPe(:,kk)*two
               
               ! - parenthesis
               dA3dPe_U(:) = zero
               Do k = 1,3
                  dA3dPe_U(k)
     &                 = dRdxi(cp,1)*( dJdPe(k,kk)*SUM(AI(3,:)*UxA2(:))
     &                   + SUM(dA1dPe(:,k,kk)*UA2xA2(:))
     &                   + SUM(dA2dPe(:,k,kk)*A1xUA2(:)) )
     &                 + dRdxi(cp,2)*( dJdPe(k,kk)*SUM( AI(3,:)*A1xU(:))
     &                   + SUM(dA1dPe(:,k,kk)*A1UxA2(:))
     &                   + SUM(dA2dPe(:,k,kk)*A1xA1U(:)) )
               Enddo
               dA3dPe_U(:) = -dA3dPe_U(:)/Area

               V(:) = zero
               Do k = 1,3
               V(k) = dRdxi(cp,1)*SUM(A3xU(:)*dA2dPe(:,k,kk))
     &              - dRdxi(cp,2)*SUM(A3xU(:)*dA1dPe(:,k,kk))
               Enddo
               V(:) = V(:) + dA3dPe_U(:)
               
               dfBdPeU(:,kk,1,cp)=dfBdPeU(:,kk,1,cp) + A3dA1d1/Area*V(:)
               dfBdPeU(:,kk,2,cp)=dfBdPeU(:,kk,2,cp) + A3dA2d2/Area*V(:)
               V(:) = V(:)*two
               dfBdPeU(:,kk,3,cp)=dfBdPeU(:,kk,3,cp) + A3dA1d2/Area*V(:)
            Enddo



            ! derivatives of contravariant metric coefficient
            dAAIdPe(:,:,:) = zero
            dAAEdPe(:,:,:) = zero
            Do k = 1,3
               dAAIdPe(k,1,1) = SUM( dA1dPe(:,k,cp)*AI(1,:) )*two
               dAAIdPe(k,2,2) = SUM( dA2dPe(:,k,cp)*AI(2,:) )*two
               dAAIdPe(k,1,2) = SUM( dA1dPe(:,k,cp)*AI(2,:) )
     &                        + SUM( dA2dPe(:,k,cp)*AI(1,:) )
               dAAIdPe(k,2,1) = dAAIdPe(k,1,2)
            Enddo
            Do j = 1,2
               Do i = 1,2
                  dAAEdPe(:,i,j) = 
     &                 - AAE(i,1)*AAE(1,j) * dAAIdPe(:,1,1)
     &                 - AAE(i,1)*AAE(2,j) * dAAIdPe(:,1,2)
     &                 - AAE(i,2)*AAE(1,j) * dAAIdPe(:,2,1)
     &                 - AAE(i,2)*AAE(2,j) * dAAIdPe(:,2,2)
               Enddo
            Enddo
            
            ! derivatives of material matrix
            dHdPe(:,cp,1,1) = two*dAAEdPe(:,1,1)*AAE(1,1)
            dHdPe(:,cp,2,2) = two*dAAEdPe(:,2,2)*AAE(2,2)
            dHdPe(:,cp,3,3) = 0.5d0*
     &           ((one-nu)*( dAAEdPe(:,1,1)*AAE(2,2) 
     &                       + AAE(1,1)*dAAEdPe(:,2,2) )
     &           +(one+nu)*( dAAEdPe(:,1,2)*AAE(1,2) * two ))
            dHdPe(:,cp,1,2) = 
     &           nu*( dAAEdPe(:,1,1)*AAE(2,2) + AAE(1,1)*dAAEdPe(:,2,2))
     &           + (one-nu) * two * dAAEdPe(:,1,2)*AAE(1,2)
            dHdPe(:,cp,1,3) = dAAEdPe(:,1,1)*AAE(1,2)
     &                     + AAE(1,1)*dAAEdPe(:,1,2)
            dHdPe(:,cp,2,3) = dAAEdPe(:,2,2)*AAE(1,2)
     &                     + AAE(2,2)*dAAEdPe(:,1,2)
            dHdPe(:,cp,2,1) = dHdPe(:,cp,1,2)
            dHdPe(:,cp,3,1) = dHdPe(:,cp,1,3)
            dHdPe(:,cp,3,2) = dHdPe(:,cp,2,3)
            
         Enddo
         
         dHdPe(:,:,:,:) = (h*coef) * dHdPe(:,:,:,:)
         
         
         

         
         
         
         
         
         
         
         
         
         
         
c     START LOOPS on CPs
         Do cpb = 1,NNODE
            
            Ub(:) = Uelem(:,cpb)
            Vb(:) = mBU(  :,cpb)
            Wb(:) = fBU(  :,cpb)
            Xb(:) = mHBU( :,cpb)
            Yb(:) = fHBU( :,cpb)
            
            dHdPe_mBU(:,:,:) = zero
            dHdPe_fBU(:,:,:) = zero
            Do ll = 1,3
               Do kk = 1,3
                  dHdPe_mBU(:,:,ll) = dHdPe_mBU(:,:,ll)
     &                              + dHdPe(:,:,kk,ll)*Vb(kk)
                  dHdPe_fBU(:,:,ll) = dHdPe_fBU(:,:,ll)
     &                              + dHdPe(:,:,kk,ll)*Wb(kk)*h*h/12.0D0
               Enddo
            Enddo
            
            
            
            Do cpa = 1,cpb
               
               Ua(:) = Uelem(:,cpa)
               Va(:) = mBU(  :,cpa)
               Wa(:) = fBU(  :,cpa)
               Xa(:) = mHBU( :,cpa)
               Ya(:) = fHBU( :,cpa)
               
               mBU_dHdPe_mBU(:,:) = zero
               fBU_dHdPe_fBU(:,:) = zero
               Do kk = 1,3
                  mBU_dHdPe_mBU(:,:)
     &                 = mBU_dHdPe_mBU(:,:) + Va(kk)*dHdPe_mBU(:,:,kk)
                  fBU_dHdPe_fBU(:,:)
     &                 = fBU_dHdPe_fBU(:,:) + Wa(kk)*dHdPe_fBU(:,:,kk)
               Enddo
               
               sKe = SUM( Va(:)*Xb(:) + Wa(:)*Yb(:) )
               
               
c     update sensitivity
               if (cpa==cpb) then
                  temp = -0.5d0 * DetJac
               else
                  temp = -1.0d0 * DetJac
               endif
               Do cp = 1,NNODE
                  
!     - 1. U.dBtdP.H.B.U
!     -- membrane
                  Do k = 1,3
                  gradC_elem(k,cp) = gradC_elem(k,cp) + (
     &               + dRdxi(cpa,1)*SUM( Ua(:)*dA1dPe(:,k,cp) ) * Xb(1)
     &               + dRdxi(cpa,2)*SUM( Ua(:)*dA2dPe(:,k,cp) ) * Xb(2)
     &               + dRdxi(cpa,2)*SUM( Ua(:)*dA1dPe(:,k,cp) ) * Xb(3)
     &               + dRdxi(cpa,1)*SUM( Ua(:)*dA2dPe(:,k,cp) ) * Xb(3)
     &               ) * area * temp
                  Enddo
!     -- flexion
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 +(dfBdPeU(:,cp,1,cpa)*Yb(1)
     &                 + dfBdPeU(:,cp,2,cpa)*Yb(2)
     &                 + dfBdPeU(:,cp,3,cpa)*Yb(3) ) * area * temp
                  
                  
!     - 2. U.Bt.dHdP.B.U
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 + ( mBU_dHdPe_mBU(:,cp) + fBU_dHdPe_fBU(:,cp) )
     &                   * area * temp
                  
                  
!     - 3. U.Bt.H.dBdP.U
!     -- membrane
                  Do k = 1,3
                  gradC_elem(k,cp) = gradC_elem(k,cp) + (
     &               + Xa(1) * SUM( Ub(:)*dA1dPe(:,k,cp) )*dRdxi(cpb,1)
     &               + Xa(2) * SUM( Ub(:)*dA2dPe(:,k,cp) )*dRdxi(cpb,2)
     &               + Xa(3) * SUM( Ub(:)*dA1dPe(:,k,cp) )*dRdxi(cpb,2)
     &               + Xa(3) * SUM( Ub(:)*dA2dPe(:,k,cp) )*dRdxi(cpb,1)
     &               ) * area * temp
                  Enddo
!     -- flexion
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 +(Ya(1)*dfBdPeU(:,cp,1,cpb)
     &                 + Ya(2)*dfBdPeU(:,cp,2,cpb)
     &                 + Ya(3)*dfBdPeU(:,cp,3,cpb) ) * area * temp
                  
                  
!     - 4. U.Bt.H.B.U.dJdP
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 + sKe * dJdPe(:,cp) * temp
                  
               Enddo
            Enddo
         Enddo
c     END LOOPS on CPs
         
         
         
c     Load vector
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            
            Do cpa = 1,NNODE
               if (KTypeDload == 0) then
                  temp = one
                  if (KNumFace == 6) then
                     temp =-one
                  endif
                  
                  call cross(Uelem(:,cpa),AI(2,:), UxA2(:))
                  call cross(AI(1,:),Uelem(:,cpa), A1xU(:))
                  
                  temp = temp*R(cpa)*ADLMAG(numLoad)*DetJac
                  
                  Do cp = 1,NNODE
                  Do k  = 1,3
                     gradC_elem(k,cp) = gradC_elem(k,cp)
     &                    -( SUM( UxA2(:)*dA1dPe(:,k,cp) ) 
     &                      +SUM( A1xU(:)*dA2dPe(:,k,cp) ) )* temp
                  Enddo
                  Enddo
                  
               elseif (KTypeDload == 1) then
                  VectNorm_U = Uelem(1,cpa)
               elseif (KTypeDload == 2) then
                  VectNorm_U = Uelem(2,cpa)
               elseif (KTypeDload == 3) then
                  VectNorm_U = Uelem(3,cpa)
               elseif (KTypeDload == 4) then
                  call norm(AI(1,:), MCRD, normV)
                  VectNorm_U = SUM( AI(1,:)*Uelem(:,cpa) )/normV
               elseif (KTypeDload == 5) then
                  call norm(AI(2,:), MCRD, normV)
                  VectNorm_U = SUM( AI(2,:)*Uelem(:,cpa) )/normV
               elseif (KTypeDload == 6) then
                  
                  temp = R(cpa)*Uelem(3,cpa)*ADLMAG(numLoad)*DetJac
                  Do cp = 1,NNODE
                  Do k  = 1,3
                     gradC_elem(k,cp) = gradC_elem(k,cp)
     &                    +(  AI(2,2)*dA1dPe(1,k,cp) 
     &                      - AI(1,2)*dA2dPe(1,k,cp))* temp           
     &                    -(  AI(2,1)*dA1dPe(2,k,cp) 
     &                      - AI(1,1)*dA2dPe(2,k,cp))* temp
                  Enddo
                  Enddo

               endif
               
            Enddo
            
         Enddo
c     End load traitement
         
      enddo
c     End of the loop on integration points on main surf
      
      
      
c     Build of element stiffness matrix done.
c
c     ..................................................................
c
c     Loop on load : find boundary loads
      
      
      End SUBROUTINE gradUELMAT30















































c     --
c     Construction de la sentitivite de  la matrice de rigite et du
c     second membre elementaire pour la formulation coque type 
c     Kirchhoff-Love immergee: 
c      - sensibilite par rapport aux PCs du mapping
c     --
      
      SUBROUTINE gradUELMAT30map(activeDir,activeElementMap,nb_elemMap,
     1     Uelem,NDOFEL,MCRD,NNODE,NNODEmap,nb_cp,JELEM,NBINT,COORDS,
     2     COORDSall,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,
     3     indDLoad,load_target_nbelem,JDLType,ADLMAG,gradC)
      
      use parameters
      use embeddedMapping
            
      Implicit None
      
c     Input arguments :
c     ---------------
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
 
      Double precision, intent(in) :: Uelem
      Integer, intent(in)       :: activeDir,activeElementMap,nb_elemMap
      dimension Uelem(3,NNODE),activeDir(3),activeElementMap(nb_elemMap)
      

c     Output variables :
c     ----------------
      Double precision, intent(inout) :: gradC
      dimension gradC(3,nb_cp)
            
      
c     Local variables :
c     ---------------
      
      Double precision :: gradC_elem
      dimension gradC_elem(MCRD,NNODEmap)

!     For gauss points
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(3,NBINT)
      
!     For nurbs basis functions
!     - embedded surface
      Double precision :: XI, R, dRdxi, ddRddxi, DetJac
      dimension XI(3), R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
      
!     Mapping
      ! - nurbs basis functions
      Double precision :: Rm, dRmdxi, ddRmddxi
      dimension Rm(NNODEmap), dRmdxi(NNODEmap,3), ddRmddxi(NNODEmap,6)
      ! - element infos
      Double precision :: COORDSmap
      dimension COORDSmap(MCRD,NNODEmap)
      Integer          :: sctr_map,isave
      dimension sctr_map(NNODEmap)
      
!     Composition Mapping+Surface
      Double precision :: dRRdxi, ddRRddxi
      dimension dRRdxi(NNODEmap,2),ddRRddxi(NNODEmap,3)
      
!     For curvilinear coordinate objects
      Double precision :: AI,dAI1dxi,dAI2dxi,AAE,AIxAJ,BI,dBI,VI,dVI
      dimension AI(3,3),dAI1dxi(2,3),dAI2dxi(2,3),AAE(2,2),AIxAJ(3,3),
     &     BI(3,2),dBI(3,3),VI(3,3),dVI(3,6)

!     For material matrix
      Double precision :: E, nu, h, matH, coef
      dimension matH(3,3)
      
      
!     For Membrane term
      Double precision :: Ua,Va,Xa,Ub,Vb,Xb,mBU,mHBU,dAAIdPm,dAAEdPm,
     &     dHdPm,dHdPm_mBU,mBU_dHdPm_mBU,dJdPm,sKE,dvol,Area,dA1dPm,
     &     dA2dPm
      dimension Ua(3),Va(3),Xa(3),Ub(3),Vb(3),Xb(3),mBU(3,NNODE),
     &     mHBU(3,NNODE),dAAIdPm(3,2,2),dAAEdPm(3,2,2),
     &     dHdPm(3,NNODEmap,3,3),dHdPm_mBU(3,NNODEmap,3),
     &     mBU_dHdPm_mBU(3,NNODEmap),dJdPm(3,NNODEmap),
     &     dA1dPm(3,3,NNODEmap),dA2dPm(3,3,NNODEmap)
      
!     For Bending term
      Double precision :: B1,B2,B3,dA1d1_A2,dA2d2_A2,dA1d2_A2,A1_dA1d1,
     &     A1_dA2d2,A1_dA1d2,A3dA1d1,A3dA2d2,A3dA1d2,fBU,fHBU,dHdPm_fBU,
     &     fBU_dHdPm_fBU,Wa,Ya,Wb,Yb,UxA2,A1xU,A3xU,dA3dPm_U,V,VxA2,
     &     A1xV,dfBdPmU,UxdA1d1,UxdA2d2,UxdA1d2,save1,save2,save3,save4,
     &     UA2xA2,A1xUA2,A1UxA2,A1xA1U,d_A3dA1d1_dPm,d_A3dA2d2_dPm,
     &     d_A3dA1d2_dPm,ddA1d1Pm,ddA2d2Pm,ddA1d2Pm
      dimension B1(3),B2(3),B3(3),dA1d1_A2(3),dA2d2_A2(3),dA1d2_A2(3),
     &     A1_dA1d1(3),A1_dA2d2(3),A1_dA1d2(3),fBU(3,NNODE),
     &     fHBU(3,NNODE),dHdPm_fBU(3,NNODEmap,3),
     &     fBU_dHdPm_fBU(3,NNODEmap),Wa(3),Ya(3),Wb(3),Yb(3),UxA2(3),
     &     A1xU(3),A3xU(3),dA3dPm_U(3),V(3),VxA2(3),A1xV(3),
     &     dfBdPmU(3,NNODEmap,3,NNODE),UxdA1d1(3),UxdA2d2(3),UxdA1d2(3),
     &     UA2xA2(3),A1xUA2(3),A1UxA2(3),A1xA1U(3),
     &     d_A3dA1d1_dPm(3,NNODEmap),d_A3dA2d2_dPm(3,NNODEmap),
     &     d_A3dA1d2_dPm(3,NNODEmap),ddA1d1Pm(3,3,NNODEmap),
     &     ddA2d2Pm(3,3,NNODEmap),ddA1d2Pm(3,3,NNODEmap)


!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      Double precision :: VectNorm_U,normV
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,j,k,kk,ll,cp,cpa,cpb, KTypeDload, KNumFace
      Double precision :: temp,temp1,temp2
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
            
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
         
         gradC_elem  = zero
         
c     PRELIMINARY QUANTITES
c     Computing NURBS basis functions and derivatives
         R      = zero
         dRdxi  = zero
         ddRddxi= zero
         DetJac = zero
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         DetJac = DetJac * GaussPdsCoord(1,n)
         
c     Find mapping parametric position
         XI(:)    = zero
         BI(:,:)  = zero
         dBI(:,:) = zero
         Do cp = 1,NNODE
            XI(:)    =  XI(:)   +         R(cp)*COORDS(:,cp)
            BI(:,1)  =  BI(:,1) +   dRdxi(cp,1)*COORDS(:,cp)
            BI(:,2)  =  BI(:,2) +   dRdxi(cp,2)*COORDS(:,cp)
            dBI(:,1) = dBI(:,1) + ddRddxi(cp,1)*COORDS(:,cp)
            dBI(:,2) = dBI(:,2) + ddRddxi(cp,2)*COORDS(:,cp)
            dBI(:,3) = dBI(:,3) + ddRddxi(cp,3)*COORDS(:,cp)
         Enddo
         
c     Computing NURBS basis functions and derivatives of the mapping
         ! get active element number
         call updateMapElementNumber(XI(:))
         
         
         
         IF (activeElementMap(current_map_elem) == 1) then



         call evalnurbs_mapping_w2ndDerv(XI(:),Rm(:),dRmdxi(:,:),
     &        ddRmddxi(:,:))

c     Composition of the basis functions
         dRRdxi(:,:) = zero
         Do i = 1,3
            dRRdxi(:,1) = dRRdxi(:,1) + BI(i,1)*dRmdxi(:,i)
            dRRdxi(:,2) = dRRdxi(:,2) + BI(i,2)*dRmdxi(:,i)
         Enddo
         
         ddRRddxi(:,:) = zero
         Do i = 1,3
            ddRRddxi(:,1) = ddRRddxi(:,1)
     &           + dBI(i,1)*dRmdxi(:,i) + BI(i,1)*BI(i,1)*ddRmddxi(:,i)
            ddRRddxi(:,2) = ddRRddxi(:,2) 
     &           + dBI(i,2)*dRmdxi(:,i) + BI(i,2)*BI(i,2)*ddRmddxi(:,i)
            ddRRddxi(:,3) = ddRRddxi(:,3) 
     &           + dBI(i,3)*dRmdxi(:,i) + BI(i,1)*BI(i,2)*ddRmddxi(:,i)
            Do j = i+1,3
               kk = i+j+1
               ddRRddxi(:,1) = ddRRddxi(:,1)
     &              + two*BI(i,1)*BI(j,1)*ddRmddxi(:,kk)
               ddRRddxi(:,2) = ddRRddxi(:,2)
     &              + two*BI(i,2)*BI(j,2)*ddRmddxi(:,kk)
               ddRRddxi(:,3) = ddRRddxi(:,3)
     &              +(BI(i,1)*BI(j,2) + BI(j,1)*BI(i,2))*ddRmddxi(:,kk)
            Enddo
         Enddo
         
c     Computing Curvilinear Coordinate objects
         ! extract COORDS
         If (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:,current_map_elem)
            
            Do cp = 1,NNODEmap
               COORDSmap(:,cp) = COORDSall(:,sctr_map(cp))
            Enddo
            
            isave = current_map_elem
         Endif
         
         call curvilinear(AI,dAI1dxi,dAI2dxi,AAE,Rm,dRRdxi,ddRRddxi,
     &        MCRD,NNODEmap,COORDSmap)
         
         VI(:,:) = zero
         Do k = 1,3
            Do cp = 1,NNODEmap
               VI(:,k) = VI(:,k) + dRmdxi(cp,k)*COORDSmap(:,cp)
            Enddo
         Enddo
         
         dVI(:,:) = zero
         Do k = 1,6
            Do cp = 1,NNODEmap
               dVI(:,k) = dVI(:,k) + ddRmddxi(cp,k)*COORDSmap(:,cp)
            Enddo
         Enddo
         
         call cross(AI(2,:),AI(3,:), AIxAJ(:,1))
         call cross(AI(3,:),AI(1,:), AIxAJ(:,2))
         call SurfElem(AI(1,:), AI(2,:), Area)
         
         call cross(dAI1dxi(1,:), AI(2,:), dA1d1_A2(:))
         call cross(dAI2dxi(2,:), AI(2,:), dA2d2_A2(:))
         call cross(dAI1dxi(2,:), AI(2,:), dA1d2_A2(:))
         call cross(AI(1,:), dAI1dxi(1,:), A1_dA1d1(:))
         call cross(AI(1,:), dAI2dxi(2,:), A1_dA2d2(:))
         call cross(AI(1,:), dAI1dxi(2,:), A1_dA1d2(:))
         call   dot(AI(3,:), dAI1dxi(1,:), A3dA1d1)
         call   dot(AI(3,:), dAI2dxi(2,:), A3dA2d2)
         call   dot(AI(3,:), dAI1dxi(2,:), A3dA1d2)
         
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
         matH(:,:) = (h*coef)*matH(:,:)
         
C     --
c     Preliminary derivatives
         ! - covariant basis vectors
         dA1dPm(:,:,:) = zero
         dA2dPm(:,:,:) = zero
         Do cp = 1,NNODEmap
            Do i = 1,3
               dA1dPm(i,i,cp) = dRRdxi(cp,1)
               dA2dPm(i,i,cp) = dRRdxi(cp,2)
            Enddo
         Enddo
         
         ddA1d1Pm(:,:,:) = zero
         ddA1d2Pm(:,:,:) = zero
         ddA2d2Pm(:,:,:) = zero
         Do cp = 1,NNODEmap
            Do i = 1,3
               ddA1d1Pm(i,i,cp) = ddRRddxi(cp,1)
               ddA2d2Pm(i,i,cp) = ddRRddxi(cp,2)
               ddA1d2Pm(i,i,cp) = ddRRddxi(cp,3)
            Enddo
         Enddo
         
         ! - Jacobian
         dJdPm(:,:) = zero
         Do cp = 1,NNODEmap
            Do i = 1,3
               dJdPm(i,cp) = SUM( AIxAJ(:,1)*dA1dPm(:,i,cp) )
     &                     + SUM( AIxAJ(:,2)*dA2dPm(:,i,cp) ) 
            Enddo
         Enddo
         
         ! - for bending term (a_{3}*a_{i,j})
         d_A3dA1d1_dPm(:,:) = zero
         d_A3dA2d2_dPm(:,:) = zero
         d_A3dA1d2_dPm(:,:) = zero
         Do cp = 1,NNODEmap
            Do i = 1,3
               temp = dJdPm(i,cp)*A3dA1d1
     &              + SUM( dA1d1_A2(:)*dA1dPm(:,i,cp) )
     &              + SUM( A1_dA1d1(:)*dA2dPm(:,i,cp) )
               
               d_A3dA1d1_dPm(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA1d1Pm(:,i,cp) )

               temp = dJdPm(i,cp)*A3dA2d2
     &              + SUM( dA2d2_A2(:)*dA1dPm(:,i,cp) )
     &              + SUM( A1_dA2d2(:)*dA2dPm(:,i,cp) )
               
               d_A3dA2d2_dPm(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA2d2Pm(:,i,cp) )
               
               temp = dJdPm(i,cp)*A3dA1d2
     &              + SUM( dA1d2_A2(:)*dA1dPm(:,i,cp) )
     &              + SUM( A1_dA1d2(:)*dA2dPm(:,i,cp) )
               
               d_A3dA1d2_dPm(i,cp) = -temp/Area
     &              + SUM( AI(3,:)*ddA1d2Pm(:,i,cp) )
            Enddo
         Enddo
         
         

c     Matrices and vectors
         mBU(:,:)  = zero
         mHBU(:,:) = zero
         fBU(:,:)  = zero
         fHBU(:,:) = zero
         dHdPm(:,:,:,:) = zero
         Do cp = 1,NNODE
            
            ! membrane strain
            temp1 = SUM( AI(1,:)*Uelem(:,cp) )
            temp2 = SUM( AI(2,:)*Uelem(:,cp) )
            mBU(1,cp) = dRdxi(cp,1) * temp1
            mBU(2,cp) = dRdxi(cp,2) * temp2
            mBU(3,cp) = dRdxi(cp,1) * temp2 + dRdxi(cp,2) * temp1
            
            ! membrane stress
            mHBU(1,cp) = SUM( matH(:,1) * mBU(:,cp) )
            mHBU(2,cp) = SUM( matH(:,2) * mBU(:,cp) )
            mHBU(3,cp) = SUM( matH(:,3) * mBU(:,cp) )
            
            ! bending strain
            B1(:) = -ddRddxi(cp,1)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA1d1_A2(:) + dRdxi(cp,2)*A1_dA1d1(:)
     &           + A3dA1d1*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            B2(:) = -ddRddxi(cp,2)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA2d2_A2(:) + dRdxi(cp,2)*A1_dA2d2(:)
     &           + A3dA2d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            B3(:) = -ddRddxi(cp,3)*AI(3,:) +
     &           ( dRdxi(cp,1)*dA1d2_A2(:) + dRdxi(cp,2)*A1_dA1d2(:)
     &           + A3dA1d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &                     + dRdxi(cp,2)*AIXAJ(:,2) )
     &           )/Area
            fBU(1,cp) = SUM( B1(:)*Uelem(:,cp) )
            fBU(2,cp) = SUM( B2(:)*Uelem(:,cp) )
            fBU(3,cp) = SUM( B3(:)*Uelem(:,cp) ) * 2.0d0
            
            ! bending stress
            fHBU(1,cp) = h*h/12.0D0 * SUM( matH(:,1) * fBU(:,cp) )
            fHBU(2,cp) = h*h/12.0D0 * SUM( matH(:,2) * fBU(:,cp) )
            fHBU(3,cp) = h*h/12.0D0 * SUM( matH(:,3) * fBU(:,cp) )
            
            
            
            ! derivations of bending strain
            call cross(Uelem(:,cp),AI(2,:), UxA2(:))
            call cross(AI(1,:),Uelem(:,cp), A1xU(:))
            call cross(Uelem(:,cp),dAI1dxi(1,:), UxdA1d1(:))
            call cross(Uelem(:,cp),dAI2dxi(2,:), UxdA2d2(:))
            call cross(Uelem(:,cp),dAI1dxi(2,:), UxdA1d2(:))
            call cross(UxA2(:),AI(2,:), UA2xA2(:) )
            call cross(AI(1,:),UxA2(:), A1xUA2(:) )
            call cross(A1xU(:),AI(2,:), A1UxA2(:) )
            call cross(AI(1,:),A1xU(:), A1xA1U(:) )
            
            
            save1 = SUM( Uelem(:,cp)*(
     &           dRdxi(cp,1)*dA1d1_A2(:) + dRdxi(cp,2)*A1_dA1d1(:)
     &           + A3dA1d1*( dRdxi(cp,1)*AIxAJ(:,1)
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ) )/Area/Area
            save2 = SUM( Uelem(:,cp)*( 
     &           dRdxi(cp,1)*dA2d2_A2(:) + dRdxi(cp,2)*A1_dA2d2(:)
     &           + A3dA2d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ))/Area/Area
            save3 = SUM( Uelem(:,cp)*( 
     &           dRdxi(cp,1)*dA1d2_A2(:) + dRdxi(cp,2)*A1_dA1d2(:)
     &           + A3dA1d2*( dRdxi(cp,1)*AIxAJ(:,1) 
     &             + dRdxi(cp,2)*AIXAJ(:,2)) ))/Area/Area * two
            
            save4 = one/Area * SUM( Uelem(:,cp)*(
     &           dRdxi(cp,1)*AIxAJ(:,1) + dRdxi(cp,2)*AIXAJ(:,2) ) )
            
            
            Do kk = 1,NNODEmap
               ! 1st term
               Do k = 1,3
               dA3dPm_U(k) = dJdPm(k,kk)*SUM( AI(3,:)*Uelem(:,cp) )
     &              + SUM( UxA2(:)*dA1dPm(:,k,kk) )
     &              + SUM( A1xU(:)*dA2dPm(:,k,kk) )
               Enddo
               dA3dPm_U(:) = -dA3dPm_U(:)/Area
               
               dfBdPmU(:,kk,1,cp) = -ddRddxi(cp,1)*dA3dPm_U(:)
               dfBdPmU(:,kk,2,cp) = -ddRddxi(cp,2)*dA3dPm_U(:)
               dfBdPmU(:,kk,3,cp) = -ddRddxi(cp,3)*dA3dPm_U(:)*two
               
               ! 2nd term dJ**(-1)dP
               dfBdPmU(:,kk,1,cp)=dfBdPmU(:,kk,1,cp) - dJdPm(:,kk)*save1
               dfBdPmU(:,kk,2,cp)=dfBdPmU(:,kk,2,cp) - dJdPm(:,kk)*save2
               dfBdPmU(:,kk,3,cp)=dfBdPmU(:,kk,3,cp) - dJdPm(:,kk)*save3
               
               ! 3rd term
               Do i = 1,3
               dfBdPmU(i,kk,1,cp) = dfBdPmU(i,kk,1,cp) + one/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA1d1Pm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA1d1Pm(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA1d1(:) * dA2dPm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA1d1(:) * dA1dPm(:,i,kk) ))
               dfBdPmU(i,kk,2,cp) = dfBdPmU(i,kk,2,cp) + one/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA2d2Pm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA2d2Pm(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA2d2(:) * dA2dPm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA2d2(:) * dA1dPm(:,i,kk) ))
               dfBdPmU(i,kk,3,cp) = dfBdPmU(i,kk,3,cp) + two/Area*(
     &              - dRdxi(cp,1)*SUM( UxA2(:) *  ddA1d2Pm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( A1xU(:) *  ddA1d2Pm(:,i,kk) )
     &              + dRdxi(cp,1)*SUM( UxdA1d2(:) * dA2dPm(:,i,kk) )
     &              - dRdxi(cp,2)*SUM( UxdA1d2(:) * dA1dPm(:,i,kk) ))
               Enddo
               
               
               ! 4th term
               ! - scalar product
               dfBdPmU(:,kk,1,cp) = dfBdPmU(:,kk,1,cp) 
     &                            + save4*d_A3dA1d1_dPm(:,kk)
               dfBdPmU(:,kk,2,cp) = dfBdPmU(:,kk,2,cp) 
     &                            + save4*d_A3dA2d2_dPm(:,kk)
               dfBdPmU(:,kk,3,cp) = dfBdPmU(:,kk,3,cp) 
     &                            + save4*d_A3dA1d2_dPm(:,kk)*two
               
               ! - parenthesis
               dA3dPm_U(:) = zero
               Do k = 1,3
                  dA3dPm_U(k)
     &                 = dRdxi(cp,1)*( dJdPm(k,kk)*SUM(AI(3,:)*UxA2(:))
     &                   + SUM(dA1dPm(:,k,kk)*UA2xA2(:))
     &                   + SUM(dA2dPm(:,k,kk)*A1xUA2(:)) )
     &                 + dRdxi(cp,2)*( dJdPm(k,kk)*SUM( AI(3,:)*A1xU(:))
     &                   + SUM(dA1dPm(:,k,kk)*A1UxA2(:))
     &                   + SUM(dA2dPm(:,k,kk)*A1xA1U(:)) )
               Enddo
               dA3dPm_U(:) = -dA3dPm_U(:)/Area

               V(:) = zero
               Do k = 1,3
               V(k) = dRdxi(cp,1)*SUM(A3xU(:)*dA2dPm(:,k,kk))
     &              - dRdxi(cp,2)*SUM(A3xU(:)*dA1dPm(:,k,kk))
               Enddo
               V(:) = V(:) + dA3dPm_U(:)
               
               dfBdPmU(:,kk,1,cp)=dfBdPmU(:,kk,1,cp) + A3dA1d1/Area*V(:)
               dfBdPmU(:,kk,2,cp)=dfBdPmU(:,kk,2,cp) + A3dA2d2/Area*V(:)
               V(:) = V(:)*two
               dfBdPmU(:,kk,3,cp)=dfBdPmU(:,kk,3,cp) + A3dA1d2/Area*V(:)
            Enddo
         Enddo
         
         
         
         Do kk = 1,NNODEmap

            ! derivatives of contravariant metric coefficient
            dAAIdPm(:,:,:) = zero
            dAAEdPm(:,:,:) = zero
            Do k = 1,3
               dAAIdPm(k,1,1) = SUM( dA1dPm(:,k,kk)*AI(1,:) )*two
               dAAIdPm(k,2,2) = SUM( dA2dPm(:,k,kk)*AI(2,:) )*two
               dAAIdPm(k,1,2) = SUM( dA1dPm(:,k,kk)*AI(2,:) )
     &                        + SUM( dA2dPm(:,k,kk)*AI(1,:) )
               dAAIdPm(k,2,1) = dAAIdPm(k,1,2)
            Enddo
            Do j = 1,2
               Do i = 1,2
                  dAAEdPm(:,i,j) = 
     &                 - AAE(i,1)*AAE(1,j) * dAAIdPm(:,1,1)
     &                 - AAE(i,1)*AAE(2,j) * dAAIdPm(:,1,2)
     &                 - AAE(i,2)*AAE(1,j) * dAAIdPm(:,2,1)
     &                 - AAE(i,2)*AAE(2,j) * dAAIdPm(:,2,2)
               Enddo
            Enddo
            
            ! derivatives of material matrix
            dHdPm(:,kk,1,1) = two*dAAEdPm(:,1,1)*AAE(1,1)
            dHdPm(:,kk,2,2) = two*dAAEdPm(:,2,2)*AAE(2,2)
            dHdPm(:,kk,3,3) = 0.5d0*
     &           ((one-nu)*( dAAEdPm(:,1,1)*AAE(2,2) 
     &                       + AAE(1,1)*dAAEdPm(:,2,2) )
     &           +(one+nu)*( dAAEdPm(:,1,2)*AAE(1,2) * two ))
            dHdPm(:,kk,1,2) = 
     &           nu*( dAAEdPm(:,1,1)*AAE(2,2) + AAE(1,1)*dAAEdPm(:,2,2))
     &           + (one-nu) * two * dAAEdPm(:,1,2)*AAE(1,2)
            dHdPm(:,kk,1,3) = dAAEdPm(:,1,1)*AAE(1,2)
     &                     + AAE(1,1)*dAAEdPm(:,1,2)
            dHdPm(:,kk,2,3) = dAAEdPm(:,2,2)*AAE(1,2)
     &                     + AAE(2,2)*dAAEdPm(:,1,2)
            dHdPm(:,kk,2,1) = dHdPm(:,kk,1,2)
            dHdPm(:,kk,3,1) = dHdPm(:,kk,1,3)
            dHdPm(:,kk,3,2) = dHdPm(:,kk,2,3)
            
         Enddo
         
         dHdPm(:,:,:,:) = (h*coef) * dHdPm(:,:,:,:)
         
         
         

         
         
         
         
         
         
         
         
         
         
         
c     START LOOPS on CPs
         Do cpb = 1,NNODE
            
            Ub(:) = Uelem(:,cpb)
            Vb(:) = mBU(  :,cpb)
            Wb(:) = fBU(  :,cpb)
            Xb(:) = mHBU( :,cpb)
            Yb(:) = fHBU( :,cpb)
            
            dHdPm_mBU(:,:,:) = zero
            dHdPm_fBU(:,:,:) = zero
            Do ll = 1,3
               Do kk = 1,3
                  dHdPm_mBU(:,:,ll) = dHdPm_mBU(:,:,ll)
     &                              + dHdPm(:,:,kk,ll)*Vb(kk)
                  dHdPm_fBU(:,:,ll) = dHdPm_fBU(:,:,ll)
     &                              + dHdPm(:,:,kk,ll)*Wb(kk)*h*h/12.0D0
               Enddo
            Enddo
            
            
            
            Do cpa = 1,cpb
               
               Ua(:) = Uelem(:,cpa)
               Va(:) = mBU(  :,cpa)
               Wa(:) = fBU(  :,cpa)
               Xa(:) = mHBU( :,cpa)
               Ya(:) = fHBU( :,cpa)
               
               mBU_dHdPm_mBU(:,:) = zero
               fBU_dHdPm_fBU(:,:) = zero
               Do kk = 1,3
                  mBU_dHdPm_mBU(:,:)
     &                 = mBU_dHdPm_mBU(:,:) + Va(kk)*dHdPm_mBU(:,:,kk)
                  fBU_dHdPm_fBU(:,:)
     &                 = fBU_dHdPm_fBU(:,:) + Wa(kk)*dHdPm_fBU(:,:,kk)
               Enddo
               
               sKe = SUM( Va(:)*Xb(:) + Wa(:)*Yb(:) )
               
               
c     update sensitivity
               if (cpa==cpb) then
                  temp = -0.5d0 * DetJac
               else
                  temp = -1.0d0 * DetJac
               endif
               Do cp = 1,NNODEmap
                  
!     - 1. U.dBtdP.H.B.U
!     -- membrane
                  Do k = 1,3
                  gradC_elem(k,cp) = gradC_elem(k,cp) + (
     &               + dRdxi(cpa,1)*SUM( Ua(:)*dA1dPm(:,k,cp) ) * Xb(1)
     &               + dRdxi(cpa,2)*SUM( Ua(:)*dA2dPm(:,k,cp) ) * Xb(2)
     &               + dRdxi(cpa,2)*SUM( Ua(:)*dA1dPm(:,k,cp) ) * Xb(3)
     &               + dRdxi(cpa,1)*SUM( Ua(:)*dA2dPm(:,k,cp) ) * Xb(3)
     &               ) * area * temp
                  Enddo
!     -- flexion
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 +(dfBdPmU(:,cp,1,cpa)*Yb(1)
     &                 + dfBdPmU(:,cp,2,cpa)*Yb(2)
     &                 + dfBdPmU(:,cp,3,cpa)*Yb(3) ) * area * temp
                  
                  
!     - 2. U.Bt.dHdP.B.U
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 + ( mBU_dHdPm_mBU(:,cp) + fBU_dHdPm_fBU(:,cp) )
     &                   * area * temp
                  
                  
!     - 3. U.Bt.H.dBdP.U
!     -- membrane
                  Do k = 1,3
                  gradC_elem(k,cp) = gradC_elem(k,cp) + (
     &               + Xa(1) * SUM( Ub(:)*dA1dPm(:,k,cp) )*dRdxi(cpb,1)
     &               + Xa(2) * SUM( Ub(:)*dA2dPm(:,k,cp) )*dRdxi(cpb,2)
     &               + Xa(3) * SUM( Ub(:)*dA1dPm(:,k,cp) )*dRdxi(cpb,2)
     &               + Xa(3) * SUM( Ub(:)*dA2dPm(:,k,cp) )*dRdxi(cpb,1)
     &               ) * area * temp
                  Enddo
!     -- flexion
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 +(Ya(1)*dfBdPmU(:,cp,1,cpb)
     &                 + Ya(2)*dfBdPmU(:,cp,2,cpb)
     &                 + Ya(3)*dfBdPmU(:,cp,3,cpb) ) * area * temp
                  
                  
!     - 4. U.Bt.H.B.U.dJdP
                  gradC_elem(:,cp) = gradC_elem(:,cp)
     &                 + sKe * dJdPm(:,cp) * temp
                  
               Enddo
            Enddo
         Enddo
c     END LOOPS on CPs
         
         
         
c     Load vector
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            
            Do cpa = 1,NNODE
               if (KTypeDload == 0) then
                  temp = one
                  if (KNumFace == 6) then
                     temp =-one
                  endif
                  
                  call cross(Uelem(:,cpa),AI(2,:), UxA2(:))
                  call cross(AI(1,:),Uelem(:,cpa), A1xU(:))
                  
                  temp = temp*R(cpa)*ADLMAG(numLoad)*DetJac
                  
                  Do cp = 1,NNODEmap
                  Do k  = 1,3
                     gradC_elem(k,cp) = gradC_elem(k,cp)
     &                    -( SUM( UxA2(:)*dA1dPm(:,k,cp) ) 
     &                      +SUM( A1xU(:)*dA2dPm(:,k,cp) ) )* temp
                  Enddo
                  Enddo
                  
               elseif (KTypeDload == 1) then
                  VectNorm_U = Uelem(1,cpa)
               elseif (KTypeDload == 2) then
                  VectNorm_U = Uelem(2,cpa)
               elseif (KTypeDload == 3) then
                  VectNorm_U = Uelem(3,cpa)
               elseif (KTypeDload == 4) then
                  call norm(AI(1,:), MCRD, normV)
                  VectNorm_U = SUM( AI(1,:)*Uelem(:,cpa) )/normV
               elseif (KTypeDload == 5) then
                  call norm(AI(2,:), MCRD, normV)
                  VectNorm_U = SUM( AI(2,:)*Uelem(:,cpa) )/normV
               elseif (KTypeDload == 6) then
                  
                  temp = R(cpa)*Uelem(3,cpa)*ADLMAG(numLoad)*DetJac
                  Do cp = 1,NNODEmap
                  Do k  = 1,3
                     gradC_elem(k,cp) = gradC_elem(k,cp)
     &                    +(  AI(2,2)*dA1dPm(1,k,cp) 
     &                      - AI(1,2)*dA2dPm(1,k,cp))* temp           
     &                    -(  AI(2,1)*dA1dPm(2,k,cp) 
     &                      - AI(1,1)*dA2dPm(2,k,cp))* temp
                  Enddo
                  Enddo

               endif
               
            Enddo
            
         Enddo
c     End load traitement
         
         
c     Update global gradient
         Do cp = 1,NNODEmap
            kk = sctr_map(cp)
            gradC(:,kk) = gradC(:,kk) + gradC_elem(:,cp)
         Enddo


         ENDIF ! test current mapping elem is active

      Enddo
c     End of the loop on integration points on main surf
      
      
      
c     Build of element stiffness matrix done.
c
c     ..................................................................
c
c     Loop on load : find boundary loads
      
      
      End SUBROUTINE gradUELMAT30map






















      
      SUBROUTINE gradUELMAT30gloAdj(
     1     activeDir,activeElementMap,nb_elemMap,Uelem,UAelem,
     2     NADJ,NDOFEL,MCRD,NNODE,NNODEmap,nb_cp,JELEM,NBINT,COORDS,
     3     COORDSall,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,nb_load,
     4     indDLoad,load_target_nbelem,JDLType,ADLMAG,computeWint,
     5     computeWext,gradWint,gradWext)
      
      use parameters
      use embeddedMapping
      
      implicit none
      
c     Input arguments :
c     ---------------
!     Embedded surface
      Integer, intent(in) :: NADJ,NDOFEL,MCRD,NNODE,JELEM,NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES, PROPS
      dimension COORDS(3,NNODE),MATERIAL_PROPERTIES(2),PROPS(JPROPS)
      
!     Global Mapping
      Integer, intent(in) :: NNODEmap, nb_cp
      Double precision, intent(in) :: COORDSall
      dimension COORDSall(3,nb_cp)
      
      Integer, intent(in) :: indDLoad,load_target_nbelem,JDLType,nb_load
      Double precision, intent(in) :: ADLMAG
      dimension ADLMAG(nb_load),indDLoad(SUM(load_target_nbelem)),
     &     load_target_nbelem(nb_load),JDLType(nb_load)
      
      Double precision, intent(in) :: Uelem,UAelem
      Integer, intent(in)          :: activeDir,activeElementMap,
     &     nb_elemMap
      dimension Uelem(3,NNODE),UAelem(3,NNODE,NADJ),activeDir(3),
     &     activeElementMap(nb_elemMap)

      Logical, intent(in) :: computeWint,computeWext
      
c     Output variables :
c     ----------------
      Double precision, intent(inout) :: gradWint,gradWext
      dimension gradWint(NADJ,3,nb_cp),gradWext(NADJ,3,nb_cp)
      
      
c     Local variables :
c     ---------------
      
!     For gauss points
      Integer :: NbPtInt
      Double precision :: GaussPdsCoord
      dimension GaussPdsCoord(3,NBINT)
      
!     For nurbs basis functions
!     - embedded surface
      Double precision :: R, dRdxi, ddRddxi, DetJac
      dimension R(NNODE), dRdxi(NNODE,2), ddRddxi(NNODE,3)
!     - global mapping
      Double precision :: Rm, dRmdxi, ddRmddxi
      dimension Rm(NNODEmap), dRmdxi(NNODEmap,3), ddRmddxi(NNODEmap,6)
!     - composition 
      Double precision :: dRRdxi, ddRRddxi
      dimension dRRdxi(NNODEmap,2),ddRRddxi(NNODEmap,3)
      
!     For curvilinear coordinate objects
!     - embedded surface
      Double precision :: XI,BI,dBI
      dimension XI(3),BI(3,2),dBI(3,3)
!     - global mapping
      Double precision :: COORDSmap
      dimension COORDSmap(3,NNODEmap)
      Integer          :: sctr_map,isave
      dimension sctr_map(NNODEmap)
!     - composition
      Double precision :: AI, dAIdxi,AAI,AAE, AIxAJ, Area, Det
      dimension AI(3,3), dAIdxi(3,3), AAI(3,3), AAE(3,3),AIxAJ(3,3)

!     For material matrix
      Integer :: voigt
      dimension voigt(3,2)
      Double precision :: E, nu, lambda, mu, h, matH, coef
      dimension matH(3,3)
      
!     For disp/strain/stress fields (state and adjoint)
      Double precision :: dUdxi,dUAdxi,ddUddxi,ddUAddxi,
     &     strainMem,stressMem,strainMemAdj,stressMemAdj,
     &     strainBnd,stressBnd,strainBndAdj,stressBndAdj,
     &     work,coef1,coef2,UA
      dimension dUdxi(3,2),dUAdxi(3,2,NADJ),ddUddxi(3,3),
     &     ddUAddxi(3,3,NADJ),strainMem(3),stressMem(3),
     &     strainMemAdj(3,nadj),stressMemAdj(3,nadj),strainBnd(3),
     &     stressBnd(3),strainBndAdj(3,nadj),stressBndAdj(3,nadj),
     &     work(nadj),UA(3,nadj)

!     For derivatives (membrane)
      Double precision :: dAAIdP,dAAEdP,dJdP,dEAdP_N,dEdP_NA,dCdP,
     &     dNdP_EA
      dimension dAAIdP(3,2,2),dAAEdP(3,2,2),dJdP(3),
     &     dEAdP_N(3,nadj),dEdP_NA(3,nadj),dCdP(3),dNdP_EA(3,nadj)

!     For derivatives (bending)
      Double precision :: dKAdP_M,dKdP_MA,dMdP_KA,temp,vect,vdA3dP,
     &     vectsave
      dimension dKAdP_M(3,nadj),dKdP_MA(3,nadj),dMdP_KA(3,nadj),vect(3),
     &     vdA3dP(3),vectsave(3)
      
!     Elementary sensitivities
      Double precision :: gradWint_elem,gradWext_elem
      dimension gradWint_elem(NADJ,3,NNODEmap),
     &          gradWext_elem(NADJ,3,NNODEmap)
      
!     For loads
      Integer :: nb_load_bnd,nb_load_srf,ind_load_loc,numLoad
      dimension ind_load_loc(nb_load)
      Double precision :: VectNorm,normV,UA_V,UA_dVdP
      dimension VectNorm(3),UA_V(nadj),UA_dVdP(3)
      
!      For loops
      Integer ntens
      Integer n,k1,k2,i,j,k,l,iA,ij,kl,kk,ll,cp, KTypeDload, KNumFace
      Double precision :: temp1,temp2
      
!     External functions
      Double precision, external :: DOTPROD, SCATRIPLEPROD
      INTERFACE
         FUNCTION CROSSPROD(u,v)
         IMPLICIT NONE
         double precision, dimension(3) :: CROSSPROD
         double precision, dimension(3), intent(in) :: u,v
         END FUNCTION CROSSPROD
      END INTERFACE 
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      ntens = 3                 ! size of stiffness tensor
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      if (NbPtInt*NbPtInt < NBINT) NbPtInt = NbPtInt + 1
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
c     Material behaviour
      h = PROPS(3)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      mu     = E/two/(one+nu)
      lambda = E*nu/(one+nu)/(one-two*nu)
      lambda = two*lambda*mu/(lambda+two*mu)
      
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
      isave = 0
      
c     Loop on integration points on main surface
      Do n = 1,NBINT
         
c     PRELIMINARY QUANTITES
c     Computing NURBS basis functions and derivatives
         R      = zero
         dRdxi  = zero
         ddRddxi= zero
         DetJac = zero
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,GaussPdsCoord(2:,n))
         DetJac = DetJac * GaussPdsCoord(1,n)
         
c     Find mapping parametric position
         XI(:)    = zero
         BI(:,:)  = zero
         dBI(:,:) = zero
         Do cp = 1,NNODE
            XI(:)    =  XI(:)   +         R(cp)*COORDS(:,cp)
            BI(:,1)  =  BI(:,1) +   dRdxi(cp,1)*COORDS(:,cp)
            BI(:,2)  =  BI(:,2) +   dRdxi(cp,2)*COORDS(:,cp)
            dBI(:,1) = dBI(:,1) + ddRddxi(cp,1)*COORDS(:,cp)
            dBI(:,2) = dBI(:,2) + ddRddxi(cp,2)*COORDS(:,cp)
            dBI(:,3) = dBI(:,3) + ddRddxi(cp,3)*COORDS(:,cp)
         Enddo
         
c     Go to global mapping level
c     - basis functions
         call updateMapElementNumber(XI(:))
         IF (activeElementMap(current_map_elem) == 1) then
         call evalnurbs_mapping_w2ndDerv(XI(:),Rm(:),dRmdxi(:,:),
     &        ddRmddxi(:,:))
c     - extract COORDS
         If (isave /= current_map_elem) then
            sctr_map(:) = IEN_map(:,current_map_elem)
            Do cp = 1,NNODEmap
               COORDSmap(:,cp) = COORDSall(:,sctr_map(cp))
            Enddo
            isave = current_map_elem
         Endif
         
c     Composition
         dRRdxi(:,:) = zero
         Do i = 1,3
            dRRdxi(:,1) = dRRdxi(:,1) + BI(i,1)*dRmdxi(:,i)
            dRRdxi(:,2) = dRRdxi(:,2) + BI(i,2)*dRmdxi(:,i)
         Enddo
         
         ddRRddxi(:,:) = zero
         Do i = 1,3
            ddRRddxi(:,1) = ddRRddxi(:,1)
     &           + dBI(i,1)*dRmdxi(:,i) + BI(i,1)*BI(i,1)*ddRmddxi(:,i)
            ddRRddxi(:,2) = ddRRddxi(:,2) 
     &           + dBI(i,2)*dRmdxi(:,i) + BI(i,2)*BI(i,2)*ddRmddxi(:,i)
            ddRRddxi(:,3) = ddRRddxi(:,3) 
     &           + dBI(i,3)*dRmdxi(:,i) + BI(i,1)*BI(i,2)*ddRmddxi(:,i)
            Do j = i+1,3
               kk = i+j+1
               ddRRddxi(:,1) = ddRRddxi(:,1)
     &              + two*BI(i,1)*BI(j,1)*ddRmddxi(:,kk)
               ddRRddxi(:,2) = ddRRddxi(:,2)
     &              + two*BI(i,2)*BI(j,2)*ddRmddxi(:,kk)
               ddRRddxi(:,3) = ddRRddxi(:,3)
     &              +(BI(i,1)*BI(j,2) + BI(j,1)*BI(i,2))*ddRmddxi(:,kk)
            Enddo
         Enddo
         
c     Computing Curvilinear Coordinate objects
c     - Covariant basis vectors and derivatives
         AI(:,:) = zero
         Do i = 1,2
            Do cp = 1,NNODEmap
               AI(:,i) = AI(:,i) + dRRdxi(cp,i)*COORDSmap(:,cp)
            Enddo
         Enddo
         
         call cross(AI(:,1),AI(:,2), AIxAJ(:,3))
         call norm( AIxAJ(:,3),3, Area)
         AI(:,3) = AIxAJ(:,3)/Area
         call cross(AI(:,2),AI(:,3), AIxAJ(:,1))
         call cross(AI(:,3),AI(:,1), AIxAJ(:,2))
         
         dAIdxi(:,:) = zero
         Do i = 1,ntens
            Do cp = 1,NNODEmap
               dAIdxi(:,i) = dAIdxi(:,i) +ddRRddxi(cp,i)*COORDSmap(:,cp)
            Enddo
         Enddo
c     - Computing metrics
         AAI(:,:) = zero
         AAE(:,:) = zero
         AAI(3,3) = one
         AAE(3,3) = one
         Do j = 1,2
            Do i = 1,2
               call dot(AI(:,i), AI(:,j), AAI(i,j))
            Enddo
         Enddo
         call MatrixInv(AAE(:2,:2), AAI(:2,:2), det, 2)
         
c     Computing material matrix
         voigt(:,1) = (/ 1,2,1 /)
         voigt(:,2) = (/ 1,2,2 /)
         
         matH(:,:) = zero
         Do kl = 1,ntens
            k=voigt(kl,1); l=voigt(kl,2)
            Do ij = 1,kl
               i=voigt(ij,1); j=voigt(ij,2)
               matH(ij,kl) = lambda*AAE(i,j)*AAE(k,l)
     &              + mu*(AAE(i,k)*AAE(j,l) + AAE(i,l)*AAE(j,k))
            Enddo
         Enddo
         Do kl = 2,ntens
            matH(kl:,kl-1) = matH(kl-1,kl:)
         Enddo
         matH(:,:) = h*matH(:,:)
         
c     Computing disp and adjoint derivatives
c     - 1st derivatives
         dUdxi(:,:) = zero
         Do i = 1,2
            Do cp = 1,NNODE
               dUdxi(:,i) = dUdxi(:,i) + dRdxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         
         dUAdxi(:,:,:) = zero
         Do iA = 1,nadj
         Do i = 1,2
            Do cp = 1,NNODE
               dUAdxi(:,i,iA) 
     &              = dUAdxi(:,i,iA) + dRdxi(cp,i)*UAelem(:,cp,iA)
            Enddo
         Enddo
         Enddo
c     - 2nd derivatives
         ddUddxi(:,:) = zero
         Do i = 1,3
            Do cp = 1,NNODE
               ddUddxi(:,i) = ddUddxi(:,i) + ddRddxi(cp,i)*Uelem(:,cp)
            Enddo
         Enddo
         
         ddUAddxi(:,:,:) = zero
         Do iA = 1,nadj
         Do i = 1,3
            Do cp = 1,NNODE
               ddUAddxi(:,i,iA) 
     &              = ddUAddxi(:,i,iA) + ddRddxi(cp,i)*UAelem(:,cp,iA)
            Enddo
         Enddo
         Enddo


         If (computeWint) then
c     Computing strain and stress (state and adjoint)
c     - Membrane
         strainMem(:) = zero
         stressMem(:) = zero
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUdxi(:,i), strainMem(ij))
            Else
               call dot(AI(:,i),dUdxi(:,j), coef1)
               call dot(AI(:,j),dUdxi(:,i), coef2)
               strainMem(ij) = coef1 + coef2
            Endif
         Enddo
         call MulVect(matH,strainMem,stressMem,ntens,ntens)
         
         strainMemAdj(:,:) = zero
         stressMemAdj(:,:) = zero
         Do iA = 1,nadj
         Do ij = 1,ntens
            i=voigt(ij,1); j=voigt(ij,2)
            If (i==j) then
               call dot(AI(:,i),dUAdxi(:,i,iA), strainMemAdj(ij,iA))
            Else
               call dot(AI(:,i),dUAdxi(:,j,iA), coef1)
               call dot(AI(:,j),dUAdxi(:,i,iA), coef2)
               strainMemAdj(ij,iA) = coef1 + coef2
            Endif
         Enddo 
         call MulVect(matH,strainMemAdj(:,iA),stressMemAdj(:,iA),ntens,
     &        ntens)       
         Enddo
         
c     - Bending
         matH(:,:) = h**two/12.0D0 * matH(:,:)
         
         strainBnd(:) = zero
         stressBnd(:) = zero
         temp = SUM(dUdxi(:,1)*AIxAJ(:,1) + dUdxi(:,2)*AIxAJ(:,2))/Area
         Do ij = 1,ntens            
            strainBnd(ij) =
     &           - SUM( ddUddxi(:,ij)*AI(:,3) )
     &           + one/Area*( 
     &             ScaTripleProd(dUdxi(:,1), dAIdxi(:,ij), AI(:,2))
     &             + ScaTripleProd(dUdxi(:,2), AI(:,1), dAIdxi(:,ij)) )
     &           + SUM(AI(:,3)*dAIdxi(:,ij))*temp
         Enddo
         strainBnd(3) = two*strainBnd(3)
         call MulVect(matH,strainBnd,stressBnd,ntens,ntens)
         
         strainBndAdj(:,:) = zero
         stressBndAdj(:,:) = zero
         Do iA = 1,nadj
         temp = SUM(dUAdxi(:,1,iA)*AIxAJ(:,1)+dUAdxi(:,2,iA)*AIxAJ(:,2))
     &           /Area
         Do ij = 1,ntens
            strainBndAdj(ij,iA) =
     &           - SUM( ddUAddxi(:,ij,iA)*AI(:,3) )
     &           + one/Area*( 
     &             ScaTripleProd(dUAdxi(:,1,iA), dAIdxi(:,ij), AI(:,2))
     &           + ScaTripleProd(dUAdxi(:,2,iA), AI(:,1), dAIdxi(:,ij)))
     &           + SUM(AI(:,3)*dAIdxi(:,ij))*temp
         Enddo 
         strainBndAdj(3,iA) = two*strainBndAdj(3,iA)
         call MulVect(matH,strainBndAdj(:,iA),stressBndAdj(:,iA),ntens,
     &        ntens)       
         Enddo
         
         
c     Computing local work
         work(:)  = zero
         Do ij = 1,ntens
            work(:) = work(:) 
     &           + strainMemAdj(ij,:)*stressMem(ij) 
     &           + strainBndAdj(ij,:)*stressBnd(ij)
         Enddo

         Endif  ! test computeWint is True
         
c     --
c     Derivatives
         If (computeWint) then
         gradWint_elem  = zero
         Do cp = 1,NNODEmap
            
            ! 1. derivatives of the jacobian
            dJdP(:) = zero
            dJdP(:) = AIxAJ(:,1)*dRRdxi(cp,1) + AIxAJ(:,2)*dRRdxi(cp,2)
            Do iA = 1,nadj
               gradWint_elem(iA,:,cp)
     &              = gradWint_elem(iA,:,cp) - work(iA)*dJdP(:)*detJac
            Enddo
            
            
            ! 2. derivatives of the adjoint membrane strain 
            !    (with dble prod. by membrane stress)
            dEAdP_N(:,:) = zero
            Do iA = 1,nadj
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               If (i==j) then
                  dEAdP_N(:,iA) = dEAdP_N(:,iA) 
     &                 + stressMem(ij)*dUAdxi(:,i,iA)*dRRdxi(cp,i)
               Else
                  dEAdP_N(:,iA) = dEAdP_N(:,iA) 
     &                 + stressMem(ij)*dUAdxi(:,i,iA)*dRRdxi(cp,j)
     &                 + stressMem(ij)*dUAdxi(:,j,iA)*dRRdxi(cp,i)
               Endif
            Enddo
            gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &           - dEAdP_N(:,iA)*Area*detJac
            Enddo
            
            
            ! 3. derivatives of the membrane stress 
            !    (with dble prod. by adjoint membrane strain)
            dNdP_EA(:,:) = zero
            !    - derivatives of covariant metrics
            dAAIdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  dAAIdP(:,i,j) = dRRdxi(cp,i)*AI(:,j)
     &                 + AI(:,i)*dRRdxi(cp,j)
               Enddo
            Enddo
            !    - derivatives of contravariant metrics
            dAAEdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  Do l = 1,2
                     Do k = 1,2
                        dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                       - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
            !    - first subterm (material dot derivative strain)
            dEdP_NA(:,:) = zero
            Do iA = 1,nadj
            Do ij = 1,ntens
               i=voigt(ij,1); j=voigt(ij,2)
               If (i==j) then
                  dEdP_NA(:,iA) = dEdP_NA(:,iA) 
     &                 + stressMemAdj(ij,iA)*dUdxi(:,i)*dRRdxi(cp,i)
               Else
                  dEdP_NA(:,iA) = dEdP_NA(:,iA) 
     &                 + stressMemAdj(ij,iA)*dUdxi(:,i)*dRRdxi(cp,j)
     &                 + stressMemAdj(ij,iA)*dUdxi(:,j)*dRRdxi(cp,i)
               Endif
            Enddo
            Enddo
            dNdP_EA(:,:) = dEdP_NA(:,:)
            
            !    - second subterm (derivative material tensor)
            Do kl = 1,ntens
               k=voigt(kl,1); l=voigt(kl,2)
               
               Do ij = 1,ntens
                  i=voigt(ij,1); j=voigt(ij,2)
                  
                  dCdP(:) =
     &                   lambda*dAAEdP(:,i,j)*AAE(k,l)
     &                 + lambda*AAE(i,j)*dAAEdP(:,k,l)
     &                 + mu*dAAEdP(:,i,k)*AAE(j,l)
     &                 + mu*AAE(i,k)*dAAEdP(:,j,l)
     &                 + mu*dAAEdP(:,i,l)*AAE(j,k)
     &                 + mu*AAE(i,l)*dAAEdP(:,j,k)
                  
                  Do iA = 1,nadj
                     dNdP_EA(:,iA) = dNdP_EA(:,iA) 
     &                    + h*dCdP(:)*strainMem(ij)*strainMemAdj(kl,iA)
                  Enddo
               Enddo
            Enddo

            Do iA = 1,nadj
               gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &              - dNdP_EA(:,iA)*Area*detJac
            Enddo
            
            
            ! 4. derivatives of the adjoint bending strain 
            !    (with dble prod. by bending stress)
            dKAdP_M(:,:) = zero
            Do iA = 1,nadj
            temp = SUM(dUAdxi(:,1,iA)*AIxAJ(:,1))
     &              + SUM(dUAdxi(:,2,iA)*AIxAJ(:,2))
            
            vect(:)   = CROSSPROD(dUAdxi(:,1,iA),AI(:,2)) 
     &                + CROSSPROD(AI(:,1),dUAdxi(:,2,iA))
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            vectsave(:) = vdA3dP(:)
            
            coef = one
            Do ij = 1,ntens
               
               if (ij==3) coef=two
               
               ! 1st term
               vect(:)   = ddUAddxi(:,ij,iA)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               dKAdP_M(:,iA) =dKAdP_M(:,iA)-vdA3dP(:)*stressBnd(ij)*coef
               
               ! 2nd term
               vect(:) =  
     &              -dJdP(:)/Area/Area*(
     &               ScaTripleProd(dUAdxi(:,1,iA),dAIdxi(:,ij),AI(:,2))
     &              +ScaTripleProd(dUAdxi(:,2,iA),AI(:,1),dAIdxi(:,ij)))
     &              +one/Area*(
     &               CROSSPROD(dUAdxi(:,1,iA),dAIdxi(:,ij))*dRRdxi(cp,2)
     &               -CROSSPROD(dUAdxi(:,1,iA),AI(:,2))*ddRRddxi(cp,ij))
     &              +one/Area*(
     &               CROSSPROD(dAIdxi(:,ij),dUAdxi(:,2,iA))*dRRdxi(cp,1)
     &               -CROSSPROD(AI(:,1),dUAdxi(:,2,iA))*ddRRddxi(cp,ij))
               
               dKAdP_M(:,iA) = dKAdP_M(:,iA) +vect(:)*stressBnd(ij)*coef
               
               ! 3rd term
               vect(:)   = dAIdxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               vect(:) = 
     &             -DOTPROD(AI(:,3),dAIdxi(:,ij))*dJdP(:)/Area/Area*temp
     &              +(vdA3dP(:)/Area +AI(:,3)*ddRRddxi(cp,ij)/Area)*temp
     &              +DOTPROD(AI(:,3),dAIdxi(:,ij))/Area*(
     &                  CROSSPROD(AI(:,3),dUAdxi(:,1,iA))*dRRdxi(cp,2)
     &                + CROSSPROD(dUAdxi(:,2,iA),AI(:,3))*dRRdxi(cp,1)
     &                + vectsave(:) )
               
               dKAdP_M(:,iA) = dKAdP_M(:,iA) +vect(:)*stressBnd(ij)*coef
               
               
            Enddo
            gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &           - dKAdP_M(:,iA)*Area*detJac
            Enddo
            
            
            ! 5. derivatives of the bending stress 
            !    (with dble prod. by adjoint bending strain)
            dMdP_KA(:,:) = zero
            !    - derivatives of covariant metrics
            dAAIdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  dAAIdP(:,i,j) = dRRdxi(cp,i)*AI(:,j)
     &                 + AI(:,i)*dRRdxi(cp,j)
               Enddo
            Enddo
            !    - derivatives of contravariant metrics
            dAAEdP(:,:,:) = zero
            Do j = 1,2
               Do i = 1,2
                  Do l = 1,2
                     Do k = 1,2
                        dAAEdP(:,i,j) = dAAEdP(:,i,j) 
     &                       - AAE(i,k)*AAE(l,j)*dAAIdP(:,k,l)
                     Enddo
                  Enddo
               Enddo
            Enddo
            
            !    - first subterm (material dot derivative strain)
            dKdP_MA(:,:) = zero
            temp = SUM(dUdxi(:,1)*AIxAJ(:,1))+SUM(dUdxi(:,2)*AIxAJ(:,2))
            
            vect(:)   = CROSSPROD(dUdxi(:,1),AI(:,2)) 
     &                + CROSSPROD(AI(:,1),dUdxi(:,2))
            vdA3dP(:) = 
     &           -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &           - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &           - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
            vdA3dP(:) = vdA3dP(:)/Area
            vectsave(:) = vdA3dP(:)
            
            coef = one
            Do ij = 1,ntens
               
               if (ij==3) coef=two

               ! 1st term
               vect(:)   = ddUddxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA)
     &                 - vdA3dP(:)*stressBndAdj(ij,iA)*coef
               Enddo
               
               ! 2nd term
               vect(:) =  
     &              -dJdP(:)/Area/Area*(
     &               ScaTripleProd(dUdxi(:,1),dAIdxi(:,ij),AI(:,2))
     &              +ScaTripleProd(dUdxi(:,2),AI(:,1),dAIdxi(:,ij)))
     &              +one/Area*(
     &                CROSSPROD(dUdxi(:,1),dAIdxi(:,ij))*dRRdxi(cp,2)
     &               -CROSSPROD(dUdxi(:,1),AI(:,2))*ddRRddxi(cp,ij))
     &              +one/Area*(
     &                CROSSPROD(dAIdxi(:,ij),dUdxi(:,2))*dRRdxi(cp,1)
     &               -CROSSPROD(AI(:,1),dUdxi(:,2))*ddRRddxi(cp,ij))
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA) 
     &                 + vect(:)*stressBndAdj(ij,iA)*coef
               Enddo
               
               ! 3rd term
               vect(:)   = dAIdxi(:,ij)
               vdA3dP(:) = 
     &              -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &              - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &              - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
               vdA3dP(:) = vdA3dP(:)/Area
               
               vect(:) = 
     &              (-dJdP(:)/Area/Area*DOTPROD(AI(:,3),dAIdxi(:,ij))
     &              + vdA3dP(:)/Area +AI(:,3)*ddRRddxi(cp,ij)/Area)*temp
     &              +DOTPROD(AI(:,3),dAIdxi(:,ij))/Area*(
     &                  CROSSPROD(AI(:,3),dUdxi(:,1))*dRRdxi(cp,2)
     &                + CROSSPROD(dUdxi(:,2),AI(:,3))*dRRdxi(cp,1)
     &                + vectsave(:) )
               
               Do iA = 1,nadj
                  dKdP_MA(:,iA) = dKdP_MA(:,iA) 
     &                 + vect(:)*stressBndAdj(ij,iA)*coef
               Enddo
            Enddo
            dMdP_KA(:,:) = dKdP_MA(:,:)
            
            !    - second subterm (derivative material tensor)
            Do kl = 1,ntens
               k=voigt(kl,1); l=voigt(kl,2)
               
               Do ij = 1,ntens
                  i=voigt(ij,1); j=voigt(ij,2)
                  
                  dCdP(:) =
     &                   lambda*dAAEdP(:,i,j)*AAE(k,l)
     &                 + lambda*AAE(i,j)*dAAEdP(:,k,l)
     &                 + mu*dAAEdP(:,i,k)*AAE(j,l)
     &                 + mu*AAE(i,k)*dAAEdP(:,j,l)
     &                 + mu*dAAEdP(:,i,l)*AAE(j,k)
     &                 + mu*AAE(i,l)*dAAEdP(:,j,k)
                  
                  Do iA = 1,nadj
                     dMdP_KA(:,iA) = dMdP_KA(:,iA) 
     &                    + h*h*h/12.0d0*dCdP(:)
     &                      *strainBnd(ij)*strainBndAdj(kl,iA)
                  Enddo
               Enddo
            Enddo
            
            Do iA = 1,nadj
               gradWint_elem(iA,:,cp) = gradWint_elem(iA,:,cp) 
     &              - dMdP_KA(:,iA)*Area*detJac
            Enddo
            
         Enddo
         Endif ! test computeWint is True
         
         
         
c     Load vector
         If (computeWext) then
         gradWext_elem  = zero
         If (nb_load_srf>0) then
            UA(:,:) = zero
            Do iA = 1,nadj
               Do cp = 1,NNODE
                  UA(:,iA) = UA(:,iA) + R(cp)*UAelem(:,cp,iA)
               Enddo
            Enddo
         Endif
         
         Do i = 1,nb_load_srf
            numLoad = ind_load_loc(nb_load+1-i)
            call LectCle(JDLType(numLoad),KNumFace,KTypeDload)
            
            coef = one
            if (KTypeDload == 0) then
               if (KNumFace == 5) then
                  VectNorm(:) = AI(:,3)
               else
                  VectNorm(:) =-AI(:,3)
                  coef = -one
               endif
            elseif (KTypeDload == 1) then
               VectNorm(:) = (/ one, zero, zero/)
            elseif (KTypeDload == 2) then
               VectNorm(:) = (/ zero, one, zero/)
            elseif (KTypeDload == 3) then
               VectNorm(:) = (/ zero, zero, one/)
            elseif (KTypeDload == 4) then
               print*,'Warning: incomplete sensitivity in the load term'
               ! --> the derivative needs to be completed
               call norm(AI(1,:), MCRD, normV)
               VectNorm(:) = AI(:,1)/normV
            elseif (KTypeDload == 5) then
               print*,'Warning: incomplete sensitivity in the load term'
               ! --> the derivative needs to be completed
               call norm(AI(2,:), MCRD, normV)
               VectNorm(:) = AI(:,2)/normV
            elseif (KTypeDload == 6) then
               ! constant vertical load
               VectNorm(:) = (/ zero, zero, AI(3,3) /)
            endif
            
            UA_V(:) = zero
            Do iA = 1,nadj
               UA_V(iA) = DOTPROD(UA(:,iA),VectNorm(:))
            Enddo
            UA_V(:) = UA_V(:)*ADLMAG(numLoad)*detJac
            
            Do cp = 1,NNODEmap
               dJdP(:) = zero
               dJdP(:) = AIxAJ(:,1)*dRRdxi(cp,1)+AIxAJ(:,2)*dRRdxi(cp,2)

               Do iA = 1,nadj
                  UA_dVdP(:) = zero
                  if (KTypeDload == 0) then
                     vect(:)   = UA(:,iA)
                     vdA3dP(:) = 
     &                    -   DOTPROD(vect(:),AI(:,3))*dJdP(:)
     &                    - CROSSPROD(vect(:),AI(:,2))*dRRdxi(cp,1)
     &                    - CROSSPROD(AI(:,1),vect(:))*dRRdxi(cp,2)
                     vdA3dP(:) = vdA3dP(:)/Area
                     UA_dVdP(:) = vdA3dP(:)*ADLMAG(numLoad)*detJac*coef
                  elseif (KTypeDload == 6) then
                     vect(:) = (/ zero, zero, one/AI(3,3) /)
                     dJdP(:) = CROSSPROD(AI(:,2),vect(:))*dRRdxi(cp,1)
     &                       + CROSSPROD(vect(:),AI(:,1))*dRRdxi(cp,2)
                  Endif
                  
                  gradWext_elem(iA,:,cp) = gradWext_elem(iA,:,cp) 
     &                 + UA_V(iA)*dJdP(:)
     &                 + UA_dVdP(:)*Area
               Enddo
            Enddo
            
         Enddo
         Endif ! test computeWext is True
         
c     End load traitement
         
         
c     Update global gradient
         Do cp = 1,NNODEmap
            kk = sctr_map(cp)
            gradWext(:,:,kk) = gradWext(:,:,kk) + gradWext_elem(:,:,cp)
            gradWint(:,:,kk) = gradWint(:,:,kk) + gradWint_elem(:,:,cp)
         Enddo
         
         ENDIF ! test current mapping elem is active
         
      Enddo
c     End of the loop on integration points on main surf

c
c     ..................................................................
c
      
      End SUBROUTINE gradUELMAT30gloAdj
