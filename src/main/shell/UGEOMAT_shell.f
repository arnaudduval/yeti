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
c     Construction elementaire de la matrice de raideur geomtrique pour
c     la formulation coque type Kirchhoff-Love
c     --
      
      SUBROUTINE UGEOMAT3(NDOFEL,MCRD,NNODE,JELEM,NBINT,NumPatch,
     1     COORDS,TENSOR,MATERIAL_PROPERTIES,PROPS,JPROPS,Uelem,
     2     GEOMATRIXkl)
      
      
c     Input arguments :
c     ---------------
      Integer, intent(in) :: NDOFEL,MCRD,NNODE,JELEM,NumPatch,
     &     NBINT,JPROPS
      Character(len=*), intent(in) :: TENSOR
      Double precision, intent(in) :: COORDS, MATERIAL_PROPERTIES,PROPS,
     &     Uelem
      dimension COORDS(MCRD,NNODE), MATERIAL_PROPERTIES(2), PROPS(10),
     &     Uelem(MCRD,NNODE)
           
      
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: GEOMATRIXkl
      dimension GEOMATRIXkl(NDOFEL,NDOFEL)
      
      
      
c     Local variables :
c     ---------------
      
!     Parameters
      Double precision zero, one, two
      parameter (zero=0.d0, one=1.d0,two=2.d0)
      
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
      
!     For membrane forces and bending moments
      Double precision :: dvol,Area,temp1,temp2,temp12,temp21,res,
     &     memStrain,memForces,bndStrain,bndMoments,dRidxi,dRjdxi,Ui,
     &     ddRjddxi, ddRiddxi
      dimension memStrain(3),memForces(3),bndStrain(3),bndMoments(3),
     &     dRidxi(2),dRjdxi(2),Ui(MCRD),ddRjddxi(3), ddRiddxi(3)
      
      
!     --
!     For second derivatives of curvature
      double precision :: a1,a2,a3, a1xa2,a2xa3,a3xa1, area_inv, matRES
      dimension a1(3),a2(3),a3(3),a1xa2(3),a2xa3(3),a3xa1(3),matRES(3,3)
      
!     1. first term and second term
      double precision :: At_a2,At_a1T, Ab_a3_a2xa3,Ab_a3_a3xa1
      dimension At_a2(3,3),At_a1T(3,3),Ab_a3_a2xa3(3,3),Ab_a3_a3xa1(3,3)
      
      double precision :: matAA_t2_b3_23, matAAT_t2_b3_23
      double precision :: matAA_t1_b3_31, matAAT_t1_b3_31
      dimension matAA_t2_b3_23(3,3), matAAT_t2_b3_23(3,3)
      dimension matAA_t1_b3_31(3,3), matAAT_t1_b3_31(3,3)
      
!     2. third
!     - subterm 1
      double precision :: dRidRj_11,dRidRj_22,dRidRj_12,dRidRj_21,
     &     dRidRj_1221
      double precision :: akk, At_akk, akka3, At_a1xa2
      dimension akk(3,3), At_akk(3,3), At_a1xa2(3,3)
      
      double precision ::  matAA_tt22_3b_2323, matAA_tt11_3b_3131,
     &     matAA_tt21_3b_2331, matAA_tt12_3b_2331, Att_22, Att_11, 
     &     Att_21, Att_12, Ab_2323, Ab_3131, Ab_2331, Ab_3123
      dimension matAA_tt22_3b_2323(3,3),  matAA_tt11_3b_3131(3,3)
      dimension matAA_tt21_3b_2331(3,3),  matAA_tt12_3b_2331(3,3)
      dimension Att_22(3,3), Att_11(3,3), Att_21(3,3), Att_12(3,3),
     &      Ab_2323(3,3), Ab_3131(3,3), Ab_2331(3,3), Ab_3123(3,3)
!     - subterm 2
      double precision :: matAA_b23_2kk, matAA_b31_kk1
      double precision :: matAA_b23_kk1_b31_2kk, matAAT_b23_kk1_b31_2kk
      dimension matAA_b23_2kk(3,3), matAA_b31_kk1(3,3)
      dimension matAA_b23_kk1_b31_2kk(3,3), matAAT_b23_kk1_b31_2kk(3,3)
      
      double precision :: a2xakk, akkxa1, Ab_232kk, Ab_31kk1, Ab_232kkT,
     &     Ab_31kk1T, Ab_23kk1, Ab_312kkT
      dimension a2xakk(3), akkxa1(3), Ab_232kk(3,3), Ab_31kk1(3,3),
     &     Ab_232kkT(3,3), Ab_31kk1T(3,3), Ab_23kk1(3,3),
     &     Ab_312kkT(3,3)
      
      
!     Indices
      Integer :: n,nodi,nodj,dof_i,dof_j,k
      
      
      
C     ------------------------------------------------------------------

c     Initialization :
c     --------------
      NbPtInt = int(NBINT**(1.0/2.0)) ! number of Gauss points per dir.
      
c     Defining Gauss points coordinates and Gauss weights
      call Gauss(NbPtInt,2,GaussPdsCoord,0)
      
c     Stiffness matrix and force vector are initialized to zero
      GEOMATRIXkl(:,:) = zero
      
c     Material behaviour
      h = PROPS(2)
      E = MATERIAL_PROPERTIES(1)
      nu= MATERIAL_PROPERTIES(2)
      coef = E/(one-nu*nu)
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
         dRdxi  = zero
         ddRddxi= zero
         DetJac = zero
         call nurbsbasis(R,dRdxi,ddRddxi,DetJac,NNODE,
     &        GaussPdsCoord(2:,n),JELEM,NumPatch)
         
         
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
         
         
c     Computing membrane forces
         matH(:,:) = h*matH(:,:)
!         memStrain(:) = zero
!         Do nodi = 1,NNODE
!            dRidxi(:) = dRdxi(nodi,:)
!            Ui(:) = Uelem(:,nodi)
!            
!            call dot(AI(1,:), Ui(:), temp1)
!            call dot(AI(2,:), Ui(:), temp2)
!            memStrain(1) = memStrain(1) + dRidxi(1)*temp1
!            memStrain(2) = memStrain(2) + dRidxi(2)*temp2
!            memStrain(3) = memStrain(3) 
!     &           + dRidxi(2)*temp1 + dRidxi(1)*temp2
!         Enddo
         call uStrainMem_shell(Uelem,NNODE,MCRD,AI,dRdxi,memStrain)
         call MulVect(matH, memStrain, memForces, MCRD,MCRD)
         
c     Computing bending moments
         matH(:,:) = h*h/12.0D0 * matH(:,:)
         call uStrainBnd_shell(Uelem,NNODE,MCRD,AI,dAI1dxi,dAI2dxi,
     &        dRdxi,ddRddxi,bndStrain)
         call MulVect(matH, bndStrain, bndMoments, MCRD,MCRD)
         
c     
c     Assembling AMATRIX (membrane part)
         call SurfElem(AI(1,:), AI(2,:), Area)
         dvol = GaussPdsCoord(1,n)*DetJac
         memForces(:) = memForces(:)*Area*dvol
         Do nodj = 1,NNODE
            dRjdxi(:) = dRdxi(nodj,:)
            
            temp1 = memForces(1)*dRjdxi(1)
            temp2 = memForces(2)*dRjdxi(2)
            temp12= memForces(3)*dRjdxi(1)
            temp21= memForces(3)*dRjdxi(2)
            
            dof_j = (nodj-1)*MCRD
            Do nodi = 1,NNODE  !nodj,NNODE
               dRidxi(:) = dRdxi(nodi,:)
               
               res = temp1*dRidxi(1) + temp2*dRidxi(2)
     &              + temp12*dRidxi(2) + temp21*dRidxi(1)
               
               dof_i = (nodi-1)*MCRD
               Do k=1,MCRD
                  GEOMATRIXkl(dof_i+k,dof_j+k)
     &                 = GEOMATRIXkl(dof_i+k,dof_j+k) + res
               Enddo
            Enddo
         Enddo
         
c     
c     Assembling AMATRIX (curvature part)
         
         a1 = AI(1,:)
         a2 = AI(2,:)
         a3 = AI(3,:)
         call cross(a1, a2, a1xa2)
         call cross(a2, a3, a2xa3)
         call cross(a3, a1, a3xa1)
         
!     First and second term
         call matA_tilde (a2, At_a2)
         call matA_tildeT(a1, At_a1T)
         call matA_bar(a3,a2xa3, Ab_a3_a2xa3)
         call matA_bar(a3,a3xa1, Ab_a3_a3xa1)
         
         matAA_t2_b3_23(:,:) = At_a2( :,:) - Ab_a3_a2xa3(:,:)
         matAA_t1_b3_31(:,:) = At_a1T(:,:) - Ab_a3_a3xa1(:,:)
         call transposeMat(matAA_t2_b3_23, matAAT_t2_b3_23, 3)
         call transposeMat(matAA_t1_b3_31, matAAT_t1_b3_31, 3)
            
!     Third term
!     - subterm 1
         call matA_tilde(a1xa2, At_a1xa2)
         
         call matA_dbletidle( a2(:),a2(:), Att_22)
         call matA_dbletidle( a1(:),a1(:), Att_11)
         call matA_dbletidle( a2(:),a1(:), Att_21)
         call matA_dbletidle( a1(:),a2(:), Att_12)
         call matA_bar(a2xa3(:),a2xa3(:), Ab_2323)
         call matA_bar(a3xa1(:),a3xa1(:), Ab_3131)
         call matA_bar(a2xa3(:),a3xa1(:), Ab_2331)
         call matA_bar(a3xa1(:),a2xa3(:), Ab_3123)
         
         !call transposeMat(Ab_2331, Ab_2331T, 3)
         matAA_tt22_3b_2323(:,:) = Att_22(:,:) - 3.d0*Ab_2323(:,:)
         matAA_tt11_3b_3131(:,:) = Att_11(:,:) - 3.d0*Ab_3131(:,:)
         matAA_tt21_3b_2331(:,:) =-Att_21(:,:) - 3.d0*Ab_2331(:,:)
         matAA_tt12_3b_2331(:,:) =-Att_12(:,:) - 3.d0*Ab_3123(:,:)
         !call transposeMat(matAA_tt12_3b_2331, matAAT_tt12_3b_2331,3)
         
         
         akk(:,1) = dAI1dxi(1,:)
         akk(:,2) = dAI2dxi(2,:)
         akk(:,3) = dAI1dxi(2,:)
         
         area_inv = one/area
         bndMoments(3) = two*bndMoments(3)
         bndMoments(:) = bndMoments(:)*dvol
         Do k = 1,3
            
            call matA_tilde(akk(:,k), At_akk(:,:))
            call dot(akk(:,k),a3(:), akka3)
            
            
            call cross(a2(:),akk(:,k), a2xakk(:))
            call cross(akk(:,k),a1(:), akkxa1(:))
            call matA_bar( a2xa3(:),a2xakk(:), Ab_232kk)
            call matA_bar( a3xa1(:),akkxa1(:), Ab_31kk1)
            call matA_bar( a2xa3(:),akkxa1(:), Ab_23kk1)
            call matA_barT(a3xa1(:),a2xakk(:), Ab_312kkT)
            call transposeMat(Ab_232kk, Ab_232kkT, 3)
            call transposeMat(Ab_31kk1, Ab_31kk1T, 3)
                        
            matAA_b23_2kk(:,:) = Ab_232kk(:,:) + Ab_232kkT(:,:)
            matAA_b31_kk1(:,:) = Ab_31kk1(:,:) + Ab_31kk1T(:,:)
            matAA_b23_kk1_b31_2kk(:,:) = Ab_23kk1(:,:) + Ab_312kkT(:,:)
            call transposeMat(matAA_b23_kk1_b31_2kk,
     &           matAAT_b23_kk1_b31_2kk, 3)
            
            
            Do nodj = 1,NNODE
               dRjdxi(:)   = dRdxi(nodj,:)
               ddRjddxi(:) = ddRddxi(nodj,:)
               
               dof_j = (nodj-1)*MCRD
               Do nodi = 1,NNODE !nodj  !nodj,NNODE
                  dRidxi(:) = dRdxi(nodi,:)
                  ddRiddxi(:) = ddRddxi(nodi,:)
                  
                  dof_i = (nodi-1)*MCRD
                  
                  dRidRj_11 = dRidxi(1)*dRjdxi(1)
                  dRidRj_22 = dRidxi(2)*dRjdxi(2)
                  dRidRj_12 = dRidxi(1)*dRjdxi(2)
                  dRidRj_21 = dRidxi(2)*dRjdxi(1)
                  dRidRj_1221 = dRidRj_12 - dRidRj_21
                  
                  
!     --
!     1. First term:  \V{a}_{\alpha,\beta,r}\cdot\V{a}_{3,s}
                  matRES(:,:) = 
     &                   (ddRiddxi(k)*dRjdxi(1)) * matAA_t2_b3_23(:,:)
     &                 + (ddRiddxi(k)*dRjdxi(2)) * matAA_t1_b3_31(:,:)
!     --
!     2. Second term: \V{a}_{\alpha,\beta,s}\cdot\V{a}_{3,r}
                  matRES(:,:) = matRES(:,:)
     &                 + (dRidxi(1)*ddRjddxi(k)) * matAAT_t2_b3_23(:,:)
     &                 + (dRidxi(2)*ddRjddxi(k)) * matAAT_t1_b3_31(:,:)
!     
!     --
!     3. Third term: \V{a}_{\alpha,\beta}\cdot\V{a}_{3,rs}
!     -
!     Subterm 1
                  matRES(:,:) = matRES(:,:)
     &                 + dRidRj_1221 * At_akk(:,:)
     &                 -(dRidRj_1221 * At_a1xa2(:,:)
     &                 + dRidRj_11 * matAA_tt22_3b_2323(:,:)
     &                 + dRidRj_22 * matAA_tt11_3b_3131(:,:)
     &                 + dRidRj_12 * matAA_tt21_3b_2331(:,:)
     &                 + dRidRj_21 * matAA_tt12_3b_2331(:,:)
     &                 )*area_inv*akka3

!     -
!     Subterm 2
                  matRES(:,:) = matRES(:,:)
     &                 -(dRidRj_11 * matAA_b23_2kk(:,:)
     &                 + dRidRj_22 * matAA_b31_kk1(:,:)
     &                 + dRidRj_12 * matAA_b23_kk1_b31_2kk(:,:)
     &                 + dRidRj_21 *matAAT_b23_kk1_b31_2kk(:,:)
     &                 )*area_inv
                  
                  
                  Do k2 = 1,MCRD
                     Do k1 = 1,MCRD
                        GEOMATRIXkl(dof_i+k1,dof_j+k2)
     &                       = GEOMATRIXkl(dof_i+k1,dof_j+k2)
     &                       - matRES(k1,k2)*bndMoments(k)
                     Enddo
                  Enddo
                  
               Enddo
            Enddo
         Enddo
         
      Enddo
!     End loop on integration points
      
c     Symmetry: lower triangular part
!      Do dof_j = 1,NDOFEL-1 !2,NDOFEL
!         Do dof_i = dof_j+1,NDOFEL !1,dof_j-1
!            GEOMATRIXkl(dof_i,dof_j) = GEOMATRIXkl(dof_j,dof_i)
!         Enddo
!      Enddo
      
      
c     Fin calcul .......................................................
      
      End SUBROUTINE UGEOMAT3
