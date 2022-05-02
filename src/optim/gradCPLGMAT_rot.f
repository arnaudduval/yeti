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

      SUBROUTINE gradCPLGrot(Xi,RI,dRIdxi,BI,dim,dimInterface,MCRD,
     &     NNODE,NNODEinterface,COORDS,SOL,lmbda, dCdPI,dCdP)
      
      use parameters
      
      Implicit None
      
c     Input arguments :
c     ---------------
      
      ! patch
      Integer,          intent(in) :: dim,MCRD,NNODE
      Double precision, intent(in) :: COORDS,SOL,lmbda
      dimension COORDS(3,NNODE),SOL(3,NNODE),lmbda(3)
      
      ! interface
      Integer,         intent(in) :: dimInterface,NNODEinterface
      Double precision,intent(in) :: Xi,RI,dRIdxi,BI
      dimension Xi(3),RI(NNODEinterface),
     &     dRIdxi(NNODEinterface,dimInterface),BI(3,dimInterface)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: dCdPI,dCdP
      dimension dCdPI(3,NNODEinterface),dCdP(3,NNODE)
      
c     Local variables :
c     ---------------
      
!     Nurbs basis fcts
      Double precision :: R,dRdxi,ddRddxi,dUdxi,ddUddxi,ddUddiPI
      dimension R(NNODE),dRdxi(NNODE,3),ddRddxi(NNODE,6),
     &     dUdxi(3,dim),ddUddxi(3,6),ddUddiPI(3,dim,dim,NNODEinterface)
      
!     Rotation field
      Double precision :: PHI,phi1,phi2,vect,dA3dP,dA3dPI,A3_A1,A2_A3,
     &     vect1,vect2,normJ,invJ,dJdP,dJdPI,dPhi1dP,dPhi2dP,dPhi1dPI,
     &     dPhi2dPI,dPHIdP_T,dPHIdPI_T,coef,coef1,coef2,temp1,temp2
      dimension PHI(3),vect(3),dA3dP(3,MCRD,NNODE),
     &     dA3dPI(3,dim,NNODEinterface),A3_A1(3),A2_A3(3),
     &     vect1(3),vect2(3),dJdP(MCRD,NNODE),dJdPI(dim,NNODEinterface),
     &     dPhi1dP(MCRD,NNODE),dPhi2dP(MCRD,NNODE),
     &     dPhi1dPI(dim,NNODEinterface),dPhi2dPI(dim,NNODEinterface),
     &     dPHIdP_T(MCRD,NNODE),dPHIdPI_T(dim,NNODEinterface)
      
!     Curivilinear basis vectos
      Double precision :: AI,dAIdxi,dAIdP,dAIdPI
      dimension AI(3,3),dAIdxi(3,6),
     &     dAIdP( 3,dim,MCRD,NNODE),
     &     dAIdPI(3,dim, dim,NNODEinterface)
      
!     Jacobian
      Double precision :: T,dTdP,dTdPI
      dimension T(3),dTdP(3,MCRD,NNODE),dTdPI(3,dim,NNODEinterface)
      
!     Assembly
      Integer :: i,j,k,kk

C     ------------------------------------------------------------------
      
C     Computation :
c     -----------
      
      call evalnurbs_w2ndDerv(XI(:),R(:),dRdxi(:,:),ddRddxi(:,:))
      
c     Covariant vectors
      AI(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            AI(:,j) = AI(:,j) + dRdxi(i,j)*COORDS(:,i)
         Enddo
      Enddo
      
      dAIdxi(:,:) = zero
      Do j = 1,6
         Do i = 1,NNODE
            dAIdxi(:,j) = dAIdxi(:,j) + ddRddxi(i,j)*COORDS(:,i)
         Enddo
      Enddo
            
!     - derivative versus physical CPs
      dAIdP(:,:,:,:) = zero
      Do k = 1,NNODE
         Do i = 1,MCRD
            Do j = 1,dim
               dAIdP(i,j,i,k) = dRdxi(k,j)
            Enddo
         Enddo
      Enddo
      
!     - derivative versus interface CPs
      dAIdPI(:,:,:,:) = zero
      Do k = 1,NNODEinterface
         Do i = 1,dim
            dAIdPI(:,i, i,k) = RI(k)*dAIdxi(:,i)
            
            Do j =1,dim
            If (j /= i) then
               kk = i+j+1
               dAIdPI(:,j, i,k) = RI(k)*dAIdxi(:,kk)
               dAIdPI(:,i, j,k) = RI(k)*dAIdxi(:,kk)
            Endif
            Enddo
            
         Enddo
      Enddo
      
      
c     Tangent vectors
      T(:) = zero
      Do i = 1,dim
         T(:) = T(:) + BI(i,1)*AI(:,i)
      Enddo
      
!     - derivative versus physical CPs
      dTdP(:,:,:)  = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            Do i = 1,dim
               dTdP(:,k,kk) = dTdP(:,k,kk) + BI(i,1)*dAIdP(:,i,k,kk)
            Enddo
         Enddo
      Enddo
      
!     - derivative versus interface CPs
      dTdPI(:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            dTdPI(:,k,kk) = dRIdxi(kk,1)*AI(:,k)
            Do i = 1,dim
               dTdPI(:,k,kk) = dTdPI(:,k,kk) + BI(i,1)*dAIdPI(:,i,k,kk)
            Enddo
         Enddo
      Enddo
      
      
c     Rotation field
      call cross(AI(:,1),AI(:,2),vect(:))
      normJ= sqrt( vect(1)**two + vect(2)**two + vect(3)**two )
      invJ = one/normJ
      AI(:,3) = invJ*vect(:)
      
      vect(:) = zero
      Do kk = 1,NNODE
         vect(:) = vect(:) + dRdxi(kk,2)*SOL(:,kk)
      Enddo
      call dot(vect(:),AI(:,3),phi1)
      phi1 = invJ*phi1
      
      vect(:) = zero
      Do kk = 1,NNODE
         vect(:) = vect(:) + dRdxi(kk,1)*SOL(:,kk)
      Enddo
      call dot(vect(:),AI(:,3),phi2)
      phi2 =-invJ*phi2
      
      PHI(:) = phi1*AI(:,1) + phi2*AI(:,2)
      
!     - derivative of J versus CPs
      call cross(AI(:,3),AI(:,1), A3_A1)
      call cross(AI(:,2),AI(:,3), A2_A3)
      
      dJdP(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            call dot(A2_A3(:),dAIdP(:,1,k,kk),temp1)
            call dot(A3_A1(:),dAIdP(:,2,k,kk),temp2)
            dJdP(k,kk) = temp1 + temp2
         Enddo
      Enddo
      
      dJdPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(A2_A3(:),dAIdPI(:,1,k,kk),temp1)
            call dot(A3_A1(:),dAIdPI(:,2,k,kk),temp2)
            dJdPI(k,kk) = temp1 + temp2
         Enddo
      Enddo
      
!     - derivative of A3 versus CPs
      dA3dP(:,:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            call cross(dAIdP(:,1,k,kk),AI(:,2),vect1(:))
            call cross(AI(:,1),dAIdP(:,2,k,kk),vect2(:))
            dA3dP(:,k,kk) = -invJ*AI(:,3)*dJdP(k,kk)
     &           + invJ*vect1(:) + invJ*vect2(:)
         Enddo
      Enddo
      
      dA3dPI(:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call cross(dAIdPI(:,1,k,kk),AI(:,2),vect1(:))
            call cross(AI(:,1),dAIdPI(:,2,k,kk),vect2(:))
            dA3dPI(:,k,kk) = -invJ*AI(:,3)*dJdPI(k,kk)
     &           + invJ*vect1(:) + invJ*vect2(:)
         Enddo
      Enddo
      
!     - derivative of Disp versus interface CPs
      dUdxi(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            dUdxi(:,j) = dUdxi(:,j) + dRdxi(i,j)*SOL(:,i)
         Enddo
      Enddo
      
      ddUddxi(:,:) = zero
      Do i = 1,NNODE
         ddUddxi(:,1) = ddUddxi(:,1) + ddRddxi(i,1)*SOL(:,i)
         ddUddxi(:,2) = ddUddxi(:,2) + ddRddxi(i,2)*SOL(:,i)
         ddUddxi(:,4) = ddUddxi(:,4) + ddRddxi(i,4)*SOL(:,i)
      Enddo
      
      ddUddiPI(:,:,:,:) = zero
      Do k = 1,NNODEinterface
         Do i = 1,dim
            ddUddiPI(:,i, i,k) = RI(k)*ddUddxi(:,i)
            Do j =1,dim
            If (j /= i) then
               kk = i+j+1
               ddUddiPI(:,j, i,k) = RI(k)*ddUddxi(:,kk)
               ddUddiPI(:,i, j,k) = RI(k)*ddUddxi(:,kk)
            Endif
            Enddo
         Enddo
      Enddo
      
!     - derivative of angle versus CPs
      dPhi1dP(:,:) = zero
      dPhi2dP(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            call dot(dUdxi(:,2),dA3dP(:,k,kk),coef2)
            dPhi1dP(k,kk) = invJ*( coef2 - phi1*dJdP(k,kk) )
            
            call dot(dUdxi(:,1),dA3dP(:,k,kk),coef1)
            dPhi2dP(k,kk) = invJ*(-coef1 - phi2*dJdP(k,kk) )
         Enddo
      Enddo

      dPhi1dPI(:,:) = zero
      dPhi2dPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(        dUdxi(:,2), dA3dPI(:,k,kk), coef2)
            call dot(ddUddiPI(:,2,k,kk),        AI(:,3), temp2)
            dPhi1dPI(k,kk) = invJ*( temp2+coef2 - phi1*dJdPI(k,kk) )
            
            call dot(       dUdxi(:,1),  dA3dPI(:,k,kk), coef1)
            call dot(ddUddiPI(:,1,k,kk),        AI(:,3), temp1)
            dPhi2dPI(k,kk) = invJ*(-temp1-coef1 - phi2*dJdPI(k,kk) )
         Enddo
      Enddo
      
      dPHIdP_T(:,:) = zero
      call dot(AI(:,1),T(:),coef1)
      call dot(AI(:,2),T(:),coef2)
      Do kk = 1,NNODE
         Do k = 1,MCRD
            call dot(dAIdP(:,1,k,kk),T(:),temp1)
            call dot(dAIdP(:,2,k,kk),T(:),temp2)
            dPHIdP_T(k,kk) = coef1*dPhi1dP(k,kk) + coef2*dPhi2dP(k,kk)
     &           + temp1*Phi1 + temp2*Phi2
         Enddo
      Enddo
      
      dPHIdPI_T(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(dAIdPI(:,1,k,kk),T(:),temp1)
            call dot(dAIdPI(:,2,k,kk),T(:),temp2)
            dPHIdPI_T(k,kk) = coef1*dPhi1dPI(k,kk) +coef2*dPhi2dPI(k,kk)
     &           + temp1*Phi1 + temp2*Phi2
         Enddo
      Enddo
      
      
c     Assembly
      dCdP(:,:)  = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            call dot(PHI(:),dTdP(:,k,kk),coef)
            dCdP(k,kk) = lmbda(1)*dPHIdP_T(k,kk) + lmbda(1)*coef
         Enddo
      Enddo
      
      dCdPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(PHI(:),dTdPI(:,k,kk),coef)
            dCdPI(k,kk) = lmbda(1)*dPHIdPI_T(k,kk) + lmbda(1)*coef
         Enddo
      Enddo
      
            
c     
c     ..................................................................
c
      
      END SUBROUTINE gradCPLGrot






































      SUBROUTINE gradCPLGrotEMBD(Xi,RI,dRIdxi,CI,dimmap,dim,
     &     dimInterface,MCRD,NNODEmap,NNODE,NNODEinterface,nb_cp,COORDS,
     &     COORDSall,SOL,lmbda, dCdPI,dCdPe,dCdPm)
      
      use parameters
      use embeddedMapping
      
      Implicit None
      
c     Input arguments :
c     ---------------
      
      ! mapping
      Integer,          intent(in) :: dimmap,NNODEmap,nb_cp
      Double precision, intent(in) :: COORDSall
      dimension COORDSall(3,nb_cp)
      
      ! patch
      Integer,          intent(in) :: dim,MCRD,NNODE
      Double precision, intent(in) :: COORDS,SOL,lmbda
      dimension COORDS(3,NNODE),SOL(3,NNODE),lmbda(3)
      
      ! interface
      Integer,         intent(in) :: dimInterface,NNODEinterface
      Double precision,intent(in) :: Xi,RI,dRIdxi,CI
      dimension Xi(3),RI(NNODEinterface),
     &     dRIdxi(NNODEinterface,dimInterface),CI(3,dimInterface)
      
c     Output variables :
c     ----------------
      Double precision, intent(out) :: dCdPI,dCdPe,dCdPm
      dimension dCdPI(3,NNODEinterface),dCdPe(3,NNODE),dCdPm(3,NNODEmap)
      
c     Local variables :
c     ---------------
      
!     Embedded entities
!     - Nurbs basis fcts
      Double precision :: R,dRdxi,ddRddxi,dUdxi,ddUddxi,ddUddiPI
      dimension R(NNODE),dRdxi(NNODE,3),ddRddxi(NNODE,6),
     &     dUdxi(3,dim),ddUddxi(3,6),ddUddiPI(3,dim,dim,NNODEinterface)
!     - vectors
      Double precision :: SI,dSIdxi,dSIdPI,dSIdPe
      dimension SI(3,dim),dSIdxi(3,dim,dim),
     &     dSIdPI(3,dim,dim,NNODEinterface),dSIdPe(3,dim,dimmap,NNODE)
      
!     Mapping
!     - Nurbs basis fcts
      Double precision :: vectV, Rm, dRmdxi, ddRmddxi, COORDSmap
      dimension vectV(3), Rm(NNODEmap), dRmdxi(NNODEmap,3),
     &     ddRmddxi(NNODEmap,6),COORDSmap(3,NNODEmap)
!     - vectors
      Double precision :: VI,dVIdxi,dVIdPm,dVIdPe,dVIdPI
      dimension VI(3,dimmap),dVIdxi(3,dimmap,dimmap),
     &     dVIdPm(3,dimmap,MCRD,NNODEmap),dVIdPe(3,dimmap,dimmap,NNODE),
     &     dVIdPI(3,dimmap,dim,NNODEinterface)
      
!     Curivilinear basis vectors
      Double precision :: AI,dAIdxi,dAIdPm,dAIdPe,dAIdPI,dA3dPm,dA3dPe,
     &     dA3dPI,A3_A1,A2_A3
      dimension AI(3,3),dAIdxi(3,dim,dim),
     &     dAIdPm(3,dim,MCRD,NNODEmap),dAIdPe(3,dim,dimmap,NNODE),
     &     dAIdPI(3,dim, dim,NNODEinterface),
     &     dA3dPm(3,MCRD,NNODEmap),dA3dPe(3,dimmap,NNODE),
     &     dA3dPI(3,dim,NNODEinterface),A3_A1(3),A2_A3(3)
      

!     Rotation field
      Double precision :: PHI,phi1,phi2,vect,vect1,vect2,normJ,invJ,
     &     dJdPm,dJdPe,dJdPI,dPhi1dPm,dPhi2dPm,dPhi1dPe,dPhi2dPe,
     &     dPhi1dPI,dPhi2dPI,dPHIdPm_T,dPHIdPe_T,dPHIdPI_T,coef,coef1,
     &     coef2,temp1,temp2
      dimension PHI(3),vect(3),vect1(3),vect2(3),dJdPm(MCRD,NNODEmap),
     &     dJdPe(dimmap,NNODE),dJdPI(dim,NNODEinterface),
     &     dPhi1dPm(  MCRD,NNODEmap),dPhi2dPm(  MCRD,NNODEmap),
     &     dPhi1dPe(dimmap,NNODE),   dPhi2dPe(dimmap,NNODE),
     &     dPhi1dPI(dim,NNODEinterface),dPhi2dPI(dim,NNODEinterface),
     &     dPHIdPm_T(MCRD,NNODEmap),dPHIdPe_T(dimmap,NNODE),
     &     dPHIdPI_T(dim,NNODEinterface)
      
!     Jacobian
      Double precision :: T,dTdPm,dTdPe,dTdPI
      dimension T(3),dTdPm(3,MCRD,NNODEmap),dTdPe(3,dimmap,NNODE),
     &     dTdPI(3,dim,NNODEinterface)
      
!     Assembly
      Integer :: i,j,k,l,kk

C     ------------------------------------------------------------------
      
C     Computation :
c     -----------
      
c     Embedded Entities
      call evalnurbs_w2ndDerv(XI(:),R(:),dRdxi(:,:),ddRddxi(:,:))
      vectV(:) = zero
      Do i = 1,NNODE
         vectV(:MCRD) = vectV(:MCRD) + R(i)*COORDS(:,i)
      Enddo
      
      SI(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            SI(:,j) = SI(:,j) + dRdxi(i,j)*COORDS(:,i)
         Enddo
      Enddo
      
      dSIdxi(:,:,:) = zero
      Do k = 1,dim
         Do i = 1,NNODE
            dSIdxi(:,k,k) = dSIdxi(:,k,k) + ddRddxi(i,k)*COORDS(:,i)
         Enddo
         Do j = k+1,dim
            kk = k+j+1
            Do i = 1,NNODE
               dSIdxi(:,j,k) = dSIdxi(:,j,k) + ddRddxi(i,kk)*COORDS(:,i)
            Enddo
            dSIdxi(:,k,j) = dSIdxi(:,j,k)
         Enddo
      Enddo

      dSIdPe(:,:,:,:) = zero
      Do kk = 1,NNODE
         Do i = 1,dimmap
            Do j = 1,dim
               dSIdPe(i,j,i,kk) = dRdxi(kk,j)
            Enddo
         Enddo
      Enddo
      
      dSIdPI(:,:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do i = 1,dim
            Do j = 1,dim
               dSIdPI(:,j,i,kk) = RI(kk)*dSIdxi(:,j,i)
            Enddo
         Enddo
      Enddo
      
      
c     Mapping
      call updateMapElementNumber(VectV(:))
      call evalnurbs_mapping_w2ndDerv(VectV(:),Rm(:),dRmdxi(:,:),
     &     ddRmddxi(:,:))
      
      Do i = 1,NNODEmap
         COORDSmap(:,i) = COORDSall(:,IEN_map(i,current_map_elem))
      Enddo
      
      VI(:,:) = zero
      Do j = 1,dimmap
         Do i = 1,NNODEmap
            VI(:,j) = VI(:,j) + dRmdxi(i,j)*COORDSmap(:,i)
         Enddo
      Enddo
      
      dVIdxi(:,:,:) = zero
      Do k = 1,dimmap
         Do i = 1,NNODEmap
            dVIdxi(:,k,k) = dVIdxi(:,k,k) + ddRmddxi(i,k)*COORDSmap(:,i)
         Enddo
         Do j = k+1,dimmap
            kk = k+j+1
            Do i = 1,NNODEmap
               dVIdxi(:,j,k)=dVIdxi(:,j,k)+ddRmddxi(i,kk)*COORDSmap(:,i)
            Enddo
            dVIdxi(:,k,j) = dVIdxi(:,j,k)
         Enddo
      Enddo
      
      dVIdPm(:,:,:,:) = zero
      Do kk = 1,NNODEmap
         Do i = 1,MCRD
            Do j = 1,dimmap
               dVIdPm(i,j,i,kk) = dRmdxi(kk,j)
            Enddo
         Enddo
      Enddo
      
      dVIdPe(:,:,:,:) = zero
      Do kk = 1,NNODE
         Do i = 1,dimmap
            Do j = 1,dimmap
               dVIdPe(:,j,i,kk) = R(kk)*dVIdxi(:,j,i)
            Enddo
         Enddo
      Enddo
      
      dVIdPI(:,:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do i = 1,dim
            Do j =  1,dimmap
               vectV(:) = zero
               Do k = 1,dimmap
                  vectV(:) = vectV(:) + SI(l,i)*dVIdxi(:,j,l)
               Enddo
               dVIdPI(:,j,i,kk) = RI(kk)*vectV(:)
            Enddo
         Enddo
      Enddo
      
      
      
c     Covariant vectors
      AI(:,:) = zero
      Do j = 1,dim
         Do i = 1,dimmap
            AI(:,j) = AI(:,j) + SI(i,j)*VI(:,i)
         Enddo
      Enddo
      
      dAIdxi(:,:,:) = zero
      Do j = 1,dim
         Do i = j,dim
            Do k = 1,dimmap
               dAIdxi(:,i,j) = dAIdxi(:,i,j) + dSIdxi(k,i,j)*VI(:,k)
               
               vectV(:) = zero
               Do l = 1,dimmap
                  vectV(:) = vectV(:) + SI(l,j)*dVIdxi(:,k,l)
               Enddo
               dAIdxi(:,i,j) = dAIdxi(:,i,j) + SI(k,i)*vectV(:)
            Enddo
            dAIdxi(:,j,i) = dAIdxi(:,i,j)
         Enddo
      Enddo
      
!     - derivative versus mapping CPs
      dAIdPm(:,:,:,:) = zero
      Do kk = 1,NNODEmap
         Do i = 1,MCRD
            Do j = 1,dim
               vectV(:) = zero
               Do l = 1,dimmap
                  vectV(:) = vectV(:) + SI(l,j)*dVIdPm(:,l,i,kk)
               Enddo
               dAIdPm(:,j,i,k) = vectV(:)
            Enddo
         Enddo
      Enddo
      
!     - derivative versus embedded CPs
      dAIdPe(:,:,:,:) = zero
      Do kk = 1,NNODE
         Do i = 1,dimmap
            Do j = 1,dim
               vectV(:) = zero
               Do l = 1,dimmap
                  vectV(:) = vectV(:) + dSIdPe(l,j,i,kk)*VI(:,l)
                  vectV(:) = vectV(:) + SI(l,j)*dVIdPe(:,l,i,kk)
               Enddo
               dAIdPe(:,j,i,kk) = vectV(:)
            Enddo
         Enddo
      Enddo
      
!     - derivative versus interface CPs
      dAIdPI(:,:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do i = 1,dim
            Do j = 1,dim
               vectV(:) = zero
               Do l = 1,dimmap
                  vectV(:) = vectV(:) + dSIdPI(l,j,i,kk)*VI(:,l)
                  vectV(:) = vectV(:) + SI(l,j)*dVIdPI(:,l,i,kk)
               Enddo
               dAIdPI(:,j,i,kk) = vectV(:)
            Enddo
         Enddo
      Enddo
      
      
      
c     Tangent vectors
      T(:) = zero
      Do i = 1,dim
         T(:) = T(:) + CI(i,1)*AI(:,i)
      Enddo
      
      
c     - derivative versus mapping CPs
      dTdPm(:,:,:)  = zero
      Do kk = 1,NNODEmap
         Do i = 1,MCRD
            vectV(:) = zero
            Do l = 1,dim
               vectV(:) = vectV(:) + CI(l,1)*dAIdPm(:,l,i,kk)
            Enddo
            dTdPm(:,i,kk) = vectV(:)
         Enddo
      Enddo
      
c     - derivative versus embedded CPs
      dTdPe(:,:,:)  = zero
      Do kk = 1,NNODE
         Do i = 1,dimmap
            vectV(:) = zero
            Do l = 1,dim
               vectV(:) = vectV(:) + CI(l,1)*dAIdPe(:,l,i,kk)
            Enddo
            dTdPe(:,i,kk) = vectV(:)
         Enddo
      Enddo
      
c     - derivative versus interface CPs
      dTdPI(:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do i = 1,dim
            vectV(:) = zero
            Do l = 1,dim
               vectV(:) = vectV(:) + CI(l,1)*dAIdPI(:,l,i,kk)
            Enddo
            dTdPI(:,i,kk) = dRIdxi(kk,1)*AI(:,i) + vectV(:)
         Enddo
      Enddo
      
      
      
c     Rotation field
      call cross(AI(:,1),AI(:,2),vect(:))
      normJ= sqrt( vect(1)**two + vect(2)**two + vect(3)**two )
      invJ = one/normJ
      AI(:,3) = invJ*vect(:)
      
      vect(:) = zero
      Do kk = 1,NNODE
         vect(:) = vect(:) + dRdxi(kk,2)*SOL(:,kk)
      Enddo
      call dot(vect(:),AI(:,3),phi1)
      phi1 = invJ*phi1
      
      vect(:) = zero
      Do kk = 1,NNODE
         vect(:) = vect(:) + dRdxi(kk,1)*SOL(:,kk)
      Enddo
      call dot(vect(:),AI(:,3),phi2)
      phi2 =-invJ*phi2
      
      PHI(:) = phi1*AI(:,1) + phi2*AI(:,2)
      
!     - derivative of J versus CPs
      call cross(AI(:,3),AI(:,1), A3_A1)
      call cross(AI(:,2),AI(:,3), A2_A3)
      
      dJdPm(:,:) = zero
      Do kk = 1,NNODEmap
         Do k = 1,MCRD
            call dot(A2_A3(:),dAIdPm(:,1,k,kk),temp1)
            call dot(A3_A1(:),dAIdPm(:,2,k,kk),temp2)
            dJdPm(k,kk) = temp1 + temp2
         Enddo
      Enddo
      
      dJdPe(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,dimmap
            call dot(A2_A3(:),dAIdPe(:,1,k,kk),temp1)
            call dot(A3_A1(:),dAIdPe(:,2,k,kk),temp2)
            dJdPe(k,kk) = temp1 + temp2
         Enddo
      Enddo

      dJdPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(A2_A3(:),dAIdPI(:,1,k,kk),temp1)
            call dot(A3_A1(:),dAIdPI(:,2,k,kk),temp2)
            dJdPI(k,kk) = temp1 + temp2
         Enddo
      Enddo
      
!     - derivative of A3 versus CPs
      dA3dPm(:,:,:) = zero
      Do kk = 1,NNODEmap
         Do k = 1,MCRD
            call cross(dAIdPm(:,1,k,kk),AI(:,2),vect1(:))
            call cross(AI(:,1),dAIdPm(:,2,k,kk),vect2(:))
            dA3dPm(:,k,kk) = -invJ*AI(:,3)*dJdPm(k,kk)
     &           + invJ*vect1(:) + invJ*vect2(:)
         Enddo
      Enddo
      
      dA3dPe(:,:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,dimmap
            call cross(dAIdPe(:,1,k,kk),AI(:,2),vect1(:))
            call cross(AI(:,1),dAIdPe(:,2,k,kk),vect2(:))
            dA3dPe(:,k,kk) = -invJ*AI(:,3)*dJdPe(k,kk)
     &           + invJ*vect1(:) + invJ*vect2(:)
         Enddo
      Enddo
      
      dA3dPI(:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call cross(dAIdPI(:,1,k,kk),AI(:,2),vect1(:))
            call cross(AI(:,1),dAIdPI(:,2,k,kk),vect2(:))
            dA3dPI(:,k,kk) = -invJ*AI(:,3)*dJdPI(k,kk)
     &           + invJ*vect1(:) + invJ*vect2(:)
         Enddo
      Enddo
      
!     - derivative of Disp versus interface CPs
      dUdxi(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            dUdxi(:,j) = dUdxi(:,j) + dRdxi(i,j)*SOL(:,i)
         Enddo
      Enddo
      
      ddUddxi(:,:) = zero
      Do i = 1,NNODE
         ddUddxi(:,1) = ddUddxi(:,1) + ddRddxi(i,1)*SOL(:,i)
         ddUddxi(:,2) = ddUddxi(:,2) + ddRddxi(i,2)*SOL(:,i)
         ddUddxi(:,4) = ddUddxi(:,4) + ddRddxi(i,4)*SOL(:,i)
      Enddo
      
      ddUddiPI(:,:,:,:) = zero
      Do k = 1,NNODEinterface
         Do i = 1,dim
            ddUddiPI(:,i, i,k) = RI(k)*ddUddxi(:,i)
            Do j =1,dim
            If (j /= i) then
               kk = i+j+1
               ddUddiPI(:,j, i,k) = RI(k)*ddUddxi(:,kk)
               ddUddiPI(:,i, j,k) = RI(k)*ddUddxi(:,kk)
            Endif
            Enddo
         Enddo
      Enddo
      
!     - derivative of angle versus CPs
      dPhi1dPm(:,:) = zero
      dPhi2dPm(:,:) = zero
      Do kk = 1,NNODEmap
         Do k = 1,MCRD
            call dot(dUdxi(:,2),dA3dPm(:,k,kk),coef2)
            dPhi1dPm(k,kk) = invJ*( coef2 - phi1*dJdPm(k,kk) )
            
            call dot(dUdxi(:,1),dA3dPm(:,k,kk),coef1)
            dPhi2dPm(k,kk) = invJ*(-coef1 - phi2*dJdPm(k,kk) )
         Enddo
      Enddo

      dPhi1dPe(:,:) = zero
      dPhi2dPe(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,dimmap
            call dot(dUdxi(:,2),dA3dPe(:,k,kk),coef2)
            dPhi1dPe(k,kk) = invJ*( coef2 - phi1*dJdPe(k,kk) )
            
            call dot(dUdxi(:,1),dA3dPe(:,k,kk),coef1)
            dPhi2dPe(k,kk) = invJ*(-coef1 - phi2*dJdPe(k,kk) )
         Enddo
      Enddo

      dPhi1dPI(:,:) = zero
      dPhi2dPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(        dUdxi(:,2), dA3dPI(:,k,kk), coef2)
            call dot(ddUddiPI(:,2,k,kk),        AI(:,3), temp2)
            dPhi1dPI(k,kk) = invJ*( temp2+coef2 - phi1*dJdPI(k,kk) )
            
            call dot(       dUdxi(:,1),  dA3dPI(:,k,kk), coef1)
            call dot(ddUddiPI(:,1,k,kk),        AI(:,3), temp1)
            dPhi2dPI(k,kk) = invJ*(-temp1-coef1 - phi2*dJdPI(k,kk) )
         Enddo
      Enddo
      
      dPHIdPm_T(:,:) = zero
      call dot(AI(:,1),T(:),coef1)
      call dot(AI(:,2),T(:),coef2)
      Do kk = 1,NNODEmap
         Do k = 1,MCRD
            call dot(dAIdPm(:,1,k,kk),T(:),temp1)
            call dot(dAIdPm(:,2,k,kk),T(:),temp2)
            dPHIdPm_T(k,kk) = coef1*dPhi1dPm(k,kk) +coef2*dPhi2dPm(k,kk)
     &           + temp1*Phi1 + temp2*Phi2
         Enddo
      Enddo
      
      dPHIdPe_T(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,dimmap
            call dot(dAIdPe(:,1,k,kk),T(:),temp1)
            call dot(dAIdPe(:,2,k,kk),T(:),temp2)
            dPHIdPe_T(k,kk) = coef1*dPhi1dPe(k,kk) +coef2*dPhi2dPe(k,kk)
     &           + temp1*Phi1 + temp2*Phi2
         Enddo
      Enddo      
      
      dPHIdPI_T(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(dAIdPI(:,1,k,kk),T(:),temp1)
            call dot(dAIdPI(:,2,k,kk),T(:),temp2)
            dPHIdPI_T(k,kk) = coef1*dPhi1dPI(k,kk) +coef2*dPhi2dPI(k,kk)
     &           + temp1*Phi1 + temp2*Phi2
         Enddo
      Enddo
      
      
c     Assembly
      dCdPm(:,:) = zero
      Do kk = 1,NNODEmap
         Do k = 1,MCRD
            call dot(PHI(:),dTdPm(:,k,kk),coef)
            dCdPm(k,kk) = lmbda(1)*dPHIdPm_T(k,kk) + lmbda(1)*coef
         Enddo
      Enddo
      
      dCdPe(:,:) = zero
      Do kk = 1,NNODE
         Do k = 1,dimmap
            call dot(PHI(:),dTdPe(:,k,kk),coef)
            dCdPe(k,kk) = lmbda(1)*dPHIdPe_T(k,kk) + lmbda(1)*coef
         Enddo
      Enddo

      dCdPI(:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            call dot(PHI(:),dTdPI(:,k,kk),coef)
            dCdPI(k,kk) = lmbda(1)*dPHIdPI_T(k,kk) + lmbda(1)*coef
         Enddo
      Enddo
      
            
c     
c     ..................................................................
c
      
      END SUBROUTINE gradCPLGrotEMBD

