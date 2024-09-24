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

      SUBROUTINE gradCPLGdisp(Xi,RI,dRIdxi,BI,dim,dimInterface,MCRD,
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
      Double precision :: R,dRdxi,ddRddxi,coef,U,dUdxi,dUdPI
      dimension R(NNODE),dRdxi(NNODE,3),ddRddxi(NNODE,6),U(3),
     &     dUdxi(3,dim),dUdPI(dim,NNODEinterface)
      
!     Curivilinear basis vectos
      Double precision :: AI,dAIdxi,dAIdP,dAIdPI
      dimension AI(3,dim),dAIdxi(3,6),
     &     dAIdP( 3,dim,MCRD,NNODE),
     &     dAIdPI(3,dim, dim,NNODEinterface)
      
!     Jacobian
      Double precision :: TI,T,normT,invT,vect1,vect2,dTIdP,dTIdPI,dTdP,
     &     dTdPI,dJTdP,dJTdPI
      dimension TI(3,3),T(3),vect1(3),vect2(3),
     &     dTIdP( 3,dimInterface,MCRD,NNODE),
     &     dTIdPI(3,dimInterface, dim,NNODEinterface),
     &     dTdP(3,MCRD,NNODE),dTdPI(3,dim,NNODEinterface),
     &     dJTdP(MCRD,NNODE),dJTdPI(dim,NNODEinterface)
      
!     Assembly
      Integer :: i,j,k,kk

C     ------------------------------------------------------------------
      
C     Computation :
c     -----------
      
      call evalnurbs_w2ndDerv(XI(:),R(:),dRdxi(:,:),ddRddxi(:,:))
      
!     Evaluate displacement and derv
      U(:) = zero
      Do i = 1,NNODE
         U(:) = U(:) + R(i)*SOL(:,i)
      Enddo
      
      dUdxi(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            dUdxi(:,j) = dUdxi(:,j) + dRdxi(i,j)*SOL(:,i)
         Enddo
      Enddo
      
c     - derivative versus internal CPs (L.dUdPi)
      dUdPI(:,:) = zero
      Do k = 1,dim
         call dot(lmbda(:),dUdxi(:,k),coef)
         dUdPI(k,:) = coef*RI(:)
      Enddo
      
      
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
      TI(:,:) = zero
      Do j = 1,dimInterface
         Do i = 1,dim
            TI(:,j) = TI(:,j) + BI(i,j)*AI(:,i)
         Enddo
      Enddo
      
c     - derivative versus physical CPs
      dTIdP(:,:,:,:)  = zero
      Do kk = 1,NNODE
         Do k = 1,MCRD
            Do j = 1,dimInterface
               Do i = 1,dim
                  dTIdP(:,j,k,kk) = dTIdP(:,j,k,kk) 
     &                            + BI(i,j)*dAIdP(:,i,k,kk)
               Enddo
            Enddo
         Enddo
      Enddo
      
c     - derivative versus interface CPs
      dTIdPI(:,:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do k = 1,dim
            Do j = 1,dimInterface
               dTIdPI(:,j,k,kk) = dRIdxi(kk,j)*AI(:,k)
               Do i = 1,dim
                  dTIdPI(:,j,k,kk) = dTIdPI(:,j,k,kk) 
     &                             + BI(i,j)*dAIdPI(:,i,k,kk)
               Enddo
            Enddo
         Enddo
      Enddo
      
      
c     Jacobian
      T(:) = zero
      dTdP(:,:,:)  = zero
      dTdPI(:,:,:) = zero
      If     (dimInterface == 1) then
         T(:) = TI(:,1)
         dTdP( :,:,:) = dTIdP( :,1,:,:)
         dTdPI(:,:,:) = dTIdPI(:,1,:,:)
      Elseif (dimInterface == 2) then
         call cross(TI(:,1),TI(:,2),T(:))
         
         Do kk = 1,NNODE
            Do k = 1,MCRD
               call cross(dTIdP(:,1,k,kk),TI(:,2),vect1)
               call cross(TI(:,1),dTIdP(:,2,k,kk),vect2)
               dTdP(:,k,kk) = vect1(:) + vect2(:)
            Enddo
         Enddo
         Do kk = 1,NNODEinterface
            Do k = 1,dim
               call cross(dTIdPI(:,1,k,kk),TI(:,2),vect1)
               call cross(TI(:,1),dTIdPI(:,2,k,kk),vect2)
               dTdPI(:,k,kk) = vect1(:) + vect2(:)
            Enddo
         Enddo
      Endif
      
c     - derivative versus physical CPs
      call norm(T(:),3,normT)
      invT = one/normT
      
      dJTdP(:,:) = zero
      Do k = 1,NNODE
         Do i = 1,MCRD
            dJTdP(i,k) = invT * SUM( T(:)*dTdP(:,i,k) )
         Enddo
      Enddo
            
c     - derivative versus interface CPs
      dJTdPI(:,:) = zero
      Do k = 1,NNODEinterface
         Do i = 1,dim
            dJTdPI(i,k) = invT * SUM( T(:)*dTdPI(:,i,k) )
         Enddo
      Enddo
      
c     Assembly
      dCdP(:,:)  = zero
      dCdPI(:,:) = zero
      
      call dot(lmbda(:),U(:),coef)
      dCdP(:MCRD,:) = coef*dJTdP(:,:)
      dCdPI(:dim,:) = dUdPI(:,:)*normT + coef*dJTdPI(:,:)
      !dCdPI(:dim,:) = dUdPI(:,:)*normT
      
c     
c     ..................................................................
c
      
      END SUBROUTINE gradCPLGdisp













































      SUBROUTINE gradCPLGdispEMBD(Xi,RI,dRIdxi,CI,dimmap,dim,
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
      Double precision :: R,dRdxi,ddRddxi,coef,U,dUdxi,dUdPI
      dimension R(NNODE),dRdxi(NNODE,3),ddRddxi(NNODE,6),U(3),
     &     dUdxi(3,dim),dUdPI(dim,NNODEinterface)
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
      Double precision :: AI,dAIdxi,dAIdPm,dAIdPe,dAIdPI
      dimension AI(3,dim),dAIdxi(3,dim,dim),
     &     dAIdPm(3,dim,MCRD,NNODEmap),dAIdPe(3,dim,dimmap,NNODE),
     &     dAIdPI(3,dim, dim,NNODEinterface)
      
!     Jacobian
      Double precision :: TI,T,normT,invT,vect1,vect2
      dimension TI(3,3),T(3),vect1(3),vect2(3)
      Double precision :: dTIdPm,dTIdPe,dTIdPI,dTdPm,dTdPe,dTdPI,dJTdPm,
     &     dJTdPe,dJTdPI
      dimension dTIdPm(3,dimInterface,MCRD,NNODEmap),
     &     dTIdPe(3,dimInterface,dimmap,NNODE),
     &     dTIdPI(3,dimInterface,   dim,NNODEinterface),
     &     dTdPm(3,MCRD,NNODEmap),dTdPe(3,dimmap,NNODE),
     &     dTdPI(3,dim,NNODEinterface),
     &     dJTdPm(MCRD,NNODEmap),dJTdPe(dimmap,NNODE),
     &     dJTdPI( dim,NNODEinterface)
      
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
      
      
      
!     Evaluate displacement and derv
      U(:) = zero
      Do i = 1,NNODE
         U(:) = U(:) + R(i)*SOL(:,i)
      Enddo
      
      dUdxi(:,:) = zero
      Do j = 1,dim
         Do i = 1,NNODE
            dUdxi(:,j) = dUdxi(:,j) + dRdxi(i,j)*SOL(:,i)
         Enddo
      Enddo
      
c     - derivative versus internal CPs (L.dUdPi)
      dUdPI(:,:) = zero
      Do k = 1,dim
         call dot(lmbda(:),dUdxi(:,k),coef)
         dUdPI(k,:) = coef*RI(:)
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
      TI(:,:) = zero
      Do j = 1,dimInterface
         Do i = 1,dim
            TI(:,j) = TI(:,j) + CI(i,j)*AI(:,i)
         Enddo
      Enddo
      
c     - derivative versus mapping CPs
      dTIdPm(:,:,:,:)  = zero
      Do kk = 1,NNODEmap
         Do i = 1,MCRD
            Do j = 1,dimInterface
               vectV(:) = zero
               Do l = 1,dim
                  vectV(:) = vectV(:) + CI(l,j)*dAIdPm(:,l,i,kk)
               Enddo
               dTIdPm(:,j,i,kk) = vectV(:)
            Enddo
         Enddo
      Enddo
      
c     - derivative versus embedded CPs
      dTIdPe(:,:,:,:)  = zero
      Do kk = 1,NNODE
         Do i = 1,dimmap
            Do j = 1,dimInterface
               vectV(:) = zero
               Do l = 1,dim
                  vectV(:) = vectV(:) + CI(l,j)*dAIdPe(:,l,i,kk)
               Enddo
               dTIdPe(:,j,i,kk) = vectV(:)
            Enddo
         Enddo
      Enddo
      
c     - derivative versus interface CPs
      dTIdPI(:,:,:,:) = zero
      Do kk = 1,NNODEinterface
         Do i = 1,dim
            Do j = 1,dimInterface
               vectV(:) = zero
               Do l = 1,dim
                  vectV(:) = vectV(:) + CI(l,j)*dAIdPI(:,l,i,kk)
               Enddo
               dTIdPI(:,j,i,kk) = dRIdxi(kk,j)*AI(:,i) + vectV(:)
            Enddo
         Enddo
      Enddo
      
      
c     Jacobian
      T(:) = zero
      dTdPm(:,:,:) = zero
      dTdPe(:,:,:) = zero
      dTdPI(:,:,:) = zero
      If     (dimInterface == 1) then
         T(:) = TI(:,1)
         dTdPm(:,:,:) = dTIdPm(:,1,:,:)
         dTdPe(:,:,:) = dTIdPe(:,1,:,:)
         dTdPI(:,:,:) = dTIdPI(:,1,:,:)
      Elseif (dimInterface == 2) then
         call cross(TI(:,1),TI(:,2),T(:))
         
         Do kk = 1,NNODEmap
            Do k = 1,MCRD
               call cross(dTIdPm(:,1,k,kk),TI(:,2),vect1)
               call cross(TI(:,1),dTIdPm(:,2,k,kk),vect2)
               dTdPm(:,k,kk) = vect1(:) + vect2(:)
            Enddo
         Enddo
         Do kk = 1,NNODE
            Do k = 1,dimmap
               call cross(dTIdPe(:,1,k,kk),TI(:,2),vect1)
               call cross(TI(:,1),dTIdPe(:,2,k,kk),vect2)
               dTdPe(:,k,kk) = vect1(:) + vect2(:)
            Enddo
         Enddo
         Do kk = 1,NNODEinterface
            Do k = 1,dim
               call cross(dTIdPI(:,1,k,kk),TI(:,2),vect1)
               call cross(TI(:,1),dTIdPI(:,2,k,kk),vect2)
               dTdPI(:,k,kk) = vect1(:) + vect2(:)
            Enddo
         Enddo
      Endif
      
c     - derivative versus physical CPs
      call norm(T(:),3,normT)
      invT = one/normT

      dJTdPm(:,:) = zero
      Do k = 1,NNODEmap
         Do i = 1,MCRD
            dJTdPm(i,k) = invT * SUM( T(:)*dTdPm(:,i,k) )
         Enddo
      Enddo
      
      dJTdPe(:,:) = zero
      Do k = 1,NNODE
         Do i = 1,dimmap
            dJTdPe(i,k) = invT * SUM( T(:)*dTdPe(:,i,k) )
         Enddo
      Enddo
            
c     - derivative versus interface CPs
      dJTdPI(:,:) = zero
      Do k = 1,NNODEinterface
         Do i = 1,dim
            dJTdPI(i,k) = invT * SUM( T(:)*dTdPI(:,i,k) )
         Enddo
      Enddo
      
c     Assembly
      dCdPm(:,:) = zero
      dCdPe(:,:) = zero
      dCdPI(:,:) = zero
      
      call dot(lmbda(:),U(:),coef)
      dCdPm(  :MCRD,:) = coef*dJTdPm(:,:)
      dCdPe(:dimmap,:) = coef*dJTdPe(:,:)
      dCdPI(   :dim,:) = dUdPI(:,:)*normT + coef*dJTdPI(:,:)
      !dCdPI(   :dim,:) = dUdPI(:,:)*normT !+ coef*dJTdPI(:,:)
      
c     
c     ..................................................................
c
      
      END SUBROUTINE gradCPLGdispEMBD
