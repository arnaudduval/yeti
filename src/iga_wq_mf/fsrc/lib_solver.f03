! module CGsolver

!     type CGtype
!         procedure(), pointer, nopass :: matrixfree => NULL()
!         procedure(), pointer, nopass :: fastdiag => NULL()
!     end type CGtype

! contains

!     subroutine set_MF(obj, MFfun)

!         implicit none
!         ! Input / output data
!         ! -------------------
!         type(CGtype), pointer :: obj
!         external :: MFfun
!         obj%matrixfree => MFfun

!     end subroutine set_MF

!     subroutine set_FD(obj, FDfun)

!         implicit none
!         ! Input / output data
!         ! -------------------
!         type(CGtype), pointer :: obj
!         external :: FDfun
!         obj%fastdiag => FDfun

!     end subroutine set_FD

!     subroutine PBiCGSTAB(obj, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
!                         indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
!                         data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                         data_W_u, data_W_v, data_W_w, U_u, U_v, U_w, D_u, D_v, D_w, &
!                         nbIterPCG, threshold, b, x, RelRes)

!         implicit none
!         ! Input / output data
!         ! -------------------
!         type(CGtype), pointer :: obj
!         integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
!         integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
!         dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
!         integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
!         dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
!         double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
!         dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

!         integer, intent(in) :: indi_u, indi_v, indi_w
!         dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
!         integer, intent(in) :: indj_u, indj_v, indj_w
!         dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
!         double precision, intent(in) :: data_W_u, data_W_v, data_W_w
!         dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

!         double precision, intent(in) :: U_u, U_v, U_w, D_u, D_v, D_w
!         dimension ::    U_u(nr_u, nr_u), U_v(nr_v, nr_v), U_w(nr_w, nr_w), &
!                         D_u(nr_u), D_v(nr_v), D_w(nr_w)


!         integer, intent(in) :: nbIterPCG
!         double precision, intent(in) :: threshold, b
!         dimension :: b(nr_total)
        
!         double precision, intent(out) :: x, RelRes
!         dimension :: x(nr_total), RelRes(nbIterPCG+1)

!         ! Local data
!         ! -----------
!         ! Conjugate gradient algorithm
!         double precision :: rsold, rsnew, alpha, omega, beta, normb
!         double precision :: r, rhat, p, s, ptilde, Aptilde, Astilde, stilde
!         dimension ::    r(nr_total), rhat(nr_total), p(nr_total), s(nr_total), &
!                         ptilde(nr_total), Aptilde(nr_total), Astilde(nr_total), stilde(nr_total)
!         integer :: iter

!         x = 0.d0; r = b; rhat = r; p = r
!         rsold = dot_product(r, rhat); normb = maxval(abs(r))
!         RelRes = 0.d0; RelRes(1) = 1.d0

!         do iter = 1, nbIterPCG
!             call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, p, ptilde)
!             call mf_wq_get_cu_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                         nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
!                         data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                         data_W_u, data_W_v, data_W_w, ptilde, Aptilde)
!             alpha = rsold/dot_product(Aptilde, rhat)
!             s = r - alpha*Aptilde
            
!             call fd_interpolation_3d(nr_total, nr_u, nr_v, nr_w, U_u, U_v, U_w, s, stilde)
!             call mf_wq_get_cu_3D(coefs, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
!                         nnz_u, nnz_v, nnz_w, indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
!                         data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
!                         data_W_u, data_W_v, data_W_w, stilde, Astilde)
!             omega = dot_product(Astilde, s)/dot_product(Astilde, Astilde)
!             x = x + alpha*ptilde + omega*stilde
!             r = s - omega*Astilde    
            
!             RelRes(iter+1) = maxval(abs(r))/normb
!             if (RelRes(iter+1).le.threshold) exit
    
!             rsnew = dot_product(r, rhat)
!             beta = (alpha/omega)*(rsnew/rsold)
!             p = r + beta*(p - omega*Aptilde)
!             rsold = rsnew
!         end do

!     end subroutine 

! end module CGsolver