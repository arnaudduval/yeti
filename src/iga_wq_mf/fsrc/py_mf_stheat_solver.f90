! ====================================================
! This modules aims to compute the dot product between a matrix and a vector 
! exploiting the tensor-product structure of those matrices.
! Moreover, it uses weighted quadrature in order to reduce the number of quadrature points.
! It is implemented conjugated gradient algorithms to solve interpolation problems
! ====================================================

subroutine mf_stcapacity_2d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, &
                        data_W_t, invJ, detJ, detG, prop, array_in, array_out)

    use matrixfreestheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 2
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    double precision, intent(in) :: invJ, detJ, detG, prop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(stthermomat) :: mat
    type(basis_data) :: basisdata
    call init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                        indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_B_u, data_B_v, data_B_t, &
                        data_W_u, data_W_v, data_W_t)
    call getcsrc2dense(basisdata)
    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(prop), prop)
    call mf_u_partialt_v(mat, basisdata, nr_total, array_in, array_out)

end subroutine mf_stcapacity_2d

subroutine mf_stcapacity_3d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, &
                        data_W_w, data_W_t, invJ, detJ, detG, prop, array_in, array_out)

    use matrixfreestheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, &
                            nc_w, nr_t, nc_t, nnz_u, nnz_v, nnz_w, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    double precision, intent(in) :: invJ, detJ, detG, prop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), prop(nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data 
    ! ---------- 
    type(stthermomat) :: mat
    type(basis_data) :: basisdata
    call init_4basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                        data_W_u, data_W_v, data_W_w, data_W_t)
    call getcsrc2dense(basisdata)
    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(prop), prop)
    call mf_u_partialt_v(mat, basisdata, nr_total, array_in, array_out)

end subroutine mf_stcapacity_3d

subroutine mf_stconductivity_2d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, &
                        data_W_t, invJ, detJ, detG, prop, array_in, array_out)

    use matrixfreestheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 2
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    double precision, intent(in) :: invJ, detJ, detG, prop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), &
                detG(nc_tm), prop(dimen_sp, dimen_sp, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(basis_data) :: basisdata
    call init_3basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                        indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_B_u, data_B_v, data_B_t, &
                        data_W_u, data_W_v, data_W_t)
    call getcsrc2dense(basisdata)
    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_gradx_u_gradx_v(mat, basisdata, nr_total, array_in, array_out)

end subroutine mf_stconductivity_2d

subroutine mf_stconductivity_3d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                        indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, &
                        data_W_w, data_W_t, invJ, detJ, detG, prop, array_in, array_out)

    use matrixfreestheat
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                            nnz_u, nnz_v, nnz_w, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    double precision, intent(in) :: invJ, detJ, detG, prop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), &
                detG(nc_tm), prop(dimen_sp, dimen_sp, nc_total)
    double precision, intent(in) :: array_in
    dimension :: array_in(nr_total)

    double precision, intent(out) :: array_out
    dimension :: array_out(nr_total)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(basis_data) :: basisdata
    call init_4basisdata(basisdata, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                        data_W_u, data_W_v, data_W_w, data_W_t)
    call getcsrc2dense(basisdata)

    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_conductivityprop(mat, nc_total, prop)
    call mf_gradx_u_gradx_v(mat, basisdata, nr_total, array_in, array_out)

end subroutine mf_stconductivity_3d

subroutine solver_picardspacetime_heat_2d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                                indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, &
                                data_W_t, table, invJ, detJ, detG, Cprop, Kprop, Fext, &
                                iterations, threshold, linprecond, linsolver, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreestheat
    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 2
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    logical, intent(in) :: table
    dimension :: table(dimen_sp+1, 2) 

    double precision, intent(in) :: invJ, detJ, detG, Cprop, Kprop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), &
                Cprop(nc_total), Kprop(dimen_sp, dimen_sp, nc_total)
    character(len=10), intent(in) :: linprecond, linsolver
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen_sp+1), nbRestarts=1
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
    integer :: dimen 

    dimen = dimen_sp + 1
    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                        indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_B_u, data_B_v, data_B_t, &
                        data_W_u, data_W_v, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.

    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(Cprop), Cprop)
    call setup_conductivityprop(mat, size(Kprop, dim=3), Kprop)
    solv%matrixfreetype = 1
    nc_list = (/nc_u, nc_v, nc_t/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen_sp, mat%Kmean(1:dimen_sp))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        solv%scalarleft = mat%Cmean
        if (solv%applyfd) then
            call space_eigendecomposition(redsyst)
            call time_schurdecomposition(redsyst)
        end if
        call initialize_solver(solv, globsyst, redsyst)
        if (linsolver.eq.'BICG') then
            call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
        else if (linsolver.eq.'GMRES') then
            call PGMRES(solv, mat, nr_total,  nbRestarts, iterations, threshold, Fext, x, residual)
        else
            stop 'Unknown method'
        end if
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_picardspacetime_heat_2d

subroutine solver_newtonspacetime_heat_2d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                                indi_t, indj_t, data_B_u, data_B_v, data_B_t, data_W_u, data_W_v, &
                                data_W_t, table, invJ, detJ, detG, Cprop, Cdersprop, Kprop, Kdersprop, &
                                Fext, iterations, threshold, linprecond, linsolver, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreestheat
    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 2
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    logical, intent(in) :: table
    dimension :: table(dimen_sp+1, 2) 

    double precision, intent(in) :: invJ, detJ, detG, Cprop, Cdersprop, Kprop, Kdersprop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), &
                Cprop(nc_total), Cdersprop(nc_total), &
                Kprop(dimen_sp, dimen_sp, nc_total), Kdersprop(dimen_sp, nc_total)
    character(len=10), intent(in) :: linprecond, linsolver
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen_sp+1), nbRestarts=1
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
    integer :: dimen

    dimen = dimen_sp + 1
    call init_3basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_t, nc_t, nnz_u, nnz_v, nnz_t, &
                        indi_u, indj_u, indi_v, indj_v, indi_t, indj_t, data_B_u, data_B_v, data_B_t, &
                        data_W_u, data_W_v, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.

    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(Cprop), Cprop)
    call setup_capacityDersprop(mat, size(Cdersprop), Cdersprop)
    call setup_conductivityprop(mat, size(Kprop, dim=3), Kprop)
    call setup_conductivityDersprop(mat, size(Kdersprop, dim=2), Kdersprop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v, nc_t/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen_sp, mat%Kmean(1:dimen_sp))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        solv%scalarleft = mat%Cmean
        if (solv%applyfd) then
            call space_eigendecomposition(redsyst)
            call time_schurdecomposition(redsyst)
        end if
        call initialize_solver(solv, globsyst, redsyst)
        if (linsolver.eq.'BICG') then
            call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
        else if (linsolver.eq.'GMRES') then
            call PGMRES(solv, mat, nr_total, nbRestarts, iterations, threshold, Fext, x, residual)
        else
            stop 'Unknown method'
        end if
        
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_newtonspacetime_heat_2d

subroutine solver_picardspacetime_heat_3d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, &
                                data_W_t, table, invJ, detJ, detG, Cprop, Kprop, Fext, iterations, threshold, linprecond, &
                                linsolver, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreestheat
    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, &
                            nc_t, nnz_u, nnz_v, nnz_w, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    logical, intent(in) :: table
    dimension :: table(dimen_sp+1, 2) 

    double precision, intent(in) :: invJ, detJ, detG, Cprop, Kprop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), &
                Cprop(nc_total), Kprop(dimen_sp, dimen_sp, nc_total)
    character(len=10), intent(in) :: linprecond, linsolver
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen_sp+1), nbRestarts=1
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
    integer :: dimen 

    dimen = dimen_sp + 1
    call init_4basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                        data_W_u, data_W_v, data_W_w, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.
    
    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(Cprop), Cprop)
    call setup_conductivityprop(mat, size(Kprop, dim=3), Kprop)
    solv%matrixfreetype = 1
    nc_list = (/nc_u, nc_v, nc_w, nc_t/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen_sp, mat%Kmean(1:dimen_sp))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        solv%scalarleft = mat%Cmean
        if (solv%applyfd) then
            call space_eigendecomposition(redsyst)
            call time_schurdecomposition(redsyst)
        end if
        call initialize_solver(solv, globsyst, redsyst)
        if (linsolver.eq.'BICG') then
            call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
        else if (linsolver.eq.'GMRES') then
            call PGMRES(solv, mat, nr_total,  nbRestarts, iterations, threshold, Fext, x, residual)
        else
            stop 'Unknown method'
        end if
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_picardspacetime_heat_3d

subroutine solver_newtonspacetime_heat_3d(nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                                nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, data_W_u, data_W_v, data_W_w, &
                                data_W_t, table, invJ, detJ, detG, Cprop, Cdersprop, Kprop, Kdersprop, &
                                Fext, iterations, threshold, linprecond, linsolver, x, residual)
    !! Precontionned bi-conjugate gradient to solve transient heat problems
    !! It solves Ann un = bn, where Ann is (thetadt*Knn + Cnn) and bn = Fn - And ud
    !! bn is compute beforehand (In python or fortran).
    !! IN CSR FORMAT

    use matrixfreestheat
    use stheatsolver
    use structured_data
    implicit none 
    ! Input / output data
    ! -------------------
    integer, parameter :: dimen_sp = 3
    integer, intent(in) :: nr_total, nc_total, nc_sp, nc_tm, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, &
                            nc_t, nnz_u, nnz_v, nnz_w, nnz_t    
    integer, intent(in) :: indi_u, indi_v, indi_w, indi_t
    dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1), indi_t(nr_t+1)
    integer, intent(in) ::  indj_u, indj_v, indj_w, indj_t
    dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w), indj_t(nnz_t)
    double precision, intent(in) :: data_B_u, data_W_u, data_B_v, data_W_v, data_B_w, data_W_w, data_B_t, data_W_t
    dimension ::    data_B_u(nnz_u, 2), data_W_u(nnz_u, 4), &
                    data_B_v(nnz_v, 2), data_W_v(nnz_v, 4), &
                    data_B_w(nnz_w, 2), data_W_w(nnz_w, 4), &
                    data_B_t(nnz_t, 2), data_W_t(nnz_t, 4)
    logical, intent(in) :: table
    dimension :: table(dimen_sp+1, 2) 

    double precision, intent(in) :: invJ, detJ, detG, Cprop, Cdersprop, Kprop, Kdersprop
    dimension :: invJ(dimen_sp, dimen_sp, nc_sp), detJ(nc_sp), detG(nc_tm), &
                Cprop(nc_total), Cdersprop(nc_total), &
                Kprop(dimen_sp, dimen_sp, nc_total), Kdersprop(dimen_sp, nc_total)
    character(len=10), intent(in) :: linprecond, linsolver
    integer, intent(in) :: iterations    
    double precision, intent(in) :: threshold

    double precision, intent(in) :: Fext
    dimension :: Fext(nr_total)
    
    double precision, intent(out) :: x, residual
    dimension :: x(nr_total), residual(iterations+1)

    ! Local data
    ! ----------
    type(stthermomat) :: mat
    type(stcgsolver) :: solv
    type(basis_data), target :: globsyst
    type(reduced_system), target :: redsyst
    integer :: nc_list(dimen_sp+1), nbRestarts=1
    double precision, allocatable, dimension(:, :) :: univMcoefs, univKcoefs
    integer :: dimen

    dimen = dimen_sp + 1
    call init_4basisdata(globsyst, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nr_t, nc_t, &
                        nnz_u, nnz_v, nnz_w, nnz_t, indi_u, indj_u, indi_v, indj_v, &
                        indi_w, indj_w, indi_t, indj_t, data_B_u, data_B_v, data_B_w, data_B_t, &
                        data_W_u, data_W_v, data_W_w, data_W_t)
    call copybasisdata(globsyst, redsyst%basisdata)
    call update_reducedsystem(redsyst, dimen, table)
    call getcsrc2dense(redsyst%basisdata)
    call getcsrc2dense(globsyst)
    redsyst%isspacetime = .true.

    call setup_geometry(mat, dimen_sp, size(detJ), size(detG), invJ, detJ, detG)
    call setup_capacityprop(mat, size(Cprop), Cprop)
    call setup_capacityDersprop(mat, size(Cdersprop), Cdersprop)
    call setup_conductivityprop(mat, size(Kprop, dim=3), Kprop)
    call setup_conductivityDersprop(mat, size(Kdersprop, dim=2), Kdersprop)
    solv%matrixfreetype = 2
    nc_list = (/nc_u, nc_v, nc_w, nc_t/)

    if ((linprecond.eq.'WP').or.(linprecond.eq.'JMC').or.(linprecond.eq.'C').or.(linprecond.eq.'TDC')) then

        if (linprecond.eq.'WP') then
            solv%applyfd = .false.
        end if

        if (linprecond.eq.'JMC') then 
            call compute_variablesmean(mat, nc_list)
            call setup_meancoefs(redsyst, dimen_sp, mat%Kmean(1:dimen_sp))
        end if

        if (linprecond.eq.'TDC') then
            allocate(univMcoefs(dimen, maxval(nc_list)), univKcoefs(dimen, maxval(nc_list)))
            call compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
            call setup_univariatecoefs(redsyst, dimen, maxval(nc_list), univMcoefs, univKcoefs)
        end if

        solv%scalarleft = mat%Cmean
        if (solv%applyfd) then
            call space_eigendecomposition(redsyst)
            call time_schurdecomposition(redsyst)
        end if
        call initialize_solver(solv, globsyst, redsyst)
        if (linsolver.eq.'BICG') then
            call PBiCGSTAB(solv, mat, nr_total, iterations, threshold, Fext, x, residual)
        else if (linsolver.eq.'GMRES') then
            call PGMRES(solv, mat, nr_total,  nbRestarts, iterations, threshold, Fext, x, residual)
        else
            stop 'Unknown method'
        end if
    else 
        stop 'Unknown method' 
    end if

end subroutine solver_newtonspacetime_heat_3d