module matrixfreeplasticity

    implicit none
    type :: mecamat
    
        integer :: dimen, nvoigt
        logical :: isLumped = .false., isElastic = .true.
        double precision :: scalars(2) = (/1.d0, 1.d0/)
        double precision, dimension(:), pointer :: detJ=>null(), Mprop=>null(), Hprop=>null()
        double precision, dimension(:, :), pointer :: CepArgs=>null(), NN=>null(), BB=>null()
        double precision, dimension(:, :, :), pointer :: invJ=>null()
        double precision, dimension(:, :, :), allocatable :: JJjj, JJnn, JJbb
        double precision, dimension(:, :), allocatable :: Smean
        double precision, dimension(:), allocatable :: Mmean
        integer :: ncols_sp
    
    end type mecamat

contains

    subroutine setup_geometry(mat, dimen, nnz, invJ, detJ)
        !! Points to the data of the inverse and determinant of the Jacobian. 
        !! It also computes and saves inv(JJ) inv(JJ).transpose

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: dimen, nnz
        double precision, target, intent(in) :: invJ, detJ
        dimension :: invJ(dimen, dimen, nnz), detJ(nnz)

        mat%dimen = dimen
        mat%nvoigt = dimen*(dimen+1)/2
        mat%invJ => invJ
        mat%detJ => detJ
        mat%ncols_sp = nnz

        allocate(mat%Smean(mat%dimen, mat%dimen), mat%Mmean(mat%dimen))
        mat%Smean = 1.d0; mat%Mmean = 1.d0

    end subroutine setup_geometry

    subroutine setup_massprop(mat, nnz, prop)

        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Mprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_massprop

    subroutine setup_thmchcoupledprop(mat, nnz, prop)
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat
        integer, intent(in) :: nnz
        double precision, target, intent(in) ::  prop
        dimension :: prop(nnz)

        mat%Hprop => prop
        mat%ncols_sp = nnz

    end subroutine setup_thmchcoupledprop

    subroutine setup_mechanicalArguments(mat, nbrows, mechArgs)
        !! Points to data of the mechanical behavior 

        implicit none
        ! Input / output data
        ! -------------------
        double precision, parameter :: threshold = 1.d-12
        type(mecamat) :: mat
        integer, intent(in) :: nbrows
        double precision, target, intent(in) :: mechArgs
        dimension :: mechArgs(nbrows, mat%ncols_sp)

        ! Local data
        ! ----------
        double precision :: tensor(mat%dimen, mat%dimen)
        integer :: i

        if (nbrows.gt.2) mat%isElastic = .false.
        if (mat%isElastic) then
            mat%CepArgs => mechArgs(:2, :)
        else
            mat%CepArgs => mechArgs(:4, :)    
            mat%NN      => mechArgs(5:5+mat%nvoigt, :)
            mat%BB      => mechArgs(5+mat%nvoigt:5+2*mat%nvoigt, :)

            if (.not.allocated(mat%JJnn)) allocate(mat%JJnn(mat%dimen, mat%dimen, mat%ncols_sp))
            if (.not.allocated(mat%JJbb)) allocate(mat%JJbb(mat%dimen, mat%dimen, mat%ncols_sp))
            mat%JJnn = 0.d0; mat%JJbb = 0.d0
            if (all(abs(mat%NN).lt.threshold)) return
            do i = 1, mat%ncols_sp
                call array2symtensor(mat%dimen, mat%nvoigt, mat%NN(:, i), tensor)
                mat%JJnn(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
                call array2symtensor(mat%dimen, mat%nvoigt, mat%BB(:, i), tensor)
                mat%JJbb(:, :, i) = matmul(mat%invJ(:, :, i), tensor)
            end do
        end if

    end subroutine setup_mechanicalArguments

    subroutine setup_jacobienjacobien(mat)
        
        implicit none
        ! Input / output data
        ! -------------------
        type(mecamat) :: mat

        ! Local data
        ! ----------
        integer :: i
        
        if (.not.allocated(mat%JJjj)) allocate(mat%JJjj(mat%dimen, mat%dimen, mat%ncols_sp))
        mat%JJjj = 0.d0
        do i = 1, mat%ncols_sp
            mat%JJjj(:, :, i) = matmul(mat%invJ(:, :, i), transpose(mat%invJ(:, :, i)))
        end do

    end subroutine setup_jacobienjacobien

    ! subroutine compute_separationvariables(mat, nc_list, univMcoefs, univKcoefs)
        
    !     use separatevariables
    !     implicit none 
    !     ! Input / output data
    !     ! -------------------
    !     type(mecamat) :: mat
    !     integer, intent(in) :: nc_list
    !     dimension :: nc_list(mat%dimen)

    !     double precision, intent(out) :: univMcoefs(mat%dimen, mat%dimen, maxval(nc_list)), &
    !                                     univKcoefs(mat%dimen, mat%dimen, maxval(nc_list))

    !     ! Local data
    !     ! ----------
    !     type(sepoperator) :: oper
    !     integer :: i, j, k, gp
    !     integer, allocatable, dimension(:) :: nc_list_t
    !     logical, allocatable, dimension(:) :: update
    !     double precision, allocatable, dimension(:, :, :) :: CC
    !     double precision :: DD, NN, TNN, BB, TBB
    !     dimension :: DD(mat%dimen, mat%dimen), NN(mat%nvoigt), TNN(mat%dimen, mat%dimen), &
    !                 BB(mat%nvoigt), TBB(mat%dimen, mat%dimen)

    !     if (.not.associated(mat%Mprop)) then
    !         allocate(CC(mat%dimen, mat%dimen, mat%ncols_sp), update(mat%dimen), nc_list_t(mat%dimen))
    !         update = .true.; nc_list_t = nc_list
    !         call initialize_operator(oper, mat%dimen, nc_list_t, update)

    !         do i = 1, mat%dimen
    !             do gp = 1, mat%ncols_sp
    !                 NN = mat%NN(:, gp); BB = mat%BB(:, gp)
    !                 call array2symtensor(mat%dimen, size(NN), NN, TNN)
    !                 call array2symtensor(mat%dimen, size(BB), BB, TBB)
                    
    !                 ! Elastic
    !                 DD = 0.d0
    !                 DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
    !                 do k = 1, mat%dimen
    !                     DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
    !                 end do

    !                 ! Plastic
    !                 if (mat%isElastic.eqv..false.) then
    !                     do j = 1, mat%dimen
    !                         do k = 1, mat%dimen
    !                             DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
    !                                         + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
    !                         end do
    !                     end do
    !                 end if
    !                 CC(:, :, gp) = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
    !             end do

    !             if (mat%dimen.eq.2) then
    !                 call separatevariables_2d(oper, CC)
    !             else if (mat%dimen.eq.3) then
    !                 call separatevariables_3d(oper, CC)
    !             end if
    !             univMcoefs(i, :, :) = oper%univmasscoefs; univKcoefs(i, :, :) = oper%univstiffcoefs
    !         end do
    !     else
    !         allocate(CC(mat%dimen+1, mat%dimen+1, mat%ncols_sp), update(mat%dimen+1), nc_list_t(mat%dimen+1))
    !         update = .true.; update(mat%dimen+1) = .false.
    !         nc_list_t(:mat%dimen) = nc_list; nc_list_t(mat%dimen+1) = 1
    !         call initialize_operator(oper, mat%dimen+1, nc_list_t, update)

    !         do i = 1, mat%dimen
    !             do gp = 1, mat%ncols_sp
    !                 NN = mat%NN(:, gp); BB = mat%BB(:, gp)
    !                 call array2symtensor(mat%dimen, size(NN), NN, TNN)
    !                 call array2symtensor(mat%dimen, size(BB), BB, TBB)
                    
    !                 ! Elastic
    !                 DD = 0.d0
    !                 DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
    !                 do k = 1, mat%dimen
    !                     DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
    !                 end do
                    
    !                 ! Plastic
    !                 if (mat%isElastic.eqv..false.) then
    !                     do j = 1, mat%dimen
    !                         do k = 1, mat%dimen
    !                             DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
    !                                         + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
    !                         end do
    !                     end do
    !                 end if
    !                 CC(:mat%dimen, :mat%dimen, gp) = matmul(mat%invJ(:, :, gp), &
    !                                                 matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
    !                 CC(mat%dimen+1, mat%dimen+1, gp) = mat%Mprop(gp)*mat%detJ(gp)
    !             end do

    !             if (mat%dimen.eq.2) then
    !                 call separatevariables_3d(oper, CC)
    !             else if (mat%dimen.eq.3) then
    !                 call separatevariables_4d(oper, CC)
    !             end if
    !             univMcoefs(i, :, :) = oper%univmasscoefs(:mat%dimen, :); univKcoefs(i, :, :) = oper%univstiffcoefs(:mat%dimen, :)
    !         end do
            
    !     end if
    ! end subroutine compute_separationvariables

    ! subroutine compute_mean(mat, nclist)
    !     !! Computes the average of the material properties (for the moment it only considers elastic materials)

    !     implicit none 
    !     ! Input / output data
    !     ! -------------------
    !     type(mecamat) :: mat
    !     integer, intent(in) :: nclist
    !     dimension :: nclist(mat%dimen)

    !     ! Local data
    !     ! ----------
    !     integer :: i, j, k, c, gp, pos, ind(3)
    !     integer, dimension(:), allocatable :: sample
    !     integer, dimension(:, :), allocatable :: indlist
    !     double precision :: DD, TNN, NN, TBB, BB, Mmean
    !     dimension :: DD(mat%dimen, mat%dimen), TNN(mat%dimen, mat%dimen), NN(mat%nvoigt), &
    !                 TBB(mat%dimen, mat%dimen), BB(mat%nvoigt)
    !     double precision, allocatable, dimension(:, :, :) :: Scoefs
    !     double precision, allocatable, dimension(:) :: Mcoefs

    !     if (product(nclist).ne.mat%ncols_sp) stop 'Size problem'
    !     allocate(indlist(mat%dimen, 3), sample(3**mat%dimen))
    !     do i = 1, mat%dimen 
    !         pos = int((nclist(i) + 1)/2); ind = (/1, pos, nclist(i)/)
    !         indlist(i, :) = ind
    !     end do
    
    !     ! Select a set of coefficients
    !     c = 1
    !     if (mat%dimen.eq.2) then
    !         do j = 1, 3
    !             do i = 1, 3
    !                 gp = indlist(1, i) + (indlist(2, j) - 1)*nclist(1)
    !                 sample(c) = gp
    !                 c = c + 1
    !             end do
    !         end do
    !     else if (mat%dimen.eq.3) then
    !         do k = 1, 3
    !             do j = 1, 3
    !                 do i = 1, 3
    !                     gp = indlist(1, i) + (indlist(2, j) - 1)*nclist(1) + (indlist(3, k) - 1)*nclist(1)*nclist(2)
    !                     sample(c) = gp
    !                     c = c + 1
    !                 end do
    !             end do
    !         end do
    !     else
    !         stop 'Try 2 or 3 dimensions'
    !     end if
        
    !     allocate(Scoefs(mat%dimen, mat%dimen, size(sample)))
    !     do i = 1, mat%dimen
    !         do c = 1, size(sample)
    !             gp = sample(c)
    !             NN = mat%NN(:, gp); BB = mat%NN(:, gp)
    !             call array2symtensor(mat%dimen, size(NN), NN, TNN)
    !             call array2symtensor(mat%dimen, size(BB), BB, TBB)

    !             ! Elastic
    !             DD = 0.d0
    !             DD(i, i) = DD(i, i) + mat%CepArgs(1, gp) + mat%CepArgs(2, gp)
    !             do k = 1, mat%dimen
    !                 DD(k, k) = DD(k, k) + mat%CepArgs(2, gp)
    !             end do
                
    !             ! Plastic
    !             if (mat%isElastic.eqv..false.) then
    !                 do j = 1, mat%dimen
    !                     do k = 1, mat%dimen
    !                         DD(j, k) = DD(j, k) + mat%CepArgs(3, gp)*TNN(i, j)*TNN(i, k) &
    !                                     + mat%CepArgs(4, gp)*(TBB(i, j)*TNN(i, k) - TNN(i, j)*TBB(i, k))
    !                     end do
    !                 end do
    !             end if

    !             Scoefs(:, :, c) = matmul(mat%invJ(:, :, gp), matmul(DD, transpose(mat%invJ(:, :, gp))))*mat%detJ(gp)
    !         end do
    
    !         do j = 1, mat%dimen
    !             if (mat%dimen.eq.2) then
    !                 call trapezoidal_rule_2d(3, 3, Scoefs(j, j, :), mat%Smean(i, j))
    !             else if (mat%dimen.eq.3) then
    !                 call trapezoidal_rule_3d(3, 3, 3, Scoefs(j, j, :), mat%Smean(i, j))
    !             end if
    !         end do   
    !     end do

    !     if (associated(mat%Mprop)) then
    !         allocate(Mcoefs(size(sample)))
    !         do c = 1, size(sample)
    !             gp = sample(c)
    !             Mcoefs(c) = mat%Mprop(gp)*mat%detJ(gp)
    !         end do
    !         if (mat%dimen.eq.2) then
    !             call trapezoidal_rule_2d(3, 3, Mcoefs, Mmean)
    !         else if (mat%dimen.eq.3) then
    !             call trapezoidal_rule_3d(3, 3, 3, Mcoefs, Mmean)
    !         end if 
    !         mat%Mmean = Mmean
    !     end if

    ! end subroutine compute_mean

    subroutine mf_tu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i
        double precision :: tmp_in, array_tmp, array_tmp2
        dimension :: tmp_in(nr_total), array_tmp(nc_total), array_tmp2(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        do i = 1, dimen
            tmp_in = array_in(i, :); if (mat%isLumped) tmp_in = 1.d0

            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                tmp_in, array_tmp)

            array_tmp = array_tmp*mat%Mprop*mat%detJ

            call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                array_tmp, array_tmp2)

            array_out(i, :) = array_tmp2; if (mat%isLumped) array_out(i, :) = array_tmp2*array_in(i, :)

        end do

    end subroutine mf_tu_tv_2d

    subroutine mf_tu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i
        double precision :: tmp_in, array_tmp, array_tmp2
        dimension :: tmp_in(nr_total), array_tmp(nc_total), array_tmp2(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'

        do i = 1, dimen
            tmp_in = array_in(i, :); if (mat%isLumped) tmp_in = 1.d0

            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                                nnz_u, indi_T_u, indj_T_u, data_BT_u(:, 1), & 
                                nnz_v, indi_T_v, indj_T_v, data_BT_v(:, 1), &
                                nnz_w, indi_T_w, indj_T_w, data_BT_w(:, 1), & 
                                tmp_in, array_tmp)

            array_tmp = array_tmp*mat%Mprop*mat%detJ

            call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, 1), &
                                nnz_v, indi_v, indj_v, data_W_v(:, 1), &
                                nnz_w, indi_w, indj_w, data_W_w(:, 1), &
                                array_tmp, array_tmp2)
            array_out(i, :) = array_tmp2; if (mat%isLumped) array_out(i, :) = array_tmp2*array_in(i, :)
        end do
            
    end subroutine mf_tu_tv_3d

    subroutine mf_gradtu_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                            data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 2D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, j, k, l, alpha, beta, zeta, nbCepArgs = 4
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision, allocatable, dimension(:) :: t4, t5, t6
        double precision, allocatable, dimension(:, :) :: kt1
        double precision :: t1, t2, t3, t7, t8, t9
        dimension :: t1(nc_total), t2(nc_total), t3(nc_total), &
                        t7(nc_total), t8(nr_total), t9(nr_total)

        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        array_out = 0.d0

        if (mat%isElastic) then 
            nbCepArgs = 2
        else
            allocate(t4(nc_total), t5(nc_total), t6(nc_total))
        end if
        allocate(kt1(nbCepArgs, nc_total))

        do j = 1, dimen
            do l = 1, dimen
                beta = 1; beta(l) = 2
                call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in(j, :), t1) 

                do i = 1, nbCepArgs
                    kt1(i, :) = mat%CepArgs(i, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(l, j, :)
                if (mat%isElastic) then
                    t4 = kt1(3, :)*mat%JJnn(l, j, :)
                    t5 = kt1(4, :)*mat%JJbb(l, j, :)
                    t6 = kt1(4, :)*mat%JJnn(l, j, :)
                end if

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(l, i, :)
                    t9 = 0.d0

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta  = beta + (alpha - 1)*2
                        
                        t7 = t2*mat%invJ(k, i, :) + t3*mat%invJ(k, j, :)
                        if (mat%isElastic.eqv..false.) then
                            t7 = t7 + t4*mat%JJnn(k, i, :) - t5*mat%JJnn(k, i, :) + t6*mat%JJbb(k, i, :)
                        end if
                        if (i.eq.j) t7 = t7 + kt1(2, :)*mat%JJjj(k, l, :)
                        call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                t7, t8)
                        t9 = t9 + t8
                    end do

                    array_out(i, :) = array_out(i, :) + t9
                end do
            end do
        end do
            
    end subroutine mf_gradtu_gradtv_2d

    subroutine mf_gradtu_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, j, k, l, alpha, beta, zeta, nbCepArgs = 4
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision, allocatable, dimension(:) :: t4, t5, t6
        double precision, allocatable, dimension(:, :) :: kt1
        double precision :: t1, t2, t3, t7, t8, t9
        dimension :: t1(nc_total), t2(nc_total), t3(nc_total), &
                        t7(nc_total), t8(nc_total), t9(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        array_out = 0.d0
        
        if (mat%isElastic) then 
            nbCepArgs = 2
        else
            allocate(t4(nc_total), t5(nc_total), t6(nc_total))
        end if
        allocate(kt1(nbCepArgs, nc_total))

        do j = 1, dimen
            do l = 1, dimen
                beta = 1; beta(l) = 2
                call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in(j, :), t1) 

                do i = 1, nbCepArgs
                    kt1(i, :) = mat%CepArgs(i, :)*t1*mat%detJ
                end do

                t2 = kt1(1, :)*mat%invJ(l, j, :)
                if (mat%isElastic) then
                    t4 = kt1(3, :)*mat%JJnn(l, j, :)
                    t5 = kt1(4, :)*mat%JJbb(l, j, :)
                    t6 = kt1(4, :)*mat%JJnn(l, j, :)
                end if

                do i = 1, dimen
                    t3 = kt1(2, :)*mat%invJ(l, i, :)
                    t9 = 0.d0

                    do k = 1, dimen
                        alpha = 1; alpha(k) = 2
                        zeta  = beta + (alpha - 1)*2
                        
                        t7 = t2*mat%invJ(k, i, :) + t3*mat%invJ(k, j, :)
                        if (mat%isElastic.eqv..false.) then
                            t7 = t7 + t4*mat%JJnn(k, i, :) - t5*mat%JJnn(k, i, :) + t6*mat%JJbb(k, i, :)
                        end if
                        if (i.eq.j) t7 = t7 + kt1(2, :)*mat%JJjj(k, l, :)
                        call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t8, t9)
                        t9 = t9 + t8
                    end do

                    array_out(i, :) = array_out(i, :) + t9
                end do
            end do
        end do
            
    end subroutine mf_gradtu_gradtv_3d

    subroutine mf_tutv_gradtugradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        !! Computes S.u in 2D where S is stiffness matrix
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(dimen, nr_total)

        call mf_tu_tv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradtu_gradtv_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                        indi_T_u, indj_T_u, indi_T_v, indj_T_v, &
                        data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                        data_W_u, data_W_v, array_in, array_tmp)

        array_out = array_out + mat%scalars(2)*array_tmp

    end subroutine mf_tutv_gradtugradtv_2d

    subroutine mf_tutv_gradtugradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)

        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(dimen, nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ---------------
        double precision :: array_tmp
        dimension :: array_tmp(dimen, nr_total)

        call mf_tu_tv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                    indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                    data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                    data_W_u, data_W_v, data_W_w, array_in, array_out)

        array_out = mat%scalars(1)*array_out

        call mf_gradtu_gradtv_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                                indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                                data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                                data_W_u, data_W_v, data_W_w, array_in, array_tmp)
        
        array_out = array_out + mat%scalars(2)*array_tmp

    end subroutine mf_tutv_gradtugradtv_3d

    subroutine mf_gradu_v_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, data_BT_u, data_BT_v, indi_u, indj_u, indi_v, indj_v, &
                            data_W_u, data_W_v, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v

        integer, intent(in) :: indi_T_u, indi_T_v
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1)
        integer, intent(in) :: indj_T_u, indj_T_v
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v)
        double precision, intent(in) :: data_BT_u, data_BT_v
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2)

        integer, intent(in) :: indi_u, indi_v
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1)
        integer, intent(in) :: indj_u, indj_v
        dimension :: indj_u(nnz_u), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)
        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        array_out = 0.d0     
        do l = 1, dimen
            beta = 1; beta(l) = 2
            call sumfacto2d_spM(nc_u, nr_u, nc_v, nr_v, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        array_in, t1) 
            t1 = t1*mat%Hprop*mat%detJ
            do i = 1, dimen
                alpha = 1; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(l, i, :)
                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), t2, t3)
                array_out(i, :) = array_out(i, :) + t3
            end do
        end do
            
    end subroutine mf_gradu_v_2d

    subroutine mf_gradu_v_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_T_u, indj_T_u, indi_T_v, indj_T_v, indi_T_w, indj_T_w, &
                            data_BT_u, data_BT_v, data_BT_w, indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, &
                            data_W_u, data_W_v, data_W_w, array_in, array_out)
        implicit none 
        ! Input / output data 
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w

        integer, intent(in) :: indi_T_u, indi_T_v, indi_T_w
        dimension :: indi_T_u(nc_u+1), indi_T_v(nc_v+1), indi_T_w(nc_w+1)
        integer, intent(in) :: indj_T_u, indj_T_v, indj_T_w
        dimension :: indj_T_u(nnz_u), indj_T_v(nnz_v), indj_T_w(nnz_w)
        double precision, intent(in) :: data_BT_u, data_BT_v, data_BT_w
        dimension :: data_BT_u(nnz_u, 2), data_BT_v(nnz_v, 2), data_BT_w(nnz_w, 2)

        integer, intent(in) :: indi_u, indi_v, indi_w
        dimension :: indi_u(nr_u+1), indi_v(nr_v+1), indi_w(nr_w+1)
        integer, intent(in) :: indj_u, indj_v, indj_w
        dimension :: indj_u(nnz_u), indj_v(nnz_v), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)

        double precision, intent(in) :: array_in
        dimension :: array_in(nr_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data 
        ! ----------
        integer :: i, l, alpha, beta, zeta
        dimension :: alpha(dimen), beta(dimen), zeta(dimen)
        double precision :: t1, t2, t3
        dimension :: t1(nc_total), t2(nc_total), t3(nr_total)

        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        array_out = 0.d0     
        do l = 1, dimen
            beta = 1; beta(l) = 2
            call sumfacto3d_spM(nc_u, nr_u, nc_v, nr_v, nc_w, nr_w, &
                        nnz_u, indi_T_u, indj_T_u, data_BT_u(:, beta(1)), & 
                        nnz_v, indi_T_v, indj_T_v, data_BT_v(:, beta(2)), & 
                        nnz_w, indi_T_w, indj_T_w, data_BT_w(:, beta(3)), & 
                        array_in, t1) 
            t1 = t1*mat%Hprop*mat%detJ
            do i = 1, dimen
                alpha = 1; zeta = beta + (alpha - 1)*2
                t2 = t1*mat%invJ(l, i, :)
                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, & 
                            nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                            nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                            nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), t2, t3)
                array_out(i, :) = array_out(i, :) + t3
            end do
        end do
            
    end subroutine mf_gradu_v_3d

    subroutine intforce_2d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v, &
                            indi_u, indj_u, indi_v, indj_v, data_W_u, data_W_v, stress, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 2
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nnz_u, nnz_v
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v)
        double precision, intent(in) :: data_W_u, data_W_v
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4)
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ----------
        double precision :: Tstress, t1, t2
        dimension :: Tstress(dimen, dimen), t1(dimen, dimen, nc_total), t2(nr_total)
        integer :: i, k, alpha(dimen), zeta(dimen)
        
        if (nr_total.ne.nr_u*nr_v) stop 'Size problem'
        
        do i = 1, nc_total
            call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
            t1(:, :, i) = matmul(mat%invJ(:, :, i), Tstress)*mat%detJ(i)
        end do
        
        array_out = 0.d0
        do i = 1, dimen
            do k = 1, dimen
                alpha = 1; alpha(k) = 2
                zeta  = 1 + (alpha - 1)*2

                call sumfacto2d_spM(nr_u, nc_u, nr_v, nc_v, &
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                t1(k, i, :), t2)
                array_out(i, :) = array_out(i, :) + t2
            end do
        end do

    end subroutine intforce_2d

    subroutine intforce_3d(mat, nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w, &
                            indi_u, indj_u, indi_v, indj_v, indi_w, indj_w, data_W_u, data_W_v, data_W_w, stress, array_out)
        !! Computes internal force vector in 3D 
        !! IN CSR FORMAT

        implicit none 
        ! Input / output data
        ! -------------------
        integer, parameter :: dimen = 3
        type(mecamat) :: mat
        integer, intent(in) :: nr_total, nc_total, nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, nnz_u, nnz_v, nnz_w
        integer, intent(in) ::  indi_u, indj_u, indi_v, indj_v, indi_w, indj_w
        dimension ::    indi_u(nr_u+1), indj_u(nnz_u), &
                        indi_v(nr_v+1), indj_v(nnz_v), &
                        indi_w(nr_w+1), indj_w(nnz_w)
        double precision, intent(in) :: data_W_u, data_W_v, data_W_w
        dimension :: data_W_u(nnz_u, 4), data_W_v(nnz_v, 4), data_W_w(nnz_w, 4)
        double precision, intent(in) :: stress
        dimension :: stress(mat%nvoigt, nc_total)

        double precision, intent(out) :: array_out
        dimension :: array_out(dimen, nr_total)

        ! Local data
        ! ----------
        double precision :: Tstress, t1, t2
        dimension :: Tstress(dimen, dimen), t1(dimen, dimen, nc_total), t2(nr_total)
        integer :: i, k, alpha(dimen), zeta(dimen)
        
        if (nr_total.ne.nr_u*nr_v*nr_w) stop 'Size problem'
        
        do i = 1, nc_total
            call array2symtensor(mat%dimen, mat%nvoigt, stress(:, i), Tstress)
            t1(:, :, i) = matmul(mat%invJ(:, :, i), Tstress)*mat%detJ(i)
        end do
        
        array_out = 0.d0
        do i = 1, dimen
            do k = 1, dimen
                alpha = 1; alpha(k) = 2
                zeta  = 1 + (alpha - 1)*2

                call sumfacto3d_spM(nr_u, nc_u, nr_v, nc_v, nr_w, nc_w, &
                                nnz_u, indi_u, indj_u, data_W_u(:, zeta(1)), &
                                nnz_v, indi_v, indj_v, data_W_v(:, zeta(2)), &
                                nnz_w, indi_w, indj_w, data_W_w(:, zeta(3)), &
                                t1(k, i, :), t2)
                array_out(i, :) = array_out(i, :) + t2
            end do
        end do

    end subroutine intforce_3d

end module matrixfreeplasticity