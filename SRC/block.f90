
module block_t

    use variables


    implicit none


    type :: block_st

        integer :: smr_level
        integer :: ista_b, iend_b
        real(8) :: smr_dx

        real(8),allocatable :: qf(:,:)    ! Primitive
        real(8),allocatable :: Qc(:,:)    ! Conservative
        real(8),allocatable :: Fplus(:,:) ! Numerical flux
        real(8),allocatable :: xm_b(:)    ! Location

        real(8),allocatable :: godunov_state(:,:)

    end type block_st

    real(8),allocatable :: dx_block(:)

    integer :: ista_rank, iend_rank
    integer :: nx_all_0, nx_all_n, nx_all

    integer :: bsta, bend
    integer :: block_all_0, block_all_n

    type(block_st),allocatable :: b_obj(:)

contains


    subroutine set_one_block_size(dx_comp, nx_comp, eos_obj)

        use parallel
        use module_params, only:eos_st

        implicit none


        integer :: i
        integer :: ista_b, iend_b

        real(8), intent(in) :: dx_comp ! space size [m]
        integer, intent(in) :: nx_comp ! number of space step
        type(eos_st),allocatable,intent(inout) :: eos_obj(:)


        nx_all_0 = 1       ! cell center min
        nx_all_n = nx_comp ! cell center max
        nx_all   = nx_all_n - nx_all_0 + 1
        
        call set_parallel(ista_b, iend_b, nx_all_0, nx_all_n)

        ! print *, 'myrank = ', myrank, 'ista_b = ', ista_b, 'iend_b = ', iend_b

        ista_rank = ista_b
        iend_rank = iend_b

        bsta = 1
        bend = 1
        call set_block_obj(bsta, bend)

        b_obj(1)%smr_level = 1
        b_obj(1)%smr_dx = dx_comp
        b_obj(1)%ista_b = ista_b
        b_obj(1)%iend_b = iend_b

        allocate( b_obj(1)%xm_b(b_obj(1)%ista_b:b_obj(1)%iend_b) )

        do i = b_obj(1)%ista_b, b_obj(1)%iend_b
            b_obj(1)%xm_b(i) = b_obj(1)%smr_dx * i - b_obj(1)%smr_dx/2.d0
        end do

        allocate( eos_obj(ista_rank:iend_rank) )

        do i = ista_rank, iend_rank
            eos_obj(i)%x_m = b_obj(1)%xm_b(i)
            eos_obj(i)%dx  = b_obj(1)%smr_dx
        end do


    end subroutine set_one_block_size


    subroutine set_block_size(dx_comp, nx_comp, ncell_per_block, eos_obj, &
                              lev_smr_0, num_smr_0, lev_smr_n, num_smr_n)
        
        use parallel
        use module_params, only:eos_st

        implicit none

        real(8), intent(in) :: dx_comp ! space size [m]
        integer, intent(in) :: nx_comp ! number of space step
        integer, intent(in) :: ncell_per_block
        integer, intent(in) :: lev_smr_0, num_smr_0, lev_smr_n, num_smr_n
        type(eos_st),allocatable,intent(inout) :: eos_obj(:)

        integer :: n_rank      ! counter in nprocs
        integer :: count_block ! counts for send dx
        integer :: tmp_displs_recv
        integer,allocatable :: count_recv(:), displs_recv(:) ! counts for recv dx
        real(8),allocatable :: tmp_dx_block(:)

        integer :: block_num_comp
        integer :: block_num_all

        integer :: i, ib, ib_smrlev, ib_smrnum
        integer :: tmp_ib_smr, tmp_ib_smr_n
        real(8) :: length_smr_0, length_smr_n
        real(8) :: length_smr, length_smr_tmp
        real(8) :: xm_b_lefti, dx_lefti

        real(8) :: IS_left_xm, IR_left_xm
        real(8) :: IS_left_dx, IR_left_dx


        nx_all_0 = 1       - int(ncell_per_block*(lev_smr_0-1)*num_smr_0) ! cell center min
        nx_all_n = nx_comp + int(ncell_per_block*(lev_smr_n-1)*num_smr_n) ! cell center max
        nx_all   = nx_all_n - nx_all_0 + 1

        if ( myrank == 0 ) then
            if ( mod(nx_all, ncell_per_block) /= 0 ) then
                print *, 'Num. of cell should be a multiple of 10'
                call fin_parallel
                stop
            else if ( mod(nx_comp, ncell_per_block) /= 0 ) then
                print *, 'Num. of cell should be a multiple of 10'
                call fin_parallel
                stop
            end if
        end if

        block_num_comp = int(nx_comp/ncell_per_block)
        block_num_all  = block_num_comp + int((lev_smr_0-1)*num_smr_0) &
                                        + int((lev_smr_n-1)*num_smr_n)


        block_all_0 = 1
        block_all_n = block_num_all


        call set_parallel_block(bsta, bend, block_all_0, block_all_n)

        ! print *, 'myrank = ', myrank, 'bsta = ', bsta, 'bend = ', bend


        call set_block_obj(bsta, bend)

        allocate( tmp_dx_block(bsta:bend) )
        allocate( dx_block(block_all_0:block_all_n) )

        count_block = bend - bsta + 1

        do ib = bsta, bend

            b_obj(ib)%smr_level = 1

            b_obj(ib)%smr_dx = dx_comp
            tmp_dx_block(ib) = dx_comp

        end do


        tmp_ib_smr = 0

        length_smr_0 = 0.d0
        length_smr_n = 0.d0

        do ib_smrlev = 2, lev_smr_0

            do ib_smrnum = 1, num_smr_0

                tmp_ib_smr = tmp_ib_smr + 1

                do ib = bsta, bend

                    if (ib == tmp_ib_smr) then

                        b_obj(ib)%smr_level = lev_smr_0-ib_smrlev+2
                        b_obj(ib)%smr_dx    = dx_comp * 2**(lev_smr_0-ib_smrlev+1)
                        tmp_dx_block(ib)    = b_obj(ib)%smr_dx

                        length_smr_0 = length_smr_0 + b_obj(ib)%smr_dx*ncell_per_block

                        ! print *, myrank, ib ,tmp_ib_smr
                        ! print *, length_smr_0

                    end if

                end do

            end do

        end do

        tmp_ib_smr = 0
        tmp_ib_smr_n = 0

        do ib_smrlev = 2, lev_smr_n

            do ib_smrnum = 1, num_smr_n

                tmp_ib_smr = tmp_ib_smr + 1
                tmp_ib_smr_n = block_all_n-tmp_ib_smr + 1

                ! print *, block_all_n, tmp_ib_smr, tmp_ib_smr_n

                do ib = bsta, bend

                    if (ib == tmp_ib_smr_n) then

                        b_obj(ib)%smr_level = lev_smr_n-ib_smrlev+2
                        b_obj(ib)%smr_dx    = dx_comp * 2**(lev_smr_n-ib_smrlev+1)
                        tmp_dx_block(ib)    = b_obj(ib)%smr_dx

                        length_smr_n = length_smr_n + b_obj(ib)%smr_dx*ncell_per_block

                        ! print *, ib ,tmp_ib_smr_n
                        ! print *, length_smr_n

                    end if

                end do

            end do

        end do

        call MPI_ALLREDUCE(length_smr_0, length_smr_tmp, 1 ,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        length_smr = length_smr_tmp


        allocate( count_recv(0:nprocs-1) )
        allocate( displs_recv(0:nprocs-1) )

        call MPI_ALLGATHER(count_block, 1, MPI_INTEGER, &
                           count_recv,  1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        tmp_displs_recv = 0
        do n_rank = 0, nprocs-1
            displs_recv(n_rank) = tmp_displs_recv
            tmp_displs_recv = tmp_displs_recv + count_recv(n_rank)
        end do

        call MPI_ALLGATHERV(tmp_dx_block, count_block,             MPI_REAL8, &
                                dx_block, count_recv, displs_recv, MPI_REAL8, MPI_COMM_WORLD, ierr)


        do ib = bsta, bend

            b_obj(ib)%ista_b = nx_all_0 + ncell_per_block*(ib-1)
            b_obj(ib)%iend_b = b_obj(ib)%ista_b + ncell_per_block - 1

            allocate( b_obj(ib)%xm_b(b_obj(ib)%ista_b:b_obj(ib)%iend_b) )


            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                if (ib == bsta .and. i == b_obj(ib)%ista_b) then

                    b_obj(ib)%xm_b(i) = 0.d0
                    xm_b_lefti = 0.d0
                    dx_lefti   = 0.d0

                else

                    b_obj(ib)%xm_b(i) = xm_b_lefti + b_obj(ib)%smr_dx
                    xm_b_lefti        = xm_b_lefti + b_obj(ib)%smr_dx

                    if (dx_lefti > b_obj(ib)%smr_dx) then
                        b_obj(ib)%xm_b(i) = b_obj(ib)%xm_b(i) + b_obj(ib)%smr_dx/2.d0
                        xm_b_lefti        = xm_b_lefti        + b_obj(ib)%smr_dx/2.d0
                    end if

                    if (dx_lefti < b_obj(ib)%smr_dx) then
                        b_obj(ib)%xm_b(i) = b_obj(ib)%xm_b(i) - b_obj(ib)%smr_dx/4.d0
                        xm_b_lefti        = xm_b_lefti        - b_obj(ib)%smr_dx/4.d0
                    end if

                end if
                dx_lefti = b_obj(ib)%smr_dx

            end do

        end do


        if ( myrank == 0 ) then

            length_smr = length_smr - b_obj(bsta)%smr_dx/2.d0

            do ib = bsta, bend
                do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                    b_obj(ib)%xm_b(i) = b_obj(ib)%xm_b(i) - length_smr

                end do
            end do

            IS_left_xm = b_obj(bend)%xm_b(b_obj(bend)%iend_b)
            IS_left_dx = b_obj(bend)%smr_dx

            call MPI_SEND(IS_left_xm, 1, MPI_REAL8, IUP, 1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(IS_left_dx, 1, MPI_REAL8, IUP, 1, MPI_COMM_WORLD, ierr)

        else if ( 1 <= myrank .and. myrank <= nprocs-2 ) then

            call MPI_RECV(IR_left_xm, 1, MPI_REAL8, IDOWN, 1, MPI_COMM_WORLD, ISTATUS, ierr)
            call MPI_RECV(IR_left_dx, 1, MPI_REAL8, IDOWN, 1, MPI_COMM_WORLD, ISTATUS, ierr)

            xm_b_lefti = IR_left_xm + (IR_left_dx + b_obj(bsta)%smr_dx)/2.d0

            do ib = bsta, bend
                do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                    b_obj(ib)%xm_b(i) = xm_b_lefti + b_obj(ib)%xm_b(i)

                end do
            end do

            IS_left_xm = b_obj(bend)%xm_b(b_obj(bend)%iend_b)
            IS_left_dx = b_obj(bend)%smr_dx

            call MPI_SEND(IS_left_xm, 1, MPI_REAL8, IUP, 1, MPI_COMM_WORLD, ierr)
            call MPI_SEND(IS_left_dx, 1, MPI_REAL8, IUP, 1, MPI_COMM_WORLD, ierr)

        else if ( myrank == nprocs-1 ) then

            call MPI_RECV(IR_left_xm, 1, MPI_REAL8, IDOWN, 1, MPI_COMM_WORLD, ISTATUS, ierr)
            call MPI_RECV(IR_left_dx, 1, MPI_REAL8, IDOWN, 1, MPI_COMM_WORLD, ISTATUS, ierr)

            xm_b_lefti = IR_left_xm + (IR_left_dx + b_obj(bsta)%smr_dx)/2.d0

            do ib = bsta, bend
                do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                    b_obj(ib)%xm_b(i) = xm_b_lefti + b_obj(ib)%xm_b(i)

                end do
            end do

        end if

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                b_obj(ib)%xm_b(i) = b_obj(ib)%xm_b(i)

            end do
        end do


        ! print *, myrank, 'length_smr', length_smr


        ista_rank = b_obj(bsta)%ista_b
        iend_rank = b_obj(bend)%iend_b

        ! print *, myrank, ista_rank, iend_rank

        allocate( eos_obj(ista_rank:iend_rank) )

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                ! print *, myrank, ib, i, b_obj(ib)%smr_dx, b_obj(ib)%xm_b(i)
                eos_obj(i)%x_m = b_obj(ib)%xm_b(i)
                eos_obj(i)%dx  = b_obj(ib)%smr_dx
            end do
        end do


    end subroutine set_block_size


    subroutine set_block_obj(bsta, bend)

        implicit none


        integer,intent(in) :: bsta, bend

        allocate( b_obj(bsta:bend) )

    end subroutine set_block_obj


    subroutine set_block_params(ista_rank, iend_rank, nspe, nequ, lvir, eos_obj)

        use module_params, only:eos_st

        implicit none


        integer :: ib, i
        integer,intent(in) :: nspe, nequ, lvir
        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(in) :: eos_obj(ista_rank:iend_rank)


        do ib = bsta, bend

            allocate( b_obj(ib)%Qc(nequ,b_obj(ib)%ista_b-lvir:b_obj(ib)%iend_b+lvir) )
            allocate( b_obj(ib)%qf(nequ,b_obj(ib)%ista_b-lvir:b_obj(ib)%iend_b+lvir) )
            allocate( b_obj(ib)%Fplus(nequ,b_obj(ib)%ista_b-1:b_obj(ib)%iend_b) )

            allocate( b_obj(ib)%godunov_state(nprim,b_obj(ib)%ista_b-1:b_obj(ib)%iend_b) )

            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                b_obj(ib)%qf(iqdens,i) = eos_obj(i)%rho
                b_obj(ib)%qf(iqxvel,i) = eos_obj(i)%u
                b_obj(ib)%qf(iqtemp,i) = eos_obj(i)%T
                b_obj(ib)%qf(nprim+1:nequ,i) = eos_obj(i)%y_sp(:)

                b_obj(ib)%Qc(iudens,i) = eos_obj(i)%rho
                b_obj(ib)%Qc(iumomx,i) = eos_obj(i)%rho*eos_obj(i)%u
                b_obj(ib)%Qc(iuener,i) = eos_obj(i)%rho*eos_obj(i)%eTotal
                b_obj(ib)%Qc(ncons+1:nequ,i) = eos_obj(i)%rho*eos_obj(i)%y_sp(:)

            end do

        end do

    end subroutine set_block_params


    subroutine set_eos_params(ista_rank, iend_rank, nspe, nequ, lvir, &
                              bsta, bend, b_obj, eos_obj)

        use module_params, only:eos_st
        use thermo, only:Yi2Xi, pressure, specific_energy

        implicit none


        integer :: ib, i
        integer,intent(in) :: nspe, nequ, lvir
        integer,intent(in) :: ista_rank, iend_rank
        integer,intent(in) :: bsta, bend
        type(block_st),intent(in)  :: b_obj(bsta:bend)
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)


        do ib = bsta, bend

            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                eos_obj(i)%rho     = b_obj(ib)%qf(iqdens,i)
                eos_obj(i)%u       = b_obj(ib)%qf(iqxvel,i)
                eos_obj(i)%T       = b_obj(ib)%qf(iqtemp,i)
                eos_obj(i)%y_sp(:) = b_obj(ib)%qf(nprim+1:nequ,i)
                call Yi2Xi(nspe, eos_obj(i)%y_sp(:), eos_obj(i)%x_sp(:))
                eos_obj(i)%p       = pressure(eos_obj(i)%rho, &
                                     eos_obj(i)%T, nspe, eos_obj(i)%x_sp(:))

                eos_obj(i)%eKinetic  = 0.5d0*eos_obj(i)%u**2
                eos_obj(i)%eInternal = specific_energy(eos_obj(i)%rho, eos_obj(i)%p, &
                                       eos_obj(i)%T, nspe, eos_obj(i)%x_sp(:))
                eos_obj(i)%eTotal    = b_obj(ib)%Qc(iuener,i)/eos_obj(i)%rho

            end do

        end do

    end subroutine set_eos_params


    subroutine interpolation_block(nequ, lvir)

        implicit none


        integer,intent(in) :: nequ, lvir

        integer :: ib, ilvir, i_c, n
        real(8) :: Qc_ave_ib1(nequ), Qc_ave_ib2(nequ)
        real(8) :: dQc_ep_m(nequ), dQc_ep_p(nequ)
        real(8) :: dQc_ep(nequ)


        do ib = bsta, bend-1

            if ( b_obj(ib)%smr_dx == b_obj(ib+1)%smr_dx ) then

                do ilvir = 1, lvir

                    b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b-ilvir) = b_obj(ib  )%Qc(:,b_obj(ib  )%iend_b-ilvir+1)
                    b_obj(ib  )%Qc(:,b_obj(ib  )%iend_b+ilvir) = b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+ilvir-1)

                end do

            else if ( b_obj(ib)%smr_dx > b_obj(ib+1)%smr_dx ) then

                do ilvir = 1, lvir
                    b_obj(ib)%Qc(:,b_obj(ib)%iend_b+ilvir) = &
                                   ( b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+int(2*(ilvir-1)  )) &
                                   + b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+int(2*(ilvir-1)+1)) ) / 2.d0
                end do

                do ilvir = 1, lvir

                    i_c = -int((ilvir-1) / 2)

                    dQc_ep_m = b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c  ) &
                             - b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c-1)
                    dQc_ep_p = b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c+1) &
                             - b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c  )

                    do n = 1, nequ
                        dQc_ep(n) = 0.5*(dsign(1.d0, dQc_ep_m(n))+dsign(1.d0, dQc_ep_p(n))) &
                                       *(dmin1(dabs(dQc_ep_m(n)), dabs(dQc_ep_p(n))))
                    end do

                    if( mod(ilvir,2) == 1 ) then
                        b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b-ilvir) = b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c) &
                                                                     + 0.25d0*dQc_ep
                    else if( mod(ilvir,2) == 0 ) then
                        b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b-ilvir) = b_obj(ib)%Qc(:,b_obj(ib)%iend_b+i_c) &
                                                                     - 0.25d0*dQc_ep
                    end if

                end do


            else if ( b_obj(ib)%smr_dx < b_obj(ib+1)%smr_dx ) then

                do ilvir = 1, lvir

                    b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b-ilvir) = &
                                   ( b_obj(ib)%Qc(:,b_obj(ib)%iend_b-int(2*(ilvir-1)  )  ) &
                                   + b_obj(ib)%Qc(:,b_obj(ib)%iend_b-int(2*(ilvir-1)+1))) / 2.d0
                end do

                do ilvir = 1, lvir

                    i_c = int((ilvir-1) / 2)

                    dQc_ep_m = b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c  ) &
                             - b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c-1)
                    dQc_ep_p = b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c+1) &
                             - b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c  )

                    do n = 1, nequ
                        dQc_ep(n) = 0.5*(dsign(1.d0, dQc_ep_m(n))+dsign(1.d0, dQc_ep_p(n))) &
                                       *(dmin1(dabs(dQc_ep_m(n)), dabs(dQc_ep_p(n))))
                    end do

                    if( mod(ilvir,2) == 1 ) then
                        b_obj(ib)%Qc(:,b_obj(ib)%iend_b+ilvir) = b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c) &
                                                                 - 0.25d0*dQc_ep
                    else if( mod(ilvir,2) == 0 ) then
                        b_obj(ib)%Qc(:,b_obj(ib)%iend_b+ilvir) = b_obj(ib+1)%Qc(:,b_obj(ib+1)%ista_b+i_c) &
                                                                 + 0.25d0*dQc_ep
                    end if

                end do

            end if

        end do

    end subroutine interpolation_block


    subroutine communicate_parallel_block(nequ, lvir)

        use parallel

        implicit none


        integer,intent(in) :: nequ, lvir

        integer :: ilvir, i_c, n
        real(8) :: dQc_ep_m(nequ), dQc_ep_p(nequ)
        real(8) :: dQc_ep(nequ)
        real(8) :: Qc_minus_send(nequ, lvir), Qc_plus_send(nequ, lvir)
        real(8) :: Qc_minus_recv(nequ, lvir), Qc_plus_recv(nequ, lvir)

        integer :: ISEND_Minus, ISEND_Plus
        integer :: IRECV_Minus, IRECV_Plus


        Qc_plus_send  = 0.d0
        Qc_plus_recv  = 0.d0

        if (b_obj(bend)%smr_dx == dx_block(bend+1)) then

            do ilvir = 1, lvir

                Qc_plus_send(:,ilvir) = b_obj(bend)%Qc(:,b_obj(bend)%iend_b-ilvir+1)

            end do

        else if (b_obj(bend)%smr_dx < dx_block(bend+1)) then

            do ilvir = 1, lvir
                
                Qc_plus_send(:,ilvir) = ( b_obj(bend)%Qc(:,b_obj(bend)%iend_b-int(2*(ilvir-1)  )) &
                                        + b_obj(bend)%Qc(:,b_obj(bend)%iend_b-int(2*(ilvir-1)+1))) / 2.d0

            end do

        else if (b_obj(bend)%smr_dx > dx_block(bend+1)) then

            do ilvir = 1, lvir
            
                i_c = -int((ilvir-1) / 2)
            
                dQc_ep_m = b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c  ) &
                         - b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c-1)
                dQc_ep_p = b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c+1) &
                         - b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c  )
            
                do n = 1, nequ
                    dQc_ep(n) = 0.5*(dsign(1.d0, dQc_ep_m(n))+dsign(1.d0, dQc_ep_p(n))) &
                                   *(dmin1(dabs(dQc_ep_m(n)), dabs(dQc_ep_p(n))))
                end do
            
                if( mod(ilvir,2) == 1 ) then
                    Qc_plus_send(:,ilvir) = b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c) + 0.25d0*dQc_ep
                else if( mod(ilvir,2) == 0 ) then
                    Qc_plus_send(:,ilvir) = b_obj(bend)%Qc(:,b_obj(bend)%iend_b+i_c) - 0.25d0*dQc_ep
                end if

            end do

        end if


        Qc_minus_send = 0.d0
        Qc_minus_recv = 0.d0

        if (dx_block(bsta-1) == b_obj(bsta)%smr_dx) then

            do ilvir = 1, lvir

                Qc_minus_send(:,ilvir) = b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+ilvir-1)

            end do

        else if (dx_block(bsta-1) > b_obj(bsta)%smr_dx) then

            do ilvir = 1, lvir

                Qc_minus_send(:,ilvir) = ( b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+int(2*(ilvir-1)  )) &
                                         + b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+int(2*(ilvir-1)+1)) ) / 2.d0

            end do

        else if (dx_block(bsta-1) < b_obj(bsta)%smr_dx) then

            do ilvir = 1, lvir
            
                i_c = int((ilvir-1) / 2)
            
                dQc_ep_m = b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c  ) &
                         - b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c-1)
                dQc_ep_p = b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c+1) &
                         - b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c  )
            
                do n = 1, nequ
                    dQc_ep(n) = 0.5*(dsign(1.d0, dQc_ep_m(n))+dsign(1.d0, dQc_ep_p(n))) &
                                   *(dmin1(dabs(dQc_ep_m(n)), dabs(dQc_ep_p(n))))
                end do
            
                if( mod(ilvir,2) == 1 ) then
                    Qc_minus_send(:,ilvir) = b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c) - 0.25d0*dQc_ep
                else if( mod(ilvir,2) == 0 ) then
                    Qc_minus_send(:,ilvir) = b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+i_c) + 0.25d0*dQc_ep
                end if
            
            end do

        end if

        do ilvir = 1, lvir

            call MPI_ISEND(Qc_plus_send(:,ilvir),nequ,MPI_REAL8,IUP,  1,MPI_COMM_WORLD,ISEND_Plus,ierr)
            call MPI_IRECV(Qc_plus_recv(:,ilvir),nequ,MPI_REAL8,IDOWN,1,MPI_COMM_WORLD,IRECV_Plus,ierr)
            call MPI_WAIT(IRECV_Plus,ISTATUS,ierr)
            call MPI_WAIT(ISEND_Plus,ISTATUS,ierr)


            call MPI_ISEND(Qc_minus_send(:,ilvir),nequ,MPI_REAL8,IDOWN,1,MPI_COMM_WORLD,ISEND_Minus,ierr)
            call MPI_IRECV(Qc_minus_recv(:,ilvir),nequ,MPI_REAL8,IUP,  1,MPI_COMM_WORLD,IRECV_Minus,ierr)
            call MPI_WAIT(IRECV_Minus,ISTATUS,ierr)
            call MPI_WAIT(ISEND_Minus,ISTATUS,ierr)

        end do


        do ilvir = 1, lvir

            b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b-ilvir) = Qc_plus_recv(:,ilvir)

            b_obj(bend)%Qc(:,b_obj(bend)%iend_b+ilvir) = Qc_minus_recv(:,ilvir)

        end do

    end subroutine communicate_parallel_block


end module block_t


