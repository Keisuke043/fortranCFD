module update

    use block_t, only:block_st
    use module_params, only:eos_st
    use timeintegration
    use source

    implicit none

contains


    subroutine calc_dt(bsta, bend, b_obj, dx, nspe, nequ, cfl, dt)

        implicit none


        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: dx, cfl
        real(8),intent(inout) :: dt

        integer :: i, ib
        real(8) :: cvel, aupc
        real(8) :: rho, temp, pres, a_sound, x_sp(nspe)
        real(8) :: dt_tmp

        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        cvel = 0.d0

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                call Yi2Xi(nspe, b_obj(ib)%qf(nprim+1:nequ,i), x_sp)

                rho     = b_obj(ib)%qf(iqdens,i)
                temp    = b_obj(ib)%qf(iqtemp,i)
                pres    = pressure(rho, temp, nspe, x_sp)
                a_sound = sonic_speed(rho, pres, temp, nspe, x_sp)
                aupc    = dmax1(dabs(b_obj(ib)%qf(iqxvel,i) - a_sound), &
                                dabs(b_obj(ib)%qf(iqxvel,i)), &
                                dabs(b_obj(ib)%qf(iqxvel,i) + a_sound))

                if (aupc > cvel) then
                    cvel = aupc
                end if

            end do
        end do

        dt = (cfl*dx)/cvel

        call MPI_ALLREDUCE(dt, dt_tmp, 1 ,MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
        dt = dt_tmp


    end subroutine calc_dt


    subroutine solve_fluid(ista_rank, iend_rank, bsta, bend, b_obj, &
                           lvir, dt, nspe, nequ, block_all_0, block_all_n)

        use module_params, only: select_time

        implicit none


        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: dt

        integer,intent(in) :: ista_rank, iend_rank, lvir
        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        if (select_time == 0) then

            call explicit_euler(bsta, bend, b_obj, lvir, dt, nspe, nequ, &
                                block_all_0, block_all_n)

        else if (select_time == 1) then

            call TVDRK3(ista_rank, iend_rank, bsta, bend, b_obj, &
                        lvir, dt, nspe, nequ, block_all_0, block_all_n)

        else if (select_time == 2) then

            call RK4(ista_rank, iend_rank, bsta, bend, b_obj, &
                     lvir, dt, nspe, nequ, block_all_0, block_all_n)

        end if


    end subroutine solve_fluid


    subroutine solve_source(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                            lvir, dt, nspe, nequ, block_all_0, block_all_n, time)

        use module_params, only: do_reaction, flag_Tw, flag_igsource, flag_radiation

        implicit none


        integer,intent(in) :: nspe, nequ, lvir
        real(8),intent(in) :: dt, time

        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        if (do_reaction == 1) then
        
            call reaction_macks(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                                lvir, dt, nspe, nequ, block_all_0, block_all_n)
        
        end if


        if (flag_radiation == 1) then
        
            call heat_transfer_radiation(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                                   lvir, dt, nspe, nequ, block_all_0, block_all_n)

        end if


        if (flag_Tw == 1) then
        
            call heat_transfer_mfr(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                                   lvir, dt, nspe, nequ, block_all_0, block_all_n)
        
        end if


        if (flag_igsource /= 0) then
        
            call energy_source(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                               lvir, dt, nspe, nequ, block_all_0, block_all_n, time, flag_igsource)
        
        end if


    end subroutine solve_source


end module update

