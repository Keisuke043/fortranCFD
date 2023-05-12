program main


    use thermo
    use source
    use update
    use parallel
    use hdf5_mpi
    use block_t
    use read_hdf5
    use module_params 
    use module_boundary
    use radiation

    implicit none


    real(8) :: time = 0.d0   ! time [s]
    integer :: nt   = 0      ! number of time step
    integer :: nout = 0      ! number of output
    integer :: nt_rstr   = 0 ! number of restart time step
    integer :: nout_rstr = 0 ! number of restart output
    real(8) :: time_rstr = 0.d0
    real(8) :: time_rstr_out = 0.d0
    real(8) :: t1, t2
    real(8) :: comp_time = 0.d0
    real(8) :: comp_time_rstr
    real(8) :: comp_time_rstr_out

    call init_solvers()
    call set_thermo(chem_dir)
    call set_trans(chem_dir)
    call main_sim()

contains


    subroutine main_sim

        implicit none

        integer :: i, ib
        integer :: ilvir
        
        call cpu_time( t1 )

        call init_parallel()

        call read_inputfiles(myrank)

        if ( flag_lbound == 1 .or. flag_rbound == 1 ) then
            call set_flow_boundary(flow_bound_file, nspe)
        end if

        if ( flag_radiation /= 0) then
            call set_radiation_params(nspe, radiation_file)
        end if

        if ( flag_Tw == 1 ) then
            call set_wallT_boundary(wallT_bound_file)
        end if

        if ( flag_igsource /= 0) then
            call set_igsource(igsource_file, is_spherical)
        end if

        if (flag_block == 0) then

            call set_one_block_size(dx_comp, nx_comp, eos_obj)

        else if (flag_block == 1) then

            call set_block_size(dx_comp, nx_comp, ncell_per_block, eos_obj, &
                                lev_smr_0, num_smr_0, lev_smr_n, num_smr_n)

        end if

        print*, myrank,    nx_all_0, '/', ista_rank, iend_rank,'/', nx_all_n,    'cell' 
        print*, myrank, block_all_0, '/',      bsta,      bend,'/', block_all_n, 'block'

        call initialize_variables(myrank, ista_rank, iend_rank)

        if (flag_restart == 1) then
            call read_hdf5_restart(ista_rank, iend_rank, nspe, eos_obj, &
                                   nx_all_0, nx_all_n, readpath, &
                                   nout_rstr, nt_rstr, time_rstr, comp_time_rstr)
            nout = nout + nout_rstr
            nt   = nt   + nt_rstr
            time = time + time_rstr
        end if

        call set_block_params(ista_rank, iend_rank, nspe, nequ, lvir, eos_obj)

        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)

        call interpolation_block(nequ, lvir)

        if (nprocs /= 1) call communicate_parallel_block(nequ, lvir)

        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)

        do ib = bsta, bend
            do ilvir = 1, lvir
                i = b_obj(ib)%ista_b-ilvir
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)

                i = b_obj(ib)%iend_b+ilvir
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do

        call init_output(savedir, nout)

        if (flag_restart == 0) then
            call write_hdf5_b(nout, nt, comp_time, comp_time_rstr, dt, &
                              time, time_rstr, ista_rank, iend_rank, lvir, &
                              nspe, eos_obj, nx_all_0, nx_all_n, chem_dir)
        end if

        if (flag_check_init == 1) then
            call fin_output
            call fin_parallel
            stop
        end if

        do while (time < tend)

            if (flag_cfl == 1) then

                call calc_dt(bsta, bend, b_obj, dx_comp, nspe, nequ, cfl, dt)

            end if

            nt = nt + 1
            time = time+dt


            call solve_fluid(ista_rank, iend_rank, bsta, bend, b_obj, &
                             lvir, dt, nspe, nequ, block_all_0, block_all_n)

            call solve_source(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, &
                              lvir, dt, nspe, nequ, block_all_0, block_all_n, time)


            call update_boundary(bsta, bend, b_obj, lvir, &
                                 nspe, nequ, block_all_0, block_all_n)

            call interpolation_block(nequ, lvir)

            if (nprocs /= 1) call communicate_parallel_block(nequ, lvir)

            call update_boundary(bsta, bend, b_obj, lvir, &
                                 nspe, nequ, block_all_0, block_all_n)

            do ib = bsta, bend
                do i = b_obj(ib)%ista_b-lvir, b_obj(ib)%iend_b+lvir
                    b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
                end do
            end do


            if (int(time/out_timestep) == nout) then
                
                call set_eos_params(ista_rank, iend_rank, nspe, nequ, &
                                    lvir, bsta, bend, b_obj, eos_obj)
                call cpu_time( t2 )

                comp_time = t2-t1+comp_time_rstr
                comp_time_rstr_out = t2-t1
                time_rstr_out = time - time_rstr

                call write_hdf5_b(nout, nt, comp_time_rstr_out, comp_time, dt, &
                                  time_rstr_out, time, ista_rank, iend_rank, lvir, &
                                  nspe, eos_obj, nx_all_0, nx_all_n, chem_dir)

                if ( mod(nout, out_num_hdf) == 0 ) then

                    call fin_output

                    ! if ( time < tend - 2*dt ) then

                    if ( time < tend) then
                        call init_output(savedir, nout)
                    end if
                end if
            end if

            if (myrank == 0) write(*,*) nout, nt, 'dt= ', dt, 't = ', time-time_rstr, time

        end do

        call fin_output
        call fin_parallel

    end subroutine main_sim

end program


