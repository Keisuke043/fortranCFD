module source

    use thermo
    use ODEs
    use radiation
    use mfr
    use ignitionsource
    use module_boundary
    use eos
    use module_params, only:eos_st
    use parallel, only: myrank, nprocs

    implicit none

contains

    subroutine reaction_macks(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, lvir, dt, &
                              nspe, nequ, block_all_0, block_all_n)
                                

        implicit none


        integer :: i, ib, k_sp

        real(8) :: rho, temp, y_sp(nspe)
        real(8) :: yiT(nspe+1)

        real(8),intent(in) :: dt
        integer,intent(in) :: nspe, nequ, lvir

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)

        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)


        do ib = bsta, bend

            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                eos_obj(i)%y_sp0(:) = b_obj(ib)%qf(nprim+1:nequ,i)


                rho  = b_obj(ib)%qf(iqdens,i)
                temp = b_obj(ib)%qf(iqtemp,i)
                do k_sp = 1, nspe
                    y_sp(k_sp) = b_obj(ib)%qf(nprim+k_sp,i)
                end do

                call adjust_sum(nspe, y_sp)

                yiT(1:nspe) = y_sp
                yiT(nspe+1) = temp

                call macks(dt, nspe+1, yiT, rho)

                y_sp = yiT(1:nspe)
                temp = yiT(nspe+1)

                call adjust_sum(nspe, y_sp)

                b_obj(ib)%qf(iqtemp,i) = temp
                do k_sp = 1, nspe
                    b_obj(ib)%qf(nprim+k_sp,i) = y_sp(k_sp)
                end do

                eos_obj(i)%HRR = heat_release_rate(nspe, eos_obj(i)%y_sp0(:), b_obj(ib)%qf(nprim+1:nequ,i), &
                                                                               b_obj(ib)%qf(iqdens,i), dt)
                eos_obj(i)%HRR_dot = heat_release_rate_omegadot(nspe, b_obj(ib)%qf(nprim+1:nequ,i), &
                                                   b_obj(ib)%qf(iqdens,i), b_obj(ib)%qf(iqtemp,i))
            end do

        end do

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                b_obj(ib)%Qc(:,i) = qfToQc(b_obj(ib)%qf(:,i), nspe, nequ)
            end do
        end do

    end subroutine reaction_macks



    subroutine heat_transfer_radiation(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, lvir, dt, &
                                       nspe, nequ, block_all_0, block_all_n)

        implicit none
        
        
        integer :: i, ib
        real(8) :: rho, temp, pres
        real(8) :: x_sp(nspe), y_sp(nspe)
        real(8) :: Qrad
        
        integer,intent(in) :: nspe, nequ, lvir
        real(8),intent(in) :: dt
        
        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)
        
        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)
        
        
        do ib = bsta, bend
        
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
        
                Qrad = 0.d0
        
                rho  = b_obj(ib)%qf(iqdens,i)
                temp = b_obj(ib)%qf(iqtemp,i)
                y_sp = b_obj(ib)%qf(nprim+1:nequ,i)
                call Yi2Xi(nspe, y_sp, x_sp)
                pres = pressure(rho, temp, nspe, x_sp)
                call radiation_heat_loss(pres, temp, nspe, x_sp, Qrad)
        
                eos_obj(i)%Qrad = Qrad
        
                b_obj(ib)%Qc(iuener,i) = b_obj(ib)%Qc(iuener,i) + dt*Qrad
        
            end do
        
        end do
        
        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do
    
    end subroutine heat_transfer_radiation


    subroutine heat_transfer_mfr(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, lvir, dt, &
                                 nspe, nequ, block_all_0, block_all_n)

        implicit none


        integer :: i, ib
        real(8) :: rho, temp, pres
        real(8) :: x_sp(nspe), y_sp(nspe)
        real(8) :: Qw
        
        integer,intent(in) :: nspe, nequ, lvir
        real(8),intent(in) :: dt

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)

        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)


        do ib = bsta, bend

            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                Qw = 0.d0

                rho  = b_obj(ib)%qf(iqdens,i)
                temp = b_obj(ib)%qf(iqtemp,i)
                y_sp = b_obj(ib)%qf(nprim+1:nequ,i)
                call Yi2Xi(nspe, y_sp, x_sp)
                pres = pressure(rho, temp, nspe, x_sp)
                call heat_transfer(rho, pres, temp, nspe, x_sp, b_obj(ib)%xm_b(i), Qw)

                eos_obj(i)%qTw = Qw

                b_obj(ib)%Qc(iuener,i) = b_obj(ib)%Qc(iuener,i) + dt*Qw

            end do

        end do

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do

    end subroutine heat_transfer_mfr


    subroutine energy_source(ista_rank, iend_rank, bsta, bend, b_obj, eos_obj, lvir, dt, &
                             nspe, nequ, block_all_0, block_all_n, time, flag_igsource, igsource_type)

        implicit none


        integer :: i, ib
        real(8) :: time, xm, Eadd

        integer,intent(in) :: flag_igsource, igsource_type
        integer,intent(in) :: nspe, nequ, lvir
        real(8),intent(in) :: dt

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)

        integer,intent(in) :: ista_rank, iend_rank
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)


        do ib = bsta, bend

            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                Eadd = 0.d0

                xm = b_obj(ib)%xm_b(i)

                if (igsource_type == 1) then
                    call add_igsource(time, xm, Eadd, flag_igsource)
                else if (igsource_type == 2) then
                    call add_igsource_WZhang(time, xm, Eadd, flag_igsource)
                end if

                b_obj(ib)%Qc(iuener,i) = b_obj(ib)%Qc(iuener,i) + dt*Eadd

            end do

        end do

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do

    end subroutine energy_source


end module source


