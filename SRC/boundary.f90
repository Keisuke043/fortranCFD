module module_boundary
   
    use block_t, only:block_st
    use eos
    use flow_boundary

    implicit none


contains

    subroutine update_boundary(bsta, bend, b_obj, lvir, nspe, nequ, block_all_0, block_all_n)

        use parallel, only: myrank, nprocs
        use module_params, only: flag_lbound, flag_rbound

        implicit none


        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        integer,intent(in) :: lvir, nspe, nequ
        type(block_st),intent(inout) :: b_obj(bsta:bend)

        integer :: i, ilvir
        real(8) :: rho_in, p_in
        real(8) :: E_internal_in


        if ( bsta == block_all_0 ) then

            if ( flag_lbound == 0) then
                call lbound_wall(bsta, bend, b_obj, lvir)
            else if ( flag_lbound == 1 ) then
                call lbound_inlet(bsta, bend, b_obj, lvir, nspe, nequ) 
            end if

            do ilvir = 1, lvir

                i = b_obj(block_all_0)%ista_b-ilvir
                b_obj(block_all_0)%qf(:,i) = &
                QcToqf(b_obj(block_all_0)%Qc(:,i), nspe, nequ)

            end do

        end if

        if ( bend == block_all_n ) then

            if ( flag_rbound == 0 ) then
                call rbound_wall(bsta, bend, b_obj, lvir)
            else if ( flag_rbound == 1 ) then
                call rbound_outlet(bsta, bend, b_obj, lvir, nspe, nequ)
            end if

            do ilvir = 1, lvir

                i = b_obj(block_all_n)%iend_b+ilvir
                b_obj(block_all_n)%qf(:,i) = &
                QcToqf(b_obj(block_all_n)%Qc(:,i), nspe, nequ)

            end do

        end if


    end subroutine update_boundary


    subroutine lbound_inlet(bsta, bend, b_obj, lvir, nspe, nequ)

        implicit none


        integer :: ilvir, k_sp
        real(8) :: rho_in, p_in
        real(8) :: Einternal_in

        integer,intent(in) :: lvir, nspe, nequ
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ilvir = 1, lvir

            rho_in = b_obj(bsta)%Qc(iudens,b_obj(bsta)%ista_b+ilvir-1)
            p_in   = pressure(rho_in, T_IN, nspe, x_sp_IN)

            b_obj(bsta)%Qc(iudens,b_obj(bsta)%ista_b-ilvir) = rho_in
            b_obj(bsta)%Qc(iumomx,b_obj(bsta)%ista_b-ilvir) = rho_u_IN

            u_in = rho_u_IN/rho_in
            Einternal_in = specific_energy(rho_in, p_in, T_IN, nspe, x_sp_IN)
            b_obj(bsta)%Qc(iuener,b_obj(bsta)%ista_b-ilvir) = rho_in*(Einternal_in + 0.5d0*u_in**2.d0)

            do k_sp = 1, nspe
                b_obj(bsta)%Qc(ncons+k_sp,b_obj(bsta)%ista_b-ilvir) = rho_in*y_sp_IN(k_sp)
            end do

        end do


    end subroutine lbound_inlet


    subroutine rbound_outlet(bsta, bend, b_obj, lvir, nspe, nequ)

        implicit none


        integer :: ilvir, k_sp
        real(8) :: u_out, rho_out, T_out
        real(8) :: Einternal_out
        real(8) :: x_sp_out(nspe), y_sp_out(nspe)

        integer,intent(in) :: lvir, nspe, nequ
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ilvir = 1, lvir

            b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b+ilvir) = - &
            b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b-ilvir+1)

            u_out         = b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b-ilvir+1) / &
                            b_obj(bend)%Qc(iudens,b_obj(bend)%iend_b-ilvir+1)
            Einternal_out = b_obj(bend)%Qc(iuener,b_obj(bend)%iend_b-ilvir+1) / &
                            b_obj(bend)%Qc(iudens,b_obj(bend)%iend_b-ilvir+1) - 0.5d0*u_out**2.d0
            y_sp_out      = b_obj(bend)%Qc(ncons+1:nequ,b_obj(bend)%iend_b-ilvir+1) / &
                            b_obj(bend)%Qc(iudens,b_obj(bend)%iend_b-ilvir+1)
            call Yi2Xi(nspe, y_sp_out, x_sp_out)
            T_out         = temperature_con(Einternal_out, nspe, x_sp_out)
            rho_out       = density(p_OUT, T_out, nspe, x_sp_out)

            b_obj(bend)%Qc(iudens,b_obj(bend)%iend_b+ilvir) = rho_out
            b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b+ilvir) = rho_out*u_out

            Einternal_out = specific_energy(rho_out, p_OUT, T_out, nspe, x_sp_out)
            b_obj(bend)%Qc(iuener,b_obj(bend)%iend_b+ilvir) = rho_out*(Einternal_out + 0.5d0*u_out**2.d0)

            do k_sp = 1, nspe
                b_obj(bend)%Qc(ncons+k_sp,b_obj(bend)%iend_b+ilvir) = rho_out*y_sp_out(k_sp)
            end do

        end do


    end subroutine rbound_outlet


    subroutine lbound_wall(bsta, bend, b_obj, lvir)

        implicit none


        integer :: ilvir
        integer,intent(in) :: lvir
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ilvir = 1, lvir

            b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b-ilvir) = &
            b_obj(bsta)%Qc(:,b_obj(bsta)%ista_b+ilvir-1)

            b_obj(bsta)%Qc(iumomx,b_obj(bsta)%ista_b-ilvir) = - &
            b_obj(bsta)%Qc(iumomx,b_obj(bsta)%ista_b+ilvir-1)

        end do


    end subroutine lbound_wall


    subroutine rbound_wall(bsta, bend, b_obj, lvir)

        implicit none


        integer :: ilvir
        integer,intent(in) :: lvir
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ilvir = 1, lvir

            b_obj(bend)%Qc(:,b_obj(bend)%iend_b+ilvir) = &
            b_obj(bend)%Qc(:,b_obj(bend)%iend_b-ilvir+1)

            b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b+ilvir) = - &
            b_obj(bend)%Qc(iumomx,b_obj(bend)%iend_b-ilvir+1)

        end do


    end subroutine rbound_wall


end module module_boundary

