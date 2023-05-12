module flow_boundary

    use thermo

    implicit none


    real(8) :: T_IN, u_cm_in, u_in, rho_u_IN
    real(8),allocatable :: x_sp_IN(:), y_sp_IN(:)
    real(8) :: p_atm_OUT, p_OUT


contains

    subroutine set_flow_boundary(flow_bound_file, nspe)

        implicit none


        integer,intent(in) :: nspe
        character(len=128),intent(in) :: flow_bound_file

        integer :: k_sp
        integer :: nspe_in
        integer :: k_sp_in
        real(8) :: sp_mole_in
        character(len=128) :: sp_name_in


        allocate(x_sp_IN(nspe), y_sp_IN(nspe))
        x_sp_IN(:) = 0.d0
        
        open(19,file=trim(flow_bound_file), status='old')

        read(19,*) T_IN
        read(19,*) u_cm_in

        read(19,*) nspe_in
        do k_sp_in = 1, nspe_in
            read(19,*) sp_name_in, sp_mole_in
            do k_sp = 1, nspe
                if (trim(sp_name_in) == species_name(k_sp)) then
                    x_sp_IN(k_sp) = sp_mole_in
                end if
            end do
        end do

        read(19,*) p_atm_OUT

        u_in  = u_cm_in  *1.d-2    ! cm/s -> m/s
        p_OUT = p_atm_OUT*101325d0 ! atm  -> Pa

        call Xi2Yi(nspe, x_sp_IN, y_sp_IN)
        rho_u_IN = density(p_OUT, T_IN, nspe, x_sp_IN) * u_in

        close(19)

    end subroutine set_flow_boundary

end module flow_boundary




