module mfr

    use trans
    
    implicit none


    ! mfr wall temperature variables
    real(8) :: Nu
    real(8) :: d_mm, d_m
    real(8) :: Tw_min, Tw_max
    real(8) :: alpha, beta 


    ! namelist wall temperature variables
    namelist /mfrwalltemp/ Nu, d_mm, Tw_min, Tw_max, alpha, beta

contains


    subroutine set_wallT_boundary(wallT_bound_file)

        implicit none


        character(len=128),intent(in) :: wallT_bound_file
        integer :: lun


        ! read mfr wall temperature params
        open(newunit=lun, file=trim(wallT_bound_file), status="old", action="read")
        read(unit=lun, nml=mfrwalltemp)
        close(unit=lun)

        d_m = d_mm * 1.0d-3 ! mm -> m


    end subroutine set_wallT_boundary


    subroutine heat_transfer(rho, p, T, nspe, x_sp, x_m, Qw)

        implicit none


        integer,intent(in)  :: nspe
        real(8),intent(in)  :: rho, p, T
        real(8),intent(in)  :: x_sp(nspe)
        real(8),intent(in)  :: x_m
        real(8),intent(out) :: Qw

        real(8) :: lambda 


        lambda = thermal_conductivity(rho, p, T, nspe, x_sp)
        Qw = - 4.d0*lambda*Nu*(T-Tw(x_m)) / (d_m**2)


    end subroutine heat_transfer


    function Tw(x)

        implicit none


        real(8),intent(in) :: x

        real(8) :: Tw
        real(8) :: A, B


        A = 0.5d0*(Tw_max - Tw_min)
        B = 0.5d0*(Tw_max + Tw_min)
        Tw = A * dtanh(alpha*(x-beta)) + B


    end function Tw


end module mfr


