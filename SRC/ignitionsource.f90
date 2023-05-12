module ignitionsource

    implicit none


    ! energy source variables
    real(8) :: Eig   ! mJ
    real(8) :: rig   ! micro m
    real(8) :: tig   ! micro sec
    real(8) :: tig_0 ! micro sec

    real(8) :: Eig_j ! J
    real(8) :: rig_m ! m
    real(8) :: tig_s ! s
    real(8) :: tig_0_s  ! s
    real(8) :: volumeig ! m^3

    real(8),parameter :: pi = 3.14159265359d0

    integer :: alpha
    real(8) :: pstep = 6
    real(8) :: gcoef
    real(8) :: kcoef


    ! namelist energy source variables
    namelist /probEsource/ Eig, rig, tig, tig_0


contains


    subroutine set_igsource(igsource_file, is_spherical)

        implicit none


        character(len=128),intent(in) :: igsource_file
        integer,intent(in) :: is_spherical
        integer :: lun


        ! read energy source params
        open(newunit=lun, file=trim(igsource_file), status="old", action="read")
        read(unit=lun, nml=probEsource)
        close(unit=lun)

        Eig_j = Eig * 1.d-3  ! mJ      -> J
        rig_m = rig * 1.d-6  ! micro m -> m
        tig_s = tig * 1.d-6  ! micro s -> s
        tig_0_s = tig_0 * 1.d-6  ! micro s -> s
        volumeig = 4.d0 * pi * rig_m**3 / 3.d0

        if (is_spherical == 0) then

            ! coefficients for planar geometries
            alpha = 0
            kcoef = 1

        else if (is_spherical == 1) then

            ! coefficients for spherical geometries
            alpha = 2
            kcoef = 4.d0*pi

        else

            ! coefficients for cylindrical geometries
            alpha = 1
            kcoef = 2.d0*pi

        end if

        gcoef = dgamma(1+(alpha+1)/pstep)**int(pstep/(alpha+1))


    end subroutine set_igsource


    subroutine add_igsource(time, rm, Eadd)

        ! A.Frendi et al. CST(1990)

        implicit none


        real(8),intent(in)  :: time, rm
        real(8),intent(out) :: Eadd

        real(8) :: coef


        if ( tig_0_s < time .and. time < tig_0_s + tig_s ) then

            ! coef = Eig_j / (kcoef/(alpha+1) * rig_m**(alpha+1) * tig_s)
            ! Eadd = coef * dexp(-gcoef * (rm/rig_m)**pstep)
            coef = Eig_j / (volumeig * tig_s)
            Eadd = coef * dexp(-(pi/4.d0) * (rm/rig_m)**6)

        else

            Eadd = 0.d0

        end if

    end subroutine add_igsource


    subroutine add_igsource_WZhang(time, rm, Eadd)

        ! W.Zhang et al. CNF(2012)

        implicit none


        real(8),intent(in)  :: time, rm
        real(8),intent(out) :: Eadd

        real(8) :: coef


        if ( tig_0_s < time .and. time < tig_0_s + tig_s ) then

            coef = Eig_j / (pi**1.5 * rig_m**3 * tig_s)
            Eadd = coef * dexp(-(rm/rig_m)**2)

        else

            Eadd = 0.d0

        end if

    end subroutine add_igsource_WZhang


end module ignitionsource


