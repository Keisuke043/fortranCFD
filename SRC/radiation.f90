
module radiation

    ! This module was firstly added by T. Akiba on Jul. 15, 2022.
    ! The coding was refered to the one by A. Tsunoda in Cantera


    use thermo
    use module_params
    
    implicit none
    

    integer :: nradsp                                   ! the number of radiative species
    real(8), allocatable :: PMACs(:, :)                 ! coefficients array of radiative parameters
    integer, allocatable :: radsp_index(:)              ! index list of radiative species
    character(len=128), allocatable :: radsp_name(:)    ! name list of radiative species


contains

    subroutine set_radiation_params(nspe, radparam_file)

        ! This subroutine reads the parameters of radiative absorption coefficients.
        ! Input : 
        !   - nspe          : number of chemical species in calculation
        !   - radparam_file : file name of radiation parameter

        implicit none

        integer, intent(in) :: nspe
        character(len=128), intent(in) :: radparam_file

        integer :: k_sp, k_rsp, fi          ! index for do loop
        integer :: lun                      ! file number indicator
        character(len=128) :: radsp_n       ! dammy radiative species name 
        real(8) :: p1, p2, p3, p4, p5, p6   ! dammy for radiative specied parameter

        ! First file reading for counting the number of radiative species
        open(newunit=lun, file=trim(radparam_file), status="old", action="read")
        nradsp = 0
        do 
            read(lun, *, end=101) radsp_n, p1, p2, p3, p4, p5, p6
            nradsp = nradsp + 1 ! count up the number of radiative species
        end do

        101 continue
        rewind(lun)
        
        ! Allocate arrays based on number of radiative species
        allocate(PMACs(6, nradsp))
        allocate(radsp_index(nradsp))
        allocate(radsp_name(nradsp))

        ! preparing the radiative heat loss parametes
        do fi = 1, nradsp
            read(lun, *) radsp_name(fi), PMACs(1, fi), PMACs(2, fi), PMACs(3, fi),&
                                         PMACs(4, fi), PMACs(5, fi), PMACs(6, fi)
        end do
        close(lun)

        ! Making the index list of radiative species
        do k_sp = 1, nspe
            do k_rsp = 1, nradsp
                if (species_name(k_sp) == radsp_name(k_rsp)) then
                    radsp_index(k_rsp) = k_sp
                end if
            end do
        end do

    end subroutine set_radiation_params


    subroutine radiation_heat_loss(p, T, nspe, x_sp, Qrad)

        ! This subroutine evaluate the radiative heat loss under Optically Thin assumption.
        ! Input:
        !   - p     : pressure
        !   - T     : temperature
        !   - nspe  : number of chemical species in calculation
        !   - x_sp  : mole fraction of all species
        ! Output:
        !   - Qrad : radiative heat loss with unit of energy

        implicit none

        integer,intent(in)  :: nspe
        real(8),intent(in)  :: p, T
        real(8),intent(in)  :: x_sp(nspe)
        real(8),intent(out) :: Qrad

        integer :: k_rsp
        real(8) :: stefan_boltzmann = 5.6703744191844294d-08 ! [W / (m^2 K^4)]
        real(8) :: T1, T2, T3, T4, T5, T0, T04
        real(8) :: Kp
        real(8) :: kpi

        ! preparing the power of temperature
        T1 = T
        T2 = T1 * T1
        T3 = T2 * T1
        T4 = T3 * T1
        T5 = T4 * T1
        T0 = 300.d0
        T04 = T0 * T0 * T0 * T0

        ! loop for radiative species
        Kp = 0.d0
        do k_rsp = 1, nradsp
            kpi = PMACs(1, k_rsp)      + PMACs(2, k_rsp) * T1 + PMACs(3, k_rsp) * T2 &
                + PMACs(4, k_rsp) * T3 + PMACs(5, k_rsp) * T4 + PMACs(6, k_rsp) * T5

            Kp = Kp + kpi * x_sp(radsp_index(k_rsp)) * (p / 101325.d0)
            ! print *, species_name(radsp_index(k_rsp)), kpi
        end do
        ! print *, Kp
    
        ! Returning the radiation heat loss
        Qrad = -4.d0 * Kp * stefan_boltzmann * (T4 - T04)

    end subroutine radiation_heat_loss

end module radiation

