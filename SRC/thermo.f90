MODULE VAR_PLOG
  PARAMETER (LENI_PLOG=200000, LENR_PLOG=99000000)
  DIMENSION I_PLOG(LENI_PLOG)
  REAL*8  R_PLOG(LENR_PLOG) !ADD FOR PLOG BY GXL,  MUST be REAL*8
  COMMON /pressure/PinATM
END MODULE VAR_PLOG

module thermo
    implicit none

    integer                    :: n_species
    integer                    :: n_con_spe
    real(8)                    :: universal_gas_constant
    
    ! CKLIB
    integer       ,allocatable :: ICKWRK(:)
    real(8)       ,allocatable :: RCKWRK(:)
    character(16) ,allocatable :: CCKWRK(:)
    real(8)                    :: RU, RUC, PA
    character(16) ,allocatable :: species_name(:)
    real(8)       ,allocatable :: molar_weight(:)
    real(8)       ,allocatable :: std_heat_formation(:)

contains
    subroutine set_thermo(chem_dir)
        

        ! local parameters
        character(len=256),intent(in) :: chem_dir
        character(len=256) :: filename_chem
        integer,parameter  :: LINC = 25
        integer,parameter  :: LOUT = 6
        ! local variables
        integer                 :: LENIWK
        integer                 :: LENRWK
        integer                 :: LENCWK
        integer                 :: MM
        integer                 :: KK
        integer                 :: II
        integer                 :: NFIT
        logical                 :: KERR

        filename_chem = trim(chem_dir)//'/chem.bin'

        ! open the CKLIB link file
        open(LINC, file=trim(filename_chem), form='UNFORMATTED', status='UNKNOWN')
        call CKLEN( LINC, LOUT, LENIWK, LENRWK, LENCWK)
        allocate(ICKWRK(LENIWK))
        allocate(RCKWRK(LENRWK))
        allocate(CCKWRK(LENCWK))
        
        call CKINIT(LENIWK, LENRWK, LENCWK, LINC, LOUT, ICKWRK, RCKWRK, CCKWRK)
        call CKINDX(ICKWRK, RCKWRK, MM, KK, II, NFIT)
        n_species = KK
        allocate(species_name(KK))
        allocate(molar_weight(KK))
        allocate(std_heat_formation(KK))


        call CKSYMS(CCKWRK, LOUT, species_name, KERR)
        call CKWT(ICKWRK, RCKWRK, molar_weight)
        call CKHMS(298.15d0, ICKWRK, RCKWRK, std_heat_formation)
        
        call CKRP(ICKWRK, RCKWRK, RU, RUC, PA)
        universal_gas_constant = RU * 1.d-4 ! ergs/(mole*K) -> J/(Kmole*K)
        close(LINC)
    end subroutine set_thermo

    subroutine Yi2Xi( n, Yi, Xi )
        integer ,intent(in)  :: n
        real(8) ,intent(in)  :: Yi(n)
        real(8) ,intent(out) :: Xi(n)
        call CKYTX( Yi, ICKWRK, RCKWRK, Xi )
    end subroutine Yi2Xi

    subroutine Xi2Yi( n, Xi, Yi )
        integer,intent(in)  :: n
        real(8),intent(in)  :: Xi(n)
        real(8),intent(out) :: Yi(n)
        call CKXTY( Xi, ICKWRK, RCKWRK, Yi )
    end subroutine Xi2Yi

    function heat_release_rate( n, Yi0, Yi, rho, dt )

        integer,intent(in) :: n
        real(8),intent(in) :: Yi0(n) 
        real(8),intent(in) :: Yi(n)
        real(8),intent(in) :: dt
        real(8),intent(in) :: rho
        real(8) :: rho_cgs
        real(8) :: heat_release_rate_cgs
        real(8) :: heat_release_rate

        ! heat_release_rate = -1.d-7*sum(molar_weight(1:n)*std_heat_formation(1:n) &
        !                                                 *(Yi(1:n)-Yi0(1:n))/dt)

        !>>>>> HRR modified by A.H>>>>>!
        rho_cgs = rho * 1.d-3  ! SI -> cgs PʊZ
        heat_release_rate_cgs = -sum(rho_cgs*std_heat_formation(1:n) &
                                            *(Yi(1:n) - Yi0(1:n)) / dt)
        heat_release_rate = heat_release_rate_cgs * 1.d-1 ! cgs -> SI PʊZ

    end function heat_release_rate

    function heat_release_rate_omegadot( n, Yi, rho, temp)

        integer,intent(in) :: n
        real(8),intent(in) :: Yi(n)
        real(8),intent(in) :: rho
        real(8),intent(in) :: temp

        integer :: i
        real(8) :: cdot(n), ddot(n)
        real(8) :: rho_cgs
        real(8) :: heat_release_rate_omegadot_cgs
        real(8) :: heat_release_rate_omegadot

        !real(8) :: c(n), d(n)
        !call ckcdyr(rho, temp, Yi, ICKWRK, RCKWRK, cdot, ddot)
        !do i = 1, n
        !    c(i) = molar_weight(i) * cdot(i) / rho
        !    d(i) = molar_weight(i) * ddot(i) / rho
        !end do
        !heat_release_rate_omegadot = -1.d-7*sum(molar_weight(1:n)* &
        !                                    std_heat_formation(1:n)*(c(1:n)-d(1:n)))
                                                                        
        !>>>>> HRR modified by A.H>>>>>!
        rho_cgs = rho * 1.d-3  ! SI -> cgs PʊZ
        call ckcdyr(rho_cgs, temp, Yi, ICKWRK, RCKWRK, cdot, ddot)
        heat_release_rate_omegadot_cgs = -sum(molar_weight(1:n)*std_heat_formation(1:n) &
                                                               *(cdot(1:n) - ddot(1:n)))
        heat_release_rate_omegadot = heat_release_rate_omegadot_cgs * 1.d-1 ! cgs -> SI PʊZ

    end function heat_release_rate_omegadot

    subroutine adjust_sum( n, vec )
        integer,intent(in)    :: n
        real(8),intent(inout) :: vec(n)
        ! Local variables
        real(8)               :: vec_sum
        integer               :: i

        vec_sum = 0.d0
        do i = 1, n
          vec(i) = dmax1( vec(i), 0.d0 )
          vec_sum = vec_sum + vec(i)
        end do
        vec(1:n) = vec(1:n)/vec_sum
    end subroutine adjust_sum

    function density( p, T, n, Xi )
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: density
        ! Local variables
        real(kind=8)            :: p_cgs
        real(kind=8)            :: rho_cgs

        p_cgs = p * 10.d0 ! SI -> cgs PʊZ
        call CKRHOX( p_cgs, T, Xi, ICKWRK, RCKWRK, rho_cgs )
        density = rho_cgs * 1.d3 ! cgs -> SI PʊZ
    end function density

    function pressure( rho, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: pressure
        ! Local variables
        real(kind=8)       :: p_cgs
        real(kind=8)       :: rho_cgs

        rho_cgs = rho * 1.d-3 ! SI -> cgs PʊZ
        call CKPX( rho_cgs, T, Xi, ICKWRK, RCKWRK, p_cgs )
        pressure = p_cgs * 0.1d0 ! cgs -> SI PʊZ
    end function pressure

    function temperature( rho, p, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: temperature
        ! Local variables
        real(kind=8)            :: wtm
        integer                 :: i

        wtm = 0.d0
        do i = 1, n
          wtm = wtm + Xi(i)*molar_weight(i)
        end do
        temperature = p/(universal_gas_constant/wtm*rho)
    end function temperature

    function pressure_con( rho, specificE, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: specificE
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: pressure_con
        ! Local variables
        real(8)            :: T

        T = temperature_con( specificE, n, Xi )
        pressure_con = pressure( rho, T, n, Xi )
    end function pressure_con

    function temperature_con( specificE, n, Xi )
        real(8),intent(in) :: specificE
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: temperature_con
        ! Local parameters
        real(kind=8) ,parameter :: TEMP_THRE_REDI = 1.d-9
        real(kind=8) ,parameter :: TEMP_THRE_ERRO = 1.d-6
        integer      ,parameter :: NR_ITR_MAX     = 10
        real(kind=8) ,parameter :: TEMP_MIN       = 20.d0
        real(kind=8) ,parameter :: TEMP_MAX       = 6000.d0
        ! Local variables
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: spE_cgs
        real(kind=8)            :: UBML
        real(kind=8)            :: CVBML
        integer                 :: NR_iteration
        real(kind=8)            :: func_x
        real(kind=8)            :: Temp, dTemp

        wtm = 0.d0
        do i = 1, n
          wtm = wtm + Xi(i)*molar_weight(i)
        end do
        spE_cgs = specificE * 1.d4 * wtm ! SI -> cgs PʊZ

        ! Temp = current temperature
        Temp    = 1000.d0

        func_x       = 1.d30
        dTemp        = 1.d30
        NR_iteration = 0
        do while( dabs(func_x) >= dabs(spE_cgs)*TEMP_THRE_REDI &
&       .AND. dabs(dTemp) >= Temp*TEMP_THRE_ERRO          &
&       .AND. NR_iteration < NR_ITR_MAX                  )
            call CKUBML( Temp, Xi, ICKWRK, RCKWRK, UBML )
            call CKCVBL( Temp, Xi, ICKWRK, RCKWRK, CVBML )
            func_x = UBML - spE_cgs
            dTemp  = -func_x/CVBML
            Temp   = Temp + dTemp

            Temp = dmax1( TEMP_MIN, dmin1( TEMP_MAX, Temp ) )
            NR_iteration = NR_iteration + 1
        end do
        temperature_con = Temp
    end function temperature_con

    function temperature_con2( specificE, n, Xi, Temp )
        real(8),intent(in) :: specificE
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: temperature_con2
        ! Local parameters
        real(kind=8) ,parameter :: TEMP_THRE_REDI = 1.d-9
        real(kind=8) ,parameter :: TEMP_THRE_ERRO = 1.d-6
        integer      ,parameter :: NR_ITR_MAX     = 10
        real(kind=8) ,parameter :: TEMP_MIN       = 20.d0
        real(kind=8) ,parameter :: TEMP_MAX       = 6000.d0
        ! Local variables
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: spE_cgs
        real(kind=8)            :: UBML
        real(kind=8)            :: CVBML
        integer                 :: NR_iteration
        real(kind=8)            :: func_x
        real(kind=8)            :: Temp, dTemp

        wtm = 0.d0
        do i = 1, n
            wtm = wtm + Xi(i)*molar_weight(i)
        end do
        spE_cgs = specificE * 1.d4 * wtm ! SI -> cgs PʊZ

        func_x       = 1.d30
        dTemp        = 1.d30
        NR_iteration = 0
        do while( dabs(func_x) >= dabs(spE_cgs)*TEMP_THRE_REDI &
&         .AND. dabs(dTemp) >= Temp*TEMP_THRE_ERRO          &
&         .AND. NR_iteration < NR_ITR_MAX                  )
            call CKUBML( Temp, Xi, ICKWRK, RCKWRK, UBML )
            call CKCVBL( Temp, Xi, ICKWRK, RCKWRK, CVBML )
            func_x = UBML - spE_cgs
            dTemp  = -func_x/CVBML
            Temp   = Temp + dTemp

            Temp = dmax1( TEMP_MIN, dmin1( TEMP_MAX, Temp ) )
            NR_iteration = NR_iteration + 1
        end do
        temperature_con2 = Temp
    end function temperature_con2

    function Cv( rho, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: Cv
        ! Local variables
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: CVBML

        wtm = 0.d0
        do i = 1, n
            wtm = wtm + Xi(i)*molar_weight(i)
        end do
        call CKCVBL( T, Xi, ICKWRK, RCKWRK, CVBML )

        Cv = CVBML * 1.d-4 / wtm ! cgs -> SI PʊZ
    end function Cv

    function Cp( rho, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: Cp
        ! Local variables
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: CPBML

        wtm = 0.d0
        do i = 1, n
            wtm = wtm + Xi(i)*molar_weight(i)
        end do
        call CKCPBL( T, Xi, ICKWRK, RCKWRK, CPBML )

        Cp = CPBML * 1.d-4 / wtm ! cgs -> SI PʊZ
    end function Cp

    function heat_capacity_ratio( rho, p, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: heat_capacity_ratio

        heat_capacity_ratio = Cp( rho, T, n, Xi ) &
&                           / Cv( rho, T, n, Xi )
    end function heat_capacity_ratio

    function specific_energy( rho, p, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: specific_energy
        ! Local variables 
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: UBML

        ! print *,  rho, p, T, n, Xi
        wtm = 0.d0
        do i = 1, n
            wtm = wtm + Xi(i)*molar_weight(i)
        end do
        call CKUBML( T, Xi, ICKWRK, RCKWRK, UBML )
        specific_energy = UBML * 1.d-4 / wtm ! cgs -> SI PʊZ
    end function specific_energy

    subroutine specific_Hi( rho, p, T, n, Hi )
        real(8),intent(in)  :: rho
        real(8),intent(in)  :: p
        real(8),intent(in)  :: T
        integer,intent(in)  :: n
        real(8),intent(out) :: Hi(n)

        call CKHMS( T, ICKWRK, RCKWRK, Hi )

        Hi(1:n) = Hi(1:n) * 1.d-4 ! cgs -> SI PʊZ
    end subroutine specific_Hi

    function specific_enthalpy( rho, p, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: specific_enthalpy
        ! Local variables 
        real(kind=8)            :: wtm
        integer                 :: i
        real(kind=8)            :: HBML

        wtm = 0.d0
        do i = 1, n
            wtm = wtm + Xi(i)*molar_weight(i)
        end do
        call CKHBML( T, Xi, ICKWRK, RCKWRK, HBML )

        specific_enthalpy = HBML * 1.d-4 / wtm ! cgs -> SI PʊZ
    end function specific_enthalpy

    function sonic_speed2( rho, p, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: sonic_speed2

        sonic_speed2 = heat_capacity_ratio( rho, p, T, n, Xi ) * p / rho
    end function sonic_speed2

    function sonic_speed( rho, p, T, n, Xi )
        real(8),intent(in) :: rho
        real(8),intent(in) :: p
        real(8),intent(in) :: T
        integer,intent(in) :: n
        real(8),intent(in) :: Xi(n)
        real(8)            :: sonic_speed

        sonic_speed = sqrt( sonic_speed2( rho, p, T, n, Xi ) )
    end function sonic_speed
end module thermo
