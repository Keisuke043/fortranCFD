module eos

    use thermo
    use variables

    implicit none


contains

    function qfToQc(qf, nspe, nequ)
        
        implicit none

    
        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: qf(nequ)

        integer :: k_sp
        real(8) :: pres, E_internal
        real(8) :: x_sp(nspe), y_sp(nspe)

        real(8) :: Qc(nequ)
        real(8) :: qfToQc(nequ)


        Qc(iudens) = qf(iqdens) ! rho
        Qc(iumomx) = qf(iqdens)*qf(iqxvel) ! rho*u

        y_sp = qf(nprim+1:nequ)
        call adjust_sum(nspe, y_sp)

        do k_sp = 1, nspe
            Qc(ncons+k_sp) = qf(iqdens)*y_sp(k_sp) ! rho*Yi
        end do

        call Yi2Xi(nspe, qf(nprim+1:nequ), x_sp)
        call adjust_sum(nspe, x_sp)

        pres = pressure(qf(iqdens), qf(iqtemp), nspe, x_sp)
        E_internal = specific_energy(qf(iqdens), pres, qf(iqtemp), nspe, x_sp)
        Qc(iuener) = qf(iqdens) * (E_internal + 0.5d0*qf(iqxvel)**2.d0) ! rho*e

        qfToQc = Qc

    end function qfToQc


    function QcToqf(Qc, nspe, nequ)

        implicit none


        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: Qc(nequ)

        integer :: k_sp
        real(8) :: E_internal
        real(8) :: x_sp(nspe)

        real(8) :: qf(nequ)
        real(8) :: QcToqf(nequ)


        qf = 0.d0

        qf(iqdens) = Qc(iudens) ! rho
        qf(iqxvel) = Qc(iumomx)/Qc(iudens) ! rho*u/rho
        do k_sp = 1, nspe
            qf(nprim+k_sp) = Qc(ncons+k_sp)/Qc(iudens) ! rho*Yi/rho
        end do

        call adjust_sum(nspe, qf(nprim+1:nequ))

        call Yi2Xi(nspe, qf(nprim+1:nequ), x_sp)
        call adjust_sum(nspe, x_sp)

        E_internal = Qc(iuener)/Qc(iudens) - 0.5d0*qf(iqxvel)**2.d0
        qf(iqtemp) = temperature_con(E_internal, nspe, x_sp)

        QcToqf = qf

    end function QcToqf


end module eos

