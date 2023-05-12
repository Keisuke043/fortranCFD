
module weno


    implicit none


contains


    subroutine weno3z(s, sl, sr)

        implicit none


        real(8),intent(in)  :: s(4)
        real(8),intent(out) :: sl, sr

        real(8) :: vl(2)
        real(8) :: vr(2)
        real(8) :: beta(2)
        real(8) :: alpha(2)
        real(8) :: alpha1
        real(8) :: tau

        real(8),parameter :: eps = 1.0d-6


        beta(2) = (s(1) - s(2)) * (s(1) - s(2))
        beta(1) = (s(2) - s(3)) * (s(2) - s(3))

        tau = dabs(beta(2)-beta(1))

        beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
        beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        alpha(2) = 2.d0 * beta(2)
        alpha(1) = 1.d0 * beta(1)
        alpha1 = 1.d0 / (alpha(2) + alpha(1))

        vl(2) = -s(1) + 3.d0 * s(2)
        vl(1) = s(2) + s(3)

        sl = 0.5d0 * alpha1 * (alpha(2) * vl(2) + alpha(1) * vl(1))


        beta(2) = (s(4) - s(3)) * (s(4) - s(3))
        beta(1) = (s(3) - s(2)) * (s(3) - s(2))

        tau = dabs(beta(2)-beta(1))

        beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
        beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        alpha(2) = 2.d0 * beta(2)
        alpha(1) = 1.d0 * beta(1)
        alpha1 = 1.d0 / (alpha(2) + alpha(1))

        vr(2) = 3.d0 * s(3) - s(4)
        vr(1) = s(2) + s(3)

        sr = 0.5d0 * alpha1 * (alpha(2) * vr(2) + alpha(1) * vr(1))


    end subroutine weno3z


    subroutine weno5(s, sl, sr, flag_weno)

        implicit none


        integer,intent(in)  :: flag_weno

        real(8),intent(in)  :: s(6)
        real(8),intent(out) :: sl, sr

        real(8) :: vl(3)
        real(8) :: vr(3)
        real(8) :: beta(3)
        real(8) :: alpha(3)
        real(8) :: alpha1
        real(8) :: tau

        real(8),parameter :: eps = 1.0d-6


        beta(3) = (13.d0 / 12.d0) * (s(1) - 2.d0 * s(2) + s(3))**2 + &
                  0.25d0 * (s(1) - 4.d0 * s(2) + 3.d0 * s(3))**2

        beta(2) = (13.d0 / 12.d0) * (s(2) - 2.d0 * s(3) + s(4))**2 + &
                  0.25d0 * (s(2) - s(4)) * (s(2) - s(4))

        beta(1) = (13.d0 / 12.d0) * (s(3) - 2.d0 * s(4) + s(5))**2 + &
                  0.25d0 * (3.d0 * s(3) - 4.d0 * s(4) + s(5))**2

        if ( flag_weno == 2 ) then

            ! 5th-order weno-js
            beta(3) = 1.0 / ((eps + beta(3)) * (eps + beta(3)))
            beta(2) = 1.0 / ((eps + beta(2)) * (eps + beta(2)))
            beta(1) = 1.0 / ((eps + beta(1)) * (eps + beta(1)))

        else if ( flag_weno == 3 ) then

            ! 5th-order weno-z
            tau = dabs(beta(3)-beta(1))
            beta(3) = 1.d0 + (tau / (eps + beta(3))) * (tau / (eps + beta(3)))
            beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
            beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        end if

        alpha(3) = beta(3)
        alpha(2) = 6.d0 * beta(2)
        alpha(1) = 3.d0 * beta(1)
        alpha1 = 1.d0 / (alpha(3) + alpha(2) + alpha(1))

        vl(3) = 2.d0 * s(1) - 7.d0 * s(2) + 11.d0 * s(3)
        vl(2) = -s(2) + 5.d0 * s(3) + 2.d0 * s(4)
        vl(1) = 2.d0 * s(3) + 5.d0 * s(4) - s(5)
        
        sl = (1.d0 / 6.d0) * alpha1 * &
             (alpha(3) * vl(3) + alpha(2) * vl(2) + alpha(1) * vl(1))


        beta(3) = (13.d0 / 12.d0) * (s(6) - 2.d0 * s(5) + s(4))**2 + &
                  0.25d0 * (s(6) - 4.d0 * s(5) + 3.d0 * s(4))**2

        beta(2) = (13.d0 / 12.d0) * (s(5) - 2.d0 * s(4) + s(3))**2 + &
                  0.25d0 * (s(5) - s(3)) * (s(5) - s(3))

        beta(1) = (13.d0 / 12.d0) * (s(4) - 2.d0 * s(3) + s(2))**2 + &
                  0.25d0 * (3.d0 * s(4) - 4.d0 * s(3) + s(2))**2

        if ( flag_weno == 2 ) then

            ! 5th-order weno-js
            beta(3) = 1.0 / ((eps + beta(3)) * (eps + beta(3)))
            beta(2) = 1.0 / ((eps + beta(2)) * (eps + beta(2)))
            beta(1) = 1.0 / ((eps + beta(1)) * (eps + beta(1)))

        else if ( flag_weno == 3 ) then

            ! 5th-order weno-z
            tau = dabs(beta(3)-beta(1))
            beta(3) = 1.d0 + (tau / (eps + beta(3))) * (tau / (eps + beta(3)))
            beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
            beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        end if

        alpha(3) = beta(3)
        alpha(2) = 6.d0 * beta(2)
        alpha(1) = 3.d0 * beta(1)
        alpha1 = 1.d0 / (alpha(3) + alpha(2) + alpha(1))

        vr(3) = 11.d0 * s(4) - 7.d0 * s(5) + 2.d0 * s(6)
        vr(2) = -s(5) + 5.d0 * s(4) + 2.d0 * s(3)
        vr(1) = 2.d0 * s(4) + 5.d0 * s(3) - s(2)

        sr = (1.d0 / 6.d0) * alpha1 * &
             (alpha(3) * vr(3) + alpha(2) * vr(2) + alpha(1) * vr(1))


    end subroutine weno5


    subroutine weno7z(s, sl, sr)

        implicit none


        real(8),intent(in)  :: s(8)
        real(8),intent(out) :: sl, sr

        real(8) :: vl(4)
        real(8) :: vr(4)
        real(8) :: beta(4)
        real(8) :: alpha(4)
        real(8) :: alpha1
        real(8) :: tau
        real(8) :: weno7_face_wghts_1(4)

        real(8),parameter :: eps = 1.0d-6


        weno7_face_wghts_1 = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)


        beta(4) = s(1) * (547.d0 * s(1) - 3882.d0 * s(2) + 4642.d0 * s(3) - 1854.d0 * s(4)) + &
                  s(2) * (7043.d0 * s(2) - 17246.d0 * s(3) + 7042.d0 * s(4)) + &
                  s(3) * (11003.d0 * s(3) - 9402.d0 * s(4)) + 2107.d0 * (s(4) * s(4))

        beta(3) = s(2) * (267.d0 * s(2) - 1642.d0 * s(3) + 1602.d0 * s(4) - 494.d0 * s(5)) + &
                  s(3) * (2843.d0 * s(3) - 5966.d0 * s(4) + 1922.d0 * s(5)) + &
                  s(4) * (3443.d0 * s(4) - 2522.d0 * s(5)) + 547.d0 * (s(5) * s(5))

        beta(2) = s(3) * (547.d0 * s(3) - 2522.d0 * s(4) + 1922.d0 * s(5) - 494.d0 * s(6)) + &
                  s(4) * (3443.d0 * s(4) - 5966.d0 * s(5) + 1602.d0 * s(6)) + &
                  s(5) * (2843.d0 * s(5) - 1642.d0 * s(6)) + 267.d0 * (s(6) * s(6))

        beta(1) = s(4) * (2107.d0 * s(4) - 9402.d0 * s(5) + 7042.d0 * s(6) - 1854.d0 * s(7)) + &
                  s(5) * (11003.d0 * s(5) - 17246.d0 * s(6) + 4642.d0 * s(7)) + &
                  s(6) * (7043.d0 * s(6) - 3882.d0 * s(7)) + 547.d0 * (s(7) * s(7))

        tau = dabs(beta(4)-beta(1))

        beta(4) = 1.d0 + (tau / (eps + beta(4))) * (tau / (eps + beta(4)))
        beta(3) = 1.d0 + (tau / (eps + beta(3))) * (tau / (eps + beta(3)))
        beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
        beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        alpha(4) = weno7_face_wghts_1(1) * beta(4)
        alpha(3) = weno7_face_wghts_1(2) * beta(3)
        alpha(2) = weno7_face_wghts_1(3) * beta(2)
        alpha(1) = weno7_face_wghts_1(4) * beta(1)
        alpha1 = 1.d0 / (alpha(4) + alpha(3) + alpha(2) + alpha(1))

        vl(4) = -3.d0 * s(1) + 13.d0 * s(2) - 23.d0 * s(3) + 25.d0 * s(4)
        vl(3) =  1.d0 * s(2) -  5.d0 * s(3) + 13.d0 * s(4) +  3.d0 * s(5)
        vl(2) = -1.d0 * s(3) +  7.d0 * s(4) +  7.d0 * s(5) -  1.d0 * s(6)
        vl(1) =  3.d0 * s(4) + 13.d0 * s(5) -  5.d0 * s(6) +  1.d0 * s(7)

        sl = (1.d0 / 12.d0) * alpha1 * &
             (alpha(4) * vl(4) + alpha(3) * vl(3) + alpha(2) * vl(2) + alpha(1) * vl(1))


        beta(4) = s(8) * (547.d0 * s(8) - 3882.d0 * s(7) + 4642.d0 * s(6) - 1854.d0 * s(5)) + &
                  s(7) * (7043.d0 * s(7) - 17246.d0 * s(6) + 7042.d0 * s(5)) + &
                  s(6) * (11003.d0 * s(6) - 9402.d0 * s(5)) + 2107.d0 * (s(5) * s(5))

        beta(3) = s(7) * (267.d0 * s(7) - 1642.d0 * s(6) + 1602.d0 * s(5) - 494.d0 * s(4)) + &
                  s(6) * (2843.d0 * s(6) - 5966.d0 * s(5) + 1922.d0 * s(4)) + &
                  s(5) * (3443.d0 * s(5) - 2522.d0 * s(4)) + 547.d0 * (s(4) * s(4))

        beta(2) = s(6) * (547.d0 * s(6) - 2522.d0 * s(5) + 1922.d0 * s(4) - 494.d0 * s(3)) + &
                  s(5) * (3443.d0 * s(5) - 5966.d0 * s(4) + 1602.d0 * s(3)) + &
                  s(4) * (2843.d0 * s(4) - 1642.d0 * s(3)) + 267.d0 * (s(3) * s(3))

        beta(1) = s(5) * (2107.d0 * s(5) - 9402.d0 * s(4) + 7042.d0 * s(3) - 1854.d0 * s(2)) + &
                  s(4) * (11003.d0 * s(4) - 17246.d0 * s(3) + 4642.d0 * s(2)) + &
                  s(3) * (7043.d0 * s(3) - 3882.d0 * s(2)) + 547.d0 * (s(2) * s(2))

        tau = dabs(beta(4)-beta(1))

        beta(4) = 1.d0 + (tau / (eps + beta(4))) * (tau / (eps + beta(4)))
        beta(3) = 1.d0 + (tau / (eps + beta(3))) * (tau / (eps + beta(3)))
        beta(2) = 1.d0 + (tau / (eps + beta(2))) * (tau / (eps + beta(2)))
        beta(1) = 1.d0 + (tau / (eps + beta(1))) * (tau / (eps + beta(1)))

        alpha(4) = weno7_face_wghts_1(1) * beta(4)
        alpha(3) = weno7_face_wghts_1(2) * beta(3)
        alpha(2) = weno7_face_wghts_1(3) * beta(2)
        alpha(1) = weno7_face_wghts_1(4) * beta(1)
        alpha1 = 1.d0 / (alpha(4) + alpha(3) + alpha(2) + alpha(1))

        vr(4) = -3.d0 * s(8) + 13.d0 * s(7) - 23.d0 * s(6) + 25.d0 * s(5)
        vr(3) =  1.d0 * s(7) -  5.d0 * s(6) + 13.d0 * s(5) +  3.d0 * s(4)
        vr(2) = -1.d0 * s(6) +  7.d0 * s(5) +  7.d0 * s(4) -  1.d0 * s(3)
        vr(1) =  3.d0 * s(5) + 13.d0 * s(4) -  5.d0 * s(3) +  1.d0 * s(2)

        sr = (1.d0 / 12.d0) * alpha1 * &
             (alpha(4) * vr(4) + alpha(3) * vr(3) + alpha(2) * vr(2) + alpha(1) * vr(1))


    end subroutine weno7z


end module weno


