
module muscl

    implicit none


contains

    subroutine muscl3(s, sl, sr, flag_muscl)

        implicit none


        integer,intent(in)  :: flag_muscl
        real(8),intent(in)  :: s(4)
        real(8),intent(out) :: sl, sr

        real(8) :: dplus_i, dminus_i, dplus_ip, dminus_ip
        real(8) :: limiter_m_i,  limiter_p_i
        real(8) :: limiter_m_ip, limiter_p_ip

        real(8),parameter :: k_muscl = 1.d0/3.d0   ! Accuracy: 2up:k=0.d0, 2:k=1.d0, 3:k = 1.d0/3.d0
        real(8) :: b_muscl


        b_muscl = (3.d0-k_muscl)/(1.d0-k_muscl)

        ! print *, k_muscl, b_muscl

        dminus_i  = s(2) - s(1)
        dplus_i   = s(3) - s(2)
        dminus_ip = s(3) - s(2)
        dplus_ip  = s(4) - s(3)

        if ( flag_muscl == 1 ) then

            sl = s(2) + (1.d0/4.d0) * &
                 ((1.d0-k_muscl)*dminus_i + (1.d0+k_muscl)*dplus_i)

            sr = s(3) - (1.d0/4.d0) * &
                 ((1.d0-k_muscl)*dplus_ip + (1.d0+k_muscl)*dminus_ip)

        else if ( flag_muscl == 2 .or. flag_muscl == 3 .or. flag_muscl == 4 ) then

            if ( flag_muscl == 2 ) then

                limiter_p_i  = minmod(dplus_i,   b_muscl*dminus_i)
                limiter_m_i  = minmod(dminus_i,  b_muscl*dplus_i)
                limiter_p_ip = minmod(dplus_ip,  b_muscl*dminus_ip)
                limiter_m_ip = minmod(dminus_ip, b_muscl*dplus_ip)

            else if ( flag_muscl == 3 ) then

                limiter_p_i  = superbee(dplus_i,   b_muscl*dminus_i)
                limiter_m_i  = superbee(dminus_i,  b_muscl*dplus_i)
                limiter_p_ip = superbee(dplus_ip,  b_muscl*dminus_ip)
                limiter_m_ip = superbee(dminus_ip, b_muscl*dplus_ip)

            else if ( flag_muscl == 4 ) then

                limiter_p_i  = mc(dplus_i,   b_muscl*dminus_i)
                limiter_m_i  = mc(dminus_i,  b_muscl*dplus_i)
                limiter_p_ip = mc(dplus_ip,  b_muscl*dminus_ip)
                limiter_m_ip = mc(dminus_ip, b_muscl*dplus_ip)

            end if

            sl = s(2) + (1.d0/4.d0) * &
                 ((1.d0-k_muscl)*limiter_m_i + (1.d0+k_muscl)*limiter_p_i)

            sr = s(3) - (1.d0/4.d0) * &
                 ((1.d0-k_muscl)*limiter_p_ip + (1.d0+k_muscl)*limiter_m_ip)

        end if

    end subroutine muscl3


    function minmod(x, y)

        ! minmod flux limiter

        implicit none


        real(8),intent(in) :: x, y
        real(8) :: minmod


        minmod = dsign(1.d0,x) * dmax1(0.d0, dmin1(dabs(x), dsign(1.d0,x) * y))


    end function minmod


    function superbee(x, y)

        ! superbee flux limiter

        implicit none


        real(8),intent(in) :: x, y
        real(8) :: superbee


        ! superbee = dsign(1.d0,y) * dmax1( 0.d0, &
        !            dmin1(dsign(1.d0,y) * b * x, dabs(y)), &
        !            dmin1(dsign(1.d0,y) * x, b * dabs(y)) )
        ! superbee = dsign(1.d0,x) * dmax1( 0.d0, &
        !            dmin1(     dabs(x), 2.d0*dsign(1.d0,x) * y), &
        !            dmin1(2.d0*dabs(x),      dsign(1.d0,x) * y) )
        ! superbee = dmax1( 0.d0, &
        !            dmin1(     dabs(x), 2.d0*y), &
        !            dmin1(2.d0*dabs(x),      y) )
        superbee = dmax1( 0.d0, &
                   dmin1(     x, 2.d0*y), &
                   dmin1(2.d0*x,      y) )


    end function superbee


    function mc(x, y)

        ! monotonized central (MC) flux limiter

        implicit none


        real(8),intent(in) :: x, y
        real(8) :: mc


        ! mc = dsign(1.d0,x) * dmax1(0.d0, & 
        !      dmin1(2.d0*dabs(x), dmin1(2.d0*dsign(1.d0,x)*y, 0.5d0*(dabs(x)+dsign(1.d0,x)*y))))
        ! mc = dmax1(0.d0, & 
        !      dmin1(2.d0*dabs(x), dmin1(2.d0*y, 0.5d0*(dabs(x)+y))))
        mc = dmax1(0.d0, dmin1(2.d0*x, 2.d0*y, 0.5d0*(x+y)))


    end function mc


endmodule muscl


