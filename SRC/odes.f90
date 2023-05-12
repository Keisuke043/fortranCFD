MODULE ODEs
  IMPLICIT NONE

CONTAINS
  SUBROUTINE MACKS(dtcfd, n, yeq, rho)
    USE thermo
    USE module_params, only: atol, rtol
    IMPLICIT NONE
    INTEGER :: i, j, k, l, m, n
    REAL(8) :: yeq(n), y1(n), ysum, c(n), d(n), yt(n)
    REAL(8) :: t_now, dtcfd
    REAL(8) :: u(n),e0,e1,cv_t,temp1,dtemp,rmax
    REAL(8) :: beta,dt1,dtmax,ymax,tau(n)
    REAL(8) :: y0(n), c0(n), tau0(n), d0(n)
    REAL(8) :: lambda
    REAL(8) :: rho, rho_cgs
    REAL(8) :: expo, alpha
    REAL(8) :: rmaxold
    ! REAL(8) :: atol, rtol
    
    ! atol = 1.d-10
    ! rtol = 1.d-4
    rho_cgs = rho*1d-3 

    ! do i = 1, n
    !     print *,yeq(i) 
    ! end do
    ysum = 0d0
    DO i = 1,n-1
       yeq(i) = dmax1(yeq(i),0.d0)
       ysum = ysum+yeq(i)
    END DO
    DO i = 1,n-1
       yeq(i) = yeq(i)/ysum
       ! if (yeq(i) /= yeq(i)) then
       !     print *,i, yeq(i), rmax, dt1, t_now
       !     stop
       ! end if
    END DO
    
    CALL CKUMS(yeq(n), ICKWRK, RCKWRK, u)
    e0=SUM(u(1:n-1)*yeq(1:n-1))

    t_now = 0.0d0
    dtmax = dtcfd
    ! print *, dtcfd, l, t_now,  dt1, n, yeq, rho
    DO
       l = 0
       rmaxold = 1d100
       DO
          l = l + 1
          CALL rhs(n, yeq, c, d, rho_cgs)
          DO i=1,n-1
             tau(i) = (yeq(i)+1d-100)/(d(i)+1d-50)
          ENDDO
          IF(l==1)THEN
             dt1 = dmin1(dtmax, dtcfd-t_now)
             y0 = yeq
             c0 = c
             d0 = d
             tau0 = tau
          ENDIF
200       CONTINUE
          rmax = 0.0d0
          DO i = 1, n-1
             alpha = dt1/tau(i)
             IF ( alpha > 1d2 ) THEN
                lambda = 0.0d0
             ELSE IF (alpha < 1d-2) THEN
                lambda = 0.5d0
             ELSE
                expo = dexp(-alpha)
                lambda = (1.0d0-(1.0d0+alpha)*expo)/((1.0d0-expo)*alpha)
             ENDIF
             beta = tau(i)*dt1/(tau(i)+(1.0d0-lambda)*dt1)
             y1(i) = yeq(i) + beta*(lambda*(c0(i)-y0(i)/tau0(i)) &
                   + (1.0d0-lambda)*c(i) + y0(i)/dt1 &
                   - (1.0d0/dt1+(1.0d0-lambda)/tau(i))*yeq(i))
             ! rmax = rmax + dabs(y1(i)-yeq(i))/(atol+rtol*dmax1(y1(i),yeq(i)))
             rmax = dmax1(rmax, dabs(y1(i)-yeq(i))/(atol+rtol*dmax1(y1(i),yeq(i))))
             yeq(i) = y1(i)
          ENDDO

          ! IF ( rmax - rmaxold > 1.0d0) THEN
          ! print *, rmax, dt1, t_now
          IF ( rmax - rmaxold > 100.0d0) THEN
             yeq = y0
             tau = tau0
             c = c0
             d = d0
             dtmax = dtmax*0.25d0
             dt1 = dmin1(dtmax, dtcfd-t_now)
             rmaxold = 1d100
             GOTO 200
          ENDIF
          rmaxold = rmax
          temp1=yeq(n)
          DO k = 1, 10
             CALL CKUMS(temp1, ICKWRK, RCKwRK, u )
             CALL ckcvbs(temp1, yeq, ICKWRK, RCKWRK, cv_t)
             e1 = -e0
             DO i=1,n-1
                e1 = e1 + u(i)*yeq(i)
             ENDDO
             !dtemp = -e1/cv(rho, yeq(n), n, )
             dtemp = -e1/cv_t
             temp1 = temp1+dtemp
             IF (DABS(dtemp)<1.e-6) EXIT
          ENDDO
          ! k=10
          yeq(n) = temp1
          ! print *, t_now,  dt1, l, rmax
          IF (l==1.AND.rmax<1.0d0.AND.dtmax < dtcfd) THEN
             dtmax = 2.0d0*dtmax
             ! dtmax = 2.1d0*dtmax
             dt1 = dmin1(dtmax, dtcfd-t_now)
          ENDIF
          IF (rmax<1d0) EXIT
       END DO
       ! print *, 'y0', y0
       ! print *, 'yeq', yeq

       ! print *, l, t_now,  dt1
       t_now = t_now + dt1
       IF(DABS(t_now-dtcfd)<1d-15) EXIT
       IF(t_now>dtcfd) EXIT
    ENDDO
  END SUBROUTINE MACKS

  SUBROUTINE rhs(neq, yeq, c, d, rho)
    USE thermo
    IMPLICIT NONE
    INTEGER i
    INTEGER :: neq
    REAL(8) yeq(neq), c(neq), d(neq)
    REAL(8) u_v(n_species), cv_v, avu_c, avu_d, cdot(n_species), ddot(n_species)
    REAL(8) temp, rho

    temp = yeq(neq)
    CALL ckcdyr( rho, temp, yeq, ICKWRK, RCKWRK, cdot, ddot )
    CALL ckums( temp, ICKWRK, RCKWRK, u_v )
    CALL ckcvbs( temp, yeq, ICKWRK, RCKWRK, cv_v )
    avu_c = 0.0d0
    avu_d = 0.0d0
    DO i = 1,n_species
       avu_c = avu_c+u_v(i)*molar_weight(i)*cdot(i)
       avu_d = avu_d+u_v(i)*molar_weight(i)*ddot(i)
       c(i) = molar_weight(i)*cdot(i)/rho
       d(i) = molar_weight(i)*ddot(i)/rho
    END DO
    c(neq) = -avu_c/rho/cv_v
    d(neq) = -avu_d/rho/cv_v
  END SUBROUTINE rhs

END MODULE ODEs
