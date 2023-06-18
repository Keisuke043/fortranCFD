module flux

    use thermo
    use trans
    use muscl
    use weno
    use module_boundary
    use eos
    use block_t

    implicit none

    !!! FVS (Flux vector splitting) !!!
    ! real(8),allocatable :: R(:,:), R_inv(:,:), Lam(:,:), Lam_abs(:,:)

contains

    subroutine calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)


        implicit none

        integer :: i, ib, k_sp
        integer,intent(in) :: nspe, nequ, lvir
        integer,intent(in) :: bsta, bend
        type(block_st) :: b_obj(bsta:bend)


        real(8) :: rho_i, rho_ip
        real(8) :: u_i, u_ip
        real(8) :: y_sp_i(nspe), y_sp_ip(nspe)
        real(8) :: x_sp_i(nspe), x_sp_ip(nspe)
        real(8) :: eTotal_i, eTotal_ip
        real(8) :: eKinetic_i, eKinetic_ip
        real(8) :: eInternal_i, eInternal_ip
        real(8) :: T_i, T_ip
        real(8) :: p_i, p_ip

        real(8) :: mu_i, mu_ip
        real(8) :: kap_i, kap_ip
        real(8) :: dm_i(nspe), dm_ip(nspe)
        real(8) :: dif_i(nspe), dif_ip(nspe)
        real(8) :: rho_ave, T_ave, u_ave
        real(8) :: dm(nspe), hs(nspe), hsum
        real(8) :: mu_ave, kap_ave, dif_ave, tau_xx


        do ib = bsta, bend
            do i = b_obj(ib)%ista_b-1, b_obj(ib)%iend_b

                rho_i  = b_obj(ib)%Qc(iudens,i)
                rho_ip = b_obj(ib)%Qc(iudens,i+1)
                u_i    = b_obj(ib)%Qc(iumomx,i)  /rho_i
                u_ip   = b_obj(ib)%Qc(iumomx,i+1)/rho_ip

                do k_sp = 1, nspe
                    y_sp_i(k_sp)  = b_obj(ib)%Qc(ncons+k_sp,i)  /rho_i
                    y_sp_ip(k_sp) = b_obj(ib)%Qc(ncons+k_sp,i+1)/rho_ip
                end do

                call Yi2Xi(nspe, y_sp_i(:),  x_sp_i(:))
                call Yi2Xi(nspe, y_sp_ip(:), x_sp_ip(:))
                eTotal_i     = b_obj(ib)%Qc(iuener,i)  /rho_i
                eTotal_ip    = b_obj(ib)%Qc(iuener,i+1)/rho_ip
                eKinetic_i   = 0.5d0*(u_i**2d0)
                eKinetic_ip  = 0.5d0*(u_ip**2d0)
                eInternal_i  = eTotal_i  - eKinetic_i
                eInternal_ip = eTotal_ip - eKinetic_ip
                T_i    = temperature_con(eInternal_i, nspe,x_sp_i)
                T_ip   = temperature_con(eInternal_ip,nspe,x_sp_ip)
                p_i    = pressure(rho_i, T_i, nspe,x_sp_i)
                p_ip   = pressure(rho_ip,T_ip,nspe,x_sp_ip)
                mu_i   = viscosity_ea(rho_i, p_i ,T_i, nspe,x_sp_i)
                mu_ip  = viscosity_ea(rho_ip,p_ip,T_ip,nspe,x_sp_ip)
                kap_i  = thermal_conductivity(rho_i, p_i, T_i, nspe,x_sp_i)
                kap_ip = thermal_conductivity(rho_ip,p_ip,T_ip,nspe,x_sp_ip)
                call set_diffusivity(rho_i, p_i, T_i, nspe,x_sp_i, dm_i)
                call set_diffusivity(rho_ip,p_ip,T_ip,nspe,x_sp_ip,dm_ip)

                do k_sp = 1, nspe

                    dif_i(k_sp)  = (1.0d0-x_sp_i(k_sp)) /(1.0d0-y_sp_i(k_sp))*dm_i(k_sp)
                    dif_ip(k_sp) = (1.0d0-x_sp_ip(k_sp))/(1.0d0-y_sp_ip(k_sp))*dm_ip(k_sp)

                    if (y_sp_i(k_sp) == 1.0d0) then
                        dif_i(k_sp) = 0.0d0
                    end if
                    if (y_sp_ip(k_sp) == 1.0d0) then
                        dif_ip(k_sp) = 0.0d0
                    end if

                end do

                rho_ave = 0.5d0*(rho_i+rho_ip)
                u_ave   = 0.5d0*(u_i+u_ip)
                T_ave   = 0.5d0*(T_i+T_ip)
                mu_ave  = 0.5d0*(mu_i+mu_ip)
                kap_ave = 0.5d0*(kap_i+kap_ip)
                tau_xx  = -2.0d0/3.0d0*mu_ave*(2.0d0*((u_ip-u_i)/b_obj(ib)%smr_dx))

                call ckhms(T_ave,ickwrk,rckwrk,hs) ! enthalpy
                hs = hs*1.0d-4
                hsum = 0.0d0

                do k_sp = 1, nspe
                    dif_ave = 0.5d0*(dif_i(k_sp)+dif_ip(k_sp))
                    hsum = hsum - hs(k_sp)*dif_ave*((y_sp_ip(k_sp)-y_sp_i(k_sp))/b_obj(ib)%smr_dx)
                end do

                b_obj(ib)%Fplus(iudens,i) = b_obj(ib)%Fplus(iudens,i) + 0.0d0
                b_obj(ib)%Fplus(iumomx,i) = b_obj(ib)%Fplus(iumomx,i) + tau_xx
                b_obj(ib)%Fplus(iuener,i) = b_obj(ib)%Fplus(iuener,i) + &
                          tau_xx*u_ave + rho_ave*hsum - kap_ave*((T_ip-T_i)/b_obj(ib)%smr_dx)

                do k_sp = 1, nspe

                    dif_ave = 0.5d0*(dif_i(k_sp)+dif_ip(k_sp))

                    b_obj(ib)%Fplus(ncons+k_sp,i) = b_obj(ib)%Fplus(ncons+k_sp,i) - &
                              rho_ave*dif_ave*((y_sp_ip(k_sp)-y_sp_i(k_sp))/b_obj(ib)%smr_dx)

                end do

            end do
        end do

    end subroutine calc_flux_tran


    subroutine calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)

        use module_params, only: dx_comp, flag_acu, flag_muscl, flag_weno, is_spherical

        implicit none


        integer :: ib, i, n, k_sp

        real(8) :: uL, uR       ! velocity
        real(8) :: rhoL, rhoR   ! density
        real(8) :: tempL, tempR ! temperature
        real(8) :: presL, presR ! pressure
        real(8) :: hL, hR       ! enthalpy
        real(8) :: eL, eR       ! total energy
        real(8) :: sL, sR, sM   ! wave speed
        real(8) :: aL, aR       ! speed of sound
        real(8) :: x_spL(nspe), x_spR(nspe) ! mole fraction
        real(8) :: fL(nequ), fR(nequ)    ! intercell flux
        real(8) :: usL(nequ), usR(nequ)
        real(8) :: vec_hllc(nequ)

        real(8) :: lambdaL, lambdaR
        real(8) :: alpha
        
        real(8) :: rt   ! Roe Average
        real(8) :: r_roe, u_roe, h_roe
        real(8) :: a_roe, p_roe

        real(8) :: qf_muscl(4)
        real(8) :: qf_weno3(4)
        real(8) :: qf_weno5(6)
        real(8) :: qf_weno7(8)

        real(8) :: qfL(nequ), qfR(nequ)
        real(8) :: QcL(nequ), QcR(nequ)

        real(8) :: flux_advec(nequ)

        integer,intent(in) :: lvir, nspe, nequ
        integer,intent(in) :: select_riemann

        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ib = bsta, bend

            do i = b_obj(ib)%ista_b-1, b_obj(ib)%iend_b

                if ( flag_acu == 1 ) then

                    do n = 1, nequ

                        qf_muscl(1) = b_obj(ib)%qf(n,i-1)
                        qf_muscl(2) = b_obj(ib)%qf(n,i)
                        qf_muscl(3) = b_obj(ib)%qf(n,i+1)
                        qf_muscl(4) = b_obj(ib)%qf(n,i+2)

                        call muscl3(qf_muscl, qfL(n), qfR(n), flag_muscl)

                    end do

                else if ( flag_acu == 2 .and. flag_weno == 1 ) then

                    do n = 1, nequ

                        qf_weno3(1) = b_obj(ib)%qf(n,i-1)
                        qf_weno3(2) = b_obj(ib)%qf(n,i)
                        qf_weno3(3) = b_obj(ib)%qf(n,i+1)
                        qf_weno3(4) = b_obj(ib)%qf(n,i+2)

                        call weno3z(qf_weno3, qfL(n), qfR(n))

                    end do

                else if ( flag_acu == 2 .and. (flag_weno == 2 .or. flag_weno == 3) ) then

                    do n = 1, nequ

                        qf_weno5(1) = b_obj(ib)%qf(n,i-2)
                        qf_weno5(2) = b_obj(ib)%qf(n,i-1)
                        qf_weno5(3) = b_obj(ib)%qf(n,i)
                        qf_weno5(4) = b_obj(ib)%qf(n,i+1)
                        qf_weno5(5) = b_obj(ib)%qf(n,i+2)
                        qf_weno5(6) = b_obj(ib)%qf(n,i+3)

                        call weno5(qf_weno5, qfL(n), qfR(n), flag_weno)

                    end do

                else if ( flag_acu == 2 .and. flag_weno == 4 ) then

                    do n = 1, nequ

                        qf_weno7(1) = b_obj(ib)%qf(n,i-3)
                        qf_weno7(2) = b_obj(ib)%qf(n,i-2)
                        qf_weno7(3) = b_obj(ib)%qf(n,i-1)
                        qf_weno7(4) = b_obj(ib)%qf(n,i)
                        qf_weno7(5) = b_obj(ib)%qf(n,i+1)
                        qf_weno7(6) = b_obj(ib)%qf(n,i+2)
                        qf_weno7(7) = b_obj(ib)%qf(n,i+3)
                        qf_weno7(8) = b_obj(ib)%qf(n,i+3)

                        call weno7z(qf_weno7, qfL(n), qfR(n))

                    end do

                end if

                if ( flag_acu == 0 .or. dx_comp /= b_obj(ib)%smr_dx ) then

                    qfL = b_obj(ib)%qf(:,i)
                    qfR = b_obj(ib)%qf(:,i+1)

                end if

                QcL = qfToQc(qfL, nspe, nequ)
                QcR = qfToQc(qfR, nspe, nequ)

                call Yi2Xi(nspe, qfL(nprim+1:nequ), x_spl)
                call Yi2Xi(nspe, qfR(nprim+1:nequ), x_spr)
                call adjust_sum(nspe, x_spl)
                call adjust_sum(nspe, x_spr)

                rhoL = qfL(iqdens)
                rhoR = qfR(iqdens)
                uL = qfL(iqxvel)
                uR = qfR(iqxvel)
                tempL = qfL(iqtemp)
                tempR = qfR(iqtemp)
                presL = pressure(rhoL, tempL, nspe, x_spL)
                presR = pressure(rhoR, tempR, nspe, x_spR)
                aL = sonic_speed(rhoL, presL, tempL, nspe, x_spL)
                aR = sonic_speed(rhoR, presR, tempR, nspe, x_spR)
                eL = QcL(iuener)
                eR = QcR(iuener)
                hL = (eL + presL) / rhoL
                hR = (eR + presR) / rhoR

                fL(iudens) = rhoL*uL ! rho*u
                fR(iudens) = rhoR*uR ! rho*u
                if (is_spherical == 1) then
                    fL(iumomx) = rhoL*uL*uL ! rho*u*u
                    fR(iumomx) = rhoR*uR*uR ! rho*u*u
                else 
                    fL(iumomx) = rhoL*uL*uL + presL ! rho*u*u + p
                    fR(iumomx) = rhoR*uR*uR + presR ! rho*u*u + p
                end if
                fL(iuener) = (eL + presL)*uL ! (e+p)*u
                fR(iuener) = (eR + presR)*uR ! (e+p)*u
                do k_sp = 1, nspe
                    fL(ncons+k_sp) = rhoL*qfL(nprim+k_sp)*uL ! rho*Yi*u
                    fR(ncons+k_sp) = rhoR*qfR(nprim+k_sp)*uR ! rho*Yi*u
                end do


                if ( select_riemann == 0 ) then

                    ! Central Scheme (LxF)

                    lambdaL = dmax1(dabs(uL-aL), dabs(uL), dabs(uL+aL))
                    lambdaR = dmax1(dabs(uR-aR), dabs(uR), dabs(uR+aR))
                    alpha = dmax1(lambdaL, lambdaR)
                    flux_advec = 0.5d0 * (fL + fR - alpha*(QcR-QcL))

                else if ( select_riemann == 1 ) then

                    ! HLL flux

                    rt = dsqrt(rhoR / rhoL) ! Roe Average
                    r_roe = dsqrt(rhoR * rhoL)
                    u_roe = (uL + rt * uR) / (1.d0 + rt)
                    h_roe = (hL + rt * hR) / (1.d0 + rt)
                    a_roe = (aL + rt * aR) / (1.d0 + rt)

                    sl = dmin1(ul - al, u_roe - a_roe)
                    sr = dmax1(ur + ar, u_roe + a_roe)

                    if (0.d0 <= sL) then
                        flux_advec = fL
                    else if (sL <= 0.d0 .and. 0.d0 <= sR) then
                        flux_advec = (sR*fL-sL*fR+sL*sR*(QcR-QcL))/(sR-sL)
                    else if (sR <= 0.d0) then
                        flux_advec = fR
                    end if

                else if ( select_riemann == 2 ) then

                    ! HLLC flux

                    rt = dsqrt(rhoR / rhoL) ! Roe Average
                    r_roe = dsqrt(rhoR * rhoL)
                    u_roe = (uL + rt * uR) / (1.d0 + rt)
                    h_roe = (hL + rt * hR) / (1.d0 + rt)
                    a_roe = (aL + rt * aR) / (1.d0 + rt)

                    sL = dmin1(uL - aL, u_roe - a_roe)
                    sR = dmax1(uR + aR, u_roe + a_roe)
                    sM = (presL - presR + rhoR * uR * (sR - uR) - rhoL * uL * (sL - uL)) &
                                            / (rhoR * (sR - uR) - rhoL * (sL - uL))

                    if (0.d0 <= sL) then
                        flux_advec = fL
                    else if (sL <= 0.d0 .and. 0.d0 <= sM) then
                        vec_hllc(iudens) = 1.d0
                        vec_hllc(iumomx) = sM
                        vec_hllc(iuener) = eL/rhoL + (sM-uL)*(sM+presL/(rhoL*(sL-uL)))
                        do k_sp = 1, nspe
                            vec_hllc(ncons+k_sp) = QcL(ncons+k_sp)/rhoL
                        end do
                        usL = rhoL * (sL-uL)/(sL-sM) * vec_hllc
                        flux_advec = fL + sL * (usL - QcL)
                    else if (sM <= 0.d0 .and. 0.d0 <= sR) then
                        vec_hllc(iudens) = 1.d0
                        vec_hllc(iumomx) = sM
                        vec_hllc(iuener) = eR/rhoR + (sM-uR)*(sM+presR/(rhoR*(sR-uR)))
                        do k_sp = 1, nspe
                            vec_hllc(ncons+k_sp) = QcR(ncons+k_sp)/rhoR
                        end do
                        usR = rhoR * (sR-uR)/(sR-sM) * vec_hllc
                        flux_advec = fR + sR * (usR - QcR)
                    else if (sR <= 0.d0) then
                        flux_advec = fR
                    end if

                end if

                b_obj(ib)%Fplus(:,i) = b_obj(ib)%Fplus(:,i) + flux_advec

                rt = dsqrt(rhoR / rhoL) ! Roe Average
                r_roe = dsqrt(rhoR * rhoL)
                u_roe = (uL + rt * uR) / (1.d0 + rt)
                p_roe = (presL + rt * presR) / (1.d0 + rt)

                b_obj(ib)%godunov_state(iqdens,i) = r_roe
                b_obj(ib)%godunov_state(iqxvel,i) = u_roe
                b_obj(ib)%godunov_state(iqpres,i) = p_roe

            end do

        end do

    end subroutine calc_flux_advection


    ! subroutine calc_flux_central(qf, Qc, dt, nx_0, nx_n, lvir, lnx_0, lnx_n, dx, nspe, nequ)
    !         Fplus(:,j) = 0.5d0*(FL+FR-cwave*(QcR(:,j)-QcL(:,j)))
    ! end subroutine calc_flux_central

    ! subroutine calc_flux_Lax(qf, Qc, nx_0, nx_n, lvir, dx, dt, lnx_0, lnx_n, nspe, nequ, Fplus)
    !         Fplus(:,j) = 0.5d0*((FL+FR)-(dx/dt)*(QcR(:,j)-QcL(:,j)))
    ! end subroutine calc_flux_Lax


end module flux


