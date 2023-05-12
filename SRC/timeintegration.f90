module timeintegration

    use flux
    use parallel
    use block_t, only:block_st

    implicit none

contains


    subroutine explicit_euler(bsta, bend, b_obj, lvir, dt, nspe, nequ, &
                              block_all_0, block_all_n)

        use module_params, only: do_hydro, do_tran, select_riemann, is_spherical

        implicit none


        integer :: ib, i
        real(8) :: dx
        real(8) :: rm, Al, Ar, dV

        integer,intent(in) :: nspe, nequ, lvir
        real(8),intent(in) :: dt

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        do ib = bsta, bend
            b_obj(ib)%Fplus(:,:) = 0.d0
        end do

        if ( do_tran == 1 ) then

            call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)

        end if

        if ( do_hydro == 1 ) then

            call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)

        end if

        if (is_spherical == 0) then

            do ib = bsta, bend
                dx = b_obj(ib)%smr_dx
                do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                    b_obj(ib)%Qc(:,i) = b_obj(ib)%Qc(:,i) - & 
                              (dt/dx) * (b_obj(ib)%Fplus(:,i) - b_obj(ib)%Fplus(:,i-1))
                end do
            end do

        else if (is_spherical == 1) then

            do ib = bsta, bend
                dx = b_obj(ib)%smr_dx
                do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                    rm = b_obj(ib)%xm_b(i)
                    Al = (rm-(dx/2.d0))**2
                    Ar = (rm+(dx/2.d0))**2
                    dV = rm**2 * dx

                    b_obj(ib)%Qc(:,i) = b_obj(ib)%Qc(:,i) - & 
                              (dt/dV) * (Ar*b_obj(ib)%Fplus(:,i) - Al*b_obj(ib)%Fplus(:,i-1))

                    b_obj(ib)%Qc(iumomx,i) = b_obj(ib)%Qc(iumomx,i) - & 
                              (dt/dx) * (b_obj(ib)%godunov_state(iqpres,i) - b_obj(ib)%godunov_state(iqpres,i-1))
                    
                end do
            end do

        end if

        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)

        call interpolation_block(nequ, lvir)

        if (nprocs /= 1) call communicate_parallel_block(nequ, lvir)

        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b-lvir, b_obj(ib)%iend_b+lvir
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do

    end subroutine explicit_euler


    subroutine TVDRK3(ista_rank, iend_rank, bsta, bend, b_obj, &
                      lvir, dt, nspe, nequ, block_all_0, block_all_n)

        use module_params, only: do_hydro, do_tran, select_riemann, is_spherical

        implicit none


        integer :: istage, ib, i
        real(8) :: dx
        real(8) :: rm, Al, Ar, dV

        real(8) :: h1, h2
        real(8),allocatable :: Qc_RK_0(:,:)

        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: dt

        integer,intent(in) :: ista_rank, iend_rank, lvir
        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        allocate ( Qc_RK_0(nequ,ista_rank:iend_rank) )

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                Qc_RK_0(:,i) = b_obj(ib)%Qc(:,i)
            end do
        end do


        do istage = 1, 3

            do ib = bsta, bend
                b_obj(ib)%Fplus(:,:) = 0.d0
            end do

            if ( do_tran == 1 ) then

                call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)

            end if

            if ( do_hydro == 1 ) then

                call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)

            end if

            if ( istage == 1 ) then
                h1  = 1.0d0
                h2  = 1.0d0
                do ib = bsta, bend
                    b_obj(ib)%Qc(:,:) = 0.d0
                end do
            else if ( istage == 2 ) then
                h1  = 3.d0/4.d0
                h2  = 1.d0/4.d0
            else if ( istage == 3 ) then
                h1  = 1.0d0/3.0d0
                h2  = 2.0d0/3.0d0
            end if


            if (is_spherical == 0) then

                do ib = bsta, bend
                    dx = b_obj(ib)%smr_dx

                    do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                        b_obj(ib)%Qc(:,i) = h1*Qc_RK_0(:,i) + h2*(b_obj(ib)%Qc(:,i) - &
                                  (dt/dx) * (b_obj(ib)%Fplus(:,i) - b_obj(ib)%Fplus(:,i-1)))
                    end do
                end do

            else if (is_spherical == 1) then

                do ib = bsta, bend
                    dx = b_obj(ib)%smr_dx

                    do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b

                        rm = b_obj(ib)%xm_b(i)
                        Al = (rm-(dx/2.d0))**2
                        Ar = (rm+(dx/2.d0))**2
                        dV = rm**2 * dx

                        b_obj(ib)%Qc(:,i) = h1*Qc_RK_0(:,i) + h2*(b_obj(ib)%Qc(:,i) - &
                                  (dt/dV) * (Ar*b_obj(ib)%Fplus(:,i) - Al*b_obj(ib)%Fplus(:,i-1)))

                        b_obj(ib)%Qc(iumomx,i) = b_obj(ib)%Qc(iumomx,i) - &
                                  (dt/dx) * (b_obj(ib)%godunov_state(iqpres,i) - b_obj(ib)%godunov_state(iqpres,i-1))

                    end do
                end do

            end if


            call update_boundary(bsta, bend, b_obj, lvir, &
                                 nspe, nequ, block_all_0, block_all_n)

            call interpolation_block(nequ, lvir)

            if (nprocs /= 1) call communicate_parallel_block(nequ, lvir)

            call update_boundary(bsta, bend, b_obj, lvir, &
                                 nspe, nequ, block_all_0, block_all_n)

            do ib = bsta, bend
                do i = b_obj(ib)%ista_b-lvir, b_obj(ib)%iend_b+lvir
                    b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
                end do
            end do

        end do

    end subroutine TVDRK3


    subroutine RK4(ista_rank, iend_rank, bsta, bend, b_obj, &
                   lvir, dt, nspe, nequ, block_all_0, block_all_n)

        use module_params, only: do_hydro, do_tran, select_riemann, is_spherical

        implicit none


        integer :: i, ib
        real(8),allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
        real(8),allocatable :: Qc_RK_0(:,:)

        integer,intent(in) :: nspe, nequ
        real(8),intent(in) :: dt

        integer,intent(in) :: block_all_0, block_all_n
        integer,intent(in) :: bsta, bend
        integer,intent(in) :: ista_rank, iend_rank, lvir
        type(block_st),intent(inout) :: b_obj(bsta:bend)


        allocate( k1(nequ,ista_rank:iend_rank) )
        allocate( k2(nequ,ista_rank:iend_rank) )
        allocate( k3(nequ,ista_rank:iend_rank) )
        allocate( k4(nequ,ista_rank:iend_rank) )

        allocate ( Qc_RK_0(nequ,ista_rank:iend_rank) )

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                Qc_RK_0(:,i) = b_obj(ib)%Qc(:,i)
            end do
        end do

        do ib = bsta, bend
            b_obj(ib)%Fplus(:,:) = 0.d0
        end do

        if ( do_tran == 1 ) then
            call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)
        end if
        if ( do_hydro == 1 ) then
            call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)
        end if

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                k1(:,i) = - (dt/b_obj(ib)%smr_dx)*(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
                b_obj(ib)%Qc(:,i) = b_obj(ib)%Qc(:,i) - 0.5d0*(0.5d0*dt/b_obj(ib)%smr_dx)&
                                      *(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
            end do
        end do
        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)


        do ib = bsta, bend
            b_obj(ib)%Fplus(:,:) = 0.d0
        end do

        if ( do_tran == 1 ) then
            call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)
        end if
        if ( do_hydro == 1 ) then
            call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)
        end if
        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                k2(:,i) = - (dt/b_obj(ib)%smr_dx)*(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
                b_obj(ib)%Qc(:,i) = b_obj(ib)%Qc(:,i) - 0.5d0*(0.5d0*dt/b_obj(ib)%smr_dx)&
                                      *(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
            end do
        end do
        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)


        do ib = bsta, bend
            b_obj(ib)%Fplus(:,:) = 0.d0
        end do

        if( do_tran == 1 ) then
            call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)
        end if
        if ( do_hydro == 1 ) then
            call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)
        end if
        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                k3(:,i) = - (dt/b_obj(ib)%smr_dx)*(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
                b_obj(ib)%Qc(:,i) = b_obj(ib)%Qc(:,i) - (dt/b_obj(ib)%smr_dx)&
                                      *(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
            end do
        end do
        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)


        do ib = bsta, bend
            b_obj(ib)%Fplus(:,:) = 0.d0
        end do

        if( do_tran == 1 ) then
            call calc_flux_tran(bsta, bend, b_obj, lvir, nspe, nequ)
        end if
        if ( do_hydro == 1 ) then
            call calc_flux_advection(bsta, bend, b_obj, lvir, nspe, nequ, select_riemann)
        end if

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                k4(:,i) = - (dt/b_obj(ib)%smr_dx)*(b_obj(ib)%Fplus(:,i)-b_obj(ib)%Fplus(:,i-1))
            end do
        end do

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b, b_obj(ib)%iend_b
                b_obj(ib)%Qc(:,i) = Qc_RK_0(:,i) + (1.d0/6.d0)*(k1(:,i)+2.d0*k2(:,i)+2.d0*k3(:,i)+k4(:,i))
            end do
        end do

        call update_boundary(bsta, bend, b_obj, lvir, &
                             nspe, nequ, block_all_0, block_all_n)

        deallocate(k1, k2, k3, k4)
        deallocate(Qc_RK_0)

        do ib = bsta, bend
            do i = b_obj(ib)%ista_b-lvir, b_obj(ib)%iend_b+lvir
                b_obj(ib)%qf(:,i) = QcToqf(b_obj(ib)%Qc(:,i), nspe, nequ)
            end do
        end do


    end subroutine RK4


end module timeintegration

