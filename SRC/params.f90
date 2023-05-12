module module_params
    
    use thermo

    implicit none


    !!! solver flags
    integer,save :: do_hydro
    integer,save :: do_tran
    integer,save :: do_reaction
    integer,save :: select_time
    integer,save :: select_riemann
    integer,save :: flag_cfl
    integer,save :: flag_acu ! 0:1st order, 1:higher order
    integer,save :: flag_muscl
    integer,save :: flag_weno
    real(8) :: atol, rtol

    character(len=256),save :: chem_dir
    character(len=128),save :: prob_file
    character(len=128),save :: init_file
    character(len=128),save :: flow_bound_file
    character(len=128),save :: radiation_file
    character(len=128),save :: wallT_bound_file
    character(len=128),save :: igsource_file

    real(8),save :: dt            ! time step size [s]
    real(8),save :: tend          ! end time [s]
    real(8),save :: cfl           ! Courant number
    real(8),save :: out_timestep  ! time interval of output file
    integer,save :: out_num_hdf   ! number of step in output hdf file
    real(8),save :: dx_comp       ! space size [m]
    real(8),save :: x_length_comp ! domain length [m]
    integer,save :: nx_comp       ! number of space step
    integer,save :: lvir          ! number of virtual boundary cell
    integer :: is_spherical = 0   ! 0:cartesian, 1:spherical

    integer,save :: flag_block
    integer,save :: ncell_per_block
    integer,save :: lev_smr_0, num_smr_0
    integer,save :: lev_smr_n, num_smr_n

    integer,save :: flag_lbound, flag_rbound
    integer,save :: flag_radiation
    integer,save :: flag_Tw
    integer,save :: flag_Tini_Tw
    integer,save :: flag_igsource
    integer,save :: flag_const_rhou

    integer,save :: flag_check_init
    integer,save :: flag_restart
    character(len=128),save :: readdir    ! directory where restart is stored 
    character(len=128),save :: readfile   ! read file to restart simulation
    character(len=128),save :: readpath   ! read path to restart simulation
    character(len=128),save :: savedir    ! directory to save hdf5 results

    integer,save :: nspe  ! number of species
    integer,save :: nequ  ! number of equations

    real(8) :: lVel , rVel
    real(8) :: lTemp, rTemp
    real(8) :: lPres, rPres
    character(len=128) :: sp_name
    real(8) :: sp_mole
    integer :: lk_sp, rk_sp
    integer :: lnsp, rnsp
    real(8),allocatable :: lXi(:), rXi(:)

    integer :: flag_errorfunc_T
    real(8) :: center_efunc_T, radius_efunc_T
    real(8) :: smoother_efunc_T
    integer :: flag_errorfunc_p
    real(8) :: center_efunc_p, radius_efunc_p
    real(8) :: smoother_efunc_p
    integer :: flag_errorfunc_v
    real(8) :: center_efunc_v, radius_efunc_v
    real(8) :: smoother_efunc_v

    integer :: flag_errorfunc_sp
    real(8) :: center_efunc_sp, radius_efunc_sp
    real(8) :: smoother_efunc_sp
    character(len=128) :: sp_name_other
    integer :: flag_equilXi
    character(len=128),save :: readpath_equilXi ! read equil Xi path

    type :: eos_st

        real(8) :: dx, x_m
        real(8) :: rho, u, T, p
        real(8) :: eTotal, eKinetic, eInternal
        real(8) :: HRR, HRR_dot
        real(8) :: Tw, qTw
        real(8) :: Qrad ! Radiation heat loss 

        real(8),allocatable :: x_sp(:), y_sp(:)
        real(8),allocatable :: y_sp0(:)

    end type eos_st

    type(eos_st),allocatable :: eos_obj(:)


    ! namelist solver
    namelist /params/ do_hydro, do_tran, do_reaction, &
                      select_time, select_riemann, &
                      flag_cfl, flag_acu, &
                      flag_muscl, flag_weno, &
                      atol, rtol, chem_dir, &
                      prob_file, init_file

    ! namelist time and spatial size
    namelist /problem/ dt, tend, &
                       out_timestep, out_num_hdf, &
                       cfl, dx_comp, x_length_comp, &
                       is_spherical

    ! namelist blocks
    namelist /blocks/ flag_block, ncell_per_block, &
                      lev_smr_0, lev_smr_n, &
                      num_smr_0, num_smr_n

    ! namelist boundary
    namelist /boundary/ flag_lbound, flag_rbound, &
                        flow_bound_file

    ! namelist radiation
    namelist /radiationloss/ flag_radiation, radiation_file

    ! namelist mfr
    namelist /wallMFR/ flag_Tw, flag_Tini_Tw, &
                       wallT_bound_file

    ! namelist energy source
    namelist /energysource/ flag_igsource, igsource_file

    ! namelist input and output
    namelist /params_io/ flag_check_init, &
                         flag_restart, readfile, &
                         readdir, savedir

contains


    subroutine init_solvers()

        implicit none

        integer :: lun
        character(len=128) :: infile
        

        infile = 'solver.inp'
        
        open(newunit=lun, file=trim(infile), status="old", action="read")
        read(unit=lun, nml=params)
        close(unit=lun)

        if ( flag_acu == 0 ) then
            lvir = 1
        else if ( flag_acu == 1 ) then
            lvir = 2
        else if (flag_acu == 2 .and. flag_weno == 1) then
            lvir = 2
        else if (flag_acu == 2 .and. flag_weno == 2) then
            lvir = 3
        else if (flag_acu == 2 .and. flag_weno == 3) then
            lvir = 3
        else if (flag_acu == 2 .and. flag_weno == 4) then
            lvir = 4
        end if

    end subroutine init_solvers


    subroutine read_inputfiles(myrank)
        
        use mymkdir
        
        implicit none

    
        integer,intent(in) :: myrank
        integer :: k_sp

        integer :: lun1, lun2, lun3, lun4, lun5, lun6, lun7


        open(newunit=lun1, file=trim(prob_file), status="old", action="read")
        read(unit=lun1, nml=problem)
        close(unit=lun1)

        open(newunit=lun2, file=trim(prob_file), status="old", action="read")
        read(unit=lun2, nml=blocks)
        close(unit=lun2)

        open(newunit=lun3, file=trim(prob_file), status="old", action="read")
        read(unit=lun3, nml=boundary)
        close(unit=lun3)

        open(newunit=lun4, file=trim(prob_file), status="old", action="read")
        read(unit=lun4, nml=radiationloss)
        close(unit=lun4)

        open(newunit=lun5, file=trim(prob_file), status="old", action="read")
        read(unit=lun5, nml=wallMFR)
        close(unit=lun5)

        open(newunit=lun6, file=trim(prob_file), status="old", action="read")
        read(unit=lun6, nml=energysource)
        close(unit=lun6)

        open(newunit=lun7, file=trim(prob_file), status="old", action="read")
        read(unit=lun7, nml=params_io)
        close(unit=lun7)


        ! make directory
        call makedirs(savedir)
        write (readpath,'("./",a,"/",a)') trim(readdir), trim(readfile)


        if (lev_smr_0 < 1) then
            lev_smr_0 = 1
        else if (lev_smr_n < 1) then
            lev_smr_n = 1
        end if

        nx_comp = nint(x_length_comp/dx_comp)
        nspe = n_species 
        nequ = 3+nspe

        open(18,file=trim(init_file), status='old')

        ! Read initial conditions
        read(18,*) lTemp, rTemp
        read(18,*) flag_errorfunc_T
        read(18,*) center_efunc_T
        read(18,*) smoother_efunc_T
        read(18,*) radius_efunc_T

        read(18,*) lPres, rPres
        read(18,*) flag_errorfunc_p
        read(18,*) center_efunc_p
        read(18,*) smoother_efunc_p
        read(18,*) radius_efunc_p
        lPres = lPres*101325.d0 ! atm  -> Pa
        rPres = rPres*101325.d0 ! atm  -> Pa

        read(18,*) flag_const_rhou
        read(18,*) lVel , rVel
        read(18,*) flag_errorfunc_v
        read(18,*) center_efunc_v
        read(18,*) smoother_efunc_v
        read(18,*) radius_efunc_v
        lVel = lVel*1.0d-2   ! cm/s -> m/s
        rVel = rVel*1.0d-2   ! cm/s -> m/s

        allocate(lXi(nspe), rXi(nspe))

        read(18,*) lnsp

        lXi = 0.d0
        do lk_sp = 1, lnsp
            read(18,*) sp_name, sp_mole
            ! print *, sp_name, sp_mole
            do k_sp = 1, nspe
                if (trim(sp_name) == species_name(k_sp)) then
                    lXi(k_sp) = sp_mole
                end if
            end do
        end do
        call adjust_sum(nspe, lXi)

        read(18,*) rnsp

        rXi = 0.d0
        do rk_sp = 1, rnsp
            read(18,*) sp_name, sp_mole
            ! print *, sp_name, sp_mole
            do k_sp = 1, nspe
                if (trim(sp_name) == species_name(k_sp)) then
                    rXi(k_sp) = sp_mole
                end if
            end do
        end do
        call adjust_sum(nspe, rXi)

        read(18,*) flag_errorfunc_sp
        read(18,*) center_efunc_sp
        read(18,*) smoother_efunc_sp
        read(18,*) radius_efunc_sp

        read(18,*) flag_equilXi
        read(18,*) sp_name_other
        read(18,*) readpath_equilXi

        close(18)

    end subroutine read_inputfiles


    subroutine initialize_variables(myrank, ista_rank, iend_rank)

        use mfr, only: Tw
        use ignitionsource
        use flow_boundary

        implicit none


        integer,intent(in) :: myrank
        integer,intent(in) :: ista_rank, iend_rank
        integer :: i, k_sp

        integer :: k_sp_equil, nspe_equil
        real(8) :: Xi_low = 0.d0
        real(8) :: Xi_sp_rest
        real(8),allocatable :: Xi_equil(:)

        real(8) :: x_m        ! location [m]
        real(8) :: rho_u_FILL ! initial constant rhou [kg/m^2-s]
        real(8) :: coef, Eadd ! initial energy addition [J/m^3]
        real(8) :: eInternal_vol ! internal energy [J/m^3]


        if ( myrank == 0 ) then
            print *, 'number of grid points', nx_comp
            print *, 'number of species  ', nspe
            print *, 'number of equations', nequ
        end if

        if (flag_equilXi == 1) then

            allocate(Xi_equil(nspe))

            open(21,file=trim(readpath_equilXi), status='old')

            read(21,*) nspe_equil

            do k_sp_equil = 1, nspe_equil
                read(21,*) sp_name, sp_mole
                do k_sp = 1, nspe
                    if (trim(sp_name) == species_name(k_sp)) then
                        Xi_equil(k_sp) = sp_mole
                    end if
                end do
            end do

            close(21)

        end if

        do i = ista_rank, iend_rank

            allocate(eos_obj(i)%x_sp(nspe))
            allocate(eos_obj(i)%y_sp(nspe))
            allocate(eos_obj(i)%y_sp0(nspe))

            eos_obj(i)%x_sp(:) = 0.0d0

            x_m = eos_obj(i)%x_m

            if (flag_errorfunc_T == 0) then

                if ( x_m <= center_efunc_T ) then
                    eos_obj(i)%T = lTemp
                else
                    eos_obj(i)%T = rTemp
                end if

            else if (flag_errorfunc_T == 1) then

                eos_obj(i)%T = lTemp + 0.5d0*(rTemp-lTemp) &
                *(erf((x_m-center_efunc_T)/smoother_efunc_T)+1.0d0)

            else if (flag_errorfunc_T == 2) then

                eos_obj(i)%T = lTemp + 0.5d0*(rTemp-lTemp) &
                *(1.d0 - erf((x_m-center_efunc_T)/smoother_efunc_T))

            else if (flag_errorfunc_T == 3) then

                eos_obj(i)%T = lTemp + 0.5d0*(rTemp-lTemp) &
                *(1.d0 - erf((dsqrt((x_m-center_efunc_T)**2)-radius_efunc_T) &
                                                          /smoother_efunc_T))

            else if (flag_errorfunc_T == 4) then ! X.Chen et al. PCI (2021)

                eos_obj(i)%T = (lTemp - rTemp) * dexp(-(x_m/radius_efunc_T)**2) + rTemp

            else if (flag_errorfunc_T == 5) then ! M. Tao et al. Frontiers in Mechanical Engineering (2020)

                eos_obj(i)%T = lTemp
                if (x_m <= 0.03d0) then
                    eos_obj(i)%T = lTemp + (rTemp - lTemp)/0.03d0 * x_m
                else
                    eos_obj(i)%T = rTemp
                end if

            end if


            if (flag_errorfunc_p == 0) then

                if ( x_m <= center_efunc_p ) then
                    eos_obj(i)%p = lPres
                else
                    eos_obj(i)%p = rPres
                end if

            else if (flag_errorfunc_p == 1) then

                eos_obj(i)%p = lPres + 0.5d0*(rPres-lPres) &
                *(erf((x_m-center_efunc_p)/smoother_efunc_p)+1.0d0)

            else if (flag_errorfunc_p == 2) then

                eos_obj(i)%p = lPres + 0.5d0*(rPres-lPres) &
                *(1.d0 - erf((x_m-center_efunc_p)/smoother_efunc_p))

            else if (flag_errorfunc_p == 3) then

                eos_obj(i)%p = lPres + 0.5d0*(rPres-lPres) &
                *(1.d0 - erf((dsqrt((x_m-center_efunc_p)**2)-radius_efunc_p) &
                                                          /smoother_efunc_p))
            end if


            if (flag_errorfunc_sp == 0) then

                if ( x_m <= center_efunc_sp ) then
                    eos_obj(i)%x_sp = lXi
                else
                    eos_obj(i)%x_sp = rXi
                end if

                call adjust_sum(nspe, eos_obj(i)%x_sp)

            else if (flag_errorfunc_sp == 1) then

                do k_sp = 1, nspe
                    eos_obj(i)%x_sp(k_sp) = Xi_low + 0.5d0*(lXi(k_sp)-Xi_low) &
                    *(erf((x_m-center_efunc_sp)/smoother_efunc_sp)+1.0d0)
                end do

            else if (flag_errorfunc_sp == 2) then

                do k_sp = 1, nspe
                    eos_obj(i)%x_sp(k_sp) = Xi_low + 0.5d0*(lXi(k_sp)-Xi_low) &
                    *(1.d0 - erf((x_m-center_efunc_sp)/smoother_efunc_sp))
                end do

            else if (flag_errorfunc_sp == 3) then

                do k_sp = 1, nspe
                    eos_obj(i)%x_sp(k_sp) = Xi_low + 0.5d0*(lXi(k_sp)-Xi_low) &
                    *(1.d0 - erf((dsqrt((x_m-center_efunc_sp)**2)-radius_efunc_sp) &
                                                               /smoother_efunc_sp))
                end do

            end if

            if (flag_equilXi == 0) then

                do k_sp = 1, nspe
                    if (trim(sp_name_other) == species_name(k_sp)) then
                        eos_obj(i)%x_sp(k_sp) = eos_obj(i)%x_sp(k_sp) + 1.d0 - sum(eos_obj(i)%x_sp(:))
                    end if
                end do

            else if (flag_equilXi == 1) then

                Xi_sp_rest = 1.d0 - sum(eos_obj(i)%x_sp(:))
                do k_sp = 1, nspe
                    eos_obj(i)%x_sp = eos_obj(i)%x_sp + Xi_sp_rest*Xi_equil
                end do

            end if
 
            call adjust_sum(nspe, eos_obj(i)%x_sp)
            call Xi2Yi(nspe, eos_obj(i)%x_sp(:), eos_obj(i)%y_sp(:))
            call Xi2Yi(nspe, eos_obj(i)%x_sp(:), eos_obj(i)%y_sp0(:))


            if ( flag_Tini_Tw == 1 ) then

                eos_obj(i)%p = p_OUT
                eos_obj(i)%T = Tw(x_m)  ! Initial temperature == mfr wall temp profile

            end if

            if ( flag_Tw == 1 ) then

                eos_obj(i)%Tw = Tw(x_m) ! Set mfr wall temperature profile

            end if

            if (tig_0_s < 1.0d-9 .and. tig_s < 1.0d-9 .and. flag_igsource /= 0) then

                if (flag_igsource == 1) then

                    coef = Eig_j / volumeig
                    Eadd = coef * dexp(-(pi / 4.d0) * (x_m / rig_m)**6)

                else if (flag_igsource == 2) then

                    coef = Eig_j / (pi**1.5 * rig_m**3 * tig_s)
                    Eadd = coef * dexp(-(x_m/rig_m)**2)

                end if

                eos_obj(i)%rho = density(eos_obj(i)%p, eos_obj(i)%T, nspe, eos_obj(i)%x_sp)
                eos_obj(i)%eInternal = specific_energy(eos_obj(i)%rho, eos_obj(i)%p, &
                                                       eos_obj(i)%T, nspe, eos_obj(i)%x_sp)
                eInternal_vol = eos_obj(i)%rho*eos_obj(i)%eInternal + Eadd
                eos_obj(i)%eInternal = eInternal_vol/eos_obj(i)%rho
                eos_obj(i)%T = temperature_con(eos_obj(i)%eInternal, nspe, eos_obj(i)%x_sp)

            end if


            eos_obj(i)%rho = density(eos_obj(i)%p, eos_obj(i)%T, nspe, eos_obj(i)%x_sp)

            if ( flag_const_rhou == 0 ) then

                if (flag_errorfunc_v == 0) then

                    if ( x_m <= center_efunc_v ) then
                        eos_obj(i)%u = lVel
                    else
                        eos_obj(i)%u = rVel
                    end if

                else if (flag_errorfunc_v == 1) then

                    eos_obj(i)%u = lVel + 0.5d0*(rVel-lVel) &
                    *(erf((x_m-center_efunc_v)/smoother_efunc_v)+1.0d0)

                else if (flag_errorfunc_v == 2) then

                    eos_obj(i)%u = lVel + 0.5d0*(rVel-lVel) &
                    *(1.d0 - erf((x_m-center_efunc_v)/smoother_efunc_v))

                else if (flag_errorfunc_v == 3) then

                    eos_obj(i)%u = lVel + 0.5d0*(rVel-lVel) &
                    *(1.d0 - erf((dsqrt((x_m-center_efunc_v)**2)-radius_efunc_v) &
                                                              /smoother_efunc_v))
                end if

            else if ( flag_const_rhou == 1 ) then

                ! rho_u_FILL = lVel*density(lPres, lTemp, nspe, lXi)
                rho_u_FILL = rho_u_IN
                eos_obj(i)%u = rho_u_FILL/eos_obj(i)%rho

            end if

            eos_obj(i)%eInternal = specific_energy(eos_obj(i)%rho, eos_obj(i)%p, &
                                                   eos_obj(i)%T, nspe, eos_obj(i)%x_sp)


            eos_obj(i)%eKinetic = 0.5d0*eos_obj(i)%u**2
            eos_obj(i)%eTotal = eos_obj(i)%eInternal + eos_obj(i)%eKinetic

            eos_obj(i)%Qrad = 0.d0
            eos_obj(i)%qTw = 0.d0
            eos_obj(i)%HRR = 0.d0
            eos_obj(i)%HRR_dot = 0.d0

        end do

        if (tig_s < 1.0d-9) then
            flag_igsource = 0
        end if

    end subroutine initialize_variables

end module module_params


