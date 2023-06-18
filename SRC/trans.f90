module trans
  use iso_c_binding
  use thermo

  implicit none

  ! +++++ TRANLIB ì‹Æ—p”z—ñ +++++
  integer     ,allocatable :: IMCWRK(:)
  real(kind=8),allocatable :: RMCWRK(:)
  ! species bundle
  integer, save :: nmax_sbgrp
  integer, save, allocatable:: numsgrp(:)
  integer, save, allocatable:: indxrps(:)
  integer, save, allocatable:: indxgrp(:)
  integer, save, allocatable:: idspgrp(:,:)

  contains

    ! subroutine set_trans(chem_dir) bind(c)
    subroutine set_trans(chem_dir)


      ! +++++ Local parameters +++++
      character(len=128),intent(in) :: chem_dir
      character(len=128) :: filename_tran
      integer           ,parameter :: LINKMC        = 35
      integer           ,parameter :: LOUT          = 6
      ! +++++ Local variables +++++
      integer                      :: st_open
      integer                      :: LENIMC
      integer                      :: LENRMC
      integer                      :: iflag
      ! ------------------------------------------------------------------------

      write (filename_tran, '(a,"/",a)') trim(chem_dir), 'tran.bin'

      ! +++ open the TRANLIB link file +++
      open(LINKMC, file=trim(filename_tran), form='UNFORMATTED', status='UNKNOWN')
      call MCLEN( LINKMC, LOUT, LENIMC, LENRMC, iflag )
      allocate( IMCWRK(LENIMC) )
      allocate( RMCWRK(LENRMC) )
      call MCINIT( LINKMC, LOUT, LENIMC, LENRMC, IMCWRK, RMCWRK, iflag )
    end subroutine set_trans


    function viscosity( rho, p, T, n, Xi )
      real(8),intent(in) :: rho
      real(8),intent(in) :: p
      real(8),intent(in) :: T
      integer     ,intent(in) :: n
      real(8),intent(in) :: Xi(n)
      real(8)            :: viscosity

      ! +++++ Local variables +++++
      real(kind=8)            :: vismix
      ! ------------------------------------------------------------------------
      call MCAVIS( T, Xi, RMCWRK, vismix )

      viscosity = vismix * 0.1d0 ! cgs -> SI ’PˆÊŠ·ŽZ
    end function viscosity


    function thermal_conductivity( rho, p, T, n, Xi )
      real(8),intent(in) :: rho
      real(8),intent(in) :: p
      real(8),intent(in) :: T
      integer,intent(in) :: n
      real(8),intent(in) :: Xi(n)
      real(8)            :: thermal_conductivity

      ! +++++ Local variables +++++
      real(kind=8)            :: conmix
      ! ------------------------------------------------------------------------
      call MCACON( T, Xi, RMCWRK, conmix )

      thermal_conductivity = conmix * 1.d-5 ! cgs -> SI ’PˆÊŠ·ŽZ
    end function thermal_conductivity


    subroutine set_diffusivity( rho, p, T, n, Xi, Di )
      real(8),intent(in)  :: rho
      real(8),intent(in)  :: p
      real(8),intent(in)  :: T
      integer     ,intent(in)  :: n
      real(8),intent(in)  :: Xi(n)
      real(8),intent(out) :: Di(n)

      ! +++++ Local variables +++++
      real(kind=8)             :: p_cgs
      real(kind=8)             :: d_cgs(n_species)
      real(kind=8)             :: dt_cgs(n_species)
      real(kind=8)             :: cond
      ! ------------------------------------------------------------------------
      ! print *, n_species
      p_cgs = p * 10.d0 ! SI -> cgs ’PˆÊŠ·ŽZ

      call MCADIF( p_cgs, T, Xi, RMCWRK, d_cgs )

      Di(1:n) = d_cgs(1:n)  * 1.d-4 ! cgs -> SI ’PˆÊŠ·ŽZ
    end subroutine set_diffusivity

    ! subroutine sbinput() bind(c)
    !   implicit none
    !   integer i,j,idum,jdum,kkdum

    !   open(10,file='spbdl.list',form='formatted')

    !   read(10,*) nmax_sbgrp
    !   allocate(numsgrp(nmax_sbgrp))
    !   do i = 1,nmax_sbgrp
    !     read(10,*) idum,numsgrp(i)
    !   enddo 
    !   allocate(indxrps(nmax_sbgrp))
    !   do i = 1,nmax_sbgrp
    !     read(10,*) idum,indxrps(i)
    !   enddo 
    !   read(10,*) kkdum
    !   if (kkdum/=nspecies) then
    !     write(*,*) 'wrong species# in sbinit',kkdum
    !     stop
    !   endif
    !   allocate(indxgrp(nspecies))
    !   do i = 1,nspecies
    !     read(10,*) idum,indxgrp(i)
    !   enddo
    !   allocate(idspgrp(nspecies,nmax_sbgrp))
    !   do i = 1,nmax_sbgrp
    !     do j = 1,numsgrp(i)
    !       read(10,*) jdum,idum,idspgrp(j,i)
    !     enddo
    !   enddo

    !   close(10)
    ! endsubroutine sbinput


    function viscosity_ea( rho, p, T, n, Xi )
      real(8),intent(in) :: rho
      real(8),intent(in) :: p
      real(8),intent(in) :: T
      integer     ,intent(in) :: n
      real(8),intent(in) :: Xi(n)
      real(8)            :: viscosity_ea

      ! +++++ Local variables +++++
      real(kind=8)            :: vismix
      ! ------------------------------------------------------------------------
      call MCAVIS_EA( T, Xi, RMCWRK, vismix )

      viscosity_ea = vismix * 0.1d0 ! cgs -> SI ’PˆÊŠ·ŽZ
    end function viscosity_ea

    ! subroutine set_diffusivity_sb( rho, p, T, n, Xi, Yi, Di )
    !   real(8),intent(in)  :: rho
    !   real(8),intent(in)  :: p
    !   real(8),intent(in)  :: T
    !   integer     ,intent(in)  :: n
    !   real(8),intent(in)  :: Xi(n), Yi(n)
    !   real(8),intent(out) :: Di(n)

    !   ! +++++ Local variables +++++
    !   real(kind=8)             :: p_cgs
    !   real(kind=8)             :: d_cgs(NSPECIES)
    !   real(kind=8)             :: dt_cgs(NSPECIES)
    !   real(kind=8)             :: cond
    !   real(kind=8),allocatable :: dmij(:,:)
    !   allocate(dmij(nmax_sbgrp,nmax_sbgrp))
    !   ! ------------------------------------------------------------------------
    !   p_cgs = p * 10.d0 ! SI -> cgs ’PˆÊŠ·ŽZ

    !   call sb_mcsdif( p_cgs, T, Xi, Yi, nspecies, nmax_sbgrp, RMCWRK, &
    !                 & indxrps,indxgrp,numsgrp,idspgrp,dmij,d_cgs )

    !   Di(1:n) = d_cgs(1:n)  * 1.d-4 ! cgs -> SI ’PˆÊŠ·ŽZ
    !   deallocate(dmij)
    ! end subroutine set_diffusivity_sb

end module trans
