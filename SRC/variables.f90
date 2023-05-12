module variables

    implicit none

    ! conserved quantities of fluid
    integer,parameter :: ncons = 3
    integer,parameter :: iudens = 1
    integer,parameter :: iumomx = 2
    integer,parameter :: iuener = 3

    ! primitive quantities of fluid
    integer,parameter :: nprim = 3
    integer,parameter :: iqdens = 1
    integer,parameter :: iqxvel = 2
    integer,parameter :: iqtemp = 3
    integer,parameter :: iqpres = 3

    ! characteristic waves
    ! integer, parameter :: nwaves = 3

end module variables

