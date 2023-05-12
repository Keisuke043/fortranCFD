
module parallel

        implicit none

        include 'mpif.h'

        ! MPI variable
        integer :: ISTATUS(MPI_STATUS_SIZE)

        integer :: nprocs, myrank
        integer :: IUP, IDOWN
        integer :: nx_total  ! number of total cell

        integer :: ierr

contains


    subroutine init_parallel

        implicit none


        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    end subroutine init_parallel


    subroutine fin_parallel

        implicit none

        
        call MPI_FINALIZE(ierr)
    
    end subroutine fin_parallel


    subroutine set_parallel_block(bsta, bend, nblock_0, nblock_n)

        implicit none

        integer,intent(in) :: nblock_0, nblock_n
        integer,intent(inout) :: bsta, bend


        call block_range(nblock_0, nblock_n, myrank, nprocs, bsta, bend)

        if ( myrank == 0        ) bsta = nblock_0
        if ( myrank == nprocs-1 ) bend = nblock_n

        IUP   = myrank + 1
        IDOWN = myrank - 1

        if ( myrank == nprocs-1 ) IUP   = MPI_PROC_NULL
        if ( myrank == 0        ) IDOWN = MPI_PROC_NULL

        ! print *, 'myrank = ', myrank, 'bsta = ', bsta, 'bend = ', bend

    end subroutine set_parallel_block


    subroutine block_range(n1, n2, myrank, nprocs, bsta, bend)

        implicit none

        integer :: n1, n2
        integer :: bsta, bend
        integer :: iwork1, iwork2
        integer :: myrank, nprocs


        iwork1 = (N2-N1+1) / nprocs
        iwork2 = mod(  N2-N1+1,   nprocs )
        bsta = myrank*iwork1 + N1 + min(myrank, iwork2)
        bend = bsta + iwork1 - 1

        if (iwork2 > myrank) bend = bend + 1

    end subroutine block_range



    subroutine set_parallel(ista_b, iend_b, nx_all_0, nx_all_n)

        implicit none

        integer :: ista_b, iend_b
        integer :: nx_all_0, nx_all_n

        
        call para_range(nx_all_0, nx_all_n, myrank, nprocs, ista_b, iend_b)

        if ( myrank == 0        ) ista_b = nx_all_0
        if ( myrank == nprocs-1 ) iend_b = nx_all_n

        IUP   = myrank + 1
        IDOWN = myrank - 1

        if ( myrank == nprocs-1 ) IUP   = MPI_PROC_NULL
        if ( myrank == 0        ) IDOWN = MPI_PROC_NULL

        ! print *, 'myrank = ', myrank, 'nprocs = ', nprocs
        ! print *, 'myrank = ', myrank, 'ista_b = ', ista_b, 'iend_b = ', iend_b

    end subroutine set_parallel


    subroutine para_range(n1, n2, myrank, nprocs, ista_b, iend_b)

        implicit none


        integer :: n1, n2
        integer :: ista_b, iend_b
        integer :: iwork1, iwork2
        integer :: myrank, nprocs


        iwork1 = int( (N2-N1+1) / nprocs )
        iwork2 = mod(  N2-N1+1,   nprocs )
        ista_b = myrank*iwork1 + N1 + min(myrank, iwork2)
        iend_b = ista_b + iwork1 - 1
        
        if (iwork2>myrank) iend_b = iend_b + 1


    end subroutine para_range


end module parallel


