module mymkdir

    implicit none


contains

    subroutine makedirs(outdir)
    
        implicit none


        character(len=*), intent(in) :: outdir
        character(len=256) command

    
        write(command,*) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', &
                          trim(outdir), '; fi'

        call system(command)
    
    end subroutine makedirs

end module mymkdir

