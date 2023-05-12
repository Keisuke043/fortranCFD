
module read_hdf5

    USE HDF5 ! This module contains all necessary modules
    USE parallel
    use thermo

    implicit none

contains


    subroutine read_hdf5_restart(ista_rank, iend_rank, nspe, eos_obj, &
                                 nx_all_0, nx_all_n, readpath, &
                                 nout_rstr, nt_rstr, time_rstr, cputime_rstr)

        use module_params, only:species_name, eos_st

        implicit none
        

        integer :: access

        integer(HID_T) :: file_id       ! File identifier
        integer(HID_T) :: plist_id      ! Property list identifier
        integer(HID_T) :: step_group_id ! Group identifier
        integer(HID_T) :: dset_id       ! Dataset identifier
        integer(HID_T) :: dspace_id     ! Dataspace identifier
        integer        :: hdferror      ! Error flags
        character(len=128) :: groupname ! Group name
        character(len=128) :: dsetname  ! Dataset name
        character(len=128) :: filepath  ! Restart file path
        integer(SIZE_T)    :: filepath_size
        integer(HSIZE_T), dimension(1) :: dims, maxdims ! Dimensions of dataspace


        integer :: i, k_sp
        integer,intent(in) :: ista_rank, iend_rank
        integer,intent(in) :: nspe
        integer,intent(in) :: nx_all_0, nx_all_n
        character(len=128),intent(in) :: readpath
        type(eos_st),intent(inout) :: eos_obj(ista_rank:iend_rank)
        integer,intent(inout) :: nout_rstr
        integer,intent(inout) :: nt_rstr
        real(8),intent(inout) :: time_rstr
        real(8),intent(inout) :: cputime_rstr


        integer :: storage_type1      ! Type of storage for links in group
        integer :: nlinks_group_num   ! Number of links in group
        integer :: max_corder1        ! Current maximum creation order value for group
        logical :: mounted1           ! Whether group has a file mounted on it
        integer :: nmembers1

        integer :: storage_type2      ! Type of storage for links in group
        integer :: nlinks_dataset_num ! Number of links in group
        integer :: max_corder2        ! Current maximum creation order value for group
        logical :: mounted2           ! Whether group has a file mounted on it
        integer :: nmembers2

        integer :: idx                ! Index for loop
        integer :: idx_stepgroup      ! Index of stepgroup member object
        integer :: obj_type_group     ! Group object type

        integer :: i_var, loop_var
        integer :: var_idx_0
        character(len=128),allocatable :: varname(:)
        real(8),allocatable            :: var(:,:)


        real(8),allocatable :: attr_var_time(:)
        real(8),allocatable :: attr_var_cputime(:)
        integer,allocatable :: attr_var_nout(:)
        integer,allocatable :: attr_var_nt(:)
        character(len=128)  :: aname     ! Attribute name
        integer(HID_T)      :: attr_id   ! Attribute identifier
        integer(HID_T)      :: atype_id  ! Attribute Datatype identifier
        integer,parameter   :: nadim = 1 ! Attribure rank
        integer(HSIZE_T), dimension(nadim) :: adims ! Attribute dimension
        integer(HSIZE_T), dimension(nadim) :: amaxdims  

        if (access(readpath, " ") /= 0) then
            print *, 'There is no file in ', trim(readpath)
            call fin_parallel
            stop
        end if

        ! Initialize FORTRAN interface.

        call h5open_f(hdferror)


        ! Setup file access property list with parallel I/O access.

        call h5pcreate_f (H5P_FILE_ACCESS_F, plist_id, hdferror)
        call h5pset_fapl_mpio_f (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferror)


        ! Open an existing file collectively.

        call h5fopen_f (readpath, H5F_ACC_RDWR_F, file_id, hdferror, &
                        access_prp = plist_id)
        call h5pclose_f (plist_id, hdferror)


        ! Returns the length of the filename

        call h5fget_name_f(file_id, filepath, filepath_size, hdferror)
        if (myrank == 0) print*, 'filepath:  ', trim(filepath)


        ! Retrieves information about the group specified by group_id

        call h5gget_info_f(file_id, storage_type1, nlinks_group_num, &
                           max_corder1, hdferror, mounted1)


        ! Retrieves the number of group members

        call h5gn_members_f(file_id, './', nmembers1, hdferror)


        ! Retrieve name and type of the group member identified by its index

        do idx = 0,  nlinks_group_num-1
            call h5gget_obj_info_idx_f(file_id, '/', idx, groupname, &
                                       obj_type_group, hdferror)
            if (myrank == 0) print *, 'idx:', idx, 'obj_name:', trim(groupname), &
                                      'obj_type:', obj_type_group
        end do
        idx_stepgroup = nlinks_group_num-1
        call h5gget_obj_info_idx_f(file_id, '/', idx_stepgroup, &
                                   groupname, obj_type_group, hdferror)
        if (myrank == 0) print *, 'idx:', idx_stepgroup, trim(groupname)


        ! Open an existing group in the specified file.

        call h5gopen_f(file_id, groupname, step_group_id, hdferror)


        ! Retrieves information about the group specified by group_id

        call h5gget_info_f(step_group_id, storage_type2, nlinks_dataset_num, &
                           max_corder2, hdferror, mounted2)


        ! Retrieves the number of group members

        call h5gn_members_f(file_id, groupname, nmembers2, hdferror)


        ! Retrieve name and type of the group member identified by its index

        loop_var = nlinks_dataset_num
        allocate( varname(loop_var) )
        do idx = 0,  nlinks_dataset_num-1
            call h5gget_obj_info_idx_f(file_id, groupname, idx, &
                                       varname(idx+1), obj_type_group, hdferror)
            if (myrank == 0) print *, 'idx:', idx, 'obj_name:', &
                                       trim(varname(idx+1)), &
                                      'obj_type:', obj_type_group
        end do


        ! Retrieves dataspace dimension size and maximum size.

        dsetname = trim(varname(1))
        call h5dopen_f(step_group_id, dsetname, dset_id, hdferror)
        call h5dget_space_f(dset_id, dspace_id, hdferror)
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferror)
        call h5dclose_f(dset_id, hdferror)
        allocate( var(loop_var, dims(1)) )


        do i_var = 1, loop_var

            ! Open an existing dataset.

            dsetname = trim(varname(i_var))
            call h5dopen_f(step_group_id, dsetname, dset_id, hdferror)


            ! Retrieves dataspace dimension size and maximum size.

            call h5dget_space_f(dset_id, dspace_id, hdferror)
            call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferror)


            ! Create property list

            call h5pcreate_f (H5P_DATASET_CREATE_F, plist_id, hdferror)


            ! Create property list for collective dataset write.

            call h5pcreate_f (H5P_DATASET_XFER_F, plist_id, hdferror)
            call h5pset_dxpl_mpio_f (plist_id, H5FD_MPIO_COLLECTIVE_F, hdferror)


            ! Read the dataset

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, var(i_var,:), dims, hdferror)


            ! Close the dataset.

            call h5dclose_f(dset_id, hdferror)


            ! Close the property list.

            call h5pclose_f (plist_id, hdferror)

        end do


        ! read time in attribute 
        aname = 'time from 0s'

        ! Open attribute.
        CALL h5aopen_f(step_group_id, aname, attr_id, hdferror)

        ! Get the datatype.
        call h5aget_type_f(attr_id, atype_id, hdferror)

        ! Get dataspace and allocate memory for read buffer.
        call h5aget_space_f(attr_id, dspace_id, hdferror)
        call h5sget_simple_extent_dims_f(dspace_id, adims, amaxdims, hdferror)

        allocate( attr_var_time(1:adims(1)) )

        ! Read the attribute data.
        call h5aread_f(attr_id, atype_id, attr_var_time, adims, hdferror)

        ! Close attribute.
        call h5aclose_f(attr_id, hdferror)


        ! read cpu time in attribute 
        aname = 'cpu time total'

        ! Open attribute.
        CALL h5aopen_f(step_group_id, aname, attr_id, hdferror)

        ! Get the datatype.
        call h5aget_type_f(attr_id, atype_id, hdferror)

        ! Get dataspace and allocate memory for read buffer.
        call h5aget_space_f(attr_id, dspace_id, hdferror)
        call h5sget_simple_extent_dims_f(dspace_id, adims, amaxdims, hdferror)

        allocate( attr_var_cputime(1:adims(1)) )

        ! Read the attribute data.
        call h5aread_f(attr_id, atype_id, attr_var_cputime, adims, hdferror)

        ! Close attribute.
        call h5aclose_f(attr_id, hdferror)


        ! read number of output in attribute 
        aname = 'number of output'

        ! Open attribute.
        CALL h5aopen_f(step_group_id, aname, attr_id, hdferror)

        ! Get the datatype.
        call h5aget_type_f(attr_id, atype_id, hdferror)

        ! Get dataspace and allocate memory for read buffer.
        call h5aget_space_f(attr_id, dspace_id, hdferror)
        call h5sget_simple_extent_dims_f(dspace_id, adims, amaxdims, hdferror)

        allocate( attr_var_nout(1:adims(1)) )


        ! Read the attribute data.
        call h5aread_f(attr_id, atype_id, attr_var_nout, adims, hdferror)

        ! Close attribute.
        call h5aclose_f(attr_id, hdferror)

        ! read number of output in attribute 
        aname = 'number of time step'

        ! Open attribute.
        CALL h5aopen_f(step_group_id, aname, attr_id, hdferror)

        ! Get the datatype.
        call h5aget_type_f(attr_id, atype_id, hdferror)

        ! Get dataspace and allocate memory for read buffer.
        call h5aget_space_f(attr_id, dspace_id, hdferror)
        call h5sget_simple_extent_dims_f(dspace_id, adims, amaxdims, hdferror)

        allocate( attr_var_nt(1:adims(1)) )

        ! Read the attribute data.
        call h5aread_f(attr_id, atype_id, attr_var_nt, adims, hdferror)

        ! Close attribute.
        call h5aclose_f(attr_id, hdferror)


        time_rstr = attr_var_time(1)
        cputime_rstr = attr_var_cputime(1)
        nout_rstr = attr_var_nout(1)+1
        nt_rstr = attr_var_nt(1)+1


        ! Close the group.
        call h5gclose_f(step_group_id, hdferror)

        ! Close the file.
        call h5fclose_f(file_id, hdferror)

        ! Close FORTRAN interface.
        call h5close_f(hdferror)


        do i_var = 1, loop_var
            if ('idx' == trim(varname(i_var))) then
                var_idx_0 = int(var(i_var, 1)-1)
                if (myrank == 0) print *, 'var_idx_0', var_idx_0
            end if
        end do

        do i_var = 1, loop_var
            if ('x_m' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    if (eos_obj(i)%x_m < var(i_var, i-var_idx_0)-1d-7 .or. &
                        var(i_var, i-var_idx_0)+1d-7 < eos_obj(i)%x_m) then
                        print *, 'x_m in the input file does not match &
                                &that in the hdf5 read file.'
                        print *, i, var(i_var, i-var_idx_0), eos_obj(i)%x_m
                        ! call fin_output
                        ! call fin_parallel
                        stop
                    end if
                end do
            else if ('idx' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    if ( i /= int(var(i_var, i-var_idx_0)) ) then
                        print *, 'idx in the input file does not match &
                                 &that in the hdf5 read file.'
                        print *, i, var(i_var, i-var_idx_0)
                        ! call fin_output
                        ! call fin_parallel
                        stop
                    end if
                end do
            else if (' dx' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    if (eos_obj(i)%dx < var(i_var, i-var_idx_0)-1d-7 .or. &
                        var(i_var, i-var_idx_0)+1d-7 < eos_obj(i)%dx) then
                        print *, 'dx in the input file does not match &
                                 &that in the hdf5 read file.'
                        print *, i, var(i_var, i-var_idx_0), eos_obj(i)%dx
                        ! call fin_output
                        ! call fin_parallel
                        stop
                    end if
                end do
            end if
        end do
            
        do i_var = 1, loop_var

            if (' dx' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%dx=var(i_var, i-var_idx_0)
                end do
            end if

            if ('rho' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%rho=var(i_var, i-var_idx_0)
                end do
            else if ('u_x' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%u=var(i_var, i-var_idx_0)
                end do
            else if ('  T' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%T=var(i_var, i-var_idx_0)
                end do
            else if ('  p' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%p=var(i_var, i-var_idx_0)
                end do
            else if (' et' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%eTotal=var(i_var, i-var_idx_0)
                end do
            else if (' ek' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%eKinetic=var(i_var, i-var_idx_0)
                end do
            else if (' ei' == trim(varname(i_var))) then
                do i = ista_rank, iend_rank
                    eos_obj(i)%eInternal=var(i_var, i-var_idx_0)
                end do

            end if
            do k_sp = 1, nspe
                if (species_name(k_sp) == trim(varname(i_var))) then
                    do i = ista_rank, iend_rank
                        eos_obj(i)%x_sp(k_sp)=var(i_var, i-var_idx_0)
                    end do
                end if
            end do

        end do
        do i = ista_rank, iend_rank
            call Xi2Yi(nspe, eos_obj(i)%x_sp(:), eos_obj(i)%y_sp(:))
        end do


        deallocate(varname, var)


    end subroutine read_hdf5_restart


end module read_hdf5

