
module hdf5_mpi


    USE HDF5 ! This module contains all necessary modules
    USE parallel
    USE module_params, only: species_name


    implicit none

    character(len=128) :: filename_h5 = 'results'
    character(len=128) :: savepath

    integer(HID_T) :: file_id  ! File identifier
    integer(HID_T) :: plist_id ! Property list identifier

    integer        :: hdferror    ! Error flags

contains


    subroutine init_output(savedir, nout)

        implicit none

        integer :: nout
        character :: savedir*128    ! directory name


        write (savepath,'("./",a,"/",a,"_nout",i6.6,".h5")') &
                        trim(savedir), trim(filename_h5), nout


        ! Initialize HDF5 library and Fortran interfaces.

        call h5open_f (hdferror)


        ! Setup file access property list with parallel I/O access.

        call h5pcreate_f (H5P_FILE_ACCESS_F, plist_id, hdferror)
        call h5pset_fapl_mpio_f (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferror)


        ! Create the file collectively.

        call h5fcreate_f (savepath, H5F_ACC_TRUNC_F, file_id, hdferror, access_prp = plist_id)
        call h5pclose_f (plist_id, hdferror)

        return

    end subroutine init_output


    subroutine fin_output

        implicit none

        
        ! Close the file.
        
        call h5fclose_f (file_id, hdferror)


        ! Close HDF5 library and Fortran interfaces.

        call h5close_f (hdferror)

        return

    end subroutine fin_output


    subroutine write_hdf5_b(nout, nt, comp_time, comp_time_total, dt, &
                            time, time_total, ista_rank, iend_rank, lvir, & 
                            nspe, eos_obj, nx_all_0, nx_all_n, chem_dir)

        use module_params, only: eos_st

        implicit none


        integer,intent(inout) :: nout
        integer,intent(in) :: nt
        real(8),intent(in) :: comp_time, comp_time_total
        real(8),intent(in) :: dt, time, time_total
        integer,intent(in) :: ista_rank, iend_rank, lvir, nspe
        integer,intent(in) :: nx_all_0, nx_all_n
        type(eos_st),intent(in) :: eos_obj(ista_rank:iend_rank)
        character(len=128),intent(in) :: chem_dir


        character(len=128) :: stepname      ! Group name
        character(len=128) :: dsetname      ! Dataset name
        integer(HID_T)     :: step_group_id ! Group identifier
        integer(HID_T)     :: dset_id       ! Dataset identifier
        integer(HID_T)     :: dspace_id     ! Dataspace identifier in file
        integer(HID_T)     :: memdspace_id  ! Memory dataspace identifier

        integer,parameter :: ndim = 1
        integer(HSIZE_T), dimension(ndim) :: subset_dims ! Dataset dimensions of subset
        integer(HSIZE_T), dimension(ndim) :: dims        ! Dataset dimensions 

        integer(HSIZE_T), dimension(ndim) :: offset ! Hyperslab offset
        integer(HSIZE_T), dimension(ndim) :: counts ! Hyperslab size
        integer(HSIZE_T), dimension(ndim) :: stride ! Hyperslab stride
        integer(HSIZE_T), dimension(ndim) :: blocks ! Hyperslab block size


        integer :: i_attr, loop_attr(1), loop_real(1), loop_int(1)
        integer :: n_attr = 8
        integer :: n_attr_real = 5
        character(len=32),allocatable :: varname_a(:)
        real(8),allocatable           :: var_real_a(:)
        integer,allocatable           :: var_int_a(:)
        character(len=128) :: aname     ! Attribute name
        integer(HID_T)     :: attr_id   ! Attribute identifier
        integer(HID_T)     :: atype_id  ! Attribute Datatype identifier
        integer(HID_T)     :: aspace_id ! Attribute Dataspace identifier

        integer,parameter :: nadim = 1                    ! Attribure rank
        integer(HSIZE_T), dimension(nadim)   :: adims     ! Attribute dimension
        integer(HSIZE_T), dimension(nadim)   :: data_dims ! Attribute dimension
        integer, dimension(nadim)            :: attr_var_int   ! Attribute variable
        real(8), dimension(nadim)            :: attr_var_real  ! Attribute variable
        integer, dimension(nadim)            :: attr_int  ! Attribute integer variable
        character(len=128), dimension(nadim) :: attr_char ! Attribute character
        integer(SIZE_T)                      :: attrlen   ! Length of the attribute string


        integer :: i, k_sp
        integer :: i_var, loop_var(2)
        integer :: n_wospe = 15
        character(len=32),allocatable :: varname(:)
        real(8),allocatable           :: var(:,:)


        subset_dims = (/iend_rank-ista_rank+1/) ! Dataset dimensions of subset
        dims = (/nx_all_n-nx_all_0+1/)          ! Dataset dimensions

        allocate(var(n_wospe+nspe,ista_rank:iend_rank))
        var(1,:) = eos_obj(:)%rho
        var(2,:) = eos_obj(:)%u
        var(3,:) = eos_obj(:)%T
        var(4,:) = eos_obj(:)%p
        var(5,:) = eos_obj(:)%eTotal
        var(6,:) = eos_obj(:)%eKinetic
        var(7,:) = eos_obj(:)%eInternal
        var(8,:) = eos_obj(:)%HRR
        var(9,:) = eos_obj(:)%HRR_dot
        var(10,:) = eos_obj(:)%Qrad
        var(11,:) = eos_obj(:)%Tw
        var(12,:) = eos_obj(:)%qTw
        var(13,:) = eos_obj(:)%x_m
        var(14,:) = eos_obj(:)%dx
        do i = ista_rank, iend_rank
            var(n_wospe,i) = i
            do k_sp = 1, nspe
                var(n_wospe+k_sp,i) = eos_obj(i)%x_sp(k_sp)
            end do
        end do


        allocate(varname(1:n_wospe+nspe))
        varname(1:n_wospe) = (/'rho','u_x','  T','  p',' et',' ek',' ei','HRR','HRD','Qrd', ' Tw','QTw',&
                               'x_m',' dx','idx'/)
        do k_sp = 1, nspe
            varname(k_sp+n_wospe) = species_name(k_sp)
        end do

        offset(:) = (/ista_rank - nx_all_0/)
        counts(:) = (/1/)
        stride(:) = (/1/)
        blocks(:) = subset_dims


        ! Create the timestep group for the specific I/O step.

        write (stepname, '("nout",i6.6,"_time",e11.5)') nout, time_total
        call h5gcreate_f (file_id, stepname, step_group_id, hdferror)


        loop_var=shape(var)

        do i_var = 1, loop_var(1)

            ! Create dataspace for the dataset. Setting maximum size to be the current size.

            call h5screate_simple_f (ndim, dims, dspace_id, hdferror)


            ! Create property list

            call h5pcreate_f (H5P_DATASET_CREATE_F, plist_id, hdferror)


            ! Create the dataset.

            dsetname = trim(varname(i_var))
            call h5dcreate_f (step_group_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                              dset_id, hdferror, plist_id)


            ! Close dataspaces.

            call h5sclose_f (dspace_id, hdferror)


            ! Select hyperslab in the file.

            call h5dget_space_f (dset_id, dspace_id, hdferror)
            call h5sselect_hyperslab_f (dspace_id, H5S_SELECT_SET_F, offset, counts, &
                                        hdferror, stride, blocks)

            ! Specify memory dataspace.

            call h5screate_simple_f (ndim, subset_dims, memdspace_id, hdferror)


            ! Create property list for collective dataset write.

            call h5pcreate_f (H5P_DATASET_XFER_F, plist_id, hdferror)
            call h5pset_dxpl_mpio_f (plist_id, H5FD_MPIO_COLLECTIVE_F, hdferror)


            ! Write the dataset collectively.

            call h5dwrite_f (dset_id, H5T_NATIVE_DOUBLE, var(i_var,ista_rank:iend_rank),  &
                             subset_dims, ierr, file_space_id = dspace_id,  &
                             mem_space_id = memdspace_id, xfer_prp = plist_id )


            ! Close the dataset.

            call h5dclose_f (dset_id, hdferror)


            ! Close the property list.

            call h5pclose_f (plist_id, hdferror)


            ! Close dataspaces.

            call h5sclose_f (memdspace_id, hdferror)
            call h5sclose_f (dspace_id, hdferror)

        end do

        adims = (/1/)

        allocate(varname_a(1:n_attr))
        varname_a(1) = 'dt'
        varname_a(2) = 'time from 0s'
        varname_a(3) = 'time from restart'
        varname_a(4) = 'cpu time total'
        varname_a(5) = 'cpu time restart'
        varname_a(6) = 'number of output'
        varname_a(7) = 'number of time step'
        varname_a(8) = 'chem directory name'

        allocate(var_real_a(1:n_attr_real))
        var_real_a(1) = dt
        var_real_a(2) = time_total
        var_real_a(3) = time
        var_real_a(4) = comp_time_total
        var_real_a(5) = comp_time

        allocate(var_int_a(1:2))
        var_int_a(1) = nout
        var_int_a(2) = nt

        loop_attr=shape(varname_a)
        loop_real=shape(var_real_a)
        loop_int =shape(var_int_a)

        do i_attr = 1, loop_attr(1)-1

            aname = varname_a(i_attr)
            
      
            ! Create scalar data space for the attribute.
            call h5screate_simple_f(nadim, adims, aspace_id, hdferror)

            if (i_attr <= n_attr) then

                if (i_attr <= loop_real(1)) then

                    ! Initialize attribute's data
                    attr_var_real(1) = var_real_a(i_attr)

                    ! Create dataset attribute.
                    call h5acreate_f(step_group_id, aname, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferror)

                    ! Write the attribute data.
                    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_var_real, adims, hdferror)

                else if (loop_real(1) < i_attr .and. i_attr < loop_attr(1)) then

                    attr_var_int(1) = var_int_a(i_attr-loop_real(1))

                    ! Create dataset attribute.
                    call h5acreate_f(step_group_id, aname, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferror)

                    ! Write the attribute data.
                    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_var_int, adims, hdferror)

                end if

            else if (i_attr == loop_attr(1)) then

                ! Initialize attribute's data
                attrlen = 32
                attr_char(1) = trim(chem_dir)

                ! Create datatype for the attribute.
                CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferror)
                CALL h5tset_size_f(atype_id, attrlen, hdferror)

                ! Create dataset attribute.
                call h5acreate_f(step_group_id, aname, atype_id, aspace_id, attr_id, hdferror)
            
                ! Write the attribute data.
                call h5awrite_f(attr_id, atype_id, attr_char(1), adims, hdferror)

            end if
             
            ! Close the attribute.
            call h5aclose_f(attr_id, hdferror)

            ! Close dataspaces.
            call h5sclose_f (aspace_id, hdferror)

        end do

        ! Close the timestep group.

        call h5gclose_f (step_group_id, hdferror)

        deallocate(varname, var)

        nout = nout + 1
       
        return

    end subroutine write_hdf5_b


end module hdf5_mpi

