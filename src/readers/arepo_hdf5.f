***********************************************************************
      subroutine read_arepo_hdf5_npart(iter, files_per_snap, 
     & flag_filter, flag_machfield)
***********************************************************************
*       reads the number of gas particles in the simulation
***********************************************************************
        use hdf5
        use particle_data  ! contains parti
        implicit none 

        integer :: iter, files_per_snap, flag_filter, flag_machfield
        character*3 iter_string
        character*1 ifile_string
        integer :: ifile
        character*200 fil1, fil2

        integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, 
     & file_space_id
        integer(hid_t) :: memtype_id
        integer :: status
        integer, dimension(6) :: numpart_thisfile
        integer(hsize_t), dimension(1) :: dims = 6

        ! default values for optional parameters
        integer(hid_t) :: xfer_prp
        !parameter (xfer_prp = h5p_default_f)  ! default transfer property list

        parti = 0
        do ifile = 0, files_per_snap - 1 
          write(iter_string, '(i3.3)') iter
          fil1 = './simulation/snapshot_' // iter_string
          if (files_per_snap .eq. 1) then
            fil2 = fil1
          else
            write(ifile_string, '(i1.1)') ifile
            fil2 = trim(adjustl(fil1)) // '.' // ifile_string
          end if
          fil2 = trim(adjustl(fil2)) // '.hdf5'

          ! open the hdf5 file in read-only mode
          call h5fopen_f(fil2, h5f_acc_rdonly_f, file_id, status)
          if (status /= 0) then
            print *, "error opening file: ", fil2
            stop
          end if

          ! open the header group
          call h5gopen_f(file_id, "/header", group_id, status)
          if (status /= 0) then
            print *, "error opening group: /header"
            call h5fclose_f(file_id, status)
            stop
          end if
      
          ! open the numpart_thisfile attribute
          call h5aopen_f(group_id, "numpart_thisfile",
     &                   attr_id, status)
          if (status /= 0) then
            print *, "error opening dataset: numpart_thisfile"
            call h5gclose_f(group_id, status)
            call h5fclose_f(file_id, status)
            stop
          end if

          ! get the datatype of the attribute to confirm it is h5t_native_integer
          call h5aget_type_f(attr_id, memtype_id, status)
          if (status /= 0) then
            print *, "error getting attribute type"
            call h5aclose_f(attr_id, status)
            call h5gclose_f(group_id, status)
            call h5fclose_f(file_id, status)
            stop
          end if
      
          ! read the attribute into the numpart_thisfile array
          call h5aread_f(attr_id, memtype_id, numpart_thisfile,
     &                   dims, status)

          if (status /= 0) then
            print *, "error reading attribute: numpart_thisfile"
            call h5aclose_f(attr_id, status)
            call h5gclose_f(group_id, status)
            call h5fclose_f(file_id, status)
            stop
          end if

          ! print the contents of numpart_thisfile
          ! print *, "numpart_thisfile: ", numpart_thisfile

          ! close the dataset, group, and file
          call h5aclose_f(attr_id, status)
          call h5gclose_f(group_id, status)
          call h5fclose_f(file_id, status)

          ! we just read parttype0
          parti = parti + numpart_thisfile(1)
        end do 

        write(*,*) 'total number of particles: ', parti
        return
      end
***********************************************************************

***********************************************************************
      subroutine read_arepo_hdf5(iter, files_per_snap,
     &             flag_filter, flag_machfield, low2)
***********************************************************************
*       reads the gas particle data of the simulation
***********************************************************************
      use hdf5
      use particle_data
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      character*3 iter_string
      character*1 ifile_string
      integer ifile
      character*200 fil1,fil2

      integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, 
     & file_space_id
      integer(hid_t) :: memtype_id
      integer :: status
      integer, dimension(6) :: numpart_thisfile
      integer(hsize_t), dimension(1) :: dims1d
      integer(hsize_t), dimension(2) :: dims2d
      
      integer basint1,basint2,basint3,basint4,i
      real*8 bas81,bas82,bas83,bas84,bas85,bas86
      integer low1,low2
      real*4,allocatable::scr4(:)
      real*4,allocatable::scr42(:,:)

      low2=0
      do ifile=0,files_per_snap-1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*     reading data
      write(iter_string, '(i3.3)') iter
      fil1 = './simulation/snapshot_' // iter_string
      if (files_per_snap .eq. 1) then
        fil2 = fil1
      else
        write(ifile_string, '(i1.1)') ifile
        fil2 = trim(adjustl(fil1)) // '.' // ifile_string
      end if
      fil2 = trim(adjustl(fil2)) // '.hdf5'

      ! open the hdf5 file in read-only mode
      write(*,*) 'reading iteration file: ',iter,' ',
     &              trim(adjustl(fil2))

      call h5fopen_f(fil2, h5f_acc_rdonly_f, file_id, status)
      if (status /= 0) then
        print *, "error opening file: ", fil2
        stop
      end if

      dims1d(1) = 6
      call h5gopen_f(file_id, "/header", group_id, status)
      call h5aopen_f(group_id, "numpart_thisfile",
     &                   attr_id, status)
      call h5aget_type_f(attr_id, memtype_id, status)
      call h5aread_f(attr_id, memtype_id, numpart_thisfile,
     &                   dims1d, status)

      write(*,*) numpart_thisfile(1), 'gas particles'
      low1=low2+1
      low2=low1+numpart_thisfile(1)-1

      call h5aclose_f(attr_id, status)
      call h5gclose_f(group_id, status)

      dims1d(1) = numpart_thisfile(1)
      dims2d(1) = numpart_thisfile(1)
      dims2d(2) = 3

      call h5gopen_f(file_id, "/parttype0",
     &                   group_id, status)
      if (status /= 0) then
        print *, "error opening group: /parttype0"
        call h5fclose_f(file_id, status)
        stop
      end if
      
      allocate(scr42(3,numpart_thisfile(1)))

      write(*,*) 'reading positions ...'
      call h5dopen_f(group_id, "coordinates", attr_id, status)
      call h5dget_type_f(attr_id, memtype_id, status)
      call h5dread_f(attr_id, memtype_id, scr42, dims2d, status)
      rxpa(low1:low2)=scr42(1,1:numpart_thisfile(1))
      rypa(low1:low2)=scr42(2,1:numpart_thisfile(1))
      rzpa(low1:low2)=scr42(3,1:numpart_thisfile(1))
      call h5dclose_f(attr_id, status)

      write(*,*) 'reading velocities ...'
      call h5dopen_f(group_id, "velocities", attr_id, status)
      call h5dget_type_f(attr_id, memtype_id, status)
      call h5dread_f(attr_id, memtype_id, scr42, dims2d, status)
      u2dm(low1:low2)=scr42(1,1:numpart_thisfile(1))
      u3dm(low1:low2)=scr42(2,1:numpart_thisfile(1))
      u4dm(low1:low2)=scr42(3,1:numpart_thisfile(1))
      call h5dclose_f(attr_id, status)

      deallocate(scr42)

      allocate(scr4(numpart_thisfile(1)))
      write(*,*) 'reading masses ...'
      call h5dopen_f(group_id, "masses", attr_id, status)
      call h5dget_type_f(attr_id, memtype_id, status)
      call h5dread_f(attr_id, memtype_id, scr4, dims1d, status)
      masap(low1:low2)=scr4(1:numpart_thisfile(1))
      call h5dclose_f(attr_id, status)

      write(*,*) 'reading effective kernel length ...' ! we read volume
      call h5dopen_f(group_id, "volume", attr_id, status)
      call h5dget_type_f(attr_id, memtype_id, status)
      call h5dread_f(attr_id, memtype_id, scr4, dims1d, status)
!$omp parallel do shared(numpart_thisfile, scr4, kernel, low1), 
!$omp+            private(i), default(none)
      do i=1,numpart_thisfile(1)
        kernel(low1+i-1) = (0.2387*scr4(i))**0.33333
      end do
      call h5dclose_f(attr_id, status)
      deallocate(scr4)

#if use_filter == 1
      if (flag_filter.eq.1) then
        allocate(scr4(numpart_thisfile(1)))
        if (flag_machfield.eq.0) then
          write(*,*) 'warning! in arepo, we always read mach!'
        end if 
        write(*,*) 'reading mach ...'
        call h5dopen_f(group_id, "machnumber", attr_id, status)
        call h5dget_type_f(attr_id, memtype_id, status)
        call h5dread_f(attr_id, memtype_id, scr4, dims1d, status)
        abvc(low1:low2)=scr4(1:numpart_thisfile(1))  
        call h5dclose_f(attr_id, status)         
        deallocate(scr4)      
      end if
#endif

#if weight_scheme == 2
      allocate(scr4(numpart_thisfile(1)))
      write(*,*) 'reading density ...'
      call h5dopen_f(group_id, "density", attr_id, status)
      call h5dget_type_f(attr_id, memtype_id, status)
      call h5dread_f(attr_id, memtype_id, scr4, dims1d, status)
      vol(low1:low2)=scr4(1:numpart_thisfile(1))
      call h5dclose_f(attr_id, status)
      deallocate(scr4)
#endif
        
      call h5gclose_f(group_id, status)
      call h5fclose_f(file_id, status)

      end do !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return 
      end 
***********************************************************************