***********************************************************************
      subroutine read_gadget_unformatted_npart(iter, files_per_snap,
     &             flag_filter, flag_machfield)
***********************************************************************
*       reads the number of gas particles in the simulation
***********************************************************************
      use gadget_read
      use particle_data  ! contains parti
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      character*3 iter_string
      character*7 ifile_string
      integer ifile
      character*200 fil1,fil2

      integer npart_gadget(6), nall(6), blocksize
      real*8 massarr(6)
      integer basint1,basint2,basint3,basint4
      real*8 bas81,bas82,bas83,bas84,bas85,bas86

      parti=0
      do ifile=0,files_per_snap-1
        write(iter_string,'(i3.3)') iter
        fil1='./simulation/'
        if (files_per_snap.eq.1) then
          fil2=trim(adjustl(fil1))//'snap_'//iter_string
        else
          fil1=trim(adjustl(fil1))//'snapdir_'//iter_string//
     &                            '/snap_'//iter_string
          write(ifile_string,'(i7.1)') ifile
          fil2=trim(adjustl(fil1))//'.'//trim(adjustl(ifile_string))
        end if

        call read_head(fil2,npart_gadget,massarr,bas81,bas82,
     &                  basint1,basint2,nall,basint3,basint4,bas83,
     &                  bas84,bas85,bas86,blocksize)

        parti=parti+npart_gadget(1)
      end do 

      return 
      end 
***********************************************************************

***********************************************************************
      subroutine read_gadget_unformatted(iter, files_per_snap,
     &             flag_filter, flag_machfield, low2)
***********************************************************************
*       reads the gas particle data of the simulation
***********************************************************************
      use gadget_read 
      use particle_data
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      character*3 iter_string
      character*7 ifile_string
      integer ifile
      character*200 fil1,fil2

      integer npart_gadget(6), nall(6), blocksize
      real*8 massarr(6)
      integer basint1,basint2,basint3,basint4,i
      real*8 bas81,bas82,bas83,bas84,bas85,bas86
      integer low1,low2
      real*4,allocatable::scr4(:)
      real*4,allocatable::scr42(:,:)

      low2=0
      do ifile=0,files_per_snap-1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*     reading data
      write(iter_string,'(i3.3)') iter
      fil1='./simulation/'
      if (files_per_snap.eq.1) then
        fil2=trim(adjustl(fil1))//'snap_'//iter_string
      else
        fil1=trim(adjustl(fil1))//'snapdir_'//iter_string//
     &                           '/snap_'//iter_string
        write(ifile_string,'(i7.1)') ifile
        fil2=trim(adjustl(fil1))//'.'//trim(adjustl(ifile_string))
      end if

      write(*,*) 'reading iteration file: ',iter,' ',
     &              trim(adjustl(fil2))

      call read_head(fil2,npart_gadget,massarr,bas81,bas82,
     &                  basint1,basint2,nall,basint3,basint4,bas83,
     &                  bas84,bas85,bas86,blocksize)
      write(*,*) npart_gadget(1), 'gas particles'
      low1=low2+1
      low2=low1+npart_gadget(1)-1
      
      allocate(scr42(3,sum(npart_gadget(1:6))))
      write(*,*) 'reading positions ...'
      call read_float3(fil2,'POS ',scr42,blocksize)
      write(*,*) ' found for ',(blocksize-8)/12,' particles'
      rxpa(low1:low2)=scr42(1,1:npart_gadget(1))
      rypa(low1:low2)=scr42(2,1:npart_gadget(1))
      rzpa(low1:low2)=scr42(3,1:npart_gadget(1))

      write(*,*) 'reading velocities ...'
      call read_float3(fil2,'VEL ',scr42,blocksize)
      write(*,*) ' found for ',(blocksize-8)/12,' particles'
      u2dm(low1:low2)=scr42(1,1:npart_gadget(1))
      u3dm(low1:low2)=scr42(2,1:npart_gadget(1))
      u4dm(low1:low2)=scr42(3,1:npart_gadget(1))

      deallocate(scr42)

      allocate(scr4(sum(npart_gadget(1:6))))
      write(*,*) 'reading masses ...'
      call read_float(fil2,'MASS',scr4,blocksize)
      write(*,*) ' found for ',(blocksize-8)/4,' particles'
      masap(low1:low2)=scr4(1:npart_gadget(1))

      write(*,*) 'reading kernel length ...'
      call read_float(fil2,'HSML',scr4,blocksize)
      write(*,*) ' found for ',(blocksize-8)/4,' particles'
      kernel(low1:low2)=scr4(1:npart_gadget(1)) 
      deallocate(scr4)       

#if use_filter == 1
      if (flag_filter.eq.1) then
        allocate(scr4(sum(npart_gadget(1:6))))
        if (flag_machfield.eq.0) then
          write(*,*) 'reading abvc ...'
          call read_float(fil2,'ABVC',scr4,blocksize)
          write(*,*) ' found for ',(blocksize-8)/4,' particles'
          abvc(low1:low2)=scr4(1:npart_gadget(1))  
        else 
          write(*,*) 'reading mach ...'
          call read_float(fil2,'MACH',scr4,blocksize)
          write(*,*) ' found for ',(blocksize-8)/4,' particles'
          abvc(low1:low2)=scr4(1:npart_gadget(1)) 
        end if            
        deallocate(scr4)      
      end if
#endif

#if weight_scheme == 2 || weight_filter == 2
      allocate(scr4(sum(npart_gadget(1:6))))
      write(*,*) 'reading density ...'
      call read_float(fil2,'RHO ',scr4,blocksize)
      write(*,*) ' found for ',(blocksize-8)/4,' particles'
      vol(low1:low2)=scr4(1:npart_gadget(1))
      deallocate(scr4)
#endif
#if weight_scheme == 3
      allocate(scr4(sum(npart_gadget(1:6))))
      write(*,*) 'reading density ...'
      call read_float(fil2,'RHO ',scr4,blocksize)
      write(*,*) ' found for ',(blocksize-8)/4,' particles'
      emissivity(low1:low2)=scr4(1:npart_gadget(1))
      write(*,*) 'reading internal energy ...'
      call read_float(fil2,'U   ',scr4,blocksize)
      write(*,*) ' found for ',(blocksize-8)/4,' particles'
!$omp parallel do shared(npart_gadget, scr4, emissivity, low1), 
!$omp+            private(i), default(none)
      do i=1,npart_gadget(1)
         ! calculate the weight as rho^2*sqrt(T)
         ! as this is only used as weight, we don't convert values to actual temperatures!
         emissivity(low1+i-1) = emissivity(low1+i-1)*
     &        emissivity(low1+i-1)*sqrt(scr4(i))
      end do
      deallocate(scr4)
#endif

      end do !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return 
      end 
***********************************************************************
