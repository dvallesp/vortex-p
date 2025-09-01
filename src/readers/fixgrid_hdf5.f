***********************************************************************
       subroutine read_fixgrid_hdf5_grid(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         nl_vortex)
***********************************************************************
*      Reads the grid structure information from the simulation.
*       For the case of a fix-grid, it is just a trivial initialization.
***********************************************************************

       use particle_data
       implicit none

       include 'vortex_parameters.dat'

       integer iter, flag_filter, nx, ny, nz, nl
       integer npatch(0:nlevels)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)
       integer pare(npalev), nl_vortex

       character*200 fil1
       character*5 iter_string

       integer ir,low1,low2,i,j,k,ix,jy,kz,n1,n2,n3,ip,irr
       real t,map,zeta,bas1,bas2,bas3,bas4

       npatch(:)=0 
       patchnx(:)=0
       patchny(:)=0
       patchnz(:)=0
       patchx(:)=0
       patchy(:)=0
       patchz(:)=0
       patchrx(:)=0.0
       patchry(:)=0.0
       patchrz(:)=0.0
       pare(:)=0

       nl = 0

*      This is just needed for being able to use the python readers later on...
       call write_trivial_gridsfile(iter,nx,ny,nz,nl,npatch,
     &            patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,pare)

       return 
      end


***********************************************************************
       subroutine read_fixgrid_hdf5_fields(iter, flag_filter, nx, ny,
     &         nz, nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         mach_thr)
***********************************************************************
*      Reads the grid fields from the simulation
***********************************************************************
       use hdf5
       use particle_data
       implicit none

       include 'vortex_parameters.dat'

       integer iter, flag_filter, nx, ny, nz, nl
       integer npatch(0:nlevels)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)
       integer pare(npalev)
       real mach_thr

       character*200 fil1
       character*5 iter_string

       integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, 
     & file_space_id
       integer(hid_t) :: memtype_id
       integer :: status
       integer, dimension(6) :: numpart_thisfile
       integer(hsize_t), dimension(1) :: dims1d
       integer(hsize_t), dimension(2) :: dims2d
       integer(hsize_t), dimension(3) :: dims3d

       integer ir,low1,low2,i,j,k,ix,jy,kz,n1,n2,n3,ip,irr
       real t,map,zeta,bas1,bas2,bas3,bas4,basx,basy
       real*4, allocatable :: scr4 (:,:,:)
       !integer, allocatable :: scr4int(:,:,:)

       ! Read data 
       real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       common /veloc/ u2,u3,u4,u12,u13,u14

#if use_filter == 1
       real*4 mach0(1:nmax,1:nmay,1:nmaz)
       real*4 mach1(1:namrx,1:namry,1:namrz,npalev)
       integer*1 shock0(1:nmax,1:nmay,1:nmaz)
       integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
       common /shocked/ shock0,shock1
#endif

#if weight_filter == 1
       real dens0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real dens1(namrx,namry,namrz,npalev)
       common /densi/ dens0,dens1
#elif weight_filter == 2
       real emis0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real emis1(namrx,namry,namrz,npalev)
       common /emiss/ emis0,emis1
#endif 

       ! I assume input data is simple precision.
       real*4,allocatable :: scr3d(:,:,:)

!      Initialize the arrays before reading 

!$omp parallel do shared(u2,u3,u4,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
       do kz=0,nz+1 
       do jy=0,ny+1
       do ix=0,nx+1
        u2(ix,jy,kz) = 0.0 
        u3(ix,jy,kz) = 0.0
        u4(ix,jy,kz) = 0.0
       end do
       end do
       end do

#if use_filter == 1
!$omp parallel do shared(shock0,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
       do kz=1,nz 
       do jy=1,ny
       do ix=1,nx
        shock0(ix,jy,kz) = 0
       end do
       end do
       end do

#if weight_filter == 1
!$omp parallel do shared(dens0,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
       do kz=0,nz+1 
       do jy=0,ny+1
       do ix=0,nx+1
        dens0(ix,jy,kz) = 0.0 
       end do
       end do
       end do
#elif weight_filter == 2
!$omp parallel do shared(emis0,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
       do kz=0,nz+1 
       do jy=0,ny+1
       do ix=0,nx+1
        emis0(ix,jy,kz) = 0.0 
       end do
       end do
       end do
#endif
#endif

*      FILENAME: let's assume the filename convention is 
*        simulation/snapXXX, with XXX being the iteration number
       write(iter_string, '(I3.3)') iter
       fil1 = 'simulation/snap'//iter_string
       fil1 = trim(adjustl(fil1))//'.hdf5'
       write(*,*) 'reading data from file:',fil1

*      Open the hdf5 file in read-only mode
       call h5fopen_f(fil1, h5f_acc_rdonly_f, file_id, status)
       if (status /= 0) then
         print *, "error opening file: ", fil1
         stop
       end if

       allocate(scr3d(nx,ny,nz))
       dims3d(1) = nx
       dims3d(2) = ny
       dims3d(3) = nz

*      Read velocities (always)
       write(*,*) 'Reading velocities'

       ! file_id because the arrays are at the root level of the hdf5 file;
       ! else change it by group_id
       call h5dopen_f(file_id, "vx", attr_id, status)
       call h5dget_type_f(attr_id, memtype_id, status)
       call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
       call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
       u2(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
       call h5dclose_f(attr_id, status)

       call h5dopen_f(file_id, "vy", attr_id, status)
       call h5dget_type_f(attr_id, memtype_id, status)
       call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
       call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
       u3(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
       call h5dclose_f(attr_id, status)
       
       call h5dopen_f(file_id, "vz", attr_id, status)
       call h5dget_type_f(attr_id, memtype_id, status)
       call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
       call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
       u4(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
       call h5dclose_f(attr_id, status)

#if use_filter == 1
       if (flag_filter.eq.1) then
*        Read Mach number (if necessary)
         write(*,*) 'Reading Mach number'
         call h5dopen_f(file_id, "mach", attr_id, status)
         call h5dget_type_f(attr_id, memtype_id, status)
         call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
         call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
         mach0(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
         call h5dclose_f(attr_id, status)

#if weight_filter == 1
*        Read density (if necessary: filter + mass-weighting)
         write(*,*) 'Reading density'
         call h5dopen_f(file_id, "density", attr_id, status)
         call h5dget_type_f(attr_id, memtype_id, status)
         call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
         call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
         dens0(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
         call h5dclose_f(attr_id, status)
#elif weight_filter == 2
*        Read density and internal energy (if necessary: filter + emissivity-weighting)
         write(*,*) 'Reading density'
         call h5dopen_f(file_id, "density", attr_id, status)
         call h5dget_type_f(attr_id, memtype_id, status)
         call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
         call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
         emis0(1:nx,1:ny,1:nz) = scr3d(1:nx,1:ny,1:nz)
         call h5dclose_f(attr_id, status)
         write(*,*) 'Reading internal energy'
         call h5dopen_f(file_id, "InternalEnergy", attr_id, status)
         call h5dget_type_f(attr_id, memtype_id, status)
         call h5dread_f(attr_id, memtype_id, scr3d, dims3d, status)
         call CtoF_order_3d_real_parallel(nx,ny,nz,scr3d)
!$omp parallel do shared(emis0,scr3d,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
         do kz=0,nz+1 
         do jy=0,ny+1
         do ix=0,nx+1
            emis0(ix,jy,kz) = emis(ix,iy,iz)*emis(ix,iy,iz)*
     &           sqrt(scr3d(ix,iy,iz)) 
         end do
         end do
         end do
         call h5dclose_f(attr_id, status)
#endif

       end if
#endif


       write(*,*) 'At level', 0
       call p_minmax_ir(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,0,basx,basy)
       write(*,*) 'vx min,max',basx,basy
       call p_minmax_ir(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,0,basx,basy)
       write(*,*) 'vy min,max',basx,basy
       call p_minmax_ir(u4,u14,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,0,basx,basy)
       write(*,*) 'vz min,max',basx,basy
#if use_filter == 1
       if (flag_filter.eq.1) then
        call p_minmax_ir(mach0,mach1,0,0,nx,ny,nz,nl,patchnx,patchny,
     &                   patchnz,npatch,0,basx,basy)
        write(*,*) 'mach min,max',basx,basy
        k = sum(int(shock0(1:nx,1:ny,1:nz), kind=4))
        write(*,*) 'number of shocked cells:', k
#if weight_filter == 1
        call p_minmax_ir(dens0,dens1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                   patchnz,npatch,0,basx,basy)
        write(*,*) 'dens min,max',basx,basy
#elif weight_filter == 2
        call p_minmax_ir(emis0,emis1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                   patchnz,npatch,0,basx,basy)
        write(*,*) 'emis min,max',basx,basy
#endif
        write(*,*)
       end if
#endif

       return 
      end

***********************************************************************
      subroutine CtoF_order_3d_real_parallel(nx, ny, nz, arr)
***********************************************************************
*       Converts a 3D array from C to Fortran order, since this is 
*        not done by default in the HDF5 library.
***********************************************************************
        implicit none
    
        integer nx, ny, nz
        real arr(nx,ny,nz)

        integer i,j,k 
        real arrcopy(nx,ny,nz)

!$omp parallel do shared(arr,arrcopy,nx,ny,nz),
!$omp+ private(i,j,k), default(none)
        do k=1,nz
        do j=1,ny
        do i=1,nx
          arrcopy(i,j,k) = arr(i,j,k)
        end do
        end do
        end do

!$omp parallel do shared(arr,arrcopy,nx,ny,nz),
!$omp+ private(i,j,k), default(none)
        do k=1,nz
        do j=1,ny
        do i=1,nx
          arr(i,j,k) = arrcopy(k,j,i)
        end do
        end do
        end do

        
        return
      end
***********************************************************************


***********************************************************************
       subroutine write_trivial_gridsfile(iter,nx,ny,nz,nl,npatch,
     &            patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,pare)
***********************************************************************
*      This is needed for non-Masclet grids to use the python readers 
*      in the post-processing
***********************************************************************
        implicit none 
        include 'vortex_parameters.dat'

        integer iter,nx,ny,nz,nl
        integer npatch(0:nlevels)
        integer patchnx(npalev),patchny(npalev),patchnz(npalev)
        integer patchx(npalev),patchy(npalev),patchz(npalev)
        real patchrx(npalev),patchry(npalev),patchrz(npalev)
        integer pare(npalev)
        
        character*5 iter_string
        character*200 filenom
        integer i,ir,low1,low2

        write(iter_string, '(I5.5)') iter
        filenom = 'output_files/grids'//trim(adjustl(iter_string))

        open (23,file=filenom,status='unknown')
        write(23,*) iter,' ',0.0,' ',nl,' ',0.0,' ',0.0
        write(23,*) 0.0
        ir=0
        write(23,*) ir,0,0,nx,ny,nz
        do ir=1,nl
         write(23,*) ir,' ',npatch(ir),' ',0,' ',0,' ',0
         write(23,*) '----------------- within level=',ir,' -----------'
         low1=sum(npatch(0:ir-1))+1
         low2=sum(npatch(0:ir))
         do i=low1,low2
          write(23,*) patchnx(i),' ',patchny(i),' ',patchnz(i)
          write(23,*) patchx(i),' ',patchy(i),' ',patchz(i)
          write(23,*) patchrx(i),' ',patchry(i),' ',patchrz(i)
          write(23,*) pare(i)
         end do
        end do
       close(23)

       end
***********************************************************************
