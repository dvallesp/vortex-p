***********************************************************************
       subroutine read_masclet_grid(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         nl_vortex)
***********************************************************************
*     reads the grid information from the simulation
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

       write(iter_string, '(I5.5)') iter
       fil1 = 'simulation/grids'//iter_string

       write(*,*) 'reading grid from file:',fil1
       open(33, file=fil1, status='unknown', action='read') 

        read(33,*) ir, t, nl, map 
        read(33,*) zeta 
        read(33,*) n1,n2 ! garbage 
 
        if (ir.ne.iter) then
         stop 'error: iteration number in file and in code do not match'
        end if

        if (nl.gt.nl_vortex) then 
          write(*,*) 'Warning! The calculations will be brought up to' 
          write(*,*) ' level: ', nl_vortex
          nl = nl_vortex  
        endif

        do ir=1,nl 
          read(33,*) irr, npatch(ir)
          read(33,*)
          write(*,*) 'level',ir,'number of patches:',npatch(ir)
          low1 = sum(npatch(0:ir-1))+1
          low2 = sum(npatch(0:ir))
          
          do ip=low1,low2 
            read(33,*) patchnx(ip), patchny(ip), patchnz(ip)
            read(33,*) patchx(ip), patchy(ip), patchz(ip)
            read(33,*) patchrx(ip), patchry(ip), patchrz(ip)
            read(33,*) pare(ip)
          end do

          write(*,*) 'ir, patchx', ir, minval(patchx(low1:low2)),
     &                                 maxval(patchx(low1:low2))
          write(*,*) 'ir, patchy', ir, minval(patchy(low1:low2)),
     &                                 maxval(patchy(low1:low2))
          write(*,*) 'ir, patchz', ir, minval(patchz(low1:low2)),
     &                                 maxval(patchz(low1:low2))
          write(*,*) 'ir, patchnx', ir, minval(patchnx(low1:low2)),
     &                                  maxval(patchnx(low1:low2))
          write(*,*) 'ir, patchny', ir, minval(patchny(low1:low2)),
     &                                  maxval(patchny(low1:low2))
          write(*,*) 'ir, patchnz', ir, minval(patchnz(low1:low2)),
     &                                  maxval(patchnz(low1:low2))
          write(*,*) 'ir, patchrx', ir, minval(patchrx(low1:low2)),
     &                                  maxval(patchrx(low1:low2))
          write(*,*) 'ir, patchry', ir, minval(patchry(low1:low2)),
     &                                  maxval(patchry(low1:low2))
          write(*,*) 'ir, patchrz', ir, minval(patchrz(low1:low2)),
     &                                  maxval(patchrz(low1:low2))
          write(*,*)

        end do
       close(33)

       return 
      end

***********************************************************************
       subroutine read_masclet_fields(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         mach_thr)
***********************************************************************
*     reads the grid information from the simulation
***********************************************************************

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
       real emissivity0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real emissivity1(namrx,namry,namrz,npalev)
       common /emiss/ emissivity0,emissivity1
#endif 

       integer is_mascletB 
       is_mascletB = 0

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
!$omp parallel do shared(emissivity0,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
       do kz=0,nz+1 
       do jy=0,ny+1
       do ix=0,nx+1
        emissivity0(ix,jy,kz) = 0.0 
       end do
       end do
      end do
#endif
#endif

       low1 = 1 
       low2 = sum(npatch)

!$omp parallel do shared(u12,u13,u14,low1,low2)
!$omp+ private(ip), default(none)
       do ip = low1,low2 
        u12(:,:,:,ip) = 0.0 
        u13(:,:,:,ip) = 0.0
        u14(:,:,:,ip) = 0.0
       end do

#if use_filter == 1
!$omp parallel do shared(shock1,low1,low2)
!$omp+ private(ip), default(none)
       do ip = low1,low2 
        shock1(:,:,:,ip) = 0.0 
       end do

#if weight_filter == 1
!$omp parallel do shared(dens1,low1,low2)
!$omp+ private(ip), default(none)
       do ip = low1,low2 
        dens1(:,:,:,ip) = 0.0 
       end do
#elif weight_filter == 2
!$omp parallel do shared(emissivity1,low1,low2)
!$omp+ private(ip), default(none)
       do ip = low1,low2 
        emissivity1(:,:,:,ip) = 0.0 
       end do
#endif
#endif

       write(iter_string, '(I5.5)') iter
       fil1 = 'simulation/clus'//iter_string

       write(*,*) 'reading data from file:',fil1
       open(31, file=fil1, status='unknown', action='read', 
     &      form='unformatted') 
        read(31) ! ignore heading

        allocate(scr4(nx,ny,nz))

#if weight_filter == 1
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        dens0(1:nx,1:ny,1:nz) = 1.0 + scr4(1:nx,1:ny,1:nz)
#elif weight_filter == 2
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        emissivity0(1:nx,1:ny,1:nz) = scr4(1:nx,1:ny,1:nz)
#else
        read(31)
#endif 

        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        u2(1:nx,1:ny,1:nz) = scr4(1:nx,1:ny,1:nz)
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        u3(1:nx,1:ny,1:nz) = scr4(1:nx,1:ny,1:nz)
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        u4(1:nx,1:ny,1:nz) = scr4(1:nx,1:ny,1:nz)
        
        read(31) ! pressure
        read(31) ! pot 
        read(31) ! opot 
#if weight_filter == 2
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!$omp parallel do shared(emissivity0,scr4,nx,ny,nz), 
!$omp+ private(ix,jy,kz), default(none)
        do kz=0,nz+1 
        do jy=0,ny+1
        do ix=0,nx+1
           emissivity0(ix,jy,kz) = emissivity0(ix,jy,kz)*
     &          emissivity0(ix,jy,kz)*sqrt(scr4(ix,jy,kz)
        end do
        end do
        end do
#else
        read(31) ! temp 
#endif
        read(31) ! metal 
        read(31) ! cr0 
        if (is_mascletB.eq.1) then
         read(31) ! Bx 
         read(31) ! By
         read(31) ! Bz
        end if

        deallocate(scr4)

        low1 = 1 
        low2 = sum(npatch)
        do ip = low1, low2 
          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)
          allocate(scr4(n1,n2,n3))

#if weight_filter == 1
          read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          dens1(1:n1,1:n2,1:n3,ip) = 1.0 + scr4(1:n1,1:n2,1:n3)
#elif weight_filter == 2
          read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          emissivity1(1:n1,1:n2,1:n3,ip) = scr4(1:n1,1:n2,1:n3)
#else
          read(31)
#endif
          read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          u12(1:n1,1:n2,1:n3,ip) = scr4(1:n1,1:n2,1:n3)
          read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          u13(1:n1,1:n2,1:n3,ip) = scr4(1:n1,1:n2,1:n3)
          read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          u14(1:n1,1:n2,1:n3,ip) = scr4(1:n1,1:n2,1:n3)

          read(31) ! pressure
          read(31) ! pot
          read(31) ! opot
#if weight_filter == 2
        read(31) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!$omp parallel do shared(emissivity1,scr4,n1,n2,n3), 
!$omp+ private(ix,jy,kz), default(none)
          do kz=0,n1+1
          do jy=0,n2+1
          do ix=0,n3+1
             emissivity1(ix,jy,kz,ip) = emissivity1(ix,jy,kz,ip)*
     &            emissivity1(ix,jy,kz,ip)*sqrt(scr4(ix,jy,kz)
          end do
          end do
          end do
#else
          read(31) ! temp 
#endif
          read(31) ! temp
          read(31) ! metal
          read(31) ! cr0
          read(31) ! solap 

          if (is_mascletB.eq.1) then
           read(31) ! Bx 
           read(31) ! By
           read(31) ! Bz
          end if

          deallocate(scr4)

        end do
       close(31)

! if necessary, read machnum to tag shocked cells
#if use_filter == 1
       if (flag_filter.eq.1) then 
        fil1 = 'shocks/MachNum_'//iter_string
        write(*,*) 'reading Mach number from file:',fil1
        open(35, file=fil1, status='unknown', action='read',
     &       form='unformatted')

        allocate(scr4(nx,ny,nz))
        read(35) (((scr4(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        mach0(1:nx,1:ny,1:nz) = scr4(1:nx,1:ny,1:nz)
        deallocate(scr4)

        low1 = 1 
        low2 = sum(npatch)
        do ip = low1, low2 
          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)
          allocate(scr4(n1,n2,n3))
          read(35) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
          mach1(1:n1,1:n2,1:n3,ip) = scr4(1:n1,1:n2,1:n3)
          deallocate(scr4)
        end do
        close(35)

        shock0(1:nx,1:ny,1:nz) = 0
        where (mach0(1:nx,1:ny,1:nz).gt.mach_thr) 
         shock0(1:nx,1:ny,1:nz) = 1
        end where

        do ip = low1, low2 
          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)
          shock1(1:n1,1:n2,1:n3,ip) = 0
          where (mach1(1:n1,1:n2,1:n3,ip).gt.mach_thr) 
           shock1(1:n1,1:n2,1:n3,ip) = 1
          end where
        end do
       end if
#endif

      write(*,*) 'At level', 0
#if weight_filter == 1
      call p_minmax_ir(dens0,dens1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                 patchnz,npatch,0,basx,basy)
      write(*,*) 'dens min,max',basx,basy
#elif weight_filter == 2
      call p_minmax_ir(emis0,emis1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                 patchnz,npatch,0,basx,basy)
      write(*,*) 'emis min,max',basx,basy
#endif
      call p_minmax_ir(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                 npatch,0,basx,basy)
      write(*,*) 'vx min,max',basx,basy
      call p_minmax_ir(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                 npatch,0,basx,basy)
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
        write(*,*)
      end if
#endif


      do ir=1,nl
        write(*,*) 'At level', ir
#if weight_filter == 1
       call p_minmax_ir(dens0,dens1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch,ir,basx,basy)
       write(*,*) 'dens min,max',basx,basy
#elif weight_filter == 2
       call p_minmax_ir(emis0,emis1,1,0,nx,ny,nz,nl,patchnx,patchny,
     &                  patchnz,npatch,ir,basx,basy)
       write(*,*) 'emis min,max',basx,basy
#endif
       call p_minmax_ir(u2,u12,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,ir,basx,basy)
       write(*,*) 'vx min,max',basx,basy
       call p_minmax_ir(u3,u13,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,ir,basx,basy)
       write(*,*) 'vy min,max',basx,basy
       call p_minmax_ir(u4,u14,1,1,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                  npatch,ir,basx,basy)
       write(*,*) 'vz min,max',basx,basy
#if use_filter == 1
       if (flag_filter.eq.1) then
        call p_minmax_ir(mach0,mach1,0,0,nx,ny,nz,nl,patchnx,patchny,
     &                   patchnz,npatch,ir,basx,basy)
        write(*,*) 'mach min,max',basx,basy
        k = 0 
        low1 = sum(npatch(0:ir-1))+1
        low2 = sum(npatch(0:ir))
        do ip=low1,low2
         n1 = patchnx(ip)
         n2 = patchny(ip)
         n3 = patchnz(ip)
         k = k + sum(int(shock1(1:n1,1:n2,1:n3,ip), kind=4))
C         write(*,*) ip,k,minval(shock1(1:n1,1:n2,1:n3,ip)),
C     &              maxval(shock1(1:n1,1:n2,1:n3,ip))
        end do
        write(*,*) 'number of shocked cells:', k
        write(*,*)
       end if
#endif
      end do ! ir=1,nl

       return 
      end

************************************************************************
      subroutine read_masclet_fix_grid_level(nl, nl_vortex, npatch, 
     &            patchnx, patchny, patchnz, patchx, patchy, patchz,
     &            patchrx, patchry, patchrz, pare)
************************************************************************
*     Allows to use a coarser peak resolution for the calculations
************************************************************************
       implicit none

       include 'vortex_parameters.dat'

       integer nl, nl_vortex

       integer npatch(0:nlevels)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)
       integer pare(npalev)

       integer low1

       if (nl.gt.nl_vortex) then 
        write(*,*) 'Warning! The calculations will be brought up to' 
        write(*,*) ' level: ', nl_vortex
        nl = nl_vortex

        low1 = sum(npatch(0:nl))+1

        npatch(nl+1:) = 0

        patchnx(low1:) = 0
        patchny(low1:) = 0
        patchnz(low1:) = 0
        patchx(low1:) = 0
        patchy(low1:) = 0
        patchz(low1:) = 0
        patchrx(low1:) = 0.0
        patchry(low1:) = 0.0
        patchrz(low1:) = 0.0
        pare(low1:) = 0
       end if

       return 
      end
