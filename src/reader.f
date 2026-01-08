#if reader==0
#include "readers/gadget_unformatted.f"
#endif

#if reader==1 
#include "readers/arepo_hdf5.f"
#endif

#if reader == 2 
#include "readers/masclet.f"
#endif 

#if reader == 3 
#include "readers/fixgrid_hdf5.f"
#endif

#if input_is_grid == 0
***********************************************************************
       subroutine read_particles(iter,files_per_snap,nx,ny,nz,t,zeta,
     &            nl_particle_grid,refine_thr,parchlim,borgrid,
     &            npatch,pare,patchnx,patchny,patchnz,
     &            patchx,patchy,patchz,patchrx,patchry,patchrz,lado0,
     &            npart,
     &            flag_filter,kneighbours,div_thr,abvc_thr,
     &            flag_machfield,mach_thr,flag_mass)
***********************************************************************
*     reads the gas particles of the simulation, builds a set of amr
*     grids and interpolates a continuous velocity field.
*     this subroutine may be changed for different simulation codes.
*     in particular, after mesh-building and interpolation, in any
*     case we'll need to feed the main code with:
****  grids info:
*     npatch: number of patches per refinement level
*     patchnx, patchny, patchnz: cell extensions of each refinement patch
*     patchx, patchy, patchz: grid coordinates of the leftmost cell of each patch
*     patchrx, patchry, patchrz: origin position of each patch (position of the leftmost "mother" cell)
*     pare: coarser patch a given patch is embedded in
****  clus info:
*     u2, u3, u4: initial velocity field (base level, i.e coarse grid)
*     u12, u13, u14: initial velocity field (refinement patches)
*     cr0amr: whether a cell is refined (=0) or it isn't (=1)
***********************************************************************

       use particle_data
       implicit none

       include 'vortex_parameters.dat'

       integer nx,ny,nz,iter,ndxyz,low1,low2,files_per_snap
       real t,aaa,bbb,ccc,map,zeta,lado0
       integer i,j,k,ix,nl,ir,irr,n1,n2,n3,nl_particle_grid
       integer refine_thr,parchlim,borgrid,kneighbours
       real div_thr,abvc_thr
       integer flag_machfield,flag_mass
       real mach_thr

       integer flag_filter

       integer npatch(0:nlevels),pare(npalev)
       integer npart(0:nlevels)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)

       character*200 fil1,fil2

!      this might be only for testing
       !real u1(1:nmax,1:nmay,1:nmaz)
       !real u11(1:namrx,1:namry,1:namrz,npalev)
       !common /dens/ u1,u11

       real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       common /veloc/ u2,u3,u4,u12,u13,u14

#if use_filter == 1
       integer*1 shock0(1:nmax,1:nmay,1:nmaz)
       integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
       common /shocked/ shock0,shock1
#endif 

       integer cr0amr(1:nmax,1:nmay,1:nmaz)
       integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
       common /cr0/ cr0amr, cr0amr1

*      ---parallel---
       integer num,omp_get_num_threads,numor, flag_parallel
       common /procesadores/ num

       real ddxl,ddxr,ddyl,ddyr,ddzl,ddzr
       real cio_xc,cio_yc,cio_zc
       common /dom_decomp/ ddxl,ddxr,ddyl,ddyr,ddzl,ddzr

       integer flag_verbose
       integer fl_gr_kernl,fl_gr_den,fl_gr_vel
       integer fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &         fl_gr_div,fl_gr_curl
       integer fl_p_err,fl_p_res
       integer fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb
       real fl_smooth_filtlen
       common /flags/ flag_verbose,fl_gr_kernl,fl_gr_den,fl_gr_vel,
     &        fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &        fl_gr_div,fl_gr_curl,fl_p_err,fl_p_res,
     &        fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb,
     &        fl_smooth_filtlen

       character*3 iter_string
       character*1 ifile_string
       integer ifile

       ! scratch variables for restricting the domain
       real*4,allocatable::scr42(:,:)
       integer,allocatable::elim(:)
       ! end scratch variables for restricting the domain

#if use_filter == 1
       real visc0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real visc1(namrx,namry,namrz,npalev)
#else 
       ! dummy variables
       real visc0, visc1 
#endif


       real xmin,ymin,zmin,xmax,ymax,zmax

************************************************************************
**************** change for other simulation codes *********************
************************************************************************
       ! first, get the number of particles to be read in the snapshot 
#if reader == 0
       call read_gadget_unformatted_npart(iter, files_per_snap,
     &                                    flag_filter, flag_machfield)
#endif
#if reader == 1
       call read_arepo_hdf5_npart(iter, files_per_snap,
     &                                    flag_filter, flag_machfield)
#endif
************************************************************************
************************************************************************
************************************************************************

       ! deallocate the particle arrays if they were allocated
       if (allocated(rxpa)) then 
          deallocate(rxpa,rypa,rzpa)
          deallocate(u2dm,u3dm,u4dm)
          deallocate(masap,kernel)
#if use_filter == 1
          deallocate(abvc)
#endif

#if weight_scheme == 2 || weight_filter == 2
          deallocate(vol)
#endif
#if weight_scheme == 3 || weight_filter == 2
          deallocate(emissivity)
#endif
        end if

       ! allocate the particle arrays
       allocate(rxpa(parti),rypa(parti),rzpa(parti))
       allocate(u2dm(parti),u3dm(parti),u4dm(parti))
       allocate(masap(parti),kernel(parti))
#if use_filter == 1
       allocate(abvc(parti))
#endif

#if weight_scheme == 2 || weight_filter == 2
       allocate(vol(parti))
#endif
#if weight_scheme == 3 || weight_filter == 2
       allocate(emissivity(parti))
#endif

       npart(:)=0

************************************************************************
**************** change for other simulation codes *********************
************************************************************************
       ! first, get the number of particles to be read in the snapshot 
#if reader == 0
       call read_gadget_unformatted(iter, files_per_snap,
     &                        flag_filter, flag_machfield, low2)
#endif
#if reader == 1
       call read_arepo_hdf5(iter, files_per_snap,
     &                        flag_filter, flag_machfield, low2)
#endif
************************************************************************
************************************************************************
************************************************************************

#if weight_scheme == 2 || weight_filter == 2
!$omp parallel do shared(vol, masap, parti), private(i), default(none)
       do i=1,parti
         vol(i)=masap(i)/vol(i)
      end do
#endif 
       
       npart(0)=low2 !retrocompatibility with general reader

       low1=1
       low2=sum(npart(0:nlevels))
       write(*,*)
       write(*,*) 'input particles:',low2
       xmin=minval(rxpa(low1:low2))
       xmax=maxval(rxpa(low1:low2))
       ymin=minval(rypa(low1:low2))
       ymax=maxval(rypa(low1:low2))
       zmin=minval(rzpa(low1:low2))
       zmax=maxval(rzpa(low1:low2))

       write(*,*) 'rxpa =',xmin,xmax
       write(*,*) 'rypa =',ymin,ymax
       write(*,*) 'rzpa =',zmin,zmax
       write(*,*) 'u2dm =',minval(u2dm(low1:low2)),
     &                     maxval(u2dm(low1:low2))
       write(*,*) 'u3dm =',minval(u3dm(low1:low2)),
     &                     maxval(u3dm(low1:low2))
       write(*,*) 'u4dm =',minval(u4dm(low1:low2)),
     &                     maxval(u4dm(low1:low2))
       write(*,*) 'masap=',minval(masap(low1:low2)),
     &                     maxval(masap(low1:low2))
       write(*,*) 'kernel length=',minval(kernel(low1:low2)),
     &                             maxval(kernel(low1:low2))
      
#if use_filter == 1
      if (flag_filter.eq.1) then
        if (flag_machfield.eq.0) then
         write(*,*) 'abvc=',minval(abvc(low1:low2)),
     &                      maxval(abvc(low1:low2))
        else 
         write(*,*) 'mach=',minval(abvc(low1:low2)),
     &                      maxval(abvc(low1:low2))
        end if
       end if
#endif
#if weight_scheme == 3 || weight_filter == 2
       write(*,*) 'emis=',minval(emissivity(low1:low2)),
     &      maxval(emissivity(low1:low2))
#endif

       if (xmin.lt.ddxl.or.xmax.gt.ddxr.or.
     &     ymin.lt.ddyl.or.ymax.gt.ddyr.or.
     &     zmin.lt.ddzl.or.zmax.gt.ddzr) then

        write(*,*)
        write(*,*) 'warning: particles outside the domain'
        allocate(elim(low1:low2))
        
        j=0
!$omp parallel do shared(rxpa,rypa,rzpa,low1,low2,ddxl,ddxr,ddyl,ddyr,
!$omp+            ddzl,ddzr,elim), 
!$omp+            private(i), 
!$omp+            reduction(+: j),
!$omp+            default(none)
        do i=low1,low2 
          elim(i)=0
          if (rxpa(i).lt.ddxl.or.rxpa(i).gt.ddxr.or.
     &        rypa(i).lt.ddyl.or.rypa(i).gt.ddyr.or.
     &        rzpa(i).lt.ddzl.or.rzpa(i).gt.ddzr) then
            j=j+1
            elim(i)=1
          end if
        end do

        ! new number of particles gets updated here
        parti=parti-j
        
        write(*,*) 'particles outside the domain:',j
        allocate(scr42(9,parti))

        j=0
        do i=low1,low2 
          if (elim(i).eq.0) then
            j=j+1
            scr42(1,j)=rxpa(i)
            scr42(2,j)=rypa(i)
            scr42(3,j)=rzpa(i)
            scr42(4,j)=u2dm(i)
            scr42(5,j)=u3dm(i)
            scr42(6,j)=u4dm(i)
            scr42(7,j)=masap(i)
            scr42(8,j)=kernel(i)
#if use_filter == 1
            if (flag_filter.eq.1) scr42(9,j)=abvc(i)
#endif
          end if
        end do
        deallocate(elim)

        deallocate(rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,kernel)
#if use_filter == 1
        deallocate(abvc)
#endif

        allocate(rxpa(parti),rypa(parti),rzpa(parti))
        allocate(u2dm(parti),u3dm(parti),u4dm(parti))
        allocate(masap(parti),kernel(parti))
#if use_filter == 1
        allocate(abvc(parti))
#endif

!$omp parallel do shared(scr42,rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,
!$omp+                   kernel,abvc,parti,flag_filter), 
!$omp+            private(i), 
!$omp+            default(none)
        do i=1,parti 
          rxpa(i)=scr42(1,i)
          rypa(i)=scr42(2,i)
          rzpa(i)=scr42(3,i)
          u2dm(i)=scr42(4,i)
          u3dm(i)=scr42(5,i)
          u4dm(i)=scr42(6,i)
          masap(i)=scr42(7,i)
          kernel(i)=scr42(8,i)
#if use_filter == 1
          if (flag_filter.eq.1) abvc(i)=scr42(9,i)
#endif
        end do

        deallocate(scr42)

        low2=parti
        npart(0)=parti !correct the number of particles!

       end if

       cio_xc=0.5*(ddxl+ddxr)
       cio_yc=0.5*(ddyl+ddyr)
       cio_zc=0.5*(ddzl+ddzr)
       if (abs(cio_xc).gt.1.e-6*max(ddxl,ddxr).or.
     &     abs(cio_yc).gt.1.e-6*max(ddyl,ddyr).or.
     &     abs(cio_zc).gt.1.e-6*max(ddzl,ddzr)) then

!$omp parallel do shared(rxpa,rypa,rzpa,low1,low2,cio_xc,cio_yc,cio_zc),
!$omp+            private(i), default(none)
        do i=low1,low2 
          rxpa(i)=rxpa(i)-cio_xc
          rypa(i)=rypa(i)-cio_yc
          rzpa(i)=rzpa(i)-cio_zc
        end do

       end if

       write(*,*)
       write(*,*) 'after recentering and domain decomposition',low2
       xmin=minval(rxpa(low1:low2))
       xmax=maxval(rxpa(low1:low2))
       ymin=minval(rypa(low1:low2))
       ymax=maxval(rypa(low1:low2))
       zmin=minval(rzpa(low1:low2))
       zmax=maxval(rzpa(low1:low2))

       write(*,*) 'rxpa =',xmin,xmax
       write(*,*) 'rypa =',ymin,ymax
       write(*,*) 'rzpa =',zmin,zmax
       write(*,*) 'u2dm =',minval(u2dm(low1:low2)),
     &                     maxval(u2dm(low1:low2))
       write(*,*) 'u3dm =',minval(u3dm(low1:low2)),
     &                     maxval(u3dm(low1:low2))
       write(*,*) 'u4dm =',minval(u4dm(low1:low2)),
     &                     maxval(u4dm(low1:low2))
       write(*,*) 'masap=',minval(masap(low1:low2)),
     &                     maxval(masap(low1:low2))
       write(*,*) 'kernel length=',minval(kernel(low1:low2)),
     &                             maxval(kernel(low1:low2))
       
#if use_filter == 1
       if (flag_filter.eq.1) then
        if (flag_machfield.eq.0) then
         write(*,*) 'abvc=',minval(abvc(low1:low2)),
     &                      maxval(abvc(low1:low2))
        else 
         write(*,*) 'mach=',minval(abvc(low1:low2)),
     &                      maxval(abvc(low1:low2))
        end if
       end if
#endif
#if weight_scheme == 3 || weight_filter == 2
       write(*,*) 'emis=',minval(emissivity(low1:low2)),
     &      maxval(emissivity(low1:low2))
#endif
       

       write(*,*) 'routine create mesh ------------------------------'
       npatch(0:ir)=0
       if (parchlim.gt.0) then
        call create_mesh(nx,ny,nz,nl_particle_grid,npatch,
     &            pare,patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,
     &            npart,lado0,refine_thr,parchlim,borgrid)
       else
        call create_mesh_octree(nx,ny,nz,nl_particle_grid,npatch,
     &            pare,patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,
     &            npart,lado0,refine_thr,parchlim,borgrid)
       end if
       
       nl=nl_particle_grid
       do ir=1,nl_particle_grid
        if (npatch(ir).eq.0) then 
          nl=ir-1
          exit
        end if
       end do

       call gridamr(nx,ny,nz,nl,npatch,
     &                   patchnx,patchny,patchnz,
     &                   patchx,patchy,patchz,
     &                   patchrx,patchry,patchrz,pare)
       write(*,*) 'end mesh creation --------------------------------'

       write(*,*) 'routine interpolate velocity ---------------------'
       call interpolate_velocities(nx,ny,nz,nl,npatch,pare,
     &            patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,
     &            npart,lado0,flag_filter,kneighbours,
     &            visc0,visc1,flag_machfield,flag_mass)

#if output_particles == 1
       if (fl_p_err.eq.1) then
        write(*,*) 'locating particles onto the grid'
        call place_particles(nx,ny,nz,nl,npatch,patchnx,patchny,
     &             patchnz,patchrx,patchry,patchrz,pare,
     &             npart,lado0,parchlim)

        call error_particles(nx,ny,nz,nl,npatch,patchnx,patchny,
     &             patchnz,patchrx,patchry,patchrz,pare,
     &             npart,lado0)
        deallocate(lihal, lihal_ix, lihal_jy, lihal_kz)
       end if
#endif
       write(*,*) 'end velocity interpolation -----------------------'

#if use_filter == 1
       if (flag_filter.eq.1) then 
*     all patches are extended with one extra cell per direction
        call extend_var(nx,ny,nz,nl,npatch,pare,patchnx,patchny,patchnz,
     &                  patchx,patchy,patchz,patchrx,patchry,patchrz)

        call identify_shocks(iter,nx,ny,nz,nl,npatch,pare,patchnx,
     &                       patchny,patchnz,patchrx,patchry,patchrz,
     &                       patchx,patchy,patchz,lado0,visc0,visc1,
     &                       div_thr,abvc_thr,flag_machfield,mach_thr)
       end if
#endif

!      if we do not want to output the particles, we deallocate them
!      here
#if output_particles == 0
        deallocate(rxpa,rypa,rzpa,u2dm,u3dm,u4dm,masap,kernel)
#if use_filter == 1
        deallocate(abvc)
#endif
#endif


       return
       end
#endif

#if use_filter == 1
***********************************************************************
       subroutine identify_shocks(iter,nx,ny,nz,nl,npatch,pare,patchnx,
     &                            patchny,patchnz,patchrx,patchry,
     &                            patchrz,patchx,patchy,patchz,lado0,
     &                            visc0,visc1,div_thr,abvc_thr,
     &                            flag_machfield,mach_thr)
***********************************************************************
*      for the multiscale filter, an indication of (strong) shocked
*       cells is required. this routine does that job by using a
*       combination of velocity divergence and artificial bulk 
*       viscosity constant.
***********************************************************************

      implicit none 
      include 'vortex_parameters.dat'
      integer iter,nx,ny,nz,nl 
      integer npatch(0:nlevels),pare(npalev),patchnx(npalev),
     &        patchny(npalev),patchnz(npalev)
      real patchrx(npalev),patchry(npalev),patchrz(npalev)
      integer patchx(npalev),patchy(npalev),patchz(npalev)
      real lado0
      real visc0(0:nmax+1,0:nmay+1,0:nmaz+1)
      real visc1(namrx,namry,namrz,npalev)
      real div_thr,abvc_thr
      integer flag_machfield
      real mach_thr

      real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      common /veloc/ u2,u3,u4,u12,u13,u14

      real diver0(0:nmax+1,0:nmay+1,0:nmaz+1)
      real diver(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
      common /divergence/ diver0, diver

      integer*1 shock0(1:nmax,1:nmay,1:nmaz)
      integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
      common /shocked/ shock0,shock1

      integer i,j,k,ipatch,n1,n2,n3,low1,low2,ir

      integer flag_verbose
      integer fl_gr_kernl,fl_gr_den,fl_gr_vel
      integer fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &         fl_gr_div,fl_gr_curl
      integer fl_p_err,fl_p_res
      integer fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb
      real fl_smooth_filtlen
      common /flags/ flag_verbose,fl_gr_kernl,fl_gr_den,fl_gr_vel,
     &        fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &        fl_gr_div,fl_gr_curl,fl_p_err,fl_p_res,
     &        fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb,
     &        fl_smooth_filtlen

      character*5 iter_string 
      write(iter_string, '(i5.5)') iter

      if (flag_machfield.eq.0) then
        write(*,*) 'identifying shocks with div_thr=',div_thr,
     &             ' and abvc_thr=',abvc_thr
       call diver_fina(nx,ny,nz,nl,npatch,pare,patchnx,patchny,patchnz,
     &                 patchx,patchy,patchz,patchrx,patchry,patchrz)

c       open(55,file='output_files/diver_unfiltered_'//iter_string,
c     &      status='unknown', form='unformatted')
c
c        write(55) (((diver0(i,j,k),i=1,nx),j=1,ny),k=1,nz)
c        do ipatch=1,sum(npatch(0:nl))
c         n1=patchnx(ipatch)
c         n2=patchny(ipatch)
c         n3=patchnz(ipatch)
c         write(55) (((diver(i,j,k,ipatch),i=1,n1),j=1,n2),k=1,n3)
c        end do
c
c       close(55)

       do k=1,nz 
       do j=1,ny
       do i=1,nx 
         shock0(i,j,k) = 0
         if (diver0(i,j,k).lt.div_thr.and.visc0(i,j,k).gt.abvc_thr) then
           shock0(i,j,k) = 1
         end if
       end do 
       end do 
       end do

       do ir=1,nl
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))
        do ipatch=low1,low2 
          n1=patchnx(ipatch)
          n2=patchny(ipatch)
          n3=patchnz(ipatch)
          do k=1,n3 
          do j=1,n2
          do i=1,n1 
            shock1(i,j,k,ipatch) = 0
            if (diver(i,j,k,ipatch).lt.div_thr.and.
     &          visc1(i,j,k,ipatch).gt.abvc_thr) then
              shock1(i,j,k,ipatch) = 1
            end if
          end do 
          end do 
          end do
        end do
       end do

      else 
        write(*,*) 'identifying shocks with mach_thr=',mach_thr

        do k=1,nz 
        do j=1,ny
        do i=1,nx 
          shock0(i,j,k) = 0
          if (visc0(i,j,k).gt.mach_thr) then
            shock0(i,j,k) = 1
          end if
        end do 
        end do 
        end do
  
        do ir=1,nl
          low1=sum(npatch(0:ir-1))+1
          low2=sum(npatch(0:ir))
          do ipatch=low1,low2 
            n1=patchnx(ipatch)
            n2=patchny(ipatch)
            n3=patchnz(ipatch)
            do k=1,n3 
            do j=1,n2
            do i=1,n1 
              shock1(i,j,k,ipatch) = 0
              if (visc1(i,j,k,ipatch).gt.mach_thr) then
                shock1(i,j,k,ipatch) = 1
              end if
            end do 
            end do 
            end do
          end do
        end do
      
      end if

#if output_filter == 1
      if (fl_filt_shock.eq.1) then
       call write_shocked(nx,ny,nz,iter,nl,npatch,patchnx,patchny,
     &                    patchnz,shock0,shock1)
      end if
#endif

      return 
      end 
#endif 


#if input_is_grid == 1
***********************************************************************
       subroutine read_grid(iter,files_per_snap,nx,ny,nz,t,zeta,
     &            nl_vortex,refine_thr,parchlim,borgrid,
     &            npatch,pare,patchnx,patchny,patchnz,
     &            patchx,patchy,patchz,patchrx,patchry,patchrz,lado0,
     &            npart,
     &            flag_filter,kneighbours,div_thr,abvc_thr,
     &            flag_machfield,mach_thr,flag_mass)
***********************************************************************
*     reads the gas particles of the simulation, builds a set of amr
*     grids and interpolates a continuous velocity field.
*     this subroutine may be changed for different simulation codes.
*     in particular, after mesh-building and interpolation, in any
*     case we'll need to feed the main code with:
****  grids info:
*     npatch: number of patches per refinement level
*     patchnx, patchny, patchnz: cell extensions of each refinement patch
*     patchx, patchy, patchz: grid coordinates of the leftmost cell of each patch
*     patchrx, patchry, patchrz: origin position of each patch (position of the leftmost "mother" cell)
*     pare: coarser patch a given patch is embedded in
****  clus info:
*     u2, u3, u4: initial velocity field (base level, i.e coarse grid)
*     u12, u13, u14: initial velocity field (refinement patches)
*     cr0amr: whether a cell is refined (=0) or it isn't (=1)
***********************************************************************

       use particle_data
       implicit none

       include 'vortex_parameters.dat'

       integer nx,ny,nz,iter,ndxyz,low1,low2,files_per_snap
       real t,aaa,bbb,ccc,map,zeta,lado0
       integer i,j,k,ix,nl,ir,irr,n1,n2,n3,nl_vortex
       integer refine_thr,parchlim,borgrid,kneighbours
       real div_thr,abvc_thr
       integer flag_machfield,flag_mass
       real mach_thr

       integer flag_filter

       integer npatch(0:nlevels),pare(npalev)
       integer npart(0:nlevels)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)

       character*200 fil1,fil2

       real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       common /veloc/ u2,u3,u4,u12,u13,u14

#if use_filter == 1
       integer*1 shock0(1:nmax,1:nmay,1:nmaz)
       integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
       common /shocked/ shock0,shock1
#endif 

       integer cr0amr(1:nmax,1:nmay,1:nmaz)
       integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
       common /cr0/ cr0amr, cr0amr1

*      ---parallel---
       integer num,omp_get_num_threads,numor, flag_parallel
       common /procesadores/ num

       integer flag_verbose
       integer fl_gr_kernl,fl_gr_den,fl_gr_vel
       integer fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &         fl_gr_div,fl_gr_curl
       integer fl_p_err,fl_p_res
       integer fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb
       real fl_smooth_filtlen
       common /flags/ flag_verbose,fl_gr_kernl,fl_gr_den,fl_gr_vel,
     &        fl_gr_vcomp,fl_gr_vsol,fl_gr_spot,fl_gr_vpot,
     &        fl_gr_div,fl_gr_curl,fl_p_err,fl_p_res,
     &        fl_filt_mach,fl_filt_shock,fl_filt_len,fl_filt_vturb,
     &        fl_smooth_filtlen

       character*3 iter_string
       character*1 ifile_string
       integer ifile

       ! first, get the number of particles to be read in the snapshot 
       parti = 0
       npart(:)=0

************************************************************************
**************** change for other simulation codes *********************
************************************************************************
       ! first, get the number of particles to be read in the snapshot 
#if reader == 2 
        call read_masclet_grid(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         nl_vortex)
        call read_masclet_fix_grid_level(nl, nl_vortex, npatch, 
     &            patchnx, patchny, patchnz, patchx, patchy, patchz,
     &            patchrx, patchry, patchrz, pare)
        call read_masclet_fields(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         mach_thr)
#elif reader == 3 
        call read_fixgrid_hdf5_grid(iter, flag_filter, nx, ny, nz, 
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         nl_vortex)
        call read_fixgrid_hdf5_fields(iter, flag_filter, nx, ny, nz,
     &         nl, npatch, patchnx, patchny, patchnz, 
     &         patchx, patchy, patchz,patchrx, patchry, patchrz, pare,
     &         mach_thr)
#endif
************************************************************************
************************************************************************
************************************************************************

       nl=nl_vortex
       do ir=1,nl_vortex
        if (npatch(ir).eq.0) then 
          nl=ir-1
          exit
        end if
       end do

       call compute_cr0amr(nx,ny,nz,nl,npatch,patchnx,patchny,patchnz,
     &                     patchx,patchy,patchz,patchrx,patchry,patchrz,
     &                     pare)

       call gridamr(nx,ny,nz,nl,npatch,patchnx,patchny,patchnz,
     &              patchx,patchy,patchz,patchrx,patchry,patchrz,
     &              pare)

!#if use_filter == 1
!       if (flag_filter.eq.1) then 
*     all patches are extended with one extra cell per direction
        call extend_var(nx,ny,nz,nl,npatch,pare,patchnx,patchny,patchnz,
     &                  patchx,patchy,patchz,patchrx,patchry,patchrz)
!       end if
!#endif

       return
       end
#endif

