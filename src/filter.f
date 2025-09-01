************************************************************************ 
      subroutine fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &                     patchnz,var0,var)
************************************************************************
*     Fills the gaps left by computing quantities 1 of each 27 cells 
************************************************************************
      implicit none 
      include 'vortex_parameters.dat'

      integer nx,ny,nz,bor0,borl 
      integer npatch(0:nlevels)
      integer patchnx(npalev),patchny(npalev),patchnz(npalev)
      real, intent(inout) :: var0(1-bor0:nmax+bor0,1-bor0:nmay+bor0,
     &                            1-bor0:nmaz+bor0)
      real, intent(inout) :: var(1-borl:namrx+borl,1-borl:namry+borl,
     &                           1-borl:namrz+borl,npalev)


      integer i,j,k,i1,j1,k1,i2,j2,k2,ip,n1,n2,n3,ix,jy,kz
      real bas1,bas2,bas3
      real w111,w211,w121,w221,w112,w212,w122,w222

      ! Fill gaps in the bulk by interpolation
!$omp parallel do shared(nx,ny,nz,var0),
!$omp+            private(i,j,k,i1,j1,k1,i2,j2,k2,bas1,bas2,bas3,
!$omp+                     w111,w211,w121,w221,w112,w212,w122,w222), 
!$omp+            default(none)
      do k=1,nz 
      do j=1,ny
      do i=1,nx 
        if ((mod(i,3).eq.1.or.i.eq.nx).and.
     &      (mod(j,3).eq.1.or.j.eq.ny).and.
     &      (mod(k,3).eq.1.or.k.eq.nz)) cycle

        i1 = ((i-1)/3)*3+1
        j1 = ((j-1)/3)*3+1
        k1 = ((k-1)/3)*3+1
        i2 = i1 + 3 
        j2 = j1 + 3
        k2 = k1 + 3
        if (i2.gt.nx) i2 = nx
        if (j2.gt.ny) j2 = ny
        if (k2.gt.nz) k2 = nz

        if (i1.ne.i2) then
          bas1 = float(i-i1)/float(i2-i1)
        else 
          bas1 = 0.0
        end if
        if (j1.ne.j2) then
          bas2 = float(j-j1)/float(j2-j1)
        else 
          bas2 = 0.0
        end if
        if (k1.ne.k2) then
          bas3 = float(k-k1)/float(k2-k1)
        else 
          bas3 = 0.0
        end if

        w111 = (1.0-bas1)*(1.0-bas2)*(1.0-bas3)
        w211 = bas1*(1.0-bas2)*(1.0-bas3)
        w121 = (1.0-bas1)*bas2*(1.0-bas3)
        w221 = bas1*bas2*(1.0-bas3)
        w112 = (1.0-bas1)*(1.0-bas2)*bas3
        w212 = bas1*(1.0-bas2)*bas3
        w122 = (1.0-bas1)*bas2*bas3
        w222 = bas1*bas2*bas3

        var0(i,j,k) = var0(i1,j1,k1)*w111 + 
     &                var0(i2,j1,k1)*w211 +
     &                var0(i1,j2,k1)*w121 +
     &                var0(i2,j2,k1)*w221 +
     &                var0(i1,j1,k2)*w112 +
     &                var0(i2,j1,k2)*w212 +
     &                var0(i1,j2,k2)*w122 +
     &                var0(i2,j2,k2)*w222

      end do 
      end do 
      end do

!$omp parallel do shared(npatch,patchnx,patchny,patchnz,var),
!$omp+            private(ip,n1,n2,n3,ix,jy,kz,i1,j1,k1,i2,j2,k2,
!$omp+                    bas1,bas2,bas3,w111,w211,w121,w221,w112,
!$omp+                    w212,w122,w222), 
!$omp+            default(none)
      do ip=1,sum(npatch)
        n1 = patchnx(ip)
        n2 = patchny(ip)
        n3 = patchnz(ip)

        do kz=1,n3
        do jy=1,n2 
        do ix=1,n1
          if ((mod(ix,3).eq.1.or.ix.eq.n1).and.
     &        (mod(jy,3).eq.1.or.jy.eq.n2).and.
     &        (mod(kz,3).eq.1.or.kz.eq.n3)) cycle 

          i1 = ((ix-1)/3)*3+1
          j1 = ((jy-1)/3)*3+1
          k1 = ((kz-1)/3)*3+1
          i2 = i1 + 3 
          j2 = j1 + 3
          k2 = k1 + 3
          if (i2.gt.n1) i2 = n1
          if (j2.gt.n2) j2 = n2
          if (k2.gt.n3) k2 = n3

          if (i1.ne.i2) then
            bas1 = float(ix-i1)/float(i2-i1)
          else 
            bas1 = 0.0
          end if
          if (j1.ne.j2) then
            bas2 = float(jy-j1)/float(j2-j1)
          else 
            bas2 = 0.0
          end if
          if (k1.ne.k2) then
            bas3 = float(kz-k1)/float(k2-k1)
          else 
            bas3 = 0.0
          end if

          w111 = (1.0-bas1)*(1.0-bas2)*(1.0-bas3)
          w211 = bas1*(1.0-bas2)*(1.0-bas3)
          w121 = (1.0-bas1)*bas2*(1.0-bas3)
          w221 = bas1*bas2*(1.0-bas3)
          w112 = (1.0-bas1)*(1.0-bas2)*bas3
          w212 = bas1*(1.0-bas2)*bas3
          w122 = (1.0-bas1)*bas2*bas3
          w222 = bas1*bas2*bas3
  
          var(ix,jy,kz,ip) = var(i1,j1,k1,ip)*w111 + 
     &                       var(i2,j1,k1,ip)*w211 +
     &                       var(i1,j2,k1,ip)*w121 +
     &                       var(i2,j2,k1,ip)*w221 +
     &                       var(i1,j1,k2,ip)*w112 +
     &                       var(i2,j1,k2,ip)*w212 +
     &                       var(i1,j2,k2,ip)*w122 +
     &                       var(i2,j2,k2,ip)*w222

        end do 
        end do 
        end do
      end do


      return 
      end 

************************************************************************
      subroutine compute_bulk(nx,ny,nz,nl,lado0,
     &   l,thisx,thisy,thisz,
     &   npatch,patchnx,patchny,patchnz,patchrx,patchry,patchrz,dx,
     &   solap,dens0,dens1,
     &   bulk_vx,bulk_vy,bulk_vz,marca_shock,cellcount)
************************************************************************
*     Computes the bulk velocity in a given radius
*     (used in the multiscale filtering)
************************************************************************

      implicit none
      include 'vortex_parameters.dat'

      integer nx, ny, nz, nl
      real lado0,l, thisx, thisy, thisz
      integer npatch(0:nlevels)
      integer patchnx(npalev), patchny(npalev), patchnz(npalev)
      real patchrx(npalev), patchry(npalev), patchrz(npalev)
      real dx
      integer solap(1:namrx,1:namry,1:namrz,npalev)
      real dens0(1:nmax,1:nmay,1:nmaz)
      real dens1(1:namrx,1:namry,1:namrz,npalev)
      ! Intent: out
      real bulk_vx,bulk_vy,bulk_vz
      integer marca_shock,cellcount

*     global variables
      integer cr0amr(1:nmax,1:nmay,1:nmaz)
      integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
      common /cr0/ cr0amr, cr0amr1

      real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      common /veloc/ u2,u3,u4,u12,u13,u14

      integer*1 shock0(1:nmax,1:nmay,1:nmaz)
      integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
      common /shocked/ shock0, shock1

      real rx(-2:namrx+3,npalev)
      real ry(-2:namrx+3,npalev)
      real rz(-2:namrx+3,npalev)
      real rmx(-2:namrx+3,npalev)
      real rmy(-2:namrx+3,npalev)
      real rmz(-2:namrx+3,npalev)
      common /minigrids/ rx,ry,rz,rmx,rmy,rmz

      real  radx(0:nmax+1),radmx(0:nmax+1),
     &      rady(0:nmay+1),radmy(0:nmay+1),
     &      radz(0:nmaz+1),radmz(0:nmaz+1)
      common /grid/ radx,radmx,rady,radmy,radz,radmz


*     Local variables 
      integer i,j,k,ii,jj,kk,mini,maxi,minj,maxj,mink,maxk,basint
      integer ir,ipatch,low1,low2,n1,n2,n3
      real bas1,bas2,bas3,bas4,l2,dist2_k,dist2_kj,dist2_kji
      real sph_xl,sph_xr,sph_yl,sph_yr,sph_zl,sph_zr,dxpa
      real p_xl,p_xr,p_yl,p_yr,p_zl,p_zr

      bulk_vx = 0.0
      bulk_vy = 0.0
      bulk_vz = 0.0
      marca_shock = 0
      cellcount = 0

      sph_xl = thisx - l
      sph_xr = thisx + l
      sph_yl = thisy - l
      sph_yr = thisy + l
      sph_zl = thisz - l
      sph_zr = thisz + l

      bas1 = 0.0
      bas2 = 0.0
      bas3 = 0.0
      bas4 = 0.0
      basint = 0
      l2 = l**2

      mini = int(((thisx - l) / lado0 + 0.5) * nx) + 1
      maxi = int(((thisx + l) / lado0 + 0.5) * nx) + 1
      minj = int(((thisy - l) / lado0 + 0.5) * ny) + 1
      maxj = int(((thisy + l) / lado0 + 0.5) * ny) + 1
      mink = int(((thisz - l) / lado0 + 0.5) * nz) + 1
      maxk = int(((thisz + l) / lado0 + 0.5) * nz) + 1
      if (mini.lt.1) mini=1
      if (maxi.gt.nx) maxi=nx
      if (minj.lt.1) minj=1
      if (maxj.gt.ny) maxj=ny
      if (mink.lt.1) mink=1
      if (maxk.gt.nz) maxk=nz

      outer0_c: do kk=mink,maxk
        dist2_k = (radz(kk)-thisz)**2
        if (dist2_k.gt.l2) cycle
        do jj=minj,maxj
          dist2_kj = dist2_k + (rady(jj)-thisy)**2
          if (dist2_kj.gt.l2) cycle
          do ii=mini,maxi
            if (cr0amr(ii,jj,kk).eq.0.and.nl.gt.0) cycle
            dist2_kji = dist2_kj + (radx(ii)-thisx)**2
            if (dist2_kji.gt.l2) cycle
            bas1 = bas1 + dens0(ii,jj,kk)
            bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
            bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
            bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
            basint = basint + 1
            if (shock0(ii,jj,kk).eq.1) marca_shock = 1
          end do
        end do
      end do outer0_c

      do ir = 1, nl
        low1 = sum(npatch(0:ir-1)) + 1
        low2 = sum(npatch(0:ir))
        dxpa = dx / (2.0**ir)
        do ipatch = low1, low2
          n1 = patchnx(ipatch)
          p_xl = patchrx(ipatch) - dxpa
          p_xr = p_xl + n1*dxpa
          if (p_xr.lt.sph_xl.or.p_xl.gt.sph_xr) cycle
          n2 = patchny(ipatch)
          p_yl = patchry(ipatch) - dxpa
          p_yr = p_yl + n2*dxpa
          if (p_yr.lt.sph_yl.or.p_yl.gt.sph_yr) cycle
          n3 = patchnz(ipatch)
          p_zl = patchrz(ipatch) - dxpa
          p_zr = p_zl + n3*dxpa
          if (p_zr.lt.sph_zl.or.p_zl.gt.sph_zr) cycle

          mini = floor((sph_xl-p_xl)/dxpa) + 1
          maxi = floor((sph_xr-p_xl)/dxpa) + 1
          minj = floor((sph_yl-p_yl)/dxpa) + 1
          maxj = floor((sph_yr-p_yl)/dxpa) + 1
          mink = floor((sph_zl-p_zl)/dxpa) + 1
          maxk = floor((sph_zr-p_zl)/dxpa) + 1

          if (mini.lt.1) mini=1
          if (maxi.gt.n1) maxi=n1
          if (minj.lt.1) minj=1
          if (maxj.gt.n2) maxj=n2
          if (mink.lt.1) mink=1
          if (maxk.gt.n3) maxk=n3
          if (mini.gt.n1.or.maxi.lt.1) cycle
          if (minj.gt.n2.or.maxj.lt.1) cycle
          if (mink.gt.n3.or.maxk.lt.1) cycle
          
          do kk=mink,maxk 
            dist2_k = (rz(kk,ipatch)-thisz)**2
            if (dist2_k.gt.l2) cycle
            do jj=minj,maxj
              dist2_kj = dist2_k + (ry(jj,ipatch)-thisy)**2
              if (dist2_kj.gt.l2) cycle
              do ii=mini,maxi
                if (solap(ii,jj,kk,ipatch).eq.0) cycle
                if (cr0amr1(ii,jj,kk,ipatch).eq.0.and.nl.gt.ir) cycle
                dist2_kji = dist2_kj + (rx(ii,ipatch)-thisx)**2
                if (dist2_kji.gt.l2) cycle
                bas1=bas1 + dens1(ii,jj,kk,ipatch)
                bas2=bas2 + dens1(ii,jj,kk,ipatch)*u12(ii,jj,kk,ipatch)
                bas3=bas3 + dens1(ii,jj,kk,ipatch)*u13(ii,jj,kk,ipatch)
                bas4=bas4 + dens1(ii,jj,kk,ipatch)*u14(ii,jj,kk,ipatch)
                basint = basint + 1
                if (shock1(ii,jj,kk,ipatch).eq.1) marca_shock = 1
              end do
            end do
          end do

        end do ! ipatch 
      end do ! ir

      if (basint.ne.0) then 
        bulk_vx = bas2 / bas1
        bulk_vy = bas3 / bas1
        bulk_vz = bas4 / bas1
      else
        bulk_vx = 0.0
        bulk_vy = 0.0
        bulk_vz = 0.0
      end if 
      cellcount = basint

      return 
      end 

************************************************************************ 
      subroutine kern_smooth(r,h,kval)
************************************************************************ 
      implicit none

      real r,h,kval
      real q

      q = r/h

      if (q.lt.0.5) then
        kval = 1.0 - 6.0*q*q + 6.0*q*q*q
      else if (q.lt.1.0) then
        kval = 2.0*(1.0-q)**3
      else 
        kval = 0.0
      end if

      return
      end
************************************************************************ 


************************************************************************ 
      subroutine smooth_l(nx,ny,nz,nl,npatch,patchnx,patchny,
     &                     patchnz,patchrx,patchry,patchrz,l0,l1,solap,
     &                     fl_smooth_filtlen)
************************************************************************
      implicit none 
      include 'vortex_parameters.dat'

      integer nx,ny,nz,nl,bor0,borl 
      integer npatch(0:nlevels)
      integer patchnx(npalev),patchny(npalev),patchnz(npalev)
      real patchrx(npalev),patchry(npalev),patchrz(npalev)
      real l0(1:nmax,1:nmay,1:nmaz)
      real l1(1:namrx,1:namry,1:namrz,npalev)
      integer solap(1:namrx,1:namry,1:namrz,npalev)
      real fl_smooth_filtlen


      real rx(-2:namrx+3,npalev)
      real ry(-2:namrx+3,npalev)
      real rz(-2:namrx+3,npalev)
      real rmx(-2:namrx+3,npalev)
      real rmy(-2:namrx+3,npalev)
      real rmz(-2:namrx+3,npalev)
      common /minigrids/ rx,ry,rz,rmx,rmy,rmz

      real  radx(0:nmax+1),radmx(0:nmax+1),
     &        rady(0:nmay+1),radmy(0:nmay+1),
     &        radz(0:nmaz+1),radmz(0:nmaz+1)
      common /grid/   radx,radmx,rady,radmy,radz,radmz

      integer cr0amr(1:nmax,1:nmay,1:nmaz)
      integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
      common /cr0/ cr0amr, cr0amr1

      real dx,dy,dz 
      common /espaciado/ dx,dy,dz

      real l0new(1:nmax,1:nmay,1:nmaz)
      real l1new(1:namrx,1:namry,1:namrz,npalev)
      real w0(1:nmax,1:nmay,1:nmaz)
      real w1(1:namrx,1:namry,1:namrz,npalev)

      
      integer i,j,k,ip,low1,low2,ir,n1,n2,n3,ii,jj,kk,jp,nn1,nn2,nn3 
      integer irr,mini,maxi,minj,maxj,mink,maxk,low3,low4
      real dxpa,dypa,dzpa,val,basevol,vol,l,l2,xl,yl,zl,xr,yr,zr
      real dista,kval,x1,y1,z1,x2,y2,z2,xxl,xxr,yyl,yyr,zzl,zzr
      real xldom,yldom,zldom,bas1,bas2,exp1,exp2,bas

      basevol = dx*dy*dz

!$omp parallel do shared(nx,ny,nz,l0new,w0),
!$omp+            private(i,j,k), default(none)
      do k=1,nz 
      do j=1,ny 
      do i=1,nx 
        l0new(i,j,k) = 0.0
        w0(i,j,k) = 0.0
      end do 
      end do 
      end do

!$omp parallel do shared(npatch,patchnx,patchny,patchnz,l1new,w1),
!$omp+            private(ip,n1,n2,n3,i,j,k),
!$omp+            default(none)
      do ip=1,sum(npatch)
        n1 = patchnx(ip)
        n2 = patchny(ip)
        n3 = patchnz(ip)
        do k=1,n3
        do j=1,n2 
        do i=1,n1 
          l1new(i,j,k,ip) = 0.0
          w1(i,j,k,ip) = 0.0
        end do 
        end do 
        end do
      end do ! ip

      xldom = -(nx/2)*dx
      yldom = -(ny/2)*dy
      zldom = -(nz/2)*dz

!$omp parallel do shared(nx,ny,nz,cr0amr,l0,basevol,radx,rady,radz,dx,
!$omp+                   dy,dz,nl,npatch,patchnx,patchny,patchnz,
!$omp+                   patchrx,patchry,patchrz,rx,ry,rz,xldom,yldom,
!$omp+                   zldom,l0new,w0,l1new,w1),
!$omp+            private(i,j,k,l,l2,val,vol,x1,y1,z1,xl,xr,yl,yr,zl,zr,
!$omp+                    mini,maxi,minj,maxj,mink,maxk,ii,jj,kk,x2,y2,
!$omp+                    z2,dista,kval,irr,low1,low2,dxpa,jp,nn1,nn2,
!$omp+                    nn3,xxl,xxr,yyl,yyr,zzl,zzr),
CX !$omp+            reduction(+:l0new,w0,l1new,w1),
!$omp+            default(none)
      do k=1,nz 
      do j=1,ny 
      do i=1,nx 
        if (cr0amr(i,j,k).eq.0) cycle
        l = l0(i,j,k)
        l2 = l**2
        val = log(l)
        vol = basevol 

        x1 = radx(i)
        y1 = rady(j)
        z1 = radz(k)

        xl = x1 - l
        xr = x1 + l
        yl = y1 - l
        yr = y1 + l
        zl = z1 - l
        zr = z1 + l
        mini = int((xl-xldom)/dx) + 1
        maxi = int((xr-xldom)/dx) + 1
        minj = int((yl-yldom)/dy) + 1
        maxj = int((yr-yldom)/dy) + 1
        mink = int((zl-zldom)/dz) + 1
        maxk = int((zr-zldom)/dz) + 1
        if (mini.lt.1) mini=1
        if (maxi.gt.nx) maxi=nx
        if (minj.lt.1) minj=1
        if (maxj.gt.ny) maxj=ny
        if (mink.lt.1) mink=1
        if (maxk.gt.nz) maxk=nz

        do kk=mink,maxk 
          if (.not.(mod(kk,3).eq.1.or.kk.eq.nz)) cycle
        do jj=minj,maxj
          if (.not.(mod(jj,3).eq.1.or.jj.eq.ny)) cycle
        do ii=mini,maxi 
          if (.not.(mod(ii,3).eq.1.or.ii.eq.nx)) cycle
          x2 = radx(ii)
          y2 = rady(jj)
          z2 = radz(kk)
          dista = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
          if (dista.gt.l2) cycle 
          dista = sqrt(dista)
          call kern_smooth(dista,l,kval)

!$omp atomic
          l0new(ii,jj,kk) = l0new(ii,jj,kk) + kval*vol*val 
!$omp atomic
          w0(ii,jj,kk) = w0(ii,jj,kk) + kval*vol
        end do 
        end do 
        end do 

        do irr=1,nl
          low1 = sum(npatch(0:irr-1)) + 1
          low2 = sum(npatch(0:irr))
          dxpa = dx / (2.0**irr)
          do jp = low1,low2 
            nn1 = patchnx(jp)
            xxl = patchrx(jp) - dxpa
            xxr = xxl + nn1*dxpa
            if (xxr.lt.xl.or.xxl.gt.xr) cycle

            nn2 = patchny(jp)
            yyl = patchry(jp) - dxpa
            yyr = yyl + nn2*dxpa
            if (yyr.lt.yl.or.yyl.gt.yr) cycle

            nn3 = patchnz(jp)
            zzl = patchrz(jp) - dxpa
            zzr = zzl + nn3*dxpa
            if (zzr.lt.zl.or.zzl.gt.zr) cycle

            mini = floor((xl-xxl)/dxpa) + 1
            maxi = floor((xr-xxl)/dxpa) + 1
            minj = floor((yl-yyl)/dxpa) + 1
            maxj = floor((yr-yyl)/dxpa) + 1
            mink = floor((zl-zzl)/dxpa) + 1
            maxk = floor((zr-zzl)/dxpa) + 1
            if (mini.lt.1) mini=1
            if (maxi.gt.nn1) maxi=nn1
            if (minj.lt.1) minj=1
            if (maxj.gt.nn2) maxj=nn2
            if (mink.lt.1) mink=1
            if (maxk.gt.nn3) maxk=nn3

            do kk=mink,maxk
              if (.not.(mod(kk,3).eq.1.or.kk.eq.nn3)) cycle
            do jj=minj,maxj
              if (.not.(mod(jj,3).eq.1.or.jj.eq.nn2)) cycle
            do ii=mini,maxi
              if (.not.(mod(ii,3).eq.1.or.ii.eq.nn1)) cycle
              x2 = rx(ii,jp)
              y2 = ry(jj,jp)
              z2 = rz(kk,jp)
              dista = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
              if (dista.gt.l2) cycle 
              dista = sqrt(dista)
              call kern_smooth(dista,l,kval)

!$omp atomic
              l1new(ii,jj,kk,jp) = l1new(ii,jj,kk,jp) + kval*vol*val 
!$omp atomic
              w1(ii,jj,kk,jp) = w1(ii,jj,kk,jp) + kval*vol
            end do 
            end do
            end do
          end do
        end do
        
      end do 
      end do 
      end do

      ! Now the same for the cells coming from refinement patches
      do ir=1,nl
        low3 = sum(npatch(0:ir-1)) + 1
        low4 = sum(npatch(0:ir))
        vol = basevol / (2.0**ir)**3
!$omp parallel do shared(low3,low4,patchnx,patchny,patchnz,cr0amr1,
!$omp+                   solap,l1,rx,ry,rz,dx,dy,dz,radx,rady,radz,vol,
!$omp+                   nl,npatch,patchrx,patchry,patchrz,xldom,yldom,
!$omp+                   zldom,nx,ny,nz,l0new,w0,l1new,w1),
!$omp+            private(ip,n1,n2,n3,i,j,k,l,l2,val,x1,y1,z1,xl,xr,yl,
!$omp+                    yr,zl,zr,mini,maxi,minj,maxj,mink,maxk,ii,jj,
!$omp+                    kk,x2,y2,z2,dista,kval,irr,low1,low2,dxpa,
!$omp+                    jp,nn1,nn2,nn3,xxl,xxr,yyl,yyr,zzl,zzr),
CX !$omp+            reduction(+:l0new,w0,l1new,w1),
!$omp+            default(none)
        do ip = low3, low4 
          !write(*,*) 'smooth ip=',ip
          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)
          do k=1,n3
          do j=1,n2 
          do i=1,n1 
            if (cr0amr1(i,j,k,ip).eq.0) cycle
            if (solap(i,j,k,ip).eq.0) cycle
            l = l1(i,j,k,ip)
            l2 = l**2
            val = log(l)
    
            x1 = rx(i,ip)
            y1 = ry(j,ip)
            z1 = rz(k,ip)
    
            xl = x1 - l
            xr = x1 + l
            yl = y1 - l
            yr = y1 + l
            zl = z1 - l
            zr = z1 + l
            mini = int((xl-xldom)/dx) + 1
            maxi = int((xr-xldom)/dx) + 1
            minj = int((yl-yldom)/dy) + 1
            maxj = int((yr-yldom)/dy) + 1
            mink = int((zl-zldom)/dz) + 1
            maxk = int((zr-zldom)/dz) + 1
            if (mini.lt.1) mini=1
            if (maxi.gt.n1) maxi=nx
            if (minj.lt.1) minj=1
            if (maxj.gt.n2) maxj=ny
            if (mink.lt.1) mink=1
            if (maxk.gt.n3) maxk=nz
    
            do kk=mink,maxk
              if (.not.(mod(kk,3).eq.1.or.kk.eq.nz)) cycle 
            do jj=minj,maxj
              if (.not.(mod(jj,3).eq.1.or.jj.eq.ny)) cycle
            do ii=mini,maxi 
              if (.not.(mod(ii,3).eq.1.or.ii.eq.nx)) cycle
              x2 = radx(ii)
              y2 = rady(jj)
              z2 = radz(kk)
              dista = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
              if (dista.gt.l2) cycle 
              dista = sqrt(dista)
              call kern_smooth(dista,l,kval)
    
!$omp atomic
              l0new(ii,jj,kk) = l0new(ii,jj,kk) + kval*vol*val 
!$omp atomic
              w0(ii,jj,kk) = w0(ii,jj,kk) + kval*vol
            end do 
            end do 
            end do 
    
            do irr=1,nl
              low1 = sum(npatch(0:irr-1)) + 1
              low2 = sum(npatch(0:irr))
              dxpa = dx / (2.0**irr)
              do jp = low1,low2 
                nn1 = patchnx(jp)
                xxl = patchrx(jp) - dxpa
                xxr = xxl + nn1*dxpa
                if (xxr.lt.xl.or.xxl.gt.xr) cycle
    
                nn2 = patchny(jp)
                yyl = patchry(jp) - dxpa
                yyr = yyl + nn2*dxpa
                if (yyr.lt.yl.or.yyl.gt.yr) cycle
    
                nn3 = patchnz(jp)
                zzl = patchrz(jp) - dxpa
                zzr = zzl + nn3*dxpa
                if (zzr.lt.zl.or.zzl.gt.zr) cycle
    
                mini = floor((xl-xxl)/dxpa) + 1
                maxi = floor((xr-xxl)/dxpa) + 1
                minj = floor((yl-yyl)/dxpa) + 1
                maxj = floor((yr-yyl)/dxpa) + 1
                mink = floor((zl-zzl)/dxpa) + 1
                maxk = floor((zr-zzl)/dxpa) + 1
                if (mini.lt.1) mini=1
                if (maxi.gt.nn1) maxi=nn1
                if (minj.lt.1) minj=1
                if (maxj.gt.nn2) maxj=nn2
                if (mink.lt.1) mink=1
                if (maxk.gt.nn3) maxk=nn3
    
                do kk=mink,maxk
                  if (.not.(mod(kk,3).eq.1.or.kk.eq.nn3)) cycle
                do jj=minj,maxj
                  if (.not.(mod(jj,3).eq.1.or.jj.eq.nn2)) cycle
                do ii=mini,maxi
                  if (.not.(mod(ii,3).eq.1.or.ii.eq.nn1)) cycle
                  x2 = rx(ii,jp)
                  y2 = ry(jj,jp)
                  z2 = rz(kk,jp)
                  dista = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                  if (dista.gt.l2) cycle 
                  dista = sqrt(dista)
                  call kern_smooth(dista,l,kval)
    
!$omp atomic
                  l1new(ii,jj,kk,jp) = l1new(ii,jj,kk,jp) + kval*vol*val 
!$omp atomic
                  w1(ii,jj,kk,jp) = w1(ii,jj,kk,jp) + kval*vol
                end do 
                end do
                end do
              end do
            end do
            
          end do 
          end do 
          end do
        end do 
      end do

      ! Now we normalize the values
!$omp parallel do shared(nx,ny,nz,w0,l0new),
!$omp+            private(i,j,k),
!$omp+            default(none)
      do k=1,nz
        if (.not.(mod(k,3).eq.1.or.k.eq.nz)) cycle
      do j=1,ny
        if (.not.(mod(j,3).eq.1.or.j.eq.ny)) cycle
      do i=1,nx
        if (.not.(mod(i,3).eq.1.or.i.eq.nx)) cycle
        l0new(i,j,k) = l0new(i,j,k) / w0(i,j,k)
      end do 
      end do 
      end do

!$omp parallel do shared(npatch,patchnx,patchny,patchnz,w1,l1new),
!$omp+            private(i,j,k,n1,n2,n3,ip),
!$omp+            default(none)
      do ip=1,sum(npatch)
        n1 = patchnx(ip)
        n2 = patchny(ip)
        n3 = patchnz(ip)
        do k=1,n3
          if (.not.(mod(k,3).eq.1.or.k.eq.n3)) cycle
        do j=1,n2 
          if (.not.(mod(j,3).eq.1.or.j.eq.n2)) cycle
        do i=1,n1 
          if (.not.(mod(i,3).eq.1.or.i.eq.n1)) cycle
          l1new(i,j,k,ip) = l1new(i,j,k,ip) / w1(i,j,k,ip)
        end do 
        end do 
        end do
      end do ! ip

      bor0=0
      borl=0
      call fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &               patchnz,l0new,l1new)

      ! Now we average geometrically, if necessary, with the previous value
      exp1 = fl_smooth_filtlen
      if (exp1.lt.0.0) exp1 = 0.0
      if (exp1.gt.1.0) exp1 = 1.0
      exp2 = 1 - exp1 

!$omp parallel do shared(nx,ny,nz,l0,l0new,dx,exp1,exp2),
!$omp+            private(i,j,k,bas1,bas2,bas),
!$omp+            default(none)
      do k=1,nz
      do j=1,ny
      do i=1,nx
        bas1 = exp(l0new(i,j,k))
        bas2 = l0(i,j,k) ! this was not log 
        bas = (bas1**exp1) * (bas2**exp2)
        l0(i,j,k) = max(bas, 4.*dx)
      end do 
      end do 
      end do

      do ir=1,nl
        low1 = sum(npatch(0:ir-1)) + 1
        low2 = sum(npatch(0:ir))
        dxpa = dx / (2.0**ir)
!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,l1,dxpa,
!$omp+                   l1new,exp1,exp2),
!$omp+            private(i,j,k,n1,n2,n3,ip,bas,bas1,bas2),
!$omp+            default(none)
        do ip=low1,low2
          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)
          do k=1,n3
          do j=1,n2 
          do i=1,n1 
            bas1 = exp(l1new(i,j,k,ip))
            bas2 = l1(i,j,k,ip) ! this was not log
            bas = (bas1**exp1) * (bas2**exp2)
            l1(i,j,k,ip) = max(bas, 4.*dxpa)
          end do 
          end do 
          end do
        end do ! ip
      end do ! ir


      return 
      end

************************************************************************
      subroutine compute_bulk_given_l(nx,ny,nz,nl,npatch,patchnx,
     &       patchny,patchnz,patchrx,patchry,patchrz,l0,l1,dens0,dens1,
     &       solap,u2bulk,u3bulk,u4bulk,u12bulk,u13bulk,u14bulk)
************************************************************************

      implicit none
      include 'vortex_parameters.dat'

      integer nx,ny,nz,nl
      integer npatch(0:nlevels)
      integer patchnx(npalev),patchny(npalev),patchnz(npalev)
      real patchrx(npalev),patchry(npalev),patchrz(npalev)
      real l0(1:nmax,1:nmay,1:nmaz)
      real l1(1:namrx,1:namry,1:namrz,npalev)
      real dens0(1:nmax,1:nmay,1:nmaz)
      real dens1(1:namrx,1:namry,1:namrz,npalev)
      integer solap(1:namrx,1:namry,1:namrz,npalev)
      real u2bulk(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real u3bulk(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real u4bulk(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real u12bulk(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real u13bulk(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real u14bulk(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)

      real rx(-2:namrx+3,npalev)
      real ry(-2:namrx+3,npalev)
      real rz(-2:namrx+3,npalev)
      real rmx(-2:namrx+3,npalev)
      real rmy(-2:namrx+3,npalev)
      real rmz(-2:namrx+3,npalev)
      common /minigrids/ rx,ry,rz,rmx,rmy,rmz

      real  radx(0:nmax+1),radmx(0:nmax+1),
     &        rady(0:nmay+1),radmy(0:nmay+1),
     &        radz(0:nmaz+1),radmz(0:nmaz+1)
      common /grid/   radx,radmx,rady,radmy,radz,radmz

      real dx,dy,dz 
      common /espaciado/ dx,dy,dz

      integer i,j,k,i3,j3,k3,marca_shock,basint,ipatch,n1,n2,n3,ir
      integer bor0,borl
      real l,thisx,thisy,thisz,lado0,bas2,bas3,bas4


      lado0 = nx * dx

!$omp parallel do shared(nx,ny,nz,radx,rady,radz,lado0,
!$omp+                   npatch,patchnx,patchny,patchnz,patchrx,
!$omp+                   patchry,patchrz,dx,solap,dens0,dens1,
!$omp+                   u2bulk,u3bulk,u4bulk,l0),
!$omp+            private(i3,j3,k3,i,j,k,l,thisx,thisy,thisz,bas2,bas3,
!$omp+                    bas4,marca_shock,basint),
!$omp+            default(none), schedule(dynamic)
      do k3=1,nz+2,3
        k = k3 
        if (k3.gt.nz) k = nz
        do j3=1,ny+2,3
          j = j3
          if (j3.gt.ny) j = ny
          do i3=1,nx+2,3
            i = i3
            if (i3.gt.nx) i = nx

            l = l0(i,j,k)
            thisx = radx(i)
            thisy = rady(j)
            thisz = radz(k)

            call compute_bulk(nx,ny,nz,0,lado0,l,thisx,thisy,thisz,
     &   npatch,patchnx,patchny,patchnz,patchrx,patchry,patchrz,dx,
     &   solap,dens0,dens1,bas2,bas3,bas4,marca_shock,basint)
      
            u2bulk(i,j,k) = bas2
            u3bulk(i,j,k) = bas3
            u4bulk(i,j,k) = bas4
          end do
        end do
      end do ! do i=1,nx

!$omp parallel do shared(npatch,patchnx,patchny,patchnz,nl,
!$omp+                   rx,ry,rz,nx,ny,nz,lado0,patchrx,patchry,
!$omp+                   patchrz,dx,solap,dens0,dens1,u12bulk,u13bulk,
!$omp+                   u14bulk,l1), 
!$omp+            private(ipatch,n1,n2,n3,ir,i3,j3,k3,i,j,k,l,thisx,
!$omp+                    thisy,thisz,bas2,bas3,bas4,marca_shock,
!$omp+                    basint),
!$omp+            default(none), schedule(dynamic)
      do ipatch=1,sum(npatch)
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)

        do ir=1,nl
          if (ipatch.le.sum(npatch(0:ir))) exit
        end do

        do k3=1,n3+2,3
          k=k3
          if (k3.gt.n3) k = n3
        do j3=1,n2+2,3
          j=j3
          if (j3.gt.n2) j = n2
        do i3=1,n1+2,3
          i=i3
          if (i3.gt.n1) i = n1

          l = l1(i,j,k,ipatch)
          thisx = rx(i,ipatch)
          thisy = ry(j,ipatch)
          thisz = rz(k,ipatch)

          call compute_bulk(nx,ny,nz,ir,lado0,l,thisx,thisy,thisz,
     &   npatch,patchnx,patchny,patchnz,patchrx,patchry,patchrz,dx,
     &   solap,dens0,dens1,bas2,bas3,bas4,marca_shock,basint)

          u12bulk(i,j,k,ipatch) = bas2
          u13bulk(i,j,k,ipatch) = bas3
          u14bulk(i,j,k,ipatch) = bas4
        end do
        end do
        end do
      end do

      bor0 = 1
      borl = 1
      call fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &               patchnz,u2bulk,u12bulk)    
      call fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &               patchnz,u3bulk,u13bulk)   
      call fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &               patchnz,u4bulk,u14bulk) 

      return 
      end

************************************************************************
      subroutine multiscale_filter(nx,ny,nz,nl,npatch,pare,
     &            patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,dx,output_iter,
     &            tol,step,maxit,maxlength,flag_filter,
     &            fl_smooth_filtlen)
************************************************************************
*     Implements the multiscale filtering technique described in
*     Vazza et al. 2012 to an AMR grid (instead of a fixed grid).
************************************************************************

      implicit none

      include 'vortex_parameters.dat'

*     function parameters
      integer nx, ny, nz, nl
      integer npatch(0:nlevels), pare(npalev)
      integer patchnx(npalev), patchny(npalev), patchnz(npalev)
      integer patchx(npalev), patchy(npalev), patchz(npalev)
      real patchrx(npalev), patchry(npalev), patchrz(npalev)
      real dx
      integer output_iter, maxit
      real tol, step, maxlength
      integer flag_filter
      real fl_smooth_filtlen

*     global variables
*     original velocity: at the end, the filtered velocity will be
*     stored here
      real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      common /veloc/ u2, u3, u4, u12, u13, u14

*     filtering scales
      real l0(1:nmax,1:nmay,1:nmaz)
      real l1(1:namrx,1:namry,1:namrz,npalev)

      real rx(-2:namrx+3,npalev)
      real ry(-2:namrx+3,npalev)
      real rz(-2:namrx+3,npalev)
      real rmx(-2:namrx+3,npalev)
      real rmy(-2:namrx+3,npalev)
      real rmz(-2:namrx+3,npalev)
      common /minigrids/ rx,ry,rz,rmx,rmy,rmz

      real  radx(0:nmax+1),radmx(0:nmax+1),
     &        rady(0:nmay+1),radmy(0:nmay+1),
     &        radz(0:nmaz+1),radmz(0:nmaz+1)
      common /grid/   radx,radmx,rady,radmy,radz,radmz

      real dens0(1:nmax,1:nmay,1:nmaz)
      real dens1(1:namrx,1:namry,1:namrz,npalev)
      !common /dens/ dens0,dens1

      integer cr0amr(1:nmax,1:nmay,1:nmaz)
      integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
      common /cr0/ cr0amr, cr0amr1

      integer solap(1:namrx,1:namry,1:namrz,npalev)
      integer*1 shock0(1:nmax,1:nmay,1:nmaz)
      integer*1 shock1(1:namrx,1:namry,1:namrz,npalev)
      common /shocked/ shock0,shock1

      integer cr3amr1(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
      integer cr3amr1x(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
      integer cr3amr1y(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
      integer cr3amr1z(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
      common /cr0cell/ cr3amr1,cr3amr1x,cr3amr1y,cr3amr1z

#if weight_filter == 1
      REAL DENSI_IN0(0:nmax+1,0:nmay+1,0:nmaz+1)
      REAL DENSI_IN1(NAMRX,NAMRY,NAMRZ,NPALEV)
      COMMON /DENSI/ DENSI_IN0,DENSI_IN1
#elif weight_filter == 2
      REAL EMISS_IN0(0:nmax+1,0:nmay+1,0:nmaz+1)
      REAL EMISS_IN1(NAMRX,NAMRY,NAMRZ,NPALEV)
      COMMON /EMISS/ EMISS_IN0,EMISS_IN1
#endif 

*     Auxiliary variables
      real u2bulk(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u3bulk(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u4bulk(0:nmax+1,0:nmay+1,0:nmaz+1)
      real u12bulk(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u13bulk(0:namrx+1,0:namry+1,0:namrz+1,npalev)
      real u14bulk(0:namrx+1,0:namry+1,0:namrz+1,npalev)

*     private VARIABLES
      integer ix, jy, kz, i, j, k, low1, low2, n1, n2, n3, ir, ipatch
      integer i3, j3, k3, i1, i2, j1, j2, k1, k2, ip, ii, jj, kk
      integer marca, iter, basint, basintprev, marca_shock
      integer jpatch, nn1, nn2, nn3, cr1, cr2, cr3, bor0, borl
      integer do_smooth
      real bas1, bas2, bas3, bas4, dxpa, l, err
      real thisx, thisy, thisz, dv2, dv3, dv4, dv2prev, dv3prev, dv4prev
      real lado0,w111,w112,w121,w122,w211,w212,w221,w222
      real xc, yc, zc, dxpa_j, xl, yl, zl, xr, yr, zr
      real u(2,2,2), fuin, uw(2,2,2)
      real ubas(3,3,3), rxbas(3), rybas(3), rzbas(3)

      real basx,basy

      integer exectime, time

      do_smooth = 0
      if (fl_smooth_filtlen.gt.0.001) do_smooth = 1

      lado0 = nx * dx

*     DENS0, DENS1 PROPORTIONAL TO CELLS MASSES or VOLUME!!!
      dens0 = 1.0
      do ir=1,nl
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,dens1,ir),
!$omp+            private(i), default(none)
       do i=low1,low2
         dens1(:,:,:,i) = 1 / 8.0**ir
       end do
      end do

#if weight_filter == 1
!$omp parallel do shared(nx,ny,nz,densi_in0,dens0),
!$omp+            private(i,j,k), 
!$omp+            default(none)
      do k=1,nz 
      do j=1,ny 
      do i=1,nx 
        dens0(i,j,k) = densi_in0(i,j,k) * dens0(i,j,k)
      end do 
      end do 
      end do

      do ir=1,nl 
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,densi_in1,dens1,ir),
!$omp+            private(i), default(none)
        do i=low1,low2
          dens1(:,:,:,i) = densi_in1(:,:,:,i) * dens1(:,:,:,i) 
        end do
      end do
#elif weight_filter == 2
!$omp parallel do shared(nx,ny,nz,emiss_in0,dens0),
!$omp+            private(i,j,k), 
!$omp+            default(none)
      do k=1,nz 
      do j=1,ny 
      do i=1,nx 
        dens0(i,j,k) = emiss_in0(i,j,k) * dens0(i,j,k)
      end do 
      end do 
      end do

      do ir=1,nl 
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,emiss_in1,dens1,ir),
!$omp+            private(i), default(none)
        do i=low1,low2
          dens1(:,:,:,i) = emiss_in1(:,:,:,i) * dens1(:,:,:,i) 
        end do
      end do
#endif 

      write(*,*) 'in filter'
      call p_minmax(dens0,dens1,0,0,nx,ny,nz,nl,
     &                    patchnx,patchny,patchnz,npatch)

      call veinsgrid_all_l(nl,npatch,pare,patchnx,patchny,patchnz,
     &            patchx,patchy,patchz,patchrx,patchry,patchrz,solap)


      if (flag_filter.eq.1) then !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!! for each cell, we ought to find the optimum coherence length
      ! we first initialize the lengths
      l0 = 0.!3.0 * dx
      do ir=1,nl
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
       dxpa = dx / (2.0**ir)
!$omp parallel do shared(low1,low2,l1,dxpa),
!$omp+            private(i), default(none)
       do i=low1,low2
         l1(:,:,:,i) = 0.!3.0 * dxpa
       end do
      end do

************
* Sect. 1.1. Here we compute L for the base grid 
************

!$omp parallel do shared(nx,ny,nz,l0,radx,rady,radz,lado0,cr0amr,
!$omp+                   dens0,u2,u3,u4,nl,npatch,patchnx,patchny,
!$omp+                   patchnz,patchrx,patchry,patchrz,cr0amr1,solap,
!$omp+                   rx,ry,rz,dens1,u12,u13,u14,tol,u2bulk,u3bulk,
!$omp+                   u4bulk,step,maxit,dx,shock0,shock1,maxlength),
!$omp+            private(i,j,k,marca,iter,l,thisx,thisy,thisz,bas1,
!$omp+                    bas2,bas3,bas4,err,
!$omp+                    dv2,dv3,dv4,dv2prev,dv3prev,dv4prev,
!$omp+                    basint,basintprev,i3,j3,k3,marca_shock),
!$omp+            schedule(dynamic), default(none)
      do k3=1,nz+2,3
        k = k3 
        if (k3.gt.nz) k = nz
        do j3=1,ny+2,3
          j = j3
          if (j3.gt.ny) j = ny
          do i3=1,nx+2,3
            i = i3
            if (i3.gt.nx) i = nx

            iter = 0
            L = 3.0*dx
            thisx = radx(i)
            thisy = rady(j)
            thisz = radz(k)

            basintprev = 0
            marca = 1
            iter_while_c: do while (marca.eq.1)

              !!! 1.COMPUTE THE BULK VELOCITY at l=0
              call compute_bulk(nx,ny,nz,0,lado0,l,thisx,thisy,thisz,
     &   npatch,patchnx,patchny,patchnz,patchrx,patchry,patchrz,dx,
     &   solap,dens0,dens1,bas2,bas3,bas4,marca_shock,basint)
      
              ! If there's a shock, and we are not at the first iteration,
              !  then we stop the iteration straight away
              if (marca_shock.eq.1.and.basintprev.gt.10) then 
                marca=0 
                exit iter_while_c 
              end if

              ! At least, 10 more cells than in the previous. If not,
              !  grow even more without checking convergence
              if (basint-basintprev.lt.10.and.iter.ne.0) then 
               l = max(l*step, l+dx)
               iter = iter + 1 
               cycle iter_while_c
              else
               basintprev=basint
              end if

              !!! 2. ERROR CALCULATION
              if (iter.eq.0) then
                err = 2.0 * tol
              else
                dv2prev = u2(i,j,k) - u2bulk(i,j,k)
                dv3prev = u3(i,j,k) - u3bulk(i,j,k)
                dv4prev = u4(i,j,k) - u4bulk(i,j,k)
                dv2 = u2(i,j,k) - bas2
                dv3 = u3(i,j,k) - bas3
                dv4 = u4(i,j,k) - bas4
                err = max(abs((dv2/dv2prev-1.0)), abs(dv3/dv3prev-1.0),
     &                    abs(dv4/dv4prev-1.0))
              end if

              ! Save the new bulk velocity
              u2bulk(i,j,k) = bas2
              u3bulk(i,j,k) = bas3
              u4bulk(i,j,k) = bas4

*             stop condition: min error or max num of its
              if (err.le.tol.or.iter.gt.maxit) marca = 0
*             stop when the growing sphere touches the domain boundary
*             (so we do not worry about boundary condition)
              if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

*             the sphere grows
              if (marca.eq.1) l = max(l*step, l+dx)
              iter = iter+1

              ! Maximum length of the coherence length, to prevent 
              !  pathological cases
              if (l.gt.maxlength) then 
                l = maxlength
                marca = 0
              end if
            end do iter_while_c
            
            l0(i,j,k) = l
          end do
        end do
      end do ! do i=1,nx

************
* Sect. 1.2. Here we compute L for the refinement levels
************
      low1=1
      low2=sum(npatch)
!$omp parallel do shared(nx,ny,nz,l1,radx,rady,radz,lado0,cr0amr,
!$omp+                   dens0,u2,u3,u4,nl,npatch,patchnx,patchny,
!$omp+                   patchnz,patchrx,patchry,patchrz,cr0amr1,solap,
!$omp+                   rx,ry,rz,dens1,u12,u13,u14,tol,u12bulk,u13bulk,
!$omp+                   u14bulk,step,maxit,dx,shock0,shock1,low1,low2,
!$omp+                   maxlength),
!$omp+            private(i,j,k,marca,iter,l,thisx,thisy,thisz,bas1,
!$omp+                    bas2,bas3,bas4,err,
!$omp+                    dv2,dv3,dv4,dv2prev,dv3prev,dv4prev,
!$omp+                    ipatch,n1,n2,n3,exectime,dxpa,ir,
!$omp+                    basint,basintprev,i3,j3,k3,marca_shock),
!$omp+            schedule(dynamic), default(none)
      do ipatch=low1,low2
        exectime = time()
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        ! get the level
        do ir=1,nl
          if (ipatch.le.sum(npatch(0:ir))) exit
        end do
        !write(*,*) 'starting',ir,ipatch
        dxpa = dx/(2.0**ir)

        do k3=1,n3+2,3
         k=k3
         if (k3.gt.n3) k = n3
        do j3=1,n2+2,3
          j=j3
          if (j3.gt.n2) j = n2
        do i3=1,n1+2,3
          i=i3
          if (i3.gt.n1) i = n1

          iter = 0
          L = 3.0*dxpa
          thisx = rx(i,ipatch)
          thisy = ry(j,ipatch)
          thisz = rz(k,ipatch)

          basintprev = 0
          marca = 1
          iter_while: do while (marca.eq.1)

            !!! 1.COMPUTE THE BULK VELOCITY at l=0
            call compute_bulk(nx,ny,nz,ir,lado0,l,thisx,thisy,thisz,
     &   npatch,patchnx,patchny,patchnz,patchrx,patchry,patchrz,dx,
     &   solap,dens0,dens1,bas2,bas3,bas4,marca_shock,basint)

            ! If there's a shock, and we are not at the first iteration,
            !  then we stop the iteration straight away
            if (marca_shock.eq.1.and.basintprev.gt.10) then 
              marca=0 
              exit iter_while
            end if
              
            ! At least, 10 more cells than in the previous. If not,
            !  grow even more without checking convergence
            if (basint-basintprev.lt.10.and.iter.ne.0) then 
             l = max(l*step, l+dxpa)
             iter = iter + 1
             cycle iter_while
            else
             basintprev=basint
            end if

            !!! 2. ERROR CALCULATION
            if (iter.eq.0) then
              err = 2.0 * tol
            else
              dv2prev = u12(i,j,k,ipatch) - u12bulk(i,j,k,ipatch)
              dv3prev = u13(i,j,k,ipatch) - u13bulk(i,j,k,ipatch)
              dv4prev = u14(i,j,k,ipatch) - u14bulk(i,j,k,ipatch)
              dv2 = u12(i,j,k,ipatch) - bas2
              dv3 = u13(i,j,k,ipatch) - bas3
              dv4 = u14(i,j,k,ipatch) - bas4
              err = max(abs((dv2/dv2prev-1.0)), abs(dv3/dv3prev-1.0),
     &                    abs(dv4/dv4prev-1.0))
            end if

            ! Save the new bulk velocity
            u12bulk(i,j,k,ipatch) = bas2
            u13bulk(i,j,k,ipatch) = bas3
            u14bulk(i,j,k,ipatch) = bas4

*             stop condition: min error or max num of its
            if (err.le.tol.or.iter.gt.maxit) marca = 0
*             stop when the growing sphere touches the domain boundary
*             (so we do not worry about boundary condition)
            if (l.gt.min(thisx+0.5*lado0,0.5*lado0-thisx,
     &                     thisy+0.5*lado0,0.5*lado0-thisy,
     &                     thisz+0.5*lado0,0.5*lado0-thisz)) marca=0

*           the sphere grows
            if (marca.eq.1) l = max(l*step, l+dxpa)
            iter = iter+1

            ! Maximum length of the coherence length, to prevent
            !  pathological cases
            if (l.gt.maxlength) then 
              l = maxlength
              marca = 0
            end if
          end do iter_while
          
          l1(i,j,k,ipatch) = l
        end do
        end do
        end do ! do i=1,n1
      end do

      write(*,*) 'refinement levels done!'

************
* Sect. 2. Here we fill the gaps in L0, L1 by interpolation
************
      bor0 = 0
      borl = 0
      call fill_gaps(nx,ny,nz,bor0,borl,npatch,patchnx,patchny,
     &               patchnz,l0,l1)

************
* Sect. 3. Here we smooth the coherence length
************
      if (do_smooth.eq.1) then !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call smooth_l(nx,ny,nz,nl,npatch,patchnx,patchny,
     &              patchnz,patchrx,patchry,patchrz,l0,l1,solap,
     &              fl_smooth_filtlen)
      end if

      else if (flag_filter.eq.2) then !~ fix filtering length ~~~~~~~~~~~~~
************
* Sect. 1, 2, and 3. We just need to copy the value of L0 and L1
************

!$omp parallel do shared(nx,ny,nz,maxlength,l0),
!$omp+            private(i,j,k),
!$omp+            default(none)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            l0(i,j,k) = maxlength
          end do
        end do
      end do ! do i=1,nx

      low1 = 1 
      low2 = sum(npatch)
!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,l1,
!$omp+                   maxlength), 
!$omp+            private(ipatch,n1,n2,n3,i,j,k),
!$omp+            default(none)
      do ipatch = low1, low2
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        do k=1,n3
        do j=1,n2
        do i=1,n1
          l1(i,j,k,ipatch) = maxlength
        end do
        end do
        end do
      end do
      end if !~~~~~~~~~~~~~~~~~ (flag_filter.eq.1 or .eq.2) ~~~~~~~~~~~~~~~~

****************************************************************************
* Sect. 4. Now we compute the bulk velocity from the pre-filled L0 and L1
       call compute_bulk_given_l(nx,ny,nz,nl,npatch,patchnx,
     &       patchny,patchnz,patchrx,patchry,patchrz,l0,l1,dens0,dens1,
     &       solap,u2bulk,u3bulk,u4bulk,u12bulk,u13bulk,u14bulk)
****************************************************************************

************
* Sect. 5. Now we will the ghost cells in the base grid and the patches
************
*     Base grid: periodic 
      do kz=0,nz+1 
      do jy=0,ny+1
      do ix=0,nx+1 
        if (ix.eq.0.or.ix.eq.nx+1.or.
     &      jy.eq.0.or.jy.eq.ny+1.or.
     &      kz.eq.0.or.kz.eq.nz+1) then 
          ii = ix 
          jj = jy
          kk = kz
          if (ii.lt.1) ii = ii + nx
          if (ii.gt.nx) ii = ii - nx
          if (jj.lt.1) jj = jj + ny
          if (jj.gt.ny) jj = jj - ny
          if (kk.lt.1) kk = kk + nz
          if (kk.gt.nz) kk = kk - nz

          u2bulk(ix,jy,kz) = u2bulk(ii,jj,kk)
          u3bulk(ix,jy,kz) = u3bulk(ii,jj,kk)
          u4bulk(ix,jy,kz) = u4bulk(ii,jj,kk)
        end if
      end do 
      end do 
      end do

*     Refinement patches: look for neighbours 
      do ir=1,nl 
        low1 = sum(npatch(0:ir-1))+1
        low2 = sum(npatch(0:ir))
        dxpa = dx/(2.0**ir)

        do ipatch=low1,low2
          n1 = patchnx(ipatch)
          n2 = patchny(ipatch)
          n3 = patchnz(ipatch)

          do kz=0,n3+1 
          do jy=0,n2+1
          do ix=0,n1+1 
            if (ix.eq.0.or.ix.eq.n1+1.or.
     &          jy.eq.0.or.jy.eq.n2+1.or.
     &          kz.eq.0.or.kz.eq.n3+1) then 
              xc = rx(ix,ipatch)
              yc = ry(jy,ipatch)
              zc = rz(kz,ipatch)

              marca = 0
              dxpa_j = dxpa
              do jpatch=low1,low2 
                nn1 = patchnx(jpatch)
                nn2 = patchny(jpatch)
                nn3 = patchnz(jpatch)
                xl = patchrx(jpatch) - dxpa_j
                yl = patchry(jpatch) - dxpa_j
                zl = patchrz(jpatch) - dxpa_j
                xr = xl + nn1*dxpa_j
                yr = yl + nn2*dxpa_j
                zr = zl + nn3*dxpa_j

                if (xc.ge.xl.and.xc.le.xr.and.
     &              yc.ge.yl.and.yc.le.yr.and.
     &              zc.ge.zl.and.zc.le.zr) then 
                  ii = int((xc-xl)/dxpa_j) + 1
                  jj = int((yc-yl)/dxpa_j) + 1
                  kk = int((zc-zl)/dxpa_j) + 1

                  u12bulk(ix,jy,kz,ipatch) = u12bulk(ii,jj,kk,jpatch)
                  u13bulk(ix,jy,kz,ipatch) = u13bulk(ii,jj,kk,jpatch)
                  u14bulk(ix,jy,kz,ipatch) = u14bulk(ii,jj,kk,jpatch)

                  marca = 1 
                  !write(*,*) 'copied from sibling'
                  exit
                end if
              end do

              if (marca.eq.1) cycle 

              ! If we are here, it means that we have not found a patch
              ! at the same level. We need to go to lower levels and 
              ! interpolate...
              jpatch = cr3amr1(ix,jy,kz,ipatch)
              cr1 = cr3amr1x(ix,jy,kz,ipatch)
              cr2 = cr3amr1y(ix,jy,kz,ipatch)
              cr3 = cr3amr1z(ix,jy,kz,ipatch)

              if (jpatch.gt.0) then 
                rxbas = rx(cr1-1:cr1+1, jpatch)
                rybas = ry(cr2-1:cr2+1, jpatch)
                rzbas = rz(cr3-1:cr3+1, jpatch)

                !write(*,*) 'interpolating from fine',jpatch,cr1,cr2,cr3
                !if (xc.lt.rxbas(1)) write(*,*) 'wfx-',xc,rxbas,cr1
                !if (xc.gt.rxbas(3)) write(*,*) 'wfx+',xc,rxbas,cr1
                !if (yc.lt.rybas(1)) write(*,*) 'wfy-',yc,rybas,cr2
                !if (yc.gt.rybas(3)) write(*,*) 'wfy+',yc,rybas,cr2
                !if (zc.lt.rzbas(1)) write(*,*) 'wfz-',zc,rzbas,cr3
                !if (zc.gt.rzbas(3)) write(*,*) 'wfz+',zc,rzbas,cr3

                ! U12
                ubas(1:3, 1:3, 1:3) = u12bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1, jpatch)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u12bulk(ix,jy,kz,ipatch) = fuin

                ! U13
                ubas(1:3, 1:3, 1:3) = u13bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1, jpatch)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u13bulk(ix,jy,kz,ipatch) = fuin

                ! U14
                ubas(1:3, 1:3, 1:3) = u14bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1, jpatch)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u14bulk(ix,jy,kz,ipatch) = fuin
                !write(*,*) 'interpolated from coarser, not base'
              else ! base grid 
                rxbas = radx(cr1-1:cr1+1)
                rybas = rady(cr2-1:cr2+1)
                rzbas = radz(cr3-1:cr3+1)

                !if (xc.lt.rxbas(1)) write(*,*) 'wcx-',xc,rxbas,cr1
                !if (xc.gt.rxbas(3)) write(*,*) 'wcx+',xc,rxbas,cr1
                !if (yc.lt.rybas(1)) write(*,*) 'wcy-',yc,rybas,cr2
                !if (yc.gt.rybas(3)) write(*,*) 'wcy+',yc,rybas,cr2
                !if (zc.lt.rzbas(1)) write(*,*) 'wcz-',zc,rzbas,cr3
                !if (zc.gt.rzbas(3)) write(*,*) 'wcz+',zc,rzbas,cr3

                ! U12
                ubas(1:3, 1:3, 1:3) = u2bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u12bulk(ix,jy,kz,ipatch) = fuin

                ! U13
                ubas(1:3, 1:3, 1:3) = u3bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u13bulk(ix,jy,kz,ipatch) = fuin

                ! U14
                ubas(1:3, 1:3, 1:3) = u4bulk(cr1-1:cr1+1, cr2-1:cr2+1,
     &                                   cr3-1:cr3+1)
                call linint52d_new_real(xc,yc,zc,rxbas,rybas,rzbas,ubas,
     &                                  fuin)
                u14bulk(ix,jy,kz,ipatch) = fuin  
                !write(*,*) 'interpolated from base grid'              
              end if
            end if
          end do 
          end do 
          end do
        end do
      end do


      write(*,*) '------------------------'
      write(*,*) 'after filter'
      do ir=0,nl 
        call p_minmax_ir(l0,l1,0,0,nx,ny,nz,nl,patchnx,patchny,patchnz,
     &                 npatch,ir,basx,basy)
        write(*,*) 'L ir,min,max',ir,basx,basy
        call p_minmax_ir(abs(u2bulk),abs(u12bulk),0,0,nx,ny,nz,nl,
     &                 patchnx,patchny,patchnz,npatch,ir,basx,basy)
        write(*,*) 'abs(u2bulk) ir,min,max',ir,basx,basy
        call p_minmax_ir(abs(u3bulk),abs(u13bulk),0,0,nx,ny,nz,nl,
     &                 patchnx,patchny,patchnz,npatch,ir,basx,basy)
        write(*,*) 'abs(u3bulk) ir,min,max',ir,basx,basy
        call p_minmax_ir(abs(u4bulk),abs(u14bulk),0,0,nx,ny,nz,nl,
     &                 patchnx,patchny,patchnz,npatch,ir,basx,basy)
        write(*,*) 'abs(u4bulk) ir,min,max',ir,basx,basy
      end do

************
* Sect. 6. Now we compute the velocity fluctuation / turbulent velocities
************
      ! U2,U3,U4,U12,U13,U14 gets updated with the values of the
      ! velocity fluctuation
      low1=sum(npatch)
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,U12,U13,U14,
!$OMP+                   U12BULK,U13BULK,U14BULK,LOW1),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,I,J,K),
!$OMP+            DEFAULT(NONE)
      do ipatch=1,low1
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        do k=0,n3+1
        do j=0,n2+1
        do i=0,n1+1
          u12(i,j,k,ipatch)=u12(i,j,k,ipatch)-u12bulk(i,j,k,ipatch)
          u13(i,j,k,ipatch)=u13(i,j,k,ipatch)-u13bulk(i,j,k,ipatch)
          u14(i,j,k,ipatch)=u14(i,j,k,ipatch)-u14bulk(i,j,k,ipatch)
        end do
        end do
        end do
      end do

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,U2BULK,U3BULK,U4BULK),
!$OMP+            PRIVATE(I,J,K), DEFAULT(NONE)      
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            u2(i,j,k) = u2(i,j,k) - u2bulk(i,j,k)
            u3(i,j,k) = u3(i,j,k) - u3bulk(i,j,k)
            u4(i,j,k) = u4(i,j,k) - u4bulk(i,j,k)
          end do
        end do
      end do

#if output_filter==1
      call write_filtlen(nx,ny,nz,nl,npatch,
     &                   patchnx,patchny,patchnz,l0,l1)
#endif

      RETURN
      END

************************************************************************
      SUBROUTINE SMOOTH_BULK(U2BULK,U3BULK,U4BULK,L0,U12BULK,U13BULK,
     &                       U14BULK,L1,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                       PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                       PATCHRY,PATCHRZ,NL)
************************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      real u2bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u3bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u4bulk(1:NMAX,1:NMAY,1:NMAZ)
      real l0(1:NMAX,1:NMAY,1:NMAZ)
      real u12bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u13bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u14bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real l1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      integer nx,ny,nz
      integer npatch(0:NLEVELS)
      integer patchnx(NPALEV), patchny(NPALEV), patchnz(NPALEV)
      integer patchx(NPALEV), patchy(NPALEV), patchz(NPALEV)
      real patchrx(NPALEV), patchry(NPALEV), patchrz(NPALEV)
      integer nl
      real fl_smooth_filtlen

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

      real RX(-2:NAMRX+3,NPALEV)
      real RY(-2:NAMRX+3,NPALEV)
      real RZ(-2:NAMRX+3,NPALEV)
      real RMX(-2:NAMRX+3,NPALEV)
      real RMY(-2:NAMRX+3,NPALEV)
      real RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      integer i,j,k,ix,jy,kz,ir,ii,jj,kk,jpatch,low1,low2,ipatch
      integer kr1,kr2,kr3,kare,n1,n2,n3
      real,allocatable::arr1(:,:,:),arr2(:,:,:),arr3(:,:,:),arr4(:,:,:)
      real ubas(3,3,3),rxbas(3),rybas(3),rzbas(3),basx,basy,basz,fuin
      real gauss(-3:3,-3:3,-3:3)
      real*8 bas8,bas81,bas82,bas83,bas84

*     Initialise the gaussian filter
      bas8 = 0.d0
      do k=-3,3
      do j=-3,3
      do i=-3,3
        basx=exp(-0.5*(i**2+j**2+k**2))
        gauss(i,j,k)=basx 
        bas8=bas8+dble(basx)
      end do
      end do
      end do

      do k=-3,3
      do j=-3,3
      do i=-3,3
        gauss(i,j,k)=gauss(i,j,k)/bas8
      end do
      end do
      end do


      do ir=nl,1,-1 
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))

!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,u12bulk,
!$omp+                   u13bulk,u14bulk,l1,rx,ry,rz,cr3amr1,cr3amr1x,
!$omp+                   cr3amr1y,cr3amr1z,u2bulk,u3bulk,u4bulk,l0,
!$omp+                   radx,rady,radz,gauss),
!$omp+            private(ipatch,n1,n2,n3,arr1,arr2,arr3,arr4,i,j,k,
!$omp+                    kare,kr1,kr2,kr3,basx,basy,basz,rxbas,rybas,
!$omp+                    rzbas,ubas,fuin,ii,jj,kk,ix,jy,kz,bas81,
!$omp+                    bas82,bas83,bas84),
!$omp+            default(none), schedule(dynamic)
        do ipatch=low1,low2 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          n1=patchnx(ipatch)
          n2=patchny(ipatch)
          n3=patchnz(ipatch)
          allocate(arr1(-2:n1+3,-2:n2+3,-2:n3+3))
          allocate(arr2(-2:n1+3,-2:n2+3,-2:n3+3))
          allocate(arr3(-2:n1+3,-2:n2+3,-2:n3+3))
          allocate(arr4(-2:n1+3,-2:n2+3,-2:n3+3))

*         Backup internal values and innitialise the ghost ones         
          do k=-2,n3+3
          do j=-2,n2+3
          do i=-2,n1+3
            if (i.gt.0.and.i.le.n1.and.
     &          j.gt.0.and.j.le.n2.and.
     &          k.gt.0.and.k.le.n3) then
              arr1(i,j,k) = u12bulk(i,j,k,ipatch)
              arr2(i,j,k) = u13bulk(i,j,k,ipatch)
              arr3(i,j,k) = u14bulk(i,j,k,ipatch)
              arr4(i,j,k) = l1(i,j,k,ipatch)
            else
              kare=cr3amr1(i,j,k,ipatch)
              kr1=cr3amr1x(i,j,k,ipatch)
              kr2=cr3amr1y(i,j,k,ipatch)
              kr3=cr3amr1z(i,j,k,ipatch)

              basx=rx(i,ipatch)
              basy=ry(j,ipatch)
              basz=rz(k,ipatch)

              if (kare.gt.0) then 
                rxbas=rx(kr1-1:kr1+1,kare)
                rybas=ry(kr2-1:kr2+1,kare)
                rzbas=rz(kr3-1:kr3+1,kare)

*               U2
                ubas(1:3,1:3,1:3)=u12bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                    kr3-1:kr3+1,kare)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr1(i,j,k)=fuin

*               U3
                ubas(1:3,1:3,1:3)=u13bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                    kr3-1:kr3+1,kare)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr2(i,j,k)=fuin  
                
*               U4
                ubas(1:3,1:3,1:3)=u14bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                    kr3-1:kr3+1,kare)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr3(i,j,k)=fuin  

*               L
                ubas(1:3,1:3,1:3)=l1(kr1-1:kr1+1,kr2-1:kr2+1,
     &                               kr3-1:kr3+1,kare)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                 rzbas,ubas,fuin)
                arr4(i,j,k)=fuin  
              else 
                rxbas=radx(kr1-1:kr1+1)
                rybas=rady(kr2-1:kr2+1)
                rzbas=radz(kr3-1:kr3+1)

*               U2
                ubas(1:3,1:3,1:3)=u2bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                   kr3-1:kr3+1)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr1(i,j,k)=fuin

*               U3
                ubas(1:3,1:3,1:3)=u3bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                    kr3-1:kr3+1)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr2(i,j,k)=fuin  
                
*               U4
                ubas(1:3,1:3,1:3)=u4bulk(kr1-1:kr1+1,kr2-1:kr2+1,
     &                                    kr3-1:kr3+1)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr3(i,j,k)=fuin  

*               L
                ubas(1:3,1:3,1:3)=l0(kr1-1:kr1+1,kr2-1:kr2+1,
     &                               kr3-1:kr3+1)
                call linint52d_new_real(basx,basy,basz,rxbas,rybas,
     &                                  rzbas,ubas,fuin)
                arr4(i,j,k)=fuin 
              end if 
            end if
          end do
          end do
          end do 

*         Smooth the bulk values
          do k=1,n3
          do j=1,n2
          do i=1,n1
            u12bulk(i,j,k,ipatch)=0.0
            u13bulk(i,j,k,ipatch)=0.0
            u14bulk(i,j,k,ipatch)=0.0
            l1(i,j,k,ipatch)=0.0

            bas81=0.d0
            bas82=0.d0
            bas83=0.d0
            bas84=0.d0
            do kz=-3,3
            do jy=-3,3
            do ix=-3,3
              ii=i+ix
              jj=j+jy
              kk=k+kz
              
              bas81=bas81+dble(gauss(ix,jy,kz)*arr1(ii,jj,kk))
              bas82=bas82+dble(gauss(ix,jy,kz)*arr2(ii,jj,kk))
              bas83=bas83+dble(gauss(ix,jy,kz)*arr3(ii,jj,kk))
              bas84=bas84+dble(gauss(ix,jy,kz)*arr4(ii,jj,kk))
            end do 
            end do 
            end do

            u12bulk(i,j,k,ipatch)=bas81
            u13bulk(i,j,k,ipatch)=bas82
            u14bulk(i,j,k,ipatch)=bas83
            l1(i,j,k,ipatch)=bas84

          end do
          end do
          end do

          deallocate(arr1,arr2,arr3,arr4)
        end do !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end do 

      allocate(arr1(-2:nx+3,-2:ny+3,-2:nz+3))
      allocate(arr2(-2:nx+3,-2:ny+3,-2:nz+3))
      allocate(arr3(-2:nx+3,-2:ny+3,-2:nz+3))
      allocate(arr4(-2:nx+3,-2:ny+3,-2:nz+3))

!$omp parallel do shared(u2bulk,u3bulk,u4bulk,l0,arr1,arr2,arr3,arr4,
!$omp+                   nx,ny,nz,gauss), 
!$omp+            private(i,j,k,ix,jy,kz), 
!$omp+            default(none)
      do k=-2,nz+3
      do j=-2,ny+3
      do i=-2,nx+3
        if (i.gt.0.and.i.le.nx.and.
     &      j.gt.0.and.j.le.ny.and.
     &      k.gt.0.and.k.le.nz) then
          arr1(i,j,k) = u2bulk(i,j,k)
          arr2(i,j,k) = u3bulk(i,j,k)
          arr3(i,j,k) = u4bulk(i,j,k)
          arr4(i,j,k) = l0(i,j,k)
        else
          ix = i 
          jy = j
          kz = k
          if (ix.lt.1) ix=ix+nx 
          if (ix.gt.nx) ix=ix-nx
          if (jy.lt.1) jy=jy+ny
          if (jy.gt.ny) jy=jy-ny
          if (kz.lt.1) kz=kz+nz
          if (kz.gt.nz) kz=kz-nz
          arr1(i,j,k) = u2bulk(ix,jy,kz)
          arr2(i,j,k) = u3bulk(ix,jy,kz)
          arr3(i,j,k) = u4bulk(ix,jy,kz)
          arr4(i,j,k) = l0(ix,jy,kz)
        end if
      end do
      end do
      end do

!$omp parallel do shared(u2bulk,u3bulk,u4bulk,l0,arr1,arr2,arr3,arr4,
!$omp+                   nx,ny,nz,gauss),
!$omp+            private(i,j,k,bas81,bas82,bas83,bas84,ii,jj,kk,ix,jy,
!$omp+                    kz),
!$omp+            default(none)
      do k=1,nz
      do j=1,ny
      do i=1,nx
        u2bulk(i,j,k)=0.0
        u3bulk(i,j,k)=0.0
        u4bulk(i,j,k)=0.0
        l0(i,j,k)=0.0

        bas81=0.d0
        bas82=0.d0
        bas83=0.d0
        bas84=0.d0
        do kz=-3,3
        do jy=-3,3
        do ix=-3,3
          ii=i+ix
          jj=j+jy
          kk=k+kz
          
          bas81=bas81+dble(gauss(ix,jy,kz)*arr1(ii,jj,kk))
          bas82=bas82+dble(gauss(ix,jy,kz)*arr2(ii,jj,kk))
          bas83=bas83+dble(gauss(ix,jy,kz)*arr3(ii,jj,kk))
          bas84=bas84+dble(gauss(ix,jy,kz)*arr4(ii,jj,kk))
        end do 
        end do 
        end do

        u2bulk(i,j,k)=bas81
        u3bulk(i,j,k)=bas82
        u4bulk(i,j,k)=bas83
        l0(i,j,k)=bas84
      end do
      end do 
      end do 

      deallocate(arr1,arr2,arr3,arr4)


      RETURN 
      END
