************************************************************************
      SUBROUTINE MULTISCALE_FILTER(NX,NY,NZ,NL,NPATCH,pare,
     &            PATCHNX,PATCHNY,PATCHNZ,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz,DX,output_iter,
     &            tol,step,maxit)
************************************************************************
*     Implements the multiscale filtering technique described in
*     Vazza et al. 2012 to an AMR grid (instead of a fixed grid).
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX, NY, NZ, NL
      INTEGER NPATCH(0:NLEVELS), pare(npalev)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      real PATCHrX(NPALEV),PATCHrY(NPALEV),PATCHrZ(NPALEV)
      REAL DX
      INTEGER output_iter, maxit
      real tol, step

*     global variables
*     original velocity: at the end, the filtered velocity will be
*     stored here
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

*     filtering scales
      real L0(1:NMAX,1:NMAY,1:NMAZ)
      real L1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

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

      real dens0(1:NMAX,1:NMAY,1:NMAZ)
      real dens1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      !common /dens/ dens0,dens1

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

      integer solap(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

      INTEGER*1 SHOCK0(1:NMAX,1:NMAY,1:NMAZ)
      INTEGER*1 SHOCK1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      COMMON /SHOCKED/ SHOCK0,SHOCK1

*     Auxiliary variables
      real u2bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u3bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u4bulk(1:NMAX,1:NMAY,1:NMAZ)
      real u12bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u13bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      real u14bulk(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

*     private VARIABLES
      INTEGER IX, JY, KZ, I, J, K, LOW1, LOW2, N1, N2, N3, IR, irr
      integer ii, jj, kk, iixx, jjyy, kkzz, ipatch, jpatch, llow1, llow2
      integer nn1, nn2, nn3
      integer marca, iter, basint, basintprev
      real BAS1, BAS2, BAS3, bas4, DXPA, L, l2, err, dxpa_i
      real thisx, thisy, thisz, dv2, dv3, dv4, dv2prev, dv3prev, dv4prev
      real lado0
      integer mini, maxi, minj, maxj, mink, maxk
      real rx1, rx2, ry1, ry2, rz1, rz2, rxx1, rxx2, ryy1, ryy2,
     &     rzz1, rzz2
      real u(2,2,2), fuin, uw(2,2,2)
      character*13 filenom
      character*30 filerr5

      real basx,basy

      integer exectime, time


      lado0 = nx * dx
      write(*,*) 'in filter'
      write(*,*) minval(dens0),maxval(dens0)
      write(*,*) minval(dens1(:,:,:,1:sum(npatch(1:nl)))), 
     &           maxval(dens1(:,:,:,1:sum(npatch(1:nl))))

*     DENS0, DENS1 PROPORTIONAL TO CELLS MASSES!!!
      dens0 = 1.0
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,DENS1,IR),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=LOW1,LOW2
         DENS1(:,:,:,I) = 1 / 8.0**IR
       END DO
      END DO

      call veinsgrid_all_l(NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &            PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,solap)

!!!!! for each cell, we ought to find the optimum coherence length
      ! we first initialize the lengths
      L0 = 0.!3.0 * DX
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA = DX / (2.0**IR)
!$OMP PARALLEL DO SHARED(LOW1,LOW2,L1,DXPA),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=LOW1,LOW2
         L1(:,:,:,I) = 0.!3.0 * DXPA
       END DO
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,L0,RADX,RADY,RADZ,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,CR0AMR1,SOLAP,
!$OMP+                   RX,RY,RZ,DENS1,U12,U13,U14,TOL,U2BULK,U3BULK,
!$OMP+                   U4BULK,STEP,MAXIT,DX,SHOCK0,SHOCK1),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    dxpa,basint,basintprev),
!$OMP+            schedule(dynamic), default(none)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            marca = 0
            if (cr0amr(i,j,k).eq.1) marca = 1
            iter = 0
            L = 3.0*dx!L0(i,j,k)
            thisx = radx(i)
            thisy = rady(j)
            thisz = radz(k)

            basintprev = 0
            iter_while_c: do while (marca.eq.1)
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
                do jj=minj,maxj
                  do ii=mini,maxi
                    if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                      bas1 = bas1 + dens0(ii,jj,kk)
                      bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                      bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                      bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                      basint = basint + 1
                      if (shock0(ii,jj,kk).eq.1) then
                        if (iter.ge.1) then
                          marca = 0
                          exit outer0_c
                        end if
                      end if
                    end if
                  end do
                end do
              end do outer0_c

              if (marca.eq.0.and.basintprev.gt.10) exit iter_while_c

              if (basint-basintprev.lt.10) then 
               l = max(l*step, l+dx)
               iter = iter + 1 
               cycle iter_while_c
              else
               basintprev=basint
              end if

              if (bas1.ne.0) then
                bas2 = bas2 / bas1
                bas3 = bas3 / bas1
                bas4 = bas4 / bas1
              else
                bas2 = 0.0
                bas3 = 0.0
                bas4 = 0.0
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
            end do iter_while_c
            if (cr0amr(i,j,k).eq.1) l0(i,j,k) = l
            !if (cr0amr(i,j,k).eq.1) write(*,*) i,j,k,iter,l,err !DEBUGGING
          end do
        end do
        !write(*,*) i, 'done'
      end do ! do i=1,nx

 !!! REFINED CELLS
      LOW1=1
      LOW2=SUM(NPATCH)
!$OMP PARALLEL DO SHARED(NX,NY,NZ,L1,radx,rady,radz,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,CR0AMR1,SOLAP,
!$OMP+                   RX,RY,RZ,DENS1,U12,U13,U14,TOL,U12BULK,U13BULK,
!$OMP+                   U14BULK,STEP,MAXIT,DX,SHOCK0,SHOCK1,LOW1,LOW2),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    ipatch,n1,n2,n3,exectime,dxpa,dxpa_i,ir,
!$OMP+                    basint,basintprev),
!$OMP+            schedule(dynamic), default(none)
      DO ipatch=LOW1,LOW2
        exectime = time()
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        ! get the level
        DO IR=1,NL
          IF (IPATCH.LE.SUM(NPATCH(0:IR))) EXIT
        END DO
        write(*,*) 'starting',ir,ipatch
        dxpa_i = dx/(2.0**ir)

c        write(*,*) ipatch, sum(cr0amr1(1:n1,1:n2,1:n3,ipatch) *
c     &                         solap(1:n1,1:n2,1:n3,ipatch))
        do k=1,n3
        do j=1,n2
        do i=1,n1
          marca = 0
          if (cr0amr1(i,j,k,ipatch).eq.1.and.
     &          solap(i,j,k,ipatch).eq.1) marca = 1
          iter = 0
          L = 3.0*dxpa_i!L1(i,j,k,ipatch)
          thisx = rx(i,ipatch)
          thisy = ry(j,ipatch)
          thisz = rz(k,ipatch)

          basintprev = 0
          iter_while: do while (marca.eq.1)
            bas1 = 0.0
            bas2 = 0.0
            bas3 = 0.0
            bas4 = 0.0
            basint = 0
            !basintprev = 0
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

            outer0: do kk=mink,maxk
              do jj=minj,maxj
                do ii=mini,maxi
                  if (cr0amr(ii,jj,kk).eq.1) then
                  if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                    bas1 = bas1 + dens0(ii,jj,kk)
                    bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                    bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                    bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                    basint = basint + 1
                    if (shock0(ii,jj,kk).eq.1) then
                      if (iter.ge.1) then
                        marca=0
                        exit outer0
                      end if
                    end if
                  end if
                  end if
                end do
              end do
            end do outer0

            if (marca.eq.0.and.basintprev.gt.10) exit iter_while

            outer1: DO irr=1,ir
             LLOW1=SUM(NPATCH(0:IRR-1))+1
             LLOW2=SUM(NPATCH(0:IRR))
             dxpa = dx / (2.0 ** irr)
             DO jpatch=LLOW1,LLOW2
               nn1 = patchnx(jpatch)
               nn2 = patchny(jpatch)
               nn3 = patchnz(jpatch)

               RX1=PATCHRX(jpatch)-0.5*dxpa
               RY1=PATCHRY(jpatch)-0.5*dxpa
               RZ1=PATCHRZ(jpatch)-0.5*dxpa
               RX2=PATCHRX(jpatch)-0.5*dxpa+(nn1-1)*dxpa
               RY2=PATCHRY(jpatch)-0.5*dxpa+(nn2-1)*dxpa
               RZ2=PATCHRZ(jpatch)-0.5*dxpa+(nn3-1)*dxpa

               RXX1 = thisx - l
               RXX2 = thisx + l
               RYY1 = thisy - l
               RYY2 = thisy + l
               RZZ1 = thisz - l
               RZZ2 = thisz + l

               IF (rxx1.le.rx2.AND.rx1.le.rxx2.AND.
     &             ryy1.le.ry2.AND.ry1.le.ryy2.AND.
     &             rzz1.le.rz2.AND.rz1.le.rzz2) then

                !X
                IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=nn1
                END IF
                IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
                   mini=1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
                   mini=1
                   maxi=nn1
                END IF

                !Y
                IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=nn2
                END IF
                IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
                   minj=1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
                   minj=1
                   maxj=nn2
                END IF

                !Z
                IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=nn3
                END IF
                IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
                   mink=1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=1
                   maxk=nn3
                END IF

                 do kkzz = mink,maxk
                 do jjyy = minj,maxj
                 do iixx = mini,maxi
                   if ((cr0amr1(iixx,jjyy,kkzz,jpatch).eq.1
     &                 .or.ir.eq.irr)
     &             .and.solap(iixx,jjyy,kkzz,jpatch).eq.1) then
                   if ((rx(iixx,jpatch)-thisx)**2+
     &                 (ry(jjyy,jpatch)-thisy)**2+
     &                 (rz(kkzz,jpatch)-thisz)**2.le.l2) then
                     bas1 = bas1 + dens1(iixx,jjyy,kkzz,jpatch)
                     bas2 = bas2 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u12(iixx,jjyy,kkzz,jpatch)
                     bas3 = bas3 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u13(iixx,jjyy,kkzz,jpatch)
                     bas4 = bas4 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u14(iixx,jjyy,kkzz,jpatch)
                     basint = basint + 1
                     if (shock1(iixx,jjyy,kkzz,jpatch).eq.1) then
                       if (iter.ge.1) then
                         marca=0
                         exit outer1
                       end if
                     end if
                   end if
                   end if
                 end do
                 end do
                 end do
              END IF
             END DO
           end do outer1

            if (marca.eq.0.and.basintprev.gt.10) exit iter_while
              
            if (basint-basintprev.lt.10) then 
             l = max(l*step, l+dxpa_i)
             iter = iter + 1
             cycle iter_while
            else
             basintprev=basint
            end if

            if (bas1.ne.0) then
              bas2 = bas2 / bas1
              bas3 = bas3 / bas1
              bas4 = bas4 / bas1
            else
              bas2 = 0.0
              bas3 = 0.0
              bas4 = 0.0
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

*             the sphere grows
            if (marca.eq.1) l = max(l*step, l+dxpa_i)
            iter = iter+1
          end do iter_while
          if (cr0amr1(i,j,k,ipatch).eq.1.and.
     &        solap(i,j,k,ipatch).eq.1) l1(i,j,k,ipatch) = l

          !DEBUGGING
C            if (cr0amr1(i,j,k,ipatch).eq.1.and.
C     &          solap(i,j,k,ipatch).eq.1) write(*,*) 'amr:',ipatch,i,j,
C     &                                                k,iter,l,err
          !END DEBUGGING
        end do
        end do
        end do ! do i=1,n1
      END DO

      write(*,*) 'refinement levels done!'


      write(*,*) '------------------------'
      write(*,*) 'pre synchro'
      do ir=0,nl
      CALL P_MINMAX_IR(L0,L1,0,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                 NPATCH,ir,BASX,BASY)
      write(*,*) 'L ir,min,max',ir,BASX,BASY
      end do
      write(*,*) '------------------------'

*     refill refined and overlapping cells
      DO IR=NL,1,-1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,L1,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U12BULK,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U13BULK,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U14BULK,nl)

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!       Apparent race condition, but harmless due to previous sync
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PARE,PATCHX,
!$OMP+                   PATCHY,PATCHZ,DENS1,U12BULK,U13BULK,U14BULK,
!$OMP+                   L1,U2BULK,U3BULK,U4BULK,L0),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,JPATCH,I,J,K,II,JJ,KK,UW,U,
!$OMP+                    FUIN),
!$OMP+            DEFAULT(NONE)
        DO ipatch=LOW1,LOW2
          !WRITE(*,*) 'FINISHING PATCH', IPATCH
          N1 = PATCHNX(IPATCH)
          N2 = PATCHNY(IPATCH)
          N3 = PATCHNZ(IPATCH)
          JPATCH = PARE(IPATCH)
          DO K=1,N3,2
          DO J=1,N2,2
          DO I=1,N1,2
            II = PATCHX(ipatch) + int((I-1)/2)
            JJ = PATCHY(ipatch) + int((J-1)/2)
            KK = PATCHZ(ipatch) + int((K-1)/2)
            if (jpatch.ne.0) then
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L1(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u12bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u13bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u14bulk(II,JJ,KK,JPATCH) = FUIN
            else
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L0(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              U2bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u3bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u4bulk(II,JJ,KK) = FUIN
            end if
          END DO
          END DO
          END DO
        END DO
      END DO !IR=NL,1,-1

      write(*,*) 'Smooth bulk!'
      call SMOOTH_BULK(U2BULK,U3BULK,U4BULK,L0,U12BULK,U13BULK,
     &                 U14BULK,L1,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                 PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                 PATCHRY,PATCHRZ,NL)

**************************************************************************
      write(*,*) 'And now use the smooth L to recompute the bulk'

!$OMP PARALLEL DO SHARED(NX,NY,NZ,L0,RADX,RADY,RADZ,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,CR0AMR1,SOLAP,
!$OMP+                   RX,RY,RZ,DENS1,U12,U13,U14,TOL,U2BULK,U3BULK,
!$OMP+                   U4BULK,STEP,MAXIT,DX,SHOCK0,SHOCK1),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    dxpa,basint,basintprev),
!$OMP+            schedule(dynamic), default(none)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            marca = 0
            if (cr0amr(i,j,k).eq.1) marca = 1
            L = L0(i,j,k)
            thisx = radx(i)
            thisy = rady(j)
            thisz = radz(k)

            if (marca.eq.1) then
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

              do kk=mink,maxk
                do jj=minj,maxj
                  do ii=mini,maxi
                    if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                      bas1 = bas1 + dens0(ii,jj,kk)
                      bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                      bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                      bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                      basint = basint + 1
                    end if
                  end do
                end do
              end do

              if (bas1.ne.0) then
                bas2 = bas2 / bas1
                bas3 = bas3 / bas1
                bas4 = bas4 / bas1
              else
                bas2 = u2bulk(i,j,k)
                bas3 = u3bulk(i,j,k)
                bas4 = u4bulk(i,j,k)
              end if

              u2bulk(i,j,k) = bas2
              u3bulk(i,j,k) = bas3
              u4bulk(i,j,k) = bas4

            end if
          end do
        end do
      end do ! do i=1,nx      

 !!! REFINED CELLS
      LOW1=1
      LOW2=SUM(NPATCH)
!$OMP PARALLEL DO SHARED(NX,NY,NZ,L1,radx,rady,radz,LADO0,CR0AMR,
!$OMP+                   DENS0,U2,U3,U4,NL,NPATCH,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,CR0AMR1,SOLAP,
!$OMP+                   RX,RY,RZ,DENS1,U12,U13,U14,TOL,U12BULK,U13BULK,
!$OMP+                   U14BULK,STEP,MAXIT,DX,SHOCK0,SHOCK1,LOW1,LOW2),
!$OMP+            PRIVATE(I,J,K,MARCA,ITER,L,THISX,THISY,THISZ,BAS1,
!$OMP+                    BAS2,BAS3,BAS4,L2,II,JJ,KK,IRR,LLOW1,LLOW2,
!$OMP+                    JPATCH,NN1,NN2,NN3,IIXX,JJYY,KKZZ,ERR,
!$OMP+                    DV2,DV3,DV4,DV2PREV,DV3PREV,DV4PREV,
!$OMP+                    MINI,MAXI,MINJ,MAXJ,MINK,MAXK,rx1,rx2,ry1,
!$OMP+                    ry2,rz1,rz2,rxx1,rxx2,ryy1,ryy2,rzz1,rzz2,
!$OMP+                    ipatch,n1,n2,n3,exectime,dxpa,dxpa_i,ir,
!$OMP+                    basint,basintprev),
!$OMP+            schedule(dynamic), default(none)
      DO ipatch=LOW1,LOW2
        exectime = time()
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        ! get the level
        DO IR=1,NL
          IF (IPATCH.LE.SUM(NPATCH(0:IR))) EXIT
        END DO
        write(*,*) 'starting2',ir,ipatch
        dxpa_i = dx/(2.0**ir)

        do k=1,n3
        do j=1,n2
        do i=1,n1
          marca = 0
          if (cr0amr1(i,j,k,ipatch).eq.1.and.
     &          solap(i,j,k,ipatch).eq.1) marca = 1
          iter = 0
          L = L1(i,j,k,ipatch)
          thisx = rx(i,ipatch)
          thisy = ry(j,ipatch)
          thisz = rz(k,ipatch)

          if (marca.eq.1) then
            bas1 = 0.0
            bas2 = 0.0
            bas3 = 0.0
            bas4 = 0.0
            basint = 0
            !basintprev = 0
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

            do kk=mink,maxk
              do jj=minj,maxj
                do ii=mini,maxi
                  if (cr0amr(ii,jj,kk).eq.1) then
                  if ((radx(ii)-thisx)**2+(rady(jj)-thisy)**2+
     &                  (radz(kk)-thisz)**2.le.l2) then
                    bas1 = bas1 + dens0(ii,jj,kk)
                    bas2 = bas2 + dens0(ii,jj,kk) * u2(ii,jj,kk)
                    bas3 = bas3 + dens0(ii,jj,kk) * u3(ii,jj,kk)
                    bas4 = bas4 + dens0(ii,jj,kk) * u4(ii,jj,kk)
                    basint = basint + 1
                  end if
                  end if
                end do
              end do
            end do

            DO irr=1,ir
             LLOW1=SUM(NPATCH(0:IRR-1))+1
             LLOW2=SUM(NPATCH(0:IRR))
             dxpa = dx / (2.0 ** irr)
             DO jpatch=LLOW1,LLOW2
               nn1 = patchnx(jpatch)
               nn2 = patchny(jpatch)
               nn3 = patchnz(jpatch)

               RX1=PATCHRX(jpatch)-0.5*dxpa
               RY1=PATCHRY(jpatch)-0.5*dxpa
               RZ1=PATCHRZ(jpatch)-0.5*dxpa
               RX2=PATCHRX(jpatch)-0.5*dxpa+(nn1-1)*dxpa
               RY2=PATCHRY(jpatch)-0.5*dxpa+(nn2-1)*dxpa
               RZ2=PATCHRZ(jpatch)-0.5*dxpa+(nn3-1)*dxpa

               RXX1 = thisx - l
               RXX2 = thisx + l
               RYY1 = thisy - l
               RYY2 = thisy + l
               RZZ1 = thisz - l
               RZZ2 = thisz + l

               IF (rxx1.le.rx2.AND.rx1.le.rxx2.AND.
     &             ryy1.le.ry2.AND.ry1.le.ryy2.AND.
     &             rzz1.le.rz2.AND.rz1.le.rzz2) then

                !X
                IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
                   mini=INT(((RXX1-RX1)/DXPA)+1) + 1
                   maxi=nn1
                END IF
                IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
                   mini=1
                   maxi=INT(((RXX2-RX1)/DXPA)) + 1
                END IF
                IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
                   mini=1
                   maxi=nn1
                END IF

                !Y
                IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
                   minj=INT(((RYY1-RY1)/dxpa)+1) + 1
                   maxj=nn2
                END IF
                IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
                   minj=1
                   maxj=INT(((RYY2-RY1)/dxpa)) + 1
                END IF
                IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
                   minj=1
                   maxj=nn2
                END IF

                !Z
                IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=INT(((RZZ1-RZ1)/dxpa)+1) + 1
                   maxk=nn3
                END IF
                IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
                   mink=1
                   maxk=INT(((RZZ2-RZ1)/dxpa)) + 1
                END IF
                IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
                   mink=1
                   maxk=nn3
                END IF

                 do kkzz = mink,maxk
                 do jjyy = minj,maxj
                 do iixx = mini,maxi
                   if ((cr0amr1(iixx,jjyy,kkzz,jpatch).eq.1
     &                 .or.ir.eq.irr)
     &             .and.solap(iixx,jjyy,kkzz,jpatch).eq.1) then
                   if ((rx(iixx,jpatch)-thisx)**2+
     &                 (ry(jjyy,jpatch)-thisy)**2+
     &                 (rz(kkzz,jpatch)-thisz)**2.le.l2) then
                     bas1 = bas1 + dens1(iixx,jjyy,kkzz,jpatch)
                     bas2 = bas2 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u12(iixx,jjyy,kkzz,jpatch)
                     bas3 = bas3 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u13(iixx,jjyy,kkzz,jpatch)
                     bas4 = bas4 + dens1(iixx,jjyy,kkzz,jpatch) *
     &                              u14(iixx,jjyy,kkzz,jpatch)
                     basint = basint + 1
                   end if
                   end if
                 end do
                 end do
                 end do
              END IF
             END DO
           end do

            if (bas1.ne.0) then
              bas2 = bas2 / bas1
              bas3 = bas3 / bas1
              bas4 = bas4 / bas1
            else
              bas2 = u12bulk(i,j,k,ipatch)
              bas3 = u13bulk(i,j,k,ipatch)
              bas4 = u14bulk(i,j,k,ipatch)
            end if

            u12bulk(i,j,k,ipatch) = bas2
            u13bulk(i,j,k,ipatch) = bas3
            u14bulk(i,j,k,ipatch) = bas4

          end if
        end do
        end do
        end do ! do i=1,n1
      END DO

*     refill refined and overlapping cells
      DO IR=NL,1,-1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,L1,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U12BULK,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U13BULK,nl)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,U14BULK,nl)

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
!       Apparent race condition, but harmless due to previous sync
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PARE,PATCHX,
!$OMP+                   PATCHY,PATCHZ,DENS1,U12BULK,U13BULK,U14BULK,
!$OMP+                   L1,U2BULK,U3BULK,U4BULK,L0),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,JPATCH,I,J,K,II,JJ,KK,UW,U,
!$OMP+                    FUIN),
!$OMP+            DEFAULT(NONE)
        DO ipatch=LOW1,LOW2
          !WRITE(*,*) 'FINISHING PATCH', IPATCH
          N1 = PATCHNX(IPATCH)
          N2 = PATCHNY(IPATCH)
          N3 = PATCHNZ(IPATCH)
          JPATCH = PARE(IPATCH)
          DO K=1,N3,2
          DO J=1,N2,2
          DO I=1,N1,2
            II = PATCHX(ipatch) + int((I-1)/2)
            JJ = PATCHY(ipatch) + int((J-1)/2)
            KK = PATCHZ(ipatch) + int((K-1)/2)
            if (jpatch.ne.0) then
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L1(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u12bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u13bulk(II,JJ,KK,JPATCH) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u14bulk(II,JJ,KK,JPATCH) = FUIN
            else
              uw(1:2,1:2,1:2) = dens1(I:I+1,J:J+1,K:K+1,IPATCH)

              u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              L0(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u12bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              U2bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u13bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u3bulk(II,JJ,KK) = FUIN

              u(1:2,1:2,1:2) = u14bulk(I:I+1,J:J+1,K:K+1,IPATCH)
              call finer_to_coarser(u,uw,fuin)
              u4bulk(II,JJ,KK) = FUIN
            end if
          END DO
          END DO
          END DO
        END DO
      END DO !IR=NL,1,-1


****************************************************************************


      write(*,*) '------------------------'
      write(*,*) 'post synchro'
      do ir=0,nl 
        CALL P_MINMAX_IR(L0,L1,0,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                 NPATCH,ir,BASX,BASY)
        write(*,*) 'L ir,min,max',ir,BASX,BASY
        CALL P_MINMAX_IR(abs(u2bulk),abs(u12bulk),0,0,NX,NY,NZ,NL,
     &                 PATCHNX,PATCHNY,PATCHNZ,NPATCH,ir,BASX,BASY)
           write(*,*) 'abs(u2bulk) ir,min,max',ir,BASX,BASY
        CALL P_MINMAX_IR(abs(u3bulk),abs(u13bulk),0,0,NX,NY,NZ,NL,
     &                 PATCHNX,PATCHNY,PATCHNZ,NPATCH,ir,BASX,BASY)
            write(*,*) 'abs(u3bulk) ir,min,max',ir,BASX,BASY
        CALL P_MINMAX_IR(abs(u4bulk),abs(u14bulk),0,0,NX,NY,NZ,NL,
     &                 PATCHNX,PATCHNY,PATCHNZ,NPATCH,ir,BASX,BASY)
            write(*,*) 'abs(u4bulk) ir,min,max',ir,BASX,BASY
      end do
      write(*,*) '------------------------'
      open(99, file='output_files/velocities_after_synchro.dat', 
     &  form='unformatted',status='unknown')

        write(99) (((l0(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(99) (((u2bulk(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(99) (((u3bulk(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        write(99) (((u4bulk(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        do i=1,sum(npatch)
         n1=patchnx(i)
         n2=patchny(i)
         n3=patchnz(i)
         write(99) (((l1(ix,j,k,i),ix=1,n1),j=1,n2),k=1,n3)
         write(99) (((u12bulk(ix,j,k,i),ix=1,n1),j=1,n2),k=1,n3)
         write(99) (((u13bulk(ix,j,k,i),ix=1,n1),j=1,n2),k=1,n3)
         write(99) (((u14bulk(ix,j,k,i),ix=1,n1),j=1,n2),k=1,n3)
        end do

      close(99)

      ! U2,U3,U4,U12,U13,U14 gets updated with the values of the
      ! velocity fluctuation
      low1=sum(npatch)
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,U12,U13,U14,
!$OMP+                   U12BULK,U13BULK,U14BULK,LOW1),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,I,J,K),
!$OMP+            DEFAULT(NONE)
      DO ipatch=1,low1
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        do k=1,n3
        do j=1,n2
        do i=1,n1
          u12(i,j,k,ipatch)=u12(i,j,k,ipatch)-u12bulk(i,j,k,ipatch)
          u13(i,j,k,ipatch)=u13(i,j,k,ipatch)-u13bulk(i,j,k,ipatch)
          u14(i,j,k,ipatch)=u14(i,j,k,ipatch)-u14bulk(i,j,k,ipatch)
        end do
        end do
        end do
      end do

!$OMP PARALLEL DO SHARED(NX,NY,NZ,U2,U3,U4,U2BULK,U3BULK,U4BULK),
!$OMP+            PRIVATE(I,J,K), DEFAULT(NONE)      
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u2(i,j,k) = u2(i,j,k) - u2bulk(i,j,k)
            u3(i,j,k) = u3(i,j,k) - u3bulk(i,j,k)
            u4(i,j,k) = u4(i,j,k) - u4bulk(i,j,k)
          end do
        end do
      end do

      CALL NOMFILE_FILTLEN(output_iter,FILENOM)
      FILERR5 = './output_files/'//FILENOM
      CALL WRITE_FILTLEN(FILERR5,NX,NY,NZ,ITER,NL,NPATCH,
     &                   PATCHNX,PATCHNY,PATCHNZ,L0,L1)

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
