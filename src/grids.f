***********************************************************************
      SUBROUTINE MALLA(NX,NY,NZ,LADO)
***********************************************************************
*     Build the coarse grid: gets the cells center coordinates and
*     interface coordinates.
***********************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,I,NY,J,NZ,K,PA4
      real A,B,C,LADO

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

*     GENERAL INITIAL CONDITIONS
*     GRID LIMITS
       A=-LADO/2.0
       B=LADO/2.0

*     GRID
*     X-AXIS
      C=(B-A)/(NX-1)
      RADX(1)=A
      DO I=2,NX
        RADX(I)=RADX(I-1)+C
      END DO
*     FICTICIOUS CELLS
      RADX(0)=RADX(1)-C
      RADX(NX+1)=RADX(NX)+C

*     Y-AXIS
      C=(B-A)/(NY-1)
      RADY(1)=A
      DO J=2,NY
        RADY(J)=RADY(J-1)+C
      END DO
*     FICTICIOUS CELLS
      RADY(0)=RADY(1)-C
      RADY(NY+1)=RADY(NY)+C

*     Z-AXIS
      C=(B-A)/(NZ-1)
      RADZ(1)=A
      DO K=2,NZ
        RADZ(K)=RADZ(K-1)+C
      END DO
*     FICTICIUS CELLS
      RADZ(0)=RADZ(1)-C
      RADZ(NZ+1)=RADZ(NZ)+C


*     COORDINATE FOR INTERFACES ***************************************
      DO I=0,NX
        RADMX(I) = (RADX(I)+RADX(I+1))/2.D0
      END DO
      DO J=0,NY
        RADMY(J) = (RADY(J)+RADY(J+1))/2.D0
      END DO
      DO K=0,NZ
        RADMZ(K) = (RADZ(K)+RADZ(K+1))/2.D0
      END DO

      DX=RADX(2)-RADX(1)
      DY=RADY(2)-RADY(1)
      DZ=RADZ(2)-RADZ(1)


      RETURN

      END

************************************************************************
      SUBROUTINE GRIDAMR(NX,NY,NZ,NL,NPATCH,
     &                   PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,
     &                   PATCHRX,PATCHRY,PATCHRZ,PARE)
************************************************************************
*     Build the AMR grid: cells center and interface positions for each
*     refinement patch. Also computes the CR3AMR variables, which
*     contain the "parent" cell of a given one which is well-inside its
*     parent patch.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,NY,NZ,NL1,NL

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      real  PATCHRX(NPALEV)
      real  PATCHRY(NPALEV)
      real  PATCHRZ(NPALEV)

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

      INTEGER I,J,K,IX,JY,KZ,IR
      INTEGER N1,N2,N3,IPA
      INTEGER INMAX(3)
      real DXPA,DYPA,DZPA

      real RX(-2:NAMRX+3,NPALEV)
      real RY(-2:NAMRX+3,NPALEV)
      real RZ(-2:NAMRX+3,NPALEV)
      real RMX(-2:NAMRX+3,NPALEV)
      real RMY(-2:NAMRX+3,NPALEV)
      real RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real xl,xr,xc,yl,yr,yc,zl,zr,zc,dxpa_i

      INTEGER NP1,NP2,NP3,ir2
      INTEGER L1,L2,L3,CR1,CR2,CR3
      INTEGER LOW1,LOW2,LOW3,LOW4
      INTEGER KR1,KR2,KR3,MARCA

*      ---PARALLEL---
      INTEGER NUM
      COMMON /PROCESADORES/ NUM

      do ir=nl,2,-1 
        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))

        dxpa_i = dx / (2.0**ir)

!$omp parallel do shared(low1,low2,patchx,patchy,patchz,patchnx,patchny,
!$omp+                   patchnz,pare,cr3amr1,cr3amr1x,cr3amr1y,
!$omp+                   cr3amr1z,dx,dy,dz,rx,ry,rz,nx,ny,nz,patchrx,
!$omp+                   patchry,patchrz,npatch,ir,dxpa_i),
!$omp+            private(i,l1,l2,l3,n1,n2,n3,ipa,np1,np2,np3,ix,jy,kz,
!$omp+                    cr1,cr2,cr3,marca,ir2,dxpa,low3,low4,xc,yc,zc,
!$omp+                    j,xl,yl,zl,xr,yr,zr),
!$omp+            default(none)
        do i=low1,low2 
          l1 = patchx(i)
          l2 = patchy(i)
          l3 = patchz(i)
          n1 = patchnx(i)
          n2 = patchny(i)
          n3 = patchnz(i)

          ipa = pare(i)
          np1 = patchnx(ipa)
          np2 = patchny(ipa)
          np3 = patchnz(ipa)

          do kz = -2, n3+3
          do jy = -2, n2+3
          do ix = -2, n1+3
            if (ix.gt.0) then
              cr1 = l1-1+int((ix+1)/2)
            else 
              cr1 = l1-1+int(ix/2)
            end if
            if (jy.gt.0) then
              cr2 = l2-1+int((jy+1)/2)
            else 
              cr2 = l2-1+int(jy/2)
            end if
            if (kz.gt.0) then
              cr3 = l3-1+int((kz+1)/2)
            else 
              cr3 = l3-1+int(kz/2)
            end if

            if (cr1.ge.1.and.cr1.le.np1.and.
     &          cr2.ge.1.and.cr2.le.np2.and.
     &          cr3.ge.1.and.cr3.le.np3) then 
              cr3amr1(ix,jy,kz,i) = ipa
              cr3amr1x(ix,jy,kz,i) = cr1
              cr3amr1y(ix,jy,kz,i) = cr2
              cr3amr1z(ix,jy,kz,i) = cr3
            else ! search among uncles, or below
              marca = 0
              xc = patchrx(i) + (ix-1.5)*dxpa_i!rx(ix,i)
              yc = patchry(i) + (jy-1.5)*dxpa_i!ry(jy,i)
              zc = patchrz(i) + (kz-1.5)*dxpa_i!rz(kz,i)
              do ir2 = ir-1,1,-1
                dxpa = dx / (2.0**ir2)
                low3 = sum(npatch(0:ir2-1))+1
                low4 = sum(npatch(0:ir2)) 
                do j = low3, low4 
                  np1 = patchnx(j)
                  np2 = patchny(j)
                  np3 = patchnz(j)
                  xl = patchrx(j) - dxpa
                  xr = xl + np1*dxpa
                  if (xc.lt.xl.or.xc.gt.xr) cycle
                  yl = patchry(j) - dxpa
                  yr = yl + np2*dxpa
                  if (yc.lt.yl.or.yc.gt.yr) cycle
                  zl = patchrz(j) - dxpa
                  zr = zl + np3*dxpa
                  if (zc.lt.zl.or.zc.gt.zr) cycle
                  ! If it has not cycled, it is inside!
                  marca = 1
                  cr1 = int((xc-xl)/dxpa) + 1
                  cr2 = int((yc-yl)/dxpa) + 1
                  cr3 = int((zc-zl)/dxpa) + 1
                  cr3amr1(ix,jy,kz,i) = j
                  cr3amr1x(ix,jy,kz,i) = cr1
                  cr3amr1y(ix,jy,kz,i) = cr2
                  cr3amr1z(ix,jy,kz,i) = cr3
                  exit
                end do 
                if (marca.eq.1) exit
              end do
              if (marca.eq.0) then 
                xl = -(nx/2)*dx
                yl = -(ny/2)*dy
                zl = -(nz/2)*dz
                cr1 = int((xc-xl)/dx) + 1
                cr2 = int((yc-yl)/dy) + 1
                cr3 = int((zc-zl)/dz) + 1
                cr3amr1(ix,jy,kz,i) = 0
                cr3amr1x(ix,jy,kz,i) = cr1
                cr3amr1y(ix,jy,kz,i) = cr2
                cr3amr1z(ix,jy,kz,i) = cr3
              end if
            end if
          end do 
          end do 
          end do
        end do

      end do


      IR=1
!$OMP PARALLEL DO SHARED(IR,NL,PATCHX,PATCHY,PATCHZ,PATCHNX,PATCHNY,
!$OMP+                   PATCHNZ,CR3AMR1,NX,NY,NZ,CR3AMR1X,CR3AMR1Y,
!$OMP+                   CR3AMR1Z,PARE,NPATCH),
!$OMP+            PRIVATE(I,L1,L2,L3,IX,JY,KZ,KR1,KR2,KR3,CR1,CR2,CR3),
!$OMP+            DEFAULT(NONE)
      DO I=1,NPATCH(IR)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       DO KZ=-2,PATCHNZ(I)+3
       DO JY=-2,PATCHNY(I)+3
       DO IX=-2,PATCHNX(I)+3

         CR1=L1-1+INT((IX+1)/2)
         CR2=L2-1+INT((JY+1)/2)
         CR3=L3-1+INT((KZ+1)/2)

         IF (IX.LT.-1) CR1=INT((IX+1)/2)+L1-2
         IF (JY.LT.-1) CR2=INT((JY+1)/2)+L2-2
         IF (KZ.LT.-1) CR3=INT((KZ+1)/2)+L3-2

         KR1=CR1
         KR2=CR2
         KR3=CR3

         CR3AMR1(IX,JY,KZ,I)=0
         CR3AMR1X(IX,JY,KZ,I)=KR1
         CR3AMR1Y(IX,JY,KZ,I)=KR2
         CR3AMR1Z(IX,JY,KZ,I)=KR3

       END DO
       END DO
       END DO
      END DO

*     MINI GRIDS
      RX=0.0
      RY=0.0
      RZ=0.0
      RMX=0.0
      RMY=0.0
      RMZ=0.0
      DO IR=1,NL
       DXPA=DX/(2.0**IR)
       DYPA=DY/(2.0**IR)
       DZPA=DZ/(2.0**IR)
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)

        CALL MINIMALLA(N1,N2,N3,DXPA,DYPA,DZPA,PATCHRX(I),
     &        PATCHRY(I),PATCHRZ(I),
     &        RX(-2:N1+3,I),RY(-2:N2+3,I),RZ(-2:N3+3,I),
     &        RMX(-2:N1+3,I),RMY(-2:N2+3,I),RMZ(-2:N3+3,I))

       END DO
      END DO


      RETURN
      END

**********************************************************************
      SUBROUTINE MINIMALLA(N1,N2,N3,DX,DY,DZ,RPAX,RPAY,RPAZ,
     &                     RX,RY,RZ,RMX,RMY,RMZ)
**********************************************************************
*     Build the grid for each patch
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER N1,N2,N3,I
      real RX(-2:N1+3),RY(-2:N2+3),RZ(-2:N3+3)
      real RMX(-2:N1+3),RMY(-2:N2+3),RMZ(-2:N3+3)
      real RPAX,RPAY,RPAZ,DX,DY,DZ

*     X
      DO I=1,N1
        RX(I)=RPAX-(DX/2.0)+(I-1)*DX
      END DO

      RX(0)=RX(1)-DX
      RX(N1+1)=RX(N1)+DX
      RX(-1)=RX(0)-DX
      RX(N1+2)=RX(N1+1)+DX
      RX(-2)=RX(-1)-DX
      RX(N1+3)=RX(N1+2)+DX

      DO I=-2,N1+2
        RMX(I)=(RX(I)+RX(I+1))/2.0
      END DO

*     Y
      DO I=1,N2
        RY(I)=RPAY-(DY/2.0)+(I-1)*DY
      END DO

      RY(0)=RY(1)-DY
      RY(N2+1)=RY(N2)+DY
      RY(-1)=RY(0)-DY
      RY(N2+2)=RY(N2+1)+DY
      RY(-2)=RY(-1)-DY
      RY(N2+3)=RY(N2+2)+DY

      DO I=-2,N2+2
        RMY(I)=(RY(I)+RY(I+1))/2.0
      END DO

*     Z
      DO I=1,N3
        RZ(I)=RPAZ-(DZ/2.0)+(I-1)*DZ
      END DO

      RZ(0)=RZ(1)-DZ
      RZ(N3+1)=RZ(N3)+DZ
      RZ(-1)=RZ(0)-DZ
      RZ(N3+2)=RZ(N3+1)+DZ
      RZ(-2)=RZ(-1)-DZ
      RZ(N3+3)=RZ(N3+2)+DZ

      DO I=-2,N3+2
        RMZ(I)=(RZ(I)+RZ(I+1))/2.0
      END DO


      RETURN
      END

***********************************************************************
       subroutine extend_var(nx,ny,nz,nl,npatch,
     &            pare,patchnx,patchny,patchnz,patchx,patchy,patchz,
     &            patchrx,patchry,patchrz)
***********************************************************************
*     extend variables one cell on each side by interpolation from
*     coarser patches.
***********************************************************************

       implicit none

       include 'vortex_parameters.dat'

       integer nx,ny,nz,i,j,k,ir,ix,jy,kz
       integer nl,n1,n2,n3,l1,l2,l3
       integer ii,jj,kk,low1, low2

       real  radx(0:nmax+1),radmx(0:nmax+1),
     &         rady(0:nmay+1),radmy(0:nmay+1),
     &         radz(0:nmaz+1),radmz(0:nmaz+1)
       common /grid/   radx,radmx,rady,radmy,radz,radmz

       real u2(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u3(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u4(0:nmax+1,0:nmay+1,0:nmaz+1)
       real u12(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u13(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       real u14(0:namrx+1,0:namry+1,0:namrz+1,npalev)
       common /veloc/ u2,u3,u4,u12,u13,u14

       real rotax_0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real rotay_0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real rotaz_0(0:nmax+1,0:nmay+1,0:nmaz+1)
       real rotax_1(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       real rotay_1(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       real rotaz_1(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       common /rots/ rotax_0,rotay_0,rotaz_0,rotax_1,rotay_1,rotaz_1

       real dx,dy,dz
       common /espaciado/ dx,dy,dz

       integer npatch(0:nlevels),pare(npalev)
       integer patchnx(npalev),patchny(npalev),patchnz(npalev)
       integer patchx(npalev),patchy(npalev),patchz(npalev)
       real patchrx(npalev),patchry(npalev),patchrz(npalev)

       real dxpa,dypa,dzpa,xxx1,yyy1,zzz1
       integer cr1,cr2,cr3
       integer mark, abuelo, kr1, kr2, kr3, kare,ir_abue
       real xc,yc,zc,xl,yl,zl,xr,yr,zr 
       integer marca,ip,nn1,nn2,nn3

       real ubas(1:3,1:3,1:3),fuin,rxbas(3),rybas(3),rzbas(3)
       real aaa,bbb,ccc

       integer cr3amr1(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       integer cr3amr1x(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       integer cr3amr1y(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       integer cr3amr1z(-2:namrx+3,-2:namry+3,-2:namrz+3,npalev)
       common /cr0cell/ cr3amr1,cr3amr1x,cr3amr1y,cr3amr1z

       real rx(-2:namrx+3,npalev)
       real ry(-2:namrx+3,npalev)
       real rz(-2:namrx+3,npalev)
       real rmx(-2:namrx+3,npalev)
       real rmy(-2:namrx+3,npalev)
       real rmz(-2:namrx+3,npalev)
       common /minigrids/ rx,ry,rz,rmx,rmy,rmz


*      ---parallel---
       integer num,omp_get_num_threads,numor, flag_parallel
       common /procesadores/ num

c     the code below is useless for the particle version.
c     it is necessary, however, for the grid-based input.

#if input_is_grid == 1
*      expansion of the patches by interpolation from coarser level or, 
*       if possible, by copying the values from a sibling patch (not yet implemented)

      ir=1
      dxpa=dx/(2.0**ir)
      dypa=dy/(2.0**ir)
      dzpa=dz/(2.0**ir)
      low1=1
      low2=npatch(ir)

!$omp parallel do shared(npatch,ir,patchnx,patchny,patchnz,cr3amr1x,
!$omp+                   cr3amr1y,cr3amr1z,u2,u3,u4,u12,u13,u14,
!$omp+                   patchrx,patchry,patchrz,dxpa,dypa,dzpa,rx,ry,
!$omp+                   rz,low1,low2),
!$omp+            private(i,n1,n2,n3,ix,jy,kz,kr1,kr2,kr3,ubas,fuin,
!$omp+                    xl,yl,zl,xr,yr,zr,xc,yc,zc,marca,ip,nn1,nn2,
!$omp+                    nn3),
!$omp+            default(none)
      do i=low1,low2

       n1=patchnx(i)
       n2=patchny(i)
       n3=patchnz(i)

       do kz=0,n3+1
       do jy=0,n2+1
       do ix=0,n1+1
        ! We first look if we can copy the value from a sibling patch
        marca=0
        xc = rx(ix,i)
        yc = ry(jy,i)
        zc = rz(kz,i)
        do ip=low1,low2
          nn1 = patchnx(ip)
          xl = patchrx(ip) - dxpa 
          xr = xl + nn1*dxpa
          if (xc.lt.xl.or.xc.gt.xr) cycle
          nn2 = patchny(ip)
          yl = patchry(ip) - dxpa
          yr = yl + nn2*dxpa
          if (yc.lt.yl.or.yc.gt.yr) cycle
          nn3 = patchnz(ip)
          zl = patchrz(ip) - dxpa
          zr = zl + nn3*dxpa
          if (zc.lt.zl.or.zc.gt.zr) cycle
          ! else it is here 
          
          kr1 = int((xc-xl)/dxpa) + 1
          kr2 = int((yc-yl)/dxpa) + 1
          kr3 = int((zc-zl)/dxpa) + 1    
          if (kr1.lt.1.or.kr1.gt.nn1) cycle
          if (kr2.lt.1.or.kr2.gt.nn2) cycle
          if (kr3.lt.1.or.kr3.gt.nn3) cycle    
          u12(ix,jy,kz,i) = u12(kr1,kr2,kr3,ip)
          u13(ix,jy,kz,i) = u13(kr1,kr2,kr3,ip)
          u14(ix,jy,kz,i) = u14(kr1,kr2,kr3,ip)
          
          marca=1
          exit
        end do ! ip   
        if (marca.ne.0) cycle
        ! We have not found a sibling patch, so we have to interpolate

*      #######################################################
       if (ix.lt.1.or.ix.gt.n1.or.jy.lt.1.or.jy.gt.n2.or.
     &     kz.lt.1.or.kz.gt.n3) then
*      #######################################################

        kr1=cr3amr1x(ix,jy,kz,i)
        kr2=cr3amr1y(ix,jy,kz,i)
        kr3=cr3amr1z(ix,jy,kz,i)

        !vx
        ubas(1:3,1:3,1:3)=u2(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
        call linint52d_new(ix,jy,kz,ubas,fuin)
        u12(ix,jy,kz,i)=fuin

        !vy
        ubas(1:3,1:3,1:3)=u3(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
        call linint52d_new(ix,jy,kz,ubas,fuin)
        u13(ix,jy,kz,i)=fuin

        !vz
        ubas(1:3,1:3,1:3)=u4(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
        call linint52d_new(ix,jy,kz,ubas,fuin)
        u14(ix,jy,kz,i)=fuin

*      #######################################################
        end if
*      #######################################################

        end do
        end do
        end do
       end do


       do ir=2,nl
        dxpa=dx/(2.0**ir)
        dypa=dy/(2.0**ir)
        dzpa=dz/(2.0**ir)

        low1=sum(npatch(0:ir-1))+1
        low2=sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,cr3amr1,
!$omp+                   cr3amr1x,cr3amr1y,cr3amr1z,rx,ry,rz,u2,u3,u4,
!$omp+                   u12,u13,u14,radx,rady,radz,patchrx,patchry,
!$omp+                   patchrz,dxpa,dypa,dzpa),
!$omp+            private(i,n1,n2,n3,ix,jy,kz,kare,kr1,kr2,kr3,aaa,bbb,
!$omp+                    ccc,rxbas,rybas,rzbas,ubas,fuin,xl,xr,yl,yr,
!$omp+                    zl,zr,xc,yc,zc,nn1,nn2,nn3,marca,ip),
!$omp+            default(none)
        do i=low1,low2

        n1=patchnx(i)
        n2=patchny(i)
        n3=patchnz(i)

         do kz=0,n3+1
         do jy=0,n2+1
         do ix=0,n1+1

         if (ix.lt.1.or.ix.gt.n1.or.
     &       jy.lt.1.or.jy.gt.n2.or.
     &       kz.lt.1.or.kz.gt.n3) then
             ! We first look if we can copy the value from a sibling patch
             marca=0
             xc = rx(ix,i)
             yc = ry(jy,i)
             zc = rz(kz,i)
             do ip=low1,low2
               nn1 = patchnx(ip)
               xl = patchrx(ip) - dxpa 
               xr = xl + nn1*dxpa
               if (xc.lt.xl.or.xc.gt.xr) cycle
               nn2 = patchny(ip)
               yl = patchry(ip) - dxpa
               yr = yl + nn2*dxpa
               if (yc.lt.yl.or.yc.gt.yr) cycle
               nn3 = patchnz(ip)
               zl = patchrz(ip) - dxpa
               zr = zl + nn3*dxpa
               if (zc.lt.zl.or.zc.gt.zr) cycle
               ! else it is here 
               
               kr1 = int((xc-xl)/dxpa) + 1
               kr2 = int((yc-yl)/dxpa) + 1
               kr3 = int((zc-zl)/dxpa) + 1    
               if (kr1.lt.1.or.kr1.gt.nn1) cycle
               if (kr2.lt.1.or.kr2.gt.nn2) cycle
               if (kr3.lt.1.or.kr3.gt.nn3) cycle    
               u12(ix,jy,kz,i) = u12(kr1,kr2,kr3,ip)
               u13(ix,jy,kz,i) = u13(kr1,kr2,kr3,ip)
               u14(ix,jy,kz,i) = u14(kr1,kr2,kr3,ip)
               
               marca=1
               exit
             end do ! ip   
             if (marca.ne.0) cycle
             ! We have not found a sibling patch, so we have to interpolate

             kare=cr3amr1(ix,jy,kz,i)
             kr1=cr3amr1x(ix,jy,kz,i)
             kr2=cr3amr1y(ix,jy,kz,i)
             kr3=cr3amr1z(ix,jy,kz,i)

             aaa=rx(ix,i)
             bbb=ry(jy,i)
             ccc=rz(kz,i)

             if (kare.gt.0) then

             rxbas(1:3)=rx(kr1-1:kr1+1,kare)
             rybas(1:3)=ry(kr2-1:kr2+1,kare)
             rzbas(1:3)=rz(kr3-1:kr3+1,kare)

             !vx:
              ubas(1:3,1:3,1:3)=
     &             u12(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1,kare)
              call linint52d_new_real(aaa,bbb,ccc,
     &                             rxbas,rybas,rzbas,ubas,fuin)
              u12(ix,jy,kz,i)=fuin

              !vy:
              ubas(1:3,1:3,1:3)=
     &             u13(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1,kare)
              call linint52d_new_real(aaa,bbb,ccc,
     &                                rxbas,rybas,rzbas,ubas,fuin)
              u13(ix,jy,kz,i)=fuin

              !vz:
              ubas(1:3,1:3,1:3)=
     &             u14(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1,kare)
              call linint52d_new_real(aaa,bbb,ccc,
     &                                rxbas,rybas,rzbas,ubas,fuin)
              u14(ix,jy,kz,i)=fuin


             else
              rxbas(1:3)=radx(kr1-1:kr1+1)
              rybas(1:3)=rady(kr2-1:kr2+1)
              rzbas(1:3)=radz(kr3-1:kr3+1)

              !vx
              ubas(1:3,1:3,1:3)=
     &             u2(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
              call linint52d_new_real(aaa,bbb,ccc,
     &                                rxbas,rybas,rzbas,ubas,fuin)
              u12(ix,jy,kz,i)=fuin

              !vy
              ubas(1:3,1:3,1:3)=
     &             u3(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
              call linint52d_new_real(aaa,bbb,ccc,
     &                                rxbas,rybas,rzbas,ubas,fuin)
              u13(ix,jy,kz,i)=fuin

              !vz
              ubas(1:3,1:3,1:3)=
     &            u4(kr1-1:kr1+1,kr2-1:kr2+1,kr3-1:kr3+1)
              call linint52d_new_real(aaa,bbb,ccc,
     &                                rxbas,rybas,rzbas,ubas,fuin)
              u14(ix,jy,kz,i)=fuin

             endif
          endif
         end do
         end do
         end do

        end do
       end do
#endif 

* ir=0 (periodic boundary)
       do k=0,nz+1
       do j=0,ny+1
       do i=0,nx+1
        if (i.lt.1.or.i.gt.nx.or.
     &      j.lt.1.or.j.gt.ny.or.
     &      k.lt.1.or.k.gt.nz) then

        ix=i
        jy=j
        kz=k
        if (i.lt.1) ix=i+nx
        if (j.lt.1) jy=j+ny
        if (k.lt.1) kz=k+nz
        if (i.gt.nx) ix=i-nx
        if (j.gt.ny) jy=j-ny
        if (k.gt.nz) kz=k-nz

        u2(i,j,k)=u2(ix,jy,kz)
        u3(i,j,k)=u3(ix,jy,kz)
        u4(i,j,k)=u4(ix,jy,kz)

        end if
       end do
       end do
       end do

***    all variables have been extended one cell

       return
       end

#if input_is_grid == 1
***********************************************************************
      subroutine compute_cr0amr(nx,ny,nz,nl,npatch,patchnx,patchny,
     &                          patchnz,patchx,patchy,patchz,patchrx,
     &                          patchry,patchrz,pare)
*********************************************************************** 
*     This subroutine just computes an integer variable over the grid,
*      which is 1 for unrefined cells, 0 for refined cells.
***********************************************************************
      implicit none
      include 'vortex_parameters.dat'
      
      integer nx,ny,nz,nl 
      integer npatch(0:nlevels)
      integer patchnx(npalev), patchny(npalev), patchnz(npalev)
      integer patchx(npalev), patchy(npalev), patchz(npalev)
      real patchrx(npalev), patchry(npalev), patchrz(npalev)
      integer pare(npalev)


      integer cr0amr(1:nmax,1:nmay,1:nmaz)
      integer cr0amr1(1:namrx,1:namry,1:namrz,npalev)
      common /cr0/ cr0amr, cr0amr1

      real dx,dy,dz
      common /espaciado/ dx,dy,dz

      integer low1,low2,ir,irpa,i1,j1,k1,i2,j2,k2,n1,n2,n3,ip,ipa,i,j,k
      integer low1pa,low2pa

*     Base level 
!$omp parallel do shared(nx,ny,nz,cr0amr), private(i,j,k), default(none)
      do k=1,nz 
      do j=1,ny
      do i=1,nx 
        cr0amr(i,j,k) = 1
      end do
      end do
      end do


      irpa = 0 
      ir = 1 
      low1 = sum(npatch(0:ir-1))+1
      low2 = sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,patchx,patchy,patchz,patchnx,
!$omp+                   patchny,patchnz,cr0amr),
!$omp+            private(ip,i1,j1,k1,i2,j2,k2,n1,n2,n3),
!$omp+            default(none)
      do ip = low1, low2 
        i1 = patchx(ip)
        j1 = patchy(ip)
        k1 = patchz(ip)

        n1 = patchnx(ip)
        n2 = patchny(ip)
        n3 = patchnz(ip)

        i2 = i1 + n1/2 - 1
        j2 = j1 + n2/2 - 1
        k2 = k1 + n3/2 - 1

        do k = k1, k2 
        do j = j1, j2
        do i = i1, i2 
          cr0amr(i,j,k) = 0
        end do
        end do
        end do

      end do

*     Refinement levels
      do ir = 2, nl 
        irpa = ir - 1

        low1pa = sum(npatch(0:irpa-1))+1
        low2pa = sum(npatch(0:irpa))

!$omp parallel do shared(low1pa,low2pa,cr0amr1), private(ipa),
!$omp+            default(none)
        do ipa = low1pa, low2pa 
          cr0amr1(:,:,:,ipa) = 1
        end do

        low1 = sum(npatch(0:ir-1))+1
        low2 = sum(npatch(0:ir))

!$omp parallel do shared(low1,low2,patchx,patchy,patchz,patchnx,
!$omp+                   patchny,patchnz,pare,cr0amr1), 
!$omp+            private(ip,i1,j1,k1,n1,n2,n3,i2,j2,k2,ipa,i,j,k),
!$omp+            default(none)
        do ip = low1, low2 
          i1 = patchx(ip)
          j1 = patchy(ip)
          k1 = patchz(ip)

          n1 = patchnx(ip)
          n2 = patchny(ip)
          n3 = patchnz(ip)

          i2 = i1 + n1/2 - 1
          j2 = j1 + n2/2 - 1
          k2 = k1 + n3/2 - 1

          ipa = pare(ip)

          do k = k1, k2 
          do j = j1, j2
          do i = i1, i2 
            cr0amr1(i,j,k,ipa) = 0
          end do
          end do
          end do

        end do

        call veinsgrid_cr0amr(irpa,npatch,pare,patchnx,patchny,patchnz,
     &                       patchx,patchy,patchz,patchrx,patchry,
     &                       patchrz)

      end do ! ir = 2, nl

*     The last level 
      irpa = nl 
      low1pa = sum(npatch(0:irpa-1))+1
      low2pa = sum(npatch(0:irpa))
!$omp parallel do shared(low1pa,low2pa,cr0amr1), private(ipa),
!$omp+            default(none)
      do ipa = low1pa, low2pa 
        cr0amr1(:,:,:,ipa) = 1
      end do

      return
      end

#endif 