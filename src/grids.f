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
            cr1 = l1-1+int((ix+1)/2)
            cr2 = l2-1+int((jy+1)/2)
            cr3 = l3-1+int((kz+1)/2)

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
       SUBROUTINE EXTEND_VAR(NX,NY,NZ,NL,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ)
***********************************************************************
*     Extend variables one cell on each side by interpolation from
*     coarser patches.
***********************************************************************

       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ,I,J,K,IR,IX,JY,KZ
       INTEGER NL,N1,N2,N3,L1,L2,L3
       INTEGER II,JJ,KK,LOW1, LOW2

       real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &         RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &         RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
       COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

       real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
       COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

       real DX,DY,DZ
       COMMON /ESPACIADO/ DX,DY,DZ

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       real DXPA,DYPA,DZPA,XXX1,YYY1,ZZZ1
       INTEGER CR1,CR2,CR3
       INTEGER MARK, ABUELO, KR1, KR2, KR3, KARE,IR_ABUE

       real UBAS(1:3,1:3,1:3),FUIN,RXBAS(3),RYBAS(3),RZBAS(3)
       real AAA,BBB,CCC

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


*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

c     The code below is useless for the particle version.
c     It is necessary, however, for the grid-based input.
c     I keep it here for the future unification of the code.

cc *      EXPANSION OF THE PATCHES BY INTERPOLATION
cc 
cc       IR=1
cc       DXPA=DX/(2.0**IR)
cc       DYPA=DY/(2.0**IR)
cc       DZPA=DZ/(2.0**IR)
cc 
cc !$OMP PARALLEL DO SHARED(NPATCH,IR,PATCHNX,PATCHNY,PATCHNZ,CR3AMR1X,
cc !$OMP+                   CR3AMR1Y,CR3AMR1Z,U2,U3,U4,U12,U13,U14),
cc !$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ,KR1,KR2,KR3,UBAS,FUIN),
cc !$OMP+            DEFAULT(NONE)
cc       DO I=1,NPATCH(IR)
cc 
cc        N1=PATCHNX(I)
cc        N2=PATCHNY(I)
cc        N3=PATCHNZ(I)
cc 
cc        DO KZ=0,N3+1
cc        DO JY=0,N2+1
cc        DO IX=0,N1+1
cc 
cc *      #######################################################
cc        IF (IX.LT.1.OR.IX.GT.N1.OR.JY.LT.1.OR.JY.GT.N2.OR.
cc      &     KZ.LT.1.OR.KZ.GT.N3) THEN
cc *      #######################################################
cc 
cc         KR1=CR3AMR1X(IX,JY,KZ,I)
cc         KR2=CR3AMR1Y(IX,JY,KZ,I)
cc         KR3=CR3AMR1Z(IX,JY,KZ,I)
cc 
cc         !VX
cc         UBAS(1:3,1:3,1:3)=U2(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc         CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
cc         U12(IX,JY,KZ,I)=FUIN
cc 
cc         !VY
cc         UBAS(1:3,1:3,1:3)=U3(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc         CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
cc         U13(IX,JY,KZ,I)=FUIN
cc 
cc         !VZ
cc         UBAS(1:3,1:3,1:3)=U4(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc         CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
cc         U14(IX,JY,KZ,I)=FUIN
cc 
cc *      #######################################################
cc         END IF
cc *      #######################################################
cc 
cc         END DO
cc         END DO
cc         END DO
cc        END DO
cc 
cc 
cc        DO IR=2,NL
cc         DXPA=DX/(2.0**IR)
cc         DYPA=DY/(2.0**IR)
cc         DZPA=DZ/(2.0**IR)
cc 
cc         LOW1=SUM(NPATCH(0:IR-1))+1
cc         LOW2=SUM(NPATCH(0:IR))
cc !$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR3AMR1,
cc !$OMP+                   CR3AMR1X,CR3AMR1Y,CR3AMR1Z,RX,RY,RZ,U2,U3,U4,
cc !$OMP+                   U12,U13,U14,RADX,RADY,RADZ),
cc !$OMP+            PRIVATE(I,N1,N2,N3,IX,JY,KZ,KARE,KR1,KR2,KR3,AAA,BBB,
cc !$OMP+                    CCC,RXBAS,RYBAS,RZBAS,UBAS,FUIN),
cc !$OMP+            DEFAULT(NONE)
cc         DO I=LOW1,LOW2
cc 
cc         N1=PATCHNX(I)
cc         N2=PATCHNY(I)
cc         N3=PATCHNZ(I)
cc 
cc          DO KZ=0,N3+1
cc          DO JY=0,N2+1
cc          DO IX=0,N1+1
cc 
cc          IF (IX.LT.1.OR.IX.GT.N1.OR.
cc      &       JY.LT.1.OR.JY.GT.N2.OR.
cc      &       KZ.LT.1.OR.KZ.GT.N3) THEN
cc 
cc              KARE=CR3AMR1(IX,JY,KZ,I)
cc              KR1=CR3AMR1X(IX,JY,KZ,I)
cc              KR2=CR3AMR1Y(IX,JY,KZ,I)
cc              KR3=CR3AMR1Z(IX,JY,KZ,I)
cc 
cc              AAA=RX(IX,I)
cc              BBB=RY(JY,I)
cc              CCC=RZ(KZ,I)
cc 
cc              IF (KARE.GT.0) THEN
cc 
cc              RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
cc              RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
cc              RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)
cc 
cc              !Vx:
cc               UBAS(1:3,1:3,1:3)=
cc      &             U12(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                             RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U12(IX,JY,KZ,I)=FUIN
cc 
cc               !Vy:
cc               UBAS(1:3,1:3,1:3)=
cc      &             U13(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U13(IX,JY,KZ,I)=FUIN
cc 
cc               !Vz:
cc               UBAS(1:3,1:3,1:3)=
cc      &             U14(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U14(IX,JY,KZ,I)=FUIN
cc 
cc 
cc              ELSE
cc               RXBAS(1:3)=RADX(KR1-1:KR1+1)
cc               RYBAS(1:3)=RADY(KR2-1:KR2+1)
cc               RZBAS(1:3)=RADZ(KR3-1:KR3+1)
cc 
cc               !Vx
cc               UBAS(1:3,1:3,1:3)=
cc      &             U2(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U12(IX,JY,KZ,I)=FUIN
cc 
cc               !Vy
cc               UBAS(1:3,1:3,1:3)=
cc      &             U3(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U13(IX,JY,KZ,I)=FUIN
cc 
cc               !Vz
cc               UBAS(1:3,1:3,1:3)=
cc      &            U4(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
cc               CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
cc      &                                RXBAS,RYBAS,RZBAS,UBAS,FUIN)
cc               U14(IX,JY,KZ,I)=FUIN
cc 
cc              ENDIF
cc           ENDIF
cc          END DO
cc          END DO
cc          END DO
cc 
cc         END DO
cc        END DO

* IR=0 (periodic boundary)
       DO K=0,NZ+1
       DO J=0,NY+1
       DO I=0,NX+1
        IF (I.LT.1.OR.I.GT.NX.OR.
     &      J.LT.1.OR.J.GT.NY.OR.
     &      K.LT.1.OR.K.GT.NZ) THEN

        IX=I
        JY=J
        KZ=K
        IF (I.LT.1) IX=I+NX
        IF (J.LT.1) JY=J+NY
        IF (K.LT.1) KZ=K+NZ
        IF (I.GT.NX) IX=I-NX
        IF (J.GT.NY) JY=J-NY
        IF (K.GT.NZ) KZ=K-NZ

        U2(I,J,K)=U2(IX,JY,KZ)
        U3(I,J,K)=U3(IX,JY,KZ)
        U4(I,J,K)=U4(IX,JY,KZ)

        END IF
       END DO
       END DO
       END DO

***    ALL VARIABLES HAVE BEEN EXTENDED ONE CELL

       RETURN
       END
