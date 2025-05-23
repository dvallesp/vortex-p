************************************************************************
      SUBROUTINE POTBASE(NX,NY,NZ,KKK)
************************************************************************
      IMPLICIT NONE 

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,NY,NZ 
      REAL KKK(NMAX,NMAY,NMAZ)

#if use_fftw==0
      CALL POFFT3D(NX,NY,NZ,KKK)
#else 
      CALL POFFT3D_FFTW(NX,NY,NZ,KKK)
#endif

      END

************************************************************************
      SUBROUTINE POFFT3D(NX,NY,NZ,KKK)
************************************************************************
*     Solve Poisson's equation in Fourier space, assuming periodic
*     boundary conditions.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,NY,NZ,I,J,K

      real  U1(0:NMAX+1,0:NMAY+1,0:NMAZ+1)    ! source in poisson equation
      real  POT(0:NMAX+1,0:NMAY+1,0:NMAZ+1)   ! field to solve
      COMMON /BASE/ U1,POT

*     FFT variables
      real KKK(NMAX,NMAY,NMAZ)
      real DATA1(2*NMAX*NMAY*NMAZ)

      INTEGER I1,I2,IJK,NYZ,NXYZ,NNN(3)
      INTEGER NZ2,NZ3,NX2,NX3,NY2,NY3


      NYZ=NY*NZ
      NXYZ=NX*NY*NZ

**    FFT METHOD
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
        IJK=(I-1)*NYZ + (J-1)*NZ + K
        I1= 2*IJK -1
        I2= 2*IJK
        DATA1(I1)=U1(I,J,K)    ! real part
        DATA1(I2)=0.0          ! imaginary part
      END DO
      END DO
      END DO

      NNN(1)=NX
      NNN(2)=NY
      NNN(3)=NZ

      CALL FOURN(DATA1,NNN,3,1)

**    POISSON EQUATION IN MOMENTUM SPACE
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
         IJK=(I-1)*NYZ + (J-1)*NZ + K
         I1= 2*IJK -1
         I2= 2*IJK
         DATA1(I1)=DATA1(I1)*KKK(I,J,K)
         DATA1(I2)=DATA1(I2)*KKK(I,J,K)
      END DO
      END DO
      END DO
      CALL FOURN(DATA1,NNN,3,-1)

      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
         IJK=(I-1)*NYZ + (J-1)*NZ + K
         I1= 2*IJK -1
         POT(I,J,K)=DATA1(I1)/NXYZ
      END DO
      END DO
      END DO


*     FICTICIOUS CELLS FOR FIELD COMPUTATION
      DO K=1,NZ
      DO J=1,NY
         POT(0,J,K)=POT(NX,J,K)
         POT(NX+1,J,K)=POT(1,J,K)
      END DO
      END DO
      DO K=1,NZ
      DO I=0,NX+1
         POT(I,0,K)=POT(I,NY,K)
         POT(I,NY+1,K)=POT(I,1,K)
      END DO
      END DO
      DO J=0,NY+1
      DO I=0,NX+1
         POT(I,J,0)=POT(I,J,NZ)
         POT(I,J,NZ+1)=POT(I,J,1)
      END DO
      END DO


      RETURN
      END

#if use_fftw==1
************************************************************************
      SUBROUTINE POFFT3D_FFTW(NX,NY,NZ,KKK)
************************************************************************
*     Solve Poisson's equation in Fourier space, assuming periodic
*     boundary conditions.
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NX,NY,NZ,I,J,K

      real  U1(0:NMAX+1,0:NMAY+1,0:NMAZ+1)    ! source in poisson equation
      real  POT(0:NMAX+1,0:NMAY+1,0:NMAZ+1)   ! field to solve
      COMMON /BASE/ U1,POT

*     FFT variables
      real KKK(NMAX,NMAY,NMAZ)
      COMPLEX DATA1(NMAX,NMAY,NMAZ)

      INTEGER I1,I2,IJK,NYZ,NXYZ,NNN(3)
      INTEGER NZ2,NZ3,NX2,NX3,NY2,NY3

      INTEGER NUM, IRET
      COMMON /PROCESADORES/ NUM

      INTEGER*8 PLAN1,PLAN2
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)

*     FFTW initialization
      DATA1=(0.0,0.0)
      CALL SFFTW_INIT_THREADS(IRET)
      CALL SFFTW_PLAN_WITH_NTHREADS(NUM)

      CALL SFFTW_PLAN_DFT_3D(PLAN1,NX,NY,NZ,DATA1,DATA1,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      CALL SFFTW_PLAN_DFT_3D(PLAN2,NX,NY,NZ,DATA1,DATA1,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
*     End initialization


!$OMP PARALLEL DO SHARED(NX,NY,NZ,DATA1,U1),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)
      DO K=1,NZ 
      DO J=1,NY 
      DO I=1,NX 
        DATA1(I,J,K)=U1(I,J,K)
      END DO 
      END DO 
      END DO

      CALL SFFTW_EXECUTE_DFT(PLAN1,DATA1,DATA1)

**    POISSON EQUATION IN MOMENTUM SPACE
!$OMP PARALLEL DO SHARED(NX,NY,NZ,DATA1,KKK),
!$OMP+            PRIVATE(I,J,K),
!$OMP+            DEFAULT(NONE)                                                    
      DO K=1,NZ 
      DO J=1,NY
      DO I=1,NX
        DATA1(I,J,K)=DATA1(I,J,K)*KKK(I,J,K) 
      END DO 
      END DO 
      END DO

      CALL SFFTW_EXECUTE_DFT(PLAN2,DATA1,DATA1)

!$OMP PARALLEL DO SHARED(NX,NY,NZ,DATA1,POT),PRIVATE(I,J,K)                                                    
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        POT(I,J,K)=REAL(DATA1(I,J,K))/(NX*NY*NZ)                                     
      END DO
      END DO
      END DO

*     FICTICIOUS CELLS FOR FIELD COMPUTATION
      DO K=1,NZ
      DO J=1,NY
          POT(0,J,K)=POT(NX,J,K)
          POT(NX+1,J,K)=POT(1,J,K)
      END DO
      END DO
      DO K=1,NZ
      DO I=0,NX+1
          POT(I,0,K)=POT(I,NY,K)
          POT(I,NY+1,K)=POT(I,1,K)
      END DO
      END DO
      DO J=0,NY+1
      DO I=0,NX+1
          POT(I,J,0)=POT(I,J,NZ)
          POT(I,J,NZ+1)=POT(I,J,1)
      END DO
      END DO


      RETURN
      END
#endif

************************************************************************
      SUBROUTINE MOMENTO(DX,NX,NY,NZ,KKK)
************************************************************************
*     Computes the momentum indices (i.e., the Green function) to
*     solve Poisson's equation in Fourier space, and stores them in
*     the array KKK
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'


      INTEGER NX,NY,NZ,I,J,K
      real PI, DELTA,CONSTA,KX2,KY2,KZ2,DEL2
      INTEGER KX,KY,KZ,NX2,NY2,NZ2

      real KKK(NMAX,NMAY,NMAZ)
      real DX        !size of coarse grid cell

*
      DELTA = DX

      PI=ACOS(-1.0)

**    KERNEL SIN
**    CONSTA=(2.0*PI/(DELTA*NX))*(DELTA/2.0)

      CONSTA=PI/NX
      DEL2=(DELTA*0.5)**2

      NX2=NX/2+2
      NY2=NY/2+2
      NZ2=NZ/2+2

**  IN THIS SECTION MOMENTUM SPACE INDEX ARE COMPUTED FOR THE POISSON EQ.
      DO I=1,NX
        IF (I.GE.NX2) THEN
          KX=I-NX-1
        ELSE
          KX=I-1
        END IF
        KX2=CONSTA*KX
        KX2=(SIN(KX2))**2
      DO J=1,NY
        IF (J.GE.NY2) THEN
          KY=J-NY-1
        ELSE
          KY=J-1
        END IF
        KY2=CONSTA*KY
        KY2=(SIN(KY2))**2
      DO K=1,NZ
        IF (K.GE.NZ2) THEN
          KZ=K-NZ-1
        ELSE
          KZ=K-1
        END IF
        KZ2=CONSTA*KZ
        KZ2=(SIN(KZ2))**2

        KKK(I,J,K)=KX2+KY2+KZ2
      END DO
      END DO
      END DO

**     CENTRAL CONDITION
**     KKK(1,1,1)=0.D0
      KKK(1,1,1)=1.0E30

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
        KKK(I,J,K)=-DEL2/KKK(I,J,K)
      END DO
      END DO
      END DO


      RETURN
      END

***********************************************************************
      SUBROUTINE FOURN(DATA,NN,NDIM,ISIGN)
***********************************************************************
*     Computes the FFT (in NDIM dimensions). These function has been
*     adapted from the following reference:
*---------------------------------------------------------------------*
*     Ref.: Numerical Recipes in FORTRAN 77: Volume 1                 *
*     W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling      *
*     1992, Cambridge University Press                                *
***********************************************************************

      IMPLICIT real(A-H,O-Z)

      INCLUDE 'vortex_parameters.dat'

      DOUBLE PRECISION WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION NN(NDIM),DATA(*)
      NTOT=1
      DO 11 IDIM=1,NDIM
        NTOT=NTOT*NN(IDIM)
11    CONTINUE
      NPREV=1
      DO 18 IDIM=1,NDIM
        N=NN(IDIM)
        NREM=NTOT/(N*NPREV)
        IP1=2*NPREV
        IP2=IP1*N
        IP3=IP2*NREM
        I2REV=1
        DO 14 I2=1,IP2,IP1
          IF(I2.LT.I2REV)THEN
            DO 13 I1=I2,I2+IP1-2,2
              DO 12 I3=I1,IP3,IP2
                I3REV=I2REV+I3-I2
                TEMPR=DATA(I3)
                TEMPI=DATA(I3+1)
                DATA(I3)=DATA(I3REV)
                DATA(I3+1)=DATA(I3REV+1)
                DATA(I3REV)=TEMPR
                DATA(I3REV+1)=TEMPI
12            CONTINUE
13          CONTINUE
          ENDIF
          IBIT=IP2/2
1         IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
            I2REV=I2REV-IBIT
            IBIT=IBIT/2
          GO TO 1
          ENDIF
          I2REV=I2REV+IBIT
14      CONTINUE
        IFP1=IP1
2       IF(IFP1.LT.IP2)THEN
          IFP2=2*IFP1
          THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
          WPR=-2.D0*DSIN(0.5D0*THETA)**2
          WPI=DSIN(THETA)
          WR=1.D0
          WI=0.D0
          DO 17 I3=1,IFP1,IP1
            DO 16 I1=I3,I3+IP1-2,2
              DO 15 I2=I1,IP3,IFP2
                K1=I2
                K2=K1+IFP1
                TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                DATA(K2)=DATA(K1)-TEMPR
                DATA(K2+1)=DATA(K1+1)-TEMPI
                DATA(K1)=DATA(K1)+TEMPR
                DATA(K1+1)=DATA(K1+1)+TEMPI
15            CONTINUE
16          CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
17        CONTINUE
          IFP1=IFP2
        GO TO 2
        ENDIF
        NPREV=N*NPREV
18    CONTINUE
      RETURN
      END

************************************************************************
      SUBROUTINE POTAMR(NL,NX,NY,NZ,DX,NPATCH,PARE,
     &           PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &           PATCHRX,PATCHRY,PATCHRZ,BOR)
************************************************************************
*     Solves Poisson equation for the refinement patches, taking into
*     account the boundary conditions imposed by the coarser cells.
*     It uses a SOR method with Chebyshev acceleration procedure to
*     set the overrelaxation parameter. It uses 3 ficticious cell at
*     each boundary. Boundary conditions are enforced in the ficticious
*     boundary (cells -2 and N+3).
************************************************************************

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER N1,N2,N3,NP1,NP2,NP3,L1,L2,L3
      INTEGER I,J,K,IX,JY,KZ,I1,J2,K3,II,IR,BOR
      INTEGER NX,NY,NZ,NL

      real  U11(-1:NAMRX+2,-1:NAMRY+2,-1:NAMRZ+2,NPALEV)
      real  APOT1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /UAMR/U11,APOT1

      real  U1(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real  POT(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      COMMON /BASE/ U1,POT

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PARE(NPALEV)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV)
      REAL PATCHRY(NPALEV)
      REAL PATCHRZ(NPALEV)

*     Here pot1 is local !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real POT1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3)
*     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real BAS,ERROR,ERRMAX,BASS, PI
      real SSS,DXPA,DX,WWW,ERR,ERRTOT
      real RADIUS,SNOR,RESNOR,PRECIS
      INTEGER CR1,CR2,CR3,MARCA,LOW1,LOW2,MAXIT
      real AAA,BBB,CCC
      COMMON /SOR/ PRECIS,MAXIT

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

      real UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),FUIN
      INTEGER KARE,KR1,KR2,KR3

      INTEGER CR3AMR1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1X(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Y(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      INTEGER CR3AMR1Z(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /CR0CELL/ CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z

*      ---PARALLEL---
      INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
      COMMON /PROCESADORES/ NUM
*     ---------------------------------------------------------------------------

      integer ip,nn1,nn2,nn3,borpotcopy
      real xl,yl,zl,xr,yr,zr,xc,yc,zc

      PI=ACOS(-1.0)

      IR=1
      DXPA=DX/(2.0**IR)

!$OMP PARALLEL DO SHARED(IR,DX,DXPA,POT,APOT1,NX,NY,NZ,NPATCH,
!$OMP+         PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+         U11,PI,PRECIS,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,U1,MAXIT,BOR,
!$OMP+         PATCHRX,PATCHRY,PATCHRZ,RX,RY,RZ),
!$OMP+      PRIVATE(I,N1,N2,N3,L1,L2,L3,NP1,NP2,NP3,JY,KZ,IX,
!$OMP+         I1,J2,K3,II,SSS,BAS,BASS,SNOR,RESNOR,ERR,ERRTOT,
!$OMP+         WWW,RADIUS,ERROR,ERRMAX,CR1,CR2,CR3,POT1,MARCA,
!$OMP+         UBAS,FUIN,KR1,KR2,KR3,XL,YL,ZL,XR,YR,ZR,XC,YC,ZC,
!$OMP+         NN1,NN2,NN3),
!$OMP+      DEFAULT(NONE)
      DO I=1,NPATCH(IR)

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       NP1=NX
       NP2=NY
       NP3=NZ

       ! initialize the potential (real and ficticious cells) by interpolation
       DO KZ=-2,N3+3
       DO JY=-2,N2+3
       DO IX=-2,N1+3

        KR1=CR3AMR1X(IX,JY,KZ,I)
        KR2=CR3AMR1Y(IX,JY,KZ,I)
        KR3=CR3AMR1Z(IX,JY,KZ,I)

        UBAS(1:3,1:3,1:3)=POT(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
        CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
        POT1(IX,JY,KZ)=FUIN

       END DO
       END DO
       END DO

       ! extend the source (ficticious cells) by interpolation
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
        if(ix.lt.1.or.ix.gt.n1.or.
     &     jy.lt.1.or.jy.gt.n2.or.
     &     kz.lt.1.or.kz.gt.n3) then
          ! We first look if we can copy the value from a sibling patch
          marca=0
          xc = rx(ix,i)
          yc = ry(jy,i)
          zc = rz(kz,i)
          do ip=1,npatch(ir) 
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

            u11(ix,jy,kz,i) = u11(kr1,kr2,kr3,ip)
            
            marca=1
            exit
          end do ! ip 

          if (marca.ne.0) cycle

          ! Else, we interpolate the value from a coarser patch
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          UBAS(1:3,1:3,1:3)=U1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
          CALL LININT52D_NEW(IX,JY,KZ,UBAS,FUIN)
          U11(IX,JY,KZ,I)=FUIN
        end if

       END DO
       END DO
       END DO

       SNOR=0.0
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
         SSS=DXPA*DXPA*U11(IX,JY,KZ,I)
         SNOR=SNOR+ABS(SSS)
       END DO
       END DO
       END DO


       RADIUS=COS(PI/((N1+N2+N3+6*BOR)/3.0))
       WWW=1.0
       RESNOR=SNOR
       II=0
       MARCA=0
       DO WHILE (MARCA.EQ.0.OR.II.LT.2) ! SOR iteration
        II=II+1
        RESNOR=0.0

        ERRTOT=-1.0
        DO KZ=1-BOR,N3+BOR
        DO JY=1-BOR,N2+BOR
        DO IX=1-BOR,N1+BOR
         IF (MOD((IX+JY+KZ),2).EQ.MOD((II+1),2)) THEN
*         POISSON EQUATION
          SSS=DXPA*DXPA*U11(IX,JY,KZ,I)

          ERR=POT1(IX,JY,KZ)
          BAS=POT1(IX+1,JY,KZ)+POT1(IX-1,JY,KZ)+POT1(IX,JY+1,KZ)
     &       +POT1(IX,JY-1,KZ)+POT1(IX,JY,KZ+1)+POT1(IX,JY,KZ-1)
     &       -6.0*POT1(IX,JY,KZ) - SSS

          RESNOR=RESNOR+ABS(BAS)

          POT1(IX,JY,KZ)=POT1(IX,JY,KZ)+WWW*BAS/6.0
          IF (ERR.NE.0.0) ERR=POT1(IX,JY,KZ)/ERR - 1.0
          ERRTOT=MAX(ERRTOT,ABS(ERR))
         END IF
        END DO
        END DO
        END DO
        IF (RESNOR.LT.(PRECIS*SNOR)) MARCA=1
*
        IF (ERRTOT.LE.0.1*PRECIS.AND.II.GT.2) MARCA=1

        WWW=1.0/(1.0-0.25*WWW*RADIUS**2)
        IF (II.EQ.1) WWW=1.0/(1.0-0.5*RADIUS**2)
*
        IF (II.GT.MAXIT) MARCA=1
       END DO

      APOT1(-2:N1+3,-2:N2+3,-2:N3+3,I)=POT1(-2:N1+3,-2:N2+3,-2:N3+3)
      END DO


*     OTHER LEVELS
      DO IR=2,NL


      DXPA=DX/(2**IR)

      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(IR,DX,DXPA,POT,APOT1,NX,NY,NZ,NPATCH,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+        U11,PI,PRECIS,LOW1,LOW2,
!$OMP+        CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,PARE,
!$OMP+        RX,RY,RZ,RADX,RADY,RADZ,U1,MAXIT,BOR,PATCHRX,
!$OMP+        PATCHRY,PATCHRZ),
!$OMP+     PRIVATE(I,N1,N2,N3,L1,L2,L3,NP1,NP2,NP3,JY,KZ,IX,MARCA,
!$OMP+        I1,J2,K3,II,SSS,BAS,BASS,RESNOR,SNOR,ERR,ERRTOT,
!$OMP+        WWW,RADIUS,ERROR,ERRMAX,CR1,CR2,CR3,POT1,
!$OMP+        KARE,KR1,KR2,KR3,UBAS,RXBAS,RYBAS,RZBAS,FUIN,
!$OMP+        AAA,BBB,CCC,XL,YL,ZL,XR,YR,ZR,XC,YC,ZC,NN1,
!$OMP+        NN2,NN3),
!$OMP+     DEFAULT(NONE)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       NP1=PATCHNX(PARE(I))
       NP2=PATCHNY(PARE(I))
       NP3=PATCHNZ(PARE(I))


       DO KZ=-2,N3+3
       DO JY=-2,N2+3
       DO IX=-2,N1+3

        KARE=CR3AMR1(IX,JY,KZ,I)
        KR1=CR3AMR1X(IX,JY,KZ,I)
        KR2=CR3AMR1Y(IX,JY,KZ,I)
        KR3=CR3AMR1Z(IX,JY,KZ,I)

        AAA=RX(IX,I)
        BBB=RY(JY,I)
        CCC=RZ(KZ,I)

        IF (KARE.GT.0) THEN
          UBAS(1:3,1:3,1:3)=
     &        APOT1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)

          RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
          RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
          RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
          POT1(IX,JY,KZ)=FUIN
        ELSE
          UBAS(1:3,1:3,1:3)=
     &        POT(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)

          RXBAS(1:3)=RADX(KR1-1:KR1+1)
          RYBAS(1:3)=RADY(KR2-1:KR2+1)
          RZBAS(1:3)=RADZ(KR3-1:KR3+1)

          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
          POT1(IX,JY,KZ)=FUIN
        ENDIF

       END DO
       END DO
       END DO

       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR

        if(ix.lt.1.or.ix.gt.n1.or.
     &     jy.lt.1.or.jy.gt.n2.or.
     &     kz.lt.1.or.kz.gt.n3) then

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
            u11(ix,jy,kz,i) = u11(kr1,kr2,kr3,ip)
            
            marca=1
            exit
          end do ! ip     
          if (marca.ne.0) cycle

          ! Else, we interpolate the value from a coarser patch
          KARE=CR3AMR1(IX,JY,KZ,I)
          KR1=CR3AMR1X(IX,JY,KZ,I)
          KR2=CR3AMR1Y(IX,JY,KZ,I)
          KR3=CR3AMR1Z(IX,JY,KZ,I)

          AAA=RX(IX,I)
          BBB=RY(JY,I)
          CCC=RZ(KZ,I)

          IF (KARE.GT.0) THEN
            UBAS(1:3,1:3,1:3)=
     &        U11(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)

            RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
            RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
            RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)

            CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
            U11(IX,JY,KZ,I)=FUIN
          ELSE
            UBAS(1:3,1:3,1:3)=
     &        U1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)

            RXBAS(1:3)=RADX(KR1-1:KR1+1)
            RYBAS(1:3)=RADY(KR2-1:KR2+1)
            RZBAS(1:3)=RADZ(KR3-1:KR3+1)

            CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                   RXBAS,RYBAS,RZBAS,UBAS,FUIN)
            U11(IX,JY,KZ,I)=FUIN
          ENDIF
        end if

       END DO
       END DO
       END DO


       SNOR=0.0
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
         SSS=DXPA*DXPA*U11(IX,JY,KZ,I)
         SNOR=SNOR+ABS(SSS)
       END DO
       END DO
       END DO


       RADIUS=COS(PI/((N1+N2+N3+6*BOR)/3.0))
*
       WWW=1.0
       II=0
       RESNOR=SNOR
       MARCA=0
       DO WHILE (MARCA.EQ.0.OR.II.LT.2)
        II=II+1
        RESNOR=0.0

        ERRTOT=-1.0
        DO KZ=1-BOR,N3+BOR
        DO JY=1-BOR,N2+BOR
        DO IX=1-BOR,N1+BOR
         IF (MOD((IX+JY+KZ),2).EQ.MOD((II+1),2)) THEN
*         POISSON EQUATION
          SSS=DXPA*DXPA*U11(IX,JY,KZ,I)

          ERR=POT1(IX,JY,KZ)
          BAS=POT1(IX+1,JY,KZ)+POT1(IX-1,JY,KZ)+POT1(IX,JY+1,KZ)
     &       +POT1(IX,JY-1,KZ)+POT1(IX,JY,KZ+1)+POT1(IX,JY,KZ-1)
     &       -6.0*POT1(IX,JY,KZ) - SSS

          RESNOR=RESNOR+ABS(BAS)
          POT1(IX,JY,KZ)=POT1(IX,JY,KZ)+WWW*BAS/6.0
          IF (ERR.NE.0.0) ERR=POT1(IX,JY,KZ)/ERR - 1.0
          ERRTOT=MAX(ERRTOT,ABS(ERR))
         END IF
        END DO
        END DO
        END DO
        IF (RESNOR.LT.(PRECIS*SNOR)) MARCA=1
*
        IF (ERRTOT.LE.0.1*PRECIS.AND.II.GT.2) MARCA=1
        WWW=1.0/(1.0-0.25*WWW*RADIUS**2)
        IF (II.EQ.1) WWW=1.0/(1.0-0.5*RADIUS**2)
*
        IF (II.GT.MAXIT) MARCA=1
       END DO


      APOT1(-2:N1+3,-2:N2+3,-2:N3+3,I)=POT1(-2:N1+3,-2:N2+3,-2:N3+3)

      END DO

*     In a second step at this level, try to solve again poisson but imposing the boundary 
*     conditions and boundary sources from siblings, if any. If not, leave it as it is.

!$OMP PARALLEL DO SHARED(IR,DX,DXPA,POT,APOT1,NX,NY,NZ,NPATCH,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
!$OMP+        U11,PI,PRECIS,LOW1,LOW2,
!$OMP+        CR3AMR1,CR3AMR1X,CR3AMR1Y,CR3AMR1Z,PARE,
!$OMP+        RX,RY,RZ,RADX,RADY,RADZ,U1,MAXIT,BOR,PATCHRX,
!$OMP+        PATCHRY,PATCHRZ),
!$OMP+     PRIVATE(I,N1,N2,N3,L1,L2,L3,NP1,NP2,NP3,JY,KZ,IX,MARCA,
!$OMP+        I1,J2,K3,II,SSS,BAS,BASS,RESNOR,SNOR,ERR,ERRTOT,
!$OMP+        WWW,RADIUS,ERROR,ERRMAX,CR1,CR2,CR3,POT1,
!$OMP+        KARE,KR1,KR2,KR3,UBAS,RXBAS,RYBAS,RZBAS,FUIN,
!$OMP+        AAA,BBB,CCC,XL,YL,ZL,XR,YR,ZR,XC,YC,ZC,NN1,
!$OMP+        NN2,NN3,BORPOTCOPY),
!$OMP+     DEFAULT(NONE)
      DO I=LOW1,LOW2
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)

       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)

       NP1=PATCHNX(PARE(I))
       NP2=PATCHNY(PARE(I))
       NP3=PATCHNZ(PARE(I))


       ! Only need to play with the ficticious cells, for the other ones we 
       ! already have the values from the previous step.
       DO KZ=-2,N3+3
       DO JY=-2,N2+3
       DO IX=-2,N1+3
        if (ix.lt.1.or.ix.gt.n1.or.
     &      jy.lt.1.or.jy.gt.n2.or.
     &      kz.lt.1.or.kz.gt.n3) then

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
           pot1(ix,jy,kz) = apot1(kr1,kr2,kr3,ip)
           
           marca=1
           exit
         end do ! ip   
         if (marca.ne.0) cycle

         ! For the fix boundary cells we will admit copying from the ficticious cells 
         ! of a sibling patch!
         if (ix.eq.-2.or.ix.eq.n1+3.or.
     &       jy.eq.-2.or.jy.eq.n2+3.or.
     &       kz.eq.-2.or.kz.eq.n3+3) then 
           borpotcopy=2
           do ip=low1,low2
              nn1 = patchnx(ip)
              xl = patchrx(ip) - dxpa - borpotcopy*dxpa
              xr = xl + nn1*dxpa + 2*borpotcopy*dxpa
              if (xc.lt.xl.or.xc.gt.xr) cycle
              nn2 = patchny(ip)
              yl = patchry(ip) - dxpa - borpotcopy*dxpa
              yr = yl + nn2*dxpa + 2*borpotcopy*dxpa
              if (yc.lt.yl.or.yc.gt.yr) cycle
              nn3 = patchnz(ip)
              zl = patchrz(ip) - dxpa - borpotcopy*dxpa
              zr = zl + nn3*dxpa + 2*borpotcopy*dxpa
              if (zc.lt.zl.or.zc.gt.zr) cycle
              ! else it is here 
              
              kr1 = int((xc-xl)/dxpa) + 1 - borpotcopy
              kr2 = int((yc-yl)/dxpa) + 1 - borpotcopy
              kr3 = int((zc-zl)/dxpa) + 1 - borpotcopy
              if (kr1.lt.-1.or.kr1.gt.nn1+2) cycle
              if (kr2.lt.-1.or.kr2.gt.nn2+2) cycle
              if (kr3.lt.-1.or.kr3.gt.nn3+2) cycle    
              pot1(ix,jy,kz) = apot1(kr1,kr2,kr3,ip)
              
              marca=1
              exit
           end do ! ip   
           if (marca.ne.0) cycle
         end if

         KARE=CR3AMR1(IX,JY,KZ,I)
         KR1=CR3AMR1X(IX,JY,KZ,I)
         KR2=CR3AMR1Y(IX,JY,KZ,I)
         KR3=CR3AMR1Z(IX,JY,KZ,I)
 
         AAA=RX(IX,I)
         BBB=RY(JY,I)
         CCC=RZ(KZ,I)
 
         IF (KARE.GT.0) THEN
           UBAS(1:3,1:3,1:3)=
     &         APOT1(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1,KARE)
 
           RXBAS(1:3)=RX(KR1-1:KR1+1,KARE)
           RYBAS(1:3)=RY(KR2-1:KR2+1,KARE)
           RZBAS(1:3)=RZ(KR3-1:KR3+1,KARE)
 
           CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                    RXBAS,RYBAS,RZBAS,UBAS,FUIN)
           POT1(IX,JY,KZ)=FUIN
         ELSE
           UBAS(1:3,1:3,1:3)=
     &         POT(KR1-1:KR1+1,KR2-1:KR2+1,KR3-1:KR3+1)
 
           RXBAS(1:3)=RADX(KR1-1:KR1+1)
           RYBAS(1:3)=RADY(KR2-1:KR2+1)
           RZBAS(1:3)=RADZ(KR3-1:KR3+1)
 
           CALL LININT52D_NEW_REAL(AAA,BBB,CCC,
     &                    RXBAS,RYBAS,RZBAS,UBAS,FUIN)
           POT1(IX,JY,KZ)=FUIN
         ENDIF
        end if
       END DO
       END DO
       END DO

       ! The source is already set, no need to do it again.


       SNOR=0.0
       DO KZ=1-BOR,N3+BOR
       DO JY=1-BOR,N2+BOR
       DO IX=1-BOR,N1+BOR
         SSS=DXPA*DXPA*U11(IX,JY,KZ,I)
         SNOR=SNOR+ABS(SSS)
       END DO
       END DO
       END DO


       RADIUS=COS(PI/((N1+N2+N3+6*BOR)/3.0))
*
       WWW=1.0
       II=0
       RESNOR=SNOR
       MARCA=0
       DO WHILE (MARCA.EQ.0.OR.II.LT.2)
        II=II+1
        RESNOR=0.0

        ERRTOT=-1.0
        DO KZ=1-BOR,N3+BOR
        DO JY=1-BOR,N2+BOR
        DO IX=1-BOR,N1+BOR
         IF (MOD((IX+JY+KZ),2).EQ.MOD((II+1),2)) THEN
*         POISSON EQUATION
          SSS=DXPA*DXPA*U11(IX,JY,KZ,I)

          ERR=POT1(IX,JY,KZ)
          BAS=POT1(IX+1,JY,KZ)+POT1(IX-1,JY,KZ)+POT1(IX,JY+1,KZ)
     &       +POT1(IX,JY-1,KZ)+POT1(IX,JY,KZ+1)+POT1(IX,JY,KZ-1)
     &       -6.0*POT1(IX,JY,KZ) - SSS

          RESNOR=RESNOR+ABS(BAS)
          POT1(IX,JY,KZ)=POT1(IX,JY,KZ)+WWW*BAS/6.0
          IF (ERR.NE.0.0) ERR=POT1(IX,JY,KZ)/ERR - 1.0
          ERRTOT=MAX(ERRTOT,ABS(ERR))
         END IF
        END DO
        END DO
        END DO
        IF (RESNOR.LT.(PRECIS*SNOR)) MARCA=1
*
        IF (ERRTOT.LE.0.1*PRECIS.AND.II.GT.2) MARCA=1
        WWW=1.0/(1.0-0.25*WWW*RADIUS**2)
        IF (II.EQ.1) WWW=1.0/(1.0-0.5*RADIUS**2)
*
        IF (II.GT.MAXIT) MARCA=1
       END DO


      APOT1(-2:N1+3,-2:N2+3,-2:N3+3,I)=POT1(-2:N1+3,-2:N2+3,-2:N3+3)

      END DO ! I=LOW1,LOW2
      END DO ! IR=2,NL

      RETURN
      END
