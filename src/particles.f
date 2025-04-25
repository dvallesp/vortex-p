************************************************************************
      SUBROUTINE CREATE_MESH(NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,REFINE_THR,PARCHLIM,BORGRID)
************************************************************************
*     Creates a mesh hierarchy for the given particle distribution
************************************************************************

      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL_MESH
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0
      INTEGER REFINE_THR,PARCHLIM,BORGRID

*     COMMON VARIABLES
      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

*     LOCAL VARIABLES
C      INTEGER PLEV(NDM)
      !REAL,ALLOCATABLE::U1(:,:,:)
      !REAL,ALLOCATABLE::U11(:,:,:,:)
      INTEGER,ALLOCATABLE::CR0(:,:,:)
      INTEGER,ALLOCATABLE::CR01(:,:,:,:)
      INTEGER,ALLOCATABLE::CONTA1(:,:,:)
      INTEGER,ALLOCATABLE::CONTA11(:,:,:,:)
      REAL XL,YL,ZL,DXPA,DYPA,DZPA
      INTEGER I,IX,JY,KZ,REFINE_COUNT,BOR,MIN_PATCHSIZE
      INTEGER INI_EXTENSION,NBIS,IRPA,BORAMR,LOW1,LOW2,IPATCH,IPARE
      INTEGER INMAX(3),INMAX2(2),I1,I2,J1,J2,K1,K2,N1,N2,N3,IR,MARCA
      INTEGER NP1,NP2,NP3,BASINT,NPALEV3,II,JJ,KK

      INTEGER,ALLOCATABLE::LNPATCH(:)
      INTEGER,ALLOCATABLE::LPATCHNX(:,:),LPATCHNY(:,:),LPATCHNZ(:,:)
      INTEGER,ALLOCATABLE::LPATCHX(:,:),LPATCHY(:,:),LPATCHZ(:,:)
      REAL,ALLOCATABLE::LPATCHRX(:,:),LPATCHRY(:,:),LPATCHRZ(:,:)
      INTEGER,ALLOCATABLE::LVAL(:,:)

      REAL,ALLOCATABLE::DDD(:)
      INTEGER,ALLOCATABLE::DDDX(:),DDDY(:),DDDZ(:)
      
      NPATCH(:)=0

      IF (NL_MESH.EQ.0) THEN
       CR0AMR(1:NX,1:NY,1:NZ)=1
       RETURN
      END IF
      
!     hard-coded parameters (for now, at least)
      BORAMR=0
      INI_EXTENSION=0 !initial extension of a patch around a cell (on each direction)
      NPALEV3=(INT(NAMRX/5)**3)+1
      write(*,*) 'NPALEV3=',NPALEV3

      BOR=BORGRID
      MIN_PATCHSIZE=PARCHLIM

C      PLEV=0
C!$OMP PARALLEL DO SHARED(NPART,PLEV,MAP,MASAP),PRIVATE(I),DEFAULT(NONE)
C      DO I=1,SUM(NPART)
C       PLEV(I)=LOG(MAP/MASAP(I)+.5)/LOG(8.)
C      END DO
C      WRITE(*,*) 'Particle levels: min and max values:', MINVAL(PLEV),
C     &           MAXVAL(PLEV)

      XL=-FLOAT(NX)*DX/2.
      YL=-FLOAT(NY)*DY/2.
      ZL=-FLOAT(NZ)*DZ/2.

*     FIRST LEVEL OF REFINEMENT ========================================
      IR=1
      DXPA=DX/(2.0**IR)
      DYPA=DY/(2.0**IR)
      DZPA=DZ/(2.0**IR)
      !ALLOCATE(U1(NMAX,NMAY,NMAZ))
      ALLOCATE(CONTA1(NMAX,NMAY,NMAZ))
      ALLOCATE(CR0(NMAX,NMAY,NMAZ))

!$OMP PARALLEL DO SHARED(CONTA1,CR0,NX,NY,NZ),PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       CONTA1(IX,JY,KZ)=0
       CR0(IX,JY,KZ)=0
      END DO
      END DO
      END DO

      LOW1=SUM(NPART)
!$OMP PARALLEL DO SHARED(NPART,RXPA,RYPA,RZPA,XL,YL,ZL,DX,DY,DZ,
!$OMP+                   NX,NY,NZ,REFINE_THR,LOW1),
!$OMP+            PRIVATE(I,IX,JY,KZ), DEFAULT(NONE)
!$OMP+            REDUCTION(+: CONTA1)
      DO I=1,LOW1
       IX=INT((RXPA(I)-XL)/DX)+1
       JY=INT((RYPA(I)-YL)/DY)+1
       KZ=INT((RZPA(I)-ZL)/DZ)+1
       IF (IX.LT.1) IX=1
       IF (IX.GT.NX) IX=NX
       IF (JY.LT.1) JY=1
       IF (JY.GT.NY) JY=NY
       IF (KZ.LT.1) KZ=1
       IF (KZ.GT.NZ) KZ=NZ

       CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+1

C       IF (PLEV(I).EQ.0) THEN
C        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+1
C       ELSE
C        CONTA1(IX,JY,KZ)=CONTA1(IX,JY,KZ)+REFINE_THR
C       END IF
      END DO

!$OMP PARALLEL DO SHARED(NX,NY,NZ,BOR,CONTA1,CR0),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       IF(IX.LE.BOR.OR.IX.GE.NX-BOR+1.OR.
     &    JY.LE.BOR.OR.JY.GE.NY-BOR+1.OR.
     &    KZ.LE.BOR.OR.KZ.GE.NZ-BOR+1) THEN
         CONTA1(IX,JY,KZ)=0
       END IF
       CR0(IX,JY,KZ)=CONTA1(IX,JY,KZ)
      END DO
      END DO
      END DO

      !WRITE(*,*) 'TOTAL DM MASS: ',SUM(U1*9.18E18)
      WRITE(*,*) 'PARTICLE COUNT CHECK: ',SUM(CONTA1),SUM(NPART)
      WRITE(*,*) 'Max. number of particles in 1 cell',maxval(conta1)

      REFINE_COUNT=COUNT(CR0.GE.REFINE_THR)
      WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

      ALLOCATE(DDD(REFINE_COUNT),DDDX(REFINE_COUNT),
     &         DDDY(REFINE_COUNT),DDDZ(REFINE_COUNT))

      I=0
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       IF (CR0(IX,JY,KZ).GE.REFINE_THR) THEN
        I=I+1
        DDD(I)=FLOAT(CR0(IX,JY,KZ))
        DDDX(I)=IX
        DDDY(I)=JY
        DDDZ(I)=KZ
       END IF
      END DO
      END DO
      END DO

      CALL SORT_CELLS(REFINE_COUNT,DDD,DDDX,DDDY,DDDZ)

      IPATCH=0

      DO I=1,REFINE_COUNT  !----------------------------------------
       IX=DDDX(I)
       JY=DDDY(I)
       KZ=DDDZ(I)
       IF (CR0(IX,JY,KZ).LT.REFINE_THR) CYCLE
       !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

       I1=MAX(IX-INI_EXTENSION,BOR+1)
       I2=MIN(IX+INI_EXTENSION,NX-BOR)
       J1=MAX(JY-INI_EXTENSION,BOR+1)
       J2=MIN(JY+INI_EXTENSION,NY-BOR)
       K1=MAX(KZ-INI_EXTENSION,BOR+1)
       K2=MIN(KZ+INI_EXTENSION,NZ-BOR)

       N1=2*(I2-I1+1)
       N2=2*(J2-J1+1)
       N3=2*(K2-K1+1)
       !NBAS=MAXVAL(N1,N2,N3)
       !NBIS=MINVAL(N1,N2,N3)

       MARCA = 1
       DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MARCA=0
        IF (N1.LE.NAMRX-2.AND.I1.GT.BOR+1) THEN
         IF (COUNT(CR0(I1-1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I1=I1-1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N1.LE.NAMRX-2.AND.I2.LT.NX-BOR) THEN
         !IF (IPATCH.EQ.153) WRITE(*,*) IX,JY,KZ,I1,I2,J1,J2,K1,K2
         IF (COUNT(CR0(I2+1,J1:J2,K1:K2).GE.REFINE_THR).GT.0) THEN
          I2=I2+1
          N1=2*(I2-I1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J1.GT.BOR+1) THEN
         IF (COUNT(CR0(I1:I2,J1-1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J1=J1-1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N2.LE.NAMRY-2.AND.J2.LT.NY-BOR) THEN
         IF (COUNT(CR0(I1:I2,J2+1,K1:K2).GE.REFINE_THR).GT.0) THEN
          J2=J2+1
          N2=2*(J2-J1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K1.GT.BOR+1) THEN
         IF (COUNT(CR0(I1:I2,J1:J2,K1-1).GE.REFINE_THR).GT.0) THEN
          K1=K1-1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

        IF (N3.LE.NAMRZ-2.AND.K2.LT.NZ-BOR) THEN
         IF (COUNT(CR0(I1:I2,J1:J2,K2+1).GE.REFINE_THR).GT.0) THEN
          K2=K2+1
          N3=2*(K2-K1+1)
          MARCA=1
         END IF
        END IF

       END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       NBIS=MIN(N1,N2,N3)
       IF (NBIS.LE.MIN_PATCHSIZE) THEN
C        CR0(I1:I2,J1:J2,K1:K2)=0
        CR0(IX,JY,KZ)=0
       ELSE
        IPATCH=IPATCH+1
*       WRITE(*,*) IPATCH,N1,N2,N3,
*     &             COUNT(CONTA1(I1:I2,J1:J2,K1:K2).GE.REFINE_THR)
        CR0(I1:I2,J1:J2,K1:K2)=-1
        CONTA1(I1:I2,J1:J2,K1:K2)=0

        PATCHNX(IPATCH)=N1
        PATCHNY(IPATCH)=N2
        PATCHNZ(IPATCH)=N3

        PATCHX(IPATCH)=I1
        PATCHY(IPATCH)=J1
        PATCHZ(IPATCH)=K1

        PATCHRX(IPATCH)=RADX(I1)
        PATCHRY(IPATCH)=RADY(J1)
        PATCHRZ(IPATCH)=RADZ(K1)

        PARE(IPATCH)=0
       END IF

      END DO  !-----------------------------------------------------

      DEALLOCATE(DDD,DDDX,DDDY,DDDZ)

      NPATCH(IR)=IPATCH

!$OMP PARALLEL DO SHARED(NX,NY,NZ,CR0,CR0AMR),
!$OMP+            PRIVATE(IX,JY,KZ), DEFAULT(NONE)
      DO KZ=1,NZ
      DO JY=1,NY
      DO IX=1,NX
       IF (CR0(IX,JY,KZ).EQ.-1) THEN
        CR0AMR(IX,JY,KZ)=0
       ELSE
        CR0AMR(IX,JY,KZ)=1
       END IF
      END DO
      END DO
      END DO

      !DEALLOCATE(U1)

      WRITE(*,*) 'At l=',1,', patches:', NPATCH(IR)
      WRITE(*,*) '  --> l=',0,' cells refined:', COUNT(CR0AMR.EQ.0)
      WRITE(*,*) '    ... max num. of particles at an unrefined cell:',
     &           MAXVAL(CONTA1)
      WRITE(*,*) '    ... refinable cells not refined:',
     &           COUNT(CONTA1.GE.REFINE_THR)

      DEALLOCATE(CONTA1)
      DEALLOCATE(CR0)
*     END FIRST LEVEL OF REFINEMENT ====================================

*     START SUBSEQUENT LEVELS OF REFINEMENT ============================
      pare_levels: DO IRPA=1,NL_MESH-1 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       IF (NPATCH(IRPA).EQ.0) THEN
         WRITE(*,*) 'Mesh building stops at level: ', IRPA
         WRITE(*,*) 'There are no more candidate patches'
         EXIT pare_levels
       END IF

       DXPA=DX/(2.0**IRPA)
       DYPA=DY/(2.0**IRPA)
       DZPA=DZ/(2.0**IRPA)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       !WRITE(*,*) IRPA, LOW1,LOW2

       ALLOCATE(CR01(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))
       ALLOCATE(CONTA11(1:NAMRX,1:NAMRY,1:NAMRZ,LOW1:LOW2))

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,CONTA11),
!$OMP+            PRIVATE(IPATCH,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        DO KZ=1,NAMRZ
        DO JY=1,NAMRY
        DO IX=1,NAMRX
         CR01(IX,JY,KZ,IPATCH)=0
         CONTA11(IX,JY,KZ,IPATCH)=0
        END DO
        END DO
        END DO
       END DO

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHRX,PATCHRY,PATCHRZ,DXPA,DYPA,
!$OMP+                   DZPA,PATCHNX,PATCHNY,PATCHNZ,NPART,RXPA,RYPA,
!$OMP+                   RZPA,CONTA11,REFINE_THR,IRPA,BORAMR,CR01),
!$OMP+            PRIVATE(IPATCH,XL,YL,ZL,N1,N2,N3,I,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2 != = = = = = = = = = = = = = = = = = = = = =
        !WRITE(*,*) IPATCH, LOW2
        XL=PATCHRX(IPATCH)-DXPA
        YL=PATCHRY(IPATCH)-DYPA
        ZL=PATCHRZ(IPATCH)-DZPA

        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        DO I=1,SUM(NPART)
         IX=INT((RXPA(I)-XL)/DXPA)+1
         JY=INT((RYPA(I)-YL)/DYPA)+1
         KZ=INT((RZPA(I)-ZL)/DZPA)+1
         IF (IX.GE.1.AND.IX.LE.N1.AND.
     &       JY.GE.1.AND.JY.LE.N2.AND.
     &       KZ.GE.1.AND.KZ.LE.N3) THEN !*****************************
           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+1
          !U1(IX,JY,KZ)=U1(IX,JY,KZ)+MASAP(I)
C          IF (PLEV(I).LE.IRPA) THEN
C           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+1
C          ELSE
C           CONTA11(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)+REFINE_THR
C          END IF
         END IF !*****************************************************
        END DO

        DO KZ=1,N3
        DO JY=1,N2
        DO IX=1,N1
         IF(IX.LE.BORAMR.OR.IX.GE.N1-BORAMR+1.OR.
     &      JY.LE.BORAMR.OR.JY.GE.N2-BORAMR+1.OR.
     &      KZ.LE.BORAMR.OR.KZ.GE.N3-BORAMR+1) THEN
           CONTA11(IX,JY,KZ,IPATCH)=0
         END IF
         CR01(IX,JY,KZ,IPATCH)=CONTA11(IX,JY,KZ,IPATCH)
        END DO
        END DO
        END DO

       END DO != = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       WRITE(*,*) '  --> Max particles at a cell:',
     &            MAXVAL(CR01(:,:,:,LOW1:LOW2))

c       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
c       WRITE(*,*) 'Refinable cells BEFORE cleaning at l=',IRPA,
c     &            REFINE_COUNT

       CALL VEINSGRID_REDUCED(IRPA,NPATCH,PARE,PATCHNX,PATCHNY,
     &      PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,CR01,
     &      CONTA11,LOW1,LOW2)

       REFINE_COUNT=COUNT(CR01(:,:,:,LOW1:LOW2).GE.REFINE_THR)
       WRITE(*,*) '  --> Refinable cells AFTER cleaning:',REFINE_COUNT


       ! mesh creation at the next level
       IR=IRPA+1
       !DXPA=DX/(2.0**IR)
       !DYPA=DY/(2.0**IR)
       !DZPA=DZ/(2.0**IR)

       ALLOCATE(LNPATCH(LOW1:LOW2))
       ALLOCATE(LPATCHNX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHNZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRX(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRY(NPALEV3,LOW1:LOW2))
       ALLOCATE(LPATCHRZ(NPALEV3,LOW1:LOW2))
       ALLOCATE(LVAL(NPALEV3,LOW1:LOW2))

       LNPATCH(:)=0
       LVAL(:,:)=0

c       WRITE(*,*) 'REFINABLE CELLS:', REFINE_COUNT

!$OMP PARALLEL DO SHARED(LOW1,LOW2,CR01,REFINE_THR,NPALEV3,PATCHNX,
!$OMP+                   PATCHNY,PATCHNZ,INI_EXTENSION,BORAMR,CONTA11,
!$OMP+                   MIN_PATCHSIZE,DXPA,DYPA,DZPA,LPATCHNX,
!$OMP+                   LPATCHNY,LPATCHNZ,LPATCHX,LPATCHY,LPATCHZ,
!$OMP+                   LPATCHRX,LPATCHRY,LPATCHRZ,LVAL,PATCHRX,
!$OMP+                   PATCHRY,PATCHRZ),
!$OMP+            PRIVATE(IPARE,REFINE_COUNT,IPATCH,INMAX,IX,JY,KZ,
!$OMP+                    BASINT,NP1,NP2,NP3,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    MARCA,NBIS),
!$OMP+            DEFAULT(NONE)
       DO IPARE=LOW1,LOW2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
        IPATCH=0
        DO WHILE (REFINE_COUNT.GT.0.AND.IPATCH.LT.NPALEV3) !------------
         INMAX=MAXLOC(CR01(:,:,:,IPARE))
         IX=INMAX(1)
         JY=INMAX(2)
         KZ=INMAX(3)
         BASINT=CR01(IX,JY,KZ,IPARE)
         !IF (CONTA1(IX,JY,KZ).LT.REFINE_THR) EXIT

         NP1=PATCHNX(IPARE)
         NP2=PATCHNY(IPARE)
         NP3=PATCHNZ(IPARE)

         I1=MAX(IX-INI_EXTENSION,BORAMR+1)
         I2=MIN(IX+INI_EXTENSION,NP1-BORAMR)
         J1=MAX(JY-INI_EXTENSION,BORAMR+1)
         J2=MIN(JY+INI_EXTENSION,NP2-BORAMR)
         K1=MAX(KZ-INI_EXTENSION,BORAMR+1)
         K2=MIN(KZ+INI_EXTENSION,NP3-BORAMR)

         N1=2*(I2-I1+1)
         N2=2*(J2-J1+1)
         N3=2*(K2-K1+1)
         !NBAS=MAXVAL(N1,N2,N3)
         !NBIS=MINVAL(N1,N2,N3)

         MARCA = 1
         DO WHILE (MARCA.EQ.1) !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          MARCA=0
          IF (N1.LE.NAMRX-2.AND.I1.GT.BORAMR+1) THEN
           IF (COUNT(CR01(I1-1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I1=I1-1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N1.LE.NAMRX-2.AND.I2.LT.NP1-BORAMR) THEN
           IF (COUNT(CR01(I2+1,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            I2=I2+1
            N1=2*(I2-I1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J1.GT.BORAMR+1) THEN
           IF (COUNT(CR01(I1:I2,J1-1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J1=J1-1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N2.LE.NAMRY-2.AND.J2.LT.NP2-BORAMR) THEN
           IF (COUNT(CR01(I1:I2,J2+1,K1:K2,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            J2=J2+1
            N2=2*(J2-J1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K1.GT.BORAMR+1) THEN
           IF (COUNT(CR01(I1:I2,J1:J2,K1-1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            K1=K1-1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

          IF (N3.LE.NAMRZ-2.AND.K2.LT.NP3-BORAMR) THEN
           IF (COUNT(CR01(I1:I2,J1:J2,K2+1,IPARE).GE.REFINE_THR)
     &        .GT.0) THEN
            K2=K2+1
            N3=2*(K2-K1+1)
            MARCA=1
           END IF
          END IF

         END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         NBIS=MIN(N1,N2,N3)
         IF (NBIS.LE.MIN_PATCHSIZE) THEN
C          CR01(I1:I2,J1:J2,K1:K2,IPARE)=0
C          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
           CR01(IX,JY,KZ,IPARE)=0
         ELSE
          IPATCH=IPATCH+1
c          WRITE(*,*) 'new,pare:',IPATCH,IPARE
c          WRITE(*,*) 'N1,N2,N3,refinable:',N1,N2,N3,
c     &             COUNT(CONTA11(I1:I2,J1:J2,K1:K2,IPARE).GE.REFINE_THR)
c          write(*,*) 'x,y,z',i1,j1,k1

          CONTA11(I1:I2,J1:J2,K1:K2,IPARE)=0
          CR01(I1:I2,J1:J2,K1:K2,IPARE)=-1

          LPATCHNX(IPATCH,IPARE)=N1
          LPATCHNY(IPATCH,IPARE)=N2
          LPATCHNZ(IPATCH,IPARE)=N3

          LPATCHX(IPATCH,IPARE)=I1
          LPATCHY(IPATCH,IPARE)=J1
          LPATCHZ(IPATCH,IPARE)=K1

          ! remember that dxpa is the cellsize of the parent!!!
          LPATCHRX(IPATCH,IPARE)=PATCHRX(IPARE)+(I1-1.5)*DXPA
          LPATCHRY(IPATCH,IPARE)=PATCHRY(IPARE)+(J1-1.5)*DYPA
          LPATCHRZ(IPATCH,IPARE)=PATCHRZ(IPARE)+(K1-1.5)*DZPA

          LVAL(IPATCH,IPARE)=BASINT
         END IF

         REFINE_COUNT=COUNT(CR01(:,:,:,IPARE).GE.REFINE_THR)
         !WRITE(*,*) REFINE_COUNT
        END DO  !-------------------------------------------------------
       END DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       IPATCH=LOW2
       DO WHILE(COUNT(LVAL(:,:).GT.0).GT.0.AND.IPATCH.LT.NPALEV)
        INMAX2=MAXLOC(LVAL)
        I=INMAX2(1)
        IPARE=LOW1-1+INMAX2(2)
C        WRITE(*,*) LVAL(I,IPARE)
        LVAL(I,IPARE)=0

        IPATCH=IPATCH+1
        IF (IPATCH.GT.NPALEV) EXIT

        PATCHNX(IPATCH)=LPATCHNX(I,IPARE)
        PATCHNY(IPATCH)=LPATCHNY(I,IPARE)
        PATCHNZ(IPATCH)=LPATCHNZ(I,IPARE)

        PATCHX(IPATCH)=LPATCHX(I,IPARE)
        PATCHY(IPATCH)=LPATCHY(I,IPARE)
        PATCHZ(IPATCH)=LPATCHZ(I,IPARE)

        PATCHRX(IPATCH)=LPATCHRX(I,IPARE)
        PATCHRY(IPATCH)=LPATCHRY(I,IPARE)
        PATCHRZ(IPATCH)=LPATCHRZ(I,IPARE)

        PARE(IPATCH)=IPARE
       END DO

       NPATCH(IR)=IPATCH-SUM(NPATCH(0:IR-1))
       IF (SUM(NPATCH).GE.NPALEV) STOP 'NPALEV too small'

       DEALLOCATE(LNPATCH,LPATCHNX,LPATCHNY,LPATCHNZ,LPATCHRX,LPATCHRY,
     &            LPATCHRZ,LPATCHX,LPATCHY,LPATCHZ,LVAL)

       ! still needing to compute cr0amr1(:,:,:,low1:low2)
       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        DO KZ=1,NAMRZ
        DO JY=1,NAMRY
        DO IX=1,NAMRX
         CR0AMR1(IX,JY,KZ,IPATCH)=1
        END DO
        END DO
        END DO
       END DO

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,PATCHX,
!$OMP+                   PATCHY,PATCHZ,PARE,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,IPARE,I1,I2,J1,J2,K1,K2,N1,N2,N3,
!$OMP+                    IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2
        IPARE=PARE(IPATCH)
        I1=PATCHX(IPATCH)
        J1=PATCHY(IPATCH)
        K1=PATCHZ(IPATCH)
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        I2=I1+N1/2-1
        J2=J1+N2/2-1
        K2=K1+N3/2-1
        DO KZ=K1,K2
        DO JY=J1,J2
        DO IX=I1,I2
         CR0AMR1(IX,JY,KZ,IPARE)=0
        END DO
        END DO
        END DO
       END DO

       CALL VEINSGRID_CR0AMR(IRPA,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                       PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,
     &                       PATCHRZ)

       LOW1=SUM(NPATCH(0:IRPA-1))+1
       LOW2=SUM(NPATCH(0:IRPA))
       WRITE(*,*) 'At l=',IR,', patches:', NPATCH(IR)
       WRITE(*,*) '  --> l=',IRPA,' cells refined:',
     &            COUNT(CR0AMR1(:,:,:,LOW1:LOW2).EQ.0)

       K1=0
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1,
!$OMP+                   CONTA11),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,K2),
!$OMP+            REDUCTION(MAX:K1)
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2 
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        K2=MAXVAL(CONTA11(1:N1,1:N2,1:N3,IPATCH)*
     &            CR0AMR1(1:N1,1:N2,1:N3,IPATCH))
        K1=MAX(K1,K2)
       END DO

       WRITE(*,*) '    ... max num. of particles at an unrefined cell:',
     &               K1

       K1=0
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1,
!$OMP+                   CONTA11,REFINE_THR),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,K2),
!$OMP+            REDUCTION(+:K1)
!$OMP+            DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2 
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
 
        K2=COUNT(CONTA11(1:N1,1:N2,1:N3,IPATCH)*
     &           CR0AMR1(1:N1,1:N2,1:N3,IPATCH).GE.REFINE_THR)
        K1=K1+K2
       END DO

       WRITE(*,*) '    ... refinable cells not refined:', K1

       DEALLOCATE(CR01, CONTA11)

      END DO pare_levels !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
*     END SUBSEQUENT LEVELS OF REFINEMENT

      IF (NPATCH(IR).EQ.0) IR=IR-1 !The last non-empty level, in any case
      LOW1=SUM(NPATCH(0:IR-1))+1
      LOW2=SUM(NPATCH(0:IR))
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ)
      DO IPATCH=LOW1,LOW2
       N1=PATCHNX(IPATCH)
       N2=PATCHNY(IPATCH)
       N3=PATCHNZ(IPATCH)
       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1
        CR0AMR1(IX,JY,KZ,IPATCH)=1
       END DO
       END DO
       END DO
      END DO

      RETURN
      END

************************************************************************
      subroutine create_mesh_octree(NX,NY,NZ,NL_MESH,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,REFINE_THR,PARCHLIM,BORGRID)
************************************************************************
       
      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL_MESH
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0
      INTEGER REFINE_THR,PARCHLIM,BORGRID

*     COMMON VARIABLES
      real DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

*     Local variables
      integer nmax_parent, cellsbox, nmax_child
      integer ii, jj, kk, ix, jy, kz, i1, j1, k1, i2, j2, k2
      integer n1, n2, n3, nbis, low1, low2, ir, ipatch, irpa
      integer np1, np2, np3, i, conta2(2,2,2), nparttot
      real xl, yl, zl, dxpa, dypa, dzpa, xr, yr, zr, xm, ym, zm
      real xp, yp, zp
      integer, allocatable :: conta1(:,:,:)

      integer, allocatable :: lnpatch(:)
      integer, allocatable :: lpatchnx(:,:),lpatchny(:,:),lpatchnz(:,:)
      integer, allocatable :: lpatchx(:,:), lpatchy(:,:), lpatchz(:,:)
      real, allocatable :: lpatchrx(:,:), lpatchry(:,:), lpatchrz(:,:)
      integer lipatch

      npatch(:)=0

      if (nl_mesh.eq.0) then 
       cr0amr(1:nx,1:ny,1:nz)=1
       return 
      end if

      nmax_child = namrx
      nmax_parent = nmax_child / 2 
      if (mod(nx,nmax_parent).ne.0) then
       write(*,*) 'nx not divisible by nmax_parent'
       stop
      end if

!$omp parallel do shared(nx,ny,nz,cr0amr),
!$omp+            private(ix,jy,kz), default(none)
      do kz=1,nz
      do jy=1,ny
      do ix=1,nx
       cr0amr(ix,jy,kz) = 1
      end do
      end do
      end do

*     FIRST LEVEL OF REFINEMENT ========================================
      n1 = nx / nmax_parent 
      n2 = ny / nmax_parent
      n3 = nz / nmax_parent

      cellsbox = nmax_parent ** 3

      xl = -float(nx)*dx/2.
      yl = -float(ny)*dy/2.
      zl = -float(nz)*dz/2.

      ir = 1
      ! Caution: here dxpa are the sizes of the blocks candidates for refinement
      dxpa = nmax_parent*dx
      dypa = nmax_parent*dy
      dzpa = nmax_parent*dz
      
      allocate(conta1(1:n1,1:n2,1:n3))
      conta1(:,:,:) = 0

      nparttot = sum(npart)
!$omp parallel do shared(nparttot,rxpa,rypa,rzpa,
!$omp+                   xl,yl,zl,dxpa,dypa,dzpa,
!$omp+                   n1,n2,n3),
!$omp+            private(ix,jy,kz,i), 
!$omp+            default(none)
!$omp+            reduction(+:conta1)
      do i = 1, nparttot
       ix = int((rxpa(i)-xl)/dxpa) + 1
       jy = int((rypa(i)-yl)/dypa) + 1
       kz = int((rzpa(i)-zl)/dzpa) + 1
       if (ix.ge.1.and.ix.le.n1.and.
     &     jy.ge.1.and.jy.le.n2.and.
     &     kz.ge.1.and.kz.le.n3) then
          conta1(ix,jy,kz) = conta1(ix,jy,kz) + 1
       end if 
      end do

C      write(*,*) 'maxval minval conta1', maxval(conta1), minval(conta1)
C      write(*,*) 'sum conta1', sum(conta1)


      ! Caution: here dxpa are the cell sizes at level 0 (base grid)
      dxpa = dx / 2.0 ** (ir-1)
      dypa = dy / 2.0 ** (ir-1)
      dzpa = dz / 2.0 ** (ir-1)

      ipatch = 0

      do ix = 1, n1 
      do jy = 1, n2
      do kz = 1, n3
       if (conta1(ix,jy,kz).ge.refine_thr*cellsbox) then
        ! Accept the patch 
        ipatch = ipatch + 1
        
        patchnx(ipatch) = nmax_child
        patchny(ipatch) = nmax_child
        patchnz(ipatch) = nmax_child 

        i1 = (ix-1) * nmax_parent + 1
        j1 = (jy-1) * nmax_parent + 1
        k1 = (kz-1) * nmax_parent + 1
        i2 = i1 + nmax_parent - 1
        j2 = j1 + nmax_parent - 1
        k2 = k1 + nmax_parent - 1

        patchx(ipatch) = i1
        patchy(ipatch) = j1
        patchz(ipatch) = k1

        patchrx(ipatch) = xl + (i1 - 0.5) * dxpa
        patchry(ipatch) = yl + (j1 - 0.5) * dypa
        patchrz(ipatch) = zl + (k1 - 0.5) * dzpa

        pare(ipatch) = 0

        cr0amr(i1:i2,j1:j2,k1:k2) = 0

C        write(*,*) 'accepted patch:', ipatch, 'conta1:', i1,j1,k1,
C     &                   i2,j2,k2,patchrx(ipatch),patchry(ipatch),
C     &                   patchrz(ipatch)
       end if
      end do 
      end do 
      end do

      deallocate(conta1)

      npatch(ir) = ipatch
      
      write(*,*) 'At l=',ir,', patches:', npatch(ir)
      ix = count(cr0amr(:,:,:).eq.0)
      write(*,*) '  --> l=',ir-1,' cells refined:', ix
            write(*,*) '    ... fraction ', float(ix) / (float(nx)**3.),
     &            'wrt. the whole domain'

*     SUBSEQUENT LEVELS OF REFINEMENT ================================
      do irpa=1, nl_mesh-1 
       if (npatch(irpa).eq.0) then
         write(*,*) 'Mesh building stops at level: ', irpa
         write(*,*) 'There are no more candidate patches'
         exit
       end if

       dxpa = dx / 2.0 ** irpa
       dypa = dy / 2.0 ** irpa
       dzpa = dz / 2.0 ** irpa

       low1 = sum(npatch(0:irpa-1)) + 1
       low2 = sum(npatch(0:irpa))

       allocate(lnpatch(low1:low2))
       allocate(lpatchnx(8, low1:low2))
       allocate(lpatchny(8, low1:low2))
       allocate(lpatchnz(8, low1:low2))
       allocate(lpatchx(8, low1:low2))
       allocate(lpatchy(8, low1:low2))
       allocate(lpatchz(8, low1:low2))
       allocate(lpatchrx(8, low1:low2))
       allocate(lpatchry(8, low1:low2))
       allocate(lpatchrz(8, low1:low2))

!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,
!$omp+                   patchx,patchy,patchz,patchrx,patchry,
!$omp+                   patchrz,irpa,cr0amr1,refine_thr,npart,
!$omp+                   rxpa,rypa,rzpa,dxpa,dypa,dzpa,cellsbox,
!$omp+                   nmax_child,lpatchnx,lpatchny,lpatchnz,
!$omp+                   lpatchx,lpatchy,lpatchz,lpatchrx,
!$omp+                   lpatchry,lpatchrz,lnpatch,nmax_parent,
!$omp+                   nparttot),
!$omp+            private(ipatch,ix,jy,kz,xl,yl,zl,xr,yr,zr,
!$omp+                    np1,np2,np3,xm,ym,zm,i,conta2,xp,yp,zp,
!$omp+                    lipatch,i1,j1,k1,i2,j2,k2),
!$omp+            default(none)
       do ipatch = low1, low2 
        np1 = patchnx(ipatch)
        np2 = patchny(ipatch)
        np3 = patchnz(ipatch)

        xl = patchrx(ipatch) - dxpa
        yl = patchry(ipatch) - dypa
        zl = patchrz(ipatch) - dzpa

        xr = xl + np1 * dxpa
        yr = yl + np2 * dypa
        zr = zl + np3 * dzpa

C        write(*,*) 'ip,xl,yl,zl,xr,yr,zr:', ipatch,xl,yl,zl,xr,yr,zr

        xm = (xl + xr) / 2.0
        ym = (yl + yr) / 2.0
        zm = (zl + zr) / 2.0

        conta2(1:2,1:2,1:2) = 0

        do i = 1, nparttot
          xp = rxpa(i)
          yp = rypa(i)
          zp = rzpa(i)

          if (xp.lt.xl) cycle 
          if (xp.ge.xr) cycle
          if (yp.lt.yl) cycle
          if (yp.ge.yr) cycle
          if (zp.lt.zl) cycle
          if (zp.ge.zr) cycle

          !write(*,*) 'ipatch:', ipatch, 'xp,yp,zp:', xp,yp,zp

          ix = merge(1, 2, xp .lt. xm)
          jy = merge(1, 2, yp .lt. ym)
          kz = merge(1, 2, zp .lt. zm)

          conta2(ix,jy,kz) = conta2(ix,jy,kz) + 1
        end do

C        write(*,*) 'conta2 ipatch:', ipatch, conta2,sum(conta2)

        cr0amr1(1:np1, 1:np2, 1:np3, ipatch) = 1

        lipatch = 0
        do ix = 1, 2
        do jy = 1, 2
        do kz = 1, 2
C         write(*,*) 'ip,ix,jy,kz:', ipatch,ix, jy, kz, conta2(ix,jy,kz),
C     &                         refine_thr*cellsbox
         if (conta2(ix,jy,kz).ge.refine_thr*cellsbox) then
          lipatch = lipatch + 1

          lpatchnx(lipatch,ipatch) = nmax_child
          lpatchny(lipatch,ipatch) = nmax_child
          lpatchnz(lipatch,ipatch) = nmax_child

          i1 = (ix-1) * nmax_parent + 1
          j1 = (jy-1) * nmax_parent + 1
          k1 = (kz-1) * nmax_parent + 1

          i2 = i1 + nmax_parent - 1
          j2 = j1 + nmax_parent - 1
          k2 = k1 + nmax_parent - 1

          lpatchx(lipatch,ipatch) = i1
          lpatchy(lipatch,ipatch) = j1
          lpatchz(lipatch,ipatch) = k1
          
          lpatchrx(lipatch,ipatch) = xl + (i1 - 0.5) * dxpa
          lpatchry(lipatch,ipatch) = yl + (j1 - 0.5) * dypa
          lpatchrz(lipatch,ipatch) = zl + (k1 - 0.5) * dzpa

          cr0amr1(i1:i2,j1:j2,k1:k2,ipatch) = 0
         end if
        end do 
        end do 
        end do
        
        lnpatch(ipatch) = lipatch 
       end do

       ipatch = low2
       do i = low1, low2 
        if (lnpatch(i).eq.0) cycle
        
        do lipatch = 1, lnpatch(i)
         ipatch = ipatch + 1

         patchnx(ipatch) = lpatchnx(lipatch,i)
         patchny(ipatch) = lpatchny(lipatch,i)
         patchnz(ipatch) = lpatchnz(lipatch,i)

         patchx(ipatch) = lpatchx(lipatch,i)
         patchy(ipatch) = lpatchy(lipatch,i)
         patchz(ipatch) = lpatchz(lipatch,i)

         patchrx(ipatch) = lpatchrx(lipatch,i)
         patchry(ipatch) = lpatchry(lipatch,i)
         patchrz(ipatch) = lpatchrz(lipatch,i)

         pare(ipatch) = i
        end do

        npatch(irpa+1) = ipatch - low2 
        if (ipatch.gt.NPALEV) stop 'NPALEV too small'
       end do

       deallocate(lnpatch)
       deallocate(lpatchnx, lpatchny, lpatchnz)
       deallocate(lpatchx, lpatchy, lpatchz)
       deallocate(lpatchrx, lpatchry, lpatchrz)

       ix = count(cr0amr1(:,:,:,low1:low2).eq.0)
       write(*,*) 'At l=',irpa+1,', patches:', npatch(irpa+1)
       write(*,*) '  --> l=',irpa,' cells refined:', ix          
       write(*,*) '    ... which is a fraction of:',
     &            float(ix) / float(nmax_child**3 * npatch(irpa)),
     &            'wrt. previous level'
       write(*,*) '    ... ', float(ix) / (float(nx)**3. * 8.**irpa),
     &            'wrt. the whole domain'



      end do

      ir = irpa + 1
      if (npatch(ir).eq.0) ir = ir - 1 !The last non-empty level, in any case

      low1 = sum(npatch(0:ir-1)) + 1
      low2 = sum(npatch(0:ir))
!$omp parallel do shared(low1,low2,patchnx,patchny,patchnz,cr0amr1),
!$omp+            private(ipatch,n1,n2,n3,ix,jy,kz)
!$omp+            default(none)
      do ipatch = low1, low2
       n1 = patchnx(ipatch)
       n2 = patchny(ipatch)
       n3 = patchnz(ipatch)
       do kz=1,n3
       do jy=1,n2
       do ix=1,n1
        cr0amr1(ix,jy,kz,ipatch) = 1
       end do
       end do
       end do
      end do

      return 
      end 

************************************************************************
      SUBROUTINE PLACE_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,parchlim)
************************************************************************

      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0
      integer parchlim

*     Local variables
      INTEGER IPATCH,IX,JY,KZ,BORAMR,IP,IR,LOW1,LOW2,N1,N2,N3,MARCA
      INTEGER IRMAX,NUMPARTLOOP
      REAL BASX,BASY,BASZ,DXPA,DYPA,DZPA
      REAL XL(NPALEV),XR(NPALEV),YL(NPALEV),YR(NPALEV),ZL(NPALEV),
     &     ZR(NPALEV)

*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      ALLOCATE(LIHAL(PARTI), LIHAL_IX(PARTI), LIHAL_JY(PARTI),
     &         LIHAL_KZ(PARTI))

      BORAMR=1
      if (parchlim.lt.0) boramr=0

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DXPA=DX/2.0**IR
       DYPA=DY/2.0**IR
       DZPA=DZ/2.0**IR
       DO IPATCH=LOW1,LOW2
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        XL(IPATCH)=PATCHRX(IPATCH)+(BORAMR-1)*DXPA
        YL(IPATCH)=PATCHRY(IPATCH)+(BORAMR-1)*DYPA
        ZL(IPATCH)=PATCHRZ(IPATCH)+(BORAMR-1)*DZPA

        XR(IPATCH)=XL(IPATCH)+(N1-2*BORAMR)*DXPA
        YR(IPATCH)=YL(IPATCH)+(N2-2*BORAMR)*DYPA
        ZR(IPATCH)=ZL(IPATCH)+(N3-2*BORAMR)*DZPA
       END DO
      END DO

      NUMPARTLOOP=SUM(NPART(0:NL))
!$OMP PARALLEL DO SHARED(NPART,NL,NPATCH,XL,XR,YL,YR,ZL,ZR,RXPA,RYPA,
!$OMP+                   RZPA,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,BORAMR,
!$OMP+                   LADO0,DX,DY,DZ,NX,NY,NZ,KERNEL,NUMPARTLOOP),
!$OMP+            PRIVATE(IP,MARCA,IR,LOW1,LOW2,IPATCH,DXPA,DYPA,DZPA,
!$OMP+                    BASX,BASY,BASZ,IX,JY,KZ,IRMAX),
!$OMP+            DEFAULT(NONE)
      DO IP=1,NUMPARTLOOP
       LIHAL(IP)=-1
       MARCA=0

       ! Look for the finest patch tah contains the particle,
       ! but not with a finer resolution than the SPH kernel length.
       IF (KERNEL(IP).GT.DX) THEN
        IRMAX=0 
       ELSE
        IRMAX=INT(LOG(DX/KERNEL(IP))/LOG(2.0))+1
        IRMAX=MIN(IRMAX,NL)
       END IF

       levels: DO IR=IRMAX,1,-1
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        patches: DO IPATCH=LOW1,LOW2
         BASX = (XR(IPATCH) - RXPA(IP)) * (RXPA(IP) - XL(IPATCH))
         IF (BASX.GT.0) THEN
          BASY = (YR(IPATCH) - RYPA(IP)) * (RYPA(IP) - YL(IPATCH))
          IF (BASY.GT.0) THEN
           BASZ = (ZR(IPATCH) - RZPA(IP)) * (RZPA(IP) - ZL(IPATCH))
           IF (BASZ.GT.0) THEN
            DXPA=DX/2.0**IR
            DYPA=DY/2.0**IR
            DZPA=DZ/2.0**IR
            LIHAL(IP)=IPATCH
            LIHAL_IX(IP)=BORAMR+1+INT((RXPA(IP)-XL(IPATCH))/DXPA)
            LIHAL_JY(IP)=BORAMR+1+INT((RYPA(IP)-YL(IPATCH))/DYPA)
            LIHAL_KZ(IP)=BORAMR+1+INT((RZPA(IP)-ZL(IPATCH))/DZPA)
            MARCA=1
            EXIT levels
           END IF
          END IF
         END IF
        END DO patches
       END DO levels

       IF (MARCA.EQ.0) THEN
        LIHAL(IP)=0
        IX=INT((RXPA(IP)+0.5*LADO0)/DX)+1
        JY=INT((RYPA(IP)+0.5*LADO0)/DY)+1
        KZ=INT((RZPA(IP)+0.5*LADO0)/DZ)+1
        IF (IX.GT.NX) IX=NX
        IF (JY.GT.NY) JY=NY
        IF (KZ.GT.NZ) KZ=NZ
        LIHAL_IX(IP)=IX
        LIHAL_JY(IP)=JY
        LIHAL_KZ(IP)=KZ
       END IF


      END DO !IP

      WRITE(*,*) 'At level', 0, '-->',
     &            COUNT(LIHAL(1:SUM(NPART(0:NL))).EQ.0)
      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       WRITE(*,*) 'At level', IR, '-->', COUNT(LIHAL(1:SUM(NPART(0:NL)))
     &                  .GE.LOW1.AND.LIHAL(1:SUM(NPART(0:NL))).LE.LOW2)
      END DO
      WRITE(*,*) 'Particles not located -->',
     &            COUNT(LIHAL(1:SUM(NPART(0:NL))).EQ.-1)

      RETURN
      END

#ifdef output_particles 
#if output_particles == 1
************************************************************************
      SUBROUTINE ERROR_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0)
************************************************************************
      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0

*     Output variables
      REAL ERR(PARTI)

*     Local variables
      INTEGER IPATCH,IX,JY,KZ,LOW1,LOW2,IR,IP,II
      REAL BASX,BASY,BASZ,BASXX,BASYY,BASZZ,BAS,BAS2,DXPA,DYPA,DZPA
      REAL UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),AAA,BBB,CCC
      REAL PATCHLEV(NPALEV)

*     COMMON VARIABLES
      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      real RX(-2:NAMRX+3,NPALEV),RY(-2:NAMRX+3,NPALEV),
     &     RZ(-2:NAMRX+3,NPALEV),RMX(-2:NAMRX+3,NPALEV),
     &     RMY(-2:NAMRX+3,NPALEV),RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/ RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

      CHARACTER*5 ITER_STRING 
      CHARACTER*200 FILENOM 
      WRITE(ITER_STRING,'(I5.5)') ITER

      DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO IPATCH=LOW1,LOW2
        PATCHLEV(IPATCH)=IR
       END DO
      END DO

      LOW1=SUM(NPART(0:NL))
      DO IP=1,LOW1
       IPATCH=LIHAL(IP)
       IX=LIHAL_IX(IP)
       JY=LIHAL_JY(IP)
       KZ=LIHAL_KZ(IP)
       AAA=RXPA(IP)
       BBB=RYPA(IP)
       CCC=RZPA(IP)
       BASXX=U2DM(IP)
       BASYY=U3DM(IP)
       BASZZ=U4DM(IP)
       BAS2=BASXX**2+BASYY**2+BASZZ**2

       IF (IPATCH.GT.0) THEN
        RXBAS=RX(IX-1:IX+1,IPATCH)
        RYBAS=RY(JY-1:JY+1,IPATCH)
        RZBAS=RZ(KZ-1:KZ+1,IPATCH)

        UBAS(1:3,1:3,1:3)=U12(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH)
        CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                          BASX)

        UBAS(1:3,1:3,1:3)=U13(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH)
        CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                          BASY)

        UBAS(1:3,1:3,1:3)=U14(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH)
        CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                          BASZ)

       ELSE
        IF (IX.GT.1.AND.IX.LT.NX.AND.
     &      JY.GT.1.AND.JY.LT.NY.AND.
     &      KZ.GT.1.AND.KZ.LT.NZ) THEN
         RXBAS=RADX(IX-1:IX+1)
         RYBAS=RADY(JY-1:JY+1)
         RZBAS=RADZ(KZ-1:KZ+1)

         UBAS(1:3,1:3,1:3)=U2(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1)
         CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                           BASX)

         UBAS(1:3,1:3,1:3)=U3(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1)
         CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                           BASY)

         UBAS(1:3,1:3,1:3)=U4(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1)
         CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                           BASZ)
        ELSE
         BASX=U2(IX,JY,KZ)
         BASY=U3(IX,JY,KZ)
         BASZ=U4(IX,JY,KZ)
        END IF
       END IF

       BAS=0.0
       BAS=BAS+((BASXX**2/BAS2)*ABS((BASX-BASXX))/(ABS(BASXX)+1E-5))**2
       BAS=BAS+((BASYY**2/BAS2)*ABS((BASY-BASYY))/(ABS(BASYY)+1E-5))**2
       BAS=BAS+((BASZZ**2/BAS2)*ABS((BASZ-BASZZ))/(ABS(BASZZ)+1E-5))**2
       ERR(IP)=SQRT(BAS)

      END DO

      WRITE(*,*) 'Velocity errors:', MINVAL(ERR(1:SUM(NPART(0:NL)))),
     &                               MAXVAL(ERR(1:SUM(NPART(0:NL))))

*     Save errors to a file

      FILENOM='output_files/error-particles'//ITER_STRING
      OPEN(98,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(98) (ERR(IP),IP=1,SUM(NPART(0:NL)))
      CLOSE(98)

      RETURN
      END
#endif 
#endif

#ifdef output_particles
#if output_particles == 1
************************************************************************
      SUBROUTINE GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,
     &            VAR0,VAR1,VARPART)
************************************************************************

      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0
      REAL VAR0(NMAX,NMAY,NMAZ)
      REAL VAR1(NAMRX,NAMRY,NAMRZ,NPALEV)

*     Output variables
      REAL VARPART(PARTI)

*     Local variables
      INTEGER IPATCH,IX,JY,KZ,LOW1,LOW2,IR,IP,II,n1,n2,n3
      REAL FUIN,DXPA,DYPA,DZPA
      REAL UBAS(3,3,3),RXBAS(3),RYBAS(3),RZBAS(3),AAA,BBB,CCC
      REAL PATCHLEV(NPALEV)

*     COMMON VARIABLES
      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      real RX(-2:NAMRX+3,NPALEV),RY(-2:NAMRX+3,NPALEV),
     &     RZ(-2:NAMRX+3,NPALEV),RMX(-2:NAMRX+3,NPALEV),
     &     RMY(-2:NAMRX+3,NPALEV),RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      real  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/ RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      LOW1=SUM(NPART(0:NL))
!$OMP PARALLEL DO SHARED(NL,NPART,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ,RX,
!$OMP+                   RY,RZ,NX,NY,NZ,RADX,RADY,RADZ,VAR0,VAR1,RXPA,
!$OMP+                   RYPA,RZPA,VARPART,LOW1,patchnx,patchny,
!$omp+                   patchnz),
!$OMP+            PRIVATE(IP,IPATCH,IX,JY,KZ,AAA,BBB,CCC,RXBAS,RYBAS,
!$OMP+                    RZBAS,UBAS,FUIN,n1,n2,n3),
!$OMP+            DEFAULT(NONE)
      DO IP=1,LOW1
       IPATCH=LIHAL(IP)
       IX=LIHAL_IX(IP)
       JY=LIHAL_JY(IP)
       KZ=LIHAL_KZ(IP)
       AAA=RXPA(IP)
       BBB=RYPA(IP)
       CCC=RZPA(IP)

       IF (IPATCH.GT.0) THEN
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)

        if (ix.gt.1.and.ix.lt.n1.and.
     &      jy.gt.1.and.jy.lt.n2.and.
     &      kz.gt.1.and.kz.lt.n3) then
          RXBAS=RX(IX-1:IX+1,IPATCH)
          RYBAS=RY(JY-1:JY+1,IPATCH)
          RZBAS=RZ(KZ-1:KZ+1,IPATCH)
  
          UBAS(1:3,1:3,1:3)=VAR1(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1,IPATCH)
          CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                            FUIN)
        else
          fuin = var1(ix,jy,kz,ipatch)
        end if
       ELSE
        IF (IX.GT.1.AND.IX.LT.NX.AND.
     &      JY.GT.1.AND.JY.LT.NY.AND.
     &      KZ.GT.1.AND.KZ.LT.NZ) THEN
         RXBAS=RADX(IX-1:IX+1)
         RYBAS=RADY(JY-1:JY+1)
         RZBAS=RADZ(KZ-1:KZ+1)

         UBAS(1:3,1:3,1:3)=VAR0(IX-1:IX+1,JY-1:JY+1,KZ-1:KZ+1)
         CALL LININT52D_NEW_REAL(AAA,BBB,CCC,RXBAS,RYBAS,RZBAS,UBAS,
     &                           FUIN)
        ELSE
         FUIN=VAR0(IX,JY,KZ)
        END IF
       END IF

       VARPART(IP)=FUIN

      END DO

      RETURN
      END
#endif
#endif

************************************************************************
      SUBROUTINE KERNEL_FUNC(N,N2,W,DIST)
************************************************************************
*     DIST contains initially the distance (particle to cell), and it is
*     updated with the (unnormalised) value of the kernel
      IMPLICIT NONE
      INTEGER N,N2 ! N is the dimension of the array dist;
                   ! N2, the actual number of particles filled in
      REAL W,DIST(N) ! W is the total width of the kernel 
      REAL H

      REAL DISTS
      INTEGER I

#ifdef ikernel 
#if ikernel == 0 
!     Cubic spline kernel (M4)
      H=W/2.0
      DO I=1,N2
       DISTS=DIST(I)/H
       IF (DISTS.LE.1.0) THEN
        DIST(I)=1.0-1.5*DISTS**2*(1-0.5*DISTS)
       ELSE IF (DISTS.LE.2.0) THEN
        DIST(I)=0.25*(2-DISTS)**3
       ELSE
        DIST(I)=0.0
       END IF
      END DO
#elif ikernel == 1
!     Wendland C4 kernel
      H=W/2.0
      DO I=1,N2
       DISTS=DIST(I)/H
       IF (DISTS.LE.2.0) THEN
        DIST(I)=(1. - .5*DISTS)**6 * (35./12.*DISTS**2 + 3.*DISTS + 1.)
       ELSE
        DIST(I)=0.0
       END IF
      END DO
#elif ikernel == 2
!     Wendland C6 kernel
      H=W/2.0
      DO I=1,N2
       DISTS=DIST(I)/H
       IF (DISTS.LE.2.0) THEN
        DIST(I)=(1. - .5*DISTS)**8 * (4.*DISTS**3 + 25./4.*DISTS**2 +
     &                                4.*DISTS + 1.)
       ELSE
        DIST(I)=0.0
       END IF
      END DO
#elif ikernel == 3
!     Quartic spline kernel (M5)
      H=W/2.5
      DO I=1,N2
       DISTS=DIST(I)/H
       IF (DISTS.GT.2.5) THEN 
        DIST(I)=0.0 
       ELSE
        DIST(I)=(2.5-DISTS)**4
        IF (DISTS.LT.1.5) DIST(I)=DIST(I)-5.0*(1.5-DISTS)**4
        IF (DISTS.LT.0.5) DIST(I)=DIST(I)+10.0*(0.5-DISTS)**4 
       END IF 
      END DO 
#elif ikernel == 4
!     Quintic spline kernel (M6)
      H=W/3.0
      DO I=1,N2
       DISTS=DIST(I)/H
       IF (DISTS.GT.3.0) THEN 
        DIST(I)=0.0 
       ELSE
        DIST(I)=(3.0-DISTS)**5
        IF (DISTS.LT.2.0) DIST(I)=DIST(I)-6.0*(2.0-DISTS)**5
        IF (DISTS.LT.1.0) DIST(I)=DIST(I)+15.0*(1.0-DISTS)**5 
       END IF 
      END DO
#endif
#endif

      RETURN
      END


************************************************************************
      SUBROUTINE VEINSGRID_REDUCED(IR,NPATCH,PARE,
     &      PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &      PATCHRY,PATCHRZ,CR01,CONTA11,LOW1,LOW2)
************************************************************************
*     small fraction of patches with rare geometry were not detected
*     when overlapping

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NPALEV2

*     U11(PATCHNX,PATCHNY,PATCHNZ,NLEVEL,NPALEV)
*     PATCHNX,PATCHNY,PATCHNZ patches dimensions
*     IPATCH number of patches per level
*     NLEVELS total number of levels

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV)
      REAL PATCHRY(NPALEV)
      REAL PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK,I2,J2
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2
      INTEGER NV,A2,B2,C2,K

      INTEGER LOW1,LOW2
      INTEGER SOLAP_PATCH(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CR01(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)
      INTEGER CONTA11(NAMRX,NAMRY,NAMRZ,LOW1:LOW2)

      REAL A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL RZ1,RZZ1,RZ2,RZZ2,ORXX1,ORYY1,ORZZ1

      REAL DXPA,DYPA,DZPA
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER,ALLOCATABLE::VECINO(:,:)
      INTEGER,ALLOCATABLE::NVECI(:)

      INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
      REAL RXFIX,RYFIX,RZFIX

      REAL OVERLAP
*
       NPALEV2=MAX(100,INT(NPALEV/5))

       ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       ALLOCATE(NVECI(NPATCH(IR)))

       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

!$OMP   PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+        DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,
!$OMP+        SOLAP_PATCH,LOW1,LOW2,RXFIX,RYFIX,RZFIX),
!$OMP+  PRIVATE(I,N1,N2,N3,NV,J,NN1,NN2,NN3,
!$OMP+          RX1,RY1,RZ1,RXX1,RYY1,RZZ1,I2,
!$OMP+          IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4),
!$OMP+  DEFAULT(NONE)
       DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         SOLAP_PATCH(:,:,:,I)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA

         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA

          IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

          IG4=IG3 + NN1 - 1
          JG4=JG3 + NN2 - 1
          KG4=KG3 + NN3 - 1

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
           NV=NV+1
           VECINO(NV,I2)=J
          END IF

          END IF

         END DO
         NVECI(I2)=NV
       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) THEN
         WRITE(*,*) 'ERROR: gvecino ST second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))
         STOP
       END IF


       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
         J=VECINO(K,I2)
         J2=J-LOW1+1

         LL1=PATCHX(J)
         LL2=PATCHY(J)
         LL3=PATCHZ(J)

         NN1=PATCHNX(J)
         NN2=PATCHNY(J)
         NN3=PATCHNZ(J)

         RXX1=PATCHRX(J)-0.5*DXPA
         RYY1=PATCHRY(J)-0.5*DYPA
         RZZ1=PATCHRZ(J)-0.5*DZPA
         RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
         RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
         RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*        X
         IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=1
            CORNXX2=NN1
         END IF
         IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=N1
            CORNXX1=1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
         END IF
         IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
            CORNX1=1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=NN1
         END IF
         IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=1
            CORNX2=N1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
         END IF

*        Y
         IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=1
            CORNYY2=NN2
         END IF
         IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=N2
            CORNYY1=1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
         END IF
         IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
            CORNY1=1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=NN2
         END IF
         IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=1
            CORNY2=N2
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
         END IF

*        Z
         IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=N3
            CORNZZ1=1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
         END IF
         IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
            CORNZ1=1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=1
            CORNZ2=N3
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
         END IF

*        overlap is the fraction of volumen overlapping
*        the idea is to fix a thershold of overlapping that below
*        this therhold for instance 10%, the cells are not marked as overlap
         OVERLAP=(CORNZZ2-CORNZZ1+1)*(CORNYY2-CORNYY1+1)*
     &           (CORNXX2-CORNXX1+1)
         OVERLAP=ABS(OVERLAP)/PATCHNX(J)/PATCHNY(J)/PATCHNZ(J)

**       celdas madre del nivel inferior
         IF (OVERLAP.GT.0.1) THEN
           DO KK=CORNZZ1,CORNZZ2
           DO JJ=CORNYY1,CORNYY2
           DO II=CORNXX1,CORNXX2
            IX=II-CORNXX1+CORNX1
            JY=JJ-CORNYY1+CORNY1
            KZ=KK-CORNZZ1+CORNZ1
            IF (SOLAP_PATCH(IX,JY,KZ,I).EQ.0) THEN
*            the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i
             SOLAP_PATCH(II,JJ,KK,J)=1
             CR01(II,JJ,KK,J)=0
             CONTA11(II,JJ,KK,J)=0
            END IF
           END DO
           END DO
           END DO
          END IF
       END DO
       END DO

      DEALLOCATE(VECINO)
      DEALLOCATE(NVECI)

      RETURN
      END


************************************************************************
      SUBROUTINE VEINSGRID_CR0AMR(IR,NPATCH,PARE,
     &      PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &      PATCHRY,PATCHRZ)
************************************************************************
*     small fraction of patches with rare geometry were not detected
*     when overlapping

      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

      INTEGER NPALEV2

*     U11(PATCHNX,PATCHNY,PATCHNZ,NLEVEL,NPALEV)
*     PATCHNX,PATCHNY,PATCHNZ patches dimensions
*     IPATCH number of patches per level
*     NLEVELS total number of levels

      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV)
      INTEGER PATCHNY(NPALEV)
      INTEGER PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV)
      INTEGER PATCHY(NPALEV)
      INTEGER PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV)
      REAL PATCHRY(NPALEV)
      REAL PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER CR1,CR2,CR3,CR4,CR5,CR6
      INTEGER IR,I,J,IX,JY,KZ,II,JJ,KK,I2,J2
      INTEGER N1,N2,N3,L1,L2,L3,NN1,NN2,NN3,LL1,LL2,LL3
      INTEGER KK2,JJ2,II2,KZ2,JY2,IX2
      INTEGER NV,A2,B2,C2,K

      INTEGER LOW1,LOW2
      INTEGER SOLAP_PATCH(NAMRX,NAMRY,NAMRZ,NPATCH(IR))

      integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      common /cr0/ cr0amr, cr0amr1

      REAL A1,B1,C1,RIV1,RIV2,RIV3
      INTEGER CONTROL
      INTEGER CORNX1,CORNXX1,CORNX2,CORNXX2
      INTEGER CORNY1,CORNYY1,CORNY2,CORNYY2
      INTEGER CORNZ1,CORNZZ1,CORNZ2,CORNZZ2
      REAL RX1,RXX1,RX2,RXX2,RY1,RYY1,RY2,RYY2
      REAL RZ1,RZZ1,RZ2,RZZ2,ORXX1,ORYY1,ORZZ1

      REAL DXPA,DYPA,DZPA
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &        RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &        RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/   RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      INTEGER,ALLOCATABLE::VECINO(:,:)
      INTEGER,ALLOCATABLE::NVECI(:)

      INTEGER IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4
      REAL*4 RXFIX,RYFIX,RZFIX

      REAL*4 OVERLAP
*
       NPALEV2=MAX(100,INT(NPALEV/5))

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

       ALLOCATE(VECINO(NPALEV2,NPATCH(IR)))
       ALLOCATE(NVECI(NPATCH(IR)))

       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

*      built auxiliar grid for comparison
       RXFIX=RADX(1) - DX*0.5 + 0.5*DXPA
       RYFIX=RADY(1) - DY*0.5 + 0.5*DYPA
       RZFIX=RADZ(1) - DZ*0.5 + 0.5*DZPA

!$OMP   PARALLEL DO SHARED(IR,NPATCH,PARE,PATCHX,PATCHY,PATCHZ,
!$OMP+        PATCHNX,PATCHNY,PATCHNZ,VECINO,NVECI,
!$OMP+        DXPA,DYPA,DZPA,PATCHRX,PATCHRY,PATCHRZ,
!$OMP+        SOLAP_PATCH,LOW1,LOW2,RXFIX,RYFIX,RZFIX),
!$OMP+  PRIVATE(I,N1,N2,N3,NV,J,NN1,NN2,NN3,
!$OMP+          RX1,RY1,RZ1,RXX1,RYY1,RZZ1,I2,
!$OMP+          IG1,IG2,JG1,JG2,KG1,KG2,IG3,JG3,KG3,IG4,JG4,KG4),
!$OMP+  DEFAULT(NONE)
       DO I=LOW1,LOW2

         I2=I-LOW1+1

         NVECI(I2)=0
         VECINO(:,I2)=0

         SOLAP_PATCH(:,:,:,I2)=0

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         NV=0

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA

         IG1=INT(((RX1-RXFIX)/DXPA)+0.5) + 1
         JG1=INT(((RY1-RYFIX)/DYPA)+0.5) + 1
         KG1=INT(((RZ1-RZFIX)/DZPA)+0.5) + 1

         IG2=IG1 + N1 - 1
         JG2=JG1 + N2 - 1
         KG2=KG1 + N3 - 1

         DO J=LOW1,LOW2
          IF (J.NE.I) THEN

          NN1=PATCHNX(J)
          NN2=PATCHNY(J)
          NN3=PATCHNZ(J)

          RXX1=PATCHRX(J)-0.5*DXPA
          RYY1=PATCHRY(J)-0.5*DYPA
          RZZ1=PATCHRZ(J)-0.5*DZPA

          IG3=INT(((RXX1-RXFIX)/DXPA)+0.5) + 1
          JG3=INT(((RYY1-RYFIX)/DYPA)+0.5) + 1
          KG3=INT(((RZZ1-RZFIX)/DZPA)+0.5) + 1

          IG4=IG3 + NN1 - 1
          JG4=JG3 + NN2 - 1
          KG4=KG3 + NN3 - 1

          IF (IG1.LE.IG4.AND.IG3.LE.IG2.AND.
     &        JG1.LE.JG4.AND.JG3.LE.JG2.AND.
     &        KG1.LE.KG4.AND.KG3.LE.KG2) THEN
           NV=NV+1
           VECINO(NV,I2)=J
          END IF

          END IF

         END DO
         NVECI(I2)=NV
       END DO


       IF (MAXVAL(NVECI(1:NPATCH(IR))).GT.NPALEV2) THEN
         WRITE(*,*) 'ERROR: gvecino ST second dimension too large',
     &     MAXVAL(NVECI(1:NPATCH(IR)))
         STOP
       END IF


       DO I=LOW1,LOW2

         L1=PATCHX(I)
         L2=PATCHY(I)
         L3=PATCHZ(I)

         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)

         RX1=PATCHRX(I)-0.5*DXPA
         RY1=PATCHRY(I)-0.5*DYPA
         RZ1=PATCHRZ(I)-0.5*DZPA
         RX2=PATCHRX(I)-0.5*DXPA+(N1-1)*DXPA
         RY2=PATCHRY(I)-0.5*DYPA+(N2-1)*DYPA
         RZ2=PATCHRZ(I)-0.5*DZPA+(N3-1)*DZPA

         I2=I-LOW1+1

         DO K=1,NVECI(I2)
         J=VECINO(K,I2)
         J2=J-LOW1+1

         LL1=PATCHX(J)
         LL2=PATCHY(J)
         LL3=PATCHZ(J)

         NN1=PATCHNX(J)
         NN2=PATCHNY(J)
         NN3=PATCHNZ(J)

         RXX1=PATCHRX(J)-0.5*DXPA
         RYY1=PATCHRY(J)-0.5*DYPA
         RZZ1=PATCHRZ(J)-0.5*DZPA
         RXX2=PATCHRX(J)-0.5*DXPA+(NN1-1)*DXPA
         RYY2=PATCHRY(J)-0.5*DYPA+(NN2-1)*DYPA
         RZZ2=PATCHRZ(J)-0.5*DZPA+(NN3-1)*DZPA

*        X
         IF (RXX1.GE.RX1.AND.RXX2.LE.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=1
            CORNXX2=NN1
         END IF
         IF (RXX1.GE.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=INT(((RXX1-RX1)/DXPA)+0.5) + 1
            CORNX2=N1
            CORNXX1=1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) +1
         END IF
         IF (RXX2.LE.RX2.AND.RXX1.LT.RX1) THEN
            CORNX1=1
            CORNX2=INT(((RXX2-RX1)/DXPA)+0.5) + 1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=NN1
         END IF
         IF (RXX1.LT.RX1.AND.RXX2.GT.RX2) THEN
            CORNX1=1
            CORNX2=N1
            CORNXX1=INT(((RX1-RXX1)/DXPA)+0.5) + 1
            CORNXX2=INT(((RX2-RXX1)/DXPA)+0.5) + 1
         END IF

*        Y
         IF (RYY1.GE.RY1.AND.RYY2.LE.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=1
            CORNYY2=NN2
         END IF
         IF (RYY1.GE.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=INT(((RYY1-RY1)/DYPA)+0.5) + 1
            CORNY2=N2
            CORNYY1=1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) +1
         END IF
         IF (RYY2.LE.RY2.AND.RYY1.LT.RY1) THEN
            CORNY1=1
            CORNY2=INT(((RYY2-RY1)/DYPA)+0.5) + 1
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=NN2
         END IF
         IF (RYY1.LT.RY1.AND.RYY2.GT.RY2) THEN
            CORNY1=1
            CORNY2=N2
            CORNYY1=INT(((RY1-RYY1)/DYPA)+0.5) + 1
            CORNYY2=INT(((RY2-RYY1)/DYPA)+0.5) + 1
         END IF

*        Z
         IF (RZZ1.GE.RZ1.AND.RZZ2.LE.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.GE.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=INT(((RZZ1-RZ1)/DZPA)+0.5) + 1
            CORNZ2=N3
            CORNZZ1=1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) +1
         END IF
         IF (RZZ2.LE.RZ2.AND.RZZ1.LT.RZ1) THEN
            CORNZ1=1
            CORNZ2=INT(((RZZ2-RZ1)/DZPA)+0.5) + 1
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=NN3
         END IF
         IF (RZZ1.LT.RZ1.AND.RZZ2.GT.RZ2) THEN
            CORNZ1=1
            CORNZ2=N3
            CORNZZ1=INT(((RZ1-RZZ1)/DZPA)+0.5) + 1
            CORNZZ2=INT(((RZ2-RZZ1)/DZPA)+0.5) + 1
         END IF

**       celdas madre del nivel inferior
          DO KK=CORNZZ1,CORNZZ2
          DO JJ=CORNYY1,CORNYY2
          DO II=CORNXX1,CORNXX2
           IX=II-CORNXX1+CORNX1
           JY=JJ-CORNYY1+CORNY1
           KZ=KK-CORNZZ1+CORNZ1
           IF (SOLAP_PATCH(IX,JY,KZ,I2).EQ.0) THEN
*            the cell ii,jj,kk,ir,j overlaps ix,jy,kz,ir,i
            IF (CR0AMR1(IX,JY,KZ,I).EQ.0) THEN
             CR0AMR1(II,JJ,KK,J)=0
            ELSE IF (CR0AMR1(II,JJ,KK,J).EQ.0) THEN
             CR0AMR1(IX,JY,KZ,I)=0
            END IF
           END IF
          END DO
          END DO
          END DO
       END DO
       END DO

      DEALLOCATE(VECINO)
      DEALLOCATE(NVECI)

      RETURN
      END

***************************************************************************
      SUBROUTINE indexx(n,arr,indx)
***************************************************************************
* From Press,Teukoslky,Vetterling & Flannery, Numerical Recipes in Fortran90,
* Cambridge University Press
***************************************************************************

      INTEGER n,indx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
         write(*,*) 'NSTACK too small in indexx'
         stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


************************************************************************
      SUBROUTINE INTERPOLATE_VELOCITIES(NX,NY,NZ,NL,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,FLAG_FILTER,KNEIGHBOURS,
     &            VISC0,VISC1,FLAG_MACHFIELD,FLAG_MASS)
************************************************************************
*     Compute the velocity field on the grid
************************************************************************

*     From coretran ****************************************
      use variableKind, only: i32, r64
      use m_allocate, only: allocate
      use m_deallocate, only: deallocate
      use m_KdTree, only: KdTree, KdTreeSearch
      use dargdynamicarray_class, only: dArgDynamicArray
************************************************************

      use particle_data
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     function parameters
      INTEGER NX,NY,NZ,NL
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS),PARE(NPALEV)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      REAL LADO0
      INTEGER FLAG_FILTER,FLAG_MACHFIELD,FLAG_MASS
*     COMMON VARIABLES
      REAL DX,DY,DZ
      COMMON /ESPACIADO/ DX,DY,DZ

      REAL  RADX(0:NMAX+1),RADMX(0:NMAX+1),
     &      RADY(0:NMAY+1),RADMY(0:NMAY+1),
     &      RADZ(0:NMAZ+1),RADMZ(0:NMAZ+1)
      COMMON /GRID/  RADX,RADMX,RADY,RADMY,RADZ,RADMZ

      REAL RX(-2:NAMRX+3,NPALEV)
      REAL RY(-2:NAMRX+3,NPALEV)
      REAL RZ(-2:NAMRX+3,NPALEV)
      REAL RMX(-2:NAMRX+3,NPALEV)
      REAL RMY(-2:NAMRX+3,NPALEV)
      REAL RMZ(-2:NAMRX+3,NPALEV)
      COMMON /MINIGRIDS/ RX,RY,RZ,RMX,RMY,RMZ

      INTEGER cr0amr(1:NMAX,1:NMAY,1:NMAZ)
      INTEGER cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      COMMON /cr0/ cr0amr, cr0amr1
      INTEGER solap(NAMRX,NAMRY,NAMRZ,NPALEV)

      REAL U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

#ifdef use_filter
#if use_filter == 1
             REAL VISC0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
             REAL VISC1(NAMRX,NAMRY,NAMRZ,NPALEV)
#else 
             ! Dummy variables
             REAL VISC0, VISC1 
#endif
#endif

      !real u1(1:NMAX,1:NMAY,1:NMAZ)
      !real u11(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      !common /dens/ u1,u11

      REAL OMEGA0,ACHE,FDM,RHOB0
      COMMON /COSMO/ OMEGA0,ACHE,FDM

      REAL L0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL L1(NAMRX,NAMRY,NAMRZ,NPALEV)
#ifdef weight_filter
#if weight_filter == 1
      COMMON /DENSI/ L0,L1
#endif 
#endif

      INTEGER IX,JY,KZ,IR,I,J,K,IPATCH,LOW1,LOW2,CONTA,KNEIGHBOURS
      INTEGER N1,N2,N3,II,JJ,KK,JPATCH,I1,I2,J1,J2,K1,K2,STEP
      INTEGER NPART_TOT,II1,II2,JJ1,JJ2,KK1,KK2,IIP1,JJP1,KKP1
      REAL DXPA,DYPA,DZPA,BASX,BASY,BASZ,BAS,RBAS,BASXX,BASYY,BASZZ
      REAL MEDIOLADO0,PI,MINX,MAXX,MINY,MAXY,MINZ,MAXZ,FRAC_INT,H_KERN
      REAL*8 BAS8,BAS8X,BAS8Y,BAS8Z,BAS8M,BASMASS
      REAL,ALLOCATABLE::DIST(:)
      INTEGER,ALLOCATABLE::NEIGH(:),LB(:)
      REAL FUIN,U(2,2,2),UW(2,2,2)
      logical do_cell

      !INTEGER,ALLOCATABLE::SCRINT(:,:,:)
      integer*4 t1,t2,trate,tmax

      REAL(R64), ALLOCATABLE :: XTREE(:), YTREE(:), ZTREE(:)
      TYPE(KDTREE) :: TREE 
      TYPE (KDTREESEARCH) :: SEARCH
      TYPE (DARGDYNAMICARRAY) :: DA

      integer omp_get_thread_num
      WRITE(*,*) 'Each cell-center value is computed using at least',
     &            KNEIGHBOURS,'particles'

#ifdef ikernel 
#if ikernel == 0
       WRITE(*,*) 'Kernel: cubic spline (M4)'
#elif ikernel == 1
       WRITE(*,*) 'Kernel: Wendland C4'
#elif ikernel == 2
       WRITE(*,*) 'Kernel: Wendland C6'
#elif ikernel == 3
       WRITE(*,*) 'Kernel: quartic spline (M5)'
#elif ikernel == 4
       WRITE(*,*) 'Kernel: quintic spline (M6)'
#endif
#endif

      MEDIOLADO0=0.5*LADO0
      PI=ACOS(-1.0)
C      RHOB0=MAXVAL(MASAP)/FDM/DX**3
C      ! to get overdensity from mass in a sphere (multiply by mass in
C      ! code units, divide by radius squared)
C      CONSTA_DENS=1/(4*PI/3)/RHOB0
C      WRITE(*,*) 'Bkg density (code units), mass to overdensity const',
C     &            RHOB0,CONSTA_DENS

      CALL VEINSGRID_ALL_L(NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                     PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &                     SOLAP)

      WRITE(*,*) 'Building KDTree...'

      NPART_TOT=SUM(NPART)
      CALL ALLOCATE(XTREE, NPART_TOT)
      CALL ALLOCATE(YTREE, NPART_TOT)
      CALL ALLOCATE(ZTREE, NPART_TOT)

!$OMP PARALLEL DO SHARED(NPART_TOT,RXPA,RYPA,RZPA,XTREE,YTREE,ZTREE),
!$OMP+            PRIVATE(I),
!$OMP+            DEFAULT(NONE)
      DO I=1,NPART_TOT 
       XTREE(I)=RXPA(I)
       YTREE(I)=RYPA(I)
       ZTREE(I)=RZPA(I)
      END DO

      call system_clock(t1,trate,tmax)
      TREE = KDTREE(XTREE, YTREE, ZTREE)
      call system_clock(t2,trate,tmax)
      WRITE(*,*) 'KDTree built in ',float(t2-t1)/1.e3,' seconds'

*     Now, go (clean) cell to (clean) cell, finding the neighbouring
*     particles and interpolating the velocity field onto the cell.
*     Around each cell, we consider all the particles within a sphere 
*     of radius MAX(DXPA, LNEIGH), where LNEIGH is the distance to the
*     KNEIGHBOURS-th nearest particle. 

*     Base grid: to save computational cost, we consider three regions:
*     1) outside the domain where there are particles, we go 4 cells 
*     by 4 cells, and then will interpolate. (outside [I1,I2])
*     2) inside this domain, but outside the region defined by the 
*     5-95 percentiles around each position, we go 2 cells by 2 cells.
*     (outside [II1,II2], but inside [I1,I2])
*     3) inside this region, we go cell by cell. (inside [II1,II2])

      ALLOCATE(NEIGH(NPART_TOT))
      FRAC_INT=0.01
      
*     X      
      CALL ARGSORT(NPART_TOT,RXPA(1:NPART_TOT),NEIGH)
      MINX=RXPA(NEIGH(1))
      MAXX=RXPA(NEIGH(NPART_TOT))
      I1=INT((MINX+MEDIOLADO0)/DX)+1
      I2=INT((MAXX+MEDIOLADO0)/DX)+1
      MINX=RXPA(NEIGH(INT(FRAC_INT*NPART_TOT)))
      MAXX=RXPA(NEIGH(INT((1.-FRAC_INT)*NPART_TOT)))
      II1=INT((MINX+MEDIOLADO0)/DX)+1
      II2=INT((MAXX+MEDIOLADO0)/DX)+1
      if (mod(i1,4).eq.3.and.i1.ge.2) i1 = i1-1
      if (mod(i2,4).eq.3.and.i2.le.nx-1) i2 = i2+1 

*     Y
      CALL ARGSORT(NPART_TOT,RYPA(1:NPART_TOT),NEIGH)
      MINY=RYPA(NEIGH(1))
      MAXY=RYPA(NEIGH(NPART_TOT))
      J1=INT((MINY+MEDIOLADO0)/DY)+1
      J2=INT((MAXY+MEDIOLADO0)/DY)+1
      MINY=RYPA(NEIGH(INT(FRAC_INT*NPART_TOT)))
      MAXY=RYPA(NEIGH(INT((1.-FRAC_INT)*NPART_TOT)))
      JJ1=INT((MINY+MEDIOLADO0)/DY)+1
      JJ2=INT((MAXY+MEDIOLADO0)/DY)+1
      if (mod(j1,4).eq.3.and.j1.ge.2) j1 = j1-1
      if (mod(j2,4).eq.3.and.j2.le.ny-1) j2 = j2+1

*     Z
      CALL ARGSORT(NPART_TOT,RZPA(1:NPART_TOT),NEIGH)
      MINZ=RZPA(NEIGH(1))
      MAXZ=RZPA(NEIGH(NPART_TOT))
      K1=INT((MINZ+MEDIOLADO0)/DZ)+1
      K2=INT((MAXZ+MEDIOLADO0)/DZ)+1
      MINZ=RZPA(NEIGH(INT(FRAC_INT*NPART_TOT)))
      MAXZ=RZPA(NEIGH(INT((1.-FRAC_INT)*NPART_TOT)))
      KK1=INT((MINZ+MEDIOLADO0)/DZ)+1
      KK2=INT((MAXZ+MEDIOLADO0)/DZ)+1
      if (mod(k1,4).eq.3.and.k1.ge.2) k1 = k1-1
      if (mod(k2,4).eq.3.and.k2.le.nz-1) k2 = k2+1

      DEALLOCATE(NEIGH)

c      WRITE(*,*) I1,II1,II2,I2
c      WRITE(*,*) J1,JJ1,JJ2,J2
c      WRITE(*,*) K1,KK1,KK2,K2

************************
*     BASE GRID
************************

      WRITE(*,*) 'Base grid...',0

!$OMP PARALLEL DO SHARED(NX,NY,NZ,L0),
!$OMP+            PRIVATE(IX,JY,KZ),
!$OMP+            DEFAULT(NONE)
      DO KZ=0,NZ+1
      DO JY=0,NY+1
      DO IX=0,NX+1
       L0(IX,JY,KZ)=0.
      END DO 
      END DO 
      END DO

*     1. outside the domain where there are particles, we go 4 cells
*        by 4 cells, and then will interpolate. (outside [I1,I2]).
      STEP=4

!$OMP PARALLEL DO SHARED(NX,NY,NZ,STEP,I1,I2,J1,J2,K1,K2,RADX,RADY,
!$OMP+                   RADZ,U2DM,U3DM,U4DM,L0,U2,U3,U4,MASAP,VOL,
!$OMP+                   KNEIGHBOURS,DX,VISC0,ABVC,TREE,XTREE,
!$OMP+                   YTREE,ZTREE,PI,FLAG_MASS),
!$OMP+            PRIVATE(IX,JY,KZ,DIST,NEIGH,CONTA,H_KERN,BAS8,
!$OMP+                    BAS8X,BAS8Y,BAS8Z,BAS8M,I,BASMASS,DA,SEARCH),
!$OMP+            SCHEDULE(DYNAMIC)!, DEFAULT(NONE)
      DO KZ=1,NZ+1,STEP
      DO JY=1,NY+1,STEP
      DO IX=1,NX+1,STEP
       IF (IX.GT.INT(I1/STEP+1)*STEP.AND.
     &     IX.LT.INT(I2/STEP)*STEP.AND.
     &     JY.GT.INT(J1/STEP+1)*STEP.AND.
     &     JY.LT.INT(J2/STEP)*STEP.AND.
     &     KZ.GT.INT(K1/STEP+1)*STEP.AND.
     &     KZ.LT.INT(K2/STEP)*STEP) CYCLE

       ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
       DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &       XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &       ZQUERY=DBLE(RADZ(KZ)), K=KNEIGHBOURS)
       DIST = DA%V%VALUES 
       NEIGH = DA%I%VALUES

       IF (DIST(KNEIGHBOURS).GT.DX) THEN
        CONTA=KNEIGHBOURS 
       ELSE 
        DEALLOCATE(DIST,NEIGH)
        DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &        XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &        ZQUERY=DBLE(RADZ(KZ)), RADIUS=DBLE(DX))
        CONTA = DA%SIZE() 
        ALLOCATE(DIST(CONTA), NEIGH(CONTA))
        DIST = DA%V%VALUES 
        NEIGH = DA%I%VALUES
       END IF

       !WRITE(*,*) IX,JY,KZ,CONTA,DIST(CONTA)
       H_KERN=DIST(CONTA)
       L0(IX,JY,KZ)=H_KERN

       CALL KERNEL_FUNC(CONTA,CONTA,H_KERN,DIST)
#ifdef weight_scheme
#if weight_scheme == 1
       DO I=1,CONTA
        DIST(I)=DIST(I)*MASAP(NEIGH(I))
       END DO
#elif weight_scheme == 2
       DO I=1,CONTA 
        DIST(I)=DIST(I)*VOL(NEIGH(I))
       END DO
#endif
#endif

       BAS8=0.D0
       BAS8X=0.D0
       BAS8Y=0.D0
       BAS8Z=0.D0
       BAS8M=0.D0
       BASMASS=0.D0
       DO I=1,CONTA 
        BAS8=BAS8+DIST(I)
        BAS8X=BAS8X+DIST(I)*U2DM(NEIGH(I))
        BAS8Y=BAS8Y+DIST(I)*U3DM(NEIGH(I))
        BAS8Z=BAS8Z+DIST(I)*U4DM(NEIGH(I))
#ifdef use_filter
#if use_filter == 1
        BAS8M=BAS8M+DIST(I)*ABVC(NEIGH(I))
#endif
#endif
        !BAS8M=MAX(BAS8M,ABVC(NEIGH(I)))
        BASMASS=BASMASS+MASAP(NEIGH(I))
       END DO
       U2(IX,JY,KZ)=BAS8X/BAS8
       U3(IX,JY,KZ)=BAS8Y/BAS8
       U4(IX,JY,KZ)=BAS8Z/BAS8
#ifdef use_filter
#if use_filter == 1
       VISC0(IX,JY,KZ)=BAS8M/BAS8
#endif
#endif
       !VISC0(IX,JY,KZ)=BAS8M


       IF (FLAG_MASS.EQ.1) L0(IX,JY,KZ)=BASMASS/(4*PI/3)/H_KERN**3

       DEALLOCATE(DIST, NEIGH)
      END DO 
      END DO 
      END DO

      ! Fill the blanks by interpolation

!$OMP PARALLEL DO SHARED(NX,NY,NZ,STEP,I1,I2,J1,J2,K1,K2,RADX,RADY,
!$OMP+                   RADZ,U2,U3,U4,L0,DX,DY,DZ,VISC0),
!$OMP+            PRIVATE(IX,JY,KZ,II,JJ,KK,IIP1,JJP1,KKP1,BASX,BASY,
!$OMP+                    BASZ),
!$OMP+            DEFAULT(NONE), SCHEDULE(DYNAMIC)
      DO KZ=1,NZ
       KK=INT((KZ-1)/STEP)*STEP+1
       KKP1=KK+STEP!MOD(KK+STEP,NZ)
       BASZ=(RADZ(KZ)-RADZ(KK))/(STEP*DZ)
      DO JY=1,NY
       JJ=INT((JY-1)/STEP)*STEP+1
       JJP1=JJ+STEP!MOD(JJ+STEP,NY)
       BASY=(RADY(JY)-RADY(JJ))/(STEP*DY)
      DO IX=1,NX
       IF (IX.GE.I1.AND.IX.LE.I2.AND.
     &     JY.GE.J1.AND.JY.LE.J2.AND.
     &     KZ.GE.K1.AND.KZ.LE.K2) CYCLE
       II=INT((IX-1)/STEP)*STEP+1
       IIP1=II+STEP!MOD(II+STEP,NX)
       BASX=(RADX(IX)-RADX(II))/(STEP*DX)

       IF (IX.EQ.II.AND.JY.EQ.JJ.AND.KZ.EQ.KK) CYCLE 

       U2(IX,JY,KZ)=U2(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U2(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U2(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U2(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U2(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U2(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U2(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U2(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ

       U3(IX,JY,KZ)=U3(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U3(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U3(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U3(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U3(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U3(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U3(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U3(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ

       U4(IX,JY,KZ)=U4(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U4(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U4(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U4(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U4(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U4(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U4(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U4(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ     

       L0(IX,JY,KZ)=L0(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              L0(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              L0(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              L0(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              L0(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              L0(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              L0(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              L0(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ  
     
#ifdef use_filter
#if use_filter == 1
       VISC0(IX,JY,KZ)=VISC0(II,JJ,KK)  *(1-BASX)*(1-BASY)*(1-BASZ) +
     &                 VISC0(IIP1,JJ,KK)*  BASX  *(1-BASY)*(1-BASZ) +
     &                 VISC0(II,JJP1,KK)*(1-BASX)* BASY *(1-BASZ) +
     &                 VISC0(IIP1,JJP1,KK)* BASX * BASY *(1-BASZ) +
     &                 VISC0(II,JJ,KKP1)*(1-BASX)*(1-BASY)* BASZ  +
     &                 VISC0(IIP1,JJ,KKP1)* BASX *(1-BASY)* BASZ   +
     &                 VISC0(II,JJP1,KKP1)*(1-BASX)* BASY * BASZ   +
     &                 VISC0(IIP1,JJP1,KKP1)* BASX * BASY * BASZ  
#endif
#endif

      END DO 
      END DO 
      END DO

      ! 2. inside this domain, but outside the region defined by the
      ! 5-95 percentiles around each position, we go 2 cells by 2 cells.
      ! (outside [II1,II2], but inside [I1,I2])
      STEP=2

!$OMP PARALLEL DO SHARED(NX,NY,NZ,STEP,I1,I2,J1,J2,K1,K2,II1,II2,JJ1,
!$OMP+                   JJ2,KK1,KK2,RADX,RADY,RADZ,U2DM,U3DM,
!$OMP+                   U4DM,L0,U2,U3,U4,KNEIGHBOURS,DX,VISC0,ABVC,
!$OMP+                   TREE,XTREE,YTREE,ZTREE,FLAG_MASS,
!$OMP+                   MASAP,VOL,PI),
!$OMP+            PRIVATE(IX,JY,KZ,DIST,NEIGH,CONTA,H_KERN,BAS8,
!$OMP+                    BAS8X,BAS8Y,BAS8Z,BAS8M,I,BASMASS,DA,SEARCH),
!$OMP+            SCHEDULE(DYNAMIC)!, DEFAULT(NONE)
      DO KZ=K1,K2+1,STEP
      DO JY=J1,J2+1,STEP
      DO IX=I1,I2+1,STEP
       IF (IX.GT.INT(II1/STEP+1)*STEP.AND.
     &     IX.LT.INT(II2/STEP)*STEP.AND.
     &     JY.GT.INT(JJ1/STEP+1)*STEP.AND.
     &     JY.LT.INT(JJ2/STEP)*STEP.AND.
     &     KZ.GT.INT(KK1/STEP+1)*STEP.AND.
     &     KZ.LT.INT(KK2/STEP)*STEP) CYCLE

       ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
       DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &       XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &       ZQUERY=DBLE(RADZ(KZ)), K=KNEIGHBOURS)
       DIST = DA%V%VALUES 
       NEIGH = DA%I%VALUES

       IF (DIST(KNEIGHBOURS).GT.DX) THEN
        CONTA=KNEIGHBOURS 
       ELSE 
        DEALLOCATE(DIST,NEIGH)
        DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &        XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &        ZQUERY=DBLE(RADZ(KZ)), RADIUS=DBLE(DX))
        CONTA = DA%SIZE() 
        ALLOCATE(DIST(CONTA), NEIGH(CONTA))
        DIST = DA%V%VALUES 
        NEIGH = DA%I%VALUES
       END IF

       !WRITE(*,*) IX,JY,KZ,CONTA,DIST(CONTA)
       H_KERN=DIST(CONTA)
       L0(IX,JY,KZ)=H_KERN

       CALL KERNEL_FUNC(CONTA,CONTA,H_KERN,DIST)
#ifdef weight_scheme
#if weight_scheme == 1
      DO I=1,CONTA
       DIST(I)=DIST(I)*MASAP(NEIGH(I))
      END DO
#elif weight_scheme == 2
      DO I=1,CONTA 
       DIST(I)=DIST(I)*VOL(NEIGH(I))
      END DO
#endif
#endif

       BAS8=0.D0
       BAS8X=0.D0
       BAS8Y=0.D0
       BAS8Z=0.D0
       BAS8M=0.D0
       BASMASS=0.D0
       DO I=1,CONTA 
        BAS8=BAS8+DIST(I)
        BAS8X=BAS8X+DIST(I)*U2DM(NEIGH(I))
        BAS8Y=BAS8Y+DIST(I)*U3DM(NEIGH(I))
        BAS8Z=BAS8Z+DIST(I)*U4DM(NEIGH(I))
#ifdef use_filter
#if use_filter == 1
        BAS8M=BAS8M+DIST(I)*ABVC(NEIGH(I))
#endif
#endif
        !BAS8M=MAX(BAS8M,ABVC(NEIGH(I)))
        BASMASS=BASMASS+MASAP(NEIGH(I))
       END DO
       
       U2(IX,JY,KZ)=BAS8X/BAS8
       U3(IX,JY,KZ)=BAS8Y/BAS8
       U4(IX,JY,KZ)=BAS8Z/BAS8
#ifdef use_filter
#if use_filter == 1
       VISC0(IX,JY,KZ)=BAS8M/BAS8
#endif
#endif
       !VISC0(IX,JY,KZ)=BAS8M

       IF (FLAG_MASS.EQ.1) L0(IX,JY,KZ)=BASMASS/(4*PI/3)/H_KERN**3

       DEALLOCATE(DIST, NEIGH)
      END DO 
      END DO 
      END DO

      ! Fill the blanks by interpolation

!$OMP PARALLEL DO SHARED(NX,NY,NZ,STEP,II1,II2,JJ1,JJ2,KK1,KK2,RADX,
!$OMP+                   RADY,RADZ,U2,U3,U4,L0,DX,DY,DZ,I1,I2,J1,J2,K1,
!$OMP+                   K2,VISC0),
!$OMP+            PRIVATE(IX,JY,KZ,II,JJ,KK,IIP1,JJP1,KKP1,BASX,BASY,
!$OMP+                    BASZ),
!$OMP+            DEFAULT(NONE), SCHEDULE(DYNAMIC)
      DO KZ=K1,K2
       KK=K1+INT((KZ-K1)/STEP)*STEP
       KKP1=KK+STEP!MOD(KK+STEP,NZ)
       BASZ=(RADZ(KZ)-RADZ(KK))/(STEP*DZ)
      DO JY=J1,J2
       JJ=J1+INT((JY-J1)/STEP)*STEP
       JJP1=JJ+STEP!MOD(JJ+STEP,NY)
       BASY=(RADY(JY)-RADY(JJ))/(STEP*DY)
      DO IX=I1,I2
       IF (IX.GE.II1.AND.IX.LE.II2.AND.
     &     JY.GE.JJ1.AND.JY.LE.JJ2.AND.
     &     KZ.GE.KK1.AND.KZ.LE.KK2) CYCLE
       II=I1+INT((IX-I1)/STEP)*STEP
       IIP1=II+STEP!MOD(II+STEP,NX)
       BASX=(RADX(IX)-RADX(II))/(STEP*DX)

       IF (IX.EQ.II.AND.JY.EQ.JJ.AND.KZ.EQ.KK) CYCLE 

       U2(IX,JY,KZ)=U2(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U2(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U2(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U2(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U2(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U2(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U2(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U2(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ

       U3(IX,JY,KZ)=U3(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U3(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U3(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U3(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U3(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U3(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U3(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U3(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ

       U4(IX,JY,KZ)=U4(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              U4(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              U4(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              U4(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              U4(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              U4(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              U4(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              U4(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ     

       L0(IX,JY,KZ)=L0(II,JJ,KK)      *(1-BASX)*(1-BASY)*(1-BASZ) +
     &              L0(IIP1,JJ,KK)    *  BASX  *(1-BASY)*(1-BASZ) +
     &              L0(II,JJP1,KK)    *(1-BASX)*  BASY  *(1-BASZ) +
     &              L0(IIP1,JJP1,KK)  *  BASX  *  BASY  *(1-BASZ) +
     &              L0(II,JJ,KKP1)    *(1-BASX)*(1-BASY)*  BASZ   +
     &              L0(IIP1,JJ,KKP1)  *  BASX  *(1-BASY)*  BASZ   +
     &              L0(II,JJP1,KKP1)  *(1-BASX)*  BASY  *  BASZ   +
     &              L0(IIP1,JJP1,KKP1)*  BASX  *  BASY  *  BASZ  
     
#ifdef use_filter
#if use_filter == 1
       VISC0(IX,JY,KZ)=VISC0(II,JJ,KK)  *(1-BASX)*(1-BASY)*(1-BASZ) +
     &                 VISC0(IIP1,JJ,KK)*  BASX  *(1-BASY)*(1-BASZ) +
     &                 VISC0(II,JJP1,KK)*(1-BASX)* BASY *(1-BASZ) +
     &                 VISC0(IIP1,JJP1,KK)* BASX * BASY *(1-BASZ) +
     &                 VISC0(II,JJ,KKP1)*(1-BASX)*(1-BASY)* BASZ  +
     &                 VISC0(IIP1,JJ,KKP1)* BASX *(1-BASY)* BASZ   +
     &                 VISC0(II,JJP1,KKP1)*(1-BASX)* BASY * BASZ   +
     &                 VISC0(IIP1,JJP1,KKP1)* BASX * BASY * BASZ  
#endif
#endif

      END DO 
      END DO
      END DO

      ! 3. inside this region, we go cell by cell. (inside [II1,II2])

      ! Load balancing
      !ALLOCATE(LB())

!$OMP PARALLEL DO SHARED(NX,NY,NZ,II1,II2,JJ1,JJ2,KK1,KK2,RADX,RADY,
!$OMP+                   RADZ,U2DM,U3DM,U4DM,L0,U2,U3,U4,
!$OMP+                   KNEIGHBOURS,DX,CR0AMR,VISC0,ABVC,
!$OMP+                   TREE,XTREE,YTREE,ZTREE,FLAG_MASS,PI,MASAP,VOL),
!$OMP+            PRIVATE(IX,JY,KZ,DIST,NEIGH,CONTA,H_KERN,BAS8,
!$OMP+                    BAS8X,BAS8Y,BAS8Z,BAS8M,I,BASMASS,DA,SEARCH),
!$OMP+            SCHEDULE(DYNAMIC)!, DEFAULT(NONE)
      DO KZ=KK1,KK2
      DO JY=JJ1,JJ2
            !write(*,*) kz,jy
      DO IX=II1,II2
       IF (CR0AMR(IX,JY,KZ).EQ.1) THEN 
        ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
        DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &        XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &        ZQUERY=DBLE(RADZ(KZ)), K=KNEIGHBOURS)
        DIST = DA%V%VALUES 
        NEIGH = DA%I%VALUES

        IF (DIST(KNEIGHBOURS).GT.DX) THEN
         CONTA=KNEIGHBOURS 
        ELSE 
         DEALLOCATE(DIST,NEIGH)
         DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &         XQUERY=DBLE(RADX(IX)), YQUERY=DBLE(RADY(JY)),
     &         ZQUERY=DBLE(RADZ(KZ)), RADIUS=DBLE(DX))
         CONTA = DA%SIZE() 
         ALLOCATE(DISt(CONTA), NEIGH(CONTA))
         DIST = DA%V%VALUES 
         NEIGH = DA%I%VALUES
        END IF

        !WRITE(*,*) IX,JY,KZ,CONTA,DIST(CONTA)
        H_KERN=DIST(CONTA)
        L0(IX,JY,KZ)=H_KERN

        CALL KERNEL_FUNC(CONTA,CONTA,H_KERN,DIST)
#ifdef weight_scheme
#if weight_scheme == 1
        DO I=1,CONTA
         DIST(I)=DIST(I)*MASAP(NEIGH(I))
        END DO
#elif weight_scheme == 2
        DO I=1,CONTA 
         DIST(I)=DIST(I)*VOL(NEIGH(I))
        END DO
#endif
#endif

        BAS8=0.D0
        BAS8X=0.D0
        BAS8Y=0.D0
        BAS8Z=0.D0
        BAS8M=0.D0
        BASMASS=0.D0
        DO I=1,CONTA 
         BAS8=BAS8+DIST(I)
         BAS8X=BAS8X+DIST(I)*U2DM(NEIGH(I))
         BAS8Y=BAS8Y+DIST(I)*U3DM(NEIGH(I))
         BAS8Z=BAS8Z+DIST(I)*U4DM(NEIGH(I))
#ifdef use_filter
#if use_filter == 1
         BAS8M=BAS8M+DIST(I)*ABVC(NEIGH(I))
#endif
#endif
         !BAS8M=MAX(BAS8M,ABVC(NEIGH(I)))
         BASMASS=BASMASS+MASAP(NEIGH(I))
        END DO
      
        U2(IX,JY,KZ)=BAS8X/BAS8
        U3(IX,JY,KZ)=BAS8Y/BAS8
        U4(IX,JY,KZ)=BAS8Z/BAS8
#ifdef use_filter
#if use_filter == 1
        VISC0(IX,JY,KZ)=BAS8M/BAS8
#endif
#endif
        !VISC0(IX,JY,KZ)=BAS8M

        IF (FLAG_MASS.EQ.1) L0(IX,JY,KZ)=BASMASS/(4*PI/3)/H_KERN**3

        DEALLOCATE(DIST, NEIGH)
       END IF
      END DO 
      END DO 
      END DO

*     AMR levels 
      DO IR=1,NL 
       WRITE(*,*) 'AMR level...',IR
       DXPA=DX/(2.**IR)
       DYPA=DY/(2.**IR)
       DZPA=DZ/(2.**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,CR0AMR1,
!$OMP+                   RX,RY,RZ,U2DM,U3DM,U4DM,L1,U12,U13,U14,
!$OMP+                   KNEIGHBOURS,DXPA,DX,VISC1,ABVC,SOLAP,
!$OMP+                   TREE,XTREE,YTREE,ZTREE,FLAG_MASS,PI,MASAP,VOL),
!$OMP+            PRIVATE(IPATCH,N1,N2,N3,IX,JY,KZ,DIST,NEIGH,DA,
!$OMP+                    CONTA,H_KERN,BAS8,BAS8X,BAS8Y,BAS8Z,BAS8M,I,
!$OMP+                    BASMASS,SEARCH,DO_CELL),
!$OMP+            SCHEDULE(DYNAMIC,1)!, DEFAULT(NONE)
       DO IPATCH=LOW1,LOW2 
            !write(*,*) ir,ipatch
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)

        DO KZ=0,N3+1 
        DO JY=0,N2+1 
        DO IX=0,N1+1
         do_cell = .false.
         if (ix.eq.0.or.ix.eq.n1+1.or.
     &       jy.eq.0.or.jy.eq.n2+1.or.
     &       kz.eq.0.or.kz.eq.n3+1) then
          do_cell = .true.
         else
           if (cr0amr1(ix,jy,kz,ipatch).eq.1.and.
     &         solap(ix,jy,kz,ipatch).eq.1) then
            do_cell = .true.
           end if
         end if
         IF (do_cell) THEN
          !Q(1)=RX(IX,IPATCH)
          !Q(2)=RY(JY,IPATCH)
          !Q(3)=RZ(KZ,IPATCH)
          
          ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
          DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &          XQUERY=DBLE(RX(IX,IPATCH)), YQUERY=DBLE(RY(JY,IPATCH)),
     &          ZQUERY=DBLE(RZ(KZ,IPATCH)), K=KNEIGHBOURS)
          DIST = DA%V%VALUES 
          NEIGH = DA%I%VALUES

          IF (DIST(KNEIGHBOURS).GT.DXPA) THEN
           CONTA=KNEIGHBOURS 
          ELSE 
           DEALLOCATE(DIST,NEIGH)
           DA = SEARCH%KNEAREST(TREE, XTREE, YTREE, ZTREE, 
     &           XQUERY=DBLE(RX(IX,IPATCH)), YQUERY=DBLE(RY(JY,IPATCH)),
     &           ZQUERY=DBLE(RZ(KZ,IPATCH)), RADIUS=DBLE(DXPA))
           CONTA = DA%SIZE() 
           ALLOCATE(DIST(CONTA), NEIGH(CONTA))
           DIST = DA%V%VALUES 
           NEIGH = DA%I%VALUES
          END IF
          !WRITE(*,*) IX,JY,KZ,CONTA,DIST(CONTA)
          H_KERN=DIST(CONTA)
          !  if (isnan(h_kern).or.conta.lt.kneighbours) then 
          !   write(*,*) ipatch,ix,jy,kz,conta,h_kern,dxpa
          !  end if
          CALL KERNEL_FUNC(CONTA,CONTA,H_KERN,DIST)
#ifdef weight_scheme
#if weight_scheme == 1
          DO I=1,CONTA
           DIST(I)=DIST(I)*MASAP(NEIGH(I))
          END DO
#elif weight_scheme == 2
          DO I=1,CONTA 
            DIST(I)=DIST(I)*VOL(NEIGH(I))
          END DO
#endif
#endif
  
          BAS8=0.D0
          BAS8X=0.D0
          BAS8Y=0.D0
          BAS8Z=0.D0
          BAS8M=0.D0
          BASMASS=0.D0
          DO I=1,CONTA 
           BAS8=BAS8+DIST(I)
           BAS8X=BAS8X+DIST(I)*U2DM(NEIGH(I))
           BAS8Y=BAS8Y+DIST(I)*U3DM(NEIGH(I))
           BAS8Z=BAS8Z+DIST(I)*U4DM(NEIGH(I))
#ifdef use_filter
#if use_filter == 1
           BAS8M=BAS8M+DIST(I)*ABVC(NEIGH(I))
#endif
#endif
           !BAS8M=MAX(BAS8M,ABVC(NEIGH(I)))
           BASMASS=BASMASS+MASAP(NEIGH(I))
          END DO

          DEALLOCATE(DIST, NEIGH)
        
          U12(IX,JY,KZ,IPATCH)=BAS8X/BAS8
          U13(IX,JY,KZ,IPATCH)=BAS8Y/BAS8
          U14(IX,JY,KZ,IPATCH)=BAS8Z/BAS8

          if (ix.eq.0.or.ix.eq.n1+1.or.
     &        jy.eq.0.or.jy.eq.n2+1.or.
     &        kz.eq.0.or.kz.eq.n3+1) cycle ! these values are not saved,
      !   we actually only want u12,u13,u14 in ghost zones for taking 
      !   derivatives

#ifdef use_filter
#if use_filter == 1
          VISC1(IX,JY,KZ,IPATCH)=BAS8M/BAS8
#endif
#endif
          !VISC1(IX,JY,KZ,IPATCH)=BAS8M

          if (flag_mass.eq.0) then
            L1(IX,JY,KZ,IPATCH)=H_KERN
          else 
            L1(IX,JY,KZ,IPATCH)=BASMASS/(4*PI/3)/H_KERN**3
          end if

         END IF
        END DO 
        END DO 
        END DO
       END DO 
      END DO

      CALL TREE%DEALLOCATE()

*     refill refined and overlapping cells
      DO IR=NL,1,-1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    L1(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U12(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U13(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    U14(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
#ifdef use_filter
#if use_filter == 1
        CALL SYNC_AMR_FILTER(IR,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &    PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,
     &    VISC1(1:NAMRX,1:NAMRY,1:NAMRZ,:),NL)
#endif
#endif

        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
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
             uw(1:2,1:2,1:2) = 1.

             u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             L1(II,JJ,KK,JPATCH) = FUIN

             u(1:2,1:2,1:2) = u12(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u12(II,JJ,KK,JPATCH) = FUIN

             u(1:2,1:2,1:2) = u13(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u13(II,JJ,KK,JPATCH) = FUIN

             u(1:2,1:2,1:2) = u14(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u14(II,JJ,KK,JPATCH) = FUIN

#ifdef use_filter
#if use_filter == 1
             u(1:2,1:2,1:2) = VISC1(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             VISC1(II,JJ,KK,JPATCH) = FUIN
#endif
#endif
            else
             uw(1:2,1:2,1:2) = 1.

             u(1:2,1:2,1:2) = L1(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             L0(II,JJ,KK) = FUIN

             u(1:2,1:2,1:2) = u12(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u2(II,JJ,KK) = FUIN

             u(1:2,1:2,1:2) = u13(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u3(II,JJ,KK) = FUIN

             u(1:2,1:2,1:2) = u14(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             u4(II,JJ,KK) = FUIN

#ifdef use_filter
#if use_filter == 1
             u(1:2,1:2,1:2) = VISC1(I:I+1,J:J+1,K:K+1,IPATCH)
             call finer_to_coarser(u,uw,fuin)
             VISC0(II,JJ,KK) = FUIN
#endif
#endif
            end if
          END DO
          END DO
          END DO
        END DO
      END DO !IR=NL,1,-1

      write(*,*) 'At level', 0
      CALL P_MINMAX_IR(L0,L1,1,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                 NPATCH,0,BASX,BASY)
      write(*,*) 'L min,max',BASX,BASY
      CALL P_MINMAX_IR(U2,U12,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                 NPATCH,0,BASX,BASY)
      write(*,*) 'vx min,max',BASX,BASY
      CALL P_MINMAX_IR(U3,U13,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                 NPATCH,0,BASX,BASY)
      write(*,*) 'vy min,max',BASX,BASY
       CALL P_MINMAX_IR(U4,U14,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                  NPATCH,0,BASX,BASY)
      write(*,*) 'vz min,max',BASX,BASY
#ifdef use_filter
#if use_filter == 1
      CALL P_MINMAX_IR(VISC0,VISC1,1,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,
     &                PATCHNZ,NPATCH,0,BASX,BASY)
      IF (FLAG_MACHFIELD.EQ.0) THEN
       WRITE(*,*) 'ABVC min,max',BASX,BASY
      ELSE 
       WRITE(*,*) 'Mach min,max',BASX,BASY
      END IF
#endif
#endif


      DO IR=1,NL
       low1=sum(npatch(0:ir-1))+1
       low2=sum(npatch(0:ir))
       write(*,*) 'At level',ir
       CALL P_MINMAX_IR(L0,L1,1,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                  NPATCH,IR,BASX,BASY)
       write(*,*) 'L min,max',BASX,BASY
       CALL P_MINMAX_IR(U2,U12,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                  NPATCH,IR,BASX,BASY)
       write(*,*) 'vx min,max',BASX,BASY
       CALL P_MINMAX_IR(U3,U13,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                  NPATCH,IR,BASX,BASY)
       write(*,*) 'vy min,max',BASX,BASY
       CALL P_MINMAX_IR(U4,U14,1,1,NX,NY,NZ,NL,PATCHNX,PATCHNY,PATCHNZ,
     &                  NPATCH,IR,BASX,BASY)
       write(*,*) 'vz min,max',BASX,BASY
#ifdef use_filter
#if use_filter == 1
       CALL P_MINMAX_IR(VISC0,VISC1,1,0,NX,NY,NZ,NL,PATCHNX,PATCHNY,
     &                PATCHNZ,NPATCH,IR,BASX,BASY)
       IF (FLAG_MACHFIELD.EQ.0) THEN
        WRITE(*,*) 'ABVC min,max',BASX,BASY
       ELSE 
        WRITE(*,*) 'Mach min,max',BASX,BASY
       END IF
#endif
#endif
      END DO

#ifdef output_grid 
#if output_grid == 1
      CALL WRITE_GRID_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                          PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                          PATCHRY,PATCHRZ,PARE,CR0AMR,CR0AMR1,
     &                          SOLAP,L0,L1,VISC0,VISC1,FLAG_MACHFIELD)
#endif
#endif

      RETURN
      END
