************************************************************************
      SUBROUTINE LEE_MACH(ITER,NPATCH,PARE,PATCHNX,PATCHNY,
     &            PATCHNZ,SHOCK0,SHOCK1)
************************************************************************
*     Only used for the filter. Reads the outputs of a shock finder,
*     and identifies the shocked cells (M >= 1.3)
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     Parameters: input
      INTEGER ITER, NPATCH(0:NLEVELS), PARE(NPALEV), PATCHNX(NPALEV),
     &        PATCHNY(NPALEV), PATCHNZ(NPALEV)

*     Parameters: output
      INTEGER SHOCK0(NMAX,NMAY,NMAZ), SHOCK1(NAMRX,NAMRY,NAMRZ,NPALEV)

*     Private variables
      INTEGER I,J,K,IPATCH,N1,N2,N3
      REAL BAS, THR
      CHARACTER*20 FILNOM,FIL1
      real*4, allocatable::scr4(:,:,:)

      thr = 3.0 ! mach no. threshold (>thr --> shocked, <thr --> unshocked)

      CALL NOMFILEMACH5(ITER,FILNOM)
      FIL1='shocks/'//FILNOM
      OPEN (31,FILE=FIL1,
     &       STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

      shock0 = 0
      allocate(scr4(nmax,nmay,nmaz))
      read(31) (((scr4(i,j,k),i=1,n2),j=1,n2),k=1,n3)
      n1 = nmax
      n2 = nmay
      n3 = nmaz
      do i=1,n1
        do j=1,n2
          do k=1,n3
            if (scr4(i,j,k).ge.thr) shock0(i,j,k) = 1
          end do
        end do
      end do
      deallocate(scr4)

      allocate(scr4(namrx,namry,namrz))
      do ipatch=1,sum(npatch)
        n1 = patchnx(ipatch)
        n2 = patchny(ipatch)
        n3 = patchnz(ipatch)
        read(31) (((scr4(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        shock1(:,:,:,ipatch) = 0
        do i=1,n1
          do j=1,n2
            do k=1,n3
              if (scr4(i,j,k).ge.thr) shock1(i,j,k,ipatch) = 1
            end do
          end do
        end do
      end do
      deallocate(scr4)

      RETURN
      END

***********************************************************************
       SUBROUTINE READ_GADGET(ITER,FILES_PER_SNAP,NX,NY,NZ,T,ZETA,
     &            NL_PARTICLE_GRID,REFINE_THR,PARCHLIM,BORGRID,
     &            NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &            PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,LADO0,
     &            NPART,RXPA,RYPA,RZPA,MASAP,U2DM,U3DM,U4DM)
***********************************************************************
*     Reads the GAS particles of the simulation, builds a set of AMR
*     grids and interpolates a continuous velocity field.
*     This subroutine may be changed for different simulation codes.
*     In particular, after mesh-building and interpolation, in any
*     case we'll need to feed the main code with:
****  GRIDS info:
*     NPATCH: number of patches per refinement level
*     PATCHNX, PATCHNY, PATCHNZ: cell extensions of each refinement patch
*     PATCHX, PATCHY, PATCHZ: grid coordinates of the leftmost cell of each patch
*     PATCHRX, PATCHRY, PATCHRZ: origin position of each patch (position of the leftmost "mother" cell)
*     PARE: coarser patch a given patch is embedded in
****  CLUS info:
*     U2, U3, U4: initial velocity field (base level, i.e coarse grid)
*     U12, U13, U14: initial velocity field (refinement patches)
*     CR0AMR: whether a cell is refined (=0) or it isn't (=1)
***********************************************************************

       USE gadget_read
       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ,LOW1,LOW2,FILES_PER_SNAP
       real T,AAA,BBB,CCC,MAP,ZETA,LADO0
       INTEGER I,J,K,IX,NL,IR,IRR,N1,N2,N3,NL_PARTICLE_GRID
       INTEGER REFINE_THR,PARCHLIM,BORGRID

       INTEGER FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &         FLAG_W_VELOCITIES
       COMMON /FLAGS/ FLAG_VERBOSE, FLAG_W_DIVROT, FLAG_W_POTENTIALS,
     &                FLAG_W_VELOCITIES

       INTEGER NPATCH(0:NLEVELS),PARE(NPALEV)
       integer NPART(0:NLEVELS)
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
       real PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)

       CHARACTER*200 FIL1,FIL2

!      this might be only for testing
       !real U1(1:NMAX,1:NMAY,1:NMAZ)
       !real U11(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       !common /dens/ u1,u11

       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC/ U2,U3,U4,U12,U13,U14

!      Particles
       REAL*4 RXPA(NDM),RYPA(NDM),RZPA(NDM),
     &        U2DM(NDM),U3DM(NDM),U4DM(NDM),MASAP(NDM),
     &        KERNEL(NDM)

       integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
       integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /cr0/ cr0amr, cr0amr1

       INTEGER LIHAL(NDM),LIHAL_IX(NDM),LIHAL_JY(NDM),LIHAL_KZ(NDM)
*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       REAL DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       REAL CIO_XC,CIO_YC,CIO_ZC
       COMMON /DOM_DECOMP/ DDXL,DDXR,DDYL,DDYR,DDZL,DDZR

       CHARACTER*3 ITER_STRING
       CHARACTER*1 IFILE_STRING
       INTEGER IFILE

       ! Scratch variables for reading
       REAL*4,ALLOCATABLE::SCR4(:)
       REAL*4,ALLOCATABLE::SCR42(:,:)
       INTEGER,ALLOCATABLE::ELIM(:)
       integer npart_gadget(6),nall(6),blocksize
       real*8 massarr(6)
       integer basint1,basint2,basint3,basint4,basint5,basint6,basint7
       real*4 bas41,bas42,bas43,bas44,bas45,bas46,bas47
       real*8 bas81,bas82,bas83,bas84,bas85,bas86,bas87
       character*4 blocklabel
       ! End scratch variables for reading

       real xmin,ymin,zmin,xmax,ymax,zmax

       NPART(:)=0

       LOW2=0
       DO IFILE=0,FILES_PER_SNAP-1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*       READING DATA
        WRITE(ITER_STRING,'(I3.3)') ITER
        FIL1='./simulation/snapdir_'//ITER_STRING//'/snap_'//ITER_STRING
        IF (FILES_PER_SNAP.EQ.1) THEN
          FIL2=FIL1
        ELSE
          WRITE(IFILE_STRING,'(I1.1)') IFILE
          FIL2=TRIM(ADJUSTL(FIL1))//'.'//IFILE_STRING
        END IF

        WRITE(*,*) 'Reading iteration file: ',ITER,' ',
     &              TRIM(ADJUSTL(FIL2))

        CALL read_head(FIL2,npart_gadget,massarr,bas81,bas82,
     &                  basint1,basint2,nall,basint3,basint4,bas83,
     &                  bas84,bas85,bas86,blocksize)
        WRITE(*,*) npart_gadget(1), 'gas particles'
        LOW1=LOW2+1
        LOW2=LOW1+NPART_GADGET(1)-1
        
        ALLOCATE(SCR42(3,SUM(NPART_GADGET(1:6))))
        WRITE(*,*) 'Reading positions ...'
        CALL read_float3(FIL2,'POS ',SCR42,blocksize)
        WRITE(*,*) ' found for ',(blocksize-8)/12,' particles'
        RXPA(LOW1:LOW2)=SCR42(1,1:NPART_GADGET(1))
        RYPA(LOW1:LOW2)=SCR42(2,1:NPART_GADGET(1))
        RZPA(LOW1:LOW2)=SCR42(3,1:NPART_GADGET(1))

        WRITE(*,*) 'Reading velocities ...'
        CALL read_float3(FIL2,'VEL ',SCR42,blocksize)
        WRITE(*,*) ' found for ',(blocksize-8)/12,' particles'
        U2DM(LOW1:LOW2)=SCR42(1,1:NPART_GADGET(1))
        U3DM(LOW1:LOW2)=SCR42(2,1:NPART_GADGET(1))
        U4DM(LOW1:LOW2)=SCR42(3,1:NPART_GADGET(1))

        DEALLOCATE(SCR42)

        ALLOCATE(SCR4(SUM(NPART_GADGET(1:6))))
        WRITE(*,*) 'Reading masses ...'
        CALL read_float(FIL2,'MASS',SCR4,blocksize)
        WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
        MASAP(LOW1:LOW2)=SCR4(1:NPART_GADGET(1))

        WRITE(*,*) 'Reading filter length ...'
        CALL read_float(FIL2,'HSML',SCR4,blocksize)
        WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
        KERNEL(LOW1:LOW2)=SCR4(1:NPART_GADGET(1))        
        DEALLOCATE(SCR4)

       END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
       NPART(0)=LOW2 !Retrocompatibility with general reader

       LOW1=1
       LOW2=SUM(NPART(0:NLEVELS))
       WRITE(*,*)
       WRITE(*,*) 'Input particles:',LOW2
       xmin=minval(RXPA(LOW1:LOW2))
       xmax=maxval(RXPA(LOW1:LOW2))
       ymin=minval(RYPA(LOW1:LOW2))
       ymax=maxval(RYPA(LOW1:LOW2))
       zmin=minval(RZPA(LOW1:LOW2))
       zmax=maxval(RZPA(LOW1:LOW2))

       WRITE(*,*) 'RXPA =',xmin,xmax
       WRITE(*,*) 'RYPA =',ymin,ymax
       WRITE(*,*) 'RZPA =',zmin,zmax
       WRITE(*,*) 'U2DM =',MINVAL(U2DM(LOW1:LOW2)),
     &                     MAXVAL(U2DM(LOW1:LOW2))
       WRITE(*,*) 'U3DM =',MINVAL(U3DM(LOW1:LOW2)),
     &                     MAXVAL(U3DM(LOW1:LOW2))
       WRITE(*,*) 'U4DM =',MINVAL(U4DM(LOW1:LOW2)),
     &                     MAXVAL(U4DM(LOW1:LOW2))
       WRITE(*,*) 'MASAP=',MINVAL(MASAP(LOW1:LOW2)),
     &                     MAXVAL(MASAP(LOW1:LOW2))
       WRITE(*,*) 'KERNEL LENGTH=',MINVAL(KERNEL(LOW1:LOW2)),
     &                             MAXVAL(KERNEL(LOW1:LOW2))

       IF (XMIN.LT.DDXL.OR.XMAX.GT.DDXR.OR.
     &     YMIN.LT.DDYL.OR.YMAX.GT.DDYR.OR.
     &     ZMIN.LT.DDZL.OR.ZMAX.GT.DDZR) THEN

        WRITE(*,*)
        WRITE(*,*) 'WARNING: PARTICLES OUTSIDE THE DOMAIN'
        ALLOCATE(ELIM(LOW1:LOW2))
        
!$OMP PARALLEL DO SHARED(RXPA,RYPA,RZPA,LOW1,LOW2,DDXL,DDXR,DDYL,DDYR,
!$OMP+            DDZL,DDZR,ELIM), PRIVATE(I), DEFAULT(NONE)
        DO I=LOW1,LOW2 
          ELIM(I)=0
          IF (RXPA(I).LT.DDXL.OR.RXPA(I).GT.DDXR.OR.
     &        RYPA(I).LT.DDYL.OR.RYPA(I).GT.DDYR.OR.
     &        RZPA(I).LT.DDZL.OR.RZPA(I).GT.DDZR) THEN
            ELIM(I)=1
          END IF
        END DO

        J=0
        DO I=LOW1,LOW2 
          IF (ELIM(I).EQ.0) THEN
            J=J+1
            IF (I.EQ.J) CYCLE
            RXPA(J)=RXPA(I)
            RYPA(J)=RYPA(I)
            RZPA(J)=RZPA(I)
            U2DM(J)=U2DM(I)
            U3DM(J)=U3DM(I)
            U4DM(J)=U4DM(I)
            MASAP(J)=MASAP(I)
            KERNEL(J)=KERNEL(I)
          END IF
        END DO
        DEALLOCATE(ELIM)

        LOW2=J

       END IF

       CIO_XC=0.5*(DDXL+DDXR)
       CIO_YC=0.5*(DDYL+DDYR)
       CIO_ZC=0.5*(DDZL+DDZR)
       IF (ABS(CIO_XC).GT.1.E-6*MAX(DDXL,DDXR).OR.
     &     ABS(CIO_YC).GT.1.E-6*MAX(DDYL,DDYR).OR.
     &     ABS(CIO_ZC).GT.1.E-6*MAX(DDZL,DDZR)) THEN

!$OMP PARALLEL DO SHARED(RXPA,RYPA,RZPA,LOW1,LOW2,CIO_XC,CIO_YC,CIO_ZC),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
        DO I=LOW1,LOW2 
          RXPA(I)=RXPA(I)-CIO_XC
          RYPA(I)=RYPA(I)-CIO_YC
          RZPA(I)=RZPA(I)-CIO_ZC
        END DO

       END IF

       WRITE(*,*)
       WRITE(*,*) 'After recentering and domain decomposition',LOW2
       xmin=minval(RXPA(LOW1:LOW2))
       xmax=maxval(RXPA(LOW1:LOW2))
       ymin=minval(RYPA(LOW1:LOW2))
       ymax=maxval(RYPA(LOW1:LOW2))
       zmin=minval(RZPA(LOW1:LOW2))
       zmax=maxval(RZPA(LOW1:LOW2))

       WRITE(*,*) 'RXPA =',xmin,xmax
       WRITE(*,*) 'RYPA =',ymin,ymax
       WRITE(*,*) 'RZPA =',zmin,zmax
       WRITE(*,*) 'U2DM =',MINVAL(U2DM(LOW1:LOW2)),
     &                     MAXVAL(U2DM(LOW1:LOW2))
       WRITE(*,*) 'U3DM =',MINVAL(U3DM(LOW1:LOW2)),
     &                     MAXVAL(U3DM(LOW1:LOW2))
       WRITE(*,*) 'U4DM =',MINVAL(U4DM(LOW1:LOW2)),
     &                     MAXVAL(U4DM(LOW1:LOW2))
       WRITE(*,*) 'MASAP=',MINVAL(MASAP(LOW1:LOW2)),
     &                     MAXVAL(MASAP(LOW1:LOW2))
       WRITE(*,*) 'KERNEL LENGTH=',MINVAL(KERNEL(LOW1:LOW2)),
     &                             MAXVAL(KERNEL(LOW1:LOW2))
       

       WRITE(*,*) 'Routine create mesh ------------------------------'
       NPATCH(0:IR)=0
       CALL CREATE_MESH(NX,NY,NZ,NL_PARTICLE_GRID,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &            U4DM,MASAP,NPART,LADO0,REFINE_THR,PARCHLIM,BORGRID)
       DO IR=1,NL_PARTICLE_GRID
        IF (NPATCH(IR).EQ.0) EXIT
       END DO
       NL=IR-1

       CALL GRIDAMR(NX,NY,NZ,NL,NPATCH,
     &                   PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,
     &                   PATCHRX,PATCHRY,PATCHRZ,PARE)
       WRITE(*,*) 'End mesh creation --------------------------------'

       WRITE(*,*) 'Routine interpolate velocity ---------------------'
       CALL INTERPOLATE_VELOCITIES(NX,NY,NZ,NL,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,RXPA,RYPA,RZPA,U2DM,U3DM,
     &            U4DM,MASAP,NPART,LADO0)

       WRITE(*,*) 'Locating particles onto the grid'
       CALL PLACE_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,LIHAL_KZ)

       CALL ERROR_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,RXPA,RYPA,RZPA,
     &            U2DM,U3DM,U4DM,NPART,LADO0,LIHAL,LIHAL_IX,LIHAL_JY,
     &            LIHAL_KZ)
       WRITE(*,*) 'End velocity interpolation -----------------------'

       RETURN
       END