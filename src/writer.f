#ifdef output_grid
#if output_grid == 1
***********************************************************************
       SUBROUTINE WRITE_DIVROT(NX,NY,NZ,ITER,T,ZETA,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the divergence and each component of the rotational
*     of the velocity field to a file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

      INTEGER FLAG_VERBOSE
      INTEGER FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL
      INTEGER FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &         FL_GR_DIV,FL_GR_CURL
      INTEGER FL_P_ERR,FL_P_RES
      INTEGER FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB
      real fl_smooth_filtlen
      COMMON /FLAGS/ FLAG_VERBOSE,FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL,
     &        FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &        FL_GR_DIV,FL_GR_CURL,FL_P_ERR,FL_P_RES,
     &        FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB,
     &        fl_smooth_filtlen

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)
      real*4 scrvar1, scrvar2

      CHARACTER*5 ITER_STRING
      CHARACTER*200 FILENOM

      WRITE(ITER_STRING,'(I5.5)') ITER

      IF (FL_GR_DIV.EQ.1) THEN
       FILENOM='output_files/divv'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
        WRITE(99) (((DIVER0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO I=1,LOW1 
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         WRITE(99) (((DIVER(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

      IF (FL_GR_CURL.EQ.1) THEN
       FILENOM='output_files/curlv'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
       WRITE(99) (((ROTAX_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAY_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAZ_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       LOW1=SUM(NPATCH(0:NL))
       DO I=1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        WRITE(99) (((ROTAX_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAY_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAZ_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
       END DO
       CLOSE(99)
      END IF

      RETURN
      END


***********************************************************************
       SUBROUTINE WRITE_POTENTIALS(NX,NY,NZ,ITER,T,ZETA,NL,
     &            NPATCH,PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the scalar and the vector potentials to a file
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

      INTEGER FLAG_VERBOSE
      INTEGER FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL
      INTEGER FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &         FL_GR_DIV,FL_GR_CURL
      INTEGER FL_P_ERR,FL_P_RES
      INTEGER FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB
      real fl_smooth_filtlen
      COMMON /FLAGS/ FLAG_VERBOSE,FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL,
     &        FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &        FL_GR_DIV,FL_GR_CURL,FL_P_ERR,FL_P_RES,
     &        FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB,
     &        fl_smooth_filtlen

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)
      real*4 scrvar1, scrvar2

      CHARACTER*5 ITER_STRING
      CHARACTER*200 FILENOM

      WRITE(ITER_STRING,'(I5.5)') ITER

      IF (FL_GR_SPOT.EQ.1) THEN
       FILENOM='output_files/spot'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
        WRITE(99) (((DIVER0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO I=1,LOW1 
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         WRITE(99) (((DIVER(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

      IF (FL_GR_VPOT.EQ.1) THEN
       FILENOM='output_files/vpot'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
       WRITE(99) (((ROTAX_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAY_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAZ_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       LOW1=SUM(NPATCH(0:NL))
       DO I=1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        WRITE(99) (((ROTAX_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAY_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAZ_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
       END DO
       CLOSE(99)
      END IF

      RETURN
      END


***********************************************************************
       SUBROUTINE WRITE_VELOCITIES(NX,NY,NZ,ITER,T,ZETA,NL,
     &                             NPATCH,PATCHNX,PATCHNY,PATCHNZ)
***********************************************************************
*     Writes the total, the compressive and the rotational velocities
*     to a file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      INTEGER NX, NY, NZ, NL, ITER
      real T, ZETA
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)

*     INPUTS FROM COMMON MODULES
*     ROTA has been recycled as the rotational velocity
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

*      original, total velocity
       real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
       COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14

      INTEGER FLAG_VERBOSE
      INTEGER FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL
      INTEGER FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &         FL_GR_DIV,FL_GR_CURL
      INTEGER FL_P_ERR,FL_P_RES
      INTEGER FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB
      real fl_smooth_filtlen
      COMMON /FLAGS/ FLAG_VERBOSE,FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL,
     &        FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &        FL_GR_DIV,FL_GR_CURL,FL_P_ERR,FL_P_RES,
     &        FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB,
     &        fl_smooth_filtlen

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)
      real*4 scrvar1, scrvar2

      CHARACTER*5 ITER_STRING
      CHARACTER*200 FILENOM

      WRITE(ITER_STRING,'(I5.5)') ITER

      IF (FL_GR_VCOMP.EQ.1) THEN
       FILENOM='output_files/vcomp'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
       WRITE(99) (((U2P(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((U3P(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((U4P(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       LOW1=SUM(NPATCH(0:NL))
       DO I=1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        WRITE(99) (((U12P(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((U13P(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((U14P(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
       END DO
       CLOSE(99)
      END IF

      IF (FL_GR_VSOL.EQ.1) THEN
       FILENOM='output_files/vsol'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
       WRITE(99) (((ROTAX_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAY_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       WRITE(99) (((ROTAZ_0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       LOW1=SUM(NPATCH(0:NL))
       DO I=1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        WRITE(99) (((ROTAX_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAY_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((ROTAZ_1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
       END DO
       CLOSE(99)
      END IF

      RETURN
      END
#endif
#endif

#ifdef output_grid 
#if output_grid == 1
**********************************************************************
      SUBROUTINE WRITE_GRID_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,
     &                          PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &                          PATCHRX,PATCHRY,PATCHRZ,PARE,CR0AMR,
     &                          CR0AMR1,SOLAP,L0,L1,MACH0,MACH1,
     &                          FLAG_MACHFIELD)
***********************************************************************
*     Writes the GRIDS and CR0AMR/SOLAP variables for the created AMR
*     structure
***********************************************************************
      IMPLICIT NONE
      INCLUDE 'vortex_parameters.dat'

      INTEGER NL,NX,NY,NZ,FLAG_MACHFIELD
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      INTEGER PARE(NPALEV)

      INTEGER cr0amr(NMAX,NMAY,NMAZ)
      INTEGER cr0amr1(NAMRX,NAMRY,NAMRZ,NPALEV)
      INTEGER solap(NAMRX,NAMRY,NAMRZ,NPALEV)

      REAL L0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL L1(NAMRX,NAMRY,NAMRZ,NPALEV)

#ifdef use_filter 
#if use_filter == 1
      REAL MACH0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL MACH1(NAMRX,NAMRY,NAMRZ,NPALEV)
#else 
      REAL MACH0, MACH1 
#endif
#endif

      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

      INTEGER IX,JY,KZ,I,LOW1,LOW2,IR,IPATCH,J,K,N1,N2,N3

      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      ! runtime IO flags
      INTEGER FLAG_VERBOSE
      INTEGER FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL
      INTEGER FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &         FL_GR_DIV,FL_GR_CURL
      INTEGER FL_P_ERR,FL_P_RES
      INTEGER FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB
      real fl_smooth_filtlen
      COMMON /FLAGS/ FLAG_VERBOSE,FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL,
     &        FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &        FL_GR_DIV,FL_GR_CURL,FL_P_ERR,FL_P_RES,
     &        FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB,
     &        fl_smooth_filtlen

      CHARACTER*5 ITER_STRING
      CHARACTER*200 FILENOM
      WRITE(ITER_STRING,'(I5.5)') ITER

!     Grid description is written always
      FILENOM='output_files/grids'//ITER_STRING
      write(*,*) 'Writing grids data', FILENOM
      OPEN (23,FILE=FILENOM,STATUS='UNKNOWN')
       WRITE(23,*) ITER,' ',0.0,' ',NL,' ',0.0,' ',0.0
       WRITE(23,*) 0.0
       IR=0
       WRITE(23,*) IR,0,0,NX,NY,NZ
       DO IR=1,NL
        WRITE(23,*) IR,' ',NPATCH(IR),' ',0,' ',0,' ',0
        WRITE(23,*) '----------------- within level=',IR,' -----------'
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
         WRITE(23,*) PATCHNX(I),' ',PATCHNY(I),' ',PATCHNZ(I)
         WRITE(23,*) PATCHX(I),' ',PATCHY(I),' ',PATCHZ(I)
         WRITE(23,*) PATCHRX(I),' ',PATCHRY(I),' ',PATCHRZ(I)
         WRITE(23,*) PARE(I)
        END DO
       END DO
      CLOSE(23)

!     Grid overlaps are also always written 
      FILENOM='output_files/grid_overlaps'//ITER_STRING
      OPEN(99,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')
       WRITE(99) (((CR0AMR(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
       LOW1=SUM(NPATCH(0:NL))
       DO IPATCH=1,LOW1
        N1=PATCHNX(IPATCH)
        N2=PATCHNY(IPATCH)
        N3=PATCHNZ(IPATCH)
        WRITE(99) (((CR0AMR1(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
        WRITE(99) (((SOLAP(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
       END DO 
      CLOSE(99)

!     Write smoothing length OR density 
      IF (FL_GR_DEN.EQ.1.OR.FL_GR_KERNL.EQ.1) THEN 
       FILENOM='output_files/gridded_density'//ITER_STRING
       IF (FL_GR_KERNL.EQ.1) 
     &       FILENOM='output_files/gridded_kernel_length'//ITER_STRING
       OPEN(99,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(99) (((L0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO IPATCH=1,LOW1
         N1=PATCHNX(IPATCH)
         N2=PATCHNY(IPATCH)
         N3=PATCHNZ(IPATCH)
         WRITE(99) (((L1(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

!     Write velocities
      IF (FL_GR_VEL.EQ.1) THEN 
       FILENOM='output_files/gridded_velocity'//ITER_STRING
       OPEN(99,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(99) (((U2(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        WRITE(99) (((U3(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        WRITE(99) (((U4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO IPATCH=1,LOW1
         N1=PATCHNX(IPATCH)
         N2=PATCHNY(IPATCH)
         N3=PATCHNZ(IPATCH)
         WRITE(99) (((U12(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
         WRITE(99) (((U13(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
         WRITE(99) (((U14(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

!     Write mach number field
#ifdef use_filter
#if use_filter == 1
      IF (FL_FILT_MACH.EQ.1) THEN 
       IF (FLAG_MACHFIELD.EQ.1) THEN 
        FILENOM='output_files/gridded_mach'//ITER_STRING
       ELSE 
        FILENOM='output_files/gridded_abvc'//ITER_STRING
       END IF
       OPEN(99,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(99) (((MACH0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO IPATCH=1,LOW1
         N1=PATCHNX(IPATCH)
         N2=PATCHNY(IPATCH)
         N3=PATCHNZ(IPATCH)
         WRITE(99) (((MACH1(I,J,K,IPATCH),I=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF
#endif 
#endif

      RETURN
      END
#endif 
#endif

#ifdef output_particles 
#if output_particles == 1
*********************************************************************
      SUBROUTINE WRITE_PARTICLES(NL,NX,NY,NZ,NPATCH,PATCHNX,PATCHNY,
     &                           PATCHNZ,PATCHX,PATCHY,PATCHZ,PATCHRX,
     &                           PATCHRY,PATCHRZ,PARE,NPART,LADO0,
     &                           parchlim)
***********************************************************************
*     Writes the GRIDS and CR0AMR/SOLAP variables for the created AMR
*     structure
***********************************************************************
      use particle_data
      IMPLICIT NONE
      INCLUDE 'vortex_parameters.dat'

      INTEGER NL,NX,NY,NZ
      INTEGER NPATCH(0:NLEVELS),NPART(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV),LADO0
      INTEGER PARE(NPALEV)
      integer parchlim

      REAL SCRPART(PARTI)
      REAL SCR0(NMAX,NMAY,NMAZ)
      REAL SCR1(NAMRX,NAMRY,NAMRZ,NPALEV)

*     INPUTS FROM COMMON MODULES
*     ROTA has been recycled as the rotational velocity
      real ROTAX_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAY_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAZ_0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real ROTAX_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAY_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      real ROTAZ_1(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /ROTS/ ROTAX_0,ROTAY_0,ROTAZ_0,ROTAX_1,ROTAY_1,ROTAZ_1

      real U2P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4P(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14P(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_P/ U2P,U3P,U4P,U12P,U13P,U14P

*      original, total velocity
      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC_ORIGINAL/ U2,U3,U4,U12,U13,U14

      INTEGER IX,JY,KZ,I,LOW1,LOW2,IR,IP,NPARTTOT

      CHARACTER*5 ITER_STRING 
      CHARACTER*200 FILENOM

      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      WRITE(ITER_STRING,'(I5.5)') ITER
      FILENOM='output_files/velocity-particles'//ITER_STRING


      NPARTTOT=SUM(NPART(0:NLEVELS))
      OPEN(99,FILE=FILENOM,STATUS='UNKNOWN',FORM='UNFORMATTED')

      WRITE(99) NPARTTOT
      WRITE(99) (RXPA(I),I=1,NPARTTOT)
      WRITE(99) (RYPA(I),I=1,NPARTTOT)
      WRITE(99) (RZPA(I),I=1,NPARTTOT)
      WRITE(99) (U2DM(I),I=1,NPARTTOT)
      WRITE(99) (U3DM(I),I=1,NPARTTOT)
      WRITE(99) (U4DM(I),I=1,NPARTTOT)
      WRITE(99) (MASAP(I),I=1,NPARTTOT)

      CALL PLACE_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,parchlim)


      SCR0(1:NX,1:NY,1:NZ)=U2(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U12(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=U3(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U13(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=U4(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U14(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=U2P(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U12P(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=U3P(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U13P(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=U4P(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=U14P(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=ROTAX_0(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=ROTAX_1(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=ROTAY_0(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=ROTAY_1(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      SCR0(1:NX,1:NY,1:NZ)=ROTAZ_0(1:NX,1:NY,1:NZ)
      SCR1(1:NAMRX,1:NAMRY,1:NAMRZ,:)=ROTAZ_1(1:NAMRX,1:NAMRY,1:NAMRZ,:)
      CALL GRID_TO_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &            PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &            NPART,LADO0,SCR0,SCR1,SCRPART)
      WRITE(99) (SCRPART(I),I=1,NPARTTOT)

      CLOSE(99)

      DEALLOCATE(LIHAL, LIHAL_IX, LIHAL_JY, LIHAL_KZ)

      RETURN
      END
#endif
#endif

#ifdef output_filter 
#if output_filter == 1
**********************************************************************
       SUBROUTINE WRITE_FILTLEN(NX,NY,NZ,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ,L0,L1)
***********************************************************************
*     Writes the filter length to a separate file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      INTEGER NX, NY, NZ, NL
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      real L0(1:NMAX,1:NMAY,1:NMAZ)
      real L1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

      real U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      real U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      real U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14
*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3
      real*4, ALLOCATABLE::SCR4(:,:,:)

      ! runtime IO flags
      INTEGER FLAG_VERBOSE
      INTEGER FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL
      INTEGER FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &         FL_GR_DIV,FL_GR_CURL
      INTEGER FL_P_ERR,FL_P_RES
      INTEGER FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB
      real fl_smooth_filtlen
      COMMON /FLAGS/ FLAG_VERBOSE,FL_GR_KERNL,FL_GR_DEN,FL_GR_VEL,
     &        FL_GR_VCOMP,FL_GR_VSOL,FL_GR_SPOT,FL_GR_VPOT,
     &        FL_GR_DIV,FL_GR_CURL,FL_P_ERR,FL_P_RES,
     &        FL_FILT_MACH,FL_FILT_SHOCK,FL_FILT_LEN,FL_FILT_VTURB,
     &        fl_smooth_filtlen

      INTEGER NXBAS,NYBAS,NZBAS,ITER
      COMMON /ITERI/ NXBAS,NYBAS,NZBAS,ITER

      CHARACTER*5 ITER_STRING
      CHARACTER*200 FILENOM
      WRITE(ITER_STRING,'(I5.5)') ITER

      IF (FL_FILT_LEN.EQ.1) THEN 
       FILENOM='output_files/gridded_filtlen'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
        WRITE(99) (((L0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO I=1,LOW1
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         WRITE(99) (((L1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

      IF (FL_FILT_VTURB.EQ.1) THEN 
       FILENOM='output_files/gridded_vturb'//ITER_STRING
       OPEN(99, FILE=FILENOM, STATUS='UNKNOWN', FORM='UNFORMATTED')
        WRITE(99) (((U2(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        WRITE(99) (((U3(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        WRITE(99) (((U4(I,J,K),I=1,NX),J=1,NY),K=1,NZ)
        LOW1=SUM(NPATCH(0:NL))
        DO I=1,LOW1
         N1=PATCHNX(I)
         N2=PATCHNY(I)
         N3=PATCHNZ(I)
         WRITE(99) (((U12(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
         WRITE(99) (((U13(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
         WRITE(99) (((U14(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        END DO
       CLOSE(99)
      END IF

      RETURN
      END

**********************************************************************
       SUBROUTINE WRITE_SHOCKED(NX,NY,NZ,ITER,NL,NPATCH,
     &            PATCHNX,PATCHNY,PATCHNZ,SHOCK0,SHOCK1)
***********************************************************************
*     Writes the filter length to a separate file.
***********************************************************************
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'

*     FUNCTION ARGUMENTS
      CHARACTER*200 FILERR5
      CHARACTER*5 ITER_STRING
      INTEGER NX, NY, NZ, NL, ITER
      INTEGER NPATCH(0:NLEVELS)
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER*1 SHOCK0(1:NMAX,1:NMAY,1:NMAZ)
      INTEGER*1 SHOCK1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)

*     VARIABLES
      INTEGER IR, I, LOW1, LOW2, IX, J, K, N1, N2, N3

*     OPEN THE OUTPUT FILE
      WRITE(ITER_STRING,'(I5.5)') ITER
      FILERR5='./output_files/shocked'//ITER_STRING
      
      OPEN(25,FILE=FILERR5,STATUS='UNKNOWN',FORM='UNFORMATTED')

*     WRITE THE 'COHERENCE' LENGTH

      WRITE(25) (((SHOCK0(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

      DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO I=LOW1,LOW2
          N1=PATCHNX(I)
          N2=PATCHNY(I)
          N3=PATCHNZ(I)
          WRITE(25) (((SHOCK1(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        END DO
      END DO

      CLOSE(25)

      END
#endif
#endif
