#if reader==0
***********************************************************************
      SUBROUTINE READ_GADGET_UNFORMATTED_NPART(ITER, FILES_PER_SNAP,
     &             FLAG_FILTER, FLAG_MACHFIELD)
***********************************************************************
*       Reads the number of gas particles in the simulation
***********************************************************************
      use gadget_read
      use particle_data  ! contains PARTI
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      CHARACTER*3 ITER_STRING
      CHARACTER*7 IFILE_STRING
      INTEGER IFILE
      CHARACTER*200 FIL1,FIL2

      integer npart_gadget(6), nall(6), blocksize
      real*8 massarr(6)
      integer basint1,basint2,basint3,basint4
      real*8 bas81,bas82,bas83,bas84,bas85,bas86

      parti=0
      do ifile=0,files_per_snap-1 
        WRITE(ITER_STRING,'(I3.3)') ITER
        FIL1='./simulation/snapdir_'//ITER_STRING//'/snap_'//ITER_STRING
        IF (FILES_PER_SNAP.EQ.1) THEN
          FIL2=FIL1
        ELSE
          WRITE(IFILE_STRING,'(I7.1)') IFILE
          FIL2=TRIM(ADJUSTL(FIL1))//'.'//TRIM(ADJUSTL(IFILE_STRING))
        END IF

        CALL read_head(FIL2,npart_gadget,massarr,bas81,bas82,
     &                  basint1,basint2,nall,basint3,basint4,bas83,
     &                  bas84,bas85,bas86,blocksize)

        parti=parti+npart_gadget(1)
      end do 

      return 
      end 
***********************************************************************

***********************************************************************
      SUBROUTINE READ_GADGET_UNFORMATTED(ITER, FILES_PER_SNAP,
     &             FLAG_FILTER, FLAG_MACHFIELD, LOW2)
***********************************************************************
*       Reads the GAS particle data of the simulation
***********************************************************************
      use gadget_read 
      use particle_data
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      CHARACTER*3 ITER_STRING
      CHARACTER*7 IFILE_STRING
      INTEGER IFILE
      CHARACTER*200 FIL1,FIL2

      integer npart_gadget(6), nall(6), blocksize
      real*8 massarr(6)
      integer basint1,basint2,basint3,basint4
      real*8 bas81,bas82,bas83,bas84,bas85,bas86
      integer low1,low2
      REAL*4,ALLOCATABLE::SCR4(:)
      REAL*4,ALLOCATABLE::SCR42(:,:)

      LOW2=0
      DO IFILE=0,FILES_PER_SNAP-1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*     READING DATA
      WRITE(ITER_STRING,'(I3.3)') ITER
      FIL1='./simulation/snapdir_'//ITER_STRING//'/snap_'//ITER_STRING
      IF (FILES_PER_SNAP.EQ.1) THEN
        FIL2=FIL1
      ELSE
        WRITE(IFILE_STRING,'(I7.1)') IFILE
        FIL2=TRIM(ADJUSTL(FIL1))//'.'//TRIM(ADJUSTL(IFILE_STRING))
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

      WRITE(*,*) 'Reading kernel length ...'
      CALL read_float(FIL2,'HSML',SCR4,blocksize)
      WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
      KERNEL(LOW1:LOW2)=SCR4(1:NPART_GADGET(1)) 
      DEALLOCATE(SCR4)       

#if use_filter == 1
      IF (FLAG_FILTER.EQ.1) THEN
        ALLOCATE(SCR4(SUM(NPART_GADGET(1:6))))
        IF (FLAG_MACHFIELD.EQ.0) THEN
          WRITE(*,*) 'Reading ABVC ...'
          CALL read_float(FIL2,'ABVC',SCR4,blocksize)
          WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
          ABVC(LOW1:LOW2)=SCR4(1:NPART_GADGET(1))  
        ELSE 
          WRITE(*,*) 'Reading MACH ...'
          CALL read_float(FIL2,'MACH',SCR4,blocksize)
          WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
          ABVC(LOW1:LOW2)=SCR4(1:NPART_GADGET(1)) 
        END IF            
        DEALLOCATE(SCR4)      
      END IF
#endif

#if weight_scheme == 2
      ALLOCATE(SCR4(SUM(NPART_GADGET(1:6))))
      WRITE(*,*) 'Reading density ...'
      CALL read_float(FIL2,'RHO ',SCR4,blocksize)
      WRITE(*,*) ' found for ',(blocksize-8)/4,' particles'
      VOL(LOW1:LOW2)=SCR4(1:NPART_GADGET(1))
      DEALLOCATE(SCR4)
#endif

      END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return 
      end 
***********************************************************************
#endif

#if reader==1 
***********************************************************************
      SUBROUTINE READ_AREPO_HDF5_NPART(ITER, FILES_PER_SNAP, 
     & FLAG_FILTER, FLAG_MACHFIELD)
***********************************************************************
*       Reads the number of gas particles in the simulation
***********************************************************************
        use hdf5
        use particle_data  ! contains PARTI
        implicit none 

        integer :: iter, files_per_snap, flag_filter, flag_machfield
        CHARACTER*3 ITER_STRING
        CHARACTER*1 IFILE_STRING
        INTEGER :: IFILE
        CHARACTER*200 FIL1, FIL2

        integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, 
     & file_space_id
        INTEGER(HID_T) :: memtype_id
        integer :: status
        integer, dimension(6) :: NumPart_ThisFile
        integer(hsize_t), dimension(1) :: dims = 6

        ! Default values for optional parameters
        integer(hid_t) :: xfer_prp
        !parameter (xfer_prp = H5P_DEFAULT_F)  ! Default transfer property list

        parti = 0
        do ifile = 0, files_per_snap - 1 
          WRITE(ITER_STRING, '(I3.3)') ITER
          FIL1 = './simulation/snapshot_' // ITER_STRING
          IF (FILES_PER_SNAP .EQ. 1) THEN
            FIL2 = FIL1
          ELSE
            WRITE(IFILE_STRING, '(I1.1)') IFILE
            FIL2 = TRIM(ADJUSTL(FIL1)) // '.' // IFILE_STRING
          END IF
          FIL2 = TRIM(ADJUSTL(FIL2)) // '.hdf5'

          ! Open the HDF5 file in read-only mode
          CALL h5fopen_f(FIL2, H5F_ACC_RDONLY_F, file_id, status)
          IF (status /= 0) THEN
            PRINT *, "Error opening file: ", FIL2
            STOP
          END IF

          ! Open the Header group
          CALL h5gopen_f(file_id, "/Header", group_id, status)
          IF (status /= 0) THEN
            PRINT *, "Error opening group: /Header"
            CALL h5fclose_f(file_id, status)
            STOP
          END IF
      
          ! Open the NumPart_ThisFile attribute
          CALL h5aopen_f(group_id, "NumPart_ThisFile",
     &                   attr_id, status)
          IF (status /= 0) THEN
            PRINT *, "Error opening dataset: NumPart_ThisFile"
            CALL h5gclose_f(group_id, status)
            CALL h5fclose_f(file_id, status)
            STOP
          END IF

          ! Get the datatype of the attribute to confirm it is H5T_NATIVE_INTEGER
          CALL h5aget_type_f(attr_id, memtype_id, status)
          IF (status /= 0) THEN
            PRINT *, "Error getting attribute type"
            CALL h5aclose_f(attr_id, status)
            CALL h5gclose_f(group_id, status)
            CALL h5fclose_f(file_id, status)
            STOP
          END IF
      
          ! Read the attribute into the NumPart_ThisFile array
          CALL h5aread_f(attr_id, memtype_id, NumPart_ThisFile,
     &                   dims, status)

          IF (status /= 0) THEN
            PRINT *, "Error reading attribute: NumPart_ThisFile"
            CALL h5aclose_f(attr_id, status)
            CALL h5gclose_f(group_id, status)
            CALL h5fclose_f(file_id, status)
            STOP
          END IF

          ! Print the contents of NumPart_ThisFile
          ! PRINT *, "NumPart_ThisFile: ", NumPart_ThisFile

          ! Close the dataset, group, and file
          CALL h5aclose_f(attr_id, status)
          CALL h5gclose_f(group_id, status)
          CALL h5fclose_f(file_id, status)

          ! We just read PartType0
          parti = parti + NumPart_ThisFile(1)
        end do 

        write(*,*) 'Total number of particles: ', parti
        RETURN
      END
***********************************************************************

***********************************************************************
      SUBROUTINE READ_AREPO_HDF5(ITER, FILES_PER_SNAP,
     &             FLAG_FILTER, FLAG_MACHFIELD, LOW2)
***********************************************************************
*       Reads the GAS particle data of the simulation
***********************************************************************
      use hdf5
      use particle_data
      implicit none 

      integer iter, files_per_snap, flag_filter, flag_machfield 

      CHARACTER*3 ITER_STRING
      CHARACTER*1 IFILE_STRING
      INTEGER IFILE
      CHARACTER*200 FIL1,FIL2

      integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, 
     & file_space_id
      INTEGER(HID_T) :: memtype_id
      integer :: status
      integer, dimension(6) :: NumPart_ThisFile
      integer(hsize_t), dimension(1) :: dims1d
      integer(hsize_t), dimension(2) :: dims2d
      
      integer basint1,basint2,basint3,basint4,i
      real*8 bas81,bas82,bas83,bas84,bas85,bas86
      integer low1,low2
      REAL*4,ALLOCATABLE::SCR4(:)
      REAL*4,ALLOCATABLE::SCR42(:,:)

      LOW2=0
      DO IFILE=0,FILES_PER_SNAP-1 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*     READING DATA
      WRITE(ITER_STRING, '(I3.3)') ITER
      FIL1 = './simulation/snapshot_' // ITER_STRING
      IF (FILES_PER_SNAP .EQ. 1) THEN
        FIL2 = FIL1
      ELSE
        WRITE(IFILE_STRING, '(I1.1)') IFILE
        FIL2 = TRIM(ADJUSTL(FIL1)) // '.' // IFILE_STRING
      END IF
      FIL2 = TRIM(ADJUSTL(FIL2)) // '.hdf5'

      ! Open the HDF5 file in read-only mode
      WRITE(*,*) 'Reading iteration file: ',ITER,' ',
     &              TRIM(ADJUSTL(FIL2))

      CALL h5fopen_f(FIL2, H5F_ACC_RDONLY_F, file_id, status)
      IF (status /= 0) THEN
        PRINT *, "Error opening file: ", FIL2
        STOP
      END IF

      dims1d(1) = 6
      CALL h5gopen_f(file_id, "/Header", group_id, status)
      CALL h5aopen_f(group_id, "NumPart_ThisFile",
     &                   attr_id, status)
      CALL h5aget_type_f(attr_id, memtype_id, status)
      CALL h5aread_f(attr_id, memtype_id, NumPart_ThisFile,
     &                   dims1d, status)

      WRITE(*,*) NumPart_ThisFile(1), 'gas particles'
      LOW1=LOW2+1
      LOW2=LOW1+NumPart_ThisFile(1)-1

      CALL h5aclose_f(attr_id, status)
      CALL h5gclose_f(group_id, status)

      dims1d(1) = NumPart_ThisFile(1)
      dims2d(1) = NumPart_ThisFile(1)
      dims2d(2) = 3

      CALL h5gopen_f(file_id, "/PartType0",
     &                   group_id, status)
      if (status /= 0) then
        PRINT *, "Error opening group: /PartType0"
        CALL h5fclose_f(file_id, status)
        STOP
      end if
      
      ALLOCATE(SCR42(3,NumPart_ThisFile(1)))

      WRITE(*,*) 'Reading positions ...'
      CALL h5dopen_f(group_id, "Coordinates", attr_id, status)
      CALL h5dget_type_f(attr_id, memtype_id, status)
      CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
      RXPA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(1))
      RYPA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(1))
      RZPA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(1))
      CALL h5dclose_f(attr_id, status)

      WRITE(*,*) 'Reading velocities ...'
      CALL h5dopen_f(group_id, "Velocities", attr_id, status)
      CALL h5dget_type_f(attr_id, memtype_id, status)
      CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
      U2DM(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(1))
      U3DM(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(1))
      U4DM(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(1))
      CALL h5dclose_f(attr_id, status)

      DEALLOCATE(SCR42)

      ALLOCATE(SCR4(NumPart_ThisFile(1)))
      WRITE(*,*) 'Reading masses ...'
      CALL h5dopen_f(group_id, "Masses", attr_id, status)
      CALL h5dget_type_f(attr_id, memtype_id, status)
      CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
      MASAP(LOW1:LOW2)=SCR4(1:NumPart_ThisFile(1))
      CALL h5dclose_f(attr_id, status)

      WRITE(*,*) 'Reading effective kernel length ...' ! We read volume
      CALL h5dopen_f(group_id, "Volume", attr_id, status)
      CALL h5dget_type_f(attr_id, memtype_id, status)
      CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
!$omp parallel do shared(numPart_ThisFile, scr4, kernel, low1), 
!$omp+            private(i), default(none)
      do i=1,NumPart_ThisFile(1)
        kernel(low1+i-1) = (0.2387*SCR4(i))**0.33333
      end do
      CALL h5dclose_f(attr_id, status)
      DEALLOCATE(SCR4)

#if use_filter == 1
      IF (FLAG_FILTER.EQ.1) THEN
        ALLOCATE(SCR4(NumPart_ThisFile(1)))
        IF (FLAG_MACHFIELD.EQ.0) THEN
          WRITE(*,*) 'WARNING! In AREPO, we always read Mach!'
        END IF 
        WRITE(*,*) 'Reading MACH ...'
        CALL h5dopen_f(group_id, "Machnumber", attr_id, status)
        CALL h5dget_type_f(attr_id, memtype_id, status)
        CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
        ABVC(LOW1:LOW2)=SCR4(1:NumPart_ThisFile(1))  
        CALL h5dclose_f(attr_id, status)         
        DEALLOCATE(SCR4)      
      END IF
#endif

#if weight_scheme == 2
      ALLOCATE(SCR4(NumPart_ThisFile(1)))
      WRITE(*,*) 'Reading density ...'
      CALL h5dopen_f(group_id, "Density", attr_id, status)
      CALL h5dget_type_f(attr_id, memtype_id, status)
      CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
      VOL(LOW1:LOW2)=SCR4(1:NumPart_ThisFile(1))
      CALL h5dclose_f(attr_id, status)
      DEALLOCATE(SCR4)
#endif
        
      CALL h5gclose_f(group_id, status)
      CALL h5fclose_f(file_id, status)

      END DO !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return 
      end 
***********************************************************************
#endif


***********************************************************************
       SUBROUTINE READ_PARTICLES(ITER,FILES_PER_SNAP,NX,NY,NZ,T,ZETA,
     &            NL_PARTICLE_GRID,REFINE_THR,PARCHLIM,BORGRID,
     &            NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &            PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ,LADO0,
     &            NPART,
     &            FLAG_FILTER,KNEIGHBOURS,DIV_THR,ABVC_THR,
     &            FLAG_MACHFIELD,MACH_THR,FLAG_MASS)
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

       USE particle_data
       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'

       INTEGER NX,NY,NZ,ITER,NDXYZ,LOW1,LOW2,FILES_PER_SNAP
       real T,AAA,BBB,CCC,MAP,ZETA,LADO0
       INTEGER I,J,K,IX,NL,IR,IRR,N1,N2,N3,NL_PARTICLE_GRID
       INTEGER REFINE_THR,PARCHLIM,BORGRID,KNEIGHBOURS
       REAL DIV_THR,ABVC_THR
       INTEGER FLAG_MACHFIELD,FLAG_MASS
       REAL MACH_THR

       INTEGER FLAG_FILTER

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

#if use_filter == 1
       INTEGER*1 SHOCK0(1:NMAX,1:NMAY,1:NMAZ)
       INTEGER*1 SHOCK1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       COMMON /SHOCKED/ SHOCK0,SHOCK1
#endif 

       integer cr0amr(1:NMAX,1:NMAY,1:NMAZ)
       integer cr0amr1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
       common /cr0/ cr0amr, cr0amr1

*      ---PARALLEL---
       INTEGER NUM,OMP_GET_NUM_THREADS,NUMOR, FLAG_PARALLEL
       COMMON /PROCESADORES/ NUM

       REAL DDXL,DDXR,DDYL,DDYR,DDZL,DDZR
       REAL CIO_XC,CIO_YC,CIO_ZC
       COMMON /DOM_DECOMP/ DDXL,DDXR,DDYL,DDYR,DDZL,DDZR

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

       CHARACTER*3 ITER_STRING
       CHARACTER*1 IFILE_STRING
       INTEGER IFILE

       ! Scratch variables for restricting the domain
       REAL*4,ALLOCATABLE::SCR42(:,:)
       INTEGER,ALLOCATABLE::ELIM(:)
       ! End scratch variables for restricting the domain

#if use_filter == 1
       REAL VISC0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
       REAL VISC1(NAMRX,NAMRY,NAMRZ,NPALEV)
#else 
       ! Dummy variables
       REAL VISC0, VISC1 
#endif


       real xmin,ymin,zmin,xmax,ymax,zmax

************************************************************************
**************** CHANGE FOR OTHER SIMULATION CODES *********************
************************************************************************
       ! First, get the number of particles to be read in the snapshot 
#if reader == 0
       CALL READ_GADGET_UNFORMATTED_NPART(ITER, FILES_PER_SNAP,
     &                                    FLAG_FILTER, FLAG_MACHFIELD)
#endif
#if reader == 1
       CALL READ_AREPO_HDF5_NPART(ITER, FILES_PER_SNAP,
     &                                    FLAG_FILTER, FLAG_MACHFIELD)
#endif
************************************************************************
************************************************************************
************************************************************************

       ! Deallocate the particle arrays if they were allocated
       if (allocated(rxpa)) then 
          deallocate(rxpa,rypa,rzpa)
          deallocate(u2dm,u3dm,u4dm)
          deallocate(masap,kernel)
#if use_filter == 1
          deallocate(abvc)
#endif

#if weight_scheme == 2
          deallocate(vol)
#endif
        end if

       ! Allocate the particle arrays
       allocate(rxpa(parti),rypa(parti),rzpa(parti))
       allocate(u2dm(parti),u3dm(parti),u4dm(parti))
       allocate(masap(parti),kernel(parti))
#if use_filter == 1
       allocate(abvc(parti))
#endif

#if weight_scheme == 2
        allocate(vol(parti))
#endif

       NPART(:)=0

************************************************************************
**************** CHANGE FOR OTHER SIMULATION CODES *********************
************************************************************************
       ! First, get the number of particles to be read in the snapshot 
#if reader == 0
       CALL READ_GADGET_UNFORMATTED(ITER, FILES_PER_SNAP,
     &                        FLAG_FILTER, FLAG_MACHFIELD, LOW2)
#endif
#if reader == 1
       CALL READ_AREPO_HDF5(ITER, FILES_PER_SNAP,
     &                        FLAG_FILTER, FLAG_MACHFIELD, LOW2)
#endif
************************************************************************
************************************************************************
************************************************************************

#if weight_scheme == 2
!$OMP PARALLEL DO SHARED(VOL, MASAP, PARTI), PRIVATE(I), DEFAULT(NONE)
       DO I=1,PARTI
         VOL(I)=MASAP(I)/VOL(I)
       END DO
#endif 
       
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
      
#if use_filter == 1
      IF (FLAG_FILTER.EQ.1) THEN
        IF (FLAG_MACHFIELD.EQ.0) THEN
         WRITE(*,*) 'ABVC=',MINVAL(ABVC(LOW1:LOW2)),
     &                      MAXVAL(ABVC(LOW1:LOW2))
        ELSE 
         WRITE(*,*) 'MACH=',MINVAL(ABVC(LOW1:LOW2)),
     &                      MAXVAL(ABVC(LOW1:LOW2))
        END IF
       END IF
#endif

       IF (XMIN.LT.DDXL.OR.XMAX.GT.DDXR.OR.
     &     YMIN.LT.DDYL.OR.YMAX.GT.DDYR.OR.
     &     ZMIN.LT.DDZL.OR.ZMAX.GT.DDZR) THEN

        WRITE(*,*)
        WRITE(*,*) 'WARNING: PARTICLES OUTSIDE THE DOMAIN'
        ALLOCATE(ELIM(LOW1:LOW2))
        
        J=0
!$OMP PARALLEL DO SHARED(RXPA,RYPA,RZPA,LOW1,LOW2,DDXL,DDXR,DDYL,DDYR,
!$OMP+            DDZL,DDZR,ELIM), 
!$OMP+            PRIVATE(I), 
!$OMP+            REDUCTION(+: J),
!$OMP+            DEFAULT(NONE)
        DO I=LOW1,LOW2 
          ELIM(I)=0
          IF (RXPA(I).LT.DDXL.OR.RXPA(I).GT.DDXR.OR.
     &        RYPA(I).LT.DDYL.OR.RYPA(I).GT.DDYR.OR.
     &        RZPA(I).LT.DDZL.OR.RZPA(I).GT.DDZR) THEN
            J=J+1
            ELIM(I)=1
          END IF
        END DO

        ! NEW NUMBER OF PARTICLES GETS UPDATED HERE
        PARTI=PARTI-J
        
        WRITE(*,*) 'Particles outside the domain:',J
        ALLOCATE(SCR42(9,PARTI))

        J=0
        DO I=LOW1,LOW2 
          IF (ELIM(I).EQ.0) THEN
            J=J+1
            SCR42(1,J)=RXPA(I)
            SCR42(2,J)=RYPA(I)
            SCR42(3,J)=RZPA(I)
            SCR42(4,J)=U2DM(I)
            SCR42(5,J)=U3DM(I)
            SCR42(6,J)=U4DM(I)
            SCR42(7,J)=MASAP(I)
            SCR42(8,J)=KERNEL(I)
#if use_filter == 1
            IF (FLAG_FILTER.EQ.1) SCR42(9,J)=ABVC(I)
#endif
          END IF
        END DO
        DEALLOCATE(ELIM)

        DEALLOCATE(RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,KERNEL)
#if use_filter == 1
        DEALLOCATE(ABVC)
#endif

        ALLOCATE(RXPA(PARTI),RYPA(PARTI),RZPA(PARTI))
        ALLOCATE(U2DM(PARTI),U3DM(PARTI),U4DM(PARTI))
        ALLOCATE(MASAP(PARTI),KERNEL(PARTI))
#if use_filter == 1
        ALLOCATE(ABVC(PARTI))
#endif

!$OMP PARALLEL DO SHARED(SCR42,RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,
!$OMP+                   KERNEL,ABVC,PARTI,FLAG_FILTER), 
!$OMP+            PRIVATE(I), 
!$OMP+            DEFAULT(NONE)
        DO I=1,PARTI 
          RXPA(I)=SCR42(1,I)
          RYPA(I)=SCR42(2,I)
          RZPA(I)=SCR42(3,I)
          U2DM(I)=SCR42(4,I)
          U3DM(I)=SCR42(5,I)
          U4DM(I)=SCR42(6,I)
          MASAP(I)=SCR42(7,I)
          KERNEL(I)=SCR42(8,I)
#if use_filter == 1
          IF (FLAG_FILTER.EQ.1) ABVC(I)=SCR42(9,I)
#endif
        END DO

        DEALLOCATE(SCR42)

        LOW2=PARTI
        NPART(0)=PARTI !Correct the number of particles!

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
       
#if use_filter == 1
       IF (FLAG_FILTER.EQ.1) THEN
        IF (FLAG_MACHFIELD.EQ.0) THEN
         WRITE(*,*) 'ABVC=',MINVAL(ABVC(LOW1:LOW2)),
     &                      MAXVAL(ABVC(LOW1:LOW2))
        ELSE 
         WRITE(*,*) 'MACH=',MINVAL(ABVC(LOW1:LOW2)),
     &                      MAXVAL(ABVC(LOW1:LOW2))
        END IF
       END IF
#endif
       

       WRITE(*,*) 'Routine create mesh ------------------------------'
       NPATCH(0:IR)=0
       if (parchlim.gt.0) then
        CALL CREATE_MESH(NX,NY,NZ,NL_PARTICLE_GRID,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,REFINE_THR,PARCHLIM,BORGRID)
       else
        CALL CREATE_MESH_OCTREE(NX,NY,NZ,NL_PARTICLE_GRID,NPATCH,
     &            PARE,PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,REFINE_THR,PARCHLIM,BORGRID)
       end if
       
       NL=NL_PARTICLE_GRID
       DO IR=1,NL_PARTICLE_GRID
        IF (NPATCH(IR).EQ.0) THEN 
          NL=IR-1
          EXIT
        END IF
       END DO

       CALL GRIDAMR(NX,NY,NZ,NL,NPATCH,
     &                   PATCHNX,PATCHNY,PATCHNZ,
     &                   PATCHX,PATCHY,PATCHZ,
     &                   PATCHRX,PATCHRY,PATCHRZ,PARE)
       WRITE(*,*) 'End mesh creation --------------------------------'

       WRITE(*,*) 'Routine interpolate velocity ---------------------'
       CALL INTERPOLATE_VELOCITIES(NX,NY,NZ,NL,NPATCH,PARE,
     &            PATCHNX,PATCHNY,PATCHNZ,PATCHX,PATCHY,PATCHZ,
     &            PATCHRX,PATCHRY,PATCHRZ,
     &            NPART,LADO0,FLAG_FILTER,KNEIGHBOURS,
     &            VISC0,VISC1,FLAG_MACHFIELD,FLAG_MASS)

#if output_particles == 1
       IF (FL_P_ERR.EQ.1) THEN
        WRITE(*,*) 'Locating particles onto the grid'
        CALL PLACE_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &             PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &             NPART,LADO0,parchlim)

        CALL ERROR_PARTICLES(NX,NY,NZ,NL,NPATCH,PATCHNX,PATCHNY,
     &             PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,PARE,
     &             NPART,LADO0)
        DEALLOCATE(LIHAL, LIHAL_IX, LIHAL_JY, LIHAL_KZ)
       END IF
#endif
       WRITE(*,*) 'End velocity interpolation -----------------------'

#if use_filter == 1
       IF (FLAG_FILTER.EQ.1) THEN 
*     All patches are extended with one extra cell per direction
        CALL EXTEND_VAR(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                  PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

        CALL IDENTIFY_SHOCKS(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,
     &                       PATCHNY,PATCHNZ,PATCHRX,PATCHRY,PATCHRZ,
     &                       PATCHX,PATCHY,PATCHZ,LADO0,VISC0,VISC1,
     &                       DIV_THR,ABVC_THR,FLAG_MACHFIELD,MACH_THR)
       END IF
#endif

!      If we do not want to output the particles, we deallocate them
!      here
#if output_particles == 0
        DEALLOCATE(RXPA,RYPA,RZPA,U2DM,U3DM,U4DM,MASAP,KERNEL)
#if use_filter == 1
        DEALLOCATE(ABVC)
#endif
#endif


       RETURN
       END

#if use_filter == 1
***********************************************************************
       SUBROUTINE IDENTIFY_SHOCKS(ITER,NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,
     &                            PATCHNY,PATCHNZ,PATCHRX,PATCHRY,
     &                            PATCHRZ,PATCHX,PATCHY,PATCHZ,LADO0,
     &                            VISC0,VISC1,DIV_THR,ABVC_THR,
     &                            FLAG_MACHFIELD,MACH_THR)
***********************************************************************
*      For the multiscale filter, an indication of (strong) shocked
*       cells is required. This routine does that job by using a
*       combination of velocity divergence and artificial bulk 
*       viscosity constant.
***********************************************************************

      IMPLICIT NONE 
      INCLUDE 'vortex_parameters.dat'
      INTEGER ITER,NX,NY,NZ,NL 
      INTEGER NPATCH(0:NLEVELS),PARE(NPALEV),PATCHNX(NPALEV),
     &        PATCHNY(NPALEV),PATCHNZ(NPALEV)
      REAL PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV)
      INTEGER PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV)
      REAL LADO0
      REAL VISC0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL VISC1(NAMRX,NAMRY,NAMRZ,NPALEV)
      REAL DIV_THR,ABVC_THR
      INTEGER FLAG_MACHFIELD
      REAL MACH_THR

      REAL U2(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U3(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U4(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL U12(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U13(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      REAL U14(0:NAMRX+1,0:NAMRY+1,0:NAMRZ+1,NPALEV)
      COMMON /VELOC/ U2,U3,U4,U12,U13,U14

      REAL DIVER0(0:NMAX+1,0:NMAY+1,0:NMAZ+1)
      REAL DIVER(-2:NAMRX+3,-2:NAMRY+3,-2:NAMRZ+3,NPALEV)
      COMMON /DIVERGENCE/ DIVER0, DIVER

      INTEGER*1 SHOCK0(1:NMAX,1:NMAY,1:NMAZ)
      INTEGER*1 SHOCK1(1:NAMRX,1:NAMRY,1:NAMRZ,NPALEV)
      COMMON /SHOCKED/ SHOCK0,SHOCK1

      INTEGER I,J,K,IPATCH,N1,N2,N3,LOW1,LOW2,IR

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

      character*5 iter_string 
      write(iter_string, '(I5.5)') iter

      IF (FLAG_MACHFIELD.EQ.0) THEN
        WRITE(*,*) 'Identifying shocks with DIV_THR=',DIV_THR,
     &             ' and ABVC_THR=',ABVC_THR
       CALL DIVER_FINA(NX,NY,NZ,NL,NPATCH,PARE,PATCHNX,PATCHNY,PATCHNZ,
     &                 PATCHX,PATCHY,PATCHZ,PATCHRX,PATCHRY,PATCHRZ)

C       open(55,file='output_files/diver_unfiltered_'//iter_string,
C     &      status='unknown', form='unformatted')
C
C        write(55) (((diver0(i,j,k),i=1,nx),j=1,ny),k=1,nz)
C        do ipatch=1,sum(npatch(0:nl))
C         n1=patchnx(ipatch)
C         n2=patchny(ipatch)
C         n3=patchnz(ipatch)
C         write(55) (((diver(i,j,k,ipatch),i=1,n1),j=1,n2),k=1,n3)
C        end do
C
C       close(55)

       DO K=1,NZ 
       DO J=1,NY
       DO I=1,NX 
         SHOCK0(I,J,K) = 0
         IF (DIVER0(I,J,K).LT.DIV_THR.AND.VISC0(I,J,K).GT.ABVC_THR) THEN
           SHOCK0(I,J,K) = 1
         END IF
       END DO 
       END DO 
       END DO

       DO IR=1,NL
        LOW1=SUM(NPATCH(0:IR-1))+1
        LOW2=SUM(NPATCH(0:IR))
        DO IPATCH=LOW1,LOW2 
          N1=PATCHNX(IPATCH)
          N2=PATCHNY(IPATCH)
          N3=PATCHNZ(IPATCH)
          DO K=1,N3 
          DO J=1,N2
          DO I=1,N1 
            SHOCK1(I,J,K,IPATCH) = 0
            IF (DIVER(I,J,K,IPATCH).LT.DIV_THR.AND.
     &          VISC1(I,J,K,IPATCH).GT.ABVC_THR) THEN
              SHOCK1(I,J,K,IPATCH) = 1
            END IF
          END DO 
          END DO 
          END DO
        END DO
       END DO

      ELSE 
        WRITE(*,*) 'Identifying shocks with MACH_THR=',MACH_THR

        DO K=1,NZ 
        DO J=1,NY
        DO I=1,NX 
          SHOCK0(I,J,K) = 0
          IF (VISC0(I,J,K).GT.MACH_THR) THEN
            SHOCK0(I,J,K) = 1
          END IF
        END DO 
        END DO 
        END DO
  
        DO IR=1,NL
          LOW1=SUM(NPATCH(0:IR-1))+1
          LOW2=SUM(NPATCH(0:IR))
          DO IPATCH=LOW1,LOW2 
            N1=PATCHNX(IPATCH)
            N2=PATCHNY(IPATCH)
            N3=PATCHNZ(IPATCH)
            DO K=1,N3 
            DO J=1,N2
            DO I=1,N1 
              SHOCK1(I,J,K,IPATCH) = 0
              IF (VISC1(I,J,K,IPATCH).GT.MACH_THR) THEN
                SHOCK1(I,J,K,IPATCH) = 1
              END IF
            END DO 
            END DO 
            END DO
          END DO
        END DO
      
      END IF

#if output_filter == 1
      IF (FL_FILT_SHOCK.EQ.1) THEN
       CALL WRITE_SHOCKED(NX,NY,NZ,ITER,NL,NPATCH,PATCHNX,PATCHNY,
     &                    PATCHNZ,SHOCK0,SHOCK1)
      END IF
#endif

      RETURN 
      END 
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

#if reader == 2 
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
#elif weight_filter == 0 
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
        read(31) ! temp 
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
#elif weight_filter == 0
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
     

#endif 