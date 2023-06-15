      MODULE gadget_read

      INTEGER(kind=4),SAVE ::  lun
      
      CONTAINS
      
      SUBROUTINE find_block(filename,label,blocksize)
      
        IMPLICIT NONE
          
      !...dummy arguments
        CHARACTER(len=*), INTENT(in)     :: filename
        CHARACTER(len=4), INTENT(in)     :: label
        INTEGER(kind=4),  INTENT(out)    :: blocksize
      
        INTEGER(kind=4) ::  istat
        REAL(kind=4) ::  dummy
        CHARACTER(len=4) :: blocklabel
      
      
      !...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
      !...try to open file in little endian
          lun = 1
      
          OPEN (lun, FILE=filename,STATUS='old',FORM='unformatted',
     &          iostat=istat)
      
          READ (lun,iostat=istat) blocklabel,blocksize
      
          IF ((istat /= 0).OR.(blocksize /= 264)) THEN
            CLOSE(lun)
      
      !...try to open file in big endian
            lun = 1001
      
            OPEN (lun, FILE=filename,STATUS='old',FORM='unformatted',
     &            iostat=istat)
      
            READ (lun,iostat=istat) blocklabel,blocksize
      
            IF ((istat /= 0).OR.(blocksize /= 264)) THEN
              PRINT *,'  snapshot corrupted !',istat,blocklabel,
     &                blocksize
              PRINT *,'  I have to stop here ...'
              stop
            ELSE
              !PRINT *, ' snapshot in big endian format'
            ENDIF
          ELSE
            !PRINT *, ' snapshot in little endian format'
      
          ENDIF
      
          DO WHILE (TRIM(label) /= TRIM(blocklabel))
              READ (lun,iostat=istat) dummy
      
              READ (lun,iostat=istat) blocklabel,blocksize
      
              IF (istat /= 0) THEN
                blocklabel=label
                blocksize=-1
              ENDIF
      
          ENDDO
        
      END SUBROUTINE find_block
      
      
      SUBROUTINE read_float3(filename,label,x,blocksize)
        CHARACTER(len=*), INTENT(in)     :: filename
        CHARACTER(len=4), INTENT(in)     :: label
        REAL(kind=4),     INTENT(INOUT)  :: x(:,:)
        INTEGER(kind=4),  INTENT(INOUT)  :: blocksize
        INTEGER(kind=4) ::  istat,i,j
      
        CALL find_block(filename,label,blocksize)
        IF (blocksize > 0) READ (lun,iostat=istat) ((x(i,j),i=1,3),j=1,
     &                                              (blocksize-8)/12)
        CLOSE (lun)
      END SUBROUTINE read_float3 
      
      
      SUBROUTINE read_float(filename,label,x,blocksize)
        CHARACTER(len=*), INTENT(in)     :: filename
        CHARACTER(len=4), INTENT(in)     :: label
        REAL(kind=4),     INTENT(INOUT)  :: x(:)
        INTEGER(kind=4),  INTENT(INOUT)  :: blocksize
        INTEGER(kind=4) ::  istat,i
      
        CALL find_block(filename,label,blocksize)
        IF (blocksize > 0) READ (lun,iostat=istat) (x(i),i=1,
     &                                              (blocksize-8)/4)
        CLOSE (lun)
      END SUBROUTINE read_float 
      
      SUBROUTINE read_long(filename,label,x,blocksize)
        CHARACTER(len=*), INTENT(in)     :: filename
        CHARACTER(len=4), INTENT(in)     :: label
        INTEGER(kind=4),  INTENT(INOUT)  :: x(:)
        INTEGER(kind=4),  INTENT(INOUT)  :: blocksize
        INTEGER(kind=4) ::  istat,i
      
        CALL find_block(filename,label,blocksize)
        IF (blocksize > 0) READ (lun,iostat=istat) (x(i),i=1,
     &                                              (blocksize-8)/4)
        CLOSE (lun)
      END SUBROUTINE read_long 
      
      
      SUBROUTINE read_head(filename,npart,massarr,a,redshift,flag_sfr,
     &                     flag_feedback,nall,cooling_flag,numfiles,
     &                     boxsize,Omega0,OmegaL0,Hubblepar,blocksize)
        CHARACTER(len=*), INTENT(in)     :: filename
        REAL(kind=8),     INTENT(out)    :: massarr(6)
        REAL(kind=8),     INTENT(out)    :: a,redshift
        INTEGER(kind=4),  INTENT(out)    :: flag_sfr,flag_feedback,
     &                                      cooling_flag
        INTEGER(kind=4),  INTENT(out)    :: npart(6),nall(6),numfiles
        REAL(kind=8),     INTENT(out)    :: boxsize,Omega0,OmegaL0,
     &                                      Hubblepar
        INTEGER(kind=4),  INTENT(INOUT)  :: blocksize
      
        CALL find_block(filename,'HEAD',blocksize)
        IF (blocksize > 0) THEN 
          READ (lun,iostat=istat) npart,massarr,a,redshift,flag_sfr,
     &                            flag_feedback,nall,cooling_flag,
     &                            numfiles,boxsize,Omega0,OmegaL0,
     &                            Hubblepar
        ELSE 
         WRITE(*,*) 'Block',filename,'not found'
        END IF
        CLOSE (lun)
      END SUBROUTINE read_head 
      
      END MODULE gadget_read