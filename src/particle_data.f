      MODULE PARTICLE_DATA

        IMPLICIT NONE
        SAVE
        INTEGER PARTI

        REAL,ALLOCATABLE::RXPA(:),RYPA(:),RZPA(:)
        REAL,ALLOCATABLE::U2DM(:),U3DM(:),U4DM(:)
        REAL,ALLOCATABLE::MASAP(:),KERNEL(:)

#if use_filter == 1
        REAL,ALLOCATABLE::ABVC(:)
#else 
        ! Dummy variable 
        REAL ABVC
#endif

#if weight_scheme == 2
        REAL,ALLOCATABLE::VOL(:)
#else
        ! Dummy variable 
        REAL VOL
#endif 

        INTEGER,ALLOCATABLE::LIHAL(:)
        INTEGER,ALLOCATABLE::LIHAL_IX(:),LIHAL_JY(:),LIHAL_KZ(:)

      END MODULE PARTICLE_DATA
