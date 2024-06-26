************************************************************************
      subroutine select(k,arr_in,n,outval) ! DEPECRATED, UNUSED
* adapted from
* Numerical Recipes in Fortran 90, Vol.2
* by Press, Teukoslky et al.
************************************************************************
      !USE nrtype; USE nrutil, ONLY : assert,swap
      IMPLICIT NONE
      !INTEGER, INTENT(IN) :: k
      !REAL, DIMENSION(:), INTENT(INOUT) :: arr
      !REAL :: select
      INTEGER, INTENT(IN) :: k
      integer, intent(in) :: n
      REAL, DIMENSION(n), INTENT(IN) :: arr_in
      real, intent(out) :: outval
*     Returns the kth smallest value in the array arr. The input array
*     will be rearranged to have this value in location arr(k), with
*     all smaller elements moved to arr(1:k-1) (in arbitrary order)
*     and all larger elements in arr(k+1:) (also in arbitrary order).
      INTEGER :: i,r,j,l
      REAL :: a
      REAL, DIMENSION(size(arr_in)) :: arr
      arr=arr_in
      !call assert(k >= 1, k <= n, ’select args’)
      l=1
      r=n
      do
       if (r-l <= 1) then
        if (r-l == 1) call masked_swap_rs(arr(l),arr(r),arr(l)>arr(r))
        outval=arr(k)
        RETURN
       else
        i=(l+r)/2
        call swap_r(arr(i),arr(l+1))
        call masked_swap_rs(arr(l),arr(r),arr(l)>arr(r))
        call masked_swap_rs(arr(l+1),arr(r),arr(l+1)>arr(r))
        call masked_swap_rs(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1
        j=r
        a=arr(l+1)
        do
         do
          i=i+1
          if (arr(i) >= a) exit
         end do
         do
          j=j-1
          if (arr(j) <= a) exit
         end do
         if (j < i) exit
         call swap_r(arr(i),arr(j))
        end do
        arr(l+1)=arr(j)
        arr(j)=a
        if (j >= k) r=j-1
        if (j <= k) l=i
       end if
      end do
      END subroutine

*** swap
      SUBROUTINE swap_r(a,b)
      REAL, INTENT(INOUT) :: a,b
      REAL :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_r

      SUBROUTINE masked_swap_rs(a,b,mask)
      REAL, INTENT(INOUT) :: a,b
      LOGICAL, INTENT(IN) :: mask
      REAL :: swp
      if (mask) then
       swp=a
       a=b
       b=swp
      end if
      END SUBROUTINE masked_swap_rs

**********************************************************************
      SUBROUTINE ARGSORT(N,ARRAY,INDX)
**********************************************************************
*      Return the indices which make ARRAY sorted in increasing order.
*      Uses MERGE SORT
*      Rewritten from the idea of the algorithm (LGPL-3.0 License) in
*       https://github.com/Astrokiwi/simple_fortran_argsort/blob/master/LICENSE
**********************************************************************
      INTEGER N
      REAL*4 ARRAY(N)
      INTEGER INDX(N)

      INTEGER TEMPINT(N)

      INTEGER DELTA
      INTEGER I,J,III,K,KKK

      DO I=1,N
       INDX(I)=I
      END DO

      IF (N.EQ.1) RETURN

      DELTA=1
      DO WHILE (DELTA.LT.N)
       DO III=1,N-DELTA,DELTA*2
        I1=III
        I2=III+DELTA
        KKK=MIN(DELTA*2,N-III+1)
        K=1

        DO WHILE (I1.LT.III+DELTA.AND.I2.LT.III+KKK)
         IF (ARRAY(INDX(I1)).LT.ARRAY(INDX(I2))) THEN
          TEMPINT(K)=INDX(I1)
          I1=I1+1
          K=K+1
         ELSE
          TEMPINT(K)=INDX(I2)
          I2=I2+1
          K=K+1
         END IF
        END DO

        IF (I1.LT.III+DELTA) THEN
         TEMPINT(K:KKK)=INDX(I1:III+DELTA-1)
        ELSE
         TEMPINT(K:KKK)=INDX(I2:III+KKK-1)
        ENDIF
        INDX(III:III+KKK-1)=TEMPINT(1:KKK)
       END DO
       DELTA=DELTA*2
      END DO

      RETURN
      END SUBROUTINE      

***************************
*     PARALLEL MAX AND MIN
***************************

      SUBROUTINE P_MINMAX(ARRAY0,ARRAY1,BOR,BORAMR,NX,NY,NZ,NL,
     &                    PATCHNX,PATCHNY,PATCHNZ,NPATCH)
       IMPLICIT NONE

       INCLUDE 'vortex_parameters.dat'
       
       INTEGER NX,NY,NZ,BOR,BORAMR,NL
       INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
       INTEGER NPATCH(0:NLEVELS)
       REAL ARRAY0(1-BOR:NMAX+BOR,1-BOR:NMAY+BOR,1-BOR:NMAZ+BOR)
       REAL ARRAY1(1-BORAMR:NAMRX+BORAMR,1-BORAMR:NAMRY+BORAMR,
     &             1-BORAMR:NAMRZ+BORAMR,NPALEV)

       REAL VMIN0,VMIN1,VMAX0,VMAX1
       INTEGER N1,N2,N3,I,LOW1

       VMIN0=1.E30
       VMAX0=-1.E30
!$OMP PARALLEL DO SHARED(NX,NY,NZ,ARRAY0),
!$OMP+            PRIVATE(I),
!$OMP+            REDUCTION(MIN:VMIN0), REDUCTION(MAX:VMAX0),
!$OMP+            DEFAULT(NONE)       
       DO I=1,NX
         VMIN0=MIN(MINVAL(ARRAY0(I,1:NY,1:NZ)),VMIN0)
         VMAX0=MAX(MAXVAL(ARRAY0(I,1:NY,1:NZ)),VMAX0)
       END DO


       VMIN1=1.E30
       VMAX1=-1.E30
       LOW1=SUM(NPATCH(0:NL))
!$OMP PARALLEL DO SHARED(NPATCH,PATCHNX,PATCHNY,PATCHNZ,NL,
!$OMP+                   ARRAY1,LOW1),
!$OMP+            PRIVATE(I,N1,N2,N3),
!$OMP+            REDUCTION(MIN:VMIN1), REDUCTION(MAX:VMAX1),
!$OMP+            DEFAULT(NONE)
       DO I=1,LOW1
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        VMIN1=MIN(MINVAL(ARRAY1(1:N1,1:N2,1:N3,I)),VMIN1)
        VMAX1=MAX(MAXVAL(ARRAY1(1:N1,1:N2,1:N3,I)),VMAX1)
       END DO

       WRITE(*,*) VMIN0,VMIN1
       WRITE(*,*) VMAX0,VMAX1

      RETURN
      END SUBROUTINE


      SUBROUTINE P_MINMAX_IR(ARRAY0,ARRAY1,BOR,BORAMR,NX,NY,NZ,NL,
     &                    PATCHNX,PATCHNY,PATCHNZ,NPATCH,IR,VMIN,VMAX)
      IMPLICIT NONE

      INCLUDE 'vortex_parameters.dat'
      
      INTEGER NX,NY,NZ,BOR,BORAMR,NL,IR
      INTEGER PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV)
      INTEGER NPATCH(0:NLEVELS)
      REAL ARRAY0(1-BOR:NMAX+BOR,1-BOR:NMAY+BOR,1-BOR:NMAZ+BOR)
      REAL ARRAY1(1-BORAMR:NAMRX+BORAMR,1-BORAMR:NAMRY+BORAMR,
     &             1-BORAMR:NAMRZ+BORAMR,NPALEV)

      REAL VMIN,VMAX
      INTEGER N1,N2,N3,I,LOW1,LOW2

      IF (IR.EQ.0) THEN
       VMIN=1.E30
       VMAX=-1.E30
!$OMP PARALLEL DO SHARED(NX,NY,NZ,ARRAY0),
!$OMP+            PRIVATE(I),
!$OMP+            REDUCTION(MIN:VMIN), REDUCTION(MAX:VMAX),
!$OMP+            DEFAULT(NONE)       
       DO I=1,NX
        VMIN=MIN(MINVAL(ARRAY0(I,1:NY,1:NZ)),VMIN)
        VMAX=MAX(MAXVAL(ARRAY0(I,1:NY,1:NZ)),VMAX)
       END DO
      ELSE
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

       VMIN=1.E30
       VMAX=-1.E30
!$OMP PARALLEL DO SHARED(LOW1,LOW2,PATCHNX,PATCHNY,PATCHNZ,NL,
!$OMP+                   ARRAY1),
!$OMP+            PRIVATE(I,N1,N2,N3),
!$OMP+            REDUCTION(MIN:VMIN), REDUCTION(MAX:VMAX),
!$OMP+            DEFAULT(NONE)
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        VMIN=MIN(MINVAL(ARRAY1(1:N1,1:N2,1:N3,I)),VMIN)
        VMAX=MAX(MAXVAL(ARRAY1(1:N1,1:N2,1:N3,I)),VMAX)
       END DO
      END IF

      RETURN
      END SUBROUTINE

********************************************************************
       SUBROUTINE SORT_CELLS(KK_ENTERO,DDD,DDDX,DDDY,DDDZ)
********************************************************************
*      Sorts cells decreasingly in density (DDD)
********************************************************************

       IMPLICIT NONE

       INTEGER KK_ENTERO
       REAL DDD(KK_ENTERO)
       INTEGER DDDX(KK_ENTERO),DDDY(KK_ENTERO),DDDZ(KK_ENTERO)

       REAL DDD2(KK_ENTERO)
       INTEGER DDDX2(KK_ENTERO),DDDY2(KK_ENTERO),DDDZ2(KK_ENTERO)
       INTEGER INDICE2(KK_ENTERO)
       INTEGER I

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD2(I)=1.0/DDD(I)
        DDDX2(I)=DDDX(I)
        DDDY2(I)=DDDY(I)
        DDDZ2(I)=DDDZ(I)
       END DO

       CALL ARGSORT(KK_ENTERO,DDD2,INDICE2)

!$OMP PARALLEL DO SHARED(KK_ENTERO,DDD,DDDX,DDDY,DDDZ,DDD2,DDDX2,
!$OMP+                   DDDY2,DDDZ2,INDICE2),
!$OMP+            PRIVATE(I), DEFAULT(NONE)
       DO I=1,KK_ENTERO
        DDD(I)=1.0/DDD2(INDICE2(I))
        DDDX(I)=DDDX2(INDICE2(I))
        DDDY(I)=DDDY2(INDICE2(I))
        DDDZ(I)=DDDZ2(INDICE2(I))
       END DO

       RETURN
       END

