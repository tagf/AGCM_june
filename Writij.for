      SUBROUTINE WRITIJ (A,SCALE,TITL,IND)
         IMPLICIT REAL*8(A-H,O-Z),INTEGER (I-N)
c
	include 'dir.fi'
c
      DIMENSION A(74,46),S(74,46),ALON(18)
      INTEGER ALAT(46),RLAT
      character*4 titl
         ABS(XXX)=DABS(XXX) 
  81  FORMAT (I3,A1,18F7.1)
  82  FORMAT (I3,A1,18F7.2)
  83  FORMAT (I3,A1,18F7.3)
  84  FORMAT (I3,A1,18F7.4)
  98  FORMAT (1H0,25X,A4,'   scale=',d12.6)
      DATA ANORT/4HNNNN/
      DATA SOUTH/4HSSSS/
      DATA BLANK/4H    /
  91  FORMAT (1H0,3X,18F7.2/)
C
      lun=66
      DO 21 J=1,46
      DO 21 I=1,74
  21  S(I,J)=A(I,J)/SCALE      
C
      SMAX=-1.d20
      SMIN=1.d20
      J1=1
      IF (IND.EQ.2) J1=2
      DO 2 J=J1,46
      DO 2 I=2,73
      SMAX=DMAX1(SMAX,S(I,J))
      SMIN=DMIN1(SMIN,S(I,J))
  2   CONTINUE
      IF (SMIN.LT.0.0) SMIN=10.0*ABS(SMIN)
      GM=DMAX1(SMIN,SMAX)
C
      RLAT=90
      IF (IND.EQ.2) RLAT=88
      RLON=-180.0
      IF (IND.EQ.2) RLON=RLON+5./2.0
C
      I1=2
C
      OPEN(lun,FILE=TRIM(BaseDir)//WorkDir//'\wr#t',ACCESS='APPEND')
      DO 1000 I=1,4
      I2=I1+17
      DO 10 K=1,18
  10  ALON(K)=RLON+(I1+K-3)*5.
      DO 20 K=J1,46
  20  ALAT(K)=RLAT-(K-J1)*4
C
c      OPEN(6,FILE='F:\CLIMATE\wr#T',ACCESS='APPEND',
      WRITE (lun,98)  TITL,scale
C
      WRITE (lun,91) (ALON(K),K=1,18)
C
      DO 100 K=J1,46
      KK=46+J1-K
C
      CHAR=BLANK
      IF (ALAT(K).GT.0) CHAR=ANORT
      IF (ALAT(K).LT.0) CHAR=SOUTH
      ALAT(K)=IABS(ALAT(K))
      IF (GM.LT.1000.0) GO TO 110
      WRITE (lun,81) ALAT(K),CHAR,( S(L,KK),L=I1,I2)
      GO TO 100
 110  IF (GM.LT.100.0) GO TO 120
      WRITE (lun,82) ALAT(K),CHAR,(S(L,KK),L=I1,I2)
      GO TO 100
 120  IF (GM.LT.10.0) GO TO 130
      WRITE (lun,83) ALAT(K),CHAR,(S(L,KK),L=I1,I2)
      GO TO 100
 130  CONTINUE
      WRITE (lun,84) ALAT(K),CHAR,(S(L,KK),L=I1,I2)
 100  CONTINUE
C
      I1=I1+18
 1000 CONTINUE  !I=1,4
      close (lun)
      RETURN
      END
