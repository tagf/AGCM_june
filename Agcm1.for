C     ***************
C     ***************
      SUBROUTINE STEP
C     ***************
C     ***************
C
C
C             CONTROLS SEQUENCE OF TIME STEPS
C
C        MRCH=0  CENTERED IN SPACE AND LEAPFROG IN TIME
C        MRCH=1  CENTERED IN SPACE AND FORWARD IN TIME
C        MRCH=2  CENTERED IN SPACE AND BACKWARD IN TIME
C        MRCH=3  UP-RIGHT UNCENTERED IN SPACE AND BACKWARD IN TIME
C        MRCH=4  DOWN-LEFT UNCENTERED IN SPACE AND BACKWARD IN TIME
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C        QT ARRAY - TEMPORARY STORAGE FOR SIMULATION
C
      COMMON /QTARY / PT(74,46),UT(74,46,2),VT(74,46,2),TT(74,46,2)
     1  ,QWT(74,46,2)
C
      COMMON /QTTARY/ QTT(30636)
C
C************ END OF COMMON ******************
C
      DIMENSION Q(30636),QT(30636),QQ(74,414),QQT(74,414)
      EQUIVALENCE (P(1,1),Q(1),QQ(1,1)),(PT(1,1),QT(1),QQT(1,1))
      DATA N/30636/

C        MATSUNO STEPS
C
ccc      CALL WRITIJ(TT(1,1,2),1.d0,'T100',0)
      DO 10 J=1,414
      QQ(1,J)=QQ(73,J)
   10 QQ(74,J)=QQ(2,J)
C
      DO 20 I=1,N
   20 QT(I)=Q(I)
C
      DTS = DT
C
      DO 100 NS=1,4
      MRCH=1
      CALL COMP1
      CALL COMP2
C
      DO 30 J=1,414
      QQT(1,J)=QQT(73,J)
   30 QQT(74,J)=QQT(2,J)
C
      DO 40 I=1,N
      QTT(I)=Q(I)
      Q(I)=QT(I)
   40 QT(I)=QTT(I)
C
      MRCH=2
      IF (NS.EQ.1) MRCH=3
      IF (NS.EQ.2) MRCH=4
C
      CALL COMP1
      CALL COMP2
C
      DO 50 J=1,414
      QQT(1,J)=QQT(73,J)
   50 QQT(74,J)=QQT(2,J)
C
      DO 60 I=1,N
   60 Q(I)=QT(I)
C
  100 CONTINUE
C
C        LEAPFROG STEPS
C
      MRCH=0
C
      DT=DT+DT
C
      DO 200 NS=5,6
      DO 110 I=1,N
  110 QT(I)=QTT(I)
C
c	write (6,101) ns,mrch
c      WRITE(7,'(A5,i4)') 'MRCH0', MRCH
c      CALL WRITIJ(T(1,1,2),1.d0,'T303',0)
      CALL COMP1
      CALL COMP2
c	write (6,101) ns,mrch
c      CALL WRITIJ(T(1,1,2),1.d0,'T304',0)
C
      DO 120 J=1,414
      QQT(1,J)=QQT(73,J)
  120 QQT(74,J)=QQT(2,J)
C
      DO 130 I=1,N
      QTT(I)=Q(I)
  130 Q(I)=QT(I)
  200 CONTINUE
      DT = DTS
C
c      CALL WRITIJ(T(1,1,2),1.d0,'T3c4',0)
      CALL COMP4
c      CALL WRITIJ(T(1,1,2),1.d0,'T305',0)
C
      RETURN
      END
C     **************
C     **************
      SUBROUTINE GMP
C     **************
C     **************
C
C
C             CHECKS MASS CONSERVATION
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C************ END OF COMMON ******************
C
      DIMENSION ZM(46)
C
      DO 10 J=1,JM
   10 ZM(J)=0.0
C
      DO 20 I=2,73
      DO 20 J=1,JM
   20 ZM(J)=ZM(J)+P(I,J)
C
      DO 30 J=1,JM
   30 ZM(J)=ZM(J)/FIM*DXYP(J)
C
      W=0.0
      G=0.0
      DO 40 J=1,JM
      G=G+ZM(J)
   40 W=W+DXYP(J)
C
      G=G/W+PTROP
      DELP=PSF-G
C
C         MASS CORRECTION
C
      DO 50 J=1,JM
      DO 50 I=1,74
   50 P(I,J)=P(I,J)+DELP
C
ccc      WRITE (6,9) DELP,TAU
ccc    9 FORMAT (17H PRESSURE ADDED =,E16.8,6H TAU =,F12.4)
      RETURN
      END
C     *******************
C     *******************
      SUBROUTINE AVRX(PU)
C     *******************
C     *******************
C
C
C             PERFORMS ZONAL SMOOTHING
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C************ END OF COMMON ******************
C
      DIMENSION PU(74,46)
      DIMENSION  AVCOS(74,36),AVSIN(74,36)
      DIMENSION ANG(74),AVSP(37),AVCP(37)
      DIMENSION SMLATS (36,5), SMLATN (36,5)
      DIMENSION  PUP(74),PUM(74), ALPHA(46), NUMS(46)
         COS(XXX)=DCOS(XXX)
         SIN(XXX)=DSIN(XXX)
         SQRT(XXX)=DSQRT(XXX)
      DATA NUMS/ 5*0, 5, 8*1, 18*0, 8*1, 5, 5*0/
      DATA ALPHA/5*0.0,4.791E-2,1.433E-1,1.076E-1,8.441E-2,6.632E-2,
     1  4.913E-2,3.292E-2,1.777E-2,3.747E-3,18*0.0,3.747E-3,
     2  1.777E-2,3.292E-2,4.913E-2,6.632E-2,8.441E-2,1.076E-1,
     3  1.433E-1,4.791E-2,5*0.0/
C
C         COMPUTE COEFFICIENTS FOR FOURIER FIT
C
      DATA ISET/0/
      IF(ISET.NE.0) GO TO 5
      ISET = 1
      DEFF=1.0/DYP(45)
      DO 200 L=1,36
      PUM(L) = (37-L)*DLON
      PUP(L) = SIN(.5*PUM(L))
      PUP(L)=1.0/PUP(L)
  200 CONTINUE
      DO 300 L=2,36
      CS = 0.
      SS = 0.
      DO 210 I=2,73
  210 ANG(I) = PUM(L)*(I-1)
      DO 220 I=2,73
      AVCOS(I,L) = COS(ANG(I))
  220 AVSIN(I,L) = SIN(ANG(I))
      DO 290 I=2,73
      CS = CS + AVCOS(I,L)*AVCOS(I,L)
  290 SS = SS + AVSIN(I,L)*AVSIN(I,L)
      CS = SQRT(CS)
      SS = SQRT(SS)
      CSA=1.0/SQRT(72.0d0)
      DO 300 I=2,73
      AVCOS(I,L) = AVCOS(I,L)/CS
      AVSIN(I,L) = AVSIN(I,L)/SS
  300 CONTINUE
      DO 305 I=2,73
      AVSIN(I,1)=0.0
  305 AVCOS(I,1)=CSA*(-1.0)**(I-1)
      DO 310 J=2,5
      DO 310 N=1,36
  310 SMLATS(N,J) = 1.0E0 - DXP(J)*PUP(N)*DEFF
C
      DO 320 J=42,45
      DO 320 N=1,36
  320 SMLATN(N,J-41)=1.0E0-DXP(J)*PUP(N)*DEFF
C
    5 CONTINUE
C
      DO 400 J=6,41
      IF (DXP(J) .GE. DYP(45)) GO TO 400
C
C         OLD  AVRX, SUCCESSIVE THREE POINT SMOOTHING
C
      NM = NUMS(J)
      DO 395 N=1,NM
      PU(1,J)=PU(73,J)
      PU(74,J)=PU(2,J)
      DO 370 I=2,74
  370 PUP(I)=PU(I-1,J)-PU(I,J)
      DO 380 I=2,73
  380 PUM(I)=PUP(I)-PUP(I+1)
      DO 390 I=2,73
  390 PU(I,J)=PU(I,J)+ALPHA(J)*PUM(I)
  395 CONTINUE
  400 CONTINUE
c      CALL WRITIJ(pu(1,1),1.d7,'av01',2)
C
C         NEW  AVRX, HALF-FAST FOURIER SPECTRUM SMOOTHING
C         USED IN HIGH LATITUDES
C
C
C         SOUTH POLE
C
      DO 600 J=2,5
      DO 550 I=2,37
      PUP(I) = PU(I,J)+PU(I+36,J)
  550 PUM(I) = PU(I,J)-PU(I+36,J)
      DO 595 M=1,36,2
C
C        EVEN WAVE NUMBER
C
      IF (SMLATS(M,J) .LT. 0.0E0) GO TO 575
      DO 555 I=2,37
      AVSP(I)=AVSIN(I,M)*PUP(I)
  555 AVCP(I)=AVCOS(I,M)*PUP(I)
      BN = 0.E0
      AN = 0.0E0
      DO 560 I=2,37
      AN = AN + AVSP(I)
  560 BN = BN + AVCP(I)
      BN = SMLATS(M,J)*BN
      AN = SMLATS(M,J)*AN
      DO 570 I=2,73
  570 PU(I,J) = PU(I,J) - BN*AVCOS(I,M) - AN*AVSIN(I,M)
C
C        ODD WAVE NUMBER
C
  575 IF (SMLATS(M+1,J).LT.0.E0) GO TO 600
      DO 577 I=2,37
      AVSP(I)=AVSIN(I,M+1)*PUM(I)
  577 AVCP(I)=AVCOS(I,M+1)*PUM(I)
      BN=0.0E0
      AN=0.0E0
      DO 580 I=2,37
      AN = AN + AVSP(I)
  580 BN = BN + AVCP(I)
      AN = SMLATS(M+1,J)*AN
      BN = SMLATS(M+1,J)*BN
      DO 590 I=2,73
  590 PU(I,J) = PU(I,J)-BN*AVCOS(I,M+1)-AN*AVSIN(I,M+1)
  595 CONTINUE
  600 CONTINUE
C
C         NORTH POLE
C
      DO 700 J=42,45
      JJ=J-41
      DO 650 I=2,37
      PUP(I) = PU(I,J)+PU(I+36,J)
  650 PUM(I) = PU(I,J)-PU(I+36,J)
      DO 695 M=1,36,2
C
C         EVEN WAVE NUMBER
C
      IF (SMLATN(M,JJ) .LT. 0.0E0) GO TO 775
      DO 655 I=2,37
      AVSP(I)=AVSIN(I,M)*PUP(I)
  655 AVCP(I)=AVCOS(I,M)*PUP(I)
      BN = 0.E0
      AN = 0.0E0
      DO 660 I=2,37
      AN = AN + AVSP(I)
  660 BN = BN + AVCP(I)
      BN = SMLATN(M,JJ)*BN
      AN = SMLATN(M,JJ)*AN
      DO 670 I=2,73
  670 PU(I,J) = PU(I,J) - BN*AVCOS(I,M) - AN*AVSIN(I,M)
C
C         ODD WAVE NUMBER
C
  775 IF (SMLATN(M+1,JJ).LT.0.E0) GO TO 700
      DO 675 I=2,37
      AVSP(I)=AVSIN(I,M+1)*PUM(I)
  675 AVCP(I)=AVCOS(I,M+1)*PUM(I)
      BN=0.0E0
      AN=0.0E0
      DO 680 I=2,37
      AN = AN + AVSP(I)
  680 BN = BN + AVCP(I)
      AN = SMLATN(M+1,JJ)*AN
      BN = SMLATN(M+1,JJ)*BN
      DO 690 I=2,73
  690 PU(I,J) = PU(I,J)-BN*AVCOS(I,M+1)-AN*AVSIN(I,M+1)
  695 CONTINUE
  700 CONTINUE
C
C        FILL FOR WRAP AROUND
C
      DO 800 J=1,46
      PU(1,J)=PU(73,J)
      PU(74,J)=PU(2,J)
  800 CONTINUE
      RETURN
      END
C     *****************
C     *****************
      SUBROUTINE INTSEA
C     *****************
C     *****************
C
C
C             INTERPOLATES SEA SURFACE TEMPERATURE
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1  ,SD(74,46)
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C        MONTHLY SEA TEMPERATURES
C
      COMMON /SEA/ SS(74,46,12), SWIN(74,46,12)
C
      include 'ice.fi'
C************ END OF COMMON ******************
C
      DIMENSION   FRLEAD(18)
c      LOGICAL  ICE(74)
      DIMENSION SEAT1(74),SEAT2(74),SEAVAL(74)
      DIMENSION  IZMN(12),IHMN(12)
         FLOAT(III)=DFLOAT(III)
C     DATA FRLEAD/ .0313E0,.0241E0,.0236E0,.0316E0,.0362E0,.0588E0,
C    1             .1290E0,.2006E0,.1578E0,.0999E0,.0452E0,.0301E0,
C    1             .0313E0,.0241E0,.0236E0,.0316E0,.0362E0,.0588E0/
      DATA FRLEAD/ 18*0.E0 /
      DATA  IZMN / 31,28,31,30,31,30,31,31,30,31,30,31  /
      DATA  IHMN / 16,13,16,15,16,15,16,16,15,16,15,16  /
      DATA  MSTART / 16 /
C
C        MNTHDY = DAY OF MONTH - 1 TO 31.
C        MONTH = MONTH OF YEAR - 1 TO 12
C        SS(I,J,IPT1) = SEA TEMP. FOR MSTART(15TH) DAY OF MONTH.
C        SS(I,J,IPT2) = SEA TEMP. FOR MSTART DAY OF NEXT MONTH
C            NOTE - CLEAR OCEAN VALUES ARE IN DEG C, ICE IN KELVIN.
C        IZMN(I) = NO. OF DAYS IN MONTH I
C        IHMN(I) = NO. OF DAYS FROM MSTART DAY WHEN ICE CHANGES ARE MADE
C        MSTART = DAY OF MONTH WHEN FD,2 ARRAYS ARE VALID
C        IARY = NO. OF FIRST SEA ARRAY TO INTERPOLATE FROM
C        IPEROD = NO. OF DAY FROM MSTART OF MONTH LONG PERIOD
C        IPT1 = NO. OF MONTH IN FD ARRAY
C        IPT2 = NO. OF MONTH IN SD ARRAY
C
C
C        FIND CORRECT SEA ARRAY COUNTER AND DAY COUNTER
C
      IARY = MONTH-1
      IF( MNTHDY.GE.MSTART ) IARY = MONTH
      IF( IARY.EQ.0 ) IARY = 12
      IF( MNTHDY.GE.MSTART ) IPEROD = MNTHDY-MSTART
      IF( MNTHDY.LT.MSTART ) IPEROD = MNTHDY+IZMN(IARY)-MSTART
C
      IPT1 = IARY
      IPT2 = MOD( IARY,12 )+1
C
C        TEMPERATURES MUST HAVE BEEN READ INTO COMMON /SEA/
C
C        NOW INTERPOLATE NEW SEA TEMPERATURE
C
      FACTOR = FLOAT(IPEROD)/FLOAT(IZMN(IARY))
      DO 90 J=1,JM
      DO 70 I=2,73
      SEAT1(I) = SS(I,J,IPT1)
      SEAT2(I) = SS(I,J,IPT2)
C
      IF(SEAT1(I).LE.200. ) SEAT1(I) = SEAT1(I)+TICE !CLEAR OCEAN
      IF(SEAT2(I).LE.200. ) SEAT2(I) = SEAT2(I)+TICE !CLEAR OCEAN
   70 CONTINUE !I
      DO 80 I=2,73
      SEAVAL(I) = (SEAT2(I)-SEAT1(I))*FACTOR+SEAT1(I) !linear interpol
   80 CONTINUE  !I

C***********************
C      FOR ICE MODELLING
C***********************
	do i=2,73
       SETICE(I,J)=.FALSE.
       if (ISFTYP(I,J).eq.7.or.ISFTYP(I,J).eq.9) then
        IF (IPEROD.LT.IHMN(IARY).AND.SS(I,J,IPT1).GT.200.) !sea ice
     1   SETICE(I,J) = .TRUE.
        IF (IPEROD.GE.IHMN(IARY).AND.SS(I,J,IPT2).GT.200.)  !sea ice
     1   SETICE(I,J) = .TRUE.
	  if (FixedIceBorder) then ! фиксированная граница льда
	   if (SETICE(I,J)) then
	    if (ISFTYP(I,J).ne.9) then
           ISFTYP(I,J) = 9
           GIceComp(I,J)  = GIceCompMax
           GHIce(I,J)     = HIceMin
           GHIce2(I,J)    = 0d0
	     GIcePool(I,J)  = 0d0
	     GIce2Pool(I,J) = 0d0
	     GT(I,J)        = TICE
           GT1(I,J)       = TICE
           GT2(I,J)       = TICE
           SnoAmt1(I,J)   = 0d0
           SnoAmt2(I,J)   = 0d0
	    endif  !.ne.9
	   else !SETICE
	    if (ISFTYP(I,J).ne.7) then
           GT(I,J)     = TICE
           GHICE(I,J)  = 	0.0
           GIceComp(I,J)  = 0.0
           ISFTYP(I,J) = 7
	    endif
	   endif  !SETICE
        endif   !FixedIceBorder
	  if (ISFTYP(I,J).eq.7.and.FixedWater) then ! фиксированная T воды
	   if (.not.SETICE(I,J)) then
          GT(I,J) = SEAVAL(I)
	   else
          GT(I,J) = TICE
	   endif
        endif
c        IF (ICE(I).AND.ISFTYP(I,J).EQ.7) GT(I,J)=TICE
c        IF (.NOT.(ICE(I).OR.TDSEA.OR.SWAMP)) GT(I,J)=SEAVAL(I)
c        ISFTYP(I,J)=7
c        IF (ICE(I)) ISFTYP(I,J)=9
c        IF ((.not.ICE(I)).and.FixedWater) GT(I,J)=SEAVAL(I)
        IF (ISFTYP(I,J).EQ.7) SNOAMT(I,J)=0.0
       endif
	enddo  !i
C
c      DO 90 I=2,73	 !if no ice model
c      ICE(I)=.FALSE.
c      IF (ISFTYP(I,J).NE.7.AND.ISFTYP(I,J).NE.9) GO TO 90
c      IF (IPEROD.LT.IHMN(IARY).AND.SS(I,J,IPT1).GT.200.) ICE(I)=.TRUE.
c      IF (IPEROD.GE.IHMN(IARY).AND.SS(I,J,IPT2).GT.200.) ICE(I)=.TRUE.
c      IF (ICE(I).AND.ISFTYP(I,J).EQ.7) GT(I,J)=TICE
c      IF (.NOT.ICE(I))  GT(I,J)=SEAVAL(I)
c      ISFTYP(I,J)=7
c      IF (ICE(I)) ISFTYP(I,J)=9
cC
c      IF (ISFTYP(I,J).EQ.7) SNOAMT(I,J)=0.0
   90 CONTINUE  !J
C***********************
C
      FLEADN = FRLEAD(MONTH)
      FLEADS = FRLEAD(MONTH+6)
      RETURN
      END
C     ****************
C     ****************
      SUBROUTINE INTOZ
C     ****************
C     ****************
C
C
C             INTERPOLATES OZONE AMOUNTS
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C************ END OF COMMON ******************
C
      DIMENSION   MONTHS(12)
      DIMENSION O3(46,12)
      DIMENSION O3JAN(46,3),O3APR(46,3),O3JUL(46,3),O3OCT(46,3)
      EQUIVALENCE (O3(1,1),O3JAN(1,1)),(O3(1,4),O3APR(1,1)),
     *  (O3(1,7),O3JUL(1,1)),(O3(1,10),O3OCT(1,1))
         FLOAT(III)=DFLOAT(III)
      DATA  MONTHS/31,28,31,30,31,30,31,31,30,31,30,31  /
      DATA O3JAN/    306., 314., 321., 330., 337., 340.,
     1 338., 331., 324., 317., 311., 306., 300., 294., 289.,
     2 284., 278., 272., 266., 260., 256., 252., 250., 250.,
     3 250., 251., 253., 256., 264., 274., 285., 298., 310.,
     4 320., 334., 353., 362., 364., 366., 365., 363., 362.,
     5 360., 360., 360., 360.,
     1            296., 297., 299., 303., 310., 317., 312., 300.,
     1 297., 296., 294., 292., 289., 286., 282., 277., 272., 267.,
     2 262., 258., 256., 254., 253., 254., 255., 257., 260., 265.,
     3 276., 291., 304., 316., 330., 350., 365., 374., 378., 381.,
     4 385., 388., 390., 391., 393., 394., 395., 395.,
     1            290., 291., 293., 295., 298., 301., 299., 299.,
     1 300., 302., 300., 295., 289., 283., 278., 272., 267., 262.,
     2 258., 254., 253., 253., 253., 255., 258., 263., 270., 280.,
     3 290., 302., 316., 330., 347., 360., 374., 385., 394., 400.,
     4 408., 413., 420., 426., 434., 438., 442., 442./
      DATA O3APR/ 286., 287., 289., 292., 298., 300., 299., 301.,
     1 306., 311., 313., 305., 296., 287., 279., 271., 266., 261.,
     2 258., 257., 257., 257., 258., 260., 264., 270., 275., 283.,
     3 293., 305., 319., 331., 347., 362., 375., 387., 400., 411.,
     4 420., 424., 429., 432., 435., 438., 439., 441.,
     1            279., 280., 281., 280., 280., 285., 300., 311.,
     1 320., 328., 332., 325., 310., 298., 285., 276., 268., 262.,
     2 259., 258., 258., 260., 262., 264., 267., 271., 275., 281.,
     3 290., 301., 315., 326., 339., 350., 361., 370., 380., 386.,
     4 390., 394., 396., 399., 401., 403., 404., 405.,
     1            269., 274., 278., 280., 280., 279., 292., 311.,
     1 327., 342., 354., 350., 335., 320., 300., 287., 275., 267.,
     2 261., 258., 259., 260., 262., 264., 266., 269., 272., 277.,
     3 284., 293., 302., 312., 325., 333., 343., 355., 361., 365.,
     4 367., 369., 370., 370., 370., 369., 367., 365./
      DATA O3JUL/ 276., 280., 282., 281., 280., 282., 302., 317.,
     1 330., 360., 369., 368., 358., 340., 320., 300., 286., 275.,
     2 265., 260., 259., 261., 261., 260., 259., 261., 265., 270.,
     3 279., 285., 294., 300., 310., 320., 328., 334., 340., 344.,
     4 348., 350., 350., 348., 346., 344., 343., 342.,
     1            290., 291., 293., 295., 299., 305., 314., 327.,
     1 342., 372., 385., 384., 372., 354., 331., 310., 298., 285.,
     2 277., 272., 268., 264., 260., 258., 257., 258., 260., 263.,
     3 269., 276., 283., 290., 297., 303., 312., 320., 324., 327.,
     4 329., 330., 330., 330., 329., 327., 325., 322.,
     1            297., 298., 300., 305., 310., 318., 329., 348.,
     1 370., 391., 409., 402., 384., 364., 340., 320., 305., 292.,
     2 282., 274., 267., 262., 258., 256., 256., 257., 258., 260.,
     3 265., 270., 275., 280., 286., 291., 296., 300., 305., 309.,
     4 311., 311., 310., 308., 304., 301., 298., 297./
      DATA O3OCT/ 320., 323., 329., 335., 342., 349., 360., 380.,
     1 400., 403., 400., 380., 363., 349., 334., 317., 305., 293.,
     2 282., 273., 266., 259., 255., 252., 251., 252., 259., 257.,
     3 262., 267., 272., 276., 280., 285., 289., 294., 298., 302.,
     4 306., 310., 308., 304., 299., 295., 291., 290.,
     1            378., 382., 382., 360., 357., 361., 382., 391.,
     1 392., 382., 367., 355., 340., 329., 317., 305., 295., 286.,
     2 276., 268., 260., 255., 251., 250., 249., 251., 253., 255.,
     3 260., 266., 272., 276., 281., 286., 293., 298., 303., 308.,
     4 313., 317., 318., 318., 317., 314., 308., 300.,
     1            359., 358., 258., 358., 359., 362., 363., 361.,
     1 357., 343., 333., 325., 314., 304., 296., 288., 280., 274.,
     2 268., 262., 258., 254., 251., 250., 250., 251., 253., 256.,
     3 262., 270., 277., 285., 292., 300., 308., 317., 327., 337.,
     4 341., 343., 344., 344., 344., 342., 340., 337./
C
C        NOW CALCULATE ZONAL OZONE VALUES - INTERPEROLATE FROM MONTHLY
C          VALUES.  VALUES AT 15TH OF MONTH.
C
      I1 = MONTH
      IF( MNTHDY.LT.15 ) I1 = I1-1
      IF( I1.EQ.0 ) I1 = 12
      I2 = I1+1
      IF( I2.EQ.13 ) I2 = 1
      TEMP = 0.0
      IF( MNTHDY.LT.15 ) TEMP = MONTHS(I1)
      VALUE = (MNTHDY+TEMP-15.)/FLOAT(MONTHS(I1))
      DO 290 J=1,JM
  290 O3AMT(J) = (O3(J,I1)+VALUE*(O3(J,I2)-O3(J,I1)))/1000.
      RETURN
      END
C     ***************
C     ***************
      SUBROUTINE SDET
C     ***************
C     ***************
C
C
C             UPDATES DAY, EARTH-SUN DISTANCE, AND SOLAR DECLINATION
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C************ END OF COMMON ******************
C
      DIMENSION   MONTHS(12)
         COS(XXX)=DCOS(XXX)
         SIN(XXX)=DSIN(XXX)
      DATA  MONTHS/31,28,31,30,31,30,31,31,30,31,30,31  /
      MAXDAY=DAYPYR + 1.0E-2
      SDEDY=SDEDY+1
      IF (SDEDY .LE. MAXDAY) GO TO 211
      SDEDY=SDEDY-MAXDAY
      SDEYR=SDEYR+1.0
  211 JDYACC=0
      DO 251 L=1,12
      MONTH = L
      JDYACC=JDYACC+MONTHS(L)
      IF (SDEDY .LE. JDYACC) GO TO 241
  251 CONTINUE
  241 MNTHDY=MONTHS(L)-JDYACC+SDEDY
      DY=SDEDY
      SEASON=(DY-SOLTCE)/DAYPYR
      DIST=(DY-APHEL )/DAYPYR
C
C        SOLTCE = JUNE 22
C        APIHELION = JULY 1
C        ECCN= ORBITAL ECCENTRICITY
C
      DEC=DECMAX*COS(2.0*PI*SEASON)
      RSDIST=(1.0+ECCN*COS(2.0*PI*DIST))**2
      S0=2793.6/RSDIST
      SIND=SIN(DEC)
      COSD=COS(DEC)
      RETURN
      END
