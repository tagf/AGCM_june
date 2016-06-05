C     ****************
C     ****************
      SUBROUTINE COMP1
C     ****************
C     ****************
C
C
C             COMPUTES ADVECTION OF MOMENTUM,TEMPERATURE AND MOISTURE
C
C
C************ BEGINNING OF COMMON ************
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      REAL*8 KAPA,KAPEL,LAT
      INTEGER  SDEDY
C
      COMMON /RECOM/
     1 TAU   ,TAUI  ,TAUO  ,TAUD  ,TAUE  ,TAUH  ,TAUC  ,TOFDAY,ROT   ,
     2 DT    ,DLAT  ,DLON  ,RAD   ,RSDIST,SIND  ,COSD  ,COSR  ,SINR  ,
     3 DAYPYR,ROTPER,SDEYR ,SOLTCE,APHEL ,DECMAX,ECCN  ,GUSTY ,
     4 DAY   ,GRAV  ,RGAS  ,KAPA  ,PSF   ,PTROP ,PSL   ,CFG   ,
     5 FM    ,ED    ,PI    ,SIG1  ,SIG3  ,DSIG  ,PM    ,KAPEL ,RKAPA1,
     6 STBO  ,GWM   ,DTC3  ,FRDAY ,CLH   ,COE1  ,HICE  ,CTYICE,CNDICE,
     7 TICE  ,TCICE ,SNOWL ,COE   ,TSPD  ,PSTQ  ,QST   ,TST   ,PTRK  ,
     8 PSFHO ,CALFA ,QC    ,PC    ,QCONST,EFVCON,TCT0  ,TCT2  ,TCT4  ,
     9 TCST0 ,TCST2 ,TC01  ,TC02  ,TC03  ,TC04  ,TC23  ,TC24  ,TC34  ,
     1 BLC   ,ALOGP0,FIM   ,HRGAS ,TCST4 ,FLR   ,PS4K  ,PS8K  ,ELOG  ,
     2 S0,TOZONE    ,LAT(46)   ,DXU(46)   ,DXP(46)   ,DYU(46),DYP(46),
     3 SFCALB(9,2),SINL(46)  ,COSL(46)  ,O3AMT(46) ,TC12  ,TC14  ,
     4 DXYP(46)  ,F(46)     ,SIG(2)    ,FLEADN    ,FLEADS    ,ERROR  ,
     5 COSLN(72),SINLN(72),DXYU1(46),DXYU2(46), GMT,GMR,GMKE,
     6 JM,IM,ID,MNTHDY,SDEDY ,NCYCLE,NC3   ,MONTH ,MRCH  ,NSTEP
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)
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
      COMMON /COMP/ SM(74,46,2),SAVPU(74,46,2)
C
      COMMON /WORK/
     1   PV(74,46),FD(74,46),H(74,46),PU(74,46),
     1   FLUXT1(74),FLUXU1(74),FLUXV1(74),FLUXQ(74),
     1   FLUXQ1(74),DTDXY(46),FDXU(46),FDYU(46),
     1    FDU(74),FLUX(74),FDQ(74)
C
C************ END OF COMMON ******************
C
C****
       DIMENSION VTUR(74,46,2),PVTUR(74,46),FDTUR(74,46),GTUR(46,2),
     1           TT23(2),TT24(2)
c          DATA GTUR/0.0d0,3.0d0,3.0d0,3.0d0,3.0d0,3.0d0,2.0d0,0.0d0,
c 	1             0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,
c	2			 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,0.0d0,
c     3             0.0d0,0.0d0,0.0d0,2.0d0,2.0d0,2.0d0,2.0d0,0.0d0,
c     4             0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,
c	5			 3.0d0,3.0d0,3.5d0,3.5d0,3.5d0,1.0d0,1.0d0/
C      ZONAL MEANS TT(I,23,L),TT(I,24,L) ARE TT23(L),TT24(L)

       do 1151 j=1,17
 1151  gtur(j,1)=300. !3.
       do 1154 j=28,46
 1154  gtur(j,1)=200. !1.

       do 1152 j=1,12
 1152  gtur(j,2)=100. !2.
       do 1153 j=29,46
 1153  gtur(j,2)=100. !3.
    
       do 1150 j=1,46
       gtur(j,1)=0. !100.
 1150  gtur(j,2)=0. !100.

       TT23(1)=0.
       TT23(2)=0.
       TT24(1)=0.
       TT24(2)=0.
       DO 115 L=1,2
       DO 115 I=2,73
       TT23(L)=TT23(L)+TT(I,23,L)
 115   TT24(L)=TT24(L)+TT(I,24,L)
       DO 116 L=1,2
       TT23(L)=TT23(L)/72.
 116   TT24(L)=TT24(L)/72.
C              L=1
C             SOUTH
         DO 10 J=2,23
         DO 10 I=1,74
         TTTUR2=DABS(TT23(1)-TT(I,1,1))/26.
         TTTUR=TT(I,J,1)-TT(I,J-1,1)
         TTTUR1=0.20*GTUR(J,1)*DABS(TTTUR)
 10      VTUR(I,J,1)=-0.003*TTTUR*TTTUR1*TTTUR2
C              L=1
C             NORTH
         DO 11 J=24,46
         DO 11 I=1,74
         TTTUR2=DABS(TT24(1)-TT(I,46,1))/17.
         TTTUR=TT(I,J,1)-TT(I,J-1,1)
         TTTUR1=0.15*GTUR(J,1)*DABS(TTTUR)
 11      VTUR(I,J,1)=-0.003*TTTUR*TTTUR1*TTTUR2
C              L=2
C             SOUTH
         DO 12 J=2,23
         DO 12 I=1,74
         TTTUR2=DABS(TT23(2)-TT(I,1,2))/35.
         TTTUR=TT(I,J,2)-TT(I,J-1,2)
         TTTUR1=0.20*GTUR(J,2)*DABS(TTTUR)
 12      VTUR(I,J,2)=-0.003*TTTUR*TTTUR1*TTTUR2
C              L=2
C             NORTH
         DO 13 J=24,46
         DO 13 I=1,74
         TTTUR2=DABS(TT24(2)-TT(I,46,2))/20.
         TTTUR=TT(I,J,2)-TT(I,J-1,2)
         TTTUR1=0.10*GTUR(J,2)*DABS(TTTUR)
 13      VTUR(I,J,2)=-0.003*TTTUR*TTTUR1*TTTUR2*TTTUR2
C
C        CALL PRTABL(VTUR(1,1,1),26,16,4HTUR1,4)
C        CALL PRTABL(VTUR(1,1,2),26,16,4HTUR2,4)
C***
      JMM1=JM-1
      JMM2=JM-2
C************************************************************
C             MASS WEIGHT VARIABLES
C************************************************************
      DO 15 J=1,JM
      DO 15 I=1,74
   15 FD(I,J)=PT(I,J)*DXYP(J)
c      call PrnMap(1,fd(1,1),74,46,0.d0,1.d14,' FD1')
c      call PrnMap(1,ut(1,1,1),74,46,0.d0,1.d0,' ut1')
c      call PrnMap(1,ut(1,1,2),74,46,0.d0,1.d0,' ut3')
c      call PrnMap(1,vt(1,1,1),74,46,0.d0,1.d0,' vt1')
c      call PrnMap(1,vt(1,1,2),74,46,0.d0,1.d0,' vt3')
C
      DO 20 L=1,2
      DO 20 J=1,JM
      DO 20 I=2,73
      TT(I,J,L) = TT(I,J,L)*FD(I,J)
   20 QWT(I,J,L) = QWT(I,J,L)*FD(I,J)
C
C        DOUBLE POLAR VALUES
C
      DO 30 I=1,74
      FD(I,1)=FD(I,1)+FD(I,1)
   30 FD(I,JM)=FD(I,JM)+FD(I,JM)
C
      DO 35 J=1,JM
      DO 35 I=2,73
   35 PU(I,J)=FD(I,J)+FD(I+1,J)
C
      DO 37 J=2,JM
      DO 37 I=2,73
   37 FD(I,J)=0.25*(PU(I,J)+PU(I,J-1))
C
      DO 40 L=1,2
      DO 40 J=2,JM
      DO 40 I=2,73
      UT(I,J,L)=UT(I,J,L)*FD(I,J)
      VT(I,J,L)=VT(I,J,L)*FD(I,J)
   40 CONTINUE
C
      FC=0.5
      IF (MRCH.LE.2) FC=0.25
C
      DO 45 J=1,JM
      FDXU(J)=FC*DXU(J)
   45 FDYU(J)=FC*DYU(J)
C************************************************************
C             COMPUTE MASS FLUX
C************************************************************
      DO 430 L=1,2
C
      DO 50 J=2,JM
      DO 50 I=2,73
   50 PU(I,J)=FDYU(J)*U(I,J,L)
C
      IF (MRCH.EQ.3) GO TO 80
      IF (MRCH.EQ.4) GO TO 120
C
C        MRCH=0 OR 1 OR 2
C
      DO 60 J=2,JMM1
      DO 60 I=2,73
   60 PU(I,J)=PU(I,J)+PU(I,J+1)
C
      DO 70 J=2,JM
      DO 70 I=2,73
   70 PV(I,J) = FDXU(J)*(V(I,J,L)+V(I-1,J,L))*
     1               (P(I,J)+P(I,J-1))
C
C****
      DO 71 J=2,JM
      DO 71 I=2,73
   71 PVTUR(I,J) = FDXU(J)*(V(I,J,L)+V(I-1,J,L)+2.*VTUR(I,J,L))*
     1               (P(I,J)+P(I,J-1))
C****
      GO TO 150
C
C        MRCH = 3
C
   80 CONTINUE
C
      DO 90 J=2,JMM1
      DO 90 I=2,73
   90 PU(I,J) = PU(I,J+1)
C
      DO 110 J=2,JM
      DO 110 I=2,73
  110 PV(I,J) = FDXU(J)*V(I,J,L)*(P(I,J)+P(I,J-1))
C
C****
      DO 111 J=2,JM
      DO 111 I=2,73
  111 PVTUR(I,J) = FDXU(J)*(V(I,J,L)+VTUR(I,J,L))*(P(I,J)+P(I,J-1))
C****
      GO TO 150
C
C        MRCH = 4
C
  120 CONTINUE
C
      DO 140 J=2,JM
      DO 140 I=2,73
  140 PV(I,J) = FDXU(J)*V(I-1,J,L)*(P(I,J)+P(I,J-1))
C
C****
      DO 141 J=2,JM
      DO 141 I=2,73
  141 PVTUR(I,J) = FDXU(J)*(V(I-1,J,L)+VTUR(I,J,L))*(P(I,J)+P(I,J-1))
C****
  150 CONTINUE
C
      CALL AVRX (PU)
C
      DO 155 J=2,JMM1
      DO 155 I=2,73
  155 SAVPU(I,J,L)=PU(I,J)*DT*KAPA/(RGAS+RGAS)
C
      DO 160 J=2,JMM1
      DO 160 I=2,73
  160 PU(I,J)=PU(I,J)*(P(I,J)+P(I+1,J))
C************************************************************
C             COMPUTE MASS CONVERGENCE
C************************************************************
      DO 165 J=1,JM
      PU(1,J)=PU(73,J)
      PU(74,J)=PU(2,J)
      PV(1,J)=PV(73,J)
      PV(74,J)=PV(2,J)
C****
      PVTUR(1,J)=PVTUR(73,J)
      PVTUR(74,J)=PVTUR(2,J)
C****
  165 CONTINUE
C
      DO 170 J=2,JMM1
      DO 170 I=2,73
  170 SM(I,J,L)=0.5*(-PU(I,J)+PU(I-1,J)-PV(I,J+1)+PV(I,J))
      DO 175 J=2,JM
      DO 175 I=2,73
  175 PIV(I,J,L)=PV(I,J)
C************************************************************
C             EQUIVALENT PU AT POLES
C************************************************************
      VM1=SUMma(PV(1,2))/FIM
      VM2=SUMma(PV(1,JM))/FIM
C
      DO 185 I=2,73
      SM(I,1,L)=-0.5*VM1
      SM(I,JM,L)=0.5*VM2
  185 CONTINUE
C
      FLUXV1(2)=0.0
      DO 188 I=3,73
  188 FLUXV1(I)=PV(I,2)-VM1
      DO 190 I=3,73
  190 FLUXV1(I)=FLUXV1(I-1)+FLUXV1(I)
C
      VM1=SUM(FLUXV1)/FIM
C
      DO 220 I=2,73
  220 PU(I,1) =0.  !-(FLUXV1(I)-VM1)*3.0 *********************
C
      FLUXV1(2)=0.0
C
      DO 225 I=3,73
  225 FLUXV1(I)=PV(I,JM)-VM2
      DO 230 I=3,73
  230 FLUXV1(I) = FLUXV1(I-1)+FLUXV1(I)
C
      VM2=SUM(FLUXV1)/FIM
C
      DO 242 I=2,73
  242 PU(I,JM) =0.  !(FLUXV1(I)-VM2)*3.0 ***********************	
C
      PU(1,1)=PU(73,1)
      PU(74,1)=PU(2,1)
      PU(1,JM)=PU(73,JM)
      PU(74,JM)=PU(2,JM)
C
      DT2 = 0.5*DT
C********************************************************************
C     HORIZONTAL ADVECTION OF THERMODYNAMIC ENERGY AND MOISTURE
C********************************************************************
      DO 245 J=2,JMM1
      DO 245 I=2,74
  245 FD(I,J)=DT2*PU(I-1,J)
c	 write(6,101) (pu(i,2),i=1,5)
C
      DO 250 J=3,JMM2
      DO 250 I=2,74
  250 H(I,J)=(T(I-1,J,L)+T(I,J,L))*FD(I,J)
C
	iii=iii+1		
c	write (6,101) iii,mrch,l
c 101	 format (' Comp1 i= ',i2,' mrch= ',i2,'l= ',i2)
c      CALL WRITIJ(pu(1,1),1.d10,'PUc1',0)
c	 write(6,101) (pu(i,2),i=1,5)
c      CALL WRITIJ(h(1,1),1.d12,'h c1',0)
      DO 252 I=2,74
      H(I,2)=CVMGP(T(I-1,2,L),T(I,2,L),FD(I,2))
      H(I,JMM1)=CVMGP(T(I-1,JMM1,L),T(I,JMM1,L),FD(I,JMM1))
      H(I,2)=(H(I,2)+H(I,2))*FD(I,2)
      H(I,JMM1)=(H(I,JMM1)+H(I,JMM1))*FD(I,JMM1)
  252 CONTINUE
C
      DO 254 J=2,JMM1
      DO 254 I=2,73
  254 TT(I,J,L)=TT(I,J,L)+H(I,J)-H(I+1,J)
C
      IM1=73
      DO 285 I=2,73
      DO 260 J=2,JMM1
      FLUXQ(J) = (QW(I-1,J,L)+QW(I,J,L))
      FLUXQ1(J) = 4.0*QW(I-1,J,L)*QW(I,J,L)/
     1    (FLUXQ(J)+1.E-30)
      FDQ(J)=FD(I,J)*(QW(I-1,J,L)-QW(I,J,L))
  260 CONTINUE
C
      DO 262 J=2,JMM1
      FLUXQ(J)=FD(I,J)*CVMGP(FLUXQ(J),FLUXQ1(J),FDQ(J))
  262 CONTINUE
C
      DO 264 J=2,JMM1
      FLUXQ(J)=dMAX1(-QWT(I,J,L),dMIN1(QWT(IM1,J,L),FLUXQ(J)))
      QWT(IM1,J,L)=QWT(IM1,J,L)-FLUXQ(J)
      QWT(I,J,L)=QWT(I,J,L)+FLUXQ(J)
  264 CONTINUE
C
      DO 268 J=2,JM
  268 FD(I,J)=DT2*PV(I,J)
C
C****
      DO 269 J=2,JM
  269 FDTUR(I,J)=DT2*PVTUR(I,J)
C****
      DO 280 J=2,JM
      FLUXQ(J)=QW(I,J-1,L)+QW(I,J,L)
      FLUXQ1(J)=4.0*QW(I,J-1,L)*QW(I,J,L)/
     1  (FLUXQ(J)+1.E-30)
  280 FDQ(J)=FD(I,J)*(QW(I,J-1,L)-QW(I,J,L))
C
      DO 282 J=2,JM
  282 FLUXQ(J)=FD(I,J)*CVMGP(FLUXQ(J),FLUXQ1(J),FDQ(J))
C
      DO 284 J=2,JM
      FLUXQ(J)=dMAX1(-QWT(I,J,L),dMIN1(QWT(I,J-1,L),FLUXQ(J)))
      QWT(I,J,L)=QWT(I,J,L)+FLUXQ(J)
  284 QWT(I,J-1,L)=QWT(I,J-1,L)-FLUXQ(J)
C
      IM1=I
C
  285 CONTINUE
C
      DO 290 J=3,JMM1
      DO 290 I=2,73
  290 H(I,J)=T(I,J,L)+T(I,J-1,L)
C
      DO 292 I=2,73
      H(I,2)=CVMGP(T(I,1,L),T(I,2,L),FD(I,2))
      H(I,JM)=CVMGP(T(I,JMM1,L),T(I,JM,L),FD(I,JM))
      H(I,2)=H(I,2)+H(I,2)
  292 H(I,JM)=H(I,JM)+H(I,JM)
C
      DO 294 J=2,JM
      DO 294 I=2,73
C****
      H(I,J)=H(I,J)*FDTUR(I,J)
C****
      TT(I,J,L)=TT(I,J,L)+H(I,J)
  294 TT(I,J-1,L)=TT(I,J-1,L)-H(I,J)
C************************************************************
C             HORIZONTAL ADVECTION OF MOMENTUM
C************************************************************
C
C        PU(I,J) TO PU(I,J)+PU(I-1,J)
C
      DO 310 J=1,JM
      DO 300 I=2,73
  300 FLUXV1(I)=PU(I,J)+PU(I-1,J)
      DO 310 I=2,73
  310 PU(I,J)=FLUXV1(I)
C
C        PV(I,J) TO PV(I,J)+PV(I,J+1)
C
      DO 320 J=2,JMM1
      DO 320 I=2,73
  320 PV(I,J)=PV(I,J)+PV(I,J+1)
C
      DO 330 J=1,JM
      PV(1,J)=PV(73,J)
      PV(74,J)=PV(2,J)
      PU(1,J)=PU(73,J)
      PU(74,J)=PU(2,J)
  330 CONTINUE
C
      DT12=DT/12.
      DT24=DT/24.
      DO 360 J=2,JM
C
      DO 350 I=2,73
      FLUX(I)=DT12*(PU(I,J)+PU(I,J-1))
      FLUXU1(I)=FLUX(I)*(U(I,J,L)+U(I-1,J,L))
  350 FLUXV1(I)=FLUX(I)*(V(I,J,L)+V(I-1,J,L))
C
      FLUXU1(74)=FLUXU1(2)
      FLUXV1(74)=FLUXV1(2)
C
      DO 355 I=2,73
      UT(I,J,L)=UT(I,J,L)+FLUXU1(I)-FLUXU1(I+1)
  355 VT(I,J,L)=VT(I,J,L)+FLUXV1(I)-FLUXV1(I+1)
  360 CONTINUE
C
      DO 368 J=2,JMM1
      DO 365 I=2,73
      FLUX(I)=DT12*(PV(I,J)+PV(I+1,J))
      FLUXU1(I)=FLUX(I)*(U(I,J,L)+U(I,J+1,L))
      UT(I,J+1,L)=UT(I,J+1,L)+FLUXU1(I)
      UT(I,J,L)=UT(I,J,L)-FLUXU1(I)
      FLUXV1(I)=FLUX(I)*(V(I,J,L)+V(I,J+1,L))
      VT(I,J+1,L)=VT(I,J+1,L)+FLUXV1(I)
  365 VT(I,J,L)=VT(I,J,L)-FLUXV1(I)
  368 CONTINUE
C
      DO 400 J=2,JMM1
      DO 370 I=2,73
      FLUX(I)=DT24*(PU(I,J)+PV(I,J))
      FLUXU1(I)=FLUX(I)*(U(I,J+1,L)+U(I-1,J,L))
  370 FLUXV1(I)=FLUX(I)*(V(I,J+1,L)+V(I-1,J,L))
C
      FLUXU1(74)=FLUXU1(2)
      FLUXV1(74)=FLUXV1(2)
C
      DO 375 I=2,73
      UT(I,J+1,L)=UT(I,J+1,L)+FLUXU1(I)
      UT(I,J,L)=UT(I,J,L)-FLUXU1(I+1)
      VT(I,J+1,L)=VT(I,J+1,L)+FLUXV1(I)
  375 VT(I,J,L)=VT(I,J,L)-FLUXV1(I+1)
C
      DO 380 I=2,73
      FLUX(I)=DT24*(-PU(I,J)+PV(I,J))
      FLUXU1(I)=FLUX(I)*(U(I-1,J+1,L)+U(I,J,L))
  380 FLUXV1(I)=FLUX(I)*(V(I-1,J+1,L)+V(I,J,L))
C
      FLUXU1(74)=FLUXU1(2)
      FLUXV1(74)=FLUXV1(2)
C
      DO 385 I=2,73
      UT(I,J+1,L)=UT(I,J+1,L)+FLUXU1(I+1)
      UT(I,J,L)=UT(I,J,L)-FLUXU1(I)
      VT(I,J+1,L)=VT(I,J+1,L)+FLUXV1(I+1)
  385 VT(I,J,L)=VT(I,J,L)-FLUXV1(I)
  400 CONTINUE
  430 CONTINUE
C************************************************************
C        VERTICAL ADVECTION OF MOMENTUM
C************************************************************
      DO 455 J=1,JM
  455 DTDXY(J)=DT/DXYP(J)
      DO 460 J=1,JM
      DO 460 I=2,73
      SD(I,J)=SM(I,J,1)-SM(I,J,2)
  460 PT(I,J)=PT(I,J)+(SM(I,J,1)+SM(I,J,2))*DTDXY(J)
C
      DO 462 J=1,JM
      SD(74,J)=SD(2,J)
  462 SD(1,J)=SD(73,J)
C
      DT8=0.125*DT
C
      DO 470 J=1,JM
      DO 470 I=2,73
  470 FD(I,J)=SD(I,J)+SD(I+1,J)
C
C        DOUBLE POLAR VALUES
C
      DO 475 I=2,73
      FD(I,JM)=FD(I,JM)+FD(I,JM)
  475 FD(I,1)=FD(I,1)+FD(I,1)
C
      DO 480 J=2,JM
      DO 480 I=2,73
  480 PV(I,J)=DT8*(FD(I,J)+FD(I,J-1))
C
      DO 490 J=2,JM
      DO 490 I=2,73
      FDU(I)=PV(I,J)*(U(I,J,1)+U(I,J,2))
      UT(I,J,2)=UT(I,J,2)+FDU(I)
  490 UT(I,J,1)=UT(I,J,1)-FDU(I)
C
      DO 500 J=2,JM
      DO 500 I=2,73
      FDU(I)=PV(I,J)*(V(I,J,1)+V(I,J,2))
      VT(I,J,2)=VT(I,J,2)+FDU(I)
  500 VT(I,J,1)=VT(I,J,1)-FDU(I)
c      DO 501 J=1,JM
c      DO 501 I=2,73
c      dpt(i,j)=pt(i,j)-p(i,j)
c      dut(i,j)=ut(i,j,1)-dut(i,j)
c 501  dvt(i,j)=vt(i,j,1)-dvt(i,j)
c      write (*,*) '    comp1'
c      DO 502 J=1,JM
c      dpt(74,j)=dpt(2,j)
c      dut(74,j)=dut(2,j)
c      dvt(74,j)=dvt(2,j)
c      dpt(1,j)=dpt(73,j)
c      dut(1,j)=dut(73,j)
c 502  dvt(1,j)=dvt(73,j)
c      CALL WRITIJ(dpt,1.d0,'dpt ',0)
c      CALL WRITIJ(dut,1.d15,'dut1  ',2)
c      CALL WRITIJ(dvt,1.d15,'dvt1  ',2)
 
c      call PrnMap(1,sm(1,1,1),74,46,0.d0,1.d8,' SM1')
c      call PrnMap(1,sm(1,1,2),74,46,0.d0,1.d8,' SM3')
      RETURN
      END
