C     ***********************************************
C     ***********************************************
      SUBROUTINE AMXMN(A,ZMX,GMX,IX,JX,ZMN,GMN,IN,JN)
C     ***********************************************
C     ***********************************************
C
C
C             GIVEN ARRAY A
C             FIND THE ZONAL MAX AND MIN
C             AND THE GLOBAL MAX AND MIN AND (I,J) POINTS
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION A(74,46),ZMX(46),ZMN(46)
C
      DO 10 J=1,46
      ZMX(J)=A(2,J)
   10 ZMN(J)=A(2,J)
      DO 20 I=3,73
      DO 20 J=1,46
      ZMX(J)=dMAX1(ZMX(J),A(I,J))
   20 ZMN(J)=dMIN1(ZMN(J),A(I,J))
      JX=1
      GMX=ZMX(1)
      JN=1
      GMN=ZMN(1)
      DO 30 J=2,46
      IF (GMX.GT.ZMX(J)) GO TO 25
      GMX=ZMX(J)
      JX=J
   25 IF (GMN.LT.ZMN(J)) GO TO 30
      GMN=ZMN(J)
      JN=J
   30 CONTINUE
      IX=0
      IN=0
      DO 40 I=2,73
      IF (GMX.EQ.A(I,JX)) IX=I-1
      IF (GMN.EQ.A(I,JN)) IN=I-1
   40 CONTINUE
      RETURN
      END
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
	include 'recom.fi'
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
c      CALL WRITIJ(pu(1,1),1.d7,'0000',2)
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
  150 CONTINUE
C
      CALL AVRX (PU)
c      CALL WRITIJ(pu(1,1),1.d7,'0001',2)
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
      VM1=SUMma(FLUXV1)/FIM
C
      DO 220 I=2,73
  220 PU(I,1) = -(FLUXV1(I)-VM1)*3.0
C
      FLUXV1(2)=0.0
C
      DO 225 I=3,73
  225 FLUXV1(I)=PV(I,JM)-VM2
      DO 230 I=3,73
  230 FLUXV1(I) = FLUXV1(I-1)+FLUXV1(I)
C
      VM2=SUMma(FLUXV1)/FIM
	mx=72
!      WRITE(7,'(A5,5E25.16)') 'PU   ',FLUXV1(mx-3:mx+1)
!      WRITE(7,'(A5,i4)') 'mx   ', mx
!      STOP 'comp1'
C
      DO 242 I=2,73
  242 PU(I,JM) = (FLUXV1(I)-VM2)*3.0
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
C
      DO 250 J=3,JMM2
      DO 250 I=2,74
  250 H(I,J)=(T(I-1,J,L)+T(I,J,L))*FD(I,J)
C
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
      H(I,J)=H(I,J)*FD(I,J)
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
 
C     ****************
C     ****************
      SUBROUTINE COMP2
C     ****************
C     ****************
C
C
C             COMPUTES CORIOLIS FORCE,GEOPOTENTIAL GRADIENT AND PRESSURE
C                  GRADIENT
C
C
C***********   BEGINNING OF COMMON ************
C
	include 'recom.fi'
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
      COMMON /WORK/ PV(74,46),FD(74,46),PHI(74,46,2),
     1   H(74,46),PU(74,46),VPS(74,46,2),
     1   TEM2(46),TEM(46),TEMP(74),TEMP2(74),
     1   FDXU(46),FDYU(46),VPK(74),TP1(74),TP3(74),PL1(74),PL3(74),
     1   PLK(74),TPK1(74),TPK3(74)
C
C************ END OF COMMON ******************
C
      DIMENSION FXCDX(46)
      DIMENSION PHIPV(74,46),PHIPU(74,46)
C************************************************************
C             COMPUTE CORIOLIS FORCE
C************************************************************
      JMM1=JM-1
      JMM2=JM-2
      DT8=0.125*DT
      DO 5 J=2,JMM1
      TEM2(J)=0.25*DT8*(DXU(J)-DXU(J+1))
      TEM(J)=F(J)*DXYP(J)*DT8
    5 CONTINUE
C
      FC=0.5
      IF (MRCH.LE.2) FC=0.25
      DO 7 J=1,JM
      FDXU(J)=FC*DT*DXU(J)
    7 FDYU(J)=FC*DT*DYU(J)
      DO 20 L=1,2
C
      DO 10 I=2,73
      FD(I,1)=0.0
   10 FD(I,JM)=0.0
C
      DO 12 J=2,JM
      DO 12 I=2,73
   12 PV(I,J)=U(I,J,L)+U(I-1,J,L)
C
      DO 14 J=2,JMM1
      DO 14 I=2,73
   14 FD(I,J)=TEM(J)+TEM2(J)*(PV(I,J)+PV(I,J+1))
C
      DO 15 J=2,JM
      DO 15 I=2,73
   15 PV(I,J)=(P(I,J)+P(I,J-1))*(FD(I,J)+FD(I,J-1))
      DO 16 J=2,JM
   16 PV(74,J)=PV(2,J)
      DO 18 J=2,JM
      DO 18 I=2,73
   18 PU(I,J)=PV(I,J)+PV(I+1,J)
      DO 20 J=2,JM
      DO 20 I=2,73
      UT(I,J,L)=UT(I,J,L)+PU(I,J)*V(I,J,L)
   20 VT(I,J,L)=VT(I,J,L)-PU(I,J)*U(I,J,L)
C************************************************************
C             COMPUTATION OF GEOPOTENTIAL
C************************************************************
      DO 25 J=1,JM
      DO 25 I=2,73
   25 FD(I,J)=(SM(I,J,1)+SM(I,J,2))*KAPA
C
      DO 45 J=1,JM
      DO 30 I=2,73
      TEMP(I)=SIG1*P(I,J)
      TEMP2(I)=SIG3*P(I,J)
      PL1(I)=TEMP(I)+PTROP
      PL3(I)=TEMP2(I)+PTROP
      TP1(I)=T(I,J,1)*TEMP(I)/PL1(I)
      TP3(I)=T(I,J,2)*TEMP2(I)/PL3(I)
   30 CONTINUE
      DO 35 I=2,73
      PLK(I)=(PL1(I)/PL3(I))**KAPA
      TPK1(I)=0.5*(T(I,J,1)+PLK(I)*T(I,J,2))
      TPK3(I)=0.5*(T(I,J,2)+T(I,J,1)/PLK(I))
   35 CONTINUE
      DO 40 I=2,73
      TT(I,J,1)=TT(I,J,1)+DT*(FD(I,J)*TP1(I)-SD(I,J)*TPK1(I))
   40 TT(I,J,2)=TT(I,J,2)+DT*(FD(I,J)*TP3(I)+SD(I,J)*TPK3(I))
      DO 45 I=2,73
      VPS(I,J,1)=HRGAS*TP1(I)
      VPS(I,J,2)=HRGAS*TP3(I)
      VPK(I)=(TPK3(I)-TPK1(I))*HRGAS/KAPA
      TEMP(I)=VPS(I,J,1)+VPS(I,J,2)+PHIS(I,J)
      PHI(I,J,1)=TEMP(I)+VPK(I)
      PHI(I,J,2)=TEMP(I)-VPK(I)
   45 CONTINUE
C
      DO 46 J=1,JM
      PHI(1,J,1)=PHI(73,J,1)
      PHI(74,J,1)=PHI(2,J,1)
      PHI(1,J,2)=PHI(73,J,2)
      PHI(74,J,2)=PHI(2,J,2)
   46 CONTINUE
C
C        VPS IS HALF WHAT IT SHOULD BE - DOUBLE IT
C
      DO 47 J=1,JM
      DO 47 I=2,73
      VPS(I,J,1)=VPS(I,J,1)+VPS(I,J,1)
      VPS(I,J,2)=VPS(I,J,2)+VPS(I,J,2)
   47 CONTINUE
C*******************************************************************
C      ENERGY CONVERSION TERM IN THERMODYNAMIC ENERGY EQUATION
C*******************************************************************
      FXCO=0.5*FC*DT*KAPA/RGAS
C
      DO 50 J=1,JM
   50 FXCDX(J)=FXCO*DXU(J)
C
      DO 390 L=1,2
C************************************************************
C             GRADIENT OF GEOPOTENTIAL
C************************************************************
      DO 70 J=2,JM
      DO 70 I=2,73
      PU(I,J)=(P(I+1,J)+P(I,J))*(PHI(I+1,J,L)-PHI(I,J,L))
   70 PV(I,J)=(P(I,J)+P(I,J-1))*(PHI(I,J,L)-PHI(I,J-1,L))
C************************************************************
C             GRADIENT OF PRESSURE
C************************************************************
      DO 100 J=1,JM
      VPS(1,J,L)=VPS(73,J,L)
  100 VPS(74,J,L)=VPS(2,J,L)
C
      DO 110 J=2,JM
      DO 110 I=2,73
      PHIPU(I,J)=(VPS(I+1,J,L)+VPS(I,J,L))*(P(I+1,J)-P(I,J))
  110 PHIPV(I,J)=(VPS(I,J,L)+VPS(I,J-1,L))*(P(I,J)-P(I,J-1))
C
      DO 120 J=2,JM
      DO 120 I=2,73
      PV(I,J)=PV(I,J)+PHIPV(I,J)
  120 PU(I,J)=PU(I,J)+PHIPU(I,J)
C
      DO 130 I=2,73
  130 PU(I,1)=0.
C
      CALL AVRX (PU)
C
      DO 140 J=2,JM
      PHIPU(1,J)=PHIPU(73,J)
      PHIPU(74,J)=PHIPU(2,J)
      PHIPV(1,J)=PHIPV(73,J)
      PHIPV(74,J)=PHIPV(2,J)
      PV(1,J)=PV(73,J)
      PV(74,J)=PV(2,J)
      PU(1,J)=PU(73,J)
  140 PU(74,J)=PU(2,J)
C
      IF (MRCH.EQ.3) GO TO 230
      IF (MRCH.EQ.4) GO TO 310
C
C        MRCH = 0 OR 1 OR 2
C
      DO 170 J=2,JM
      DO 170 I=2,73
      VT(I,J,L)=VT(I,J,L)-FDXU(J)*(PV(I,J)+PV(I+1,J))
  170 UT(I,J,L)=UT(I,J,L)-FDYU(J)*(PU(I,J)+PU(I,J-1))
C
      DO 220 J=2,JM
      DO 220 I=2,73
      TEMP(I)=FXCDX(J)*(V(I,J,L)+V(I-1,J,L))*PHIPV(I,J)
      TT(I,J,L)=TT(I,J,L)+TEMP(I)
  220 TT(I,J-1,L)=TT(I,J-1,L)+TEMP(I)
C
      GO TO 350
C
C        MRCH=3.
C
  230 CONTINUE
C
      DO 250 J=2,JM
      DO 250 I=2,73
      VT(I,J,L)=VT(I,J,L)-FDXU(J)*PV(I+1,J)
  250 UT(I,J,L)=UT(I,J,L)-FDYU(J)*PU(I,J)
C
      DO 290 J=2,JM
      DO 290 I=2,73
      TEMP(I)=FXCDX(J)*V(I,J,L)*PHIPV(I,J)
      TT(I,J,L)=TT(I,J,L)+TEMP(I)
  290 TT(I,J-1,L)=TT(I,J-1,L)+TEMP(I)
C
      GO TO 350
C
C        MRCH = 4
C
  310 CONTINUE
C
      DO 330 J=2,JM
      DO 330 I=2,73
      VT(I,J,L)=VT(I,J,L)-FDXU(J)*PV(I,J)
  330 UT(I,J,L)=UT(I,J,L)-FDYU(J)*PU(I,J-1)
C
      DO 340 J=2,JM
      DO 340 I=2,73
      TEMP(I)=FXCDX(J)*V(I-1,J,L)*PHIPV(I,J)
      TT(I,J,L)=TT(I,J,L)+TEMP(I)
  340 TT(I,J-1,L)=TT(I,J-1,L)+TEMP(I)
C
  350 CONTINUE
C
      DO 360 J=2,JMM1
C
      DO 355 I=2,73
  355 PU(I,J)=SAVPU(I,J,L)*PHIPU(I,J)
C
      PU(1,J)=PU(73,J)
C
      DO 360 I=2,73
  360 TT(I,J,L)=TT(I,J,L)+PU(I,J)+PU(I-1,J)
C
  390 CONTINUE
C************************************************************
C             ADJUSTMENT AT THE POLES
C************************************************************
      DO 410 L=1,2
      PB1=SUMma(TT(1,1,L))/FIM
      PB2=SUMma(TT(1,JM,L))/FIM
C
      DO 400 I=2,73
      TT(I,1,L)=PB1
  400 TT(I,JM,L)=PB2
C
      PB1=SUMma(QWT(1,1,L))/FIM
      PB2=SUMma(QWT(1,JM,L))/FIM
C
      DO 410 I=2,73
      QWT(I,1,L)=PB1
  410 QWT(I,JM,L)=PB2
C************************************************************
C             RETURN TO UNWEIGHTED VARIABLES
C************************************************************
      DO 430 J=1,JM
      DO 430 I=2,73
  430 FD(I,J)=PT(I,J)*DXYP(J)
C
      DO 435 J=1,JM
      FD(1,J)=FD(73,J)
  435 FD(74,J)=FD(2,J)
C
      DO 440 L=1,2
      DO 440 J=1,JM
      DO 440 I=2,73
      TT(I,J,L)=TT(I,J,L)/FD(I,J)
  440 QWT(I,J,L)=QWT(I,J,L)/FD(I,J)
C
      DO 460 J=1,JM
      DO 460 I=2,73
  460 PU(I,J)=FD(I,J)+FD(I+1,J)
C
      DO 465 I=2,73
      PU(I,1)=PU(I,1)+PU(I,1)
  465 PU(I,JM)=PU(I,JM)+PU(I,JM)
C
      DO 470 J=2,JM
      DO 470 I=2,73
  470 FD(I,J)=4.0/(PU(I,J)+PU(I,J-1))
C
      DO 480 L=1,2
      DO 480 J=2,JM
      DO 480 I=2,73
      UT(I,J,L)=UT(I,J,L)*FD(I,J)
  480 VT(I,J,L)=VT(I,J,L)*FD(I,J)
      RETURN
      END
C     ****************
C     ****************
      SUBROUTINE COMP4
C     ****************
C     ****************
C
C
C              SMOOTHS LAPSE RATE
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
      COMMON /COMP/ G(74,46),B(74,46),S(74,46)
C
      DIMENSION TDSM(74),TBAR(74)
C
      JMM1=JM-1
      DO 10 J=1,JM
      DO 10 I=2,73
   10 G(I,J)=0.5*(T(I,J,2)-T(I,J,1))/P(I,J)
C
      DO 20 J=1,JM
      G(1,J)=G(73,J)
   20 G(74,J)=G(2,J)
C
      DO 30 J=1,JM
      DO 30 I=2,74
   30 B(I,J)=G(I-1,J)+G(I,J)
C
      DO 40 J=1,JM
      DO 40 I=2,73
   40 S(I,J)=B(I,J)+B(I+1,J)
C
      DO 50 J=2,JM
      DO 50 I=2,73
   50 B(I,J)=S(I,J)+S(I,J-1)
C
      DO 60 J=2,JMM1
      DO 60 I=2,73
   60 S(I,J)=(B(I,J)+B(I,J+1))*0.0625
C
C          S(I,J)= (  G(I-1,J-1)+2.*G(I,J-1)+G(I+1,J-1)
C                   +2.*G(I-1,J)+4.*G(I,J)+2.*G(I+1,J)
C                   + G(I-1,J+1)+2.*G(I,J+1)+G(I+1,J+1) )/16.
C
      DO 200 J=2,JMM1
      DO 100 I=2,73
      TDSM(I)=(G(I,J)+(S(I,J)-G(I,J))/TSPD)*P(I,J)
  100 TBAR(I)=(T(I,J,2)+T(I,J,1))*0.5
C
      DO 200 I=2,73
      T(I,J,1)=TBAR(I)-TDSM(I)
  200 T(I,J,2)=TBAR(I)+TDSM(I)
      RETURN
      END
C     ***************
C     ***************
      FUNCTION SUMma(X)
C     ***************
C     ***************
C
C
C              ZONAL ADDING ROUTINE
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION X(74),S(36),A(18),B(9),C(4),D(2)
C
      DO 10 I=2,37
   10 S(I-1)=X(I)+X(I+36)
      DO 20 I=1,18
   20 A(I)=S(I)+S(I+18)
      DO 30 I=1,9
   30 B(I)=A(I)+A(I+9)
      DO 40 I=1,4
   40 C(I)=B(I)+B(I+4)
      DO 50 I=1,2
   50 D(I)=C(I)+C(I+2)
      SUMma=D(1)+D(2)+B(9)
      RETURN
      END
