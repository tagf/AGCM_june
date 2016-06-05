C     *******************
      SUBROUTINE COMP3(J)
C     *******************
C
C
C             COMPUTES HEATING AND MOISTURE BUDGETS,CONVECTION,SURFACE
C                 WINDS,SURFACE FRICTION AND MOMENTUM REDISTRIBUTION
C
C
C***********************************************************************
C                                                                      *
C     NEW COMP3                                                        *
C     1. VERTICAL FLUX OF MOISTURE HAS BEEN ADDED                      *
C     2. LARGE SCALE PRECIPITATION HAS BEEN MODIFIED AND LARGE         *
C        SCALE EVAPORATION HAS BEEN INCLUDED                           *
C     3. CUMULUS CONVECTION PARAMETERIZATION HAS BEEN MODIFIED         *
C     4. THE RADIATION ROUTINE HAS BEEN MODIFIED AS NOTED THERE        *
C     5. ATMOSPHERIC HEAT AND MOISTURE BUDGETS HAVE BEEN MODIFIED      *
C     6. NEW SOLAR CONSTANT, 1.94 LY/MIN                               *
C     7. GROUND TEMPERATURE CALCULATION HAS BEEN MODIFIED              *
C     8. SNOWFALL, EVAPORATION AND MELTING HAVE BEEN ADDED             *
C     9. GROUND HYDROLOGY HAS BEEN MODIFIED                            *
C    10. SURFACE WIND ON U-V AND PI GRIDS HAS BEEN MODIFIED            *
C                                                                      *
C***********************************************************************
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
C        QT ARRAY - TEMPORARY STORAGE FOR SIMULATION
C
      COMMON /QTARY / PT(74,46),UT(74,46,2),VT(74,46,2),TT(74,46,2)
     1  ,QWT(74,46,2)
C
C        COMMON FOR RADIATION AND CLOUD SUB.  COMP3
C
      LOGICAL SNOW,SEAICE
C
      include 'radvar.fi'
C
      include 'comp.fi'
C
	include 'dir.fi'
C            For ocean
         COMMON /WINDatm/ TAUU(72,46),TAUV(72,46)
c         COMMON /OSA/ QSM(24,16),TXM(24,16),TYM(24,16)                  
c         COMMON /OSA/ QSN(24,16),TXN(24,16),TYN(24,16)
c      include 'paramz.fi'
      COMMON /OSA/ QSN(72,46),TXN(72,46),TYN(72,46)
      COMMON /WIND/ WIND(74,46)
c
      include 'ice.fi'
C************ END OF COMMON ******************
C
      DIMENSION H(74,3),BETADT(74,2)
      DIMENSION SPQ(74),PHIUV(74),ISFUV(74)
      DIMENSION DU13(74),DV13(74),FMCOL(74),COLMRV(74)
      DIMENSION COSF(74)
      DIMENSION STRESS(74), USP(74),VSP(74)
      DIMENSION ROS(74,2),ROUV(74)	
      DIMENSION SUBLIM2(74)		!sea ice
         SQRT(XXXX)=DSQRT(XXXX)
         ALOG10(XXX)=DLOG10(XXX)
         EXP(XXX)=DEXP(XXX)
      DATA ZL3 /2000.0/
10003 CONTINUE
C***********************************************************************
C            PICK UP STATE VARIABLES AND SURFACE CONDITIONS
C***********************************************************************
C
C     ISFTYP- 1=WOODLAND GRASS AND CULTIVATION,2=FOREST,3=STEPPE AND
C          GRASS LAND,4=STEPPE,DESERT,5=DESERT,6=TUNDRA,MOUTAIN AND
C          ARCTIC FLORA,7=WATER(I),8=LAND ICE,9=SEA ICE.
C
      DO 10004 I=2,73
      ISRFCE(I) = ISFTYP(I,J)
      WATER(I) = ISRFCE(I).EQ.7
      SEAICE(I) = ISRFCE(I).EQ.9
      ICE(I) = ISRFCE(I).EQ.8.OR.SEAICE(I)
C***  FOR ICE MODELING
      CHICE(I)  = GHICE(I,J)
      CHICE2(I) = GHICE2(I,J)
      TG1(I)    = GT1(I,J)
      TG2(I)    = GT2(I,J)
      SNR1(I)   = SnoAmt1(I,J)
      SNR2(I)   = SnoAmt2(I,J)
      CIceComp(I)  = GIceComp(I,J)
C
 	if (SEAICE(I)) then
       TG(I) = TG1(I)*CIceComp(I) + TG2(I)*(1.d0-CIceComp(I))
	else
       TG(I) = GT(I,J)
	endif
      WET(I) = GW(I,J)
      SNR(I) = SNOAMT(I,J)
      if (SEAICE(I)) SNR(I) = SNR1(I)  !add 16.03.2001
c!!!!!
      if (WATER(I)) SNR(I) = 0.
c!!!!!
c**************************
c      TG(I) = GT(I,J)
c      WET(I) = GW(I,J)
c      SNR(I) = SNOAMT(I,J)
      TGC(I) = TCICE
      SP(I)=P(I,J)
      SPQ(I)=P(I,J)*0.25
      T1(I)=T(I,J,1)
      T3(I)=T(I,J,2)
      T4(I)=TS(I,J)
      Q1(I) = QW(I,J,1)
      Q3(I) = QW(I,J,2)
10004 FLD(I) = 0.0
C
      DO 10005 I=2,73
      IF( WET(I).LT.1.0E-10 ) WET(I) = 0.0
      IF( SNR(I).GT.500. ) SNR(I) = 500.0  !2002, to stop great snow accum.
      IF (SEAICE(I))  FLD(I) = FLEADS
      IF(J.GT.23.AND.SEAICE(I)) FLD(I)=FLEADN
10005 CONTINUE
C
C*******************************************************************
C               PI AND UV GRID RELATIONSHIP
C
C                I -1                  I                   I +1
C      /-------------------------------------------------------/
C  J+1 / O                    O                    O           /
C      /                                                       /
C      /           +                    +                    + / J +1
C      /                                                       /
C   J  / O                    O                    O           /
C      /                                                       /
C      /           +                    +                    + / J
C      /                                                       /
C  J-1 / O                    O                    O           /
C      /                                                       /
C      /           +                    +                    + / J -1
C      /-------------------------------------------------------/
C       I-1                   I                   I+1
C
C            O (I,J)  IS THE PI GRID
C            + (I ,J ) IS THE UV GRID
C
C      NOTE: U AND V DO NOT EXIST AT J =1
C
C************************************************************
C             FOR SPIN
C************************************************************
      DO 10006 I=1,74
      H(I,1)=H(I,2)
      H(I,2)=H(I,3)
      BETADT(I,1)=BETADT(I,2)
      US(I,1)=US(I,2)
      VS(I,1)=VS(I,2)
      ROS(I,1)=ROS(I,2)
10006 CONTINUE
C
      IF (J.NE.JM) GO TO 10020
C************************************************************
C             NORTH POLE SURFACE WIND
C************************************************************
      DO 10007 I=2,73
      USP(I)=US(I,2)*COSLN(I-1)-VS(I,2)*SINLN(I-1)
10007 VSP(I)=US(I,2)*SINLN(I-1)+VS(I,2)*COSLN(I-1)
      USPI=SUMma(USP)
      VSPI=SUMma(VSP)
      WP=sqrt(USPI*USPI+VSPI*VSPI)/72.0
      DO 10008 I=2,73
10008 WMAGPI(I)=WP
      GO TO 11
10020 CONTINUE
C************************************************************
C             NON-POLAR SURFACE WIND
C************************************************************
      DO 10070 I=2,73
      US(I,2)=1.4*(SIG3*U(I,J+1,2)-SIG1*U(I,J+1,1))
      VS(I,2)=1.4*(SIG3*V(I,J+1,2)-SIG1*V(I,J+1,1))
10070 CONTINUE
      US(1,2)=US(73,2)
      VS(1,2)=VS(73,2)
C
      IF (J.EQ.1) GO TO 10100
C
      CALL ISFCUV (J,ISFTYP(1,J),ISFTYP(1,J-1),ISFUV)
      DO 10080 I=2,73
      PHIUV(I)=DXYU2(J)*(PHIS(I,J)+PHIS(I+1,J))+
     1         DXYU1(J)*(PHIS(I,J-1)+PHIS(I+1,J-1))
      WMAGUV(I)=sqrt(US(I,1)*US(I,1)+VS(I,1)*VS(I,1))
      CDUV(I)=0.002+CFG*dMAX1(PHIUV(I),0.0d0)
      TEMP(I)=dMIN1(0.001d0+7.d-5*WMAGUV(I),0.0025d0)
      USP(I)=US(I,1)+US(I,2)+US(I-1,1)+US(I-1,2)
      VSP(I)=VS(I,1)+VS(I,2)+VS(I-1,1)+VS(I-1,2)
10080 WMAGPI(I)=0.25*sqrt(USP(I)*USP(I)+VSP(I)*VSP(I))
      DO 10085 I=2,73
      IF (ISFUV(I).EQ.7) CDUV(I)=TEMP(I)
10085 CONTINUE
      GO TO 11
10100 CONTINUE
C************************************************************
C             SOUTH POLE SURFACE WIND
C************************************************************
C
C        FOR SPIN
C
      DO 10110 I=1,74
      CDUV(I)=0.0
      WMAGUV(I)=0.0
      US(I,1)=0.0
10110 VS(I,1)=0.0
C
      DO 10120 I=2,73
      USP(I)=US(I,2)*COSLN(I-1)+VS(I,2)*SINLN(I-1)
10120 VSP(I)=-US(I,2)*SINLN(I-1)+VS(I,2)*COSLN(I-1)
      USPI=SUMma(USP)
      VSPI=SUMma(VSP)
      WP=sqrt(USPI*USPI+VSPI*VSPI)/72.0
      DO 10130 I=2,73
10130 WMAGPI(I)=WP
   11 CONTINUE
C************************************************************
C             SET SURFACE WIND PARAMETERS
C************************************************************
C
C        WMAGPI = SURFACE WIND MAGNITUDE ON PI GRID
C        WMAGUV = SURFACE WIND MAGNITUDE ON U,V GRID
C
      DO 12 I=2,73
      WINDF(I)=dMAX1(WMAGPI(I),GUSTY)
      EDV(I)= ED*EXP(0.32E0*(TG(I)-T4(I)))
      CDPI(I)=0.002+CFG*dMAX1(PHIS(I,J),0.0d0)
   12 TEMP(I)=dMIN1(0.001d0+7.d-5*WMAGPI(I),0.0025d0)
C
      DO 13 I=2,73
      IF (WATER(I)) CDPI(I)=TEMP(I)
   13 CONTINUE
C
C        SET MAXIMUM FOR EDV
C
      DO 14 I=2,73
      TEMP(I)=100.
      IF (WATER(I)) TEMP(I)=15.
   14 CONTINUE
C
      DO 15 I=2,73
      EDV(I)=dMIN1(TEMP(I),EDV(I))/ZL3
      DRAW(I)=CDPI(I)*WINDF(I)
   15 CONTINUE
C************************************************************
C             RADIATION CONSTANTS
C************************************************************
      AREA=DXYP(J)
      TOZONE=O3AMT(J)
C
      DO 16 I=2,73
   16 COSF(I)=-COSLN(I-1)*COSR+SINLN(I-1)*SINR
C
      DO 20 I=2,73
      COSZ(I)=SINL(J)*SIND+COSL(J)*COSD*COSF(I)
   20 SCOSZ(I)=S0*COSZ(I)
C
      DO 25 I=2,73
      IF (COSZ(I).LE.0.01) SCOSZ(I)=0.0
   25 CONTINUE
C************************************************************
C             PRESSURES
C************************************************************
      DO 50 I=2,73
      COLMR(I)=PM/SP(I)
      P4(I)=SP(I)+PTROP
      P4L(I)=ALOG10(P4(I))
   50 P4K(I)=EXP(P4L(I)*KAPEL)
      DO 51 I=2,73
      PL1(I)=SPQ(I)+PTROP
      P1L(I)=ALOG10(PL1(I))
   51 PL1K(I)=EXP(P1L(I)*KAPEL)
      DO 52 I=2,73
      PL2(I)=SPQ(I)+SPQ(I)+PTROP
      P2L(I)=ALOG10(PL2(I))
   52 PL2K(I)=EXP(P2L(I)*KAPEL)
      DO 53 I=2,73
      PL3(I)=SP(I)-SPQ(I)+PTROP
      P3L(I)=ALOG10(PL3(I))
   53 PL3K(I)=EXP(P3L(I)*KAPEL)
C*********************************************************************
C         TEMPERATURES AND TEST FOR DRY-ADIABATIC INSTABILITY
C*********************************************************************
      DO 60 I=2,73
      TEM(I)=(T1(I)+T3(I))/(PL1K(I)+PL3K(I))
      THL1(I)=T1(I)/PL1K(I)
   60 THL3(I)=T3(I)/PL3K(I)
      DO 65 I=2,73
      IF (THL1(I) .GT. THL3(I)) GO TO 65
      T(I,J,1)=TEM(I)*PL1K(I)
      T(I,J,2)=TEM(I)*PL3K(I)
      T1(I)=T(I,J,1)
      T3(I)=T(I,J,2)
      THL1(I)=TEM(I)
      THL3(I)=TEM(I)
   65 CONTINUE
C************************************************************
C             MOISTURE VARIABLES
C************************************************************
      CALL QSDQV(T1,PL1,QS1,FSGAM)
C
      DO 67 I=2,73
   67 GAM1(I)=CLH*FSGAM(I)
C
      CALL QSDQV(T3,PL3,QS3,FSGAM)
C
      DO 70 I=2,73
      GAM3(I)=CLH*FSGAM(I)
      RH1(I) = Q1(I)/QS1(I)
   70 RH3(I) = Q3(I)/QS3(I)
C**********************************************************************
C       TEMPERATURE EXTRAPOLATION AND INTERPOLATION FOR RADIATION
C**********************************************************************
      DO 80 I=2,73
      TEMA(I)=(THL3(I)-THL1(I))/(PL3K(I)-PL1K(I))
      TEMB(I)=(THL1(I)*PL3K(I)-THL3(I)*PL1K(I))/(PL3K(I)-PL1K(I))
      T2(I)=(TEMA(I)*PL2K(I)+TEMB(I))*PL2K(I)
      TTROP(I)=(TEMA(I)*PTRK+TEMB(I))*PTRK
   80 CONTINUE
C
      CALL QSATV (T2,PL2,QS2)
C************************************************************
C             VERTICAL FLUX OF MOISTURE
C************************************************************
      DO 90 I=2,73
      SDP(I)=DTC3*SD(I,J)/(SP(I)*AREA)
C
C        COMPUTE Q2 BY ASSUMING Q = LINEAR FUNCTION OF QS
C
      CI1(I)=dMAX1(QS3(I)-QS2(I),1.d-6)
      CI3(I)=dMAX1(QS2(I)-QS1(I),1.d-6)
   90 Q2(I)=(CI1(I)*Q1(I)+CI3(I)*Q3(I))/(CI1(I)+CI3(I))
      DO 92 I=2,73
      IF(SD(I,J).LT.0.0.AND.RH3(I).LT.RH1(I)) Q2(I)=QS2(I)*RH3(I)
      IF(SD(I,J).GT.0.0.AND.RH3(I).GT.RH1(I)) Q2(I)=QS2(I)*RH1(I)
   92 CONTINUE
      DO 94 I=2,73
      Q2C(I)=Q2(I)
      FLUXQ(I)=SDP(I)*Q2(I)
   94 TEM(I)=0.5E0*(THL1(I)-THL3(I))*PL2K(I)/CLH
      DO 96 I=2,73
      IF (FLUXQ(I).GE.0.0) GO TO 96
      IF (Q3(I)+FLUXQ(I).GT.QS3(I))
     1  Q2C(I)=dMAX1(Q2(I),QS3(I)-TEM(I))
      IF (Q1(I)-FLUXQ(I).GT.QS1(I))
     1  Q2C(I)=dMIN1(Q2(I),QS1(I)+TEM(I))
   96 CONTINUE
      DO 98 I=2,73
      Q2(I)=Q2C(I)
      FLUXQ(I)=SDP(I)*Q2(I)
      Q1(I) = Q1(I)-FLUXQ(I)
      Q3(I) = Q3(I)+FLUXQ(I)
      RH1(I) = Q1(I)/QS1(I)
   98 RH3(I) = Q3(I)/QS3(I)
C********************************************************************
C             LARGE SCALE PRECIPITATION AND EVAPORATION
C********************************************************************
      DO 100 I=2,73
      PREC1(I) = 0.0
      PREC3(I) = 0.0
      IF (RH1(I).LE.1.0.AND.RH3(I).LE.1.0) GO TO 100
      IF( RH1(I).LE.1.0 ) GO TO 99
C
C        MOISTURE IS TO BE MOVED FROM LEVEL 1 TO LEVEL 3 - RAIN FROM
C             LEVEL 1 - PREC1 IS EVAPORATED TO LEVEL 3.
C
      EX(I) = (Q1(I)-QS1(I))/(1.0+GAM1(I))
      T1(I) = T1(I)+CLH*EX(I)
      Q1(I) = Q1(I)-EX(I)
      PREC1(I) = PREC1(I)+EX(I)
      QS1(I) = QSDQS(T1(I),PL1(I),FSGAM(I))
      GAM1(I) = CLH*FSGAM(I)
C
      EX(I) = (Q1(I)-QS1(I))/(1.0+GAM1(I))
      T1(I) = T1(I)+CLH*EX(I)
      Q1(I) = Q1(I)-EX(I)
      PREC1(I) = PREC1(I)+EX(I)
      QS1(I) = QSDQS(T1(I),PL1(I),FSGAM(I))
      GAM1(I)=CLH*FSGAM(I)
C
      THL1(I) = T1(I)/PL1K(I)
      RH1(I) = Q1(I)/QS1(I)
      Q3(I) = Q3(I)+PREC1(I)
      T3(I) = T3(I)-CLH*PREC1(I)
      THL3(I) = T3(I)/PL3K(I)
      QS3(I) = QSDQS(T3(I),PL3(I),FSGAM(I))
      RH3(I) = Q3(I)/QS3(I)
      GAM3(I)=CLH*FSGAM(I)
C
   99 IF(RH3(I).LE.1.0) GO TO 105
      EX(I) = (Q3(I)-QS3(I))/(1.0+GAM3(I))
      T3(I) = T3(I) + CLH*EX(I)
      Q3(I) = Q3(I)-EX(I)
      PREC3(I) = PREC3(I)+EX(I)
      QS3(I) = QSDQS(T3(I),PL3(I),FSGAM(I))
      GAM3(I) = CLH*FSGAM(I)
C
      EX(I) = (Q3(I)-QS3(I))/(1.0+GAM3(I))
      T3(I) = T3(I) + CLH*EX(I)
      Q3(I) = Q3(I)-EX(I)
      PREC3(I) = PREC3(I)+EX(I)
      QS3(I) = QSDQS(T3(I),PL3(I),FSGAM(I))
      GAM3(I)=CLH*FSGAM(I)
C
      THL3(I) = T3(I)/PL3K(I)
      RH3(I) = Q3(I)/QS3(I)
  105 CONTINUE
      TEMA(I)=(THL3(I)-THL1(I))/(PL3K(I)-PL1K(I))
      TEMB(I)=(THL1(I)*PL3K(I)-THL3(I)*PL1K(I))/(PL3K(I)-PL1K(I))
      T2(I)=(TEMA(I)*PL2K(I)+TEMB(I))*PL2K(I)
      TTROP(I)=(TEMA(I)*PTRK+TEMB(I))*PTRK
      QS2(I) = QSAT(T2(I),PL2(I))
      CI1(I) = dMAX1(QS3(I)-QS2(I),1.d-6)
      CI3(I) = dMAX1(QS2(I)-QS1(I),1.d-6)
      Q2(I)=(CI1(I)*Q1(I)+CI3(I)*Q3(I))/(CI1(I)+CI3(I))
      IF(SD(I,J).LT.0.0.AND.RH3(I).LT.RH1(I)) Q2(I)=QS2(I)*RH3(I)
      IF(SD(I,J).GT.0.0.AND.RH3(I).GT.RH1(I)) Q2(I)=QS2(I)*RH1(I)
  100 CONTINUE
C************************************************************
C             PREPARATION FOR CONVECTION
C************************************************************
      DO 110 I=2,73
      DTHL2(I)=0.5*(THL1(I)-THL3(I))*PL2K(I)
      SS3(I)=THL3(I)*P4K(I)
      SS2(I) = SS3(I) + DTHL2(I)
      SS1(I) = SS2(I) + DTHL2(I)
      HH3(I) = SS3(I) + CLH*Q3(I)
      HH3S(I) = SS3(I) + CLH*QS3(I)
  110 HH1S(I) = SS1(I) + CLH*QS1(I)
C************************************************************
C             MIDDLE LEVEL CONVECTION
C************************************************************
      DO 120 I=2,73
      CT1(I) = 0.
      CT3(I) = 0.
      CQ1(I) = 0.
      CQ3(I) = 0.
      Q1P(I)=Q1(I)-1.E-9
      Q3P(I)=Q3(I)-1.E-9
      TCNV(I) = 0.0
      BETADT(I,2)=0.0
      TCNV(I)=DTC3/3600.
  120 EX(I) = HH3(I) - HH1S(I)
C
      DO 125 I=2,73
      IF (EX(I).LE.0.)  GO TO 125
      TEMA(I)=(1.+GAM1(I))
      TEMB(I) = EX(I)*TCNV(I)/(EX(I) +
     1 GAM1(I)*DTHL2(I)+CLH*(Q3(I)-Q2(I)))
      DQ1(I)=GAM1(I)*EX(I)/(CLH*TEMA(I))+QS1(I)-Q2(I)
      DQ3(I)=Q2(I)-Q3(I)
      IF (TEMB(I)*DQ1(I).LT.-Q1P(I)) TEMB(I)=-Q1P(I)/DQ1(I)
      IF (TEMB(I)*DQ3(I).LT.-Q3P(I)) TEMB(I)=-Q3P(I)/DQ3(I)
      TEMB(I)=dMIN1(TEMB(I),1.0d0)
      CT1(I)=TEMB(I)*(EX(I) + TEMA(I)*DTHL2(I))/TEMA(I)
      CT3(I) = TEMB(I)*DTHL2(I)
      CQ1(I) = TEMB(I)*DQ1(I)
      CQ3(I) = TEMB(I)*DQ3(I)
      BETADT(I,2)=0.5*TEMB(I)
  125 CONTINUE
C************************************************************
C             DETERMINATION OF SURFACE TEMPERATURE
C************************************************************
      CALL QSDQV(TG,P4,QG,DQG)
C
      DO 140 I=2,73
      EVE(I)=dMIN1(1.0d0,WET(I)+WET(I))
      EVE(I)=CVMGP(EVE(I),1.d0,-SNR(I))
      EVE(I)=CVMGP(EVE(I),1.d0,QG(I)-Q3(I))
  140 CONTINUE
C
      DO 143 I=2,73
      EVACO(I)=EVE(I)*DRAW(I)
  143 Q4(I)=(EVACO(I)*QG(I)+EDV(I)*Q3(I))/(EVACO(I)+EDV(I))
      CALL QSDQV(T4,P4,QS4,FSGAM)
      DO 145 I=2,73
      GAM4(I)=CLH*FSGAM(I)
      RH4(I)=Q4(I)/QS4(I)
      T4CR(I)=THL3(I)*P4K(I)-RH4(I)*ZL3*GRAV*
     1        (KAPA*GAM4(I)-CLH*QS4(I)/T4(I))/
     2        (RGAS*(1.0+GAM4(I)))
      T4(I)=(DRAW(I)*TG(I)+EDV(I)*T4CR(I))/(DRAW(I)+EDV(I))
  145 HH4(I)=T4(I)+CLH*Q4(I)
C************************************************************
C             PENETRATING AND LOW-LEVEL CONVECTION
C************************************************************
      DO 150 I=2,73
      PCT1(I)=0.
      PCT3(I)=0.
      PCQ1(I)=0.
      PCQ3(I)=0.
  150 EX(I)=0.
      DO 155 I=2,73
      IF (HH4(I) .LE. HH3S(I) .OR. HH3(I).GT.HH1S(I)) GO TO 155
      EX(I) = HH4(I)-HH3S(I)
      IF (HH4(I) .LT. HH1S(I)) GO TO 153
      TEMB(I) = CLH*(Q4(I)-Q2(I))*EDV(I)/(EDV(I)+EVACO(I))
     *   +(SS2(I)-T4(I))*(1.E0+GAM3(I)-EDV(I)*P4K(I)/
     *    (PL3K(I)*(EDV(I)+DRAW(I))))
      IF (TEMB(I).LT.0.001E0) TEMB(I)=0.001E0
      TEMB(I)=TCNV(I)*EX(I)/TEMB(I)
      TEMA(I)=(HH4(I)-HH1S(I))/(1.0+GAM1(I))
      DQ1(I)=GAM1(I)*TEMA(I)/CLH+QS1(I)-Q2(I)
      DQ3(I)=Q2(I)-Q4(I)
      IF (TEMB(I)*DQ1(I).LT.-Q1P(I)) TEMB(I)=-Q1P(I)/DQ1(I)
      IF (TEMB(I)*DQ3(I).LT.-Q3P(I)) TEMB(I)=-Q3P(I)/DQ3(I)
      TEMB(I)=dMIN1(TEMB(I),0.5d0)
      BETADT(I,2)=TEMB(I)
      PCT1(I) = TEMB(I)*(TEMA(I)+DTHL2(I))
      PCT3(I) = TEMB(I)*(SS2(I)-T4(I))
      PCQ1(I) = TEMB(I)*DQ1(I)
      PCQ3(I) = TEMB(I)*DQ3(I)
  153 HH4(I) = HH3S(I)
      QS4(I)=QSDQS(T4(I),P4(I),FSGAM(I))
      GAM4(I)=CLH*FSGAM(I)
      RH4(I)=Q4(I)/QS4(I)
      T4(I)=T4(I)-EX(I)/(1.0+RH4(I)*GAM4(I))
      Q4(I) = (HH4(I)-T4(I))/CLH
  155 CONTINUE
C********************************************************************
C             SURFACE EVAPORATION AND SENSIBLE HEAT FLUX
C********************************************************************
      DO 160 I=2,73
      RO4(I) = 0.1E0*P4(I)/(RGAS*T4(I))
      ROS(I,2)=RO4(I)
      CEVA(I)=CDPI(I)*100.*RO4(I)*WINDF(I)
      CSEN(I)=.24E0*24.*3600.*CEVA(I)
C
C        SENSIBLE HEAT (LY/DAY) AND EVAPORATION (GM/CM**2/SEC)
C
      CEVA(I) = EVE(I)*CEVA(I)
      E4(I) = CEVA(I)*(QG(I)-Q4(I))
      E4DT(I) = E4(I)*DTC3
      F4(I) = CSEN(I)*(TG(I)-T4(I))
c##########  0.9*SENSIBLE HEAT - for seaice ################## 
c      IF (SEAICE(I)) F4(I) = 0.9*F4(I)

      E4C(I) = 0.0E0
  160 F4C(I) = 0.0E0
      DO 162 I=2,73
C***********************
C      FOR ICE MODELLING
C***********************
C     THIS IS A CHECK FOR EXISTANCE IF ICE ???
C      IF(FLD(I).LE.0.0) GO TO 162
      IF (SEAICE(I)) go to 162
C***********************
      QGC(I)=QSAT(TGC(I),P4(I))
      E4C(I) = CEVA(I)*(QGC(I)-Q4(I))
      F4C(I) = CSEN(I)*(TGC(I)-T4(I))
  162 CONTINUE
      DO 164 I=2,73
      CIMELT(I)=0.0
      RAIN(I) = 0.E0
      SNFAL(I) = 0.E0
      SMELT(I) = 0.E0
      SNOW(I) = .FALSE.
      !Total precipitation P*dtc3 [g/cm*2] 
  164 TOTALP(I) = (PREC3(I)-(CQ1(I)+CQ3(I)+PCQ1(I)+PCQ3(I)))*5.*SP(I)
     *   /9.8E0
C************************************************************
C             SNOWFALL AND SNOW EVAPORATION
C************************************************************
      DO 166 I=2,73
      TEM(I)=0.0
C***********************
C      FOR ICE MODELLING
C***********************
c      IF(WATER(I)) GO TO 166
      IF(WATER(I).and.(.not.SEAICE(I))) GO TO 166
C*********************** if land or sea ice:
      IF(T4(I).LE.TICE) SNFAL(I) = TOTALP(I)
      IF(T4(I).GT.TICE) RAIN(I) = TOTALP(I)
  166 CONTINUE  !if any point:
      DO 167 I=2,73
      SNR(I) = SNR(I) + SNFAL(I)
      SNR2(I) = SNR2(I) + SNFAL(I)   !sea ice
      SUBLIM(I)=dMIN1(E4DT(I),SNR(I))
      SUBLIM2(I)=dMIN1(E4DT(I),SNR2(I))  !sea ice
  167 CONTINUE
C
C        EVAPORATE SNOW
C
      DO 168 I=2,73
      SUBLIM(I)=CVMGP(0.d0,SUBLIM(I),-SNR(I))
      SUBLIM2(I)=CVMGP(0.d0,SUBLIM2(I),-SNR2(I))	!sea ice
      SNR(I)=SNR(I)-SUBLIM(I)
      SNR2(I)=SNR2(I)-SUBLIM2(I)	!sea ice
      E4DT(I)=E4DT(I)-SUBLIM(I)
  168 SNOW(I)=SNR(I).GT.0.0
C
C        E4DT IS RESIDUAL EVAPORATION IF SNOW WAS PRESENT
C
      DO 169 I=2,73
      ISF=ISRFCE(I)
      SNR(I)=dMAX1(SNR(I),0.0d0)
      SNR2(I)=dMAX1(SNR2(I),0.0d0)   
C***********************
C      FOR ICE MODELLING
C***********************
c      ALS(I)=SFCALB(ISF,1)+sqrt(DMIN1(SNR(I),1.0D0))
c     1    *(SFCALB(ISF,2)-SFCALB(ISF,1))
	SNR1(i) =SNR(i)		!thick sea ice

c################### 
	if (ISF.ne.9) then
       ALS(I) = SFCALB(ISF,1) + SQRT(DMIN1(SNR(I),1.0D0))
     1    *(SFCALB(ISF,2)-SFCALB(ISF,1))
	else
	 if (CHICE(I).gt.0.) then
        CIceAlb(I)  = SFCALB(9,1) + SQRT(DMIN1(SNR1(I),1.0D0))
     1    *(SFCALB(9,2)-SFCALB(9,1))
	 else
	  CIceAlb(I)  = SFCALB(7,1)
	 endif
	  if (cHIce2(I).gt.0.) then
         CIce2Alb(I) = SFCALB(9,1) + SQRT(DMIN1(SNR2(I),1.0D0))
     1    *(SFCALB(9,2)-SFCALB(9,1))
	  else
	   CIce2Alb(I) = SFCALB(7,1)
	  endif
	 ALS(I) = CIceAlb(I)*CIceComp(I) + CIce2Alb(I)*(1-CIceComp(I))
c	 CIceAlb(I)   = (1d0-CIceAlb(I))/(1d0-Albedo)
c	 CIce2Alb(I)  = (1d0-CIce2Alb(I))/(1d0-Albedo)
	endif
c################### 

C***********************
  169 CONTINUE
C
      CALL CLOUDS
      CALL VAPOR
      CALL FMAFMB
      CALL TRA
      CALL LONGW
      CALL SOLAR
C******************************************************************
C             GROUND TEMPERATURE AND GROUND HYDROLOGY
C******************************************************************
C        GROUND TEMPERATURE CALCULATED BY BULK METHOD(AFTER
C        KATAYAMA) USING ONE SKIN DEPTH. THERMAL PROPERTIES
C        FOR LAND FROM BHUMRALKAR, R-1511-ARPA(PP20-21,33).
C        THERMAL PROPERTIES FOR SNOW AND ICE FROM PRIESTLEY
C        (1959), TURBULENT TRANSFER IN THE LOWER ATMOSPHERE
C        P 100
C
C
C        BULK HEAT CAPACITY
C
      DO 170 I=2,73
      CNDTVY(I)=dMIN1(0.0012d0*(1.0d0+20.d0*WET(I)),0.006d0)
  170 HCAPTY(I)=0.6
C
C***********************
C      FOR ICE MODELLING
C***********************
c         DO 171 I=2,25
c  171    IF (WATER(I).AND.SWAMP) HCAPTY(I)=1.
      do I=2,73
       if (WATER(I).and.(.not.SEAICE(I))) HCAPTY(I)=1.
      enddo
C***********************

      DO 172 I=2,73
      IF (ICE(I))CNDTVY(I)=5.4E-3
  172 CONTINUE
C
      DO 173 I=2,73
      TEM(I)=sqrt(2.0*CNDTVY(I)*HCAPTY(I)/FRDAY)
  173 CZH(I)=TEM(I)
C
      DO 174 I=2,73
      EVAL(I)=600.  ! evaporation latent heat cal/g
      IF (WATER(I)) GO TO 174
      IF (ICE(I)) CZH(I)=8.29  !bulk ground heat cap. ly/deg
      IF (SNOW(I)) CZH(I)=1.16
      IF (ICE(I).OR.SNOW(I)) EVAL(I)=680. !sublimation
  174 CONTINUE
C
C        GROUND TEMPERATURE CHANGE
C
C***********************
C      FOR ICE MODELLING
C***********************
c      DO 175 I=2,73
c      TEMB(I)=0.0
c      IF (SEAICE(I)) TEMB(I)=CNDICE
c  175 CONTINUE
c      DO 176 I=2,73
c      CICE(I)=TEMB(I)*(TG(I)-TCICE)
c  176 DTG(I)=(S4(I)-R4(I)-EVAL(I)*E4(I)*DAY-F4(I)-CICE(I))
c     1  /(CZH(I)*TSPD+4.*STBO*TG(I)*TG(I)*TG(I)+
c     2   EVAL(I)*CEVA(I)*DQG(I)*DAY+CSEN(I)+TEMB(I))

      do i=2,73
c       if (I.eq.Ibreak.and.J.eq.Jbreak) then
c        tmpdata = 0
c       endif
       if (.not.SEAICE(I)) then
        TEMB(I) = 0.0d0
        CICE(I) = 0.0d0
       else
C       TEMB(I) = Ks*Ki/( Ks*HICE(I) + Ki*HSNOW(i) )
        TEMB(I)  = DAY*CTYICE*CTYSN/(CTYSN*CHIce(I)
	1      + CTYICE*SNR1(I)/ROSN)
        CICE(I)  = TEMB(I)*(TG1(I)-TCICE)
	  if (CHIce2(I).ge.HIceMin) then
         TEMB2(I) = DAY*CTYICE*CTYSN/(CTYSN*CHIce2(I)
	1      + CTYICE*SNR2(I)/ROSN)
         CICE2(I) = TEMB2(I)*(TG2(I)-TCICE)
	  else
         TEMB2(I) = 0.0d0
         CICE2(I) = 0.0d0
	  endif
       endif
       if (.not.SEAICE(I)) then
        DTG(I)=(S4(I)-R4(I)-EVAL(I)*E4(I)*DAY-F4(I)-CICE(I))
     1   /(CZH(I)*TSPD+4.*STBO*TG(I)*TG(I)*TG(I)+
     2   EVAL(I)*CEVA(I)*DQG(I)*DAY+CSEN(I)+TEMB(I))
       else
        DTG(I)=((S4(I)+DR4(I))
	1   - STBO*TG1(I)**4
     2   - EVAL(I)*E4(I)*DAY-F4(I)-CICE(I) )
     3   /(CZH(I)*TSPD+4.*STBO*TG1(I)*TG1(I)*TG1(I)+
     4   EVAL(I)*CEVA(I)*DQG(I)*DAY+CSEN(I)+TEMB(I))
        DTG2(I)=(( S4(I)+DR4(I))
	1   - STBO*TG2(I)**4
     2   - EVAL(I)*E4(I)*DAY-F4(I)-CICE(I) )
     3   /(CZH(I)*TSPD+4.*STBO*TG2(I)*TG2(I)*TG2(I)+
     4   EVAL(I)*CEVA(I)*DQG(I)*DAY+CSEN(I)+TEMB2(I))
        GIceCND(I,J)  = GIceCND(I,J) + DMIN1(CICE(I),0d0)/TSPD
        GIce2CND(I,J) = GIce2CND(I,J) + DMIN1(CICE2(I),0d0)/TSPD
       endif
      enddo
      DO 177 I=2,73
      IF (WATER(I)) DTG(I)=0.0
  177 CONTINUE
C***********************
C
C        SNOW MELTING
C
C***********************
C      FOR ICE MODELLING
C***********************
      DO 180 I=2,73
c       if (I.eq.Ibreak.and.J.eq.Jbreak) then
c        tmpdata = 0
c       endif
      IF (.NOT.(SNOW(I).OR.ICE(I))) GO TO 180
c      IF (TG(I)+DTG(I).LE.TICE) GO TO 180
c      HEATNG(I)=(S4(I)+DR4(I)-STBO*TICE*TICE*TICE*TICE
c     1         -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
c     2         -CSEN(I)*(TICE-T4(I))-TEMB(I)*(TICE-TCICE) ) /TSPD
c     3         -CZH(I)*(TICE-TG(I))
c      TG(I)=TICE
c      DTG(I)=0.0
c      IF (HEATNG(I).LE.0.0) GO TO 180
c      SMELT(I)=dMIN1(HEATNG(I)/SNOWL,SNR(I))
c      HEATNG(I)=HEATNG(I)-SMELT(I)*SNOWL
c      SNR(I)=SNR(I)-SMELT(I)
c      SNOW(I)=SNR(I).GT.0.0
c      IF (SNOW(I)) GO TO 180
C        EXISTING SNOW MELTED. ONLY BARE LAND OR ICE REMAINS
c      DTG(I)=HEATNG(I)/TEM(I)
c      IF (.NOT.ICE(I)) GO TO 180
C        ONLY FOR BARE ICE
c      CIMELT(I)=HEATNG(I)/SNOWL
c      DTG(I)=0.0

	if (SEAICE(I)) then
	 if (SNR1(I).lt.1.d-10) then   !толстый лед
        if (TG1(I)+DTG(I).gt.TICE) then
           HEATNG(I)=((S4(I)*(1.d0-GIceSThruPercent)	
	1       +DR4(I))
	2       -STBO*TICE*TICE*TICE*TICE
     3       -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     4       -CSEN(I)*(TICE-T4(I)) ) / TSPD
     5       -CZH(I)*(TICE-TG1(I))
	     GIcePool(I,J)   = GIcePool(I,J)
	1       + S4(I)*GIceSThruPercent/TSPD 
	     GIceHeatng(I,J) = GIceHeatng(I,J) + DMAX1(HEATNG(I),0d0)
	     TG1(I)   = TICE
	     DTG(I)   = 0d0
	  else ! TG1(I)+DTG(I) <= TICE
	   if (GIcePool(I,J).gt.0d0) then
            HEATNG(I)=((S4(I)+DR4(I))
	1		-STBO*(TG1(I)+DTG(I))**4
     2		-EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     3		-CSEN(I)*(TICE-T4(I)) ) / TSPD
     4		-CZH(I)*(TICE-TG1(I))
	        GIcePool(I,J) = DMAX1(GIcePool(I,J) + HEATNG(I), 0d0)
	      TG1(I)   = TICE
	      DTG(I)   = 0d0
	   else ! ничего не мерзнет и не тает
	   endif
	  endif
        CIMELT(I) = 0d0
	 else ! снег на толстом льду
        HEATNG(I)=((S4(I)+DR4(I))
	1         -STBO*TICE*TICE*TICE*TICE
     2         -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     3         -CSEN(I)*(TICE-T4(I))-TEMB(I)*(TICE-TCICE) ) / TSPD
     4         -CZH(I)*(TICE-TG1(I))
        if (TG1(I)+DTG(I).gt.TICE) then
         TG1(I)=TICE
         DTG(I)=0.0
	   if (HEATNG(I).gt.0d0) then
          SMELT(I)  = DMIN1(HEATNG(I)/SNOWL, SNR(I))
          HEATNG(I) = HEATNG(I) - SMELT(I)*SNOWL
          SNR1(I)   = SNR1(I) - SMELT(I)
          SNOW(I)   = SNR1(I).GT.0.0
	   endif
	  endif
	 endif ! конец толстого льда
	 if (CHIce2(I).ge.HIceMin) then  !тонкий лед
	  if (SNR2(I).lt.1d-10) then
         if (TG2(I)+DTG2(I).gt.TICE) then
          HEATNG(I)=(( S4(I)*(1.d0-GIceSThruPercent)  
	1       +DR4(I))
	2       -STBO*TICE*TICE*TICE*TICE
     3       -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     4       -CSEN(I)*(TICE-T4(I)) ) / TSPD
     5       -CZH(I)*(TICE-TG2(I))
	    GIce2Pool(I,J)   = GIce2Pool(I,J)
	1       + S4(I)*GIceSThruPercent/TSPD	
	    GIce2Heatng(I,J) = GIce2Heatng(I,J) + DMAX1(HEATNG(I),0d0)
	    TG2(I)    = TICE
	    DTG2(I)   = 0d0
	   else ! TG2(I)+DTG2(I) <= TICE
	    if (GIce2Pool(I,J).gt.0d0) then
           HEATNG(I)=((S4(I)+DR4(I))
	1		-STBO*(TG2(I)+DTG2(I))**4
     2		-EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     3		-CSEN(I)*(TICE-T4(I)) ) / TSPD
     4		-CZH(I)*(TICE-TG2(I))
	     GIce2Pool(I,J) = DMAX1(GIce2Pool(I,J) + HEATNG(I), 0d0)
	     TG2(I)    = TICE
	     DTG2(I)   = 0d0
	    else ! ничего не мерзнет и не тает
	    endif
	   endif
         CIMELT(I) = 0d0
	 else ! снег на тонком льду
         HEATNG(I)=((S4(I)+DR4(I))
	1         -STBO*TICE*TICE*TICE*TICE
     2         -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     3         -CSEN(I)*(TICE-T4(I))-TEMB(I)*(TICE-TCICE) ) / TSPD
     4         -CZH(I)*(TICE-TG2(I))
         if (TG2(I)+DTG2(I).gt.TICE) then
          TG2(I)=TICE
          DTG2(I)=0.0
	    if (HEATNG(I).gt.0d0) then
           SMELT(I)  = DMIN1(HEATNG(I)/SNOWL, SNR2(I))
           HEATNG(I) = HEATNG(I) - SMELT(I)*SNOWL
           SNR2(I)   = SNR2(I) - SMELT(I)
	    endif
	   endif
	  endif
	 else ! чистая вода в обл тонкого льда
        HEATNG(I)=((S4(I)+DR4(I))
	1         -STBO*TICE*TICE*TICE*TICE
     2         -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     3         -CSEN(I)*(TICE-T4(I))-TEMB(I)*(TICE-TCICE) ) / TSPD
     4         -CZH(I)*(TICE-TG2(I))
	  GIce2Heatng(I,J) = GIce2Heatng(I,J) + DMAX1(HEATNG(I),0d0)
	  TG2(I)    = TICE
	  DTG2(I)   = 0d0
	 endif ! конец тонкого льда
	else ! лед на суше
       HEATNG(I)=(S4(I)+DR4(I)-STBO*TICE*TICE*TICE*TICE
     1         -EVAL(I)*CEVA(I)*(QSAT(TICE,P4(I))-Q4(I))*DAY
     2         -CSEN(I)*(TICE-T4(I))-TEMB(I)*(TICE-TCICE) ) / TSPD
     3         -CZH(I)*(TICE-TG(I))
       IF (TG(I)+DTG(I).LE.TICE) GO TO 179
       TG(I)=TICE
       DTG(I)=0.0
       IF (HEATNG(I).LE.0.0) GO TO 179
       SMELT(I)=DMIN1(HEATNG(I)/SNOWL,SNR(I))
       HEATNG(I)=HEATNG(I)-SMELT(I)*SNOWL
       SNR(I)=SNR(I)-SMELT(I)
       SNOW(I)=SNR(I).GT.0.0
       IF (SNOW(I)) GO TO 179
C        EXISTING SNOW MELTED. ONLY BARE LAND OR ICE REMAINS
       DTG(I)=HEATNG(I)/TEM(I)
       IF (.NOT.ICE(I)) GO TO 179
C        ONLY FOR BARE ICE
       CIMELT(I)=HEATNG(I)/SNOWL
       DTG(I)=0.d0
  179	endif
  180  CONTINUE !конец цикла для точек снега или льда
!      OPEN(66,FILE=TRIM(BaseDir)//WorkDir//'\c3',ACCESS='APPEND')
!      WRITE (66,84) j,HEATNG 
!      WRITE (66,84) j,SNR
!      WRITE (66,84) j,SMELT
!      WRITE (66,84) j,RUNOFF
!      WRITE (66,84) j,E4DT 
!  84  FORMAT (1x,I3,(10d12.4))
C***********************
C
C        UPDATE GROUND TEMPERATURE AND WETNESS
C
      DO 188 I=2,73
      TG(I)=TG(I)+DTG(I)
C***********************
C      FOR ICE MODELLING
C***********************
      if (SEAICE(I)) then
	 TG1(I) = TG1(I) + dTG(I)
	 TG2(I) = TG2(I) + dTG2(I)
       TG(I) = TG1(I)*CIceComp(I) + TG2(I)*(1d0-CIceComp(I))
       SNR(I) = SNR1(I)*CIceComp(I) + SNR2(I)*(1d0-CIceComp(I))
	endif
C***********************
      RUNOFF(I)=WET(I)*(RAIN(I)+SMELT(I))
      WET(I)=WET(I)+(RAIN(I)+SMELT(I)-RUNOFF(I)-E4DT(I))/GWM
  188 CONTINUE
C
C        ADJUST GROUND WETNESS
C
      DO 195 I=2,73
      RUNOFF(I)=RUNOFF(I)+GWM*dMAX1(0.0d0,WET(I)-1.0d0)
      WET(I)=dMIN1(1.0d0,WET(I))
  195 CONTINUE
      DO 200 I=2,73
      IF( WET(I).LT.1.0E-10 ) WET(I) = 0.0E0
      IF (.NOT.(WATER(I).OR.ICE(I))) GO TO 200
      WET(I)=1.0
      RUNOFF(I)=0.0
  200 CONTINUE
C
      DO 210 I=2,73
      GW(I,J) = WET(I)
      GT(I,J) = TG(I)
      SNOAMT(I,J) = SNR(I)
C***********************
C      FOR ICE MODELLING
C***********************
      GHICE(I,J)  = CHICE(I)
      GHICE2(I,J) = CHICE2(I)
      GT1(I,J)    = TG1(I)
      GT2(I,J)    = TG2(I)
      GIceComp(I,J) = CIceComp(I)
      SnoAmt1(I,J)  = SNR1(I)
      SnoAmt2(I,J)  = SNR2(I)
C***********************
  210 TS(I,J)=T4(I) !surface air temperature
C************************************************************
C             TOTAL HEATING AND MOISTURE BUDGET
C************************************************************
      DO 220 I=2,73
      QW(I,J,1) = QW(I,J,1)-FLUXQ(I)-PREC1(I)+CQ1(I)+PCQ1(I)
  220 QW(I,J,2) = QW(I,J,2)-(PREC3(I)-CQ3(I)-PCQ3(I)-FLUXQ(I)-PREC1(I)
     1  -2.0*GRAV*((1.-FLD(I))*E4(I)+FLD(I)*E4C(I))* DTC3/(10.*SP(I)))
C
      DO 230 I=2,73
      H1(I)=(AS1(I)+R2(I)-R0(I))*COE1*COLMR(I)+
     1   CT1(I)+PCT1(I)+CLH*PREC1(I)
      H3(I)=(AS3(I)+(1.0-FLD(I))*(R4(I)+F4(I))+
     1   FLD(I)*(R4C(I)+F4C(I))-R2(I))
     1  *COE1*COLMR(I)+CT3(I)+PCT3(I)+CLH*(PREC3(I)-PREC1(I))
  230 CONTINUE
C
      DO 235 I=2,73
      HEAT(I)=H1(I)+H3(I)
      TEMB(I)=0.5*HEAT(I)*DXYP(J)
  235 CONTINUE
C
      TEMB(1)=TEMB(73)
      TEMB(74)=TEMB(2)
C
      DO 240 I=2,73
  240 H(I,3)=TEMB(I-1)+TEMB(I)+TEMB(I)+TEMB(I+1)
C
      IF (J.NE.1) GO TO 255
C
C        UPDATE FIRST ROW TEMPERATURES
C
      DO 250 I=2,73
      T(I,J,1)=T(I,J,1)+H1(I)
      T(I,J,2)=T(I,J,2)+H3(I)
  250 CONTINUE
C
      GO TO 290
C
  255 IF (J.NE.JM) GO TO 265
C
C        UPDATE LAST ROW TEMPERATURES
C
      DO 260 I=2,73
      T(I,J,1)=T(I,J,1)+H1(I)
  260 T(I,J,2)=T(I,J,2)+H3(I)
C
      GO TO 269
C
  265 DO 268 I=2,73
      T(I,J,1)=T(I,J,1)+0.5*(H1(I)-H3(I))
  268 T(I,J,2)=T(I,J,2)-0.5*(H1(I)-H3(I))
C
      IF (J.LT.3) GO TO 290
C
  269 DO 270 I=2,73
      TEMA(I)=0.0625*(H(I,3)+H(I,2)+H(I,2)+H(I,1))/DXYP(J-1)
      T(I,J-1,1)=T(I,J-1,1)+TEMA(I)
      T(I,J-1,2)=T(I,J-1,2)+TEMA(I)
  270 CONTINUE
C
  290 CONTINUE
C*********************************************************************
C             SURFACE FRICTION AND MOMENTUM REDISTRIBUTION
C*********************************************************************
      IF (J .EQ. 1) GO TO 440
C
      ROS(74,2)=ROS(2,2)
      ROS(74,1)=ROS(2,1)
      DO 320 I=2,73
      ROUV(I)=DXYU2(J)*(ROS(I,2)+ROS(I+1,2))
     1       +DXYU1(J)*(ROS(I,1)+ROS(I+1,1))
      FK(I)=ROUV(I)*CDUV(I)*(10.*GRAV)/(DSIG*PM)*dMAX1(WMAGUV(I),GUSTY)
  320 COLMRV(I)=4.*PM/(P(I,J)+P(I+1,J)+P(I,J-1)+P(I+1,J-1))
C
      BETADT(74,2)=BETADT(2,2)
      BETADT(74,1)=BETADT(2,1)
      DO 360 I=2,73
      BETA(I)=DXYU2(J)*(BETADT(I,2)+BETADT(I+1,2))+
     1        DXYU1(J)*(BETADT(I,1)+BETADT(I+1,1))
      FMCOL(I)=FM*COLMRV(I)*COLMRV(I)*DTC3
      TEMA(I)=dMAX1(FMCOL(I),BETA(I))
      TEMB(I)=FK(I)*COLMRV(I)*DTC3
      DU13(I)=TEMA(I)*(U(I,J,1)-U(I,J,2))
      DV13(I)=TEMA(I)*(V(I,J,1)-V(I,J,2))
  360 CONTINUE
C
      DO 370 I=2,73
c      write(21,10068) i,j,u(i,j,1),u(i,j,2),v(i,j,1),v(i,j,2)
c10068  format (1x,2i2,4f8.3)     
      U(I,J,1)=U(I,J,1)-DU13(I)
      U(I,J,2)=U(I,J,2)+DU13(I)-TEMB(I)*US(I,1)
      V(I,J,1)=V(I,J,1)-DV13(I)
      V(I,J,2)=V(I,J,2)+DV13(I)-TEMB(I)*VS(I,1)
c      write(21,10069) i,j,u(i,j,1),u(i,j,2),v(i,j,1),v(i,j,2),
c     1 DU13(I), DV13(I), US(I,1), VS(I,1), TEMB(I)
c10069  format (1x,2i2,9f8.3)     
  370 CONTINUE
  440 CONTINUE

C*********************************************************************
C             PREPARATION FOR OCEAN MODEL
C*********************************************************************
      do I=2,73                                                    
c       STRESS(I)=ROS(I,2)*CDPI(I)*DMAX1(WMAGPI(I),GUSTY)
c       TAUU(I,J)=1.d+3*STRESS(I)*USPI/24.0
c       TAUV(I,J)=1.d+3*STRESS(I)*VSPI/24.0
c       TXN(I-1,J)=TAUU(I,J)
c       TYN(I-1,J)=TAUV(I,J)
	 if (SEAICE(I)) then
        QSN(I-1,47-J)=(((S4(I)-R4(I))     !*CIce2Alb(I)
	1   -F4(I))*0.484D0 - 2.5D7*E4(I))/TSPD
     2   +QSN(I-1,47-J)
	 else
        QSN(I-1,47-J)=((S4(I)-R4(I)-F4(I))*0.484D0-2.5D7*E4(I))/TSPD
     *   +QSN(I-1,47-J)
	 endif
c       QSN(I-1,17-J)=(S4(I)+DR4(I)-F4(I)-STBO*(TG(I)**4)
c     1  -2.5D7*E4(I))*0.125D0 + QSN(I-1,17-J)
	enddo
 	if (CalcWind) then
	 do i=1,74
	  WIND(I,J) = WINDF(I)
	 enddo
	endif

	RETURN
      END
