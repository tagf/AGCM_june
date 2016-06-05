C     ************************************
C     ************************************
      SUBROUTINE QSDQV(TQS,PQS,QSAT,FSGAM)
C     ************************************
C     ************************************
C
C
C        COMPUTES SATURATION MIXING RATIO AND CHANGE IN SATURATION
C        MIXING RATIO WITH RESPECT TO TEMPERATURE. VECTOR CODED --
C        COMPUTATIONS ARE FOR AN ENTIRE LATITUDE
C
C
C     POLYNOMIAL APPROXIMATIONS FROM:
C        LOWE, P.R. AND J.M. FICKE, 1974: THE COMPUTATION OF SATURATION
C        VAPOR PRESSURE, TECH. PAPER NO. 4-74, ENVIRONMENTAL RESEARCH
C        PREDICTION FACILITY, NAVAL POSTGRADUATE SCHOOL. 11/25/75
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION TQS(74),PQS(74),QSAT(74),FSGAM(74)
      DIMENSION ESTAR(74),ESTARP(74)
      DIMENSION A0(74),A1(74),A2(74),A3(74),A4(74),A5(74)
      DIMENSION A6(74)
      DIMENSION B0(74),B1(74),B2(74),B3(74),B4(74),B5(74)
      DIMENSION B6(74)
      DIMENSION TC(74),PQEST(74)
      DATA AH0/6.107799961E0/, AH1/4.436518521E-1/, AH2/1.428945805E-2/
     *, AH3/2.650648471E-4/, AH4/3.031240396E-6/, AH5/2.034080948E-8/
     *, AH6/6.136820929E-11/
      DATA AL0/4.866786841E0/, AL1/3.152625546E-1/, AL2/8.640188586E-3/
     *, AL3/1.279669658E-4/, AL4/1.077955914E-6/, AL5/4.886796102E-9/
     *, AL6/9.296950850E-12/
      DATA BH0/4.438099984E-1/, BH1/2.857002636E-2/, BH2/7.938054040E-4/
     *, BH3/1.215215065E-5/, BH4/1.036561403E-7/, BH5/3.532421810E-10/
     *, BH6/-7.090244804E-13/
      DATA BL0/4.086240791E-1/, BL1/2.516118369E-2/, BL2/6.576862688E-4/
     *, BL3/9.325531518E-6/, BL4/7.550718726E-8/, BL5/3.303373957E-10/
     *, BL6/6.088242842E-13/
C
      DO 10 I=2,73
      TC(I) = TQS(I) - 273.155E0
   10 CONTINUE
C
      DO 20 I=2,73
      A0(I)=AH0
      A1(I)=AH1
      A2(I)=AH2
      A3(I)=AH3
      A4(I)=AH4
      A5(I)=AH5
      A6(I)=AH6
   20 CONTINUE
C
      DO 30 I=2,73
      B0(I)=BH0
      B1(I)=BH1
      B2(I)=BH2
      B3(I)=BH3
      B4(I)=BH4
      B5(I)=BH5
      B6(I)=BH6
   30 CONTINUE
C
      DO 40 I=2,73
      IF (TQS(I).GE.224.0) GO TO 40
      A0(I)=AL0
      A1(I)=AL1
      A2(I)=AL2
      A3(I)=AL3
      A4(I)=AL4
      A5(I)=AL5
      A6(I)=AL6
   40 CONTINUE
C
      DO 50 I=2,73
      IF (TQS(I).GE.224.0) GO TO 50
      B0(I)=BL0
      B1(I)=BL1
      B2(I)=BL2
      B3(I)=BL3
      B4(I)=BL4
      B5(I)=BL5
      B6(I)=BL6
   50 CONTINUE
C
      DO 60 I=2,73
      ESTAR(I)=(((((A6(I)*TC(I)+A5(I))*TC(I)+A4(I))*TC(I)
     1   +A3(I))*TC(I)+A2(I))*TC(I)+A1(I))*TC(I)+A0(I)
   60 CONTINUE
C
      DO 65 I=2,73
      ESTARP(I)=(((((B6(I)*TC(I)+B5(I))*TC(I)+B4(I))*TC(I)
     1   +B3(I))*TC(I)+B2(I))*TC(I)+B1(I))*TC(I)+B0(I)
   65 CONTINUE
C
      DO 70 I=2,73
      PQEST(I)=(PQS(I)-ESTAR(I))
      QSAT(I)=0.622*ESTAR(I)/PQEST(I)
      FSGAM(I)=0.622*ESTARP(I)*PQS(I)/(PQEST(I)*PQEST(I))
   70 CONTINUE
      RETURN
      END
C     ******************************
C     ******************************
      SUBROUTINE QSATV(TQS,PQS,QSAT)
C     ******************************
C     ******************************
C
C
C        COMPUTES SATURATION MIXING RATIO. VECTOR CODED -- COMPUTATIONS
C        ARE FOR AN ENTIRE LATITUDE
C
C
C     POLYNOMIAL APPROXIMATIONS FROM:
C        LOWE, P.R. AND J.M. FICKE, 1974: THE COMPUTATION OF SATURATION
C        VAPOR PRESSURE, TECH. PAPER NO. 4-74, ENVIRONMENTAL RESEARCH
C        PREDICTION FACILITY, NAVAL POSTGRADUATE SCHOOL. 11/25/75
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION TQS(74),PQS(74),QSAT(74)
      DIMENSION ESTAR(74),ESTARP(74)
      DIMENSION A0(74),A1(74),A2(74),A3(74),A4(74),A5(74)
      DIMENSION A6(74)
      DIMENSION TC(74)
      DATA AH0/6.107799961E0/, AH1/4.436518521E-1/, AH2/1.428945805E-2/
     *, AH3/2.650648471E-4/, AH4/3.031240396E-6/, AH5/2.034080948E-8/
     *, AH6/6.136820929E-11/
      DATA AL0/4.866786841E0/, AL1/3.152625546E-1/, AL2/8.640188586E-3/
     *, AL3/1.279669658E-4/, AL4/1.077955914E-6/, AL5/4.886796102E-9/
     *, AL6/9.296950850E-12/
C
      DO 10 I=2,73
      TC(I) = TQS(I) - 273.155E0
   10 CONTINUE
C
      DO 20 I=2,73
      A0(I)=AH0
      A1(I)=AH1
      A2(I)=AH2
      A3(I)=AH3
      A4(I)=AH4
      A5(I)=AH5
      A6(I)=AH6
   20 CONTINUE
C
      DO 40 I=2,73
      IF (TQS(I).GE.224.0) GO TO 40
      A0(I)=AL0
      A1(I)=AL1
      A2(I)=AL2
      A3(I)=AL3
      A4(I)=AL4
      A5(I)=AL5
      A6(I)=AL6
   40 CONTINUE
C
      DO 60 I=2,73
      ESTAR(I)=(((((A6(I)*TC(I)+A5(I))*TC(I)+A4(I))*TC(I)
     1   +A3(I))*TC(I)+A2(I))*TC(I)+A1(I))*TC(I)+A0(I)
   60 CONTINUE
C
      DO 70 I=2,73
      QSAT(I)=0.622*ESTAR(I)/(PQS(I)-ESTAR(I))
   70 CONTINUE
C
      RETURN
      END
C     *****************************
C     *****************************
      FUNCTION QSDQS(TQS,PQS,FSGAM)
C     *****************************
C     *****************************
C
C
C        COMPUTES SATURATION MIXING RATIO AND CHANGE IN SATURATION
C        MIXING RATIO WITH RESPECT TO TEMPERATURE. SCALAR CODED --
C        COMPUTATIONS ARE FOR A SINGLE POINT
C
C
C     POLYNOMIAL APPROXIMATIONS FROM:
C        LOWE, P.R. AND J.M. FICKE, 1974: THE COMPUTATION OF SATURATION
C        VAPOR PRESSURE, TECH. PAPER NO. 4-74, ENVIRONMENTAL RESEARCH
C        PREDICTION FACILITY, NAVAL POSTGRADUATE SCHOOL. 11/25/75
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DATA AH0/6.107799961E0/, AH1/4.436518521E-1/, AH2/1.428945805E-2/
     *, AH3/2.650648471E-4/, AH4/3.031240396E-6/, AH5/2.034080948E-8/
     *, AH6/6.136820929E-11/
      DATA AL0/4.866786841E0/, AL1/3.152625546E-1/, AL2/8.640188586E-3/
     *, AL3/1.279669658E-4/, AL4/1.077955914E-6/, AL5/4.886796102E-9/
     *, AL6/9.296950850E-12/
      DATA BH0/4.438099984E-1/, BH1/2.857002636E-2/, BH2/7.938054040E-4/
     *, BH3/1.215215065E-5/, BH4/1.036561403E-7/, BH5/3.532421810E-10/
     *, BH6/-7.090244804E-13/
      DATA BL0/4.086240791E-1/, BL1/2.516118369E-2/, BL2/6.576862688E-4/
     *, BL3/9.325531518E-6/, BL4/7.550718726E-8/, BL5/3.303373957E-10/
     *, BL6/6.088242842E-13/
C
      TC = TQS - 273.155E0
C
      IF (TQS.LT.224.0) GO TO 10
      ESTAR =(((((AH6*TC+AH5)*TC+AH4)*TC+AH3)*TC+AH2)*TC+AH1)*TC+AH0
      ESTARP=(((((BH6*TC+BH5)*TC+BH4)*TC+BH3)*TC+BH2)*TC+BH1)*TC+BH0
      GO TO 20
   10 ESTAR =(((((AL6*TC+AL5)*TC+AL4)*TC+AL3)*TC+AL2)*TC+AL1)*TC+AL0
      ESTARP=(((((BL6*TC+BL5)*TC+BL4)*TC+BL3)*TC+BL2)*TC+BL1)*TC+BL0
   20 PQEST=PQS-ESTAR
      QSDQS=0.622*ESTAR/PQEST
      FSGAM=0.622*ESTARP*PQS/(PQEST*PQEST)
      RETURN
      END
C     **********************
C     **********************
      FUNCTION QSAT(TQS,PQS)
C     **********************
C     **********************
C
C
C        COMPUTES SATURATION MIXING RATIO. SCALAR CODED -- COMPUTATIONS
C        ARE FOR A SINGLE POINT
C
C
C     POLYNOMIAL APPROXIMATIONS FROM:
C        LOWE, P.R. AND J.M. FICKE, 1974: THE COMPUTATION OF SATURATION
C        VAPOR PRESSURE, TECH. PAPER NO. 4-74, ENVIRONMENTAL RESEARCH
C        PREDICTION FACILITY, NAVAL POSTGRADUATE SCHOOL. 11/25/75
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DATA AH0/6.107799961E0/, AH1/4.436518521E-1/, AH2/1.428945805E-2/
     *, AH3/2.650648471E-4/, AH4/3.031240396E-6/, AH5/2.034080948E-8/
     *, AH6/6.136820929E-11/
      DATA AL0/4.866786841E0/, AL1/3.152625546E-1/, AL2/8.640188586E-3/
     *, AL3/1.279669658E-4/, AL4/1.077955914E-6/, AL5/4.886796102E-9/
     *, AL6/9.296950850E-12/
C
      TC = TQS - 273.155E0
C
      IF (TQS.GE.224.0)
     1 ESTAR =(((((AH6*TC+AH5)*TC+AH4)*TC+AH3)*TC+AH2)*TC+AH1)*TC+AH0
      IF (TQS.LT.224.0)
     1 ESTAR =(((((AL6*TC+AL5)*TC+AL4)*TC+AL3)*TC+AL2)*TC+AL1)*TC+AL0
      QSAT=0.622*ESTAR/(PQS-ESTAR)
      RETURN
      END
C     *************************************
C     *************************************
      SUBROUTINE ISFCUV(J,ISFJ,ISFJM,ISFUV)
C     *************************************
C     *************************************
C
C
C        IDENTIFIES SURFACE TYPE IN ORDER TO COMPUTE SURFACE DRAG
C        COEFFICIENT
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION ISFJM(74),ISFJ(74),ISFUV(74)
      LOGICAL WATER(74),HALFJ(74),ALLJ(74),HALFJM(74),ALLJM(74)
C
      DO 10 I=1,74
   10 ISFUV(I)=0
C
      IF (J.LT.2) RETURN
      ISFJM(74)=ISFJM(2)
      ISFJ (74)=ISFJ (2)
C
      DO 20 I=2,73
      HALFJ(I)=ISFJ(I).EQ.7.OR.ISFJ(I+1).EQ.7
      ALLJ(I)=ISFJ(I).EQ.7.AND.ISFJ(I+1).EQ.7
      HALFJM(I)=ISFJM(I).EQ.7.OR.ISFJM(I+1).EQ.7
      ALLJM(I)=ISFJM(I).EQ.7.AND.ISFJM(I+1).EQ.7
   20 CONTINUE
C
      DO 30 I=2,73
   30 WATER(I)=(ALLJ(I).AND.HALFJM(I)).OR.(ALLJM(I).AND.HALFJ(I))
      IF (J.EQ.24) GO TO 100
      IF (J.GT.24) GO TO 60
C
C        SOUTHERN HEMISPHERE
C
      DO 40 I=2,73
   40 WATER(I)=WATER(I).OR.ALLJ(I)
      GO TO 100
C
C        NORTHERN HEMISPHERE
C
   60 DO 70 I=2,73
   70 WATER(I)=WATER(I).OR.ALLJM(I)
C
  100 DO 150 I=2,73
      IF (WATER(I)) ISFUV(I)=7
  150 CONTINUE
      RETURN
      END
C     ****************
C     ****************
      SUBROUTINE LONGW
C     ****************
C     ****************
C
C
C             CALCULATES LONGWAVE RADIATION BUDGET
C
C
C***********************************************************************
C                                                                      *
C     NEW RADIATION SCHEME                                             *
C     1. MOISTURE CARRIED AT LEVELS 1 AND 3                            *
C     2. NEW COMPUTATION OF VERTICAL PROFILE OF EFFECTIVE WATER        *
C        VAPOR AMOUNT                                                  *
C     3. NEW METHOD FOR DETERMINING WHICH CLOUD TYPE EXISTS            *
C     4. NEW NOMENCLATURE FOR SOLAR RADIATION CALCULATION              *
C     5. NEW SOLAR CONSTANT, S0 = 1.94 LY/MIN                          *
C     6. MAGNIFACTION FACTOR AFTER RODGERS(1967)                       *
C     7. SURFACE ALBEDO A FUNCTION OF SURFACE TYPE                     *
C     8. NEW WATER VAPOR ABSORPTION FUNCTION FOR SOLAR RADIATION       *
C     9. OZONE ABSORPTION OF THE SCATTERED RADIATION                   *
C                                                                      *
C***********************************************************************
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      LOGICAL SNOW,SEAICE
C
      include 'radvar.fi'
C
      COMMON /FACTR/  EST0(74),E10(74),E20(74),E02(74),E32(74),
     1                E12(74),E42(74),E24(74),E34(74)
C
      COMMON /MOIST/ EFVT(74),EFVST(74),EFV3(74),EFV2(74),
     1               EFV1(74),EFV0(74),AK(74)
C
	include 'cover.fi'
C
	include 'trans.fi'
C
      COMMON /WORK/ TEM(74),TDR2(74),TDR4(74),TUR0(74),TUR2(74)
     X  ,BL0(74),BL1(74),BL2(74),BL3(74),BL4(74),BLG(74)
     X  ,BLST(74),DBCST(74),DBST0(74),DB12(74)
     X  ,DB01(74),DB02(74),DB23(74),DB24(74),DB34(74),DB4G(74)
     X  ,TB34(74),TB24(74),TB10(74),TB20(74),TB32(74)
     X  ,TB42(74),TBST0(74),TB02(74),TB12(74)
     X  ,UR0(74),UR2(74)
     X  ,ALSF(74),SECZ(74),SA(74),SST(74),SS(74),ALA0(74)
     X  ,SAO0(74),SAO2(74),SAO4(74),SAOG(74)
     X  ,ASO1(74),ASO3(74),ASC1(74),ASC3(74)
     X  ,SAC0(74),SAC2(74),SAC4(74),SACG(74)
     X  ,SSO4(74),SSOG(74),SSC4(74),SSCG(74)
     X  ,TWR0(74),TWR1(74),TWR2(74),TWR3(74),TWR4(74)
     X  ,ALAC(74),EFVC(74),SACRC(74),COT(74),COTL(74)
     X  ,RABAR(74),TDR4L(74),TUR0L(74),TUR2L(74)
     X  ,TDR2L(74),EFVCL(74),SACRCL(74),RSBAR(74),ALACL(74)
C
C************ END OF COMMON ******************
C
C        CALCULATE BLACKBODY FLUXES
C
      DO 160 I=2,73
      BLST(I)=STBO*TST*TST*TST*TST
      BL0(I)=STBO*TTROP(I)*TTROP(I)*TTROP(I)*TTROP(I)
      BL1(I)=STBO*T1(I)*T1(I)*T1(I)*T1(I)
      BL2(I)=STBO*T2(I)*T2(I)*T2(I)*T2(I)
      BL3(I)=STBO*T3(I)*T3(I)*T3(I)*T3(I)
      BL4(I)=STBO*T4(I)*T4(I)*T4(I)*T4(I)
      BLG(I)=STBO*TG(I)*TG(I)*TG(I)*TG(I)
      DBCST(I)=BLST(I)-BLC
      DBST0(I)=BL0(I)-BLST(I)
      DB02(I)=BL2(I)-BL0(I)
      DB01(I)=BL1(I)-BL0(I)
      DB12(I)=BL2(I)-BL1(I)
      DB23(I)=BL3(I)-BL2(I)
      DB24(I)=BL4(I)-BL2(I)
      DB34(I)=BL4(I)-BL3(I)
      DB4G(I)=BLG(I)-BL4(I)
  160 CONTINUE
C
C        CALCULATE UPWARD AND DOWNWARD COMPONENTS, CLEAR PART
C
      DO 180 I=2,73
      TB34(I)=DB34(I)*(1.0+E34(I)*TR34(I))/(1.0+E34(I))
      TB24(I)=DB24(I)*(1.0+E24(I)*TR24(I))/(1.0+E24(I))
      TB10(I)=DB01(I)*(1.0+E10(I)*TR01(I))/(1.0+E10(I))
      TB20(I)=DB02(I)*(1.0+E20(I)*TR02(I))/(1.0+E20(I))
      TB32(I)=DB23(I)*(1.0+E32(I)*TR23(I))/(1.0+E32(I))
      TB42(I)=DB24(I)*(1.0+E42(I)*TR24(I))/(1.0+E42(I))
      TBST0(I)=DBST0(I)*(1.0+EST0(I)*TRST0(I))/(1.0+EST0(I))
      TB02(I)=DB02(I)*(1.0+E02(I)*TR02(I))/(1.0+E02(I))
      TB12(I)=DB12(I)*(1.0+E12(I)*TR12(I))/(1.0+E12(I))
  180 CONTINUE
      DO 200 I=2,73
      DR0(I)=BL0(I)-BLC*T2RT0(I)-DBCST(I)*TRT0(I)-TBST0(I)
      DR2(I)=BL2(I)-BLC*T2RT2(I)-DBCST(I)*TRT2(I)-TB02(I)
     1    -0.5*DBST0(I)*(TR02(I)+TRST2(I))
      UR2(I)=BL2(I)+DB4G(I)*TR24(I)+TB42(I)
      DR4(I)=BL4(I)-BLC*T2RT4(I)-DBCST(I)*TRT4(I)-TB24(I)
     1  -0.5*(DBST0(I)*(TR04(I)+TRST4(I))+DB02(I)*(TR04(I)+TR24(I)))
      UR0(I)=BL0(I)+DB4G(I)*TR04(I)+TB20(I)+
     1    0.5*DB24(I)*(TR04(I)+TR02(I))
  200 CONTINUE
C
C        CALCULATE UPWARD AND DOWNWARD COMPONENTS, CLOUDY PART
C
      DO 220 I=2,73
      TDR2(I)=0.0
      TDR4(I)=0.0
      TUR0(I)=0.0
      TUR2(I)=0.0
      TDR4L(I)=0.0
      TUR0L(I)=0.0
      TUR2L(I)=0.0
  220 CONTINUE
      DO 310 I=2,73
      IF (JCLOUD(I).NE.1) GO TO 310
      TDR2(I)=TB12(I)
      TDR4(I)=TB24(I)+0.5*DB12(I)*(TR14(I)+TR24(I))
      TUR2(I)=UR2(I)-BL2(I)
  310 CONTINUE
      DO 320 I=2,73
      IF(JCLOUD(I).NE.2) GO TO 320
      TDR2(I)=BL2(I)-DR2(I)
      TDR4(I)=TB34(I)
      TUR0(I)=TB20(I)
  320 CONTINUE
      DO 330 I=2,73
      IF(JCLOUD(I).NE.3) GO TO 330
      TDR2(I)=BL2(I)-DR2(I)
      TDR4(I)=TB34(I)
      TUR0(I)=TB20(I)+0.5*DB23(I)*(TR03(I)+TR02(I))
      TUR2(I)=TB32(I)
  330 CONTINUE
      DO 340 I=2,73
      IF(JCLOUD(I).NE.4) GO TO 340
      TDR2(I)=0.0
      TDR4(I)=TB24(I)
      TUR0(I)=TB10(I)
      TUR2(I)=UR2(I)-BL2(I)
  340 CONTINUE
      DO 400 I=2,73
      IF(JCLOUD(I).NE.4) GO TO 400
      IF(JCLL(I).EQ.0) GO TO 400
      TDR4L(I)=TB34(I)
      TUR0L(I)=TB20(I)
      IF(JCLL(I).NE.3) GO TO 400
      TUR0L(I)=TB20(I)+0.5*DB23(I)*(TR03(I)+TR02(I))
      TUR2L(I)=TB32(I)
  400 CONTINUE
C
C        ADD CLOUDY PART TO CLEAR PART
C
      DO 500 I=2,73
      DR2(I)=DR2(I)+CLR(I)*(BL2(I)-TDR2(I)-DR2(I))
      UR2(I)=UR2(I)+CLR(I)*(BL2(I)+TUR2(I)-UR2(I))
     1   +CLRL(I)*(BL2(I)+TUR2L(I)-UR2(I))
  500 CONTINUE
      DO 600 I=2,73
      DR4(I)=DR4(I)+CLR(I)*(BL4(I)-TDR4(I)-DR4(I))
  600 DR4(I)=DR4(I)+CLRL(I)*(BL4(I)-TDR4L(I)-DR4(I))
      DO 700 I=2,73
      UR0(I)=UR0(I)+CLRL(I)*(BL0(I)+TUR0L(I)-UR0(I))
  700 UR0(I)=UR0(I)+CLR(I)*(BL0(I)+TUR0(I)-UR0(I))
C
C        CALCULATE NET UPWARD FLUX
C
      DO 800 I=2,73
      R0(I)=UR0(I)-DR0(I)
      R2(I)=UR2(I)-DR2(I)
      R4(I)=BLG(I)-DR4(I)
      R4C(I)=STBO*TGC(I)*TGC(I)*TGC(I)*TGC(I)-DR4(I)
  800 CONTINUE
      RETURN
      END
