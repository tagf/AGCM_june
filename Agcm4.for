C     ****************
C     ****************
      SUBROUTINE SOLAR
C     ****************
C     ****************
C
C
C             CALCULATES SHORT WAVE RADIATION BUDGET
C
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
      COMMON /WORK/
     1  SAO0(74),SAO2(74),SAO4(74),SAC0(74),SAC2(74),SAC4(74),
     2  AM  (74),SA  (74),SECZ(74),ALSF(74),CAO0(74),CAO2(74),
     3  CAO4(74),CCO0(74),CCO2(74),CCO4(74),SACTL(74),SACTU(74),
     4  CCTL(74),CCTU(74),EFVCU(74),EFVCL(74),YCU(74),YCL(74),
     5  ALA0(74),ALAC(74),TEMA(74),TEMB(74),TEMC(74),TC4L(74),
     6  TO0 (74),TO1 (74),TO2 (74),TO3 (74),TO4 (74),TC0 (74),
     7  TC2 (74),TC4 (74),TWT0(74),TWT1(74),TWT2(74),TWT3(74),
     8  TWT4(74),TWC (74),TOC (74),TWCL(74),TOCL(74),TWC2(74),
     9  TWC4(74),TWC4L(74),RABAR(74),RSBAR(74),ALACL(74),SS(74),
     X  SST(74),SAC4L(74),SAG(74),SSG(74),EV01(74),EV12(74),EV23(74),
     X  RACTL(74)
C
C************ END OF COMMON ******************
C
      DIMENSION RT (148)
      DIMENSION RTR(74,38),RTM(74,7),RTC(74,14)
c      integer (4) RTRint(74,39:43),RTint(148)
      integer RTRint(74,39:43),RTint(148)
      DIMENSION RTRres(74,44:47)
      EQUIVALENCE (RH1(1),RTR(1,1)),(EFVT(1),RTM(1,1))
      EQUIVALENCE (JCLOUD(1),RTRint(1,39)),(P1L(1),RTRres(1,44))
      EQUIVALENCE (ALBDOA(1),RTC(1,1))
C
         SQRT(XXX)=DSQRT(XXX)
         ALOG10(XXX)=DLOG10(XXX)
      FRACT(X,Y)=(1.0-X)/(1.0-X*Y)
C
      DO 1 I=1,74
      S4(I)=0.0
      AS1(I)=0.0
      AS3(I)=0.0
      RETOT(I)=0.0
    1 CONTINUE
C
C        FIND LONGITUDES THAT RECEIVE SOLAR RADIATION AND SPIN THAT
C        PORTION OF LATITUDE ARRAY TO ITS BEGINNING
C
      IL=1
      COSZ(1)=COSZ(73)
      COSZ(74)=COSZ(2)
      IS=2
C
      DO 2 I=2,73
      IF (COSZ(I-1).LE.0.01.AND.COSZ(I).GT.0.01) IS=I
      IF (COSZ(I).GT.0.01) IL=IL+1
    2 CONTINUE
C
      IF (IS.EQ.2) GO TO 15
C
      DO 5 J=1,38  !47
      DO 5003 I=2,73
 5003 RT(I)=RTR(I,J)
      DO 3 I=2,73
    3 RT(I+72)=RTR(I,J)
      DO 4 I=2,73
    4 RTR(I,J)=RT(I+IS-2)
    5 CONTINUE

      DO 85 J=39,43
      DO 5203 I=2,73
 5203 RTint(I)=RTRint(I,J)
      DO 83 I=2,73
   83 RTint(I+72)=RTRint(I,J)
      DO 84 I=2,73
   84 RTRint(I,J)=RTint(I+IS-2)
   85 CONTINUE

      DO 95 J=44,47
      DO 5103 I=2,73
 5103 RT(I)=RTRres(I,J)
      DO 93 I=2,73
   93 RT(I+72)=RTRres(I,J)
      DO 94 I=2,73
   94 RTRres(I,J)=RT(I+IS-2)
   95 CONTINUE
C
      DO 8 J=1,7
      DO 5006 I=2,73
 5006 RT(I)=RTM(I,J)
      DO 6 I=2,73
    6 RT(I+72)=RTM(I,J)
      DO 7 I=2,73
    7 RTM(I,J)=RT(I+IS-2)
    8 CONTINUE
C
      DO 11 J=1,14
      DO 5009 I=2,73
 5009 RT(I)=RTC(I,J)
      DO 9 I=2,73
    9 RT(I+72)=RTC(I,J)
      DO 10 I=2,73
   10 RTC(I,J)=RT(I+IS-2)
   11 CONTINUE
C
   15 IF (IL.LT.2) RETURN
C
C        SURFACE ALBEDO
C
      DO 17 I=2,IL
      TEMA(I)=ALS(I)
      TEMB(I)=0.0
   17 TEMC(I)=0.0
C
      DO 20 I=2,IL   !water albedo
      IF (ISRFCE(I).NE.7) GO TO 20
      TEMA(I)=0.99615721
      TEMB(I)=5.4077056
      TEMC(I)=9.2530638
      IF (COSZ(I).LE.0.258819) GO TO 20
      TEMA(I)=0.4655505
      TEMB(I)=1.2820763
      TEMC(I)=0.940257
      IF (COSZ(I).LT.0.707107) GO TO 20
      TEMA(I)=0.03717367
      TEMB(I)=0.01725384
      TEMC(I)=0.0
   20 CONTINUE
C
      DO 30 I=2,IL
   30 ALSF(I)=(TEMA(I)-COSZ(I)*(TEMB(I)-COSZ(I)*TEMC(I)))
C
      DO 40 I=2,IL
   40 SECZ(I)=1.0/COSZ(I)
C
C        MAGNIFICATION FACTOR AFTER RODGERS (1967)
C
      DO 50 I=2,IL
   50 AM(I)=35.*SECZ(I)/SQRT(1224.+SECZ(I)*SECZ(I))
C
      DO 60 I=2,IL
      SA(I)=0.349*SCOSZ(I)
      TEMA(I)=TOZONE*AM(I)
   60 SST(I)=0.651*SCOSZ(I)
C
      CALL TROZON (TEMA,TEMB,IL)
C
      DO 70 I=2,IL
      SS(I)=SST(I)*TEMB(I)
   70 ALA0(I)=0.085-0.247*ALOG10(COSZ(I)*1000.0/P4(I))
                    !0.045-0.247*ALOG10(COSZ(I)*1000.0/P4(I))
C
      DO 80 I=2,IL
   80 ALA0(I)=dMIN1(dMAX1(ALA0(I),0.0d0),1.0d0)
C
      DO 100 I=2,IL
      TO0(I)=(EFVT(I)-EFV0(I))*AM(I)
      TO1(I)=(EFVT(I)-EFV1(I))*AM(I)
      TO2(I)=(EFVT(I)-EFV2(I))*AM(I)
      TO3(I)=(EFVT(I)-EFV3(I))*AM(I)
      TO4(I)=EFVT(I)*AM(I)
  100 CONTINUE
C
      CALL TRWATR (TO0,TWT0,IL)
      CALL TRWATR (TO1,TWT1,IL)
      CALL TRWATR (TO2,TWT2,IL)
      CALL TRWATR (TO3,TWT3,IL)
      CALL TRWATR (TO4,TWT4,IL)
C
C        VAPOR CONTENT AND TRANSMISSION FUNCTIONS FOR CLOUDS
C
      DO 105 I=2,IL
      YCU(I)=1.0
      YCL(I)=1.0
      TOC(I)=0.0
      TWC(I)=0.0
      TOCL(I)=0.0
      TWCL(I)=0.0
  105 CONTINUE
C
C        CLOUD 1
C
      DO 110 I=2,IL
      IF (JCLOUD(I).NE.1) GO TO 110
      TOC(I)=TO0(I)
      TWC(I)=TWT0(I)
  110 CONTINUE
C
C        CLOUD 2
C
      DO 120 I=2,IL
      IF (JCLOUD(I).NE.2) GO TO 120
      TOC(I)=TO2(I)
      TWC(I)=TWT2(I)
      IF (PREC3(I).GT.0.0) GO TO 120
      YCU(I)=10.*(RH3(I)-0.9)
      YCU(I)=dMIN1(1.0d0,YCU(I))
      YCU(I)=dMAX1(0.0d0,YCU(I))
  120 CONTINUE
C
C        CLOUD 3
C
      DO 130 I=2,IL
      IF (JCLOUD(I).NE.3) GO TO 130
      TOC(I)=TO3(I)
      TWC(I)=TWT3(I)
  130 CONTINUE
C
C        CLOUD 4
C
      DO 140 I=2,IL
      IF (JCLOUD(I).NE.4) GO TO 140
      TOC(I)=TO1(I)
      TWC(I)=TWT1(I)
      IF (PREC1(I).GT.0.0) GO TO 140
      YCU(I)=10.*(RH1(I)-0.9)
      YCU(I)=dMIN1(1.0d0,YCU(I))
      YCU(I)=dMAX1(0.0d0,YCU(I))
  140 CONTINUE
C
C        LOWER LEVEL CLOUDS
C
      DO 180 I=2,IL
      IF (JCLL(I).NE.2) GO TO 180
      IF (PREC3(I).GT.0.0) GO TO 180
      YCL(I)=10.*(RH3(I)-0.9)
      YCL(I)=dMIN1(1.0d0,YCL(I))
      YCL(I)=dMAX1(0.0d0,YCL(I))
  180 CONTINUE
C
      DO 190 I=2,IL
      EV01(I)=EFV0(I)-EFV1(I)
      EV12(I)=EFV1(I)-EFV2(I)
      EV23(I)=EFV2(I)-EFV3(I)
  190 CONTINUE
C
C        EQUIVALENT WATER VAPOR FROM TOP OF
C        ATMOSPHERE TO CLOUD BASE OF UPPER LAYER
C
      DO 200 I=2,IL
  200 TEMA(I)=1.0-RWC(I)*TWC(I)
C
      CALL AMTWTR (TEMA,TEMB,IL)
C
      DO 210 I=2,IL
  210 EFVCU(I)=YCU(I)*(TEMB(I)-TOC(I))/1.66
C
      DO 220 I=2,IL
      IF (JCLOUD(I).EQ.1) EFVCU(I)=dMAX1(EFVCU(I),EV01(I))
      IF (JCLOUD(I).EQ.2) EFVCU(I)=dMAX1(EFVCU(I),EV23(I))
      IF (JCLL(I).NE.0) EFVCU(I)=dMAX1(EFVCU(I),EV12(I))
      IF (TWC(I).LE.0.) EFVCU(I)=0.0
  220 CONTINUE
C
C        VAPOR CONTENT AND TRANSMISSION FUNCTIONS
C        FOR LOWER CLOUDS
C
      DO 230 I=2,IL
      IF (JCLL(I).EQ.0) GO TO 230
      TOCL(I)=TOC(I)+1.66*EFVCU(I)
      IF (JCLL(I).EQ.3) TOCL(I)=TOCL(I)+1.66*EV23(I)
  230 CONTINUE
C
      CALL TRWATR (TOCL,TWCL,IL)
C
C        EQUIVALENT WATER VAPOR FROM TOP OF
C        ATMOSPHERE TO CLOUD BASE OF LOWER LAYER
C
      DO 300 I=2,IL
  300 TEMA(I)=1.0-RWCL(I)*TWCL(I)
C
      CALL AMTWTR (TEMA,TEMB,IL)
C
      DO 310 I=2,IL
  310 EFVCL(I)=YCL(I)*(TEMB(I)-TOCL(I))/1.66
C
      DO 320 I=2,IL
      IF (JCLL(I).EQ.2) EFVCL(I)=dMAX1(EFVCL(I),EV23(I))
      IF (TWCL(I).LE.0.0) EFVCL(I)=0.0
  320 CONTINUE
C
C        RADIATION TERMS
C        CLEAR COMPONENTS
C        ABSORBED PART
C
      DO 400 I=2,IL
      SAO0(I)=TWT0(I)
      SAO2(I)=TWT2(I)
      SAO4(I)=TWT4(I)
  400 CONTINUE
C
C        OVERCAST COMPONENTS
C        ABSORBED PART
C
      DO 450 I=2,IL
      TEMA(I)=0.0
      IF (JCLOUD(I).EQ.1) TEMA(I)=EV12(I)+EFVCU(I)
      IF (JCLOUD(I).EQ.4) TEMA(I)=EFVCU(I)
  450 CONTINUE
C
      DO 460 I=2,IL
  460 TC2(I)=TOC(I)+1.66*TEMA(I)
C
      DO 470 I=2,IL
      TEMA(I)=0.0
      IF (JCLOUD(I).EQ.1) TEMA(I)=EFV1(I)
      IF (JCLOUD(I).EQ.2) TEMA(I)=EFV3(I)
      IF (JCLOUD(I).EQ.3) TEMA(I)=EFV3(I)
      IF (JCLOUD(I).EQ.4) TEMA(I)=EFV2(I)
      IF (JCLOUD(I).EQ.3) TC2(I)=TO2(I)
  470 CONTINUE
C
      DO 475 I=2,IL
      TC4L(I)=TOCL(I)+1.66*(EFVCL(I)+EFV3(I))
  475 TC4(I)=TOC(I)+1.66*(EFVCU(I)+TEMA(I))
C
      CALL TRWATR (TC2,TWC2,IL)
      CALL TRWATR (TC4,TWC4,IL)
      CALL TRWATR (TC4L,TWC4L,IL)
C
      DO 480 I=2,IL
      CCO2(I)=1.0-ALBDOA(I)
      CCO4(I)=CCO2(I)
      IF (JCLOUD(I).EQ.2.OR.JCLOUD(I).EQ.3) CCO2(I)=1.0
  480 CONTINUE
C
      DO 490 I=2,IL
      SAC0(I)=SAO0(I)
      SAC2(I)=CCO2(I)*TWC2(I)
      SAC4(I)=CCO4(I)*TWC4(I)
  490 SACTU(I)=ALBDOA(I)*TWC(I)
C
      DO 520 I=2,IL
      TEMA(I)=ALBDAL(I)*ALBDOA(I)
      TEMB(I)=ALBDAL(I)+ALBDOA(I)-TEMA(I)
      TEMC(I)=1.0/(1.0-TEMA(I))
  520 CONTINUE
C
      DO 530 I=2,IL
      SACTL(I)=(ALBDAL(I)-TEMA(I))*TWCL(I)
      SAC4L(I)=(1.0-TEMB(I))*TWC4L(I)
  530 CONTINUE
C
      DO 540 I=2,IL
      RACTL(I)=(ALBDOA(I)-TEMA(I))*TEMC(I)
      RABAR(I)=(TEMB(I)-TEMA(I))*TEMC(I)
C
C        RSBAR=(ALBDOS+ALBDSL-2*ALBDOS*ALBDSL)/(1-ALBDOS*ALBDSL)
C        BUT ALBDOS=ALBDOA, AND ALBDSL=ALBDAL SO,
C
      RSBAR(I)=RABAR(I)
  540 CONTINUE
C
      DO 550 I=2,IL
      ALAC(I)=ALBDOS(I)+ALA0(I)-ALBDOS(I)*ALA0(I)
      ALACL(I)=RSBAR(I)+ALA0(I)-RSBAR(I)*ALA0(I)
  550 CONTINUE
C
      DO 600 I=2,IL
      AS1(I)=SA(I)*(PCL(I)*(SAO0(I)-SAO2(I))+CL(I)*(SAC0(I)-SAC2(I))-
     1 (CL1(I)+CL4(I))*SACTU(I))
  600 CONTINUE
      DO 610 I=2,IL
      AS3(I)=SA(I)*(  PCL(I)*(SAO2(I)-SAO4(I))
     1   +CL(I)*PCLL(I)*(SAC2(I)-SAC4(I))
     1  +CLL(I)*(SAC2(I)-SAC4L(I)-SACTL(I))
     2   -(CL3(I)+CL2(I))*(1.0-CL4(I))*SACTU(I))
  610 CONTINUE
C
C        ABSORBED PART OF S4
C
      DO 750 I=2,IL
      SAG(I)=SA(I)*(PCL(I)*(1.0-ALSF(I))*SAO4(I)+
     1  CL(I)*PCLL(I)*SAC4(I)*FRACT(ALS(I),ALBDOA(I))
     2 +CLL(I)*(SAC4L(I)+RACTL(I)*SACTL(I))*FRACT(ALS(I),RABAR(I))  )
  750 CONTINUE
C
C        SCATTERED PART OF S4
C
      DO 850 I=2,IL
      SSG(I)=SS(I)*(1.0-ALS(I))*(PCL(I)*FRACT(ALA0(I),ALS(I))
     1  +PCLL(I)*CL(I)*FRACT(ALAC(I),ALS(I))
     2  +CLL(I)*FRACT(ALACL(I),ALS(I)))
  850 CONTINUE
C
      DO 900 I=2,IL
  900 S4(I)=SAG(I)+SSG(I)
C
      DO 1000 I=2,IL
 1000 RETOT(I)=SCOSZ(I)-AS1(I)-AS3(I)-SA(I)-S4(I)-(SST(I)-SS(I))
     1          +SAO0(I)*SA(I)
C
C        SPIN LATITUDE ARRAY BACK TO ITS ORIGINAL POSITION
C
      IF (IS.EQ.2) RETURN
      IS=76-IS
C
      DO 2000 J=1,38  !47
      DO 1010 I=2,73
 1010 RT(I)=RTR(I,J)
      DO 1100 I=2,73
 1100 RT(I+72)=RTR(I,J)
      DO 1200 I=2,73
 1200 RTR(I,J)=RT(I+IS-2)
 2000 CONTINUE

      DO 2001 J=39,43  
      DO 1011 I=2,73
 1011 RTint(I)=RTRint(I,J)
      DO 1101 I=2,73
 1101 RTint(I+72)=RTRint(I,J)
      DO 1201 I=2,73
 1201 RTRint(I,J)=RTint(I+IS-2)
 2001 CONTINUE
      DO 2002 J=44,47  
      DO 1012 I=2,73
 1012 RT(I)=RTRres(I,J)
      DO 1102 I=2,73
 1102 RT(I+72)=RTRres(I,J)
      DO 1202 I=2,73
 1202 RTRres(I,J)=RT(I+IS-2)
 2002 CONTINUE
C
      DO 3000 J=1,7
      DO 2010 I=2,73
 2010 RT(I)=RTM(I,J)
      DO 2100 I=2,73
 2100 RT(I+72)=RTM(I,J)
      DO 2200 I=2,73
 2200 RTM(I,J)=RT(I+IS-2)
 3000 CONTINUE
C
      DO 4000 J=1,14
      DO 3010 I=2,73
 3010 RT(I)=RTC(I,J)
      DO 3100 I=2,73
 3100 RT(I+72)=RTC(I,J)
      DO 3200 I=2,73
 3200 RTC(I,J)=RT(I+IS-2)
 4000 CONTINUE
      RETURN
      END
C     *****************
C     *****************
      SUBROUTINE CLOUDS
C     *****************
C     *****************
C
C
C             DETERMINES CLOUD TYPES AND PROPERTIES FOR RADIATION
C                 CALCULATIONS
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      LOGICAL SNOW,SEAICE
C
      include 'radvar.fi'
C
	include 'cover.fi'
C
      LOGICAL CIRUS,CIRUSL
C
      COMMON /WORK/ CIRUS(74),CIRUSL(74),GREY(74),GREYL(74)
     1             ,CLC(74)
C
C************ END OF COMMON ******************
C
C
C     SCHEMATIC REPRESENTATION OF CLOUD .TYPES.
C
C LEVEL                                              P(MB) T(K)
C     ==============================================
C   T ]                                            ]    0     0
C   C ]                                            ] UNDEF  220
C     ]                                            ]
C  ST ]                                            ]  100   218
C     ]                                            ]
C   0 ]    *****                                   ]  200    T0
C     ]   *     *                                  ]
C     ]   * CL1 *                                  ]
C     ]   *     *                                  ]
C   1 ]    *****                                   ]   P1    T1(I)
C     ]    CLOUD1                                  ]
C     ]                                            ]
C     ]                                            ]
C   2 ]                   *****                    ]   P2    T2
C     ]                  *     *                   ]
C     ]                  * CL2 *          CLL      ]
C     ]                  *     *        *******    ]
C   3 ]                   *****         *******    ]   P3    T3
C     ]                   CLOUD2         CLOUD3    ]
C   CLOUD                                          ]
C   .TYPE.   1              2              3       ]
C   JCLOUD(I)                                      ]
C   4 ]                                            ]   P4(I) T4
C   G ]                                            ]   PS    TG
C     ==============================================
C     NOTE: LEVELS T AND 0 ARE DENOTED BY 0 AND T OR TROP RESPECTIVELY
C     OUTSIDE SUB. RADCAL
C
C
C     CLOUD ALBEDOS, ABSORPTIVITY AND EQUIVALENT WATER VAPOR AMOUNT
C     CLOUD PROPERTIES FROM RODGERS,C.D., 1967: THE RADIATIVE HEAT
C     BUDGET OF THE TROPOSPHERE AND LOWER STRATOSPHERE, PLANETARY
C     CIRCULATION PROJECT, DEPT. MET., MIT, REP. NO. 2, 99PP. TABLE
C     8, P. 31. UCLA CN METEOR. REP. M3-12813
C
C
C LEVEL                                              P(MB) T(K)
C     ==============================================
C   T ]                                            ]    0     0
C   C ]                                            ] UNDEF  220
C     ]                                            ]
C  ST ]                                            ]  100   218
C     ]                                            ]
C   0 ]                                            ]  200    T0
C     ]                                            ]
C     ]                                            ]
C   1 ]    *****          *****          *****     ]   P1    T1(I)
C     ]   *     *        *     *        *     *    ]
C     ]   * CL4 *        *     *        *     *    ]
C     ]   *     *        *     *        *     *    ]
C   2 ]    *****         *******         *****     ]   P2    T2
C     ]    CLOUD4        *     *                   ]
C     ]                  *     *                   ]
C     ]                  *     *        *******    ]
C   3 ]                   *****         *******    ]   P3    T3
C   CLOUD                                          ]
C   .TYPE.   4              5              6       ]
C   JCLOUD(I)                                      ]
C   4 ]                                            ]   P4(I) T4
C   G ]                                            ]   PS    TG
C     ==============================================
C
      DO 1 I=1,74
      CIRUS(I)=.FALSE.
      CIRUSL(I)=.FALSE.
      JCLOUD(I)=0
      JCLL(I)=0
      GREY(I)=1.
      GREYL(I)=1.
      ALBDOA(I)=0.0
      ABSVTY(I)=0.0
      ALBDOS(I)=0.0
      ALBDAL(I)=0.0
      ABSVYL(I)=0.0
      ALBDSL(I)=0.0
      RWC(I)=1.0
      RWCL(I)=1.0
      CL(I) = 0.0
      CLR(I) = 0.0
      CLL(I) = 0.0
      CLRL(I) = 0.0
      CL1(I) = 0.0
      CL2(I) = 0.0
      CL3(I) = 0.0
      CL4(I) = 0.0
      CLC(I)=1.0
    1 CONTINUE
C
C        DETERMINATION OF CLOUD .TYPE.
C
      NCL1=0
      NCL2=0
      NCL3=0
      NCL4=0
      DO 10 I=2,73
      IF (CT1(I).LE.0.0.AND.PCT1(I).LE.0.0) GO TO 10
      CL(I)=CLC(I)
      CL1(I)=CL(I)
      JCLOUD(I)=1
      NCL1=NCL1+1
   10 CONTINUE
C
      DO 20 I=2,73
      IF (JCLOUD(I).NE.0) GO TO 20
      IF (RH3(I).LE.0.9) GO TO 20
      CL(I)=1.0
      CL2(I)=CL(I)
      JCLOUD(I)=2
      NCL2=NCL2+1
      CIRUS(I)=(PREC3(I).LE.0.0).AND.(T2(I)+T3(I).LE.466.0)
   20 CONTINUE
C
      DO 30 I=2,73
      IF (JCLOUD(I).NE.0) GO TO 30
      IF (EX(I).LE.0.0.OR.PCT1(I).NE.0.0) GO TO 30
      CL(I)=CLC(I)
      CL3(I)=CL(I)
      JCLOUD(I)=3
      NCL3=NCL3+1
      CIRUS(I)=T3(I).LE.233.0
   30 CONTINUE
C
      DO 40 I=2,73
      IF (JCLOUD(I).EQ.1) GO TO 40
      IF (PREC1(I).EQ.0.0) GO TO 40
      CLL(I)=CL(I)
      CIRUSL(I)=CIRUS(I)
      CL(I)=1.0
      CIRUS(I)=(T1(I)+T2(I).LE.466.0).AND.(PREC1(I).LE.0.0)
      JCLL(I)=JCLOUD(I)
      JCLOUD(I)=4
      CL4(I)=CL(I)
      NCL4=NCL4+1
   40 CONTINUE
C
C        DETERMINE CLOUD PROPERTIES
C
C        CL = FRACTIONAL CLOUDINESS FOR ALL CLOUD TYPES
C        CLL =FRACTIONAL CLOUDINESS FOR LOWER CLOUD
C        CLR(I) = CL*RADIATION EMISSIVITY (.5 FOR GRAY, 1. FOR BLACK)
C        PCL(I) = 1.-CL(I)
C        PCLR(I) = 1.-CLR(I)
C        CIRUS(I) = FLAG FOR CIRUS(I) CLOUDS
C        GREY = EMISSIVITY FOR RADIATION IN THE CLOUD TYPE
C
C        GREYNESS FACTOR FOR CIRUS
C
      DO 50 I=2,73
      IF (CIRUS(I)) GREY(I)=0.5
      IF (JCLOUD(I).EQ.1.AND.(T1(I)+TTROP(I).LE.466.0)) GREY(I)=0.5
      IF (CIRUSL(I)) GREYL(I)=0.5
   50 CONTINUE
C
C        CLOUDINESS PARAMETERS FOR RADIATION
C
      DO 60 I=2,73
      CLR(I)=CL(I)*GREY(I)
      CLRL(I)=CLL(I)*GREYL(I)
      PCLR(I)=1.0-CLR(I)
      PCLRL(I)=1.0-CLRL(I)
      PCL(I)=1.0-CL(I)
      PCLL(I)=1.0-CLL(I)
   60 CONTINUE
C
      IF (NCL1+NCL2+NCL3+NCL4.LE.0) GO TO 90
C
      DO 70 I=2,73
      IF (JCLOUD(I).EQ.0) GO TO 70
      ALBDOA(I)=0.19
      ABSVTY(I)=0.04
      ALBDOS(I)=0.19
C
C        RWC(I)=1-ABSVTY/(1-ALBDOA)
C
      RWC(I)=0.95061728
      IF (CIRUS(I)) GO TO 70
      ALBDOA(I)=0.46
      ABSVTY(I)=0.2
      ALBDOS(I)=0.46
C
      RWC(I)=.62962963
      IF (JCLOUD(I).NE.1.AND.JCLOUD(I).NE.3) GO TO 70
      ALBDOA(I)=0.6
      ABSVTY(I)=0.3
      ALBDOS(I)=0.6
C
      RWC(I)=0.25
   70 CONTINUE
C
      IF (NCL4.LE.0) GO TO 90
      IF (NCL2.LE.0.AND.NCL3.LE.0) GO TO 90
      DO 80 I=2,73
      IF (JCLOUD(I).NE.4) GO TO 80
      IF (JCLL(I).LT.2) GO TO 80
      ALBDAL(I)=0.19
      ABSVYL(I)=0.04
      ALBDSL(I)=0.19
C
C        RWCL(I)=1-ABSVYL/(1-ALBDAL)
C
      RWCL(I)=0.95061728
      IF (CIRUSL(I)) GO TO 80
      ALBDAL(I)=0.46
      ABSVYL(I)=0.2
      ALBDSL(I)=0.46
C
      RWCL(I)=0.62962963
      IF (JCLL(I).NE.3) GO TO 80
      ALBDAL(I)=0.6
      ABSVYL(I)=0.3
      ALBDSL(I)=0.6
C
      RWCL(I)=0.25
   80 CONTINUE
   90 CONTINUE
      RETURN
      END
C     ******************************
C     ******************************
      SUBROUTINE TRWATR(AMTW,TWR,IL)
C     ******************************
C     ******************************
C
C
C        EMPIRICAL WATER VAPOR TRANSMISSION FUNCTION FOR SOLAR RADIATION
C        DETERMINED BY INTEGRATION OF MC CLATCHEY.S LOWTRAN DATA FOR
C        WATER, AMT=(E-5,E4) GM/CM**2, TWR IS TRANSMISSIVITY OF
C        .ABSORBED PART.
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION AA(74),BB(74),CC(74)
      DIMENSION AMTW(74),TWR(74)
      DIMENSION BND(12),A(13),B(13),C(13)
      DATA BND/6.E-5, 3.E-4, 1.E-3, 6.E-3, 3.E-2, 1.E-1, 6.E-1,
     1         4.E 0, 3.E 1, 2.E 2, 1.E 3, 6.E 3/
      DATA A/8.611328E-5, 4.567886E-4, 1.405538E-3, 3.819889E-3,
     1       1.055894E-2, 2.247760E-2, 4.225268E-2, 7.757599E-2,
     2       1.245152E-1, 1.818462E-1, 2.407783E-1, 2.935484E-1,
     3       3.265483E-1/
      DATA B/2.472894E 1, 1.370112E 1, 7.518070E 0, 3.468863E 0,
     1       1.322240E 0, 5.016671E-1, 1.438218E-1, 2.816388E-2,
     2       5.039994E-3, 8.014325E-4, 1.285009E-4, 1.981110E-5,
     3       5.975114E-6/
      DATA C/-9.923650E 4, -1.228571E 4, -1.957271E 3, -1.969240E 2,
     1       -1.619591E 1, -1.791154E 0, -1.038365E-1, -3.140233E-3,
     2       -8.296123E-5, -1.997084E-6, -5.992159E-8, -1.618743E-9,
     3       -2.016571E-10/
C
      DO 30 I=2,IL
      DO 10 J=1,12
      IF(AMTW(I).LE.BND(J)) GO TO 20
   10 CONTINUE
      J = 13
   20 AA(I)=A(J)
      BB(I)=B(J)
      CC(I)=C(J)
   30 CONTINUE
C
      DO 40 I=1,74
   40 TWR(I)=0.0
C
      DO 100 I=2,IL
      TWR(I)=1.0-(AA(I)+AMTW(I)*(BB(I)+CC(I)*AMTW(I)))/0.366
  100 CONTINUE
C
      DO 200 I=2,IL
  200 TWR(I)=dMAX1(TWR(I),0.0d0)
      RETURN
      END
C     ***************************
C     ***************************
      SUBROUTINE AMTWTR(X,AMT,IL)
C     ***************************
C     ***************************
C
C
C        INVERSE OF WATER VAPOR ABSORPTION FUNCTION FOR SOLAR RADIATION
C        DETERMINED FROM MC CLATCHEY'S LOWTRAN DATA, X IS ABSORPTIVITY
C        OF THE 'ABSORBED PART'
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION X(74),AMT(74)
      DIMENSION ABSPTN(74)
      DIMENSION AA(74),BB(74),CC(74)
      DIMENSION BND(14),A(15),B(15),C(15)
      DATA BND/1.717209E-3, 6.980809E-3, 2.018627E-2, 3.986494E-2,
     1         5.488397E-2, 8.753838E-2, 1.215694E-1, 1.609071E-1,
     2         1.891823E-1, 2.290224E-1, 2.624816E-1, 2.959839E-1,
     3         3.288028E-1, 3.551181E-1/
      DATA A/-1.172763E-6, -2.231221E-6, 4.609903E-4, 1.324878E-2,
     1       9.312623E-2, 6.289253E-1, 4.593607E 0, 2.960411E 1,
     2       1.193145E 2, 6.083588E 2, 2.595276E 3, 1.302157E 4,
     3       6.572026E4, 3.565238E5, 1.128963E6/
      DATA B/2.936614E-2, 3.198488E-2, -7.526538E-2, -1.199613E 0,
     1       -5.180382E 0, -2.326999E 1, -1.119238E 2, -5.121072E 2,
     2       -1.641023E 3, -6.703211E 3, -2.403009E 4, -1.026429E 5,
     3       -4.562541E 5, -2.215934E 6, -6.579729E 6/
      DATA C/1.720888E 1, 1.599123E 1, 2.215919E 1, 4.686455E 1,
     1       9.661239E 1, 2.486978E 2, 7.450842E 2, 2.346497E 3,
     2       5.899362E 3, 1.899839E 4, 5.678314E 4, 2.049658E 5,
     3       7.982183E 5, 3.460250 E6, 9.623603 E6/
C
      DO 10 I=2,IL
      ABSPTN(I)=0.366*X(I)
      AMT(I)=0.0
   10 CONTINUE
C
      DO 30 I=2,IL
      DO 20 J=1,14
      IF (ABSPTN(I).LE.BND(J)) GO TO 25
   20 CONTINUE
      J=15
   25 AA(I)=A(J)
      BB(I)=B(J)
      CC(I)=C(J)
   30 CONTINUE
C
      DO 100 I=2,IL
      AMT(I)=AA(I)+ABSPTN(I)*(BB(I)+CC(I)*ABSPTN(I))
  100 CONTINUE
      RETURN
      END
C     ******************************
C     ******************************
      SUBROUTINE TROZON(AMT,TROZ,IL)
C     ******************************
C     ******************************
C
C
C        EMPIRICAL OZONE TRANSMISSION FUNCTION FOR SOLAR RADIATION
C        DETERMINED BY INTEGRATION OF ACKERMAN'S DATA FOR OZONE
C        AMT =(E-5,30) ATM-CM, TROZ IS TRANSMISSIVITY OF
C        'SCATTERED PART'
C
C
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION AMT(74),TROZ(74)
      DIMENSION ABSPTN(74),AA(74),BB(74),CC(74)
      DIMENSION BND(6),A(7),B(7),C(7)
      DATA BND/1.E-3, 7.E-3, 3.E-2, 2.E-1, 1.E 0, 9.E 0/
      DATA A/2.617776E-7, 1.819089E-4, 3.240676E-3, 8.757822E-3,
     1       1.540986E-2, 2.511835E-2, 7.850732E-2/
      DATA B/1.472231E 0, 1.252737E 0, 4.741989E-1, 1.164579E-1,
     1       4.076455E-2, 2.440325E-2, 1.317787E-2/
      DATA C/-1.220810E 2, -5.785796E 1, -6.295442E 0, -2.259775E-1,
     1       -8.505830E-3, -7.831901E-4, -1.788120E-4/
      DO 10 I=2,IL
      TROZ(I)=0.0
   10 CONTINUE
C
      DO 30 I=2,IL
      DO 20 J=1,6
      IF (AMT(I).LE.BND(J)) GO TO 25
   20 CONTINUE
      J=7
   25 AA(I)=A(J)
      BB(I)=B(J)
      CC(I)=C(J)
   30 CONTINUE
C
      DO 100 I=2,IL
      TROZ(I)=1.0-(AA(I)+AMT(I)*(BB(I)+CC(I)*AMT(I)))/0.634
  100 CONTINUE
C
      DO 200 I=2,IL
      TROZ(I)=dMAX1(0.0d0,TROZ(I))
  200 CONTINUE
      RETURN
      END
C     ****************
C     ****************
      SUBROUTINE VAPOR
C     ****************
C     ****************
C
C
C             COMPUTES EFFECTIVE VAPOR CONTENT
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      LOGICAL SNOW,SEAICE
C
      include 'radvar.fi'
C
      COMMON /MOIST/ EFVT(74),EFVST(74),EFV3(74),EFV2(74),
     1               EFV1(74),EFV0(74),AK(74)
C
      COMMON /WORK/  TEM(74),TEM1(74),TEM2(74),TEM3(74),TEM4(74),
     1               P4E(74),P3E(74),P2E(74),P1E(74),PSE(74),
     2               PTE(74)
C
C************ END OF COMMON ******************
C
         ALOG10(XXX)=DLOG10(XXX)
         EXP(XXX)=DEXP(XXX)
      DO 10 I=2,73
      TEM2(I)=dMAX1(Q3(I),1.d-5)
      TEM4(I)=dMAX1(Q1(I),QST)
   10 CONTINUE
      DO 20 I=2,73
      TEM3(I)=ALOG10(TEM2(I))-QC
      TEM1(I)=ALOG10(TEM4(I))-QC
   20 CONTINUE
C
      DO 30 I=2,73
      AK(I)=((P1L(I)-PC)*TEM1(I)+(P3L(I)-PC)*TEM3(I))/
     1      ((P1L(I)-PC)*(P1L(I)-PC)+(P3L(I)-PC)*(P3L(I)-PC))
   30 TEM(I)=AK(I)+CALFA
C
      DO 40 I=2,73
      P4E(I)=EXP(TEM(I)*P4L(I)*ELOG)
      P3E(I)=EXP(TEM(I)*P3L(I)*ELOG)
      P2E(I)=EXP(TEM(I)*P2L(I)*ELOG)
      P1E(I)=EXP(TEM(I)*P1L(I)*ELOG)
      PSE(I)=EXP(TEM(I)*PC*ELOG)
      PTE(I)=EXP(TEM(I)*ALOGP0*ELOG)
   40 CONTINUE
C
      DO 100 I=2,73
  100 TEM1(I)=EFVCON*CALFA/(TEM(I)*PSE(I))
C
      DO 200 I=2,73
      EFV3(I)=TEM1(I)*(P4E(I)-P3E(I))
      EFV2(I)=TEM1(I)*(P4E(I)-P2E(I))
      EFV1(I)=TEM1(I)*(P4E(I)-P1E(I))
      EFV0(I)=TEM1(I)*(P4E(I)-PTE(I))
      EFVST(I)=TEM1(I)*(P4E(I)-PSE(I))
      EFVT(I)=EFVST(I)+EFVCON
  200 CONTINUE
      RETURN
      END
C     *****************
C     *****************
      SUBROUTINE FMAFMB
C     *****************
C     *****************
C
C
C             COMPUTES COEFFICIENTS OF TRANSMISSION FUNCTIONS
C
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
      COMMON /WORK/ AP(216),ZAK(216),TEM1(216),TEM2(216),TEM3(216),
     1              TEM4(216),P(216),Z(216),XA(216),XB(216),
     2               ZA(216),ZB(216),ZC(216),AA(216),AB(216),CA(216),
     3               CB(216),CZ(216)
C
C************ END OF COMMON ******************
C
C
C        EMPIRICAL FUNCTION OF P AND Q FOR CALCULATION OF THE LINEAR
C        INTERPOLATION FACTOR USED IN THE COMPUTATION OF THE ADJACENT
C        LAYER TRANSMISSION FUNCTIONS
C
C        NOTE - AP=ALOG10(P)  Z=ALOG10(Q) + 3
C
         ABS(XXX)=DABS(XXX)
C
      DO 10 I=2,73
      AP(I-1)=ALOGP0
      AP(71+I)=P2L(I)
      AP(143+I)=P4L(I)
      P(I-1)=PTROP
      P(71+I)=PL2(I)
      P(143+I)=P4(I)
      ZAK(I-1)=AK(I)
      ZAK(71+I)=AK(I)
      ZAK(143+I)=AK(I)
   10 CONTINUE
      DO 20 I=1,216
      Z(I)=QC+ZAK(I)*(AP(I)-PC)+3.0
   20 CONTINUE
C
      DO 30 I=1,216
      AA(I)=76.625
      AB(I)=28.39
      CA(I)=60.81
      CB(I)=22.53
   30 CONTINUE
      DO 40 I=1,216
      IF (AP(I).LT.2.7) GO TO 40
      AA(I)=61.86
      AB(I)=22.92
      CA(I)=42.59
      CB(I)=15.78
   40 CONTINUE
C
      DO 50 I=1,216
   50 CZ(I)=-0.041+0.021*Z(I)
C
      DO 60 I=1,216
   60 CZ(I)=dMIN1(CZ(I),-0.06d0)
C
      DO 70 I=1,216
      XA(I)=AA(I)-AB(I)*AP(I)
      XB(I)=-CA(I)+CB(I)*AP(I)
      ZA(I)=Z(I)-0.105*XA(I)
      ZB(I)=Z(I)-0.105*XB(I)
      ZC(I)=ABS(Z(I)+2.5)
   70 CONTINUE
C
      DO 100 I=1,216
      TEM1(I) =  -0.09*XA(I)+2.57+ZA(I)*(.233+ZA(I)*(.18+.027*ZA(I)))
      TEM3(I) = -1.66+1.76*AP(I)+Z(I)*(.3+Z(I)*(.28+.04*Z(I)))
      TEM2(I)=.01*(-0.09*XB(I)+1.42+ZB(I)*(.48+ZB(I)*(.16+.011*ZB(I)))
     1 + (0.08+(0.371-0.102*AP(I))*(Z(I)+2.1))*(ZAK(I)-3.)  )
      TEM4(I)=.01*(-0.197+0.0002*P(I)+ZC(I)*(.0812-ZC(I)*
     *   (.045-.02334*ZC(I)))
     1 + dMIN1(-0.041d0+0.021d0*Z(I),-0.06d0)*(ZAK(I)-3.0d0) )
  100 CONTINUE
C
      DO 200 I=2,73
      EST0(I)=TEM1(I-1)+TEM2(I-1)*(PTROP-PSTQ)
      E10(I)=TEM3(I-1)+TEM4(I-1)*(PL1(I)-PTROP)
      E20(I)=TEM3(I-1)+TEM4(I-1)*(PL2(I)-PTROP)
      E02(I)=TEM1(71+I)+TEM2(71+I)*(PL2(I)-PTROP)
      E12(I)=TEM1(71+I)+TEM2(71+I)*(PL2(I)-PL1(I))
      E32(I)=TEM3(71+I)+TEM4(71+I)*(PL3(I)-PL2(I))
      E42(I)=TEM3(71+I)+TEM4(71+I)*(P4(I)-PL2(I))
      E24(I)=TEM1(143+I)+TEM2(143+I)*(P4(I)-PL2(I))
      E34(I)=TEM1(143+I)+TEM2(143+I)*(P4(I)-PL3(I))
  200 CONTINUE
      RETURN
      END
C     **************
C     **************
      SUBROUTINE TRA
C     **************
C     **************
C
C
C             COMPUTES TOTAL AND MEAN TOTAL TRANSMISSION FUNCTIONS
C
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
      COMMON /WORK/ X(1080),Z(1080),TRS(1080),
     1              TA(1080),TB(1080),TC(1080),
     2              TD(219),CTCO(74)
C
C************ END OF COMMON ******************
C
C
C        EMPIRICAL FIT OF THE WEIGHTED MEAN TRANSMISSION
C        FUNCTION  TAU, KATAYAMA.S METHOD E
C
C        NOTE - Z = ALOG10(X)
C
         ALOG10(XXX)=DLOG10(XXX)
C
      DO 10 I=2,73
      X(I-1)=EFVT(I)
      X(I+71)=EFVT(I)-EFV2(I)
      X(I+143)=EFVT(I)-EFV0(I)
      X(I+215)=EFV2(I)-EFV3(I)
      X(I+287)=EFV0(I)-EFV3(I)
      X(I+359)=EFV3(I)
      X(I+431)=EFV0(I)-EFV1(I)
      X(I+503)=EFV0(I)
      X(I+575)=EFVST(I)
      X(I+647)=EFV2(I)
      X(I+719)=EFV1(I)
      X(I+791)=EFV1(I)-EFV2(I)
      X(I+863)=EFV0(I)-EFV2(I)
      X(I+935)=EFVST(I)-EFV2(I)
      X(I+1007)=EFVST(I)-EFV0(I)
   10 CONTINUE
C
      DO 50 I=1,1080
   50 Z(I)=ALOG10(X(I))
C
      DO 100 I=1,1080
      TA(I)=0.373-Z(I)*(0.2595+0.0275*Z(I))
      TB(I)=0.373-Z(I)*(0.274-0.035*Z(I))
      TC(I)=1.0/(1.0+298.7*X(I))
  100 CONTINUE
C
C        TRS = B IF X GT 1
C            = A IF X LE 1
C
      DO 110 I=1,1080
      TRS(I)=CVMGP(TA(I),TB(I),1.0-X(I))
C
C        TRS = C IF LOG(X) LT -4
C
      TRS(I)=CVMGP(TRS(I),TC(I),Z(I)+4.0)
  110 CONTINUE
C
      DO 120 I=1,1080
      TRS(I)=dMIN1(TRS(I),1.0d0)
      TRS(I)=dMAX1(TRS(I),0.0d0)
  120 CONTINUE
C
      DO 130 I=2,73
      CTCO(I)=1.09-9.E-5*P4(I)  !for CO2 transmission function
  130 CONTINUE
C
      DO 200 I=2,73
      TRT4(I)=TRS(I-1)*CTCO(I)*TCT4
      TRT2(I)=TRS(I+71)*CTCO(I)*TCT2
      TRT0(I)=TRS(I+143)*CTCO(I)*TCT0
      TR23(I)=TRS(I+215)*CTCO(I)*TC23
      TR03(I)=TRS(I+287)*CTCO(I)*TC03
      TR34(I)=TRS(I+359)*CTCO(I)*TC34
      TR01(I)=TRS(I+431)*CTCO(I)*TC01
      TR04(I)=TRS(I+503)*CTCO(I)*TC04
      TRST4(I)=TRS(I+575)*CTCO(I)*TCST4
      TR24(I)=TRS(I+647)*CTCO(I)*TC24
      TR14(I)=TRS(I+719)*CTCO(I)*TC14
      TR12(I)=TRS(I+791)*CTCO(I)*TC12
      TR02(I)=TRS(I+863)*CTCO(I)*TC02
      TRST2(I)=TRS(I+935)*CTCO(I)*TCST2
      TRST0(I)=TRS(I+1007)*CTCO(I)*TCST0
  200 CONTINUE
C
      DO 210 I=1,216
      TA(I)=0.254-Z(I)*(0.226+0.007*Z(I))
      TB(I)=0.254-Z(I)*(0.2224-0.0444*Z(I))
      TC(I)=1.0/(1.0+2.56*X(I)**0.39)
      TD(I)=0.194-Z(I)*(0.297+0.028*Z(I))
  210 CONTINUE
C
C        TRS = B IF X GT 1
C            = A IF X LE 1
C
      DO 300 I=1,216
      TRS(I)=CVMGP(TA(I),TB(I),1.0-X(I))
C
C        TRS(I) = D IF X LT 0.03
C
      TRS(I)=CVMGP(TRS(I),TD(I),X(I)-0.03)
C
C        TRS = C IF LOG(X) LT -4
C
      TRS(I)=CVMGP(TRS(I),TC(I),Z(I)+4.0)
  300 CONTINUE
C
      DO 310 I=1,216
      TRS(I)=dMIN1(TRS(I),1.0d0)
      TRS(I)=dMAX1(TRS(I),0.0d0)
  310 CONTINUE
C
      DO 400 I=2,73
      T2RT4(I)=TRS(I-1)*CTCO(I)*TCT4
C
C
C
      T2RT2(I)=TRS(I+71)*CTCO(I)*TCT2
      T2RT0(I)=TRS(I+143)*CTCO(I)*TCT0
  400 CONTINUE
      RETURN
      END
      FUNCTION CVMGP(X1,X2,X3)
         REAL*8 X1,X2,X3,CVMGP
      CVMGP=X1
      IF (X3.LT.0.) CVMGP=X2
      RETURN
      END
