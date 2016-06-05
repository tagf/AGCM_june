C     ****************
C     ****************
      SUBROUTINE FORCE
C     ****************
C     ****************
C
C
C             COMPUTES FORCING TERMS,DAILY ACCUMULATIONS AND ZONAL MEANS
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
C        ACCUMULATED VARIABLES IN COMMON
C
	include 'accum.fi'
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C        COMMON FOR RADIATION AND CLOUD SUB. RADCAL, COMP3
C
      COMMON /day_end/ end_of_day
      logical end_of_day

      LOGICAL SNOW,SEAICE
C
      include 'radvar.fi'
C
      include 'comp.fi'

      include 'ice.fi'
      include 'acc_ice.fi'
	include 'dir.fi'
cIIIIIIIIIIIIIIIIIIIIIIIIIII
      include 'StatIce.fi'
cIIIIIIIIIIIIIIIIIIIIIIIIIII
C
      common/glacc/GMTACC,GMKEAC,GMRACC,namgl 
C************ END OF COMMON ******************
c      DIMENSION AlbTemp(74,46) !for surfer
c
      DIMENSION ASCZtmp(74,46)   !test scosz
C
      DIMENSION BS(74),BA(74),ZON(72,74)
      DIMENSION AN0(74),FD(74,46),FDU(74,46)
      DIMENSION FDT(74),fdsum(74)
      DIMENSION SLP(72)
      DIMENSION SUMFD(74)
cIIIIIIIIIIIIIIIIIIIIIIIIIII
	DIMENSION xtemp(11)
	DIMENSION tmpS4(74,46),tmpR4(74,46),tmpF4(74,46), ! for surfer
     1	tmpE4(74,46),tmpBS(74,46)
cIIIIIIIIIIIIIIIIIIIIIIIIIII
!UUUUUUU For CO2 Usatuk Jl 09
      DIMENSION T4_max(72,46),T4_min(72,46),S4_max(72,46),
     1 Prec_day(72,46),Wind_day(72,46),Q4_day(72,46),
     2 CL_day(72,46),U1_day(72,46),U2_day(72,46),V1_day(72,46),
     3 V2_day(72,46),SD_day(72,46),SLP_day(72,46)
      integer step_count
!UUUUUUU
         EXP(XXX)=DEXP(XXX)
      DATA ISET /0/
      DATA CJW/0.4846d0/ !convert  ly/day to W
      GMR=0.0
C
         NAV=NAV+1
      step_count=step_count+1  !UUUUUUU For CO2 Usatuk Jl 09
       DO 7 I=1,72
       DO 7 N=1,74
  7    ZON(I,N)=0.
C
      IF (ISET.NE.0) GO TO 10
      ISET=1
      DTH=DT/3600.0
      NHIS=TAUH/DTH+0.01
      AREA=DXYP(1)
      DO 8 J=2,JM
    8 AREA=AREA+DXYP(J)
      AREA=AREA*FIM
 !UUUUUUU For CO2 Usatuk Jl 09
      T4_max=-100.
      T4_min=1000.
      S4_max=-100.
      Prec_day=0.
      Q4_day=0.
      CL_day=0.
      U1_day=0.
      U2_day=0.
      V1_day=0.
      V2_day=0.
      SD_day=0.
      SLP_day=0.
      Wind_day=0.
      step_count=0
   !   if (MNTHDY==1) then ! NumProcs = 1  - Only!
  !     OPEN (333, FILE = TRIM(BaseDir) // WorkDir //
   !    WRITE (333,999) MONTH,MNTHDY,TOFDAY
 !999    FORMAT (1X,I3,1H/,I2,1H/,F6.2)   
    !  endif    
  !UUUUUUU
   10 CONTINUE
C
C        COMPUTE ZONAL MAX AND MIN, AND GLOBAL MAX AND MIN AND
C        (I,J) POINTS FOR TEMPERATURE, ZONAL WIND VELOCITY AND
C        MERIDIONAL WIND VELOCITY AT LEVELS 1 AND 3.
C
	 if (dabs(tofday-24.0d0).le.0.1) THEN !end of day

  !    CALL WRITIJ(T(1,1,1),1.d0,'T1  ',0)
      CALL AMXMN (T(1,1,1),ZON(1,63),TMX1,ITX1,JTX1,
     1  ZON(1,64),TMN1,ITN1,JTN1)
      CALL AMXMN (T(1,1,2),ZON(1,65),TMX2,ITX2,JTX2,
     1  ZON(1,66),TMN2,ITN2,JTN2)
!      CALL WRITIJ(U(1,1,1),1.d0,'U1  ',2)
      CALL AMXMN (U(1,1,1),ZON(1,67),UMX1,IUX1,JUX1,
     1  ZON(1,68),UMN1,IUN1,JUN1)
      CALL AMXMN (U(1,1,2),ZON(1,69),UMX2,IUX2,JUX2,
     1  ZON(1,70),UMN2,IUN2,JUN2)
      CALL AMXMN (V(1,1,1),ZON(1,71),VMX1,IVX1,JVX1,
     1  ZON(1,72),VMN1,IVN1,JVN1)
      CALL AMXMN (V(1,1,2),ZON(1,73),VMX2,IVX2,JVX2,
     1  ZON(1,74),VMN2,IVN2,JVN2)
C
C        PRINT OUT GLOBAL MAX AND MIN AND (I,J) POINTS
C        FOR MODEL MONITORING
C
      ITMX1=TMX1-TICE
      ITMN1=TMN1-TICE
      ITMX2=TMX2-TICE
      ITMN2=TMN2-TICE
      IUMX1=UMX1
      IUMX2=UMX2
      IUMN1=UMN1
      IUMN2=UMN2
      IVMX1=VMX1
      IVMX2=VMX2
      IVMN1=VMN1
      IVMN2=VMN2
c      SMN8=1.E6                   !stop, if U or T are large
      if (abs(UMX1).ge.300..or.ITMN1.le.-300.) then 
         WRITE (6,2) TAU,MONTH,GMKEAC,GMTACC,GMRACC,SDEDY
         print *,UMX1,ITMN1
 2       FORMAT (1X,'** ERR **',F9.2,2X,I3,2X,3e9.2,2x,i3)
c          Play music when stop
        CALL BEEPQQ(2000, 200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(1000,200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(500, 200)
        CALL BEEPQQ(2000, 200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(1000,200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(500, 200)
c
	   stop ' TOO LARGE NUMBERS '
	endif    
c      DO 15 J=1,46 !for print snoamt
c      DO 15 I=2,73
c      IF (ISFTYP(I,J).NE.8) GO TO 15
c      SMN8=dMIN1(SNOAMT(I,J),SMN8)
c   15 CONTINUE

c      ISMN8=SMN8
c        IF (mnthdy.eq.1)  then     ! Every month output
        IF (end_of_day)  then     ! Every day output
         WRITE (7,99) ITMX1,ITX1,JTX1,ITMN1,ITN1,JTN1,
     1                ITMX2,ITX2,JTX2,ITMN2,ITN2,JTN2,
     2                IUMX1,IUX1,JUX1,IUMN1,IUN1,JUN1,
     3                IUMX2,IUX2,JUX2,IUMN2,IUN2,JUN2
	  ENDIF
cc      WRITE (6,99) ITMX1,ITX1,JTX1,ITMN1,ITN1,JTN1, !убрал печать на экран
cc     1             ITMX2,ITX2,JTX2,ITMN2,ITN2,JTN2,
cc     2             IUMX1,IUX1,JUX1,IUMN1,IUN1,JUN1,
cc     3             IUMX2,IUX2,JUX2,IUMN2,IUN2,JUN2
c     4             IVMX1,IVX1,JVX1,IVMN1,IVN1,JVN1,
c     5             IVMX2,IVX2,JVX2,IVMN2,IVN2,JVN2
c     6            ,ISMN8
	 ENDIF
   99 FORMAT (4(I3,1H(,2I2,2H)/,I4,1H(,2I2,1H)),1X,I4)
      DO 1000 J=1,JM
C***
         IF (J.EQ.1) GO TO 98
         DO 97 I=2,73
         UACC(I-1,J,1)=UACC(I-1,J,1)+U(I,J,1)
         UACC(I-1,J,2)=UACC(I-1,J,2)+U(I,J,2)
         VACC(I-1,J,1)=VACC(I-1,J,1)+V(I,J,1)
  97     VACC(I-1,J,2)=VACC(I-1,J,2)+V(I,J,2)
  98     CONTINUE
         DO 96 I=2,73
         PACC  (I-1,J  )=PACC  (I-1,J  )+P     (I,J  )
         TACC  (I-1,J,1)=TACC  (I-1,J,1)+T     (I,J,1)
         TACC  (I-1,J,2)=TACC  (I-1,J,2)+T     (I,J,2)
         QWACC (I-1,J,1)=QWACC (I-1,J,1)+QW    (I,J,1)
         QWACC (I-1,J,2)=QWACC (I-1,J,2)+QW    (I,J,2)
         GWACC (I-1,J  )=GWACC (I-1,J  )+GW    (I,J  )
         GTACC (I-1,J  )=GTACC (I-1,J  )+GT    (I,J  )
         SNOACC(I-1,J  )=SNOACC(I-1,J  )+SNOAMT(I,J  )
         TS ACC(I-1,J  )=TS ACC(I-1,J  )+TS    (I,J  )
  96     SD ACC(I-1,J  )=SD ACC(I-1,J  )+SD    (I,J  )
C***
C
      CALL COMP3 (J)
C
      DO 20 I=2,73
      APLS(I-1,J)=APLS(I-1,J)+SP(I)*PREC3(I)
      APCV(I-1,J)=APCV(I-1,J)-SP(I)*(CQ1(I)+CQ3(I)+PCQ1(I)+PCQ3(I))
C***
      APRES(I-1,J)=APRES(I-1,J)+APLS(I-1,J)+APCV(I-1,J)
   20 ACLOUD(I-1,J)=ACLOUD(I-1,J)+CL(I)
C***
      DO 30 I=2,73
      ASNF(I-1,J)=ASNF(I-1,J)+SNFAL(I)
   30 ACNO(I-1,J)=ACNO(I-1,J)+AS1(I)+AS3(I)+S4(I)-R0(I)
      DO 40 I=2,73
      BS(I)=S4(I)-(1.0-FLD(I))*R4(I)-FLD(I)*R4C(I)
   40 ACNS(I-1,J)=ACNS(I-1,J)+BS(I)
      DO 50 I=2,73
      BS(I)=BS(I)-(1.0-FLD(I))*(F4(I)+EVAL(I)*E4(I)*DAY+CICE(I))
     1            -FLD(I)*(F4C(I)+600.*E4C(I)*DAY)
   50 ACBS(I-1,J)=ACBS(I-1,J)+BS(I)
c********* for surfer files
      DO 51 I=2,73 ! W/m**2
	tmpS4(I,J)=S4(I)*CJW
	tmpR4(I,J)=R4(I)*CJW
	tmpF4(I,J)=F4(I)*CJW
	tmpE4(I,J)=EVAL(I)*E4(I)*DAY*CJW
   51	tmpBS(I,J)=BS(I)*CJW
c*********
      DO 60 I=2,73
      BA(I)=(H1(I)+H3(I))/(COE1*COLMR(I))
   60 ACBA(I-1,J)=ACBA(I-1,J)+BA(I)
      DO 70 I=2,73
      ABEA(I-1,J)=ABEA(I-1,J)+BA(I)+BS(I)+CICE(I)*(1.0-FLD(I))
   70 CONTINUE
      DO 90 I=2,73
      ASCZ(I-1,J)=ASCZ(I-1,J)+SCOSZ(I)
      ASCZtmp(I,J)=SCOSZ(I)  !test cosz
   90 ARET(I-1,J)=ARET(I-1,J)+RETOT(I)
      DO 95 I=2,73
      AT4(I-1,J)=AT4(I-1,J)+S4(I)
      hIceAcc(I-1,J)=hIceAcc(I-1,J)+GHIce(i,j) !sea ice
   95 GIceCompAcc(I-1,J)=GIceCompAcc(I-1,J)+GIceComp(i,j)
C***                    ****T4****
c      do i=2,73 !for Surfer
c      AlbTemp(I,J) =  ALs(I)
c	enddo
c
c         DO 91 I=2,73
c         TAUUAC(I-1,J)=TAUUAC(I-1,J)+TAUU(I,J)
c   91    TAUVAC(I-1,J)=TAUVAC(I-1,J)+TAUV(I,J)
C***
C
C        ZONAL AVERAGES
C
      DO 100 I=2,73
      ZON(I-1,1)=CL1(I)
      ZON(I-1,2)=CL2(I)
      ZON(I-1,3)=CL3(I)
      ZON(I-1,4)=CL4(I)
      ZON(I-1,5)=CL(I)
      ZON(I-1,6)=(1.0-FLD(I))*E4(I)+FLD(I)*E4C(I)
      ZON(I-1,7)=(1.0-FLD(I))*F4(I)+FLD(I)*F4C(I)
      ZON(I-1,8)=AS1(I)
      ZON(I-1,9)=AS3(I)
      ZON(I-1,10)=AS1(I)/SP(I)
      ZON(I-1,11)=AS3(I)/SP(I)
      ZON(I-1,12)=R0(I)
      ZON(I-1,13)=R2(I)
      ZON(I-1,14)=R4(I)
      ZON(I-1,15)=R0(I)/SP(I)
      ZON(I-1,16)=R2(I)/SP(I)
      ZON(I-1,17)=R4(I)/SP(I)
      ZON(I-1,18)=S4(I)
      ZON(I-1,19)=RETOT(I)
      ZON(I-1,20)=SCOSZ(I)
      ZON(I-1,21)=PREC1(I)*SP(I)
      ZON(I-1,22)=PREC3(I)*SP(I)
      ZON(I-1,23)=CQ1(I)*SP(I)
      ZON(I-1,24)=CQ3(I)*SP(I)
      ZON(I-1,25)=PCQ1(I)*SP(I)
      ZON(I-1,26)=PCQ3(I)*SP(I)
      ZON(I-1,27)=CT1(I)
      ZON(I-1,28)=CT3(I)
      ZON(I-1,29)=PCT1(I)
      ZON(I-1,30)=PCT3(I)
      ZON(I-1,31)=ZON(I-1,7)/SP(I)
      ZON(I-1,32)=RH1(I)
      ZON(I-1,33)=RH3(I)
      ZON(I-1,34)=T4(I)
      ZON(I-1,35)=ALS(I)
      ZON(I-1,36)=Q4(I)
      ZON(I-1,37)=TG(I)
      ZON(I-1,38)=WET(I)
      ZON(I-1,39)=SNR(I)
      ZON(I-1,40)=SD(I,J)
  100 CONTINUE
      IF (J.EQ.1) GO TO 120
      DO 110 I=2,73
      ZON(I-1,41)=U(I,J,1)
      ZON(I-1,42)=U(I,J,2)
      ZON(I-1,43)=V(I,J,1)
      ZON(I-1,44)=V(I,J,2)
      ZON(I-1,58)=BETA (I)*SP(I)*SP(I)
      ZON(I-1,59)=PIV(I,J,2)+PIV(I,J,1)
      ZON(I-1,60)=PIV(I,J,2)-PIV(I,J,1)
      ZON(I-1,61)=US(I,1)
      ZON(I-1,62)=VS(I,1)
  110 CONTINUE
  120 DO 130 I=2,73
      ZON(I-1,45)=QW(I,J,1)
      ZON(I-1,46)=QW(I,J,2)
      ZON(I-1,47)=SP(I)*(QW(I,J,1)+QW(I,J,2))
      ZON(I-1,48)=EDV(I)*2000.
      ZON(I-1,49)=SMELT(I)
      ZON(I-1,53)=SNFAL(I)
      ZON(I-1,54)=EVAL(I)*ZON(I-1,6)
      ZON(I-1,55)=CICE(I)
      ZON(I-1,56)=PREC1(I)
      ZON(I-1,57)=PREC3(I)
  130 CONTINUE
c
C*****IIIIIIIIIIIIIIIII
	if (j.eq.jpoint(1)) then
	do ii=1,10
       xtemp(ii)=xIce(ii,1)
	enddo
	i=ipoint(1)
      xtemp(1)=xtemp(1)+SP(I)*PREC3(I)-SP(I)*(CQ1(I)+CQ3(I)
	1                                            +PCQ1(I)+PCQ3(I))
      xtemp(2)=xtemp(2)+CL(I)
      xtemp(3)=xtemp(3)+SNFAL(I)
      xtemp(4)=xtemp(4)+AS1(I)+AS3(I)+S4(I)-R0(I)
      BS(I)=S4(I)-(1.0-FLD(I))*R4(I)-FLD(I)*R4C(I)
      xtemp(5)=xtemp(5)+BS(I)
      BS(I)=BS(I)-(1.0-FLD(I))*(F4(I)+EVAL(I)*E4(I)*DAY+CICE(I))
     1            -FLD(I)*(F4C(I)+600.*E4C(I)*DAY)
      xtemp(6)=xtemp(6)+BS(I)
      xtemp(7)=xtemp(7)+S4(I)
      xtemp(8)=xtemp(8)+(1.0-FLD(I))*F4(I)+FLD(I)*F4C(I)
      xtemp(9)=xtemp(9)+EVAL(I)*(1.0-FLD(I))*E4(I)+FLD(I)*E4C(I)
	xtemp(10)=xtemp(10)+T4(I)-273.
	do ii=1,10
       xIce(ii,1)=xtemp(ii)
	enddo
	endif
C*****IIIIIIIIIIIIIIIII
C         ZONALLY ACCUMULATE THE VALUES
C
      DO 150 I=1,72
      DO 150 N=53,62
  150 ZONAVG(J,N)=ZONAVG(J,N)+ZON(I,N)
      DO 200 I=1,72
      DO 200 N=1,49
  200 ZONAVG(J,N)=ZONAVG(J,N)+ZON(I,N)
C         NET RADIATION AT TOP OF ATMOSPHERE
C
c      IF (MOD(NSTEP,NHIS).NE.0) GO TO 400
      DO 300 I=2,73
  300 AN0(I)=AS1(I)+AS3(I)+S4(I)-R0(I)
      GMR=SUMma(AN0)*DXYP(J)+GMR
  400 CONTINUE
C
!UUUUUUU For CO2 Usatuk Jl 09
      DO I=2,73
       if (t4(i).gt.T4_max(i-1,j)) T4_max(i-1,j)=t4(i)
       if (t4(i).lt.T4_min(i-1,j)) T4_min(i-1,j)=t4(i)
       if (S4(i).gt.S4_max(i-1,j)) S4_max(i-1,j)=S4(i)
       Prec_day(i-1,j)=Prec_day(i-1,j)+SP(I)*PREC3(I)-SP(I)*
     1  (CQ1(I)+CQ3(I)+PCQ1(I)+PCQ3(I))
       Q4_day(i-1,j)=Q4_day(i-1,j)+Q4(I)
       Q4_AGCM(i,j)=Q4(I) ! GGGGGGGGGGG for GLDSTN
       TOTALP_AGCM(i,j)=TOTALP(I) !  prec  GGGGGGGGGGG for GLDSTN
       CL_day(i-1,j)=CL_day(i-1,j)+CL(I) 
       U1_day(i-1,j)=U1_day(i-1,j)+U(I,J,1)
       U2_day(i-1,j)=U2_day(i-1,j)+U(I,J,2)
       V1_day(i-1,j)=V1_day(i-1,j)+V(I,J,1)
       V2_day(i-1,j)=V2_day(i-1,j)+V(I,J,2)  
       SD_day(i-1,j)=SD_day(i-1,j)+SD(I,J)/DXYP(j)
       Wind_day(i-1,j)=Wind_day(i-1,j)+dsqrt(US(I,1)*US(I,1)+
     1  VS(I,1)*VS(I,1))
      END DO
    
!UUUUUUU 
 1000 CONTINUE
	nnn=nnn+1
C
C        END OF COMP3 LOOP
C
      DO 1100 N=63,74,2
      DO 1100 J=1,46
      ZONAVG(J,N)=dMAX1(ZONAVG(J,N),ZON(J,N))
 1100 ZONAVG(J,N+1)=dMIN1(ZONAVG(J,N+1),ZON(J,N+1))
      DO 2000 I=2,73
      DO 2000 J=1,46
      ZONAVG(J,50)=ZONAVG(J,50)+T(I,J,1)
      ZONAVG(J,51)=ZONAVG(J,51)+T(I,J,2)
 2000 CONTINUE
      DO 2200 J=1,JM
      DO 2100 I=2,73
      SLP(I-1)=(P(I,J)+PTROP)*
     1 EXP(PHIS(I,J)/(RGAS*(1.5*T(I,J,2)-0.5*T(I,J,1)+FLR*PHIS(I,J))))
      ASLP(I-1,J)=ASLP(I-1,J)+SLP(I-1)
      SLP_day(i-1,j)=SLP_day(i-1,j)+SLP(I-1) !UUUU For CO2 Usatuk 03.10
 2100 CONTINUE
      DO 2200 I=2,73
 2200 ZONAVG(J,52)=ZONAVG(J,52)+SLP(I-1)
 
!UUUUUUU For CO2 Usatuk Jl 09
      go to 2201
      if (end_of_day) then
!        WRITE (333,998) step_count
   !     WRITE (333,998) sdedy
       WRITE (333,999) MONTH,MNTHDY,TOFDAY, sdedy
 999    FORMAT (1X,I3,1H/,I2,1H/,F6.2,2X,I3)
       DO I=2,73
        DO J=1,46
         CL_day(i-1,j)=CL_day(i-1,j)/step_count
         if (CL_day(i-1,j).gt.1.0) CL_day(i-1,j)=1.0
          WRITE (333,998) i-1,j,(T4_max(i-1,j)-273.),
     1   (T4_min(i-1,j)-273.),CJW*S4_max(i-1,j),
     2   Prec_day(i-1,j),Q4_day(i-1,j)/step_count,
     3   CL_day(i-1,j),U1_day(i-1,j)/step_count,
     4   U2_day(i-1,j)/step_count,V1_day(i-1,j)/step_count,
     5   V2_day(i-1,j)/step_count,SD_day(i-1,j)/step_count,
     6   SLP_day(i-1,j)/step_count,Wind_day(i-1,j)/step_count
!      WRITE (333,998) i-1,j,TS(I,J)-273.
        END DO   
      END DO   
  998    FORMAT (1X,2I3,2F6.1,F7.1,F6.2,E9.2,F5.2,4F6.1,E10.2,F7.1, 
     1   F5.1)
 !     pause 'Usat'
      T4_max=-100.
      T4_min=1000.
      S4_max=-100.
      Prec_day=0.
      Q4_day=0.
      CL_day=0.
      U1_day=0.
      U2_day=0.
      V1_day=0.
      V2_day=0.
      SD_day=0.
      SLP_day=0.
      Wind_day=0.
      step_count=0
    !  if (MNTHDY==1) then ! NumProcs = 1  - Only!
    !   close (333)
   !    stop 'Month end'
   !    OPEN (333, FILE = TRIM(BaseDir) // WorkDir //
   !  1 ' /Usatuk.txt', ACCESS='APPEND')
      endif
 2201  continue        
!UUUUUUU 
C
C        GLOBAL MEAN K.E., T, NET RADIATION AT LEVEL 0
C
c      IF (MOD(NSTEP,NHIS).NE.0) GO TO 5000
      DO 3000 J=1,JM
      DO 3000 I=2,74
 3000 FD(I,J)=P(I,J)*DXYP(J)
      DO 3100 I=1,74
      FDT(I)=0.0
 3100 SUMFD(I)=0.0
      DO 3300 J=1,JM
      DO 3200 I=2,73
      FDT(I)=FD(I,J)*(T(I,J,1)+T(I,J,2))+FDT(I)
 3200 SUMFD(I)=FD(I,J)+SUMFD(I)
 3300 CONTINUE
      GMT=SUMma(FDT)/SUMma(SUMFD)*0.5-TICE
C
      DO 4000 I=2,74
      FD(I,1)=FD(I,1)+FD(I,1)
 4000 FD(I,JM)=FD(I,JM)+FD(I,JM)
C
      DO 4010 J=1,JM
      DO 4010 I=2,73
 4010 FDU(I,J)=FD(I,J)+FD(I+1,J)
      DO 4100 J=2,JM
      DO 4100 I=2,73
 4100 FD(I,J)=(FDU(I,J)+FDU(I,J-1))*(U(I,J,1)*U(I,J,1)+
     1   U(I,J,2)*U(I,J,2)+V(I,J,1)*V(I,J,1)+V(I,J,2)*V(I,J,2) )
      GMKE=0.0
      DO 4300 J=2,JM
	do 4301 i=2,73
 4301	fdsum(i)=fd(i,j)
 4300 GMKE=GMKE+SUMma(FDsum)
C
      GMKE=0.25E-3*GMKE/(4.0*GRAV*AREA)
      GMR=GMR/AREA

C***     DAYLY MEANS TGLOB,K.E.,N0
         GMTACC=GMT+GMTACC
         GMKEAC=GMKE+GMKEAC
         GMRACC=GMR+GMRACC
	   namgl=namgl+1
C***
c	call SurFile ('ALBEDO',AlbTemp,0.d0,'ALBEDO2')
c      STOP 'Normal end'
c	call SurFile ('Snow',SNOAMT,0.d0,'snow')
c	call SurFile ('S4',tmpS4,0.d0,'S4')
c	call SurFile ('R4',tmpR4,0.d0,'R4')
c	call SurFile ('F4',tmpF4,0.d0,'F4')
c	call SurFile ('E4',tmpE4,0.d0,'E4')
c	call SurFile ('BS',tmpBS,0.d0,'BS')
 5000 CONTINUE
      RETURN
      END
