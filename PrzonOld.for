C     ****************
C     ****************
      SUBROUTINE PRZON(lun)
C     ****************
C     ****************
C
C        PRINTS OUT ZONAL MEANS
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)
C
C        ACCUMULATED VARIABLES IN COMMON
C
      COMMON /ACCUM/ PACC(72,46),TACC(72,46,2),GTACC(72,46),
     1            TSACC(72,46),APRES(72,46),ACLOUD(72,46) ,ZONAVG(46,74)
      dimension zero(72*7+74,46)
      equivalence (pacc(1,1),zero(1,1))
C
C************ END OF COMMON ******************
			character*4 title
      real*8 ZONAV1(46,74)
      DIMENSION TITLE(68)
      DIMENSION AN0(46),ALFAP(46),ANS(46),BS(46),BA(46),BEA(46)
C
      DATA TITLE
     1/' CL1',' CL2',' CL3',' CL4','  CL','  E4','  F4',' AS1',
     1 ' AS3','AS1P','AS3P','  R0','  R2','  R4','R0:P',
     2 'R2:P','R4:P','  S4','RETO','COSZ','PC1P','PC3P','CQ1P',
     2 'CQ3P','PCQ1','PCQ3',' CT1',' CT3','ACT1',' PCT',
     3 'F4SP',' RH1',' RH3','  T4',' ALS','  Q4','  TG',' WET',
     3 ' SNR','  SD','  U1','  U3','  V1','  V3',' QW1',
     4 ' QW3','Q1Q3',' EDV','SMLT','  T1','  T3',' SLP','SFAL',
     4 ' LES','CICE','PRC1','PRC3','BETA',' PIV','V3V1',
     5 '  US','  VS','  N0','ALBP','  NS','  BS','  BA',' BEA'/
      DATA CJW/0.4846E0/ !transision from ly/day to W
C
      if (nav.lt.1) nav=1
C    GLOBAL T,K.E.,N0
c      GMTAC1=GMTACC/NAV*2.
c      GMKEA1=GMKEAC/NAV*2.
c      GMRAC1=GMRACC/NAV*2.
C	lun=11
	OPEN (lun, FILE='przon.OUT',ACCESS='APPEND')
      WRITE (lun,9999) TAU,MONTH,MNTHDY,TOFDAY,GMKE,GMT,GMR,SDEDY,ID
 9999 FORMAT (1X,F9.2,2X,I3,1H/,I2,1H/,F6.2,3H.00,5X,5HK.E.=,1PE15.7,
     1   3H T=,1PE15.7,4H N0=,1PE15.7,'  SDEDY=',I3,' ID=',A4)
C
      DO 40 J=1,46
      DO 40 I=1,74
 40   ZONAV1(J,I)=ZONAVG(J,I)/(72.*NAV)          !72=IM
C
         DO 41 J=1,46
         ZONAV1(J,51)=ZONAV1(J,51)-TICE
         ZONAV1(J,50)=ZONAV1(J,50)-TICE
         ZONAV1(J,34)=ZONAV1(J,34)-TICE
 41      ZONAV1(J,37)=ZONAV1(J,37)-TICE
         DO 42 J=1,46
      AN0(J)=ZONAV1(J,8)+ZONAV1(J,9)+ZONAV1(J,18)-ZONAV1(J,12)
      ANS(J)=ZONAV1(J,18)-ZONAV1(J,14)
      BS(J)=ANS(J)-ZONAV1(J,7)-ZONAV1(J,54)*day-ZONAV1(J,55) !surface balance
      BA(J)=0.
      BEA(J)=0.
      ALFAP(J)=0.
      IF (ABS(ZONAV1(J,20)).LE.1.E-10) GO TO 42
      ALFAP(J)=ZONAV1(J,19)/ZONAV1(J,20) !planetary albedo
  42  CONTINUE
C      ENERGY FLUXES W/M**2
      DO 43 J=1,46
      DO 44 I=1,3
      ZONAV1(J,6+I)=ZONAV1(J,6+I)*CJW
      ZONAV1(J,11+I)=ZONAV1(J,11+I)*CJW
  44  ZONAV1(J,17+I)=ZONAV1(J,17+I)*CJW
      ZONAV1(J,54)=ZONAV1(J,54)*CJW*day !1 day=86400 sec
      ZONAV1(J,55)=ZONAV1(J,55)*CJW
      AN0(J)=AN0(J)*CJW
      ANS(J)=ANS(J)*CJW
  43  BS(J)=BS(J)*CJW
C
C        PRINT OUT ZONAL MEANS
C
      write (lun, 98)
  98  FORMAT(30X,'********  ZONAL  MEANS ***********')
      DO 2 N=1,45,11
      I1=N
      I2=N+10
      write (lun,94)(I,TITLE(I),I=I1,I2)
      do 2 j=1,46
      write (lun,95) J,(ZONAV1(J,I),I=I1,I2),J 
  2   continue
      write (lun,94)(I,TITLE(I),I=56,66)
      do 10 j=1,46
      write (lun,95) J,(ZONAV1(J,I),I=56,62),AN0(J),ALFAP(J),
     1         ANS(J),BS(J),J
  10  continue
      close (lun)
C
  95  FORMAT (1X,I2,11E11.4,1X,I2)
  94  FORMAT (3X,11(2X,I2,1X,A4,2X))
C
      RETURN
      END
