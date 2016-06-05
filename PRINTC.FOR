      SUBROUTINE  PRINTC
C
C
C*********    BEGINNING OF COMMON
C
C
         IMPLICIT REAL*8(A-H,O-Z),INTEGER (I-N)
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
     6 JM,IM,ID,MNTHDY,SDEDY ,NCYCLE,NC3 ,MONTH ,MRCH,NSTEP,CDUM(81)
C
c      DIMENSION TITLE(9)
C
c      OPEN (6,FILE='f:\CLIMATE\AGCM_Code1\carr',	 !Default Names
c     1      FORM='FORMATTED')
      WRITE( 6,100 )
  100 FORMAT( 2H1  //10X,'***** PRINT C - ARRAY *****')
      WRITE (6,110)
 110  FORMAT (1H0,//10X,'*****  GRID DESCRIPTION')
      TEMP1 = DLON*180.0/PI
      TEMP2 = DLAT*180.0/PI
      WRITE( 6,200 ) JM,DLAT,TEMP2,IM,DLON,TEMP1
  200 FORMAT(1H0,'JM    =',I6,11X,'DLAT  = (',1PE15.7,' RAD, ',1PE15.7,
     1 ' DEG)'//' IM    =',I6,11X,'DLON  = (',1PE15.7,' RAD, ',1PE15.7,
     2 ' DEG)')
      WRITE (6,300) (J,DXU(J),DXP(J),DYU(J),DYP(J),DXYP(J),DXYU1(J),
     1  DXYU2(J),J=1,46)
 300  FORMAT (1H0/9X,' DXU',12X,'DXP',12X,'DYU',12X,'DYP',12X,'DXYP',
     1  11X,'DXYU1',10X,'DXYU2'/(1X,I3,1P7E15.7) )
      WRITE (6,350) (J,LAT(J),SINL(J),COSL(J),F(J),J=1,JM)
  350 FORMAT (1H1,/9X,' LAT',12X,'SINL',11X,'COSL',11X,'F'/
     1 (1X,I3,1P4E15.7))
      WRITE (6,360)
 360  FORMAT (1H1,/9X,' LON',5X,'COSLN',10X,'SINLN',20X,'LON',12X,
     1 'COSLN',10X,'SINLN')
      DO 365 I=1,36
      LN1=-180+(I-1)*5
      LN2=LN1+180
 365  WRITE (6,370) LN1,COSLN(I),SINLN(I),LN2,COSLN(I+36),SINLN(I+36)
 370  FORMAT (9X,I4,2E15.7,15X,I4,2E15.7)
      WRITE (6,400)
  400 FORMAT (1H1,//'*****  VERTICAL DEFINITIONS')
      WRITE (6,500) (SIG(I),I=1,2),PTROP,PSF,PSL
  500 FORMAT (1H0,' SIG   =',1P2E15.7/' PTROP =',1PE15.7,' PSF   =',
     1    E15.7,' PSL   =',  E15.7)
      WRITE (6,600)
  600 FORMAT (1H0,//'*****  TIME INTEGRATION DEFINTIIONS')
      WRITE (6,700) DT,DTC3,NCYCLE,NC3
  700 FORMAT (1H0,' DT    =',1PE15.7,' DTC3  =',1PE15.7,' NCYCLE=',I6,
     1 9X,'NC3   =',I6)
C
C     NOW WRITE THE TIME DEPENDENT VARIABLES
C
      WRITE (6,800)
  800 FORMAT (1H ,//'*****  TIME DEPENDENT VARIABLES')
      WRITE (6,900) TAU,TAUI,TAUE,NSTEP,TAUC,TAUD,TAUH,RSDIST,
     1              SIND,COSD,SDEDY,MNTHDY,MONTH
  900 FORMAT (1H0,' TAU   =',1PE15.7,' TAUI  =',1PE15.7,' TAUE  =',
     1        1PE15.7,' NSTEP =',I6/'  TAUC  =',1PE15.7,' TAUD  =',
     2        1PE15.7,' TAUH  =',1PE15.7/'  RSDIST=',1PE15.7,' SIND  ='
     3       ,1PE15.7,' COSD  =',1PE15.7,' SDEDY =',I6/
     4        '  MNTHDY',I6,9X,' MONTH ',I6)
C
C     PHYSICAL CONSTANTS SECTION
C
      WRITE( 6,1200 )
 1200 FORMAT( 1H0,//10X,'****  PHYSICAL CONSTANTS  ****')
      WRITE (6,1000) RAD,DAYPYR,ROTPER,SOLTCE,APHEL,DECMAX,
     1               ECCN,DAY
 1000 FORMAT (1H0,/'  RAD   =',1PE15.7,' DAYPYR=',1PE15.7,' ROTPER=',
     1    E15.7,' SOLTCE=',  E15.7/'  APHEL =',  E15.7,' DECMAX=',
     2    E15.7,' ECCN  =',  E15.7,' DAY   =',  E15.7)
      WRITE (6,1100) GRAV,RGAS,KAPA,A,FM,ED,PI,SIG1,SIG3,DSIG,
     1  RKAPA1,STBO,GWM,FRDAY,CLH,GUSTY,HICE,CTYICE,CNDICE,TICE,
     2  TCICE,SNOWL,COE,COE1,TSPD,PSTQ,QST,TST,PSFHO,CALFA,QC,
     3  PC,QCONST,EFVCON,TCT0,TCT2,TCT4,TCST0,TCST2,TC01
      WRITE (6,1500) TC02,TC03,
     4  TC04,TC23,TC24,TC34,TCST4,BLC,ALOGP0,HRGAS,FIM,FLEADN,FLEADS
     2  ,FLR,PS4K,PS8K,PTRK,ELOG,KAPEL,CFG
 1100 FORMAT (1H0,/'  GRAV  =',1PE15.7,' RGAS  =',1PE15.7,' KAPA  =',
     1    E15.7,' A     =',  E15.7/'  FM    =',  E15.7,' ED    =',
     2    E15.7,' PI    =',  E15.7,' SIG1  =',  E15.7/'  SIG3  =',
     3    E15.7,' DSIG  =',  E15.7,' RKAPA1=',  E15.7,' STBO  =',
     4    E15.7/'  GWM   =',  E15.7,' FRDAY =',  E15.7,' CLH   =',
     5    E15.7,' GUSTY =',  E15.7/'  HICE  =',  E15.7,' CTYICE=',
     6    E15.7,' CNDICE=',  E15.7,' TICE  =',  E15.7/'  TCICE =',
     7    E15.7,' SNOWL =',  E15.7,' COE   =',  E15.7,' COE1  =',
     8    E15.7/'  TSPD  =',  E15.7,' PSTQ  =',  E15.7,' QST   =',
     9    E15.7,' TST   =',  E15.7/'  PSFHO =',  E15.7,' CALFA =',
     1    E15.7,' QC    =',  E15.7,' PC    =',  E15.7/'  QCONST=',
     2    E15.7,' EFVCON=',  E15.7,' TCT0  =',  E15.7,' TCT2  =',
     3    E15.7/'  TCT4  =',  E15.7,' TCST0 =',  E15.7,' TCST2 =',
     4    E15.7,' TC01  =',  E15.7)
 1500 FORMAT ('  TC02  =',  E15.7,' TC03  =',
     5    E15.7,' TC04  =',  E15.7,' TC23  =',  E15.7/'  TC24  =',
     5    E15.7,' TC34  =',  E15.7,' TCST4 =',  E15.7,' BLC   =',
     7    E15.7/'  ALOGP0=',  E15.7,' HRGAS =',  E15.7,' FIM   =',
     8    E15.7,' FLEADN=',  E15.7/'  FLEADS=',  E15.7,' FLR   =',
     9    E15.7,' PS4K  =',E15.7,' PS8K  =',E15.7/'  PTRK  =',E15.7,
     9    ' ELOG  =',E15.7,' KAPEL =',E15.7,' CFG   =',E15.7)
      WRITE (6,1600) SFCALB
 1600 FORMAT (16H0SFCALB W/O SNOW,9F7.2/16H SFCALB  W  SNOW,9F7.2)
      RETURN
      END
