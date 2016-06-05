C     *****************
C     *****************
      SUBROUTINE OUTAC1(lun)
C     *****************
C     *****************
C
C
C             OUTPUTS ACCUMULATED VARIABLES HISTORY
C                      LUN=42
C		  write *4 data - to save disk space
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
     6 JM,IM,ID,MNTHDY,SDEDY ,NCYCLE,NC3   ,MONTH ,MRCH  ,NSTEP,nav
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1  ,SD(74,46)
C
C        ACCUMULATED VARIABLES IN COMMON
C
	include 'accum.fi'
C
C************ END OF COMMON ******************
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
      IF (MONTH.LT.1.OR.MONTH.GT.12) WRITE (6,9) MONTH,TAU
    9 FORMAT (14H0ERR- MONTH = ,I8,7H TAU = ,F10.2)
C
C        OUTPUT UNITS ARE 40 THRU 52
c        LUN=42
      write (LUN) C
      write (LUN)  ((ISFTYP(I,J),I=2,73),J=1,46)
      WRITE (LUN)  sngl(APLS)
      WRITE (LUN)  sngl(APCV)
      WRITE (LUN)  sngl(ASNF)
      WRITE (LUN)  sngl(ACNO)
      WRITE (LUN)  sngl(ACNS)
      WRITE (LUN)  sngl(ACBS)
      WRITE (LUN)  sngl(ACBA)
      WRITE (LUN)  sngl(ABEA)
      WRITE (LUN)  sngl(ASCZ)
      WRITE (LUN)  sngl(ARET)
      WRITE (LUN)  sngl(ASLP)
      WRITE (LUN)  sngl(AT4)
      WRITE (LUN)  sngl(ZONAVG)
C
          WRITE (LUN) ((sngl(PHIS(I,J)),I=2,73),J=1,46)
          WRITE (LUN) sngl(PACC)
          WRITE (LUN) sngl(UACC)
          WRITE (LUN) sngl(VACC)
          WRITE (LUN) sngl(TACC)
          WRITE (LUN) sngl(QWACC)
          WRITE (LUN) sngl(GWACC)
          WRITE (LUN) sngl(GTACC)
          WRITE (LUN) sngl(SNOACC)
          WRITE (LUN) sngl(TSACC)
          WRITE (LUN) sngl(SDACC)
          WRITE (LUN) sngl(TAUUAC)
          WRITE (LUN) sngl(TAUVAC)
          WRITE (LUN) sngl(APRES)
          WRITE (LUN) sngl(ACLOUD)
         write (7,111) TAU,MONTH,MNTHDY,TOFDAY ! To PROTOCOL FILE
         print 111,tau,MONTH,MNTHDY,TOFDAY
 111     FORMAT (' ****  OUTAC1 RECORDS    TAU=',F9.2,
     1           1X,I3,1H/,I2,1H/,F6.2)
      RETURN
      END
