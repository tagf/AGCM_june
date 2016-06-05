C     *****************
C     *****************
c      SUBROUTINE  INACC
C     *****************
C     *****************
C
C              INPUT MONTHLY   ACCUMULATED VARIABLES HISTORY
C
C
C************ BEGINNING OF COMMON ************
C
      REAL KAPA,KAPEL,LAT
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
     6 JM,IM,ID,MNTHDY,SDEDY,NCYCLE,NC3,MONTH,MRCH,NSTEP,NAV
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
      LOGICAL SNOW,SEAICE
C
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)
C
C        ACCUMULATED VARIABLES IN COMMON
C
      COMMON /ACCUM/ APLS(72,46),APCV(72,46),ASNF(72,46),
     1               ACNO(72,46),ACNS(72,46),ACBS(72,46),
     2               ACBA(72,46),ABEA(72,46),ASCZ(72,46),
     3               ARET(72,46),ASLP(72,46),AT4(72,46),ZONAVG(46,74)
       dimension zero(72*12+74,46)
       equivalence (apls(1,1),zero(1,1))
C
	include 'dir.fi'
C************ END OF COMMON ******************
C
      dimension arr(72,46)
	 character*30 infile
      WorkDir='id##4x5b'
C
          LUN=43
	OPEN (lun,FILE='E:\CLIMATE\'//WorkDir//'\ACC',STATUS='OLD',
     1      FORM='UNFORMATTED')
         rewind lun
      READ (LUN) C
      READ (LUN)  ((ISFTYP(I,J),I=2,73),J=1,JM)
      READ (LUN)  APLS
      READ (LUN)  APCV
      READ (LUN)  ASNF
      READ (LUN)  ACNO
      READ (LUN)  ACNS
      READ (LUN)  ACBS
      READ (LUN)  ACBA
      READ (LUN)  ABEA
      READ (LUN)  ASCZ
      READ (LUN)  ARET
      READ (LUN)  ASLP
      READ (LUN)  AT4
      READ (LUN)  ZONAVG
C
 21   DO 30 N=1,46
      DO 30 I=1,72*12+74
 30     ZERO(I,N)=ZERO(I,N)/NAV
        close (lun)
        zer=1000.
          DO 25 I=1,72
          DO 25 J=1,46
          arr(i,j)=aslp(i,j)
 25        arr(i,j)=arr(i,j)-zer

 88    print '(a,$)' , ' Input SURFER file name (no extention): '
       read (*,'(a)') surfile
       open (12, file='e:\climate\id##4x5b\'//surfile//'.dat',
     1   status='new',err=88)

      DO 2 I=1,72
      DO 2 J=1,46
 2     write (12,11) i,j,arr(i,j)
 11    FORMAT (1x,2i3,1x,e12.4)

          infile='E:\CLIMATE\'//WorkDir//'\ACC'
        PRINT 111,infile,TAU
 111     FORMAT (' INPUT ACCUMULATED DATA. FILE ',a,' TAU=',F9.2)
      stop 'Normal end'
      END
