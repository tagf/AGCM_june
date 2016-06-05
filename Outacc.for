C     *****************
C     *****************
      SUBROUTINE OUTACC(lun,nappend)
C     *****************
C     *****************
C
C
C             OUTPUTS ACCUMULATED VARIABLES HISTORY
C                      
c         nappend=0   rewind
c         nappend=1   append
C			write *4 data - disk resourses economy
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)

C        ACCUMULATED VARIABLES IN COMMON
C
	include 'accum.fi'

	 include 'dir.fi'
C
C************ END OF COMMON ******************
C
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))
C
      IF (MONTH.LT.1.OR.MONTH.GT.12) WRITE (6,9) MONTH,TAU
    9 FORMAT (14H0ERR- MONTH = ,I8,7H TAU = ,F10.2)
C
      if (nappend.eq.0) then 
	 OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\ACC',FORM='UNFORMATTED')
	 REWIND LUN
	endif
      if (nappend.eq.1) then
	 OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\ACC1',  
     1      ACCESS='append',FORM='UNFORMATTED',status='unknown'
     2      ,err=1)
	 goto 2
1       write(*,'(a)') ' open acc1 error - ok '
c first open acc1-usual open ACCESS='append' does't work, if acc1 file 
c does't exist or is empty
	  OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\ACC1',  
     1      ACCESS='SEQUENTIAL',FORM='UNFORMATTED',status='unknown')
	endif
    !!!  write (LUN) C
 2     write (lun)   
     1 TAU   ,TAUI  ,TAUO  ,TAUD  ,TAUE  ,TAUH  ,TAUC  ,TOFDAY,ROT   , 
     1 DT    ,DLAT  ,DLON  ,RAD   ,RSDIST,SIND  ,COSD  ,COSR  ,SINR  , 
     1 DAYPYR,ROTPER,SDEYR ,SOLTCE,APHEL ,DECMAX,ECCN  ,GUSTY ,        
     1 DAY   ,GRAV  ,RGAS  ,KAPA  ,PSF   ,PTROP ,PSL   ,CFG   , 
     1 FM    ,ED    ,PI    ,SIG1  ,SIG3  ,DSIG  ,PM    ,KAPEL ,RKAPA1, 
     1 STBO  ,GWM   ,DTC3  ,FRDAY ,CLH   ,COE1  ,HICE  ,CTYICE,CNDICE, 
     1 TICE  ,TCICE ,SNOWL ,COE   ,TSPD  ,PSTQ  ,QST   ,TST   ,PTRK  , 
     1PSFHO ,CALFA ,QC    ,PC    ,QCONST,EFVCON,TCT0  ,TCT2  ,TCT4  , 
     1TCST0 ,TCST2 ,TC01  ,TC02  ,TC03  ,TC04  ,TC23  ,TC24  ,TC34  , 
     1BLC   ,ALOGP0,FIM   ,HRGAS ,TCST4 ,FLR   ,PS4K  ,PS8K  ,ELOG  , 
     1S0, TOZONE    ,LAT   ,DXU   ,DXP   ,DYU,DYP, 
     1SFCALB,SINL  ,COSL  ,O3AMT ,TC12  ,TC14  , 
     1DXYP  ,F     ,SIG   ,FLEADN    ,FLEADS    ,ERROR  , 
     1COSLN,SINLN,DXYU1,DXYU2, GMT,GMR,GMKE, 
     1 JM,IM,ID,MNTHDY,SDEDY,NCYCLE,NC3,MONTH,MRCH,NSTEP,NAV

      print *,JM,IM,ID,MNTHDY,SDEDY,NCYCLE,NC3,MONTH,MRCH,NSTEP,NAV
      write (LUN)  ((ISFTYP(I,J),I=2,73),J=1,JM)
!      pause 'outacc'
      WRITE (LUN)   APLS
      WRITE (LUN)   APCV
      WRITE (LUN)   ASNF
      WRITE (LUN)   ACNO
      WRITE (LUN)   ACNS
      WRITE (LUN)   ACBS
      WRITE (LUN)   ACBA
      WRITE (LUN)   ABEA
      WRITE (LUN)   ASCZ
      WRITE (LUN)   ARET
      WRITE (LUN)   ASLP
      WRITE (LUN)   AT4
      WRITE (LUN)   ZONAVG
C
          WRITE (LUN) (( PHIS(I,J),I=2,73),J=1,46)
          WRITE (LUN)  PAC  
          WRITE (LUN)  UAC  
          WRITE (LUN)  VAC  
          WRITE (LUN)  TAC  
          WRITE (LUN)  QWAC  
          WRITE (LUN)  GWAC  
          WRITE (LUN)  GTAC  
          WRITE (LUN)  SNOAC  
          WRITE (LUN)  TSAC  
          WRITE (LUN)  SDAC  
          WRITE (LUN)  TAUUA  
          WRITE (LUN)  TAUVA  
          WRITE (LUN)  APRES
          WRITE (LUN)  ACLOUD
         CLOSE (LUN)
c         write (7,111) TAU,MONTH,MNTHDY,TOFDAY ! To PROTOCOL FILE
         print 111,tau,MONTH,MNTHDY,TOFDAY
 111     FORMAT (' ****  OUTACC RECORDS    TAU=',F9.2,
     1           1X,I3,1H/,I2,1H/,F6.2)
      
      RETURN
      END
