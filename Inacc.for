C     *****************
C     *****************
      SUBROUTINE  INACC(lun)
C     *****************
C     *****************
C
C              INPUT MONTHLY   ACCUMULATED VARIABLES HISTORY
C			 read *4 data and convert to *8 data
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))
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
	include 'dir.fi'
C************ END OF COMMON ******************
C
	 character*30 infile
	real*4 tt(72,46),tt1(46,74),tt2(72,46,2) ! temporary *4 arrays
       equivalence (tt(1,1),tt1(1,1),tt2(1,1,1))
      integer ttt(72,46)    !ISFTYP correction
      real rrr(72,46)
C
	 OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\ACC',FORM='UNFORMATTED')
      rewind lun
   !!!   READ (LUN) C
      READ (lun)   
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

      READ (LUN)  ((ISFTYP(I,J),I=2,73),J=1,JM)
      READ (LUN)   APLS
      READ (LUN)   APCV
      READ (LUN)   ASNF
      READ (LUN)   ACNO
      READ (LUN)   ACNS
      READ (LUN)   ACBS
      READ (LUN)   ACBA
      READ (LUN)   ABEA
      READ (LUN)   ASCZ
      READ (LUN)   ARET
      READ (LUN)   ASLP
      READ (LUN)   AT4
      READ (LUN)   ZONAVG
C
          READ (LUN) (( PHIS(I,J),I=2,73),J=1,46)
          READ (LUN)  PAC  
          READ (LUN)  UAC  
          READ (LUN)  VAC  
          READ (LUN)  TAC  
          READ (LUN)  QWAC  
          READ (LUN)  GWAC  
          READ (LUN)  GTAC  
          READ (LUN)  SNOAC  
          READ (LUN)  TSAC  
          READ (LUN)  SDAC  
          READ (LUN)  TAUUA  
          READ (LUN)  TAUVA  
          READ (LUN)  APRES
          READ (LUN)  ACLOUD
	close (LUN)
c**************************
!Приближение типа точек isftyp сетки AGCM к типу точек GLDSTN
      goto 40 ! нет измен
      ttt(1:72,:)=isftyp(2:73,:)
      ttt(6,4)=ttt(5,40)
      ttt(28,4)=ttt(31,4)
      ttt(29,4)=ttt(31,4)
      ttt(30,4)=ttt(31,4)
      ttt(18,5)=ttt(17,5)
      ttt(19,5)=ttt(17,5)
      ttt(20,5)=ttt(17,5)
      ttt(21,5)=ttt(17,5)
      ttt(25,5)=ttt(26,5)
      ttt(70,5)=ttt(71,5)
      ttt(23,10)=ttt(24,10)
      ttt(22,11)=ttt(23,11)
      ttt(22,12)=ttt(23,12)
      ttt(22,13)=ttt(23,13)
      ttt(22,13)=ttt(23,13)
      ttt(25,14)=ttt(26,14)
      ttt(66,14)=ttt(67,14)
      ttt(26,15)=ttt(27,15)
      ttt(26,15)=ttt(27,15)
      ttt(27,16)=ttt(28,16)
      ttt(27,16)=ttt(28,16)
      ttt(27,16)=ttt(28,16)
      ttt(43,17)=ttt(44,18)
      ttt(67,15)=ttt(68,15)
      ttt(67,16)=ttt(68,16)
      ttt(67,17)=ttt(68,17)
      ttt(23,18)=ttt(24,18)
      ttt(44,18)=ttt(45,18)
      ttt(46,18)=ttt(47,18)
      ttt(44,19)=ttt(45,19)
      ttt(29,19)=ttt(30,19)
      ttt(61,19)=ttt(62,19)
      ttt(29,20)=ttt(30,20)
      ttt(45,20)=ttt(46,20)
      ttt(62,20)=ttt(63,20)
      ttt(21,21)=ttt(22,21)
      ttt(29,21)=ttt(30,21)
      ttt(62,21)=ttt(63,21)
      ttt(64,21)=ttt(64,20)
      ttt(64,22)=ttt(64,20)
      ttt(66,22)=ttt(67,22)
      ttt(28,23)=ttt(29,23)
      ttt(45,23)=ttt(46,23)
      ttt(59,23)=ttt(60,23)
      ttt(64,23)=ttt(65,23)
      ttt(21,24)=ttt(22,24)
      ttt(45,24)=ttt(46,24)
      ttt(58,24)=ttt(57,24)
      ttt(59,24)=ttt(57,24)
      ttt(26,25)=ttt(25,25)
      ttt(57,25)=ttt(57,24)
      ttt(23,26)=ttt(25,26)
      ttt(24,26)=ttt(25,26)
      ttt(57,26)=ttt(57,27)
      ttt(20,27)=ttt(21,27)
      ttt(52,27)=ttt(51,27)
      ttt(53,27)=ttt(54,27)
      ttt(18,28)=ttt(17,28)
      ttt(19,28)=ttt(17,28)
      ttt(45,28)=ttt(46,28)
      ttt(47,28)=ttt(49,28)
      ttt(48,28)=ttt(49,28)
      ttt(53,28)=ttt(54,28)
      ttt(17,29)=ttt(18,29)
      ttt(48,29)=ttt(49,29)
      ttt(17,30)=ttt(18,30)
      ttt(34,30)=ttt(35,30)
      ttt(14,31)=ttt(15,31)
      ttt(18,31)=ttt(18,30)
      ttt(19,31)=ttt(19,30)
      ttt(20,31)=ttt(21,31)
      ttt(21,32)=ttt(22,32)
      ttt(44,32)=ttt(45,32)
      ttt(64,32)=ttt(63,32)
      ttt(36,33)=ttt(37,33)
      ttt(65,33)=ttt(66,33)
      ttt(12,34)=ttt(13,34)
      ttt(35,34)=ttt(36,34)
      ttt(63,34)=ttt(64,34)
      ttt(12,35)=ttt(13,35)
      ttt(24,35)=ttt(25,35)
      ttt(64,35)=ttt(65,35)
      ttt(25,36)=ttt(24,36)
      ttt(26,36)=ttt(27,36)
      ttt(37,36)=ttt(38,36)
      ttt(65,36)=ttt(66,36)
      ttt(11,37)=ttt(12,37)
      ttt(38,37)=ttt(39,37)
      ttt(65,37)=ttt(66,37)
      ttt(4,38)=ttt(5,38)
      ttt(36,38)=ttt(36,37)
      ttt(38,38)=ttt(39,38)
      ttt(66,38)=ttt(65,38)
      ttt(67,38)=ttt(65,38)
      ttt(68,38)=ttt(65,38)
      ttt(69,38)=ttt(70,38)
      ttt(3,39)=ttt(4,39)
      ttt(23,39)=ttt(22,39)
      ttt(28,39)=ttt(27,39)
      ttt(38,39)=ttt(39,39)
      ttt(70,39)=ttt(72,39)
      ttt(71,39)=ttt(72,39)
      ttt(2,40)=ttt(2,39)
      ttt(3,40)=ttt(4,40)
      ttt(23,40)=ttt(22,40)
      ttt(28,40)=ttt(27,40)
      ttt(38,40)=ttt(39,40)
      ttt(1,40)=ttt(72,40)
      ttt(3,41)=ttt(5,41)
      ttt(4,41)=ttt(5,41)
      ttt(6,41)=ttt(10,41)
      ttt(7,41)=ttt(10,41)
      ttt(8,41)=ttt(10,41)
      ttt(13,41)=ttt(14,41)
      ttt(23,41)=ttt(24,41)
      ttt(26,41)=ttt(27,41)
      ttt(32,41)=ttt(33,41)
      ttt(40,41)=ttt(41,41)
      ttt(41,41)=ttt(42,41)
      ttt(44,41)=ttt(43,41)
      ttt(45,41)=ttt(43,41)
      ttt(46,41)=ttt(43,41)
      ttt(47,41)=ttt(51,41)
      ttt(48,41)=ttt(51,41)
      ttt(49,41)=ttt(51,41)
      ttt(50,41)=ttt(51,41)
      ttt(14,42)=ttt(13,42)
      ttt(25,42)=ttt(26,42)
      ttt(33,42)=ttt(34,42)
      ttt(54,42)=ttt(53,42)
      ttt(58,42)=ttt(61,42)
      ttt(59,42)=ttt(61,42)
      ttt(60,42)=ttt(61,42)
      ttt(22,43)=ttt(21,43)
      ttt(23,43)=ttt(24,43)
      ttt(40,43)=ttt(39,43)
      ttt(41,43)=ttt(42,43)
      ttt(52,43)=ttt(51,43)
      ttt(53,43)=ttt(54,43)
      ttt(20,44)=ttt(19,44)
      ttt(21,44)=ttt(19,44)
      ttt(22,44)=ttt(19,44)
      ttt(23,44)=ttt(19,44)
      ttt(32,44)=ttt(35,44)
      ttt(33,44)=ttt(35,44)
      ttt(34,44)=ttt(35,44)
      isftyp(2:73,:)=ttt(1:72,:)
C
      rrr(1:72,:)=PHIS(2:73,:)
      rrr(6,4)=rrr(5,4)
      rrr(28,4)=rrr(31,4)
      rrr(29,4)=rrr(31,4)
      rrr(30,4)=rrr(31,4)
      rrr(18,5)=rrr(17,5)
      rrr(19,5)=rrr(17,5)
      rrr(20,5)=rrr(17,5)
      rrr(21,5)=rrr(17,5)
      rrr(25,5)=rrr(26,5)
      rrr(70,5)=rrr(71,5)
      rrr(23,10)=rrr(24,10)
      rrr(22,11)=rrr(23,11)
      rrr(22,12)=rrr(23,12)
      rrr(22,13)=rrr(23,13)
      rrr(22,13)=rrr(23,13)
      rrr(25,14)=rrr(26,14)
      rrr(66,14)=rrr(67,14)
      rrr(26,15)=rrr(27,15)
      rrr(26,15)=rrr(27,15)
      rrr(27,16)=rrr(28,16)
      rrr(27,16)=rrr(28,16)
      rrr(27,16)=rrr(28,16)
      rrr(43,17)=rrr(44,18)
      rrr(67,15)=rrr(68,15)
      rrr(67,16)=rrr(68,16)
      rrr(67,17)=rrr(68,17)
      rrr(23,18)=rrr(24,18)
      rrr(44,18)=rrr(45,18)
      rrr(46,18)=rrr(47,18)
      rrr(44,19)=rrr(45,19)
      rrr(29,19)=rrr(30,19)
      rrr(61,19)=rrr(62,19)
      rrr(29,20)=rrr(30,20)
      rrr(45,20)=rrr(46,20)
      rrr(62,20)=rrr(63,20)
      rrr(21,21)=rrr(22,21)
      rrr(29,21)=rrr(30,21)
      rrr(62,21)=rrr(63,21)
      rrr(64,21)=rrr(64,20)
      rrr(64,22)=rrr(64,20)
      rrr(66,22)=rrr(67,22)
      rrr(28,23)=rrr(29,23)
      rrr(45,23)=rrr(46,23)
      rrr(59,23)=rrr(60,23)
      rrr(64,23)=rrr(65,23)
      rrr(21,24)=rrr(22,24)
      rrr(45,24)=rrr(46,24)
      rrr(58,24)=rrr(57,24)
      rrr(59,24)=rrr(57,24)
      rrr(26,25)=rrr(25,25)
      rrr(57,25)=rrr(57,24)
      rrr(23,26)=rrr(25,26)
      rrr(24,26)=rrr(25,26)
      rrr(57,26)=rrr(57,27)
      rrr(20,27)=rrr(21,27)
      rrr(52,27)=rrr(51,27)
      rrr(53,27)=rrr(54,27)
      rrr(18,28)=rrr(17,28)
      rrr(19,28)=rrr(17,28)
      rrr(45,28)=rrr(46,28)
      rrr(47,28)=rrr(49,28)
      rrr(48,28)=rrr(49,28)
      rrr(53,28)=rrr(54,28)
      rrr(17,29)=rrr(18,29)
      rrr(48,29)=rrr(49,29)
      rrr(17,30)=rrr(18,30)
      rrr(34,30)=rrr(35,30)
      rrr(14,31)=rrr(15,31)
      rrr(18,31)=rrr(18,30)
      rrr(19,31)=rrr(19,30)
      rrr(20,31)=rrr(21,31)
      rrr(21,32)=rrr(22,32)
      rrr(44,32)=rrr(45,32)
      rrr(64,32)=rrr(63,32)
      rrr(36,33)=rrr(37,33)
      rrr(65,33)=rrr(66,33)
      rrr(12,34)=rrr(13,34)
      rrr(35,34)=rrr(36,34)
      rrr(63,34)=rrr(64,34)
      rrr(12,35)=rrr(13,35)
      rrr(24,35)=rrr(25,35)
      rrr(64,35)=rrr(65,35)
      rrr(25,36)=rrr(24,36)
      rrr(26,36)=rrr(27,36)
      rrr(37,36)=rrr(38,36)
      rrr(65,36)=rrr(66,36)
      rrr(11,37)=rrr(12,37)
      rrr(38,37)=rrr(39,37)
      rrr(65,37)=rrr(66,37)
      rrr(4,38)=rrr(5,38)
      rrr(36,38)=rrr(36,37)
      rrr(38,38)=rrr(39,38)
      rrr(66,38)=rrr(65,38)
      rrr(67,38)=rrr(65,38)
      rrr(68,38)=rrr(65,38)
      rrr(69,38)=rrr(70,38)
      rrr(3,39)=rrr(4,39)
      rrr(23,39)=rrr(22,39)
      rrr(28,39)=rrr(27,39)
      rrr(38,39)=rrr(39,39)
      rrr(70,39)=rrr(72,39)
      rrr(71,39)=rrr(72,39)
      rrr(2,40)=rrr(2,39)
      rrr(3,40)=rrr(4,40)
      rrr(23,40)=rrr(22,40)
      rrr(28,40)=rrr(27,40)
      rrr(38,40)=rrr(39,40)
      rrr(1,40)=rrr(72,40)
      rrr(3,41)=rrr(5,41)
      rrr(4,41)=rrr(5,41)
      rrr(6,41)=rrr(10,41)
      rrr(7,41)=rrr(10,41)
      rrr(8,41)=rrr(10,41)
      rrr(13,41)=rrr(14,41)
      rrr(23,41)=rrr(24,41)
      rrr(26,41)=rrr(27,41)
      rrr(32,41)=rrr(33,41)
      rrr(40,41)=rrr(41,41)
      rrr(41,41)=rrr(42,41)
      rrr(44,41)=rrr(43,41)
      rrr(45,41)=rrr(43,41)
      rrr(46,41)=rrr(43,41)
      rrr(47,41)=rrr(51,41)
      rrr(48,41)=rrr(51,41)
      rrr(49,41)=rrr(51,41)
      rrr(50,41)=rrr(51,41)
      rrr(14,42)=rrr(13,42)
      rrr(25,42)=rrr(26,42)
      rrr(33,42)=rrr(34,42)
      rrr(54,42)=rrr(53,42)
      rrr(58,42)=rrr(61,42)
      rrr(59,42)=rrr(61,42)
      rrr(60,42)=rrr(61,42)
      rrr(22,43)=rrr(21,43)
      rrr(23,43)=rrr(24,43)
      rrr(40,43)=rrr(39,43)
      rrr(41,43)=rrr(42,43)
      rrr(52,43)=rrr(51,43)
      rrr(53,43)=rrr(54,43)
      rrr(20,44)=rrr(19,44)
      rrr(21,44)=rrr(19,44)
      rrr(22,44)=rrr(19,44)
      rrr(23,44)=rrr(19,44)
      rrr(32,44)=rrr(35,44)
      rrr(33,44)=rrr(35,44)
      rrr(34,44)=rrr(35,44)
      PHIS(2:73,:)=rrr(1:72,:)
  40  continue
C
      DO  J=1,JM
        PHIS(1,J)=PHIS(73,J)
        PHIS(74,J)=PHIS(2,J)
        ISFTYP(1,J)=ISFTYP(73,J)
        ISFTYP(74,J)=ISFTYP(2,J)
      enddo
      
        infile=TRIM(BaseDir)//WorkDir//'\ACC'
        PRINT 111,TRIM(infile),TAU
 111     FORMAT (' INPUT ACCUM DATA. FILE ',a,' TAU=',F9.2)
      RETURN
      END
