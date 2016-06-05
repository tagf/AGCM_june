C     ****************
C     ****************
      program isftyp
C     ****************
C     ****************
C
C    change types for gldstn 2014
C             READS INITIAL CONDITIONS
C
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))
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
C        MONTHLY SEA TEMPERATURES
C
      COMMON /SEA/ SS(74,46,12), SWIN(74,46,12)
C
	include 'dir.fi'

C***  FOR ICE MODELLING
      include 'paramz.fi'
       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
      include 'ice.fi'
c**************************
      integer t(72,46)
C
C************ END OF COMMON ******************
C
         REAL*4 SNGLSS(74,46,12)
	DIMENSION SWIN1(26,16,12)
	 character*40 restfile,restfile_ice,seatfile

      OPEN (lun,FILE=restfile,FORM='UNFORMATTED')
      rewind lun
      READ (lun) ((PHIS(I,J),I=2,73),J=1,46)
   !!!   READ (lun) C
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

      READ (lun) ((ISFTYP(I,J),I=2,73),J=1,46)
      READ (lun) ((P(I,J),I=2,73),J=1,46)
      READ (lun) ((U(I,J,1),I=2,73),J=1,46)
      READ (lun) ((U(I,J,2),I=2,73),J=1,46)
      READ (lun) ((V(I,J,1),I=2,73),J=1,46)
      READ (lun) ((V(I,J,2),I=2,73),J=1,46)
      READ (lun) ((T(I,J,1),I=2,73),J=1,46)
      READ (lun) ((T(I,J,2),I=2,73),J=1,46)
      READ (lun) ((QW(I,J,1),I=2,73),J=1,46)
      READ (lun) ((QW(I,J,2),I=2,73),J=1,46)
      READ (lun) ((GW(I,J),I=2,73),J=1,46)
      READ (lun) ((GT(I,J),I=2,73),J=1,46)
      READ (lun) ((SNOAMT(I,J),I=2,73),J=1,46)
      READ (lun) ((TS(I,J),I=2,73),J=1,46)
      READ (lun) ((SD(I,J),I=2,73),J=1,46)
	close (lun)
c      SNOAMT=0.  !first time
C***  FOR ICE MODELLING
      OPEN (lun2,FILE=restfile_ice,FORM='UNFORMATTED')
      rewind lun2
      READ (lun2) TSM
      READ (lun2) TMM
      READ (lun2) TBM
      READ (lun2) HFM
      READ (lun2) HSM
      READ (lun2) UGM
      READ (lun2) VGM
      READ (lun2) UDM
      READ (lun2) VDM
      READ (lun2) W1M
      READ (lun2) W2M
      READ (lun2) GHIce
      READ (lun2) GHIce2
      READ (lun2) GIceComp
      READ (lun2) GIcePool
      READ (lun2) GIce2Pool
      READ (lun2) GIceHeatng
      READ (lun2) GIce2Heatng
      READ (lun2) GIceCND
      READ (lun2) GIce2CND
      READ (lun2) GIceUpFlow
      READ (lun2) GT1
      READ (lun2) GT2
      READ (lun2) SnoAmt1
      READ (lun2) SnoAmt2
c      READ (lun2) NewIce
      t(1:72,:)=isftyp(2:73,:)
      t(6,4)=t(5,40)
      t(28,4)=t(31,4)
      t(29,4)=t(31,4)
      t(30,4)=t(31,4)
      t(18,5)=t(17,5)
      t(19,5)=t(17,5)
      t(20,5)=t(17,5)
      t(21,5)=t(17,5)
      t(25,5)=t(26,5)
      t(70,5)=t(71,5)
      t(23,10)=t(24,10)
      t(22,11)=t(23,11)
      t(22,12)=t(23,12)
      t(22,13)=t(23,13)
      t(22,13)=t(23,13)
      t(25,14)=t(26,14)
      t(66,14)=t(67,14)
      t(26,15)=t(27,15)
      t(26,15)=t(27,15)
      t(27,16)=t(28,16)
      t(27,16)=t(28,16)
      t(27,16)=t(28,16)
      t(43,17)=t(44,18)
      t(67,15)=t(68,15)
      t(67,16)=t(68,16)
      t(67,17)=t(68,17)
      t(23,18)=t(24,18)
      t(44,18)=t(45,18)
      t(46,18)=t(47,18)
      t(44,19)=t(45,19)
      t(29,19)=t(30,19)
      t(61,19)=t(62,19)
      t(29,20)=t(30,20)
      t(45,20)=t(46,20)
      t(62,20)=t(63,20)
      t(21,21)=t(22,21)
      t(29,21)=t(30,21)
      t(62,21)=t(63,21)
      t(64,21)=t(64,20)
      t(64,22)=t(64,20)
      t(66,22)=t(67,22)
      t(28,23)=t(29,23)
      t(45,23)=t(46,23)
      t(59,23)=t(60,23)
      t(64,23)=t(65,23)
      t(21,24)=t(22,24)
      t(45,24)=t(46,24)
      t(58,24)=t(57,24)
      t(59,24)=t(57,24)
      t(26,25)=t(25,25)
      t(57,25)=t(57,24)
      t(23,26)=t(22,26)
      t(24,26)=t(22,26)
      t(57,26)=t(57,27)
      t(20,27)=t(21,27)
      t(52,27)=t(51,27)
      t(53,27)=t(54,27)
      t(18,28)=t(17,28)
      t(19,28)=t(17,28)
      t(45,28)=t(46,28)
      t(47,28)=t(49,28)
      t(48,28)=t(49,28)
      t(53,28)=t(54,28)
      t(17,29)=t(18,29)
      t(48,29)=t(49,29)
      t(17,30)=t(18,30)
      t(34,30)=t(35,30)
      t(14,31)=t(15,31)
      t(18,31)=t(18,30)
      t(19,31)=t(19,30)
      t(20,31)=t(21,31)
      t(21,32)=t(22,32)
      t(44,32)=t(45,32)
      t(64,32)=t(63,32)
      t(36,33)=t(37,33)
      t(65,33)=t(66,33)
      t(12,34)=t(13,34)
      t(35,34)=t(36,34)
      t(63,34)=t(64,34)
      t(12,35)=t(13,35)
      t(24,35)=t(25,35)
      t(64,35)=t(65,35)
      t(25,36)=t(24,36)
      t(26,36)=t(27,36)
      t(37,36)=t(38,36)
      t(65,36)=t(66,36)
      t(11,37)=t(12,37)
      t(38,37)=t(39,37)
      t(65,37)=t(66,37)
      t(4,38)=t(5,38)
      t(36,38)=t(36,37)
      t(38,38)=t(39,38)
      t(66,38)=t(65,38)
      t(67,38)=t(65,38)
      t(68,38)=t(65,38)
      t(69,38)=t(70,38)
      t(3,39)=t(4,39)
      t(23,39)=t(22,39)
      t(28,39)=t(27,39)
      t(38,39)=t(39,39)
      t(70,39)=t(69,39)
      t(71,39)=t(69,39)
      t(2,40)=t(2,39)
      t(3,40)=t(4,40)
      t(23,40)=t(22,40)
      t(28,40)=t(27,40)
      t(38,40)=t(39,40)
      t(1,40)=t(72,40)
      t(3,41)=t(5,41)
      t(4,41)=t(5,41)
      t(6,41)=t(10,41)
      t(7,41)=t(10,41)
      t(8,41)=t(10,41)
      t(13,41)=t(14,41)
      t(23,41)=t(24,41)
      t(26,41)=t(27,41)
      t(32,41)=t(33,41)
      t(40,41)=t(41,41)
      t(41,41)=t(42,41)
      t(44,41)=t(43,41)
      t(45,41)=t(43,41)
      t(46,41)=t(43,41)
      t(47,41)=t(51,41)
      t(48,41)=t(51,41)
      t(49,41)=t(51,41)
      t(50,41)=t(51,41)
      t(14,42)=t(13,42)
      t(25,42)=t(26,42)
      t(33,42)=t(34,42)
      t(54,42)=t(53,42)
      t(58,42)=t(61,42)
      t(59,42)=t(61,42)
      t(60,42)=t(61,42)
      t(22,43)=t(21,43)
      t(23,43)=t(24,43)
      t(40,43)=t(39,43)
      t(41,43)=t(42,43)
      t(52,43)=t(51,43)
      t(53,43)=t(54,43)
      t(20,44)=t(19,44)
      t(21,44)=t(19,44)
      t(22,44)=t(19,44)
      t(23,44)=t(19,44)
      t(32,44)=t(35,44)
      t(33,44)=t(35,44)
      t(34,44)=t(35,44)
      isftyp(2:73,:)=t(1:72,:)
	close (lun2)
c**************************
C
      stop 'normal end'
      END
