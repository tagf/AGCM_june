C     *****************
C     *****************
      SUBROUTINE OUTRES(restfile,lun,restfile_ice,lun2)
C     *****************
C     *****************
C
C
C             OUTPUTS RESTART CONDITIONS
C
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
     1  ,SD(74,46)
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C***  FOR ICE MODELLING
      include 'paramz.fi'
       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
      include 'ice.fi'
	include 'dir.fi'
c**************************

C************ END OF COMMON ******************
	 character*40 restfile,restfile_ice
      OPEN (lun,FILE=restfile,FORM='UNFORMATTED')
      REWIND lun
C
      WRITE (lun) ((PHIS(I,J),I=2,73),J=1,46)
  !!!    WRITE (lun) C
      write (lun)   
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

      WRITE (lun) ((ISFTYP(I,J),I=2,73),J=1,46)
      WRITE (lun) ((P(I,J),I=2,73),J=1,46)
      WRITE (lun) ((U(I,J,1),I=2,73),J=1,46)
      WRITE (lun) ((U(I,J,2),I=2,73),J=1,46)
      WRITE (lun) ((V(I,J,1),I=2,73),J=1,46)
      WRITE (lun) ((V(I,J,2),I=2,73),J=1,46)
      WRITE (lun) ((T(I,J,1),I=2,73),J=1,46)
      WRITE (lun) ((T(I,J,2),I=2,73),J=1,46)
      WRITE (lun) ((QW(I,J,1),I=2,73),J=1,46)
      WRITE (lun) ((QW(I,J,2),I=2,73),J=1,46)
      WRITE (lun) ((GW(I,J),I=2,73),J=1,46)
      WRITE (lun) ((GT(I,J),I=2,73),J=1,46)
      WRITE (lun) ((SNOAMT(I,J),I=2,73),J=1,46)
      WRITE (lun) ((TS(I,J),I=2,73),J=1,46)
      WRITE (lun) ((SD(I,J),I=2,73),J=1,46)
	close (lun)
C***  FOR ICE MODELLING
      OPEN (lun2,FILE=restfile_ice,FORM='UNFORMATTED')
      REWIND lun2
      WRITE (lun2) TSM
      WRITE (lun2) TMM
      WRITE (lun2) TBM
      WRITE (lun2) HFM
      WRITE (lun2) HSM
      WRITE (lun2) UGM
      WRITE (lun2) VGM
      WRITE (lun2) UDM
      WRITE (lun2) VDM
      WRITE (lun2) W1M
      WRITE (lun2) W2M
      WRITE (lun2) GHIce
      WRITE (lun2) GHIce2
      WRITE (lun2) GIceComp
      WRITE (lun2) GIcePool
      WRITE (lun2) GIce2Pool
      WRITE (lun2) GIceHeatng
      WRITE (lun2) GIce2Heatng
      WRITE (lun2) GIceCND
      WRITE (lun2) GIce2CND
      WRITE (lun2) GIceUpFlow
      WRITE (lun2) GT1
      WRITE (lun2) GT2
      WRITE (lun2) SnoAmt1
      WRITE (lun2) SnoAmt2
c      WRITE (lun2) NewIce
	close (lun2)
c**************************
         write (7,111) TAU,MONTH,MNTHDY,TOFDAY ! To PROTOCOL FILE
         print 111,tau,MONTH,MNTHDY,TOFDAY
 111     FORMAT (' ****  OUTREST RECORDS    TAU=',F9.2,
     1           1X,I3,1H/,I2,1H/,F6.2)
      RETURN
      END
