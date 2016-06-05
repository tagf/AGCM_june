C     ****************
C     ****************
      program difference
C     ****************
C     ****************
C
C
C             READS TWO POINTS INITIAL CONDITIONS
C                AND CALCULATES DIFFERENCE
C
C************ BEGINNING OF COMMON ************
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
     6 JM,IM,ID,MNTHDY,SDEDY ,NCYCLE,NC3   ,MONTH ,MRCH  ,NSTEP,nav
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1  ,SD(74,46)
      DIMENSION PHIS1(74,46),ISFTYP1(74,46),TS1(74,46)
     1  ,SD1(74,46)
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
      DIMENSION P1(74,46),U1(74,46,2),V1(74,46,2),T1(74,46,2)
     1  ,QW1(74,46,2),GW1(74,46),GT1(74,46),SNOAMT1(74,46)
C
	include 'dir.fi'
C************ END OF COMMON ******************
C
	 character*30 restfile
C           
            WorkDir='id##4x5d'
            restfile='E:\CLIMATE\id##4x5d\rest12'
      OPEN (9,FILE=restfile,FORM='UNFORMATTED')
      rewind 9
      READ (9) ((PHIS(I,J),I=2,73),J=1,46)
      READ (9) C
      READ (9) ((ISFTYP(I,J),I=2,73),J=1,46)
      READ (9) ((P(I,J),I=2,73),J=1,46)
      READ (9) ((U(I,J,1),I=2,73),J=1,46)
      READ (9) ((U(I,J,2),I=2,73),J=1,46)
      READ (9) ((V(I,J,1),I=2,73),J=1,46)
      READ (9) ((V(I,J,2),I=2,73),J=1,46)
      READ (9) ((T(I,J,1),I=2,73),J=1,46)
      READ (9) ((T(I,J,2),I=2,73),J=1,46)
      READ (9) ((QW(I,J,1),I=2,73),J=1,46)
      READ (9) ((QW(I,J,2),I=2,73),J=1,46)
      READ (9) ((GW(I,J),I=2,73),J=1,46)
      READ (9) ((GT(I,J),I=2,73),J=1,46)
      READ (9) ((SNOAMT(I,J),I=2,73),J=1,46)
      READ (9) ((TS(I,J),I=2,73),J=1,46)
      READ (9) ((SD(I,J),I=2,73),J=1,46)
      close (9)
C
         PRINT 111,restfile,TAU
 111     FORMAT (' INPUT FIRST INITIAL DATA. FILE ',a,' TAU=',F9.2)
C
            restfile='E:\CLIMATE\init#8\rest2'
      OPEN (9,FILE=restfile,FORM='UNFORMATTED')
      rewind 9
      READ (9) ((PHIS1(I,J),I=2,73),J=1,46)
      READ (9) C
      READ (9) ((ISFTYP1(I,J),I=2,73),J=1,46)
      READ (9) ((P1(I,J),I=2,73),J=1,46)
      READ (9) ((U1(I,J,1),I=2,73),J=1,46)
      READ (9) ((U1(I,J,2),I=2,73),J=1,46)
      READ (9) ((V1(I,J,1),I=2,73),J=1,46)
      READ (9) ((V1(I,J,2),I=2,73),J=1,46)
      READ (9) ((T1(I,J,1),I=2,73),J=1,46)
      READ (9) ((T1(I,J,2),I=2,73),J=1,46)
      READ (9) ((QW1(I,J,1),I=2,73),J=1,46)
      READ (9) ((QW1(I,J,2),I=2,73),J=1,46)
      READ (9) ((GW1(I,J),I=2,73),J=1,46)
      READ (9) ((GT1(I,J),I=2,73),J=1,46)
      READ (9) ((SNOAMT1(I,J),I=2,73),J=1,46)
      READ (9) ((TS1(I,J),I=2,73),J=1,46)
      READ (9) ((SD1(I,J),I=2,73),J=1,46)
      close (9)
C
         PRINT 112,restfile,TAU
 112     FORMAT (' INPUT SECOND INITIAL DATA. FILE ',a,' TAU=',F9.2)
C
      do 1 i=2,73
      do 1 j=1,46
      phis(i,j)=phis(i,j)-phis1(i,j)
      p(i,j)=p(i,j)-p1(i,j)
      u(i,j,1)=u(i,j,1)-u1(i,j,1)
      u(i,j,2)=u(i,j,2)-u1(i,j,2)
      v(i,j,1)=v(i,j,1)-v1(i,j,1)
      v(i,j,2)=v(i,j,2)-v1(i,j,2)
      t(i,j,1)=t(i,j,1)-t1(i,j,1)
      t(i,j,2)=t(i,j,2)-t1(i,j,2)
      qw(i,j,1)=qw(i,j,1)-qw1(i,j,1)
      qw(i,j,2)=qw(i,j,2)-qw1(i,j,2)
      gw(i,j)=gw(i,j)-gw1(i,j)
      gt(i,j)=gt(i,j)-gt1(i,j)
      snoamt(i,j)=snoamt(i,j)-snoamt1(i,j)
      ts(i,j)=ts(i,j)-ts1(i,j)
 1    sd(i,j)=sd(i,j)-sd1(i,j)
c 
c      CALL WRITIJ(PHIS,1.d0,'PHIS',0)
      CALL WRITIJ(P   ,1.d0,'P   ',0)
      CALL WRITIJ(U(1,1,1),1.d0,'U1  ',2)
      CALL WRITIJ(U(1,1,2),1.d0,'U3  ',2)
      CALL WRITIJ(V(1,1,1),1.d0,'V1  ',2)
      CALL WRITIJ(V(1,1,2),1.d0,'V3  ',2)
      CALL WRITIJ(T(1,1,1),1.d0,'T1  ',0)
      CALL WRITIJ(T(1,1,2),1.d0,'T3  ',0)
      CALL WRITIJ(QW(1,1,1),1.0d-4,'QW1 ',0)
      CALL WRITIJ(QW(1,1,2),1.0d-4,'QW3 ',0)
      CALL WRITIJ(GW  ,1.0d-4,'GW  ',0)
      CALL WRITIJ(GT  ,1.d0,'GT  ',0)
      CALL WRITIJ(SNOAMT,1.0d-4,'SNOW',0)
      CALL WRITIJ(TS    ,1.d0,'TS  ',0)
      CALL WRITIJ(SD    ,1.0d9,'SD  ',0)
c      call PrintC
      stop 'Normal end'
      END
