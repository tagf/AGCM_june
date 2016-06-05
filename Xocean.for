C
C**********************************************************************
C**********************************************************************
               SUBROUTINE OCEAN
C**********************************************************************
C**********************************************************************
C
	include 'recom.fi'
       REAL*8 K
       LOGICAL MGEO,MICE
C
      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     4            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     5            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
C
       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),TXA(Niz,Njz),TYA(Niz,Njz)
       COMMON /Y/ TSY(Niz,Njz),TMY(Niz,Njz),HFY(Niz,Njz),HSY(Niz,Njz),
     2            QSY(Niz,Njz),QBY(Niz,Njz),TXY(Niz,Njz),TYY(Niz,Njz),
     3             DTyear(Niz,Njz)
       COMMON /Z/ TSZ(Njz,13),TMZ(Njz,13),HFZ(Njz,13),HSZ(Njz,13),
     2            QSZ(Njz,13),QBZ(Njz,13),TXZ(Njz,13),TYZ(Njz,13),
     3            SQZ(Njz,13),QAZ(Njz,13),QDZ(Njz,13),FQZ(Njz,13)
C
       COMMON / QARY / P(Niz+2,Njz)  ,   U(Niz+2,Njz,2), V(Niz+2,Njz,2)
     *                ,T(Niz+2,Njz,2),  QW(Niz+2,Njz,2),GW(Niz+2,Njz)
     *               ,GT(Niz+2,Njz)  ,SNOAMT(Niz+2,Njz)
       COMMON /ATC/ TMAX(Niz,Njz),TMIN(Niz,Njz),ATCM(Niz,Njz)
       COMMON /MP/  MGEO(Niz,Njz),MICE(Niz,Njz)
       COMMON /TEST/ S4O(Niz,Njz),R4O(Niz,Njz),T4O(Niz,Njz)
C**********************************************************************
C
c       NYER = 0
c       MEC  = MONTH
c       NDAY = MNTHDY
       CALL BUFER(NYER,MEC,NDAY,1)
ccccc       PRINT 20,NYER,MEC,NDAY
 20    FORMAT(10X,3I5)
C*******************
C      PRINT 200,TSM(1,8),TMM(1,8),HFM(1,8),QSM(1,8),S4O(1,8),R4O(1,8),
C    * T4O(1,8)
200   FORMAT(1X,'TSM=',D13.5,' TMM=',D13.5,' HFM=',D13.5,' QSM=',D13.5,
     * /,1X,' S4=',D13.5,'  R4=',D13.5,'  T4=',D13.5)
       CALL DYNAMO
C*******************
C      CALL PRIN(UGM,100.,11)
C      CALL PRIN(VGM,100.,12)
C      CALL PRIN(UDM,100.,13)
C      CALL PRIN(VDM,100.,14)
C      CALL PRIN(W1M,1.D7,15)
C      CALL PRIN(W2M,1.D7,16)
C      CALL PRIN(QSM,1.,5)
C******************************
       DO 111 I=1,Niz
       DO 111 J=1,Njz
       IF(MGEO(I,J)) GOTO 111
       IF(HFM(I,J).GT.240.) HFM(I,J)=240.
       IF(HFM(I,J).LT. 10.) HFM(I,J)= 10.
C      TMIN(I,J)=100.
C      TMAX(I,J)=0.0
C      ATCM(I,J)=0.0
  111  CONTINUE
C
       CALL ADIF
       CALL TERMO
       CALL THCAL
       CALL AVERAG(MEC)
C
C***   WRITE (12) NYER,MEC,NDAY,TSM,TMM,HFM,TSA,TMA,HSA,QSA,QBA,
C*** *            TSY,TMY,QSY,QBY,DT,TSZ,HSZ,QSZ,QBZ,TMAX,TMIN,ATCM
C***   BACKSPACE 12
       DO 3 J=1,Njz
       DO 2 I=1,Niz
        IF(MICE(I,J)) TSM(I,J)=TBM(I,J)
        IF(MGEO(I,J).OR.MICE(I,J)) GOTO 2
        GT(I+1,Njz+1-J)=TSM(I,J)+273.
 2     CONTINUE
        GT( 1,Njz+1-J)=GT(Niz+1,Njz+1-J)
        GT(Niz+2,Njz+1-J)=GT( 2,Njz+1-J)
 3     CONTINUE
C
       IF(NDAY.NE.1) GOTO 1
       MC=MEC-1
       IF(MC.EQ.0) MC=12
       PRINT 10,MC
C      CALL PRIN(TSA,1.,1)
C      CALL PRIN(HSA,1.,4)
C      CALL PRIN(QSA,1.,5)
       CALL ZONAL(MC)
       CALL ATCMAP(MC)
       CALL DIF(MC)
C***   IF(MC.EQ.2.OR.MC.EQ.8) WRITE(13) HSA
       CALL INIT(2,MC)
C
       IF(MC.NE.12) GOTO 1
C
       PRINT 11,NYER
C      CALL PRIN(TSY,1.,1)
C      CALL PRIN(QSY,1.,5)
c       CALL PRIN(DTyear ,1.,7)
       CALL PRINZ(TSZ,1)
       CALL PRINZ(HSZ,4)
       CALL PRINZ(QSZ,5)
C***   WRITE(13) NYER,DTyear,QSY,ATCM
       CALL INIT(3,MC)
C
 1     CONTINUE
C
 10    FORMAT(30X,'MONTH',I3)
 11    FORMAT(/28X,'* YEAR *',I4/)
C
       RETURN
       END
