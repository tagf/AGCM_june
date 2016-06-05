C**********************************************************************
C**********************************************************************
               SUBROUTINE INIT(K,MEC)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 KIT
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /N/ TAN(Niz,Njz,12),VSN(Niz,Njz,12),QRN(Niz,Njz,12),
     2            TXN(Niz,Njz, 4),TYN(Niz,Njz, 4),DPN(Niz,Njz),
     3            TS1(Niz,Njz),TS2(Niz,Njz)
       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),TXA(Niz,Njz),TYA(Niz,Njz)
       COMMON /Y/ TSY(Niz,Njz),TMY(Niz,Njz),HFY(Niz,Njz),HSY(Niz,Njz),
     2            QSY(Niz,Njz),QBY(Niz,Njz),TXY(Niz,Njz),TYY(Niz,Njz),
     3            DTyear(Niz,Njz)
       COMMON /Z/ TSZ(Njz,13),TMZ(Njz,13),HFZ(Njz,13),HSZ(Njz,13),
     2            QSZ(Njz,13),QBZ(Njz,13),TXZ(Njz,13),TYZ(Njz,13),
     3            SQZ(Njz,13),QAZ(Njz,13),QDZ(Njz,13),FQZ(Njz,13)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,AV,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON / QARY / P(Niz+2,Njz)  ,   U(Niz+2,Njz,2), V(Niz+2,Njz,2)
     *                ,T(Niz+2,Njz,2),  QW(Niz+2,Njz,2),GW(Niz+2,Njz)
     *               ,GT(Niz+2,Njz)  ,SNOAMT(Niz+2,Njz)
       COMMON /ATC/ TMAX(Niz,Njz),TMIN(Niz,Njz),ATCM(Niz,Njz)
       COMMON /SEA/ SS(Niz+2,Njz,12),SWIN(Niz+2,Njz,12)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C**********************************************************************
	include 'dir.fi'
C**********************************************************************
C
       DIMENSION LMEC(12),tmec(12)
       DATA LMEC/31,28,31,30,31,30,31,31,30,31,30,31/
C
       IF(K.EQ.2) GOTO 7
       IF(K.EQ.3) GOTO 4
C
C============== CONSTANTS =============================================
C
	 tmec=0.d0 !темп. глуб. ок. по месяцам в одной точке
       PI=3.14159
       G=9.81
       F=7.3D-5
       R=1.D3
       R1=7.35D-5
       R2=0.47D-5
       C=4.D3
       CR=4.D6
       CG0=10.
       CE=40.
       CK0=0.3
       AV=1.00D-4
       AH=2.D3
       HB=250.
       RAD=6.37D6
       DL=2.0d0*PI/Niz
       DF=PI/(Njz-1)
       RAD2=RAD*RAD
       DL2=DL*DL
       DF2=DF*DF
       TM=2.63D6
       ST=86400.	  !=24h=1day  (sec)	  Time step
       SY=1./365.
C
       DO 1 L=1,12
 1     SA(L)=1./LMEC(L)
C
C============== OTHER DIMENSIONS ===============================
C
        SZM=0.
       DO 11 J=1,Njz
c!!!!        SJ=8.5-J
        SJ=(Njz+1)/2.0-J	!учитывает шаг сетки
c        SJ=23.5-J	 ! надо помнить при смене сетки!!!
        FI=DABS(SJ)*PI/(Njz-1)
        COSL(J)=DCOS(FI)
        SINL(J)=DSIN(FI)
        TANL(J)=SINL(J)/COSL(J)
        SINJ=SINL(J)
        IF(SINJ.LT.0.1) SINJ=0.1
        FCL(J)=2.*F*SINJ
        FML(J)=2.2
        IF(DABS(FI).LT.0.5) FML(J)=7.2*DABS(FI)/PI+1.
        SZ(J)=0.
       DO 11 I=1,Niz
        IF(MGEO(I,J)) GOTO 11
         SZ(J)=SZ(J)+1.
         SZM=SZM+COSL(J)
 11    CONTINUE
       IF(K.EQ.0) GOTO 2
       IF(K.EQ.4) GOTO 4
C
C***   READ (12) NYER,MEC1,NDAY,TSM,TMM,HFM,TSA,TMA,HSA,QSA,QBA,
C*** *            TSY,TMY,QSY,QBY,DTyear,TSZ,HSZ,QSZ,QBZ,TMAX,TMIN,ATCM
C***   BACKSPACE 12
       Open (unit=12,file=BaseDir//'ocean\ocean.dat',status='OLD')
       READ (12,1100) NYER,MEC1,NDAY
c       READ (12,1200) TSM
       READ (12,1200) ((TSM(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TSM )
c       READ (12,1200) TMM
       READ (12,1200) ((TMM(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TMM )
c       READ (12,1200) HFM
       READ (12,1200) ((HFM(I,J),J=1,16),I=1,24)
	 call ocean_stretch( HFM )
c       READ (12,1200) TSA
       READ (12,1200) ((TSA(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TSA )
c       READ (12,1200) TMA
       READ (12,1200) ((TMA(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TMA )
c       READ (12,1200) HSA
       READ (12,1200) ((HSA(I,J),J=1,16),I=1,24)
	 call ocean_stretch( HSA )
c       READ (12,1200) QSA
       READ (12,1200) ((QSA(I,J),J=1,16),I=1,24)
	 call ocean_stretch( QSA )
c       READ (12,1200) QBA
       READ (12,1200) ((QBA(I,J),J=1,16),I=1,24)
	 call ocean_stretch( QBA )
c       READ (12,1200) TSY
       READ (12,1200) ((TSY(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TSY )
c       READ (12,1200) TMY
       READ (12,1200) ((TMY(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TMY )
c       READ (12,1200) QSY
       READ (12,1200) ((QSY(I,J),J=1,16),I=1,24)
	 call ocean_stretch( QSY )
c       READ (12,1200) QBY
       READ (12,1200) ((QBY(I,J),J=1,16),I=1,24)
	 call ocean_stretch( QBY )
c       READ (12,1200) DTyear
       READ (12,1200) ((DTyear(I,J),J=1,16),I=1,24)
	 call ocean_stretch( DTyear )
C       READ (12,1300) TSZ
       READ (12,1300) ((TSZ(I,J),J=1,13),I=1,16)
C       READ (12,1300) HSZ
       READ (12,1300) ((HSZ(I,J),J=1,13),I=1,16)
C       READ (12,1300) QSZ
       READ (12,1300) ((QSZ(I,J),J=1,13),I=1,16)
C       READ (12,1300) QBZ
       READ (12,1300) ((QBZ(I,J),J=1,13),I=1,16)
C       READ (12,1200) TMAX
       READ (12,1200) ((TMAX(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TMAX )
C       READ (12,1200) TMIN
       READ (12,1200) ((TMIN(I,J),J=1,16),I=1,24)
	 call ocean_stretch( TMIN )
C       READ (12,1200) ATCM
       READ (12,1200) ((ATCM(I,J),J=1,16),I=1,24)
	 call ocean_stretch( ATCM )
       Close (unit=12)
C
1100   FORMAT (1x,i8,1x,i8,1x,i8)
1200   FORMAT (6(1x,d12.5))
1300   FORMAT (4(1x,d12.5))

       DO 111 J=1,Njz
       DO 110 I=1,Niz
c***********
c       IF(MGEO(I,J)) GOTO 110
       IF(MGEO(I,J).or.MICE(I,J)) GOTO 110
       GT(I+1,Njz+1-J)=TSM(I,J)+273.
  110  CONTINUE
       GT(1,Njz+1-J)=GT(Niz+1,Njz+1-J)
       GT(Niz+2,Njz+1-J)=GT(2,Njz+1-J)
  111  CONTINUE
C
c       CALL PRIN(TSM,1.D0,1)
c       CALL PRIN(TMM,1.D0,2)
c       CALL PRIN(HFM,1.D0,3)
C
 2     CONTINUE
c 1000   FORMAT(5X,'*** AAA =',F10.3)   not used
C============== READING DINENSIONS TBM,DPM ===============================
C
       DO 101 I=1,Niz
       DO 101 J=1,Njz
       DPM(I,J)=0.
 101   TBM(I,J)=0.
C
      Open (unit=14,file=TRIM(BaseDir)//'ocean\bosa.dat',status='OLD')
c!!!      READ (14,22) ((TBM(I,Njz+1-J),J=4,Njz-1),I=1,Niz)
c****************************** interpolating datas
      READ (14,22) ((TBM(I,17-J),J=4,15),I=1,24)
	call ocean_stretch( TBM )
c******************************
c     вместо этого берем из sst - вычисляем min за год (за 12 мес)
	do i=1,Niz
	 do j=1,Njz
	  do k1=1,12
	    tmec(k1) = ss(i+1,Njz+1-j,k1)
	   if (ss(i+1,Njz+1-j,k1).gt.271.) then	 !sea ice
	    tmec(k1) = 0.
	   endif
	  enddo

	  ssmin = 1000.d0
	  do k1=1,12
	   if (ssmin.gt.tmec(k1)) then
	    ssmin = tmec(k1)
	   endif
	  enddo

	   TBM(I,J) = ssmin
	  if (TBM(I,J).lt.1.1d0) then
	   TBM(I,J) = 1.0d0  
	  endif
	 enddo
	enddo
c******************************


c      PRINT 22,((TBM(I,Njz+1-J),J=1,Njz),I=1,Niz)
c       READ (14,23) ((DPM(I,Njz+1-J),J=3,Njz-1),I=1,Niz)
c****************************** interpolating data
      READ (14,23) ((DPM(I,17-J),J=3,15),I=1,24)
	do i=1,24
	 do j=1,16
	  DPM(i,j) = DMAX1( DPM(i,j), 30.0d0 ) ! минимум 30 см.
	 enddo
	enddo
	call ocean_stretch( DPM )
c******************************
C***  PRINT 23,((DPM(I,Njz+1-J),J=3,Njz-1),I=1,Niz)
       Close (unit=14)
 22    FORMAT(23F4.0)
 23    FORMAT(13F4.0)
C
C============== INITIAL CONDITIONS ====================================
C
       DO 3 I=1,Niz
       DO 3 J=1,Njz
       IF(MGEO(I,J)) GOTO 3
       IF(K.EQ.1)  GOTO 3
        TSM(I,J)=SS(I+1,Njz+1-J,MEC)
        IF(TSM(I,J).LT.TBM(I,J)) TSM(I,J)=TBM(I,J)
        IF(TSM(I,J).GT.  100.  ) TSM(I,J)=TBM(I,J)
        TMM(I,J)=(TSM(I,J)+TBM(I,J))*0.5
        HFM(I,J)=50.
c 33    CONTINUE
        HSM(I,J)=50.
        UGM(I,J)=0.
        VGM(I,J)=0.
        UDM(I,J)=0.
        VDM(I,J)=0.
        W1M(I,J)=0.
        W2M(I,J)=0.
 3     CONTINUE
        CALL THCAL
 4     CONTINUE
       CALL THCAL
       IF(K.EQ.1) GOTO 9
C
C============== SET 0 ANNUAL AND MONTHLY AVERAGE ==============
C
       DO 5 I=1,Niz
       DO 5 J=1,Njz
       IF(MGEO(I,J)) GOTO 5
        TSY(I,J)=0.
        TMY(I,J)=0.
        HFY(I,J)=0.
        HSY(I,J)=0.
        QSY(I,J)=0.
        QBY(I,J)=0.
        TXY(I,J)=0.
        TYY(I,J)=0.
        TMAX(I,J)=0.
        TMIN(I,J)=1000.
        DTyear(I,J)=0.
 5     CONTINUE
       DO 6 J=1,Njz
       DO 6 M=1,13
        TSZ(J,M)=0.
        TMZ(J,M)=0.
        HFZ(J,M)=0.
        HSZ(J,M)=0.
        QSZ(J,M)=0.
        QBZ(J,M)=0.
        TXZ(J,M)=0.
        TYZ(J,M)=0.
        QDZ(J,M)=0.
        SQZ(J,M)=0.
 6      FQZ(J,M)=0.
 7     CONTINUE
       DO 8 I=1,Niz
       DO 8 J=1,Njz
C
        TSA(I,J)=0.
        TMA(I,J)=0.
        HFA(I,J)=0.
        HSA(I,J)=0.
        QSA(I,J)=0.
        QBA(I,J)=0.
        TXA(I,J)=0.
        TYA(I,J)=0.
 8     CONTINUE
 9     CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE ZONAL(MEC)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),TXA(Niz,Njz),TYA(Niz,Njz)
       COMMON /Z/ TSZ(Njz,13),TMZ(Njz,13),HFZ(Njz,13),HSZ(Njz,13),
     2            QSZ(Njz,13),QBZ(Njz,13),TXZ(Njz,13),TYZ(Njz,13),
     3            SQZ(Njz,13),QAZ(Njz,13),QDZ(Njz,13),FQZ(Njz,13)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,VK,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       M=MEC
       DO 1 I=1,Niz
       DO 1 J=1,Njz
       IF(MGEO(I,J)) GOTO 1
        TSZ(J,M)=TSZ(J,M)+TSA(I,J)/SZ(J)
        TMZ(J,M)=TMZ(J,M)+TMA(I,J)/SZ(J)
C       HFZ(J,M)=HFZ(J,M)+HFA(I,J)/SZ(J)
        HSZ(J,M)=HSZ(J,M)+HSA(I,J)/SZ(J)
        QSZ(J,M)=QSZ(J,M)+QSA(I,J)/SZ(J)
        QBZ(J,M)=QBZ(J,M)+QBA(I,J)/SZ(J)
C       TXZ(J,M)=TXZ(J,M)+TXA(I,J)/SZ(J)
C       TYZ(J,M)=TYZ(J,M)+TYA(I,J)/SZ(J)
 1     CONTINUE
       IF(MEC.NE.12) GOTO 7
C
C============== ANNUAL VALUES ================================
C
       DO 2 J=1,Njz
       DO 2 N=1,12
        TSZ(J,13)=TSZ(J,13)+TSZ(J,N)/12.
        TMZ(J,13)=TMZ(J,13)+TMZ(J,N)/12.
C       HFZ(J,13)=HFZ(J,13)+HFZ(J,N)/12.
        HSZ(J,13)=HSZ(J,13)+HSZ(J,N)/12.
        QSZ(J,13)=QSZ(J,13)+QSZ(J,N)/12.
        QBZ(J,13)=QBZ(J,13)+QBZ(J,N)/12.
C       TXZ(J,13)=TXZ(J,13)+TXZ(J,N)/12.
C       TYZ(J,13)=TYZ(J,13)+TYZ(J,N)/12.
 2      CONTINUE
C
C============== ANNUAL THERMOCYCLE ===========================
C
       DO 3 N=1,12
       DO 3 J=1,Njz
 3     SQZ(J,13)=SQZ(J,13)+TMZ(J,N)*HB*SZ(J)*COSL(J)*5.E-3/12.
C
       DO 4 J=1,Njz
       DO 4 N=1,12
 4     SQZ(J,N)=TMZ(J,N)*HB*SZ(J)*COSL(J)*5.E-3-SQZ(J,13)
C
C============== ANNUAL GLOBAL CHARACTERISTICS ===============
C
        TSG=0.D0
        TMG=0.D0
        QSG=0.D0
        QBG=0.D0
C
C
       DO 5 J=1,Njz
        TSG=TSG+TSZ(J,13)*SZ(J)*COSL(J)/SZM
        TMG=TMG+TMZ(J,13)*SZ(J)*COSL(J)/SZM
        QSG=QSG+QSZ(J,13)*SZ(J)*COSL(J)/SZM
        QBG=QBG+QBZ(J,13)*SZ(J)*COSL(J)/SZM
C
 5     CONTINUE
C
C============== MERIDIONAL THERMOFLOW ==========================
C
       DO 6 J=1,Njz
       DO 6 I=1, J
 6     FQZ(J,13)=FQZ(J,13)+(QSZ(I,13)-QSG)*SZ(I)*COSL(I)
       PRINT 10,TSG,TMG,QSG,QBG
 10    FORMAT(//10X,'***GLOBAL*** TS=',F7.2,'  TM=',F7.2,
     * '  QS=',F7.1,'  QB=',F7.1//)
 7     RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE AVERAG(MEC)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),TXA(Niz,Njz),TYA(Niz,Njz)
       COMMON /Y/ TSY(Niz,Njz),TMY(Niz,Njz),HFY(Niz,Njz),HSY(Niz,Njz),
     2            QSY(Niz,Njz),QBY(Niz,Njz),TXY(Niz,Njz),TYY(Niz,Njz),
     3            DTyear(Niz,Njz)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,VK,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
       IF(MGEO(I,J)) GOTO 1
        TSA(I,J)=TSA(I,J)+TSM(I,J)*SA(MEC)
        TMA(I,J)=TMA(I,J)+TMM(I,J)*SA(MEC)
        HFA(I,J)=HFA(I,J)+HFM(I,J)*SA(MEC)
        HSA(I,J)=HSA(I,J)+HSM(I,J)*SA(MEC)
        QSA(I,J)=QSA(I,J)+QSM(I,J)*SA(MEC)
C       QBA(I,J)=QBA(I,J)+QBM(I,J)*SA(MEC)
C       TXA(I,J)=TXA(I,J)+TXM(I,J)*SA(MEC)
C       TYA(I,J)=TYA(I,J)+TYM(I,J)*SA(MEC)
        TSY(I,J)=TSY(I,J)+TSM(I,J)*SY
        TMY(I,J)=TMY(I,J)+TMM(I,J)*SY
        HFY(I,J)=HFY(I,J)+HFM(I,J)*SY
C       HSY(I,J)=HSY(I,J)+HSM(I,J)*SY
        QSY(I,J)=QSY(I,J)+QSM(I,J)*SY
C       QBY(I,J)=QBY(I,J)+QBM(I,J)*SY
C       TXY(I,J)=TXY(I,J)+TXM(I,J)*SY
C       TYY(I,J)=TYY(I,J)+TYM(I,J)*SY
 1     CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE TMCAL
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,VK,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
       IF (MGEO(I,J)) GOTO 1
       TMM(I,J)=(HLM(I,J)*TSM(I,J)+0.5*(THM(I,J)+TBM(I,J))*
     *          (HB-HLM(I,J)))/HB
 1     CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE THCAL
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,VK,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
       IF(MGEO(I,J)) GOTO 1
       TMAX=TSM(I,J)-(TSM(I,J)-TBM(I,J))/(2.*HB)
       IF(TMM(I,J).GT.TMAX) TMM(I,J)=TMAX
        THM(I,J)=TMM(I,J)*HB-TSM(I,J)*HFM(I,J)-TBM(I,J)*(HB-HFM(I,J))/2.
        THM(I,J)=THM(I,J)*2./(HB-HFM(I,J))
        HLM(I,J)=HFM(I,J)
        IF(THM(I,J).LT.TBM(I,J)) THM(I,J)=TBM(I,J)
        IF(THM(I,J).LE.TSM(I,J)) GOTO 1
        THM(I,J)=TSM(I,J)
        IF(TSM(I,J)-TBM(I,J).LT.0.01) GOTO 2
        HLM(I,J)=2.*TMM(I,J)*HB-(TSM(I,J)+TBM(I,J))*HB
        HLM(I,J)=HLM(I,J)/(TSM(I,J)-TBM(I,J))
        IF(HLM(I,J).GT.240.) HLM(I,J)=240.
        GOTO 1
 2      HLM(I,J)=HFM(I,J)
 1     CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE DIF(MEC)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),QLA(Niz,Njz),QTA(Niz,Njz)
       COMMON /Y/ TSY(Niz,Njz),TMY(Niz,Njz),HFY(Niz,Njz),HSY(Niz,Njz),
     2            QSY(Niz,Njz),QBY(Niz,Njz),TXY(Niz,Njz),TYY(Niz,Njz),
     3            DTyear(Niz,Njz)
       COMMON /SEA/ SS(Niz+2,Njz,12),SWIN(Niz+2,Njz,12)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
C**********************************************************************
C
       DIMENSION DTS(Niz,Njz)
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
        IF(MGEO(I,J)) GOTO 1
        SST=SS(I+1,Njz+1-J,MEC)
        IF(SST.GT.200.) SST=0.
        DTS(I,J)=TSA(I,J)-SST
 1     CONTINUE
       DO 2 I=1,Niz
       DO 2 J=1,Njz
        IF(MGEO(I,J))GOTO 2
        DTyear(I,J)=DTyear(I,J)+DTS(I,J)/12.
 2     CONTINUE
c       CALL PRIN(DTS,1.D0,7)
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE ATCMAP(MEC)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /A/ TSA(Niz,Njz),TMA(Niz,Njz),HFA(Niz,Njz),HSA(Niz,Njz),
     2            QSA(Niz,Njz),QBA(Niz,Njz),QLA(Niz,Njz),QTA(Niz,Njz)
       COMMON /ATC/ TMAX(Niz,Njz),TMIN(Niz,Njz),ATCM(Niz,Njz)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       DO 10 I=1,Niz
       DO 10 J=1,Njz
       IF(MGEO(I,J)) GOTO 10
        IF(TMAX(I,J).LT.TMA(I,J)) TMAX(I,J)=TMA(I,J)
        IF(TMIN(I,J).GT.TMA(I,J)) TMIN(I,J)=TMA(I,J)
 10    CONTINUE
       IF(MEC.NE.12) GOTO 30
       DO 20 I=1,Niz
       DO 20 J=1,Njz
       IF(MGEO(I,J)) GOTO 20
        ATCM(I,J)=10.*(TMAX(I,J)-TMIN(I,J))
 20    CONTINUE
c       CALL PRIN(ATCM,1.D0,8)
 30    RETURN
       END
C**********************************************************************
C**********************************************************************
               SUBROUTINE PRINZ(A,K)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)

      include 'paramz.fi'

       DIMENSION A(Njz,13)
       DIMENSION MON(13),NAM(17)
C
       DATA MON/'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     2          'JULY','AUG ','SEP ','OCT ','NOV ','DEC ','YEAR'/
       DATA NAM/'TS','TM','HF','HS','QS','QB','DTy','SQ','TX','TY',
     2          'UG','VG','UD','VD','W1','W2','FQ'/
       PRINT 10,NAM(K)
C      IF(K.NE.10) PRINT 11,((MON(M),(A(J,M),J=1,Njz)),M=1,13)
       IF(K.EQ.10) PRINT 11,MON(13),(A(J,13),J=1,Njz)
 10    FORMAT(1X,50('*'),' ANNUAL CYCLE ',A6,50('*')//
     * 12X,'90N',4X,'78N',4X,'66N',4X,'54N',4X,'42N',4X,'30N',4X,
     * '18N',4X,' 6N',4X,' 6S',4X,'18S',4X,'30S',4X,'42S',4X,'54S',
     * 4X,'66S',4X,'78X',4X,'90S'/1X,118('='))
 11    FORMAT(12(1X,A6,16F7.1/)/1X,A6,16F7.1)
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE PRIN(A,S,K)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       DIMENSION A(Niz,Njz),LX(12),LY(Njz),MA(Niz,Njz),NAM(17)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
C
       DATA NI/'I'/
       DATA NAM/'TS','TM','HF','HS','QS','QB','DTy','AT','TX','TY',
     *          'UG','VG','UD','VD','W1','W2','FQ'/
       DATA LX/180,150,120,90,60,30,0,30,60,90,120,150/,
c     *      LY/90,78,66,54,42,30,18,6,6,18,30,42,54,66,78,90/
c     данные по широтам
     1      LY/90,86,82,78,74,70,66,62,58,54,50,46,42,38,34,30,26,22
     2       ,18,14,10,6,2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62
     3       ,66,70,74,78,82,86,90/
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
        if (i.eq.3.and.j.eq.Njz-1) then
         iii = i
        end if
       IF(MGEO(I,J)) GOTO 2
        MA(I,J)=A(I,J)*S
        GOTO 1
 2      MA(I,J)=100000
 1     CONTINUE
       PRINT 11,NAM(K)
       PRINT 12
       PRINT 13,LX
       PRINT 16
       PRINT 15,(LY(I),NI,(MA(J,I),J=1,Niz),NI,LY(I),I=1,Njz)
       PRINT 16
       PRINT 13,LX
       PRINT 12
 11    FORMAT(/55X,'DIMENSION',A6/)
 12    FORMAT(1X,104('*'))
 13    FORMAT(1X,'*',I8,11I8,6X,'*')
 16    FORMAT(1X,'*',2X,98('-'),2X,'*')
 15    FORMAT(1X,'*',I2,A1,24I4,A1,I2,'*')
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               FUNCTION DIV(A,M)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       DIMENSION A(Niz,Njz),IC(4),JC(4)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
       COMMON /DV/ I,J,KL,KL1
C
       DATA IC/0,1,0,-1/,JC/1,0,-1,0/
       IP=I+IC(M)
       JP=J+JC(M)
         IF(IP.EQ.0) IP=Niz
         IF(IP.EQ.Niz+1) IP=1
       IF(MGEO(IP,JP)) GOTO 1
       B=A(IP,JP)
       GOTO 3
C
  1    CONTINUE
C
       IF(KL.NE.1) GOTO 2
       IN=I-IC(M)
       JN=J-JC(M)
         IF(IN.EQ.0)  IN=Niz
         IF(IN.EQ.Niz+1) IN= 1
       IF(MGEO(IN,JN)) GOTO 2
       B=2.*A(I,J)-A(IN,JN)
       GOTO 3
C
  2    CONTINUE
C
       IF(KL.NE.2) B=A(I,J)
       IF(KL.EQ.2) B=0.
C
  3    DIV=B
       RETURN
       END


C**********************************************************************
C**********************************************************************
               SUBROUTINE OCEAN_STRETCH(A)
C**********************************************************************
C**********************************************************************
c растяжка массивов в 3 раза. Вход:массив, размер по X (72 или 74)
      IMPLICIT REAL*8(A-H,O-Z)
c      LOGICAL MGEO,MICE
      include 'paramz.fi'
      DIMENSION A(Niz,Njz)
	do J=Njz,1,-1
	 do I=Niz,1,-1
	  A(I,J) = A( (i+2)/3, j/3+1 )
	 enddo
	enddo
      RETURN
      END
