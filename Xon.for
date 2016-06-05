C                                            VERSION 26.02.88
C *****************************
C ******* OCEAN MODEL *********           (C) Computer Center
C ********* COUPLED ***********            A.V.Ganopolsky
C ***** ATMOSPHERE-OCEAN ******
C ********* MODEL *************
C
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE TERMO
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     4            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     5            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,AV,AH,HB,RAD,
     2            RAD2,DF1,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
c
	 include 'ice.fi'
C
       ALF(X)=R1+2.*R2*X
       GF(X) =CG*X*X*X
       DF(X) =CE*FC*X*X
       DR(X) =R1*X+R2*X*X
       P(X)  =6.1*DEXP(19.9*X/(X+273.)) !// not used ?
       SUM(X,Y)=X+ST*Y
       DATA EPSR/2.D-5/
C
C============== SPACE CYCLE =================================
C
       DO 10 I=1,Niz
       DO 10 J=2,Njz-1
C***********************
C      FOR ICE MODELLING
C***********************
       IF(MGEO(I,J).OR.MICE(I,J)) GOTO 10
c       IF(MGEO(I,J)) GOTO 10
c	 TICE  = 0.1d0       ! calculated in Celsius, C := K - 273.0
c	 TCICE = -1.5d0
	 TFreeze = 0.1d0      !TICE
c       if ((I+1).eq.Ibreak.and.(47-J).eq.Jbreak) then
c        tmpdata = 0					
c       endif
C***********************
C
C============== LOCAL MODEL =================================
C
       if (I.eq.Ibreak.and.J.eq.Jbreak) then
        tmpdata = 0					
       endif
       TS=TSM(I,J)
       TM=TMM(I,J)
       TH=THM(I,J)
       TB=TBM(I,J)
       HF=HFM(I,J)
       HL=HLM(I,J)
       QS=QSM(I,J)
       VS=VSM(I,J)
       W1=W1M(I,J)
       W2=W2M(I,J)
C
C***   IF(I.EQ. 4.AND.J.EQ.6)  PRINT 234,I,J,TS,TM,TH,TB,HS,HF,HL,QS,VS
  234  FORMAT(1X,2I3,' TS,TM,TH,TB',4F5.1,' HS,HF,HL,QS,VS',5F7.1)
C
C     CURRENTS, PARAMETRES
       FC=FCL(J)
       FM=FML(J)
       AL=ALF(TS)
       V=1.25D-3*VS
       CG=CG0*FM
       CK=1.
       IF(QS.LT.0.) CK=CK0
       GN=GF(V)
       DN=DF(V)
       Q=QS/CR
       QB=AV*(TH-TB)/(HB-HL)
       HE=CG*V/(CE*FC)
       DELT=TS-TH
       IF(DELT.LT.0.1) DELT=0.1
C     DETERMINATION OF EFFECTIVE DEPTH
       IF(Q.GT.0.) HT=GN/(DN+0.5*G*AL*Q)
       IF(Q.LE.0.) HT=HE
       IF(HF.LE.HT) GOTO 2
       IF(Q.GT.0.) GOTO 1
C     CONVECTION
       QH=-CK*Q
       WE=QH/DELT
       DFH=WE+W1
       GOTO 3
C     DETRAINMENT
 1     HF=HT
       QH=0.
       DFH=0.
       GOTO 3
C     ENTRAINMENT
 2     FP=1.-HF/HE
       WE=(2.*GN*FP/(G*HF)-CK*Q*AL)/(DR(TS)-DR(TH)+EPSR)
       DFH=WE+W1
       QH=DELT*WE
C     TIME STEP
 3     DFTS=(Q-QH)/HF
       DFTM=(Q-QB+W1*(TS-TH)+W2*(TH-TB))/HB
       HF=SUM(HF,DFH)
       TS=SUM(TS,DFTS)
C***********************
C      FOR ICE MODELLING
C***********************
       NEWICE(I+1,Njz+1-J) = TS.lt.TFreeze
C***********************
       TM=SUM(TM,DFTM)
C     OTHERS CONDITIONS
       IF(TS.LT.TB)   TS=TB
       TMAX=TS-5.*(TS-TB)/HB
       IF(TM.GT.TMAX) TM=TMAX
       IF(TM.GT.TB)   GOTO 4
       QB=(TM-TB)*HB/ST
       TM=TB
 4     CONTINUE
       IF(HF.GT.240.) HF=240.
       IF(Q.GT.0..AND.HF.GT.HT) HF=HT
       HS=HF
       IF(Q.GT.0.) HS=HF/FM
       IF(TH.EQ.TS) HS=HL
C     DETERMINATION OF DIMENSIONS
       QB=QB*CR
       TSM(I,J)=TS
       TMM(I,J)=TM
       HFM(I,J)=HF
       HSM(I,J)=HS
C
       QBM(I,J)=QB
 10    CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
              SUBROUTINE DYNAMO
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
       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,AV,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
       COMMON /DV/ I,J,KL,KL1
C
C**********************************************************************
C
       DO 1 J=2,Njz-2
        SN=1.
        IF(J.GT.8) SN=-1.
        FC=FCL(J)
       DO 1 I=1,Niz
       IF(MGEO(I,J).OR.MICE(I,J)) GOTO 1
        KL=1
        HF=HFM(I,J)
        HL=HLM(I,J)
C     GEOSTROPHIC WIND COMPONENT
        UGM(I,J)=0.01*G*(DIV(DPM,1)-DIV(DPM,3))*0.5/
     *           (FC*DF*RAD)*SN
        VGM(I,J)=0.01*G*(DIV(DPM,2)-DPM(I , J))/
     *           (FC*DL*COSL(J)*RAD)*SN
C     DRIFT WIND COMPONENT
        UDM(I,J)=TYM(I,J)/(FC*HF*R)*SN
        VDM(I,J)=-TXM(I,J)/(FC*HF*R)*SN
C     VERTICAL WIND COMPONENT
        KL=2
        ROT=5.D-3*((DIV(TYM,2)-DIV(TYM,4))/(DL*RAD*COSL(J))-
     *            (DIV(TXM,3)-DIV(TXM,1))/(DF*RAD))*SN
        W=-(ROT+TXM(I,J)/(TANL(J)*RAD))/(FC*R)
        W1M(I,J)=W+VGM(I,J)*HF/(TANL(J)*RAD)*SN
        W2M(I,J)=W+.5*VGM(I,J)*(HB+HL)/(TANL(J)*RAD)*SN
1     CONTINUE
      RETURN
      END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE ADIF
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)

      include 'paramz.fi'

       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     5            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     6            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
C
       CALL ADVECT(TSM,UDM,VDM)
       CALL TMCAL
       CALL ADVECT(TSM,UGM,VGM)
       CALL ADVECT(TMM,UGM,VGM)
       CALL ADVECT(HFM,UGM,VGM)
       CALL DIFUS(TSM)
       CALL DIFUS(TMM)
       CALL DIFUS(HFM)
       CALL THCAL
C
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE ADVECT(A,U,V)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,VK,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3            TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       DIMENSION   A(Niz,Njz),AZ(Niz,Njz),U(Niz,Njz),V(Niz,Njz)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
       COMMON /DV/ I,J,KL,KL1
C
       KL=0
C
C============== X-COMPONENT ==========================================
C
       DO 111 I=1,Niz
       DO 111 J=1,Njz
       IF(MGEO(I,J)) GOTO 111
       AZ(I,J)=A(I,J)
111    CONTINUE
       DO 1 J=3,Njz
        C2=COSL(J)*COSL(J)
       DO 1 I=1,Niz
       IF(MGEO(I,J).OR.MICE(I,J)) GOTO 1
        DX= (DIV(A,2)-DIV(A,4))/(2.*DL*COSL(J))
        DX2=(DIV(A,2)+DIV(A,4)-2.*A(I,J))/(DL2*C2)+TANL(J)*DX
        DX= DX/RAD
        DX2=DX2/RAD2
        DU=-U(I,J)*ST
        AZ(I,J)=A(I,J)+DU*DX+DU*DU*DX2*0.5
        IF(AZ(I,J).GT.240) AZ(I,J)=240.
C
 1     CONTINUE
C
C============== Y-COMPONENT ==========================================
C
       DO 2 I=1,Niz
       DO 2 J=3,Njz
       IF(MGEO(I,J).OR.MICE(I,J)) GOTO 2
        DY =(DIV(AZ,3)-DIV(AZ,1))/(2.*DF)
        DY2=(DIV(AZ,3)+DIV(AZ,1)-2.*AZ(I,J))/DF2
        DY =DY/RAD
        DY2=DY2/RAD2
        DV=-V(I,J)*ST
        A(I,J)=AZ(I,J)+DV*DY+DV*DV*DY2*0.5
        IF(A(I,J).GT.240.) A(I,J)=240.
C
 2     CONTINUE
       RETURN
       END
C
C**********************************************************************
C**********************************************************************
               SUBROUTINE DIFUS(A)
C**********************************************************************
C**********************************************************************
C
       IMPLICIT REAL*8(A-H,O-Z)
       LOGICAL MGEO,MICE

      include 'paramz.fi'

       COMMON /C/ PI,G,R,R1,R2,C,CR,CG0,CE,CK0,AV,AH,HB,RAD,
     2            RAD2,DF,DF2,DL,DL2,ST,SA(12),SY,SINL(Njz),COSL(Njz),
     3           TANL(Njz),FCL(Njz),FML(Njz),SZ(Njz),SZM
       DIMENSION   A(Niz,Njz),B(Niz,Njz)
       COMMON /MP/ MGEO(Niz,Njz),MICE(Niz,Njz)
       COMMON /DV/ I,J,KL,KL1
C
       KL=0
       DO 1 J=3,Njz
        C2=COSL(J)*COSL(J)
       DO 1 I=1,Niz
        IF(MGEO(I,J).OR.MICE(I,J)) GOTO 1
        DY= (DIV(A,1)-DIV(A,3))/(2.*DF)
        DX2=(DIV(A,2)+DIV(A,4)-2.*A(I,J))/(DL2*C2)+TANL(J)*DY
        DX2=DX2/RAD2
        DY2=(DIV(A,3)+DIV(A,1)-2.*A(I,J))/DF2
        DY2=DY2/RAD2
        DLAP=DX2+DY2
        B(I,J)=A(I,J)+AH*DLAP*ST
 1     CONTINUE
       DO 2 I=1,Niz
       DO 2 J=3,Njz
       IF(MGEO(I,J).OR.MICE(I,J)) GOTO 2
       A(I,J)=B(I,J)
 2     CONTINUE
       RETURN
       END
