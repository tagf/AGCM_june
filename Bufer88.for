C
C**********************************************************************
C**********************************************************************
               SUBROUTINE BUFER(NYER,MEC,NDAY,K)
C**********************************************************************
C**********************************************************************
C
	include 'recom.fi'
       LOGICAL MGEO,MICE
C
      include 'paramz.fi'
C************** COMMON ************************************************
C
c     7 ,ARINDX(9,9),GRWET(4,7),GRCALB(6,3),FOEMIS(9,2),GREMIS(4,2),
c     8 IGRTEX(Niz+2,Njz),IGRCOL(Niz+2,Njz)
C
      COMMON / MISC / PHIS(Niz+2,Njz),ISFTYP(Niz+2,Njz),TS(Niz+2,Njz)
     1   ,SD(Niz+2,Njz)
C
      COMMON /OSA/ QSN(Niz,Njz),TXN(Niz,Njz),TYN(Niz,Njz)
C
      COMMON /WIND/ WIND(Niz+2,Njz)
C
cccc      COMMON /PSEA/ PSEA(Niz+2,Njz)
C
       COMMON /M/ TSM(Niz,Njz),TMM(Niz,Njz),THM(Niz,Njz),TBM(Niz,Njz),
     2            HFM(Niz,Njz),HSM(Niz,Njz),HLM(Niz,Njz),DPM(Niz,Njz),
     3            QSM(Niz,Njz),QBM(Niz,Njz),VSM(Niz,Njz),TXM(Niz,Njz),
     4            TYM(Niz,Njz),W1M(Niz,Njz),W2M(Niz,Njz),QAM(Niz,Njz),
     5            UGM(Niz,Njz),VGM(Niz,Njz),UDM(Niz,Njz),VDM(Niz,Njz)
C
      COMMON /MP/  MGEO(Niz,Njz),MICE(Niz,Njz)
C
C************** END COMMON ********************************************
C
C           TIME FOR OCEAN
C
       NYER=SDEYR
       MEC =MONTH
       NDAY=MNTHDY
C
C     PANAMA
C
c       ISFTYP(8,10)=5
c       PHIS(8,10)=0.
C
C              NEW MAPS
C
       DO 1 I=1,Niz
       DO 1 J=1,Njz
       IS=ISFTYP(I+1,Njz+1-J)
       MGEO(I,J)=.FALSE.
       MICE(I,J)=.FALSE.
       IF(IS.NE.7.AND.IS.NE.9) MGEO(I,J)=.TRUE.  !land
       IF(IS.EQ.9)             MICE(I,J)=.TRUE.  !sea ice
C            WIND SPEED AND HEAT FOR OCEAN
       IF(K.EQ.0) GOTO 1
       VSM(I,J)=WIND(I+1,Njz+1-J)
       QSM(I,J)=QSN(I,J)
 1     CONTINUE
c       ISFTYP(8,10)=7
       IF(K.EQ.0) GOTO 4
C
       DO 2 I=1,Niz	 ! k .ne. 0
       IM1=I-1
       IF(IM1.EQ.0) IM1=Niz
       DO 2 J=3,Njz-1
       TXM(I,Njz+1-J)=(TXN(IM1,J)+TXN(IM1,J+1)+TXN(I,J)+TXN(I,J+1))*0.25
       TYM(I,Njz+1-J)=(TYN(IM1,J)+TYN(IM1,J+1)+TYN(I,J)+TYN(I,J+1))*0.25
 2     CONTINUE
       DO 3 I=1,Niz
       IM1=I-1
       IF(IM1.EQ.0) IM1=Niz
       TXM(I,Njz-1)=(TXN(IM1,2)+TXN(I,2))*0.5
       TYM(I,Njz-1)=(TYN(IM1,2)+TYN(I,2))*0.5
       TXM(I,1)=(TXN(IM1,Njz)+TXN(I,Njz))*0.5
       TYM(I,1)=(TYN(IM1,Njz)+TYN(I,Njz))*0.5
 3     CONTINUE
 4     RETURN
       END
