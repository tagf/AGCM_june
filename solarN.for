C     ****************                                                  
C     ****************                                                  
      SUBROUTINE SOLAR(J1)                                              
C     ****************                                                  
C     ****************                                                  
C                                                                       
C   ++++ D P ++++                                                       
C             CALCULATES SHORT WAVE RADIATION BUDGET                    
C             NEW VERSION                                               
C                                                                       
C************ BEGINNING OF COMMON ************                          
	include 'recom.fi'
C                                                                       
      LOGICAL SNOW,SEAICE                                               
C                                                                       
      include 'radvar.fi'
C                                                                       
      COMMON /MOIST/ EFVT(74),EFVST(74),EFV3(74),EFV2(74),              
     1               EFV1(74),EFV0(74),AK(74)                           
C                                                                       
	include 'cover.fi'
C                                                                       
      LOGICAL CIRUS,CIRUSL       !where  used?                                       
      COMMON/NEWSOL/ CIRUS(74),CIRUSL(74),GREY(74),GREYL(74),CLC(74)    
C                                                                       
      COMMON/DSDD/ DD(74,46,5),DS(74,46,5)                              
C                                                                       
      COMMON /WORK/
     1  SAO0(74),SAO2(74),SAO4(74),SAC0(74),SAC2(74),SAC4(74),
     2  AM  (74),SA  (74),SECZ(74),ALSF(74),CAO0(74),CAO2(74),
     3  CAO4(74),CCO0(74),CCO2(74),CCO4(74),SACTL(74),SACTU(74),
     4  CCTL(74),CCTU(74),EFVCU(74),EFVCL(74),YCU(74),YCL(74),
     5  ALA0(74),ALAC(74),TEMA(74),TEMB(74),TEMC(74),TC4L(74),
     6  TO0 (74),TO1 (74),TO2 (74),TO3 (74),TO4 (74),TC0 (74),
     7  TC2 (74),TC4 (74),TWT0(74),TWT1(74),TWT2(74),TWT3(74),
     8  TWT4(74),TWC (74),TOC (74),TWCL(74),TOCL(74),TWC2(74),
     9  TWC4(74),TWC4L(74),RABAR(74),RSBAR(74),ALACL(74),SS(74),
     X  SST(74),SAC4L(74),SAG(74),SSG(74),EV01(74),EV12(74),EV23(74),
     X  RACTL(74)
C                                                                       
      COMMON /COMB/ RU,TU,R1,TT1,RL,TL,RR2,TT2,TS1,RC,TCC,R3,TT3,XS     
C                                                                       
C************ END OF COMMON ******************                          
C                                                                       
      DIMENSION D(11,5,6),G(11,5,6),W(11,5,6),DL(5,6),R(5,6),           
     1        T(5,6),RB(5,6),TB(5,6),CLD(10),WL(5,6),GL(5,6)            
      DIMENSION TS(5,6),XIN(6)                                          
      DIMENSION XX(5),YY(5),RD(6),TD(6),RR4(6),TT4(6),RH(6),TH(6),      
     1  R8(6),T8(6),RI(6),TI(6),R9(6),T9(6),RF(6),TF(6),R6(6),T6(6),    
     2  RJ(6),TJ(6),RX(6),TX(6),RG(6),TTG(6),R7(6),T7(6),RK(6),         
     3  TK(6),RZ(6),TZ(6),S0Z(6),S2(6),SS4(6),RE(6),TE(6),R5(6),T5(6),  
     4  RET(6)                                                          
C                                                                       
      EXP(XXX)=DEXP(XXX)                                                
      SQRT(XXX)=DSQRT(XXX)                                              
C                                                                       
C     DD=0.                                                             
C     DS=0.                                                             
C                                                                       
      DO 1 I=1,74                                                       
      S4(I)=0.0                                                         
      AS1(I)=0.0                                                        
      AS3(I)=0.0                                                        
      RETOT(I)=0.0                                                      
    1 CONTINUE                                                          
      DO 17 I=2,73                                                      
      TEMA(I)=ALS(I)                                                    
      TEMB(I)=0.0                                                       
   17 TEMC(I)=0.0                                                       
C                                                                       
      DO 20 I=2,73                                                      
      IF (ISRFCE(I).NE.7) GO TO 20                                      
      TEMA(I)=0.9961572                                                 
      TEMB(I)=5.4077056                                                 
      TEMC(I)=9.2530638                                                 
      IF (COSZ(I).LE.0.258819) GO TO 20                                 
      TEMA(I)=0.4655505                                                 
      TEMB(I)=1.2820763                                                 
      TEMC(I)=0.940257                                                  
      IF (COSZ(I).LT.0.707107) GO TO 20                                 
      TEMA(I)=0.0371737                                                 
      TEMB(I)=0.0172538                                                 
      TEMC(I)=0.0                                                       
   20 CONTINUE                                                          
C                                                                       
      DO 30 I=2,73                                                      
   30 ALSF(I)=TEMA(I)-COSZ(I)*(TEMB(I)-COSZ(I)*TEMC(I))                 
C     DO 400 I=2,25                                                     
C     AM(I)=35.0                                                        
C     IF (COSZ(I).LE.0.01) GO TO 400                                    
C     SECZ(I)=1.0/COSZ(I)                                               
C     AM(I)=35.*SECZ(I)/SQRT(1224.+SECZ(I)*SECZ(I))                     
C     ENDIF                                                             
C 400 CONTINUE                                                          
C                                                                       
         DO 71 M=1,6                                                    
         DO 71 J=1,5                                                    
         DO 71 I=1,11  !?????????????                                                 
         W(I,J,M)=0.                                                    
         G(I,J,M)=0.                                                    
 71      D(I,J,M)=0.                                                    
C                                                                       
      DO 1777 II=2,73                                                   
      IF (COSZ(II).LE.0.01) GO TO 1777                                  
      AS=ALSF(II)                                                       
      SC=S0*.4847593                                                    
      XS=COSZ(II)                                                       
      X1=XS                                                             
      XS=SQRT(1224.*XS*XS+1.)/35.                                       
      CLD(3)=CL1(II)                                                    
      CLD(4)=CL2(II)                                                    
      IF(CIRUSL(II).OR.(CIRUS(II).AND.CL4(II).EQ.0.0))CLD(4)=0.0        
      CLD(5)=CL2(II)-CLD(4)                                             
  500 CLD(6)=CL3(II)                                                    
      IF(CIRUSL(II).OR.(CIRUS(II).AND.CL4(II).EQ.0.0))CLD(6)=0.0        
      CLD(7)=CL3(II)-CLD(6)                                             
  501 CLD(8)=CL4(II)                                                    
      IF (CIRUS(II)) CLD(8)=0.0                                         
      CLD(9)=CL4(II)-CLD(8)                                             
      XX(1)=EFVT(II)-EFV0(II)                                           
      XX(2)=EFV0(II)-EFV1(II)                                           
      XX(3)=EFV1(II)-EFV2(II)                                           
      XX(4)=EFV2(II)-EFV3(II)                                           
      XX(5)=EFV3(II)                                                    
C        DO 71 M=1,6                                                    
C        DO 71 J=1,5                                                    
C        DO 71 I=1,11                                                   
C        W(I,J,M)=0.                                                    
C        G(I,J,M)=0.                                                    
C 71     D(I,J,M)=0.                                                    
      DO 69 J=1,5                                                       
      YY(J)=1.0                                                         
 69   CONTINUE                                                          
      IF (CLD(3).EQ.1.0) YY(2)=0.                                       
      IF (CLD(8).EQ.1.0) YY(3)=0.                                       
      IF (CLD(9).EQ.1.0) YY(3)=0.                                       
      IF (CLD(4).EQ.1.0) YY(4)=0.                                       
      IF (CLD(5).EQ.1.0) YY(4)=0.                                       
      IF (CLD(6).EQ.1.0) YY(5)=0.                                       
      IF (CLD(7).EQ.1.0) YY(5)=0.                                       
      DO 440 J=1,5                                                      
      D(2,J,1)=1.093*DD(II,J1,J)                                        
      DO 2000 M=2,6                                                     
      D(2,J,M)=.211*DD(II,J1,J)                                         
 2000 CONTINUE                                                          
 440  CONTINUE                                                          
      DO 220 J=1,5                                                      
      D(1,J,1)=.0262                                                    
      DO 2001 M=2,6                                                     
      D(1,J,M)=.0003                                                    
 2001 CONTINUE                                                          
      D(11,J,1)=0.0                                                     
      D(11,J,2)=.005*YY(J)*XX(J)                                        
      D(11,J,3)=.041*YY(J)*XX(J)                                        
      D(11,J,4)=.416*YY(J)*XX(J)                                        
      D(11,J,5)=4.75*YY(J)*XX(J)                                        
      D(11,J,6)=72.5*YY(J)*XX(J)                                        
 220  CONTINUE                                                          
      DO 330 J=1,5                                                      
      D(10,J,1)=.27*DS(II,J1,J)                                         
      DO 2002 M=2,6                                                     
      D(10,J,M)=.071*DS(II,J1,J)                                        
 2002 CONTINUE                                                          
 330  CONTINUE                                                          
C                                                                       
      DO 2003 M=1,6                                                     
      D(3,2,M)=CLD(3)*12.                                               
      D(4,4,M)=CLD(4)*7.                                                
      D(5,4,M)=CLD(5)*2.                                                
      D(6,5,M)=CLD(6)*12.                                               
      D(7,5,M)=CLD(7)*2.                                                
      D(8,3,M)=CLD(8)*7.                                                
      D(9,3,M)=CLD(9)*2.                                                
 2003 CONTINUE                                                          
C                                                                       
      G(2,1,1)=.632                                                     
      DO 2004 M=2,6                                                     
      G(2,1,M)=.573                                                     
 2004 CONTINUE                                                          
      DO 35 J=1,5                                                       
      G(1,J,1)=0.0                                                      
      DO 2005 M=2,6                                                     
      G(1,J,M)=0.0                                                      
 2005 CONTINUE                                                          
      DO 2006 M=1,6                                                     
      G(11,J,M)=0.0                                                     
 2006 CONTINUE                                                          
  35  CONTINUE                                                          
      DO 40 J=2,5                                                       
      G(10,J,1)=.706                                                    
      DO 2007 M=2,6                                                     
      G(10,J,M)=.659                                                    
 2007 CONTINUE                                                          
  40  CONTINUE                                                          
      G(3,2,1)=CLD(3)*.84                                               
      G(4,4,1)=CLD(4)*.84                                               
      G(5,4,1)=CLD(5)*.84                                               
      G(6,5,1)=CLD(6)*.84                                               
      G(7,5,1)=CLD(7)*.84                                               
      G(8,3,1)=CLD(8)*.84                                               
      G(9,3,1)=CLD(9)*.84                                               
C                                                                       
      DO 2008 M=2,6                                                     
      G(3,2,M)=CLD(3)*.6                                                
      G(4,4,M)=CLD(4)*.76                                               
      G(5,4,M)=CLD(5)*.82                                               
      G(6,5,M)=CLD(6)*.6                                                
      G(7,5,M)=CLD(7)*.82                                               
      G(8,3,M)=CLD(8)*.76                                               
      G(9,3,M)=CLD(9)*.82                                               
 2008 CONTINUE                                                          
C                                                                       
      DO 2009 M=1,6                                                     
      W(2,1,M)=1.                                                       
 2009 CONTINUE                                                          
      DO 50 J=1,5                                                       
      DO 2010 M=1,6                                                     
      W(1,J,M)=1.0                                                      
 2010 CONTINUE                                                          
      DO 2011 M=1,6                                                     
      W(11,J,M)=0.0                                                     
 2011 CONTINUE                                                          
  50  CONTINUE                                                          
      DO 60 J=2,5                                                       
      W(10,J,1)=.705                                                    
      DO 2012 M=2,6                                                     
      W(10,J,M)=.552                                                    
 2012 CONTINUE                                                          
 60   CONTINUE                                                          
      W(3,2,1)=CLD(3)                                                   
      W(4,4,1)=CLD(4)                                                   
      W(5,4,1)=CLD(5)                                                   
      W(6,5,1)=CLD(6)                                                   
      W(7,5,1)=CLD(7)                                                   
      W(8,3,1)=CLD(8)                                                   
      W(9,3,1)=CLD(9)                                                   
      DO 2013 M=2,6                                                     
      W(3,2,M)=CLD(3)*.982                                              
      W(4,4,M)=CLD(4)*.983                                              
      W(5,4,M)=CLD(5)*.99                                               
      W(6,5,M)=CLD(6)*.982                                              
      W(7,5,M)=CLD(7)*.99                                               
      W(8,3,M)=CLD(8)*.983                                              
      W(9,3,M)=CLD(9)*.99                                               
 2013 CONTINUE                                                          
      A=0.0                                                             
      DO 70 M=1,6                                                       
      DO 80 J=1,5                                                       
      DO 90 I=1,11                                                      
      A=A+D(I,J,M)                                                      
  90  CONTINUE                                                          
      DL(J,M)=A                                                         
      A=0.0                                                             
 80   CONTINUE                                                          
 70   CONTINUE                                                          
      DO 100 M=1,6                                                      
      DO 110 J=1,5                                                      
      DO 120 I=1,11                                                     
      IF (I.EQ.1) GG=SQRT(0.2D0)                                        
      IF (I.GT.1) GG=G(I,J,M)                                           
      A=A+(1.-W(I,J,M)*GG*GG)*D(I,J,M)                                  
 120  CONTINUE                                                          
      TS(J,M)=A                                                         
      A=0.0                                                             
 110  CONTINUE                                                          
 100  CONTINUE                                                          
C                                                                       
      A=0.0                                                             
      DO 130 M=1,6                                                      
      DO 140 J=1,5                                                      
      DO 150 I=1,11                                                     
      A=A+W(I,J,M)*D(I,J,M)                                             
 150  CONTINUE                                                          
      WL(J,M)=A/DL(J,M)                                                 
      A=0.0                                                             
 140  CONTINUE                                                          
 130  CONTINUE                                                          
C                                                                       
      A=0.0                                                             
      DO 160 M=1,6                                                      
      DO 170 J=1,5                                                      
      DO 180 I=1,11                                                     
      A=A+D(I,J,M)*W(I,J,M)*G(I,J,M)                                    
 180  CONTINUE                                                          
      GL(J,M)=A/(WL(J,M)*DL(J,M))                                       
      A=0.0                                                             
  170 CONTINUE                                                          
  160 CONTINUE                                                          
C                                                                       
      DO 190 M=1,6                                                      
      DO 200 J=1,5                                                      
      GG=GL(J,M)                                                        
      WW=WL(J,M)                                                        
      DDD=DL(J,M)                                                       
C     XS1=XS                                                            
      CALL REFTRA(GG,WW,DDD,R1A,T1A,RB1A,TB1A,XS)                       
C                                                                       
      R(J,M)=R1A                                                        
      T(J,M)=T1A                                                        
      RB(J,M)=RB1A                                                      
      TB(J,M)=TB1A                                                      
 200  CONTINUE                                                          
 190  CONTINUE                                                          
C                                                                       
      DO 201 M=1,6                                                      
      TS1=TS(1,M)                                                       
      RU=R(1,M)                                                         
      TU=T(1,M)                                                         
      R1=RB(1,M)                                                        
      TT1=TB(1,M)                                                       
      RL=R(2,M)                                                         
      TL=T(2,M)                                                         
      RR2=RB(2,M)                                                       
      TT2=TB(2,M)                                                       
      CALL COMBIN                                                       
      TS1=TS(1,M)+TS(2,M)                                               
      RU=RC                                                             
      TU=TCC                                                            
      R1=R3                                                             
      TT1=TT3                                                           
      RL=R(3,M)                                                         
      TL=T(3,M)                                                         
      RR2=RB(3,M)                                                       
      TT2=TB(3,M)                                                       
      CALL COMBIN                                                       
      RD(M)=RC                                                          
      TD(M)=TCC                                                         
      RR4(M)=R3                                                         
      TT4(M)=TT3                                                        
 201  CONTINUE                                                          
      DO 202 M=1,6                                                      
      TS1=TS(1,M)+TS(2,M)+TS(3,M)                                       
      RU=RD(M)                                                          
      TU=TD(M)                                                          
      R1=RR4(M)                                                         
      TT1=TT4(M)                                                        
      RL=R(4,M)                                                         
      TL=T(4,M)                                                         
      RR2=RB(4,M)                                                       
      TT2=TB(4,M)                                                       
      CALL COMBIN                                                       
      TS1=TS(1,M)+TS(2,M)+TS(3,M)+TS(4,M)                               
      RU=RC                                                             
      TU=TCC                                                            
      R1=R3                                                             
      TT1=TT3                                                           
      RL=R(5,M)                                                         
      TL=T(5,M)                                                         
      RR2=RB(5,M)                                                       
      TT2=TB(5,M)                                                       
      CALL COMBIN                                                       
      RE(M)=RC                                                          
      TE(M)=TCC                                                         
      R5(M)=R3                                                          
      T5(M)=TT3                                                         
  202 CONTINUE                                                          
C                                                                       
      DO 203 M=1,6                                                      
      TS1=TS(4,M)                                                       
      RU=R(4,M)                                                         
      TU=T(4,M)                                                         
      R1=RB(4,M)                                                        
      TT1=TB(4,M)                                                       
      RL=R(5,M)                                                         
      TL=T(5,M)                                                         
      RR2=RB(5,M)                                                       
      TT2=TB(5,M)                                                       
      CALL COMBIN                                                       
      RH(M)=RC                                                          
      TH(M)=TCC                                                         
      R8(M)=R3                                                          
      T8(M)=TT3                                                         
 203  CONTINUE                                                          
C                                                                       
      DO 204 M=1,6                                                      
      TS1=TS(2,M)                                                       
      RU=R(2,M)                                                         
      TU=T(2,M)                                                         
      R1=RB(2,M)                                                        
      TT1=TB(2,M)                                                       
      RL=R(3,M)                                                         
      TL=T(3,M)                                                         
      RR2=RB(3,M)                                                       
      TT2=TB(3,M)                                                       
      CALL COMBIN                                                       
      TS1=TS(2,M)+TS(3,M)                                               
      RU=RC                                                             
      TU=TCC                                                            
      R1=R3                                                             
      TT1=TT3                                                           
      RL=RH(M)                                                          
      TL=TH(M)                                                          
      RR2=R8(M)                                                         
      TT2=T8(M)                                                         
      CALL COMBIN                                                       
      RI(M)=RC                                                          
      TI(M)=TCC                                                         
      R9(M)=R3                                                          
      T9(M)=TT3                                                         
 204  CONTINUE                                                          
C                                                                       
      DO 205 M=1,6                                                      
      TS1=TS(5,M)                                                       
      RU=R(5,M)                                                         
      TU=T(5,M)                                                         
      R1=RB(5,M)                                                        
      TT1=TB(5,M)                                                       
      RL=R(4,M)                                                         
      TL=T(4,M)                                                         
      RR2=RB(4,M)                                                       
      TT2=TB(4,M)                                                       
      CALL COMBIN                                                       
      RF(M)=RC                                                          
      TF(M)=TCC                                                         
      R6(M)=R3                                                          
      T6(M)=TT3                                                         
 205  CONTINUE                                                          
C                                                                       
      DO 210 M=1,6                                                      
      TS1=TS(3,M)                                                       
      RU=R(3,M)                                                         
      TU=T(3,M)                                                         
      R1=RB(3,M)                                                        
      TT1=TB(3,M)                                                       
      RL=R(2,M)                                                         
      TL=T(2,M)                                                         
      RR2=RB(2,M)                                                       
      TT2=TB(2,M)                                                       
      CALL COMBIN                                                       
      TS1=TS(3,M)+TS(2,M)                                               
      RU=RC                                                             
      TU=TCC                                                            
      R1=R3                                                             
      TT1=TT3                                                           
      RL=R(1,M)                                                         
      TL=T(1,M)                                                         
      RR2=RB(1,M)                                                       
      TT2=TB(1,M)                                                       
      CALL COMBIN                                                       
      RJ(M)=RC                                                          
      TJ(M)=TCC                                                         
      RX(M)=R3                                                          
      TX(M)=TT3                                                         
 210  CONTINUE                                                          
C                                                                       
      DO 2220 M=1,6                                                     
      TS1=TS(5,M)+TS(4,M)                                               
      RU=RF(M)                                                          
      TU=TF(M)                                                          
      R1=R6(M)                                                          
      TT1=T6(M)                                                         
      RL=R(3,M)                                                         
      TL=T(3,M)                                                         
      RR2=RB(3,M)                                                       
      TT2=TB(3,M)                                                       
      CALL COMBIN                                                       
      TS1=TS(5,M)+TS(4,M)+TS(3,M)                                       
      RU=RC                                                             
      TU=TCC                                                            
      R1=R3                                                             
      TT1=TT3                                                           
      RL=R(2,M)                                                         
      TL=T(2,M)                                                         
      RR2=RB(2,M)                                                       
      TT2=TB(2,M)                                                       
      CALL COMBIN                                                       
      RG(M)=RC                                                          
      TTG(M)=TCC                                                        
      R7(M)=R3                                                          
      T7(M)=TT3                                                         
 2220 CONTINUE                                                          
C                                                                       
      DO 230 M=1,6                                                      
      TS1=TS(5,M)+TS(4,M)+TS(3,M)+TS(2,M)                               
      RU=RG(M)                                                          
      TU=TTG(M)                                                         
      R1=R7(M)                                                          
      TT1=T7(M)                                                         
      RL=R(1,M)                                                         
      TL=T(1,M)                                                         
      RR2=RB(1,M)                                                       
      TT2=TB(1,M)                                                       
      CALL COMBIN                                                       
      RK(M)=RC                                                          
      TK(M)=TCC                                                         
      RZ(M)=R3                                                          
      TZ(M)=TT3                                                         
 230  CONTINUE                                                          
C                                                                       
      XIN(1)=SC*X1*.647                                                 
      XIN(2)=SC*X1*.107                                                 
      XIN(3)=SC*X1*.104                                                 
      XIN(4)=SC*X1*.073                                                 
      XIN(5)=SC*X1*.044                                                 
      XIN(6)=SC*X1*.025                                                 
      DO 240 M=1,6                                                      
      AT=RE(M)+AS*TE(M)*T5(M)/(1.-AS*RZ(M))                             
      RET(M)=XIN(M)*AT                                                  
      A0=RI(M)+AS*TI(M)*T9(M)/(1.-AS*R7(M))                             
      A1=R9(M)+AS*T9(M)*T9(M)/(1.-AS*R7(M))                             
C                                                                       
      TS1=TS(1,M)/XS                                                    
      IF (TS1.GT.40.)TS1=40.                                            
      TS1EXP=EXP(-TS1)                                                  
C                                                                       
      DTT=TS1EXP                                                        
      DF=(TS1EXP*RB(1,M)*A0+T(1,M)-TS1EXP)/                             
     1   (1.-A1*RB(1,M))                                                
      S0Z(M)=XIN(M)*((1.-A0)*DTT+(1.-A1)*DF)                            
      A2=RH(M)+AS*TH(M)*T8(M)/(1.-AS*R6(M))                             
      A3=R8(M)+AS*T8(M)*T8(M)/(1.-AS*R6(M))                             
C                                                                       
      TS1=(TS(1,M)+TS(2,M)+TS(3,M))/XS                                  
      IF (TS1.GT.40.) TS1=40.                                           
      TS1EXP=EXP(-TS1)                                                  
C                                                                       
      DTT=TS1EXP                                                        
      DF=(TS1EXP*RX(M)*A2+TD(M)-TS1EXP)/(1.-A3*RX(M))                   
      S2(M)=XIN(M)*((1.-A2)*DTT+(1.-A3)*DF)                             
      SS4(M)=XIN(M)*(1.-AS)*TE(M)/(1.-RZ(M)*AS)                         
 240  CONTINUE                                                          
      RETO=RET(1)+RET(2)+RET(3)+RET(4)+RET(5)+RET(6)                    
      S0A=S0Z(1)+S0Z(2)+S0Z(3)+S0Z(4)+S0Z(5)+S0Z(6)                     
      S2A=S2(1)+S2(2)+S2(3)+S2(4)+S2(5)+S2(6)                           
      SS4A=SS4(1)+SS4(2)+SS4(3)+SS4(4)+SS4(5)+SS4(6)                    
C                                                                       
      S4(II)=SS4A                                                       
      AS1(II)=S0A-S2A                                                   
      AS3(II)=S2A-SS4A                                                  
      RETOT(II)=RETO                                                    
C     PRINT 2,II                                                        
      RETOT(II)=RETOT(II)*2.0628796                                     
      S4(II)=S4(II)*2.0628796                                           
      AS1(II)=AS1(II)*2.0628796                                         
      AS3(II)=AS3(II)*2.0628796                                         
C                                                                       
C     PRINT 3,SS4A,SS4(1),SS4(2),SS4(3),SS4(4),SS4(5),SS4(6)            
C    1XX(1),XX(2),XX(3),XX(4),XX(5),S4(II),AS1(II),AS3(II),RETOT(II)    
 3     FORMAT(1X,'S4=',D11.3,'SS4=',6(1X,D11.3))                        
 2     FORMAT(1X,'II=',I2)                                              
      GO TO 1696                                                        
 1696 CONTINUE                                                          
 1777 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE REFTRA(G,W,D,R,T,RB,TB,X)                              
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)                             
      EXP(XXX)=DEXP(XXX)                                                
      SQRT(XXX)=DSQRT(XXX)                                              
      IF (G.EQ.0.0) F=.2                                                
      IF (G.GT.0.0) F=G*G                                               
      TS=(1.-W*F)*D                                                     
      WS=(1.-F)*W/(1.-W*F)                                              
      GS=(G-F)/(1.-F)                                                   
  1   IF (W.EQ.1.) GO TO 341                                            
C                                                                       
      XLA=SQRT(3.*(1.-WS)*(1.-WS*GS))                                   
      U=1.5*(1.-WS*GS)/XLA                                              
      ZZ=XLA*TS                                                         
      IF(ZZ.GT.40.) ZZ=40.                                              
C                                                                       
2     ZEXP=EXP(ZZ)                                                      
3     ZZEXP=1./ZEXP                                                     
4        TZZ=TS/X                                                       
         IF (TZZ.GT.40.) TZZ=40.                                        
5        TSEXP=EXP(-TZZ)                                                
C                                                                       
      XN=(U+1.)*(U+1.)*ZEXP-(U-1.)*(U-1.)*ZZEXP                         
6     TB=4.*U/XN                                                        
7     RB=(U+1.)*(U-1.)*(ZEXP-ZZEXP)/XN                                  
8     BE=.5*WS*(3.*GS*(1.-WS)*X*X+1.)/(1.-XLA*XLA*X*X)                  
9     AL=.75*WS*X*(1.+GS*(1.-WS))/(1.-XLA*XLA*X*X)                      
10    R=(AL-BE)*TB*TSEXP+(AL+BE)*RB-AL+BE                               
11    T=(AL+BE)*TB+(AL-BE)*TSEXP*RB-(AL+BE)*TSEXP+TSEXP                 
C                                                                       
      IF(R.LT.0.0) T=TSEXP                                              
      IF(R.LT.0.0) R=0.0                                                
C                                                                       
      IF (W.LT.1.) GO TO 344                                            
C                                                                       
 341  TZZ=TS/X                                                          
      IF (TZZ.GT.40.) TZZ=40.                                           
12    TSEXP=EXP(-TZZ)                                                   
13    R=((1.-GS)*TS+0.66666666*(1.-1.5*X)*(1.-TSEXP))/(1.33333333       
     1  +(1.-GS)*TS)                                                    
      T=1.-R                                                            
14    RB=D*(1.-G)/(1.33333333+D*(1.-G))                                 
      TB=1.-RB                                                          
 344  CONTINUE                                                          
         IF (TB.LE.1.d-10) TB=0.                                        
      RETURN                                                            
      END                                                               
C                                                                       
       SUBROUTINE COMBIN                                                
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)                             
      COMMON/COMB/ RU,TU,R1,TT1,RL,TL,RR2,TT2,TS1,RC,TCC,R3,TT3,XS      
      EXP(XXX)=DEXP(XXX)                                                
      ZZ=TS1/XS                                                         
      IF (ZZ.GT.40.) ZZ=40.                                             
      ZZEXP=EXP(-ZZ)                                                    
      R1RR2=1.0/(1.0-R1*RR2)                                            
C                                                                       
      RC=RU+TT1*(RL*ZZEXP+RR2*(TU-ZZEXP))*R1RR2                         
      TCC=TL*ZZEXP+TT2*((TU-ZZEXP)+R1*RL*ZZEXP)*R1RR2                   
C                                                                       
      R3=R1+TT1*TT1*RR2*R1RR2                                           
      TT3=TT1*TT2*R1RR2                                                 
      RETURN                                                            
      END                                                               
