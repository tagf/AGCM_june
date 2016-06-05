C     *****************                                                 
C     *****************                                                 
      SUBROUTINE ACC0                                                   
C     *****************                                                 
C     *****************                                                 
C                                                                       
C                                                                       
C             ZERO ACCUMULATED VARIABLES                                
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
     1   ,SD(74,46),PIV(74,46,2)
C
C        ACCUMULATED VARIABLES IN COMMON
C
	include 'accum.fi'
      include 'acc_ice.fi'
      common/glacc/GMTACC,GMKEAC,GMRACC,namgl 
C************ END OF COMMON ******************                          
         NAV=0                                                          
C        SET ACCUMULATIONS TO ZERO
C
      ZERO=0.0
      DO 15 N=63,74,2
      DO 15 J=1,46
      ZONAVG(J,N)=-1.E20
   15 ZONAVG(J,N+1)=1.E20
      ZERO1=0.0
      zeroIce=0.0                                              
C***   GLOBAL T,K.E.,N0                                                 
         GMTACC=0.                                                      
         GMKEAC=0.                                                      
         GMRACC=0.                                                      
C***                                                                    
         print 111,tau
 111     FORMAT (' ********** ZERO ACCUMULATED DATA      TAU=',F9.2)
      RETURN                                                            
      END                                                               
