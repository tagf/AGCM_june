C *********************************************************************
C             To rewrite SST  SS(72,46,12)  UNFORMATTED
C                      
C *********************************************************************
      DIMENSION SS(72,46,12)
      DO 1 NAV=1,12
       NAV1=NAV+40
      READ (NAV1,31) I
    1 READ (NAV1,30) ((SS(I,J,NAV),I=1,72),J=1,46)
   30 FORMAT(4E20.13)
   31 FORMAT(A4)
      DO 2 K=1,12
    2 WRITE (54) ((SS(I,J,K),I=1,72),J=1,46)
      STOP
      END                                                                       
