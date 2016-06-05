      Subroutine GlobAcc(lun)
c       Writes Global values to File
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      common/glacc/GMTACC,GMKEAC,GMRACC,namgl 
*********** END OF COMMON ******************

C        GLOBAL T, K.E. ,N0 - DAY MEANS - to file GLOBAL
         GMTACC=GMTACC/namgl  
         GMKEAC=GMKEAC/namgl
	   ! gmracc (w/m**2)  
         GMRACC=GMRACC/namgl*0.484 !0.484 W/m**2= 1 ly/day 
	! To GLOBAL File:
         WRITE (lun,1) TAU,MONTH,GMKEAC,GMTACC,GMRACC,SDEDY,SDEYR,id
 1       FORMAT (1X,F9.2,2X,I3,2X,3f9.4,2x,i3,2x,f5.0,2x,a4)
      GMTACC=0.
      GMKEAC=0.
      GMRACC=0.
	namgl=0
      return
      end

