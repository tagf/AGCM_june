C     *****************
C     *****************
      SUBROUTINE OutAcIce(lun,nappend)
C     *****************
C     *****************
C
C
C             OUTPUTS ACCUMULATED ICE VARIABLES HISTORY
C                      
c         nappend=0   rewind
c         nappend=1   append
C			write *4 data - disk resourses economy
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
C        ACCUMULATED Ice VARIABLES IN COMMON
C
      include 'acc_ice.fi'
      include 'ice.fi'

	 include 'dir.fi'
C
C************ END OF COMMON ******************
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
      IF (MONTH.LT.1.OR.MONTH.GT.12) WRITE (6,9) MONTH,TAU
    9 FORMAT (14H0ERR- MONTH = ,I8,7H TAU = ,F10.2)
C
      if (nappend.eq.0) then 
	 OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\AccIce.dat',
     1	 FORM='BINARY')
	 REWIND LUN
        write (lun) sngl(hIceAcc)
        write (lun) sngl(GIceCompAcc)
       CLOSE (LUN)
	else   !nappend.eq.1
	 OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\AccIce1.dat',  
     1      ACCESS='append',FORM='BINARY',status='unknown'
     2      ,err=1)
	 goto 2
1       write(*,'(a)') ' open AccIce1 error - ok '
	  OPEN (LUN,FILE=TRIM(BaseDir)//WorkDir//'\AccIce1.dat',  
     1      ACCESS='SEQUENTIAL',FORM='BINARY',status='unknown')
c         North pole point data correction	
2      a_sredneeJM = 0.0
      do I=1,IM
       a_sredneeJM = a_sredneeJM + hIceAcc(I,JM-1)
	enddo
      do I=1,IM
       hIceAcc(I,JM)   = a_sredneeJM / IM
	enddo
      a_sredneeJM = 0.0
      do I=1,IM
       a_sredneeJM = a_sredneeJM + GIceCompAcc(I,JM-1)
	enddo
      do I=1,IM
       GIceCompAcc(I,JM)   = a_sredneeJM / IM
	enddo

        write (lun) sngl(hIceAcc/nav)
        write (lun) sngl(GIceCompAcc/nav)
         CLOSE (LUN)
      endif
c      write (lun) ((ARR(I,J),I=1,IM),J=1,JM)
c         write (7,111) TAU,MONTH,MNTHDY,TOFDAY ! To PROTOCOL FILE
c         print 111,tau,MONTH,MNTHDY,TOFDAY
c 111     FORMAT (' ****  OUTACC RECORDS    TAU=',F9.2,
c     1           1X,I3,1H/,I2,1H/,F6.2)
      RETURN
      END
