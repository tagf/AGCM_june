C     *****************
C     *****************
      SUBROUTINE  InAccIce(lun)
C     *****************
C     *****************
C
C         INPUT MONTHLY   ACCUMULATED Ice VARIABLES HISTORY
C		read *4 data and convert to *8 data
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
C        ACCUMULATED VARIABLES IN COMMON
C
      include 'acc_ice.fi'
C
	include 'dir.fi'
C************ END OF COMMON ******************
C
	 character*40 infile
	real*4 tt(72,46) ! temporary *4 array
C
        infile=TRIM(BaseDir)//WorkDir//'\AccIce.dat'
	 OPEN (LUN,FILE=infile,FORM='BINARY')
      rewind lun
      read (lun) tt
	 hIceAcc=DBLE(tt)
      read (lun) tt
 	 GIceCompAcc=DBLE(tt)
	close (LUN)
c       
        PRINT 111,TRIM(infile),TAU
 111     FORMAT (' INPUT ACCUM DATA. FILE ',a,' TAU=',F9.2)
      RETURN
      END
