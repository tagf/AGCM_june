C     *****************
C     *****************
      SUBROUTINE  GradsCtl(lun)
C     *****************
C     *****************
C
C         INPUT initial day, month, year to file AccIce1.ctl
C   	                      for Grads
C
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
	include 'dir.fi'
C************ END OF COMMON ******************
C
	 character*50 infile
	 character*3 MonTitl(12)/'jan','feb','mar','apr','may','jun',
     1 	     'jul','aug','sep','oct','nov','dec'/
C
        infile=TRIM(BaseDir)//WorkDir//'\AccIce1.ctl'
	 OPEN (LUN,FILE=infile,FORM='formatted') 
       read (lun,2)
 2      FORMAT (8/)

      WRITE(lun,4)ID,SDEDY,SDEYR,DT  !data identification
 4    FORMAT ('TITLE ID = ',A4,' *SDEDY=',I7,
     1        ' *SDEYR =',F8.2,' *DT=',F8.2)
c       
        write (lun,1) MNTHDY,MonTitl(month),nint(SDEYR) !initial time point
 1      FORMAT ('TDEF 600 LINEAR ',i2.2,a3,'00',i2.2,' 1mo')

	 close (LUN)
      RETURN
      END
