     
      subroutine SurFile(title,arr,zero,fileN)

c     To get special file for SURFER from arr(74,46) 
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))
	include 'dir.fi'
C************ END OF COMMON ******************

       real*4  arr1(74,46) !real*4 is need to print E - format (not D!)
	 DIMENSION arr(74,46)
	 integer*4 PointDate(6)
       character*(*) fileN, title

           ! To transform data - Zero point and single: 
         arr1=arr-zero

       open (12, file=TRIM(BaseDir)//WorkDir//'\'//
     1	 trim(fileN)//'.dat', status='new',err=2)    
       goto 3
  2    STOP '  ERROR: SB SurFile - File name Duplicate.'

  3	PointDate(1)=int(tau)
	PointDate(2)=month
	PointDate(3)=mnthdy
	PointDate(4)=int(tofday)
	PointDate(5)=0 !ID
	PointDate(6)=sdedy

      write (12,10) TRIM(title),PointDate 
 10   FORMAT(A,1X,I7,1Hh,2X,I3,1H/,I2,1H/,i2,3H.00,2X,
     1   4H ID=,A4,2X,' day=',I3,15H* * * * * * * *)
      DO 4 I=2,73
      DO 4 J=1,46
 4     write (12,11) i-1,j,arr1(i,j)
 11    FORMAT (1x,2i3,1x,e12.4)
      close (12)

       return
       end
