c********************************************* 
c     вывод в формате GRADS
c********************************************* 
      SUBROUTINE OutGRADS1(arr1,im1,jm1,shift,scale,filename)
C***** IM1=72 OR IM1=74  JM1=46    ******
C***** PRINT OUT TRANSFORMED ARRAY: APRINT=(A-SHIFT)/SCALE ******
C                                                                       
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))
C
	include 'dir.fi'
c
C****************** END OF COMMON ***********************

      DIMENSION ARR1(IM1,JM1)
      character*(*) FILENAME
      REAL ARR(IM,JM)

c  подготовка данных
      do J=1,JM
       do I=1,IM
	  if (IM1.eq.74) then
         ARR(I,J) = (ARR1(I+1,J) - SHIFT) / SCALE
	  else ! IM1.eq.72
         ARR(I,J) = (ARR1(I,J) - SHIFT) / SCALE
	  endif
	 enddo
	enddo

      a_sredneeJM = 0.0
      a_srednee1  = 0.0
      do I=1,IM
       a_sredneeJM = a_sredneeJM + ARR(I,JM-1)
       a_srednee1  = a_srednee1 + ARR(I,2)
	enddo
      do I=1,IM
       ARR(I,JM)   = a_sredneeJM / IM
       ARR(I,1)    = a_srednee1 / IM
	enddo

c      lun1=30 ! файл описания
c	OPEN (lun1,FILE=BaseDir//WorkDir//'\'//FILENAME//'.ctl')
      lun2=31 ! файл данных
	OPEN (lun2,FILE=TRIM(BaseDir)//WorkDir//'\'//FILENAME//'.dat',
     1      FORM='BINARY',ACCESS='APPEND')

      write (lun2) ARR

c      CLOSE(lun1)
      CLOSE(lun2)

      RETURN
      END
