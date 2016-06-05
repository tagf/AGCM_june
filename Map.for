C****************************************************************
      subroutine PrnMap(nav1,a1,im1,jm1,shift,scale,title,mname)
C***** PRINT OUT ARRAY WITH TITLE FOR MAP                  ******
C***** IM1=72 OR IM1=74  JM1=46 ** 4 CHARACTERS PER CELL   ******
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
       DIMENSION A1(IM1,JM1)
c      character*4 TITLE                                                      
      character*(*) TITLE, MNAME
      DIMENSION A(24,16),IA(24,16),IH1(24),I1(24),LAT1(16)
      DATA IH1/12*'W',' ',11*'E'/
      lun=2
c      open (lun,file='map',form='formatted')
	OPEN (lun,FILE=TRIM(BaseDir)//WorkDir//'\'//MNAME,
     1      FORM='FORMATTED',ACCESS='APPEND')
C***      IM1=74
       IF (IM1.EQ.72) GO TO 8
C***   J= 2 TO 15
       DO 9 I3=1,24
       DO 9 J3=2,15
       A(I3,J3)=0.
       DO 1  I0=1,3
       DO 1  J0=1,3
       I=3*(I3-1)+I0+1
       J=3*(J3-1)+J0-1
 1     A(I3,J3)=A(I3,J3)+A1(I,J)
       A(I3,J3)=A(I3,J3)/9.
 9     CONTINUE
C***   J=1,J=16
       DO 91 I3=1,24
       A(I3,1)=0.
       A(I3,16)=0.
       DO 911 I0=1,3
       DO 911 J0=1,2
       I=3*(I3-1)+I0+1
       A(I3,1)=A(I3,1)+A1(I,J0)
 911   A(I3,16)=A(I3,16)+A1(I,44+J0)
       A(I3,1)=A(I3,1)/6.
  91   A(I3,16)=A(I3,16)/6.
       GO TO 4
 8     CONTINUE
C***      IM1=72
C***      J=2 TO 15
       DO 3  I3=1,24
       DO 3  J3=2,15
       A(I3,J3)=0.
       DO 15 I0=1,3
       DO 15 J0=1,3
       I=3*(I3-1)+I0
       J=3*(J3-1)+J0-1
 15    A(I3,J3)=A(I3,J3)+A1(I,J)
       A(I3,J3)=A(I3,J3)/9.
 3     CONTINUE
C***   J=1,J=16
       DO 92 I3=1,24
       A(I3,1)=0.
       A(I3,16)=0.
       DO 921 I0=1,3
       DO 921 J0=1,2
       I=3*(I3-1)+I0
       A(I3,1)=A(I3,1)+A1(I,J0)
  921  A(I3,16)=A(I3,16)+A1(I,44+J0)
       A(I3,1)=A(I3,1)/6.
  92   A(I3,16)=A(I3,16)/6.
  4   DO 6 I=1,24
      DO 6 J=1,16
  6   IA(I,J)=(A(I,J)/NAV1-SHIFT)/SCALE+0.1
  5   continue
      write (lun, 10)TITLE
      write (lun, 13) MONTH,MNTHDY,TOFDAY,ID,shift,scale,nav1,nstep,tau
      DO 7 J=1,16
  7   LAT1(J)=IABS(90-12*(J-1))
      DO 2 I=1,24
  2   I1(I)=IABS(180-15*(I-1))
      write (lun, 11) (I1(I),IH1(I),I=1,24)
      write (lun, 117)
      write (lun, 12) LAT1(16),(IA(I,16),I=1,24),LAT1(16)
      write (lun, 116)
      write (lun, 12) LAT1(15),(IA(I,15),I=1,24),LAT1(15)
      write (lun, 115) 
      write (lun, 12) LAT1(14),(IA(I,14),I=1,24),LAT1(14)
      write (lun, 114) 
      write (lun, 12) LAT1(13),(IA(I,13),I=1,24),LAT1(13)
      write (lun, 113) 
      write (lun, 12) LAT1(12),(IA(I,12),I=1,24),LAT1(12)
      write (lun, 112) 
      write (lun, 12) LAT1(11),(IA(I,11),I=1,24),LAT1(11)
      write (lun, 111) 
      write (lun, 12) LAT1(10),(IA(I,10),I=1,24),LAT1(10)
      write (lun, 110)
      write (lun, 12) LAT1(9 ),(IA(I,9 ),I=1,24),LAT1(9 )
      write (lun, 109)
      write (lun, 12) LAT1(8 ),(IA(I,8 ),I=1,24),LAT1(8 )
      write (lun, 108)
      write (lun, 12) LAT1(7 ),(IA(I,7 ),I=1,24),LAT1(7)
      write (lun, 107)
      write (lun, 12) LAT1(6 ),(IA(I,6 ),I=1,24),LAT1(6)
      write (lun, 106)
      write (lun, 12) LAT1(5 ),(IA(I,5 ),I=1,24),LAT1(5)
      write (lun, 105)
      write (lun, 12) LAT1(4 ),(IA(I,4 ),I=1,24),LAT1(4)
      write (lun, 104)
      write (lun, 12) LAT1(3 ),(IA(I,3 ),I=1,24),LAT1(3)
      write (lun, 103)
      write (lun, 12) LAT1(2 ),(IA(I,2 ),I=1,24),LAT1(2)
      write (lun, 102)
      write (lun, 12) LAT1(1 ),(IA(I,1 ),I=1,24),LAT1(1)
      write (lun, 101)
 10   FORMAT (1X,5('*'),A)
 11   FORMAT (4X,24(I3,A1))
 12   FORMAT (1X,I2,'I',24I4,'I',I2)
 13   FORMAT (2X,I3,1H/,I2,1H/,F6.2,3H.00,5X,'ID=',A4,2X,'SHIFT=',F6.1,
     1      2X,'SCALE=',D8.1,2X,'NAV=',I4,2X,'NSTEP=',I7,2X,'TAU=',F9.2)
 117  FORMAT (1X,'N I',96('='),'I',' N')
 116  FORMAT (3X,'I',24X,20('='),32X,4('='),16X,'I',3X)
 114  FORMAT (3X,'I',4X,12('='),8X,4('='),4X,8('='),48X,8('='),'I')
 115  FORMAT (3X,'I',4X,20('='),8X,4('='),4X,4('='),
     1       8X,24('='),4X,16('='),'I')
 113  FORMAT(3X,'I',24X,4('='),4X,4('='),12X,12('='),24X,4('='),8X,'I')
 112  FORMAT (3X,'I',16X,4('='),8X,4('='),20X,8('='),36X,'I')
 111  FORMAT (3X,'I',20X,8('='),16X,4('='),16X,4('='),4X,
     1       4('='),4X,4('='),12X,'I')
 110  FORMAT(3X,'I',28X,8('='),8X,4('='),20X,4('='),4X,4('='),16X,'I'  )
 109  FORMAT (3X,'I',36X,4('='),8X,4('='),8X,4('='),32X,'I')
 108  FORMAT (3X,'I',28X,4('='),52X,4('='),8X,'I')
 107  FORMAT (3X,'I',36X,4('='),12X,4('='),24X,4('='),4X,4('='),4X,'I')
 106  FORMAT (3X,'I',56X,4('='),20X,12('='),4X,'I')
 105  FORMAT (3X,'I',32X,4('='),60X,'I')
 104  FORMAT (3X,'I',76X,4('='),16X,'I')
 103  FORMAT (3X,'I',8X,28('='),8X,32('='),4X,16('='),'I')
 102  FORMAT (3X,'I',8('='),28X,8('='),52X,'I')
 101  FORMAT (1X,'S I',96('='),'I',' S')

      CLOSE(lun)
      RETURN
      END
