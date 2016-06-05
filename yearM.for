      program YearMean
c
c       Reads data from Global.txt, Calculates Year means
c       and Writes Global values to File globYear.txt

	real TAU,GMKEAC,GMTACC,GMRACC,SDEYR
	real ke,t,r
	integer MONTH,SDEDY
	integer n,ny
	character*4 id
	character*8 WorkDir
	character*30 BaseDir

      BaseDir='e:\climdig\agcm_ice3\'
      WorkDir='id##00b6'

      OPEN (13, FILE=TRIM(BaseDir)//WorkDir//'\glob.txt',
	1 form='formatted',status='old')
      OPEN (14, FILE=TRIM(BaseDir)//WorkDir//'\globYear.txt',
	1 ACCESS='APPEND')
	ny=0 ! year number
 2	n=0  ! year day counter
	ke=0. ! K.e.
	t=0.  ! temperature
	r=0.  ! radiation
	do while (n.le.365)

C        GLOBAL T, K.E. ,N0 - DAY MEANS - to file GLOBAL
	! To GLOBAL File:
      read (13,10,end=1) TAU,MONTH,GMKEAC,GMTACC,GMRACC,SDEDY,SDEYR,id
 10       FORMAT (1x,F9.2,2X,I3,2X,3f9.4,2x,i3,2x,f5.0,2x,a4)
	n=n+1
	ke=ke+GMKEAC
	t=t+GMTACC
	r=r+GMRACC
	enddo
	ke=ke/365. 
	t=t/365. 
	r=r/365.
	ny=ny+1
	write (14,11) ny,ke,t,r,sdeyr,id 
 11       FORMAT (1X,I3,2X,3f9.4,2x,f5.0,2x,a4)
	goto 2 
 1    stop 'normal end'
      end

