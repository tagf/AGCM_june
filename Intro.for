
      subroutine Intro(istop,rest_off,restfile,restfile_ice)
c       Setup for environment
C************ BEGINNING OF COMMON ************
C
	include 'recom.fi'
C
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
      COMMON /SEA/ SS(74,46,12), SWIN(74,46,12)
      
c
	include 'dir.fi'
c
	include 'ice.fi'
      include 'acc_ice.fi'
C
cIIIIIIIIIIIIIIIIIIIIIIIIIII
      include 'StatIce.fi'
cIIIIIIIIIIIIIIIIIIIIIIIIIII
C************ END OF COMMON ******************
      INTEGER*2 tmphour, tmpminute, tmpsecond, tmphund
      character*4 id0,id1,IDir,Ident
      character*1 ii
      character*40 restfile
      character*60 restfile_ice
      character*8 inptfile,Irest
      logical rest_off
      integer*2 hs386,ms386
c      include'grex.fh'

      rest_off=.true.       ! do not write REST and ACC data
C
C        MAIN COMPUTATIONAL CONTROL
C**************

      OPEN (22,FILE=TRIM(BaseDir)//'\IntrData',	 !Default Names
     1      FORM='FORMATTED')
	read (22,200) IDir,Irest,Ident
200	format (1x,a4,a8,a4)

       write (*,'(a,a4,a,$)') 
	1 ' Input Working Directory name  (ENTER - ',IDir,'): '
            read (*,'(a)') id1
            if (id1.eq.' ')  id1=IDir
	      IDir=id1
            WorkDir='id##'//id1
C        Print general information file

c      OPEN (2,FILE=TRIM(BaseDir)//WorkDir//'\map', !not used
c     1      FORM='FORMATTED',access='append')
!     OPEN output files
      OPEN (7,FILE=TRIM(BaseDir)//WorkDir//'\Protocol',
     1      FORM='FORMATTED', access='append')
      OPEN (13, FILE=TRIM(BaseDir)//WorkDir//'\global.txt',
	1 ACCESS='APPEND')
      OPEN (33, FILE=TRIM(BaseDir)//WorkDir//'\GlIce.txt',
	1 ACCESS='APPEND')
      OPEN (34, FILE=TRIM(BaseDir)//WorkDir//'\IceIJ.txt',
	1 ACCESS='APPEND')
      OPEN (35, FILE=TRIM(BaseDir)//WorkDir//'\IceIJ1.txt',
	1 ACCESS='APPEND')
      OPEN (61, FILE=TRIM(BaseDir)//WorkDir//'\StatIce.txt',
	1 ACCESS='APPEND')
	!special interest points
	ipoint=36
	jpoint=23

	ipoint(1)=34
	jpoint(1)=45
         WRITE (61,126) ipoint(1),jpoint(1),id
 126       FORMAT (2x,i3,2x,i3,2x,a4)
c
c	LUN=42 - acc1
c	LUN=43 - acc
c	LUN=8  - sst
c	LUN=9  - rest
c	LUN=10 - restice
c	LUN=15 - SWIN
c	LUN=12 - surfile
c	LUN=6 -  writeij
c	LUN=34 - IceIJ
c	LUN=35 - IceIJ1
c	LUN=61 - StatIce
c	LUN=44 - AccIce1
c	LUN=45 - AccIce
c	LUN=46 - sb GradsCtl, file AccIce1.ctl
c	LUN=7 -  sb Force, Protocol
c     lun=31 - sb Grads
c     lun=22 - _info File
         WRITE (7,887)
 887     FORMAT (2X,75('*'))
        call PrnDate(6)             ! print date to screen
        call PrnDate(7)             ! print date to file
c      Current time
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund) !initial time
          write (*,210) ' Current time is ',              !Display
     1    tmphour, tmpminute, tmpsecond           
          write (7,210) ' Current time is ',              !File
     1    tmphour, tmpminute, tmpsecond 
  210      format(a,I2,':',I2.2,':',I2.2)
c             DIALOGUE:
       write(*,'(a,$)') ' Input time when to stop (hh mm): '
       read (*,'(I2,a,I2)') hs386,ii,ms386
       write (*,210) ' Time to stop: ',hs386,ms386
         iStop=60*hs386+ms386
      write(*,'(a,$)') ' To set ZERO Accumulated values (1-Yes/0-No): '
      read (*,'(i1)') iacc

      write(*,'(a,$)') ' To Write REST and ACC data (1-Yes/0-No): '
      read (*,'(i1)') ichoice
      if (ichoice.eq.1) Rest_off=.false.    ! to write REST&ACC

      write(*,'(a,$)') ' To set  Year number one (1-Yes/0-No): '
      read (*,'(i1)') ichoice
 
            inptfile=Irest   !Restart data files names
            restfile=TRIM(BaseDir)//WorkDir//'\'//inptfile
            restfile_ice=TRIM(BaseDir)//WorkDir//'\'//
	1	    TRIM(inptfile)//'ice'
C
C        INITIAL CONDITIONS
      
      CALL INPUT(restfile,9,restfile_ice,10)
      if (ichoice.eq.1) sdeyr=1.
      if (tau.gt.1.e6) tau=24.
          if (iacc.eq.1) then
            call acc0       ! to put zero accum values
            CALL OUTACC(43,0)
c           CALL OutAcIce(45,0)  !Hice, comp data file for Grads
c           call GradsCtl(46)   !Data description file for Grads

          else
            CALL INACC(43)
c            CALL InAccIce(45)
          endif

C***********************
C      FOR ICE MODELLING
C***********************
c      write(*,'(a,$)') ' Use pre-set wind in ocean model (1-Yes/0-No): '
c      read (*,'(i1)') ichoice
c      if (ichoice.eq.1) then
c	 CalcWind =.false.
c	else
	 CalcWind =.true.
c	endif
c  онстанты
      HActWater   = 3000d0  ! толщина де€тельного сло€, см
      CTWater     = 1d0     ! теплоемкость воды, кал/см**3
      CTYSN       = 7.88d-4 ! теплопроводность снега
      ROSN        = 0.32d0  ! плотность снега, г/см**3
      QIce        = 72.0d0  ! теплота плавлени€ льда, кал/см**3
      QSnow       = 26.2d0  ! теплота плавлени€ снега, кал/см**3
      QWater      = 64.0d0  ! теплота замерзани€ воды, кал/см**3
	GIceCompMin = 0.02d0  ! минимум сплоченности (not used)
	GIceCompMax = 0.98d0  ! максимум сплоченности
	GIWExchage  = 0.55d0  ! теплообмен лед-вода
	HIceMin     = 5d0     ! минимум толщины льда, см
	HIceThin    = 5d0     ! толщина, ниже которой лед считаетс€ тонким, см
	GIceSThruPercent   = 0.3d0   ! дол€  ¬ радиации, проникающей в лед
	GIcePoolMaxPercent = 0.3d0   ! макс. дол€ тепла, накапливаемого в полынь€х
	FixedWater         =  .TRUE.  ! .FALSE.темп. воды - с натуры
	FixedIceBorder     =  .TRUE. ! .FALSE.граница льдов - с натуры
	iDaysToFreeze = 2     ! 7 по прошествии этого времени вырастает новый лед
	SFCALB(9,1) =  0.75d0 !0.45d0 0.75d0  ! sea ice поправка на нехватку снега
	SFCALB(9,2) =  0.80d0  ! 0.75d0!Sea ice + snow
c	SFCALB(1,1) =  0.12 ! 0.05d0 !woodland без снега
c	SFCALB(2,1) =  0.10 ! 0.05d0 !forest без снега
c	SFCALB(3,1) =  0.13 ! 0.07d0 !steppe, grassland без снега
c	SFCALB(4,1) =  0.20 ! 0.05d0 !steppe desert без снега
c	SFCALB(5,1) =  0.25 ! 0.05d0 !desert без снега
c	SFCALB(6,1) =  0.19 !0.1d0  !тундра без снега
	GIceUpFlow  =  3.0 !3.0 1.0!5.0  ! 7.5   !15.d0 !!!
c	SNOAMT = 0.0   
     
       CALL INTSEA
c******** OCEAN
       CALL BUFER(NYER,MEC,NDAY,0)

      write(*,*) ' To reset ice conditions (1-Yes/0-No)'
      read (*,'(i1)') ichoice
      if (ichoice.eq.1) then
c******** OCEAN
       call init(0, MONTH)
c******** OCEAN
	GIceComp    = 0.0
      GHICE       = -2d0
       do j=1,46
        do i=1,74
c	   if (IceCompFlows) then
c	    do k=1,12
c	     GWFlowDays(k)    = 0
c	     GWFlowAve(i,j,k) = 0
c	    enddo
c	    GWFlowAcc(i,j) = 0
c	   else
c	    GWFlowAcc(i,j) = 1d20 ! заведомо много
c	    GIcePool(i,j)  = 0
c	    GIceComp(i,j)  = 0.02
c         endif
c         dTWater(I,J) = 0d0
c         if (ISFTYP(I,J).eq.7.or.ISFTYP(I,J).eq.9) then
c          if (SS(I,J,MONTH).gt.200d0) then  ! sea ice
          if (ISFTYP(I,J).eq.9) then  ! sea ice
            GHICE(I,J)  = 203d0
c		  SNOAMT(I,J) = 0d0	  take from input data
c           ISFTYP(I,J) = 9	  nonsence
            GIceComp(I,J) = 0.96 !GIceCompMax
          endif
c          else
         if (ISFTYP(I,J).eq.7) then	 !open water
           GHICE(I,J)  = 0d0
c           ISFTYP(I,J) = 7
c           GIceComp(I,J) = 0d0
         endif
c         else
c          GHICE(I,J) = -2d0
c          GIceComp(I,J) = 0d0
c         endif
        enddo
       enddo
c       GIceComp    = GIceCompMax
	 GT1         = GT
	 GT2         = GT
       SnoAmt1     =  0d0 !SNOAMT
       SnoAmt2     =  0d0 ! SNOAMT
       GIceCND     = 0d0
	 GIceHeatng  = 0d0
	 GIce2Heatng = 0d0
	 GIcePool    = 0d0
	 GIce2Pool   = 0d0
	 NewIce      = .FALSE.
	 iGIceFreezeCount   = 0
      else ! not resetting
c******** OCEAN
       call init(4, MONTH)
c       do j=1,46
c        do i=1,74
c         if (ISFTYP(I,J).eq.7.or.ISFTYP(I,J).eq.9) then
c          if (GHICE(I,J).lt.20d0) then
c           ISFTYP(I,J) = 7
c           GIceComp(I,J) = 0d0
c          else
c           ISFTYP(I,J) = 9
c           GIceComp(I,J) = GIceCompMax
c          endif
c	   else
c          GIceComp(I,J) = 0d0
c         endif
c	   GT1(I,J) = GT(I,J)
c	   GT2(I,J) = GT(I,J)
c        enddo
c       enddo
      endif
c***
      Ibreak = 2
      Jbreak = 2
C***********************

c********  For sea ice experiment ****************
c      do 102 kk=1,12
c      do 102 j=12,16
c      do 102 i=1,26
c       if((isftyp(i,j).ne.7).or.(ss(i,j,kk).gt.200.)) go to 102
c      ss(i,j,kk)=ss(i,j,kk)+tice
c 102  continue
c************************************************
C     not used:
c         nhis=1
c         nres=1

C      CALL INPUTO
c       FM=0.7D-5
C*******************************************************************
C          SIMULATION IDENTIFICATION
      write(*,'(a,a4,a,$)')
     1' Input SIMULATION IDENTIFICATION (4 char.) (ENTER - ',Ident,'): '
      read (*,'(a)') id0
      if (id0.eq.' ') id0=Ident
	Ident=id0
c      Call clear_text
c      call home()
           write (*,'(/1x,a,a,a)') ' ***** ', id0,' *****'
C*******************************************************************
           ID=ID0
	rewind (22)
	write (22,200) IDir,Irest,Ident	  !Default file names
	close (22)
      return
	end
