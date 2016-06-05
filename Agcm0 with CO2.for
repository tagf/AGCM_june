
C     ************
      PROGRAM AGCM
C     ************
C             2 - LEVEL AGCM
C             5x4 DEGREES VERSION
C************ BEGINNING OF COMMON ************
C
	USE MSFLIB
	USE DFLIB
	USE WinMod
C
	include 'recom.fi'
      DIMENSION C(900)
      EQUIVALENCE (TAU,C(1))

      LOGICAL SNOW,SEAICE
C
	include 'radvar.fi'
C
      COMMON /FACTR/  EST0(74),E10(74),E20(74),E02(74),E32(74),
     1                E12(74),E42(74),E24(74),E34(74)
C
      COMMON /MOIST/ EFVT(74),EFVST(74),EFV3(74),EFV2(74),
     1               EFV1(74),EFV0(74),AK(74)
C
	include 'cover.fi'
	include 'trans.fi'
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1   ,SD(74,46),PIV(74,46,2)
C
	include 'accum.fi'
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)
C
C        QT ARRAY - TEMPORARY STORAGE FOR SIMULATION
C
      COMMON /QTARY / PT(74,46),UT(74,46,2),VT(74,46,2),TT(74,46,2)
     1  ,QWT(74,46,2)
C
      COMMON /QTTARY/ QTT(74,46,9)
C
      COMMON /COMP/ DUM(74,46,4)
C
      COMMON /WORK/ DUM2(72,46,9)
C
      COMMON /SEA/ SS(74,46,12), SWIN(74,46,12)
      
         COMMON /WINDatm/ TAUU(74,46),TAUV(74,46)
      common/glacc/GMTACC,GMKEAC,GMRACC,GMT4ACC,namgl,restart
	logical  restart
c
	include 'dir.fi'
C
	include 'ice.fi'
      include 'acc_ice.fi'
C
cIIIIIIIIIIIIIIIIIIIIIIIIIII
      include 'StatIce.fi'
cIIIIIIIIIIIIIIIIIIIIIIIIIII
      COMMON /OSA/ QSN(72,46),TXN(72,46),TYN(72,46)
        !@@@@@@@@@@@@@ tarko data - accum and write data for Tarko program
         COMMON /Tarko/ PACCt(72,46),UACCt(72,46,2),VACCt(72,46,2),     1                   
     1                   SDACCt(72,46),navt
C
c        For program execution time
      INTEGER*2 tmphour, tmpminute, tmpsecond, tmphund
      logical stop_time,end_of_day,rest_off,first   !,stop_key
      character*40 restfile
      character*60 restfile_ice
	logical DT300 !small DT flag
c%%%%%%%%%%%%%!for surfer files slides

      real*4 arr(72,46)
	integer*4 date(6)
      character*9 fileName
      character*2 fileN(30)
	data fileN /'01','02','03','04','05','06','07','08','09','10',
     1      '11','12','13','14','15','16','17','18','19','20',
     2      '21','22','23','24','25','26','27','28','29','30'/
c%%%%%%%%%%%%%

c	dimension aaa(74,46) !temporal for writeij
c     CO2 concentration transmission function 
      data    co2 /1.d0/ !CO2 concentration factor
      trco(x,y)=dmin1(1.e0,1.187e0-.066e0*dlog10(dabs(x-y)*(x+y))
	1          -0.066*dlog10(co2)) !additional term 
c         ABS(XXX)=DABS(XXX)
C
c      Dialogue and Setup
	 BaseDir = 'e:\climdig\AGCM_ICE3\'  !'e:\climdig\AGCM_DLG\' 
       call Intro(istop,rest_off,restfile,restfile_ice)
	
        !@@@@@@@@@@ tarko data
      PACCt=0.
      UACCt=0.
      VACCt=0.     
      SDACCt=0.
      lunt=67
      navt=0
      OPEN(lunt,FILE=TRIM(BaseDir)//WorkDir//'\tarko',ACCESS='append')
        !@@@@@@@@@@@ tarko data
	
	 DT1=DT   !480. ! 300.  ! 600.
	 Dt300=.false.
	 restart=.false. !restart near overflow (set in sb force)
	nfile=0 !for surfer files slides

	 goto 1001
 1000	if (restart) then !slow down DT near owerflow and run
	  call Intro1(restfile,restfile_ice) !restart from beginning of month
	  restart=.false.
	  DT1=   300.
	  Dt300=.true.
	  N_ofDT1steps=0  !10 days use DT=300
      endif

c     CO2 concentration for layers
       tct0=trco(0.0d0,200.d0)
       tct2=trco(0.0d0,600.d0)
       tct4=trco(0.0d0,1000.d0)
       tcst0=trco(100.0d0,200.d0)
       tcst2=trco(100.0d0,600.d0)
       tcst4=trco(100.0d0,1000.d0)
       tc01=trco(200.0d0,400.d0)
       tc02=trco(200.0d0,600.d0)
       tc03=trco(200.0d0,800.d0)
       tc04=trco(200.0d0,1000.d0)
       tc12=trco(400.0d0,600.d0)
       tc14=trco(400.0d0,1000.d0)
       tc23=trco(600.0d0,800.d0)
       tc24=trco(600.0d0,1000.d0)
       tc34=trco(800.0d0,1000.d0)
c
	 
c	call printc
c          Initial time point
 1001   CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund) !initial time
      time0 =tmphour*3600.+ tmpminute*60.+tmpsecond+tmphund*0.01 ! (sec)

	 Ibreak = 3
	 Jbreak = 31
           stop_time=.false.
           end_of_day=.false.
	     stop_key=.false.
           first=.true.          ! for stop key
c	     nstop=0               ! stop key print counter
        write (6,*)' *****Program  Execution Began...    **************'
	write (6,*)' *****To Stop Running Program - Press ANY KEY *****'
c********** XOCEAN ***
c      CALL BUFER(NYER,MEC,NDAY,0)
c      CALL INTSEA
c      call init(0, MONTH)
c********************
c           ODT=10.*24.*3600.
       TAU=24.
       TAUE=TAU
c      Time step dependent variables
c  !!!!
	 DT=DT1   !480. ! 300.  ! 600.
       DTC3=DT*6.
       TSPD=24.*3600./DTC3
       COE1=COE*DTC3/(24.0*3600.)
	 iDaysToFreeze = 5*tspd     ! по прошествии этого времени (5 days)
	                          ! вырастает новый лед
c         TDSEA=.FALSE.
c         SWAMP=.FALSE.
c         OIL  =.FALSE.
c      call PrintC
c      write (6,91) ((isftyp(i,j),i=2,73),j=46,1,-1)
c 91   format (//(72i1))            

c	 print '(9I8)' ,(i, i=1,9) !Print surface albedo
c	do j=1,2  
c	 print '(9e8.2)', (sfcalb(i,j), i=1,9)
c	enddo
C
c        COMPLETE TIME STEP USING COMP3

      DO 1 J=1,JM
  1   CALL COMP3 (J)

c         CALL OUTRES                ! Restart records to file
c      stop 'end comp3'
c      call outacc
 2222 DTH=DT/3600.0
      NSTEP=TAU/DTH+0.1
      NHIS=TAUH/DTH+0.01
      NGMP=TAUC/DTH+0.01
      NSDET=TAUD/DTH+0.01
      NRES=TAUD/DTH+0.01
      end_of_day= (MOD(NSTEP,NSDET).EQ.0)

C        Simulation  identification print (to file WorkDir\Protocol)
         WRITE (7,888)
         WRITE (7,889) ID
         WRITE (6,889) ID
         WRITE (7,888)
 889     FORMAT (35X,'ID = ',A4)
 888     FORMAT (20X,40('*'))

      WRITE(7,804)NSTEP,NGMP,NSDET,NRES,IM,JM,NCYCLE,NC3,SDEDY,NAV
      WRITE(7,803) DT,DTH,TAUH,TAU,TAUI,TOFDAY,SDEYR,TSPD,fm,co2
      WRITE(7,805) FixedWater,FixedIceBorder,iDaysToFreeze/tspd,
	1	    SFCALB(9,1),GIceUpFlow(1,1)

      WRITE(6,804)NSTEP,NGMP,NSDET,NRES,IM,JM,NCYCLE,NC3,SDEDY,NAV
      WRITE(6,803) DT,DTH,TAUH,TAU,TAUI,TOFDAY,SDEYR,TSPD,fm,co2
 804  FORMAT (' *STEP =',I7,' *NGMP =',I7,' *NSDET=',I7,
     1   ' *NRES =',I7,' *IM   =',I7/' *JM   =',I7,' *NCYCL=',I7,
     2  ' *NC3  =',I7,' *SDEDY=',I7,' *NAV  =',I7)
 803  FORMAT (' *DT  =',F8.2,' *DTH =',F8.2,' *TAUH=',F8.2,
     1' *TAU =',F8.2,' *TAUI=',F8.2,/' *TOFDAY=',F8.2,
     2' *SDEYR =',F8.2,' *TSPD  =',F8.2,' *fm=',d12.3,/' *co2=',F8.2/)
 805  FORMAT (' *FixedWater  =',L1,' *FixedIceBorder =',L1,
     1' *iDaysToFreeze=',F6.2,/' *SFCALB =',F6.2,' *GIceUpFlow=',F6.2,/)

      IF (end_of_day) then
        CALL SDET
c********** XOCEAN ***
        CALL INTSEA       ! INTERPOLATE SEA TEMPERATURES
	 if (.not.FixedWater) then
         call ocean
         TAUU=0.d0          ! ZERO wind stress EVERY DAY
         TAUV=0.d0
         QSN=0.d0           ! ZERO SEA HEAT EVERY DAY
	  endif
c********************
        call IceDaily     ! once-per-day ice calculations
c       aaa = NEWICE
c      CALL WRITIJ(aaa,1.d0,'NEWICE',0)
        CALL INTOZ
        CALL GMP
       endif
C

 10   continue	!next step
 
c      stop_key=PEEKCHARQQ() !not used in QWin
c      stop_key=.FALSE.
      if (stop_key.and.first) then !stop_key - from sb StopRun
      write (6,*)' ****** STOP AFTER FULL DAY. WAIT ...*****'
        first=.false.
      endif 
      
      NSTEP=NSTEP+NC3
      end_of_day= (MOD(NSTEP,NSDET).EQ.0)
      TAU=(NSTEP*DT)/3600.0
      TOFDAY=DMOD(TAU-1.0,ROTPER)+1.0
      ROT=TOFDAY/ROTPER*(PI+PI)
      COSR=DCOS(ROT)
      SINR=DSIN(ROT)
C
C         INTEGRATE FORWARD
C
      CALL STEP

c      CALL WRITIJ(PHIS,1.d0,'PHIS',0)
c      CALL WRITIJ(P   ,1.d0,'P   ',0)
c      CALL WRITIJ(U(1,1,1),1.d0,'U1  ',2)
c      CALL WRITIJ(U(1,1,2),1.d0,'U3  ',2)
c      CALL WRITIJ(V(1,1,1),1.d0,'V1  ',2)
c      CALL WRITIJ(V(1,1,2),1.d0,'V3  ',2)
c      CALL WRITIJ(T(1,1,1),1.d0,'T1  ',0)
c      CALL WRITIJ(T(1,1,2),1.d0,'T3  ',0)
c      CALL WRITIJ(QW(1,1,1),1.d-4,'QW1 ',0)
c      CALL WRITIJ(QW(1,1,2),1.d-4,'QW3 ',0)
c      CALL WRITIJ(GW  ,1.d-4,'GW  ',0)
c      CALL WRITIJ(GT  ,1.d0,'GT  ',0)
c      CALL WRITIJ(SNOAMT,1.d-4,'SNOW',0)
c      CALL WRITIJ(TS    ,1.d0,'TS  ',0)
c      call PrnMap(1,p,74,46,800.d0,1.d0,' P ')
c      call PrnMap(1,gt,74,46,273.d0,1.d0,' GT ')
c      call PrnMap(1,ts,74,46,273.d0,1.d0,' TS ')
c      call PrnMap(1,u(1,1,1),74,46,0.d0,1.d0,' U1 ')
c      call PrnMap(1,u(1,1,2),74,46,0.d0,1.d0,' U3 ')
c      call PrnMap(1,v(1,1,1),74,46,0.d0,1.d0,' V1 ')
c      call PrnMap(1,v(1,1,2),74,46,0.d0,1.d0,' V3 ')
c      call PrnMap(1,t(1,1,1),74,46,273.d0,1.d0,' T1 ')
c      call PrnMap(1,t(1,1,2),74,46,273.d0,1.d0,' T3 ')
c      call PrnMap(1,qw(1,1,1),74,46,0.d0,1.d-4,' Q1 ')
c      call PrnMap(1,qw(1,1,2),74,46,0.d0,1.d-4,' Q3 ')
c      call PrnMap(1,snoamt,74,46,0.d0,1.d0,' SN')
c      call PrnMap(1,sd(1,1),74,46,0.d0,1.d8,' SDe ')
c         go to 606            !  If only one step
    
      IF (end_of_day) then
        CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)
        iTime=60*tmphour+tmpminute
c         Check when to stop
          if ((iabs(iStop-iTime).lt.2).or.stop_key) stop_time=.true.
c      if (SDEYR.ge.9.9.and.month.eq.3) stop_time=.true. !!!!!Temporal!!!!
      endif

      if (rest_off) goto 9997
      IF (end_of_day.and.
     1 (mnthdy.eq.1.or.stop_time))
     2    CALL OUTRES(restfile,9,restfile_ice,10)      ! Restart records to file

 9997 continue

      IF (end_of_day) then !10 days with small DT and then put large DT
	 if (Dt300) then
	  N_ofDT1steps=N_ofDT1steps+1
	  if (N_ofDT1steps.eq.10) then
	     DT300 =.false.
	     N_ofDT1steps=0
           DT1=480.
	     goto 1001
	  endif
	 endif
	endif
C
C         ADD FORCING TERMS
C
       CALL FORCE

	if (restart) goto 1000 ! if near overflow, to change DT 

c       CALL outac1(42)
c       STOP 'Normal end'
C     OCEAN MODEL
C###      IF (end_of_day.AND.TDSEA) CALL OCEAN
C
c      IF (mod(nstep,nhis).eq.0) then
      IF (end_of_day) then
        WRITE (6,999) MONTH,MNTHDY,TOFDAY,SDEYR,GMKE,GMT,GMR,HIceGl
 999    FORMAT (1X,I3,1H/,I2,1H/,F6.2,2X,5Hyear=,f4.0,2X,5HK.E.=,
     1   1PD10.3,3H T=,1PD10.3,4H N0=,1PD10.3,6H HIce=,1PD10.3)

c        nstop=nstop+1    ! Print announsement to screen
c         if (nstop.eq.21) then
c         write (6,*)' *****To Stop Running Program - Press ANY KEY*****'
c	   nstop=0
c         endif
      endif

C
C         OUTPUT ACCUMULATED VARIABLES HISTORY
C
      if (rest_off) goto 123
      IF (end_of_day.and.
     1 (mnthdy.eq.1.or.stop_time)) then
                 CALL OUTACC(43,0)  ! Accumulated variables for restart
                 CALL OutAcIce(45,0)
       ENDIF
 123   continue
        !@@@@@@@@ tarko data
	IF (end_of_day.and.(MOD(SDEDY,10).EQ.0)) then
       PACCt=PACCt/navt
       UACCt=UACCt/navt
       VACCt=VACCt/navt       
       SDACCt=SDACCt/navt
        write (lunt,1300) SDEDY,'ten day records'
 1300     FORMAT(1X,I3,2X,A)
          DO  I=1,72
           DO  J=1,46
            write (lunt,1190) i,j,PACCt(i,j),UACCt(i,j,1),UACCt(i,j,2),
     1		   VACCt(i,j,1),VACCt(i,j,2),SDACCt(i,j)
	     enddo
	    enddo
 1190    FORMAT (1x,2i3,1x,6e12.4)
       PACCt=0.
       UACCt=0.
       VACCt=0.     
       SDACCt=0.
       navt=0
	ENDIF
        !@@@@@@@@ tarko data

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%!for surfer files slides
      if (rest_off)       goto 9996
        goto 9996
      IF (abs(TOFDAY-8.).le.0.01.or.abs(TOFDAY-16.).le.0.01
     1           	.or.abs(TOFDAY-24.).le.0.01)  then

	nfile=nfile+1
c       do 36 J=1,JM    !SLP
c       do 36 I=2,73
c        arr(i-1,j)=sngl((P(I,J)+PTROP)*
c     1  EXP(PHIS(I,J)/(RGAS*(1.5*T(I,J,2)-0.5*T(I,J,1)+FLR*PHIS(I,J)))))
c	2   -1000.
c 36    continue
	  date(1)=int(tau)
	  date(2)=month
	  date(3)=mnthdy
	  date(4)=int(tofday)
	  date(5)=id
	  date(6)=sdedy
	! U3
        fileName='U3'//fileN(nfile)//'.dat'
        open (70, file=fileName, status='new') !,err=89)
c        goto 90
c 89    print '(a,$)' , ' ERROR: File name Duplicate.'
c       goto 88

        write (70,13) ' U3 ',date(1),'h',date(2),date(3),date(4)
     1   ,date(5),date(6),'* * * * * * * *'
 13     FORMAT(A,1X,I7,A,2X,I3,1H/,I2.2,1H/,i2.2,3H.00,2X,
     1   4H ID=,A4,2X,' day=',I3.3,A)
          DO  I=1,72
           DO  J=1,46
            write (70,119) i,j,U(i+1,j,2) !arr(i,j)
	     enddo
	    enddo
 119    FORMAT (1x,2i3,1x,e12.4)

         close (70)
	! V3
        fileName='U3'//fileN(nfile)//'y.dat'
        open (70, file=fileName, status='new') !,err=89)
c        goto 90
c 89    print '(a,$)' , ' ERROR: File name Duplicate.'
c       goto 88

        write (70,13) ' V3 ',date(1),'h',date(2),date(3),date(4)
     1   ,date(5),date(6),'* * * * * * * *'
          DO  I=1,72
           DO  J=1,46
            write (70,119) i,j,V(i+1,j,2) !arr(i,j)
	     enddo
	    enddo

         close (70)

       ENDIF
	if (nfile.eq.30) goto 606  !stop 
 9996 continue
c%%%%%%%%%%%%%%%%%%

C         OUTPUT OCEAN VARIABLES
C**      IF (end_of_day.AND.TDSEA) CALL OUTOCN
C
C         UPDATE DAY, EARTH-SUN DISTANCE, SOLAR DECLINATION
C
      IF (end_of_day) CALL SDET

C        OUTPUT MONTHLY ACCUMULATED VARIABLES HISTORY

      if (rest_off) goto 125
        IF     ((MONTH.EQ.2.OR.MONTH.EQ.8).and.        ! January and July output
     *        (mnthdy.eq.1).and.end_of_day)  then
c              CALL PRZON(11)     ! PRINT OUT MONTHLY ZONAL MEANS
              call outacc(42,1)   ! Accumulated variables
c             call maps           ! Print out maps
        endif
        IF     ((mnthdy.eq.1).and.end_of_day)  then     ! Every month output
              call OutAcIce(44,1)   ! Accumulated Ice variables
        endif
125   continue

	numgl1=namgl !force per day

      if (end_of_day)  call GlobAcc(13)       ! Print global values

c     Put zero accumulations in the end of month
       IF (end_of_day.and.mnthdy.eq.1) then
         CALL acc0
c          Time elapsed
c          Print time elapsed (h,min,sec):
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund) 
      time1 =tmphour*3600.+ tmpminute*60.+tmpsecond+tmphund*0.01
      sec=time1-time0
	     time0=time1
       if (sec.lt.0.) sec=sec+24.*3600.
	secday=sec/(tau-taue)*24.
	taue=tau
       
           ihours= floor(sec/3600.0d0)  !idint
           jmin=floor((sec-ihours*3600)/60.0d0) !idint
           sec1=sec-dfloat(ihours*3600+jmin*60)
           write(6,11) ihours,jmin,sec1,secday
           write(7,11) ihours,jmin,sec1,secday
 11        format (1x,'Time elapsed: ',i2,'h ',i2,'min ',f5.2,'sec',
	1             '   Sec/day: ',f5.2)
      endif
	
            !########### 1 January stop time
c        IF     ((MONTH.EQ.1).and.      
c     *        (mnthdy.eq.2).and.end_of_day)    stop_time= .true.
      if (stop_time) go to 606

      IF (end_of_day) then
c********** XOCEAN ***
       CALL INTSEA       ! INTERPOLATE SEA TEMPERATURES
	 if (.not.FixedWater) then
  
        call ocean
        QSN=0.d0           ! ZERO SEA HEAT EVERY DAY
	 endif
c********************
c       aaa = NEWICE
c      CALL WRITIJ(aaa,1.d0,'NEWICE',0)
c       aaa = iGIceFreezeCount
c      CALL WRITIJ(aaa,1.d0,'FrzCount',0)
        call IceDaily     ! once-per-day ice calculations
       CALL INTOZ        ! INTERPOLATE OZONE AMOUNTS
       CALL GMP          ! CHECK MASS CONSERVATION

	NIce=0	  !Output global sea ice data
	CompGl=0.0	 !compactness
	HIceGl=0.0	 !ice thickness
	SIceGl=0.0   !ice area
	do i=2,73
	do j=24,46   !North Hemisphere
	 if (isftyp(i,j).eq.9) then 
	  NIce=NIce+1
	                           !GIceComp(i,j)
	  HIceGl=HIceGl+GHIce(i,j)*GIceComp(i,j)*dxyp(j)
	  CompGl=CompGl+GIceComp(i,j)*dxyp(j)
	  SIceGl=SIceGl+dxyp(j)
	 endif
	enddo
	enddo
	if (NIce.eq.0) NIce=1
c	HIceGl=HIceGl/SIceGl 
	HIceGl=HIceGl/CompGl
	CompGl= CompGl/SIceGl
	SIceGl=SIceGl/1.d12
       WRITE (33,111) TAU,MONTH,HIceGl,10.*CompGl,SIceGl,
     1	 NIce,SDEDY,SDEYR,id
 111     FORMAT (1X,F9.2,2X,I3,2X,3f7.2,2x,i3,2x,i3,2x,f4.0,2x,a4)
	ip=34
	jp=45
         WRITE (34,112) MONTH,GHIce(ip,jp), GHIce2(ip,jp),
     1     10.*GIceComp(ip,jp),GT1(ip,jp), GT2(ip,jp),
     2     SnoAmt1(ip,jp), SnoAmt2(ip,jp),SDEDY,isftyp(ip,jp),SDEYR,id
 112     FORMAT (1X,I2,2X,7f7.2,2x,i3,2x,i1,2x,f4.0,2x,a4)
	ip1=10
	jp1=44
	ip2=58
	jp2=44
      WRITE(35,113)
	1    MONTH, GHIce(ip1,jp1), SnoAmt1(ip1,jp1),isftyp(ip1,jp1),
	2	       GHIce(ip2,jp2), SnoAmt1(ip2,jp2),isftyp(ip2,jp2),
     3           SDEDY,SDEYR,id
 113     FORMAT (1X,I2,2X,2(2f7.2,2x,i2),2x,i3,2x,f4.0,2x,a4)
cIIIIIIIIIIIIIIIIIIIIIIIIIII
       WRITE (61,126) TAU,MONTH,(xIce(i,1)/numgl1,i=1,10),SDEDY,SDEYR,id
 126       FORMAT (1X,F9.2,2X,I3,2X,10f10.4,2x,i3,2x,f5.0,2x,a4)
	xIce=0.
cIIIIIIIIIIIIIIIIIIIIIIIIIII

	 if (mnthdy.eq.1) then

        call OutGRADS1(gt,74,46,273.d0,1.d0,'monthly')
        call OutGRADS1(GHIce,74,46,0.d0,1.d0,'monthly')
        call OutGRADS1(GIceComp,74,46,0.d0,1.d0,'monthly')
c        call PrnMap(1,p,74,46,800.d0,1.d0,' P ','P.map')

        call PrnMap(1,gt,74,46,273.d0,1.d0,' GT ','GT.map')
c        call PrnMap(1,ts,74,46,273.d0,1.d0,' TS ','TS.map')
        call PrnMap(1,GHIce,74,46,0.d0,1.d0,'GHIce','GHIce.map')
c       call PrnMap(1,GIceComp,74,46,0.d0,1.d-2,'Compactness',
c     * 'COMP.map')
	 endif
      endif
c         Next day integration
           go to 10
 606   continue       !!!!!!!! END of RUN!!!!!!
c      Current time
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)
          write (6,210) ' Current time is ',      !Display
     1    tmphour, tmpminute, tmpsecond
          write (7,210) ' Current time is ',      !File
     1    tmphour, tmpminute, tmpsecond
  210      format(a,I2,':',I2.2,':',I2.2)
c          Time elapsed
c          Print time elapsed (h,min,sec):
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund) 
      time1 =tmphour*3600.+ tmpminute*60.+tmpsecond+tmphund*0.01
      sec=time1-time0
	     time0=time1
       if (sec.lt.0.) sec=sec+24.*3600.
	secday=sec/(tau-taue)*24.
       
           ihours=idint(sec/3600.0d0)
           jmin=idint((sec-ihours*3600)/60.0d0)
           sec1=sec-dfloat(ihours*3600+jmin*60)
           write(6,11) ihours,jmin,sec1,secday
           write(7,11) ihours,jmin,sec1,secday
c          Play music when stop
        CALL BEEPQQ(2000, 200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(1000,200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(500, 200)
c
      STOP !'Normal end'
      END
