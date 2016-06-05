c       INCLUDE   'FLIB.FI'

C     ************
      PROGRAM AGCM
C     ************
C             2 - LEVEL AGCM
C             5x4 DEGREES VERSION
C************ BEGINNING OF COMMON ************
C
!	USE MSFLIB
        USE IFCORE   !for stop key
C
	include 'recom.fi'
 !     DIMENSION C(900)
 !     EQUIVALENCE (TAU,C(1))

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
      dimension    RISFTYP(74,46)
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
      ! COMMON /COMP/ DUM(74,46,4) !- DUM - not used
      include 'comp.fi'
C
      COMMON /WORK/ DUM2(72,46,9)
C
      COMMON /SEA/ SS(74,46,12), SWIN(74,46,12)
      
      common/glacc/GMTACC,GMKEAC,GMRACC,namgl 
c
	include 'dir.fi'
C
	include 'ice.fi'
      include 'acc_ice.fi'
C
cIIIIIIIIIIIIIIIIIIIIIIIIIII
      include 'StatIce.fi'
cIIIIIIIIIIIIIIIIIIIIIIIIIII
	!include 'varAGCM.f90' 
	     real  TS_atm_for_oc(0:72+1,0:72+1), QS_atm_for_oc(0:72+1,0:72+1),
     1         PREC_atm_for_oc(0:72+1,0:72+1)  
        common /varsAGCM/  TS_atm_for_oc, QS_atm_for_oc, PREC_atm_for_oc
        
         COMMON /OSA/ QSN(72,46),TXN(72,46),TYN(72,46)

C
C+++
c      INCLUDE 'flib.fd'

c        For program execution time
      INTEGER*2 tmphour, tmpminute, tmpsecond, tmphund
      COMMON /day_end/ end_of_day                           
      logical stop_time,end_of_day,stop_key,rest_off,first
      logical  first_gldstn, have_TS_atm   !GLD GGGGGGGGG 
      character*40 restfile
      character*60 restfile_ice
      INTEGER day10, iTime,NF, sdedy_oc
      !GLD GGGGGGGGG 
      ! SST from GLDSTN, air temp from AGCM 
      real ts_oc_for_atm(0:72+1,0:72+1),TS_AGCM(1:74,1:46),
     1     QS_AGCM(1:74,1:46),PREC_AGCM(1:74,1:46),GT_temp(1:74,1:46)
     
      integer  k2(0:72+1,0:72+1)
      
      !GLD GGGGGGGGG 
      character*2 fileN(40)/'01','02','03','04','05','06','07','08',
     1 '09','10','11','12','13','14','15','16','17','18',
     2 '19','20','21','22','23','24','25','26','27','28','29','30',
     3 '31','32','33','34','35','36','37','38','39','40'/
c         ABS(XXX)=DABS(XXX)
C
c      Dialogue and Setup  
	 BaseDir = '..\'
!	 BaseDir = 'E:\ClimParal\AGCM-ICE3\'  !'e:\climdig\AGCM_DLG\' 
 !      OPEN (333, FILE = 'nul', status='new')
           NF=1 !For Usatjuk   UUUUUUUUU
 !          OPEN (333, FILE = 'Month'//fileN(NF)//'.txt', status='new')
            NF=NF+1
 !      OPEN (333, FILE = 'Month10-12.txt', status='new')
       call Intro(istop,rest_off,restfile,restfile_ice)
c	call printc
c          Initial time
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund) !initial time
      time0 =tmphour*3600.+ tmpminute*60.+tmpsecond+tmphund*0.01 ! (sec)

	 Ibreak = 3
	 Jbreak = 31
           stop_time=.false.
           end_of_day=.false.
	     stop_key=.false.
           first=.true.          ! for stop key
	     nstop=0               ! stop key print counter
	     first_gldstn = .true. !first call gldstn GLD GGGGGGGGG
	     have_TS_atm = .false.   !first call gldstn GLD GGGGGGGGG
        write (*,*)' *****Program  Execution Began...    **************'
	write (*,*)' *****To Stop Running Program - Press ANY KEY *****'
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
	 DT=  300.	!    480.!600. 
       DTC3=DT*6.
       TSPD=24.*3600./DTC3
       COE1=COE*DTC3/(24.0*3600.)
        day10=0
c         TDSEA=.FALSE.
c         SWAMP=.FALSE.
c         OIL  =.FALSE.
c      call PrintC
c      write (6,91) ((isftyp(i,j),i=2,73),j=46,1,-1)
c 91   format (//(72i1))            
C
c        COMPLETE TIME STEP USING COMP3

	 write (7,'(9I8)') ,(i, i=1,9) !Print surface albedo
	do j=1,2  
	 write (7, '(9e8.2)'), (sfcalb(i,j), i=1,9)
	enddo
      do i=1,im+2   !limit great snow
       do J=1,jm
        if (SNOAMT(i,J) .gt. 100.d0) SNOAMT(i,J)=100.d0
        if (SNOAMT1(i,J) .gt. 100.d0) SNOAMT1(i,J)=100.d0
       enddo
      enddo
c	RISFTYP=real(ISFTYP)
c      CALL WRITIJ(RISFTYP,1.d0,'Type',0)
      DO 1 J=1,JM
  1   CALL COMP3 (J)
c      CALL WRITIJ(SNOAMT,1.d0,'SNOW',0)  !1.d-4

c      CALL WRITIJ(WriteArr,1.d0,'snr1',0)
c       go to 2222
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
c      CALL WRITIJ(TS    ,1.d0,'TS  ',0) 
!         	call SurFile ('GT_Atm',GT,273.1E0,'GT_Atm')
!       stop 'c3'

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
      WRITE(7,803) DT,DTH,TAUH,TAU,TAUI,TOFDAY,SDEYR,TSPD,fm
      WRITE(7,805) FixedWater,FixedIceBorder,iDaysToFreeze,
	1	    SFCALB(9,1),GIceUpFlow(1,1)

      WRITE(*,804)NSTEP,NGMP,NSDET,NRES,IM,JM,NCYCLE,NC3,SDEDY,NAV
      WRITE(*,803) DT,DTH,TAUH,TAU,TAUI,TOFDAY,SDEYR,TSPD,fm
 804  FORMAT (' *STEP =',I7,' *NGMP =',I7,' *NSDET=',I7,
     1   ' *NRES =',I7,' *IM   =',I7/' *JM   =',I7,' *NCYCL=',I7,
     2  ' *NC3  =',I7,' *SDEDY=',I7,' *NAV  =',I7)
 803  FORMAT (' *DT  =',F8.2,' *DTH =',F8.2,' *TAUH=',F8.2,
     1' *TAU =',F8.2,' *TAUI=',F8.2,/' *TOFDAY=',F8.2,
     2' *SDEYR =',F8.2,' *TSPD  =',F8.2,' *fm=',d12.3/)
 805  FORMAT (' *FixedWater  =',L1,' *FixedIceBorder =',L1,
     1 ' *iDaysToFreeze=',I2,/' *SFCALB =',F6.2,' *GIceUpFlow=',F6.2,/)
     
!GLD GGGGGGGGG
      IF (end_of_day) then
        CALL SDET         !Solar
c********** XOCEAN ***
        CALL INTSEA       !INTERPOLATE SEA TEMPERATURES
!GLD GGGGGGGGG
        sdedy_oc=-1
         do while (sdedy_oc/=mod(sdedy+1,365))    
      call gldstn(first_gldstn, have_TS_atm, sdedy_oc,ts_oc_for_atm, k2)
	      first_gldstn = .false.      !not first call gldstn
         enddo
	   have_TS_atm  = .true.
	   !open (unit=222,file="ts_oc_for_atm.txt")
         !write(222,'(72F10.4)') ((ts_oc_for_atm(i,j),i=1,72),j=1,72)
         open (unit=223,file="GT1_1.txt")
         write(223,'(72F10.4)') ((GT(i,j),i=2,73),j=1,46)
!INTERPOLATE ts_oc_for_atm(72,72) on AGCM grid (72,46)
!GT(74,46) = ts_oc_for_atm after interpolation 
         call Interpol_oc_atm(ts_oc_for_atm, GT_temp, ISFTYP, k2) 
         do i=2,73
            do j=1,46
             if ((isftyp(i,j).eq.7).and.(GT_temp(i,j).gt.273.1)) then
                   GT(i,j)=GT_temp(i,j) !TTT
             endif
            enddo
         enddo     
         GT_temp=0.
    !    	call SurFile ('GT_Oc',GT,273.1E0,'GT_Oc')
 !GGGG       Stop 'GT'
         open (unit=229,file="GT1_2.txt")
         write(229,'(72F10.4)') ((GT(i,j),i=2,73),j=1,46)
!GLD GGGGGGGGG
	  if (.not.FixedWater) then
         call ocean
         QSN=0.d0           ! ZERO SEA HEAT EVERY DAY
	  endif
c********************
        call IceDaily     ! once-per-day ice calculations
        CALL INTOZ        ! ozon
        CALL GMP          ! mass check  
      endif
C
          open(111,ACCESS ='APPEND')  !global values agcm

 10   continue	!next  time step of AGCM
 
      stop_key=PEEKCHARQQ()
c      stop_key=.FALSE.
      if (stop_key.and.first) then
       write (*,*)' ****** STOP AFTER FULL DAY. WAIT ...*****'
       first=.false.
      endif 
      
      NSTEP=NSTEP+NC3    !NC3=6
      end_of_day= (MOD(NSTEP,NSDET).EQ.0)
      TAU=(NSTEP*DT)/3600.0
      TOFDAY=DMOD(TAU-1.0,ROTPER)+1.0
      ROT=TOFDAY/ROTPER*(PI+PI)
      COSR=DCOS(ROT)
      SINR=DSIN(ROT)
C
C         INTEGRATE FORWARD
c      WRITE(*,*) 'nstep', nstep
C
      CALL STEP
c       stop 'step'

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
c       if (SDEYR.ge.9.9.and.month.eq.3) stop_time=.true. !!!!!Temporal!!!!
      endif

      if (rest_off) goto 9997  !No Restart records to file
      IF (end_of_day.and.(mnthdy.eq.1.or.stop_time)) then
        !Restart AGCM records to file
       CALL OUTRES(restfile,9,restfile_ice,10) 
             open(2,file= '..\..\gldstn_F90_1980\results\w29.3') 
            rewind 2
	      ! data for GLDSTN restart
            ! write out: ts(l,i,j,k), u(l,i,j,k),tq(l,i,j),varice(l,i,j),
                                  !tice(i,j),t after any iwstp steps
            call outm(2)
            close(2)
       endif     
9997  continue
 
!UUUUUUU For CO2 Usatjuk Dec 10
      go to 9998
       if (sdedy==1.and.end_of_day) then 
            close (333)
           OPEN (333, FILE = 'Month'//fileN(NF)//'.txt', status='new')
            NF=NF+1
       endif     
       if (sdedy==91.and.end_of_day) then
            close (333)
           OPEN (333, FILE = 'Month'//fileN(NF)//'.txt', status='new')
            NF=NF+1
       endif     
       if (sdedy==182.and.end_of_day) then
            close (333)
           OPEN (333, FILE = 'Month'//fileN(NF)//'.txt', status='new')
            NF=NF+1
       endif     
       if (sdedy==274.and.end_of_day) then
            close (333)
           OPEN (333, FILE = 'Month'//fileN(NF)//'.txt', status='new')
            NF=NF+1
       endif
 9998   continue         
!UUUUUUU For CO2 Usatuk Dec 10
C
C         ADD FORCING TERMS
C
     
       CALL FORCE
!GLD GGGGGGGGG 
!TS(74,46) - surface air temp in AGCM
! Q4_AGCM(74,46) - surface air hum in AGCM
! TOTALP_AGCM(74,46) - prec in AGCM (g/cm**2)
       n_acc=n_acc+1
       TS_AGCM(2:74,1:46) = TS_AGCM(2:74,1:46)+TS(2:74,1:46)
       QS_AGCM(2:74,1:46) = QS_AGCM(2:74,1:46)+Q4_AGCM(2:74,1:46)
       PREC_AGCM(2:74,1:46) = PREC_AGCM(2:74,1:46)+
     1                                       TOTALP_AGCM(2:74,1:46)
       IF (end_of_day) then
       TS_AGCM = TS_AGCM/float(n_acc)
       QS_AGCM = QS_AGCM/float(n_acc)
       PREC_AGCM = PREC_AGCM/float(n_acc)*0.01/DTC3  !prec (m/sec)
      DO  J=1,JM
       TS_AGCM(1,J)=TS_AGCM(73,J)
       TS_AGCM(74,J)=TS_AGCM(2,J)
       QS_AGCM(1,J)=QS_AGCM(73,J)
       QS_AGCM(74,J)=QS_AGCM(2,J)
       PREC_AGCM(1,J)=PREC_AGCM(73,J)
       PREC_AGCM(74,J)=PREC_AGCM(2,J)
      enddo
!       open (unit=217,file="TQ_AGCM.txt")
    !     	call SurFile ('TS_AGCM',TS_AGCM,273.1E0,'TS_AGCM')
    !     	call SurFile ('QS_AGCM',QS_AGCM,0.,'QS_AGCM')
    !     	call SurFile ('PREC_AGCM',PREC_AGCM,0.,'PREC_AGCM')
!       write(217,'(72F10.4)') ((TS_AGCM(i,j)-273.1,i=1,72),j=1,46)
!       write(217,'(72F10.4)') ((QS_AGCM(i,j),i=1,72),j=1,46)
!INTERPOLATE TS_AGCM(2,72,46) to TS_atm_for_oc (2,72,72)
       call Interpol_atm_oc(TS_AGCM, TS_atm_for_oc)  
       TS_atm_for_oc=TS_atm_for_oc-273.1E0  ! Kelvin to celsius
       call Interpol_atm_oc(QS_AGCM, QS_atm_for_oc)  
       call Interpol_atm_oc(PREC_AGCM, PREC_atm_for_oc)  
 !      open (unit=218,file="TS_atm_for_oc.txt")
 !      write(218,'(72F10.4)') ((TS_atm_for_oc(i,j),i=1,72),j=1,72)
 !     CALL WRITIJ(PREC_AGCM,1.d-5,'pr00',0)
 !     Stop 'TS_AGCM'
       TS_AGCM = 0.
       QS_AGCM = 0.
       PREC_AGCM = 0.
       n_acc=0
      endif
 !GLD GGGGGGGGG

c       CALL outac1(42)
c       STOP 'Normal end'
C     OCEAN MODEL
C###      IF (end_of_day.AND.TDSEA) CALL OCEAN
C
c      IF (mod(nstep,nhis).eq.0) then
      IF (end_of_day) then
        WRITE (7,999) MONTH,MNTHDY,TOFDAY,SDEYR,GMKE,GMT,GMR
        WRITE (*,999) MONTH,MNTHDY,TOFDAY,SDEYR,GMKE,GMT,GMR
        WRITE (111,101) MONTH,MNTHDY,TOFDAY,SDEYR,GMKE,GMT,GMR
 999    FORMAT (1X,I3,1H/,I2,1H/,F6.2,2X,5Hyear=,f5.0,2X,5HK.E.=,
     1   1PD10.3,3H T=,1PD10.3,4H N0=,1PD10.3)
101    FORMAT (1X,2I3,5f10.3)

        nstop=nstop+1    ! Print announsement to screen
        if (nstop.eq.21) then
         write (*,*)' *****To Stop Running Program - Press ANY KEY*****'
	   nstop=0
        endif
      endif
C
C         OUTPUT ACCUMULATED VARIABLES HISTORY
C
      if (rest_off) goto 123
       IF (end_of_day.and.(mnthdy.eq.1.or.stop_time)) then
        CALL OUTACC(43,0)  ! Accumulated variables for restart
c                 CALL OutAcIce(45,0)
    !    CALL INACC(43)
       ENDIF
 123  continue

C         OUTPUT OCEAN VARIABLES
C**      IF (end_of_day.AND.TDSEA) CALL OUTOCN
C
C         UPDATE DAY, EARTH-SUN DISTANCE, SOLAR DECLINATION
C
      IF (end_of_day) CALL SDET   !solar

C        OUTPUT MONTHLY ACCUMULATED VARIABLES HISTORY

      if (rest_off) goto 125
        IF ((MONTH.EQ.2.OR.MONTH.EQ.8).and.! January and July output
     *     (mnthdy.eq.1).and.end_of_day)  then
c              CALL PRZON(11)     ! PRINT OUT MONTHLY ZONAL MEANS
              call outacc(42,1)   ! Accumulated variables
c             call maps           ! Print out maps
        endif
        IF ((mnthdy.eq.1).and.end_of_day) then ! Every month output
c          call OutAcIce(44,1)   ! Accumulated Ice variables
        endif
125   continue

	numgl1=namgl !force per day

c      if (end_of_day)  call GlobAcc(13)       ! Print global values

c     Put zero accumulations in the end of month
!       IF (end_of_day) then
       IF (end_of_day.and.mnthdy.eq.1) then
              CALL PRZON(11)     ! PRINT OUT MONTHLY ZONAL MEANS
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
           print 11,ihours,jmin,sec1,secday
           write(7,11) ihours,jmin,sec1,secday
 11        format (1x,'Time elapsed: ',i2,'h ',i2,'min ',f5.2,'sec',
	1             '   Sec/day: ',f5.2)
       endif
	
            !########### 1 January stop time
         IF     ((MONTH.EQ.1).and.      
     *   (mnthdy.eq.2).and.end_of_day.and.NF==40)    stop_time= .true.
        if (stop_time) go to 606

        IF (end_of_day) then
c********** XOCEAN ***
         CALL INTSEA       ! INTERPOLATE SEA TEMPERATURES
!GLD GGGGGGGGG
       ! pause
 
       call gldstn(first_gldstn, have_TS_atm, sdedy_oc,ts_oc_for_atm,k2)
     
   !    open (unit=224,file="ts_oc_for_atm2.txt")
    !   write(224,'(72F10.4)') ((ts_oc_for_atm(i,j),i=1,72),j=1,72)
!       open (unit=225,file="GT2_1.txt")
!       write(225,'(72F10.4)') ((GT(i,j),i=2,73),j=1,46)
        
!INTERPOLATE ts_oc_for_atm (72,72) to GT(72,46) 
       call Interpol_oc_atm(ts_oc_for_atm, GT_temp, ISFTYP, k2) 
         do i=2,73
            do j=1,46
             if ((isftyp(i,j).eq.7).and.(GT_temp(i,j).gt.273.1)) then
               if (sdedy==16) then
                write (113,*) i,j, GT(i,j)-GT_temp(i,j)!GGGGGGGGGG
               endif
               GT(i,j)=GT_temp(i,j) !TTT
            endif
            enddo
         enddo
 !        if (sdedy==16) stop 'Jn 16'     
       GT_temp=0.     
!       open (unit=227,file="GT2_2.txt")
!       write(227,'(72F10.4)') ((GT(i,j),i=2,73),j=1,46)
         
!GLD GGGGGGGGG
	   if (.not.FixedWater) then
          call ocean
          QSN=0.d0           ! ZERO SEA HEAT EVERY DAY
	   endif
c********************
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
	     HIceGl=HIceGl+GHIce(i,j)*GIceComp(i,j)*dxyp(j)
	     CompGl=CompGl+GIceComp(i,j)*dxyp(j)
	     SIceGl=SIceGl+dxyp(j)
	    endif
	   enddo
	  enddo
	  if (NIce.eq.0) NIce=1
	  HIceGl=HIceGl/CompGl
	  CompGl= CompGl/SIceGl
	  SIceGl=SIceGl/1.d12
        WRITE (33,111) TAU,MONTH,HIceGl,10.*CompGl,SIceGl,
     1	 NIce,SDEDY,SDEYR,id
 111     FORMAT (1X,F9.2,2X,I3,2X,3f7.2,2x,i3,2x,i3,2x,f3.0,2x,a4)
	  ip=34
	  jp=45
        WRITE (34,112) MONTH,GHIce(ip,jp), GHIce2(ip,jp),
     1     10.*GIceComp(ip,jp),GT1(ip,jp), GT2(ip,jp),
     2     SnoAmt1(ip,jp), SnoAmt2(ip,jp),SDEDY,isftyp(ip,jp),SDEYR,id
 112     FORMAT (1X,I2,2X,7f7.2,2x,i3,2x,i1,2x,f3.0,2x,a4)
	  ip1=10
	  jp1=44
	  ip2=58
	  jp2=44
        WRITE(35,113)
	1    MONTH, GHIce(ip1,jp1), SnoAmt1(ip1,jp1),isftyp(ip1,jp1),
	2	       GHIce(ip2,jp2), SnoAmt1(ip2,jp2),isftyp(ip2,jp2),
     3           SDEDY,SDEYR,id
 113     FORMAT (1X,I2,2X,2(2f7.2,2x,i2),2x,i3,2x,f3.0,2x,a4)
cIIIIIIIIIIIIIIIIIIIIIIIIIII
       WRITE (61,126) TAU,MONTH,(xIce(i,1)/numgl1,i=1,10),SDEDY,SDEYR,id
 126       FORMAT (1X,F9.2,2X,I3,2X,10f10.4,2x,i3,2x,f5.0,2x,a4)
	  xIce=0.
cIIIIIIIIIIIIIIIIIIIIIIIIIII

	  if (mnthdy.eq.1) then

c        call OutGRADS1(gt,74,46,273.d0,1.d0,'monthly')
c        call OutGRADS1(GHIce,74,46,0.d0,1.d0,'monthly')
c        call OutGRADS1(GIceComp,74,46,0.d0,1.d0,'monthly')
c        call PrnMap(1,p,74,46,800.d0,1.d0,' P ','P.map')

c        call PrnMap(1,gt,74,46,273.d0,1.d0,' GT ','GT.map')
c        call PrnMap(1,ts,74,46,273.d0,1.d0,' TS ','TS.map')
c        call PrnMap(1,GHIce,74,46,0.d0,1.d0,'GHIce','GHIce.map')
c       call PrnMap(1,GIceComp,74,46,0.d0,1.d-2,'Compactness',
c     * 'COMP.map')
	  endif
        day10=day10+1
       endif !if end_of_day
      go to 10

      IF (end_of_day.and.day10.eq.5) then
       CALL WRITIJ(SNOAMT,1.d0,'SNOW',0)  !1.d-4
       CALL WRITIJ(PHIS,1.d0,'PHIS',0)
       CALL WRITIJ(P   ,1.d0,'P   ',0)
       CALL WRITIJ(U(1,1,1),1.d0,'U1  ',2)
       CALL WRITIJ(U(1,1,2),1.d0,'U3  ',2)
       CALL WRITIJ(V(1,1,1),1.d0,'V1  ',2)
       CALL WRITIJ(V(1,1,2),1.d0,'V3  ',2)
       CALL WRITIJ(T(1,1,1),1.d0,'T1  ',0)
       CALL WRITIJ(T(1,1,2),1.d0,'T3  ',0)
       CALL WRITIJ(QW(1,1,1),1.d-4,'QW1 ',0)
       CALL WRITIJ(QW(1,1,2),1.d-4,'QW3 ',0)
       CALL WRITIJ(GW  ,1.d-4,'GW  ',0)
       CALL WRITIJ(GT  ,1.d0,'GT  ',0)
       CALL WRITIJ(TS    ,1.d0,'TS  ',0)
       stop  '1day'
      end if
c         Next day integration
      go to 10
 606  continue       !!!!!!!! END of RUN!!!!!!
c      Current time
      CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)
          write (*,210) ' Current time is ',      !Display
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
           print 11,ihours,jmin,sec1,secday
           write(7,11) ihours,jmin,sec1,secday
c          Play music when stop
        CALL BEEPQQ(2000, 200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(1000,200)
        CALL SLEEPQQ(100)
        CALL BEEPQQ(500, 200)
c
      STOP 'Normal end'
      END
