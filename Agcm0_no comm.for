
C     ************
      PROGRAM AGCM
C     ************
C             2 - LEVEL AGCM
C             5x4 DEGREES VERSION
C************ BEGINNING OF COMMON ************
C
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
     1 ,GT_temp1(1:74,1:46)
     
      !// dummy integer  k2(0:72+1,0:72+1)
      
      !GLD GGGGGGGGG 
c      Dialogue and Setup  
	 BaseDir = '..\'
       call Intro(istop,rest_off,restfile,restfile_ice)
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
       TAU=24.
       TAUE=TAU
c      Time step dependent variables
c  !!!!
	 DT=  300.	!    480.!600. 
       DTC3=DT*6.
       TSPD=24.*3600./DTC3
       COE1=COE*DTC3/(24.0*3600.)
        day10=0
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
!         	call SurFile ('GT_Atm',GT,273.1E0,'GT_Atm')
!       stop 'c3'

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
        CALL SDET         !Solar  atm
c********** XOCEAN ***
        CALL INTSEA       !INTERPOLATE SEA TEMPERATURES
!GLD GGGGGGGGG
        sdedy_oc=-1
        !syncronize AGCM and GLD dates before run
        print *, 'Syncronization...'
         do while (sdedy_oc/=mod(sdedy+1,365))
          call gldstn(first_gldstn, have_TS_atm,
     1           sdedy_oc,ts_oc_for_atm)
	    first_gldstn = .false.      !not first call gldstn
         enddo
	    have_TS_atm  = .true.  !if use TS_atm for GLD atm
         GT_temp1=GT
!INTERPOLATE ts_oc_for_atm(72,72) on AGCM grid (72,46)
!GT_temp(74,46) = ts_oc_for_atm after interpolation 
!         call Interpol_oc_atm(ts_oc_for_atm, GT_temp, ISFTYP, k2) !old 
         write(1446,'(74F10.4)') ((ts_oc_for_atm(i,j),i=0,73),j=73,0,-1)
        call interpolate(ts_oc_for_atm,GT_temp)     !Tagir
!        call interpolate(ts(1,:,:,kmax),t_AGCM)    !Tagir
         write(1447,'(74F10.4)') ((GT_temp(i,j),i=1,74),j=46,1,-1) !grad C
         do i=2,73
            do j=1,46
             if ((isftyp(i,j).eq.7).and.(GT(i,j).gt.273.1)) then
                   GT(i,j)=GT_temp(i,j)+273.1 
             endif
            enddo
         enddo     
    !    	call SurFile ('GT_Oc',GT,273.1E0,'GT_Oc') !GT - from GLD
    !   	call SurFile ('GT_a-o',GT_temp1-GT,0.E0,'GT_a-o4')
 !GGGG       Stop 'GT'
 !        open (unit=229,file="GT_a-o.txt")
 !        write(229,'(72F10.4)') ((GT_temp1(i,j)-GT(i,j),i=2,73),j=1,46)
 !        stop  ' GT_temp '      
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
C
      CALL STEP
    
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
 
C
C         ADD FORCING TERMS
C
       CALL FORCE
!GLD GGGGGGGGG 
!TS(74,46) - surface air temp in AGCM (Kelvin grad)
! Q4_AGCM(74,46) - surface air hum in AGCM
! TOTALP_AGCM(74,46) - prec in AGCM (g/cm**2)
       n_acc=n_acc+1
 !      TS_AGCM(2:73,1:46) = TS_AGCM(2:73,1:46)+(T(2:73,1:46,2)
 !    1  +TS(2:73,1:46))/2.
      do J=1,46
      do I=2,73
        TS_AGCM(I,J) = TS_AGCM(I,J)+TS(I,J)
   !     TS_AGCM(2:73,1:46) = TS_AGCM(2:73,1:46)+TS(2:73,1:46)
      enddo
      enddo
       QS_AGCM(2:73,1:46) = QS_AGCM(2:73,1:46)+Q4_AGCM(2:73,1:46)
       PREC_AGCM(2:73,1:46) = PREC_AGCM(2:73,1:46)+
     1                                       TOTALP_AGCM(2:73,1:46)
       IF (end_of_day) then !day mean values
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
!INTERPOLATE AGCM(72,46) to _atm_for_oc (72,72)
   !    call Interpol_atm_oc(TS_AGCM, TS_atm_for_oc)  !old
   !    call Interpol_atm_oc(QS_AGCM, QS_atm_for_oc)  
   !    call Interpol_atm_oc(PREC_AGCM, PREC_atm_for_oc)  
   !        write(1449,'(72F10.4)') ((tq(1,i,j),i=1,72),j=1,72)
      TS_AGCM=TS_AGCM-273.1
      call interpolate_back(TS_AGCM, TS_atm_for_oc) !  shift
      call interpolate_back(QS_AGCM, QS_atm_for_oc) !  shift
      call interpolate_back(PREC_AGCM, PREC_atm_for_oc) !  shift
  !     TS_atm_for_oc=TS_atm_for_oc-273.1  ! Kelvin to celsius
       write(1451,'(74F10.4)') ((TS_atm_for_oc(i,j),i=0,73),j=0,73)
       write(1452,'(74F10.4)') ((TS_AGCM(i,j),i=1,74),j=1,46)
       write(1453,'(74F10.4)') ((T(i,j,2)-273.1,i=1,74),j=1,46)
       write(1454,'(74F10.4)') ((TS(i,j)-273.1,i=1,74),j=1,46)
 !      stop 'TS_atm_for_oc'
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
       IF (end_of_day.and.mnthdy.eq.1) then
              CALL PRZON(11)     ! PRINT OUT MONTHLY ZONAL MEANS
         CALL acc0
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
         GT_temp1=GT
 
       call gldstn(first_gldstn, have_TS_atm, sdedy_oc,ts_oc_for_atm )
     
   !    open (unit=224,file="ts_oc_for_atm2.txt")
    !   write(224,'(72F10.4)') ((ts_oc_for_atm(i,j),i=1,72),j=1,72)
!       open (unit=225,file="GT2_1.txt")
!       write(225,'(72F10.4)') ((GT(i,j),i=2,73),j=1,46)
        
!INTERPOLATE ts_oc_for_atm (72,72) to GT(72,46) 
   !    call Interpol_oc_atm(ts_oc_for_atm, GT_temp, ISFTYP, k2) 
         call interpolate(ts_oc_for_atm,GT_temp)   !Tagir
        do i=2,73
            do j=1,46
             if ((isftyp(i,j).eq.7).and.(GT(i,j).gt.273.1)) then
               if (sdedy==16) then
        !       write (113,*) i,j,GT(i,j)-GT_temp(i,j)-273.1 !GGGGGGGGGG
               endif
               GT(i,j)=GT_temp(i,j)+273.1 !TTT
            endif
            enddo
         enddo
 !     if (sdedy==16) write(1448,'(72F10.4)')((GT_temp(i,j)
 !    1 ,i=2,73),j=1,46)
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

       endif !if end_of_day

c         Next day integration
      go to 10
 606  continue       !!!!!!!! END of RUN!!!!!!
c
      STOP 'Normal end'
      END
