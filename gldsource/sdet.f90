      !###############
      !***************
      SUBROUTINE SDET_OC(SCOSZ)
      !***************
      !***************
      !             UPDATES DAY, EARTH-SUN DISTANCE, AND SOLAR DECLINATION
      !************ BEGINNING OF COMMON ************
      include 'var.f90'
      !
      !************ END OF COMMON ******************
      !
      integer MONTHS(12)
      integer L
      integer MAXDAY, JDYACC,MONTH,j0,j
      real DY, SEASON, DIST, DEC,RAD,SOLTCE,APHEL,temp,RSDIST,&
            SIND,COSD,DECMAX,ECCN,DAYPYR,ROTPER,S0,fi0,fi,coe,t00
	real SCOSZ(maxj) !mean day insolation W/m**2

      DATA  MONTHS/31,28,31,30,31,30,31,31,30,31,30,31  /

      RAD= 6375000.0 !m
	SOLTCE=	173.0  !day
	APHEL=	182.0   !day
	DECMAX=	0.41015 !max declin (rad)
	ECCN=	1.7799999E-002
	MAXDAY=	365   !day
	DAYPYR=	365.0  !day
	ROTPER=	24.0  !h

      MAXDAY=DAYPYR + 1.0E-2
      SDEDY=SDEDY+1

      IF (SDEDY .GT. MAXDAY) THEN
       SDEDY=SDEDY-MAXDAY
       SDEYR=SDEYR+1.0
      END IF
      JDYACC=0
      DO L=1,12
       MONTH = L
       JDYACC=JDYACC+MONTHS(L)
       IF (SDEDY .LE. JDYACC) EXIT
      END DO
      MNTHDY=MONTHS(L)-JDYACC+SDEDY
      DY=SDEDY
      SEASON=(DY-SOLTCE)/DAYPYR
      DIST=(DY-APHEL )/DAYPYR
!
!        SOLTCE = JUNE 22
!        APIHELION = JULY 1
!        ECCN= ORBITAL ECCENTRICITY
!
      DEC=DECMAX*COS(2.0*PI*SEASON) !declin
      RSDIST=(1.0+ECCN*COS(2.0*PI*DIST))**2
      S0=1368./RSDIST   !W/m**2
      SIND=SIN(DEC)
      COSD=COS(DEC)
       !not used
	fi0=-pi/2.+abs(dec) !lat, where t00=12 or t00=0  -light day duration
	j0=1
       !not used
      do j=0,jmax-1
	 fi=sin(fi0)
	 if (fi.ge.s(j).and.fi.lt.s(j+1)) then
	  j0=j+1 !lat, where t00=12 or c=0
	  exit
	 endif
	enddo

	 SCOSZ=0.  !total day insol (W/m**2)
	 coe=0.5*rotper/pi

      do j=1,jmax
       temp=s(j)*rc(j)*sind/cosd
	 if (temp.le.-1.0) temp=-1.0  !t00 = 0 h
	 if (temp.ge. 1.0) temp= 1.0  !t00 = 12 h
	  t00=coe*acos(-temp) !solar day duration (h): 0<t00<12
        SCOSZ(j)=2.*S0*(t00*S(J)*SIND+coe*C(J)* COSD*sqrt(1.-temp*temp))/24.
       if (SCOSZ(j).LE.0.0) SCOSZ(j)=0.0
	enddo
!####################
      RETURN
      END
