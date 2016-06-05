      SUBROUTINE IceDaily()

C************ BEGINNING OF COMMON ************
	include 'recom.fi'
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
     1  ,SD(74,46)
C
C        Q ARRAY - STATE VARIABLES
C
      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)

      include 'ice.fi'
C************ END OF COMMON ******************
      real*8	IcePoolLimit, FullHeatng, Ice2PoolLimit, FullHeatng2

c      call map(1, GHICE, 74, 46, 0d0, 1d0, 'ICE THICKNESS', '\ICE0.OUT')
c      call map(1, GIcePool, 74, 46, 0d0, 1d0, 'POOL', '\POOL.OUT')
c      call map(1, GIceCND, 74, 46, 0d0, 1d0, 'CND', '\CND.OUT')
c      call map(1, GT, 74, 46, TICE, 1d0, 'GT', '\GT.OUT')
c ÷икл
      do J=1,46
       do I=2,73
         if (I.eq.Ibreak.and.J.eq.Jbreak) then
          tmpdata = 0					
         endif
	  if (NewIce(I,J)) then  ! number of days when ts<tice
	   iGIceFreezeCount(I,J) = iGIceFreezeCount(I,J) + 1
	  else
	   iGIceFreezeCount(I,J) = 0
	  endif

        if (ISFTYP(I,J).eq.7.and.(.not.FixedIceBorder)
	1     .and.(iGIceFreezeCount(I,J).ge.iDaysToFreeze)) then
c        ќбразование льда на поверхности воды (тип C)
         ISFTYP(I,J)    = 9
         GIceComp(I,J)  = GIceCompMax
         GHIce(I,J)     = HIceMin
         GHIce2(I,J)    = 0d0
	   GIcePool(I,J)  = 0d0
	   GIce2Pool(I,J) = 0d0
	   GT(I,J)        = TICE
         GT1(I,J)       = TICE
         GT2(I,J)       = TICE
         SnoAmt1(I,J)	  = 0d0
         SnoAmt2(I,J)	  = 0d0
	  else
	   if (ISFTYP(I,J).eq.9) then
c слой 1 есть
	    IcePoolLimit   = GHIce(I,J)*QIce*GIcePoolMaxPercent
	    if (GIcePool(I,J).gt.IcePoolLimit) then
	     FullHeatng = GIceHeatng(I,J) + GIceUpFlow(I,J)
	1      + GIcePool(I,J) - IcePoolLimit
	     GIcePool(I,J) = IcePoolLimit
	    else
	     FullHeatng = GIceHeatng(I,J) + GIceUpFlow(I,J)
	    endif
	    if (FullHeatng.lt.0d0) then
	     GIcePool(I,J) = GIcePool(I,J) + FullHeatng
	     FullHeatng    = 0d0
	    endif
	    if (GIcePool(I,J).lt.0d0) then
	     FullHeatng    = GIcePool(I,J)
	     GIcePool(I,J) = 0d0
	    endif
	    if (GIcePool(I,J).eq.IcePoolLimit) then
           GHIce(I,J)    = GHIce(I,J) - FullHeatng/(QIce*0.7d0)
	     GIcePool(I,J) = GHIce(I,J)*QIce*GIcePoolMaxPercent
	    else
           GHIce(I,J)    = GHIce(I,J) - FullHeatng/QIce
	1       - GIceCND(I,J)/QWater
	    endif

	    if (GHIce2(I,J).ge.HIceMin) then
c слой 2 есть
	     Ice2PoolLimit   = GHIce2(I,J)*QIce*GIcePoolMaxPercent
	     if (GIce2Pool(I,J).gt.Ice2PoolLimit) then
	      FullHeatng2 = GIce2Heatng(I,J) + GIceUpFlow(I,J)
	1       + GIce2Pool(I,J) - Ice2PoolLimit
	      GIce2Pool(I,J) = Ice2PoolLimit
	     else
	      FullHeatng2 = GIce2Heatng(I,J) + GIceUpFlow(I,J)
	     endif
	     if (FullHeatng2.lt.0d0) then
	      GIce2Pool(I,J) = GIce2Pool(I,J) + FullHeatng2
	      FullHeatng2     = 0d0
	     endif
	     if (GIce2pool(I,J).lt.0d0) then
	      FullHeatng2     = GIce2Pool(I,J)
	      GIce2Pool(I,J) = 0d0
	     endif
	     if (GIce2Pool(I,J).eq.Ice2PoolLimit) then
            GHIce2(I,J)    = GHIce2(I,J) - FullHeatng2/(QIce*0.7d0)
	      GIce2Pool(I,J) = GHIce2(I,J)*QIce*GIcePoolMaxPercent
	     else
            GHIce2(I,J)    = GHIce2(I,J) - FullHeatng2/QIce
	1        - GIce2Cnd(I,J)/QWater
	     endif
	    else ! во второй зоне - вода, подтай с краев или новый лед
	     FullHeatng2 = GIce2Heatng(I,J) + GIceUpFlow(I,J)
	     if (FullHeatng2.gt.0d0) then
	      GIceComp(I,J) = GIceComp(I,J) - (1d0-GIceComp(I,J))
	1         *GIWExchage
	2         *FullHeatng2/(GHIce(I,J)*QIce - GIcePool(I,J))
	      if (GIceComp(I,J).gt.GIceCompMax) GIceComp(I,J)=GIceCompMax 
	      if (GIceComp(I,J).lt.0) then
		   GHIce(I,J) = 0d0 ! ушло все
	       GIceComp(I,J) = 0d0
	      endif
	     endif
           if ((TS(I,J).lt.(TICE-5.0d0).and.GHIce2(I,J).lt.HIceMin) 
	1		 .or.(J.eq.46.and.GIceCND(I,J).lt.0)) then
c           if (NewIce(I,J).or.(J.eq.46.and.GIceCND(I,J).lt.0)) then
c        ќбразование льда на поверхности воды (тип A)
           GHIce2(I,J)    = HIceMin
	      GIce2Pool(I,J) = 0d0
            GT2(I,J)       = TICE
            SnoAmt2(I,J)	 = 0d0
	     endif
	    endif

c теперь посмотрим, что вышло
          if (GHIce(I,J).ge.HIceThin.and.GHIce2(I,J).gt.HIceThin) then
c          перерос тонкий лед
c          (1-C)*(H2-Hthin) == dC*(H1-Hthin)   =>
	     GIceComp(I,J) = GIceComp(I,J) + (1d0-GIceComp(I,J))
	1       *(GHIce2(I,J)-HIceThin)/(GHIce(I,J)-HIceThin)
	     GHIce2(I,J)   = HIceThin
	     if (GIceComp(I,J).gt.GIceCompMax) then
c           dC*(H1-Hthin)	== Cmax*dH1   =>
c     пошлем нафиг излишки, а то все перерастет слишком
c	      GHIce(I,J)    = GHIce(I,J) + (GHIce(I,J)-HIceThin)*
c	1       (GIceComp(I,J)/GIceCompMax-1d0)
	      GIceComp(I,J) = GIceCompMax
	     endif
	      if (GIceComp(I,J).lt.0) then
		   GHIce(I,J) = 0d0 ! ушло все
	       GIceComp(I,J) = 0d0
	      endif
	    endif

          if (GHIce(I,J).lt.GHIce2(I,J)) then
c          тонкий лед и толстый мен€ютс€ местами
	     tmpdata        = GHIce2(I,J)
	     GHIce2(I,J)    = GHIce(I,J)
	     GHIce(I,J)     = tmpdata
	     tmpdata        = GIce2Pool(I,J)
	     GIce2Pool(I,J) = GIcePool(I,J)
	     GIcePool(I,J)  = tmpdata
	     tmpdata        = SnoAmt2(I,J)
	     SnoAmt2(I,J)   = SnoAmt1(I,J)
	     SnoAmt1(I,J)   = tmpdata
	     GIceComp(I,J)  = dmin1(1d0 - GIceComp(I,J),GIceCompMax)
	     tmpdata        = GT2(I,J)
	     GT2(I,J)       = GT1(I,J)
	     GT1(I,J)       = tmpdata
	    endif

          if (GHIce(I,J).lt.GHIce2(I,J)) then
c          оба сло€ тонкие
c          C*H1+(1-C)*H2 == Cnew*Hthin + (1-Cnew)*0   =>
	     GIceComp(I,J) = (GIceComp(I,J)*GHIce(I,J) +
	1         (1d0-GIceComp(I,J))*GHIce2(I,J)) / HIceMin
	      if (GIceComp(I,J).gt.GIceCompMax) GIceComp(I,J)=GIceCompMax 
	     GHIce(I,J)     = HIceThin
	     GHIce2(I,J)    = 0d0
	     GIce2Pool(I,J) = 0d0
           SnoAmt2(I,J)	= 0d0
	     GT2(I,J)       = TICE
	    endif

          if (GHIce2(I,J).lt.HIceMin) then
	     GHIce2(I,J)    = 0d0
	     GIce2Pool(I,J) = 0d0
           SnoAmt2(I,J)	= 0d0
	     GT2(I,J)       = TICE
	    endif

	    if (GHIce(I,J).lt.HIceMin) then
           if (FixedIceBorder) then
c           надо держать минимум льда
            GHICE(I,J)  = HIceMin
	     else
c           все ста€ло, кроме сев. полюса
	      if (J.ne.46) then
             GT(I,J)     = TICE    
             GHICE(I,J)  = 0.0
             GIceComp(I,J)  = 0.0
             ISFTYP(I,J) = 7
	      else
             GHICE(I,J)     = HIceMin
	       GIceComp(I,J)  = GIceCompMax
             GHIce2(I,J)    = 0d0
	       GIcePool(I,J)  = 0d0
	       GIce2Pool(I,J) = 0d0
             GT(I,J)        = TICE
             GT1(I,J)       = TICE
             GT2(I,J)       = TICE
             SnoAmt1(I,J)	 = 0d0
             SnoAmt2(I,J)	 = 0d0
	      endif
	     endif
	    endif
	   endif
	  endif
       enddo !i
      enddo !j

c „истка всего, что надо очистить дл€ следующего шага
      GIceCND     = 0.0
      GIce2CND    = 0.0
	GIceHeatng  = 0.0
	GIce2Heatng = 0.0
	NewIce      = .FALSE.
c	GIceUpFlow  = 7.5   !15.d0 !!!

c experiment!!!
c	TCICE = TICE !!!
      RETURN
      END
      
c********************************************* 
c     счет альбедо
c********************************************* 
c      FUNCTION Albedo(I,J)

C************ BEGINNING OF COMMON ************
c	include 'recom.fi'
C
C        NON STATE AND NON TEMPORARY VARIABLES IN MISC COMMON
C
c      COMMON / MISC / PHIS(74,46),ISFTYP(74,46),TS(74,46)
c     1  ,SD(74,46)
C
C        Q ARRAY - STATE VARIABLES
C
c      COMMON / QARY / P(74,46),U(74,46,2),V(74,46,2),T(74,46,2)
c     1  ,QW(74,46,2),GW(74,46),GT(74,46),SNOAMT(74,46)

c      include 'radvar.fi'
c      include 'comp.fi'
c      include 'ice.fi'
c      include 'paramz.fi'
C************ END OF COMMON ******************

c	if (ISFTYP(I,J).ne.9) then
c       Albedo = SFCALB(ISFTYP(I,J),1) + SQRT(DMIN1(SNR(I),1.0D0))
c     1    *(SFCALB(ISFTYP(I,J),2)-SFCALB(ISFTYP(I,J),1))
c	else
c	 if (GHIce(I,J).gt.0) then
c        CIceAlb(I)  = SFCALB(9,1) + SQRT(DMIN1(SNR1(I),1.0D0))
c     1    *(SFCALB(9,2)-SFCALB(9,1))
c	 else
c	  CIceAlb(I)  = SFCALB(7,1)
c	 endif
c	 if (GHIce2(I,J).gt.0) then
c        CIce2Alb(I) = SFCALB(9,1) + SQRT(DMIN1(SNR2(I),1.0D0))
c     1    *(SFCALB(9,2)-SFCALB(9,1))
c	 else
c	  CIce2Alb(I) = SFCALB(7,1)
c	 endif
c	 Albedo = CIceAlb(I)*GIceComp(I,J) + CIce2Alb(I)*(1-GIceComp(I,J))
c	 CIceAlb(I)   = (1d0-CIceAlb(I))/(1d0-Albedo)
c	 CIce2Alb(I)  = (1d0-CIce2Alb(I))/(1d0-Albedo)
c	endif

c      RETURN
c      END
