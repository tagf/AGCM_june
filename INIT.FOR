      subroutine init
C
      REAL KAPA,KAPEL,LAT
      INTEGER  SDEDY
C
      COMMON /RECOM/
     1 TAU   ,TAUI  ,TAUO  ,TAUD  ,TAUE  ,TAUH  ,TAUC  ,TOFDAY,ROT   ,
     2 DT    ,DLAT  ,DLON  ,RAD   ,RSDIST,SIND  ,COSD  ,COSR  ,SINR  ,
     3 DAYPYR,ROTPER,SDEYR ,SOLTCE,APHEL ,DECMAX,ECCN  ,GUSTY ,
     4 DAY   ,GRAV  ,RGAS  ,KAPA  ,PSF   ,PTROP ,PSL   ,CFG   ,
     5 FM    ,ED    ,PI    ,SIG1  ,SIG3  ,DSIG  ,PM    ,KAPEL ,RKAPA1,
     6 STBO  ,GWM   ,DTC3  ,FRDAY ,CLH   ,COE1  ,HICE  ,CTYICE,CNDICE,
     7 TICE  ,TCICE ,SNOWL ,COE   ,TSPD  ,PSTQ  ,QST   ,TST   ,PTRK  ,
     8 PSFHO ,CALFA ,QC    ,PC    ,QCONST,EFVCON,TCT0  ,TCT2  ,TCT4  ,
     9 TCST0 ,TCST2 ,TC01  ,TC02  ,TC03  ,TC04  ,TC23  ,TC24  ,TC34  ,
     1 BLC   ,ALOGP0,FIM   ,HRGAS ,TCST4 ,FLR   ,PS4K  ,PS8K  ,ELOG  ,
     2 S0,TOZONE    ,LAT(46)   ,DXU(46)   ,DXP(46)   ,DYU(46),DYP(46),
     3 SFCALB(9,2),SINL(46)  ,COSL(46)  ,O3AMT(46) ,TC12  ,TC14  ,
     4 DXYP(46)  ,F(46)     ,SIG(2)    ,FLEADN    ,FLEADS    ,ERROR  ,
     5 COSLN(72),SINLN(72),DXYU1(46),DXYU2(46), GMT,GMR,GMKE,
     6 JM,IM,ID,MNTHDY,SDEDY ,NCYCLE,NC3   ,MONTH ,MRCH  ,NSTEP,nav
C
c     jm  number of north-south grid points
c     im  number of east-west grid points
c
      data jm,im /46,72/
      data nav /0/
c
c     
      data    taui /8016.0/
      data    tauo,taud /1.0,24./
      data    tauh,tauc,id,dt /6.0,24.,4hagcm,600.0/
c
      data    rad /6375000.e0/
c
      data    daypyr,rotper,soltce,aphel /365.,24.,173.,182./
c
      data  eccn,day,grav,rgas,kapa /.0178e0,86400.,9.81e0,287.,.286e0/
c
      data    psf,ptrop,psl /984.,200.,1000./
c
      data    ncycle,nc3,fm,ed /6,6,3.62425e-7,10.0/
c
      data    pi /3.1415926536e0/
      data    nstep /0/
      data    lat,dxu,dxp,dyu,dyp/46*0.0,46*0.0,46*0.0,46*0.0,46*0.0/
      data    sinl,cosl,o3amt,dxyp,f/46*0.0,46*0.0,46*0.0,46*0.0,46*0.0/
      data    sig /.25e0,.75e0/
      data    sfcalb /
     1 0.12e0,0.10e0,0.13e0,0.20e0,0.25e0,0.19e0,0.07e0,0.45e0,0.45e0,
     2 0.50e0,0.41e0,0.58e0,0.53e0,0.55e0,0.70e0,0.07e0,0.80e0,0.80e0/
c     
      trco(x,y)=amin1(1.e0,1.187e0-.066e0*alog10(abs(x-y)*(x+y)))
      decmax=23.5e0*pi/180.0
      dlat=4.0*pi/180.0
      dlon=5.0*pi/180.0
      jmm2=jm-2
      imm1=im-1
      jmm1=jm-1
      imm2=im-2
      fim=im
      hrgas=rgas/2.
      kapa=0.286e0
      cfg=0.006e0/(5000.*grav)
      rkapa1=1.0/(2.0*kapa)
      sig1=sig(1)
      sig3=sig(2)
      dsig=0.5
      stbo=1.171e-7
      gwm=15.e0
      dtc3=float(nc3)*dt
      frday=2.e0*pi/day
      clh=600.e0/.24e0
      hice=300.e0
      ctyice=5.4e-3
      cndice=day*ctyice/hice
      tice=273.1e0
      tcice=tice-1.6e0
      snowl=80.e0
      gusty=2.e0
c
      pm=psl-ptrop
      coe=grav*100./(0.5e0*pm*1000.*0.24e0)
      coe1=coe*dtc3/(24.0*3600.)
      tspd=day/dtc3
      pstq=100.
      qst=2.5e-6
      tst=218.0
      psfho=0.9e0
      calfa=1.+psfho
      qc=alog10(qst)
      pc=alog10(pstq)
      qconst=1020.4e0*qst
      efvcon=(qconst/calfa)*(pstq/1000.)**calfa
      flr=0.5e0*0.1828e0/(30.48e0*grav)
      elog=alog(10.e0)
      kapel=kapa*elog
      ps4k=400.0e0**kapa
      ps8k=800.0e0**kapa
      blc=274.3e0
      alogp0=alog10(ptrop)
      ptrk=exp(kapa*alogp0*elog)
c
       tct0=trco(0.0e0,200.e0)
       tct2=trco(0.0e0,600.e0)
       tct4=trco(0.0e0,1000.e0)
       tcst0=trco(100.0e0,200.e0)
       tcst2=trco(100.0e0,600.e0)
       tcst4=trco(100.0e0,1000.e0)
       tc01=trco(200.0e0,400.e0)
       tc02=trco(200.0e0,600.e0)
       tc03=trco(200.0e0,800.e0)
       tc04=trco(200.0e0,1000.e0)
       tc12=trco(400.0e0,600.e0)
       tc14=trco(400.0e0,1000.e0)
       tc23=trco(600.0e0,800.e0)
       tc24=trco(600.0e0,1000.e0)
       tc34=trco(800.0e0,1000.e0)
c
      sdeyr=0.0
      sdedy=365
      tau=0.0
c
      tofday=amod(tau+23.0,rotper)+1.0
      rot=tofday/rotper*(pi*pi)
      cosr=cos(rot)
      sinr=sin(rot)
      nstep=(tau*3600.0)/dt+.01e0
c
      fjm=jm
      fje=fjm/2.0+0.5e0
      do 410 j=2,jmm1
      fj=j
 410  lat(j)=dlat*(fj-fje)
      lat(1)=-pi/2.0 
      lat(jm)=pi/2.0                        
c
      do 415 j=2,jm
 415  dyu(j)=rad*(lat(j)-lat(j-1))
      dyu(1)=dyu(2)
c
      do 420 j=1,jm                          
 420  dxp(j)=rad*cos(lat(j))*dlon
      dxp(1)=0.0
      dxp(jm)=0.0
c     
      do 430 j=2,jm                                                                                     
 430  dxu(j)=0.5e0*(dxp(j)+dxp(j-1))
      dxu(1)=dxu(2)
c
      do 440 j=2,jmm1
 440  dyp(j)=0.5e0*(dyu(j+1)+dyu(j))
      dyp(1)=dyu(2)
      dyp(jm)=dyu(jm)
c
      do 445 j=2,jmm1
 445  dxyp(j)=0.5e0*(dxu(j)+dxu(j+1))*dyp(j)
      dxyp(1)=dxu(2)*dyp(1)*0.25e0
      dxyp(jm)=dxu(jm)*dyp(jm)*0.25e0
c
      do 450 j=2,jmm1
 450  f(j)=2.0*pi/day*(rad/dxyp(j))*((cos(lat(j-1))+cos(lat(j)))*
     * dxu(j)-(cos(lat(j))+cos(lat(j+1)))*dxu(j+1))/2.0
      f(jm)=2.0*pi/day*(rad/dxyp(jm))*(cos(lat(jm-1))+cos(lat(jm)))
     * *dxu(jm)/2.0
      f(1)=-f(jm)
c    
      do 510 j=1,jm
      sinl(j)=sin(lat(j))
      cosl(j)=cos(lat(j))
 510  continue
      cosl(1)=0.0
      cosl(jm)=0.0
      sinl(1)=-1.0
      sinl(jm)=1.0
      do 600 i=1,72
      alon=-pi+(i-1)*dlon
      cosln(i)=cos(alon)
 600  sinln(i)=sin(alon)
c
      do 5 k=3,45
      dxyu2(k)=0.5e0*dxyp(k)/(dxyp(k)+dxyp(k-1))
      dxyu1(k)=0.5e0*dxyp(k-1)/(dxyp(k)+dxyp(k-1))
 5    continue
      dxyu2(2)=0.5e0*dxyp(2)/(dxyp(2)+dxyp(1)+dxyp(1)) 
      dxyu2(46)=dxyp(46)/(dxyp(46)+dxyp(46)+dxyp(45))
      dxyu1(2)=dxyp(1)/(dxyp(2)+dxyp(1)+dxyp(1))
      dxyu1(46)=0.5e0*dxyp(45)/(dxyp(46)+dxyp(46)+dxyp(45))
c
      return
      end
                                                                                                                                            