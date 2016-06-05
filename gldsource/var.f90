! var.fi common block containing variables
! EMBM variables added

      implicit none
      !#####
      common /data_path/ path_source, path_results
      character*40 path_source, path_results
      !#####
      integer maxi,maxj,maxk,maxl,mpxi,mpxj,isles,mpi
! parameter ( maxi = 36 , maxj =  36, maxk =  8 , maxl = 2)
      parameter ( maxi = 72 , maxj =  72, maxk =  8 , maxl = 2)
! for streamfunction equation
      parameter ( mpxi = maxi , mpxj = maxj+1 )
! for islands
      parameter ( isles = 1, mpi = 2*(maxi+maxj))
      integer imax,jmax,kmax,lmax,ntot,intot,k1(0:maxi+1,0:maxj+1) ,ku(2,maxi,maxj),mk(maxi+1,maxj)
      common /invars/imax,jmax,kmax,lmax,ntot,intot,k1,ku,mk
      integer ips(maxj),ipf(maxj),ias(maxj),iaf(maxj),jsf
      common /lego/ips,ipf,ias,iaf,jsf
      integer lpisl(mpi,isles), ipisl(mpi,isles), jpisl(mpi,isles) ,npi(isles)
! integer lpisl(mpi,isles+1), ipisl(mpi,isles+1), jpisl(mpi,isles+1)
! 1 ,npi(isles+1)
      common /islands/lpisl,ipisl,jpisl,npi

      real dt(maxk),dphi,dphi2,ds,ds2,dz(maxk),u(3,0:maxi,0:maxj,maxk),ts(maxl,0:maxi+1,0:maxj+1,0:maxk+1),&
      t,s(0:maxj),c(0:maxj),dzu(2,maxk),tau(2,maxi,maxj),drag(2,maxi+1,maxj),dztau(2,maxi,maxj),diff(2),ec(4),&
      cn,ratm(mpxi*mpxj,mpxi+1),ub(2,0:maxi+1,0:maxj),rho(0:maxi+1,0:maxj+1,0:maxk),&
      ts1(maxl,0:maxi+1,0:maxj+1,0:maxk+1),sv(0:maxj)
	 !ts1 - T and S at the previous time step
      real cv(0:maxj),dza(maxk),dztav(2,maxi,maxj),gb(mpxi*mpxj),gap(mpxi*mpxj,2*mpxi+3),ez0,cost(maxi,maxj),&
      rh(3,0:maxi+1,0:maxj+1),gbold(mpxi*mpxj) ,sda1,sdomg,dzz,tau0(maxi,maxj),dztav0(maxi,maxj) ,&
      tau1(maxi,maxj),dztav1(maxi,maxj),tsa0(maxj),t0
      real psi(0:maxi,0:maxj)
      common /vars/dt,dphi,dphi2,ds,ds2,dz,u ,ts,t,s,c,dzu ,tau,drag,dztau ,diff,ec ,cn,ratm,ub ,rho,ts1,sv,&
      cv ,dza,dztav,gb ,gap,ez0 ,cost,rh,gbold ,sda1,sdomg,dzz,tau0,dztav0 ,tau1,dztav1,tsa0,t0
      common /holes/psi
      real rel,u1(3,0:maxi,0:maxj,maxk)
      common /relax/rel,u1
! reciprocal variables to speed up fortran
      real rc(0:maxj),rcv(0:maxj),rdphi,rds,cv2(0:maxj),rc2(0:maxj) ,rtv(maxi,maxj),rtv3(maxi,maxj),rdz(maxk),&
      rdza(maxk)
      common /recips/rc,rcv,rdphi,rds,cv2,rc2,rtv,rtv3,rdz,rdza

! variables for non-trivial islands

      real bp(maxi+1,maxj,maxk), sbp(maxi+1,maxj)
      common /press/bp,sbp

! diagnostics
      real dmax
      common /testvar/dmax
      integer limps,istepT
      common /testint/limps,istepT

! dimensional scale values
      real usc,rsc,dsc,fsc,gsc,rh0sc,rhosc,cpsc,tsc,pi
      common /dimsc/usc,rsc,dsc,fsc,gsc,rh0sc,rhosc,cpsc,tsc,pi

! EMBM
      integer ndta
!######### seasonal solar
        INTEGER SDEDY,SDEYR,MNTHDY
	  real    tq_avr(2,maxi,maxj),ice_avr(2,maxi,maxj)
	  integer i_avr
!######### seasonal solar
      common /inebm/ndta,  SDEDY,SDEYR,MNTHDY, i_avr
      real cd,tq(2,maxi,maxj),tq1(2,maxi,maxj) ,qsata(maxi,maxj),qsato(maxi,maxj) ,&
      varice(2,maxi,maxj),varice1(2,maxi,maxj),tice(maxi,maxj) ,tqa(2,maxi,maxj),solfor(maxj),&
      ghs,scf ,dtatm,ryear,coef_co2,Time_co2,delta_co2 !############
      common /ebmvar/cd,tq,tq1 ,qsata,qsato ,varice,varice1,tice ,tqa,solfor,ghs,scf,dtatm,ryear,&
      tq_avr,ice_avr,coef_co2,Time_co2,delta_co2
!############
! constants for embm
      real emo,ema,solconst ,tfreez,alpha,alphice,rfluxsc,rfluxsca ,b00,b10,b20,b01,b11,b21,b02,&
      b12,b22,b03,b13,b23,delf2x,co20
      common /ebmconsts/emo,ema,solconst ,tfreez,alpha,alphice,rfluxsc,rfluxsca ,b00,b10,b20,b01,&
      b11,b21,b02,b12,b22,b03,b13,b23,delf2x,co20
! arrays for embm
      real albcl(maxi,maxj),albedo(maxi,maxj) ,fxsw(maxi,maxj),fxplw(maxi,maxj) ,fx0a(maxi,maxj),&
      fx0o(maxi,maxj) ,fxsen(maxi,maxj),pme(maxi,maxj),pmeadj(maxi,maxj) ,pptn(maxi,maxj),&
      evap(maxi,maxj),usurf(maxi,maxj) ,fxlata(maxi,maxj),fxlato(maxi,maxj) ,fxlw(maxi,maxj) ,&
      diffa(2,2,maxj),beta(2),hatmbl(2) ,ca(maxi,maxj),co2(maxi,maxj),runoff(0:maxi,0:maxj) ,&
      qb(maxi,maxj),fwflux_qb(maxi,maxj)
       real fx0sic(maxi,maxj),fx0t(maxi,maxj) ,fx0neto(maxi,maxj),fwfxneto(maxi,maxj) ,&
       evapsic(maxi,maxj),tsfreez(maxi,maxj) ,dtha(2,maxi,maxj)
      common /ebmflux/albcl,albedo,fxsw,fxplw ,fx0a,fx0o,fxsen,pme ,pmeadj,pptn,evap,usurf,&
      fxlata,fxlato ,fxlw,diffa,beta,hatmbl,ca,co2,runoff ,qb,fwflux_qb,fx0sic,fx0t,fx0neto ,&
      fwfxneto,evapsic,tsfreez ,dtha
! arrays for runoff scheme
      integer iroff(maxi,maxj),jroff(maxi,maxj)
      common /runoff/ iroff,jroff
! constants and parameters for atmosphere
      real rhoair,rhoao,cpa ,rho0,hlv,hls,hlf,const1,const2,const3,const4,const5,rmax ,saln0,&
      rpmesca,rpmesco ,diffmod0,ppmin,ppmax
      common /atmosconsts/rhoair,rhoao,cpa ,rho0,hlv,hls,hlf,const1,const2,const3 ,const4,const5,&
      rmax,saln0,rpmesca,rpmesco ,diffmod0,ppmin,ppmax
! prescribed/diagnosed atmospheric transports and velocities
      real transptq(2,maxi,maxj) ,uatm(3,maxi,maxj)
      common /transps/transptq,uatm
! constants and parameters for sea ice
      real rsictscsf,rhoice,rho0sea,difsic ,tsic,h0sic,hmin,rhoio,rhooi,rrholf,rhmin
      common /sicconsts/ rsictscsf,rhoice,rho0sea ,difsic,tsic,h0sic,hmin,rhoio,rhooi,rrholf,rhmin
! adjustable freshwater forcing parameters
      integer nsteps_extra0
      common /fwf_int/ nsteps_extra0
      real extra0,range0,rextra0,scl_fwf
      common /fwf_real/ extra0,range0,rextra0,scl_fwf
             !##### CO2 up to 2100     
      real B1_CO2(1:13),A2_CO2(1:13)
      integer Year2100(1:13),istepCO2
      common /CO2_Data/ B1_CO2,A2_CO2,Year2100,istepCO2    
