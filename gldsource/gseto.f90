
! subroutine gseto, sets up geometry etc variable depth
! copied from v3_1 with coupled stuff added

      subroutine gseto

      include 'var.f90'

      real phix, th0, th1, z1, tv, tv1, tv2, tv3, tv4, tv5, zro(maxk), zw(0:maxk),&
       temp0, temp1, adrag, drgf, s0, s1
      real h(3,0:maxi+1,0:maxj+1)
      real tdata(maxi,maxj,maxk),sdata(maxi,maxj,maxk)
      real scl_kh, scl_kz, scl_tau, scl_adrag,arg,i0

      integer i,j,k,l,kmxdrg,jeb
      logical getj(maxi,maxj)

      common /lars/getj

      pi=4.*atan(1.0)

! dimensional scale (sc) values

      usc = 0.05   !vel (m/s)
      rsc = 6.37e6 !Earth rad (m)
      dsc = 5.e3   !depth (m)
      fsc = 2.*7.2921e-5   !Earth rotation rate (1/s)
	                    ! 7.2921e-5=2.*pi/(24*60*60)=omega
      gsc = 9.81  !g (m/s**2)
      rh0sc = 1.e3 !dens (kg/m3)
      rhosc = rh0sc*fsc*usc*rsc/gsc/dsc !=0.947  dens var
      cpsc = 3981.1  !heat cap?
      tsc = rsc/usc  !=1.27e8 sec = 1470 days = 4 year

      write(6,'(i4,3e12.4)')(k,dsc*zw(k),dsc*zro(k),dsc*dz(k) ,k=kmax,1,-1)  !k=0=>no print
! EMBM scaling for heat forcing of ocean

      rfluxsc = rsc/(dsc*usc*rh0sc*cpsc)

! EMBM reference salinity

      saln0 = 34.9

! EMBM scaling for P-E forcing of ocean  (Prec-Evap)

      rpmesco = rsc*saln0/(dsc*usc)

! parameters for setting up grid
! th is latitude, coords are sin(th), longitude phi, and z

      th0 = - pi/2.
      th1 = pi/2.
! th0 = - pi*7./18.
! th1 = pi*7./18.
      s0 = sin(th0)    !=-1
      s1 = sin(th1)    !=1
! phix = pi/3.
      phix = 2.*pi

! grid dimensions must be no greater than array dimensions in var.fi

      imax = 72  !lon
      jmax = 72  !lat
      kmax = 8   !vertic
      lmax = 2  !1 - temper, 2 - salin
! ocean positions (may be overwritten later in this routine)
! if (imax.eq.18) then
! ips = 2
! ipf = 7
! ias = 10
! iaf = 12
! jsf = 3
! else
! print*,'unknown ocean positions'
! stop
! endif

      dphi = phix/imax  !=2pi/imax lon step
      ds = (s1-s0)/jmax  !=2/jmax  sin lat step
      dphi2 = dphi*2.
      ds2 = ds*2.
      rdphi = 1.0/dphi
      rds = 1.0/ds

! set time to zero (may be overwritten if continuing run)

      t0 = 0.
      t = t0

! set timestep and initialize some 1-d arrays to zero

      print*,'ocean dt in days and A/O dt ratio'
      read(5,*)tv,ndta
      print*,tv,ndta
!####      tv=1.46  !(days)
!####      ndta=2
      tv = tv*86400./tsc  !dimless dt ocean (86400 s = 1day)
      dtatm = tv/ndta     !dimless dt atm
	 ! sdedy - day of year
	sdedy=0
	 ! SDEYR - year number
	SDEYR=0
! old nondim timestep
! print*,'dtatm, ndta'
! read(5,*)dtatm,ndta
! tv = ndta*dtatm

! variable timestep option not recommended
      do k=1,kmax
         dt(k) = tv
! initialize
         dzu(1,k) = 0.    !du/dz
         dzu(2,k) = 0.
      enddo

! set up grid
! For variable (exponential) dz use ez0 > 0, else use ez0 < 0

      ez0 = 0.1  !parameter controlling vertical grid spacing
! ez0 = - 1.0
      z1 = ez0*((1.0 + 1./ez0)**(1.0/kmax) - 1.0) !=3.5e-2
      print*,'z1',z1
      tv4 = ez0*((z1/ez0+1.)**0.5-1.)  !=1.6t-2
      tv2 = 0.
      tv1 = 0.
      zro(kmax) = -tv4
      zw(kmax) = tv2
      do k=1,kmax
         if(ez0.gt.0)then
            tv3 = ez0*((z1/ez0+1.)**k-1.)
            dz(kmax-k+1) = tv3 - tv2
            tv2 = tv3
            tv5 = ez0*((z1/ez0+1.)**(k+0.5)-1.)
            if(k.lt.kmax)dza(kmax-k) = tv5 - tv4
            tv4 = tv5
            tv1 = tv1 + dz(kmax-k+1)
! tv3 is the depth of the kth w level from the top
! tv5 is the depth of the k+1th density level from the top
         else
            dz(k) = 1.d0/kmax
            dza(k) = 1.d0/kmax
         endif
      enddo

      do k=kmax,1,-1
         if(k.gt.1)zro(k-1) = zro(k) - dza(k-1)
         zw(k-1) = zw(k) - dz(k)
      enddo
! write(6,'(i4,4e11.4)')(k,zw(k),zro(k),dz(k),dza(k),k=kmax,1,-1)
! write(6,'(i4,e11.4)')k,zw(0)

      dzz = dz(kmax)*dza(kmax-1)/2.

! efficiency array

      do k=1,kmax-1
         rdz(k) = 1.0/dz(k)
         rdza(k) = 1.0/dza(k)
      enddo
      rdz(kmax) = 1.0/dz(kmax)

! set up sin and cos factors at rho and v points (c grid)
! fix for global domain although only cv and cv2 are referred to at or beyond
! limits 24/6/2 if no flow out of N + S boundaries.

      do j=0,jmax
         sv(j) = s0 + j*ds  !s0 = sin(th0)=-1
         if(abs(1.0 - abs(sv(j))).lt.1.e-12)then
            cv(j) = 0.
            rcv(j) = 1.e12
         else
            cv(j) = sqrt(1. - sv(j)*sv(j)) !=cos
            rcv(j) = 1.0/cv(j)
         endif
         cv2(j) = cv(j)*cv(j)*rds
         s(j) = sv(j) - 0.5*ds   !sin of latitude rho points
         if(s(j).lt.-1.0) s(j) = -2.0 - s(j)
         c(j) = sqrt(1. - s(j)*s(j)) ! cos of latitude
         rc(j) = 1.0/c(j)
         rc2(j) = rc(j)*rc(j)*rdphi
! print*,'s TEST f-PLANE FUDGE'
! sv(j) = 0.5
! s(j) = 0.5
! print*,j,sv(j),s(j),c(j),cv(j)
      enddo

! set up coeffs for state equation following WS 1993

      ec(1) = - 0.0559 /1.18376
      ec(2) = 0.7968   /1.18376
      ec(3) = - 0.0063 /1.18376
      ec(4) = 3.7315e-5/1.18376

! read parameters

! print*,'temp0 temp1 rel scf'
! read(5,*)temp0,temp1,rel,scf
! print*,temp0,temp1,rel,scf

! read ocean initial temps
      print*,'temp0 temp1'
      read(5,*)temp0,temp1
!###      temp0=0.
!####	    temp1=0.
      print*,temp0,temp1

! read ocean parameters
      print*,'rel,scl_kh,scl_kz,scf'
      read(5,*)rel,scl_kh,scl_kz,scf
!###      rel=0.9
!###	scl_kh=1
!###	scl_kz=1
!###	scf=1
      print*,rel,scl_kh,scl_kz,scf

! define forcing

! read wind data

! taux,tauy at u-points  (zonal & meridional component of wind stress)
      open(96,file=trim(path_source)//'taux_u.interp')
      open(97,file=trim(path_source)//'tauy_u.interp')
! taux,tauy at v-points
      open(98,file=trim(path_source)//'taux_v.interp')
      open(99,file=trim(path_source)//'tauy_v.interp')

      do j=1,jmax
         do i=1,imax
! rotate grid to check b.c.s
! do idum=1,imax
! i=1+mod(idum+18-1,36)

            read(96,*)dztau(1,i,j)   !d(tau)/dz
            read(97,*)dztau(2,i,j)
            read(98,*)dztav(1,i,j)
            read(99,*)dztav(2,i,j)

! multiply by scaling factor scf (to drive correct gyre strengths)

            dztau(1,i,j) = scf*dztau(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztau(2,i,j) = scf*dztau(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztav(1,i,j) = scf*dztav(1,i,j)/(rh0sc*dsc*usc*fsc)/dzz
            dztav(2,i,j) = scf*dztav(2,i,j)/(rh0sc*dsc*usc*fsc)/dzz

            tau(1,i,j) = dztau(1,i,j)*dzz
            tau(2,i,j) = dztav(2,i,j)*dzz
         enddo
      enddo

      close(96)
      close(97)
      close(98)
      close(99)

! parameters for (restricted) time-dependent forcing
! set sda1 < 1e05 for steady forcing

! sda1 = 0.017
      sda1 = 0.000  !oscillatory forcing amplitude
      sdomg = 2.*pi/10.0

! seabed depth h needed BEFORE forcing if coastlines are non-trivial
! note k1(i,j) must be periodic ; k1(0,j) - k1(imax,j) = 0 and
! k1(1,j) - k1(imax+1,j) = 0

!number of wet points in the ocean and # interior wet points
      ntot = 0
      intot = 0

      open(13,file=trim(path_source)//'world.k1')
! note k1(i,j) must be periodic ; k1(0,j) - k1(imax,j) = 0 and
! k1(1,j) - k1(imax+1,j) = 0, as enforced below;

      do j=jmax+1,0,-1
         read(13,*)(k1(i,j),i=0,imax+1)
  !for debug       if (k1(i,j)==91.or.k1(i,j)==92.or.k1(i,j)==93)k1(i,j)=94
! rotate grid to check b.c.s
! read(13,*)xxx,(k1(i,j),i=19,36),(k1(i,j),i=1,18),xxx

         k1(0,j) = k1(imax,j)
         k1(imax+1,j) = k1(1,j)
         do i=0,imax+1
! boundary condition
! if(i*j*(imax+1-i)*(jmax+1-j).eq.0)then
! k1(i,j) = kmax+1
! else
! k1(i,j) = 1
! endif

            do k=1,3
              h(k,i,j) = 0.
              rh(k,i,j) = 0. !reciprocal of ocean depth at u,v,T
            enddo         !(обратная величина)            points resp.
         enddo
         write(6,'(i4,74i3)')j,(k1(i,j),i=0,imax+1)
! write(99,'(38i3)')(k1(i,j),i=0,imax+1)
      enddo

! read ips etc if possible

! read(13,*,end=200)ips,ipf,ias,iaf,jsf
!  arrays defining start and finish points of Pacific and
!  Atlantic and finish of southern ocean respectively
 200  close(13) !k1
! print*,'ocean positions ',ips,ipf,ias,iaf,jsf

! count wet cells

      do j=1,jmax
         do i=1,imax                !number of wet points in the
            if(k1(i,j).le.kmax)then !ocean and # interior wet points
               ntot = ntot + kmax - k1(i,j) + 1
               intot = intot + kmax - k1(i,j)
            endif
         enddo
      enddo

! find ocean positions semi-automatically, must start with a
! longitude i which is in the right ocean for all j, tricky in north

! stop 'testing new code'
!ias(j), iaf(j), ips(j, ipf(j)- defining start and
!finish of Pacific and Atlantic basins at each latitude j, and jsf for finish
!of southern ocean
      ias(jmax) = 49
      ips(jmax) = 16
      jsf = 1
      print*,'j   ips ipf ias iaf '
      do j=1,jmax-1
         ips(j) = ips(jmax)
         ipf(j) = ips(j)
         ias(j) = ias(jmax)
         iaf(j) = ias(j)
! if(j.eq.jmax-1)then
! ias(j) = 22
! endif
         do i=1,imax
            if(k1(ips(j)-1,j).le.kmax)ips(j) = ips(j) - 1
            if(k1(ipf(j)+1,j).le.kmax)ipf(j) = ipf(j) + 1
            if(k1(ias(j)-1,j).le.kmax)ias(j) = ias(j) - 1
            if(k1(iaf(j)+1,j).le.kmax)iaf(j) = iaf(j) + 1
            ips(j) = 1 + mod(ips(j)-1+imax,imax)
            ipf(j) = 1 + mod(ipf(j)-1+imax,imax)
            ias(j) = 1 + mod(ias(j)-1+imax,imax)
            iaf(j) = 1 + mod(iaf(j)-1+imax,imax)
         enddo  !find j border of south ocean
         if(ias(j).ge.iaf(j).and.j.le.jmax/2)jsf = j
         if(ips(j).ge.ipf(j).and.j.le.jmax/2)jsf = j
      enddo
      ips(jmax) = 1
      ipf(jmax) = 0
      ips(jmax-1) = 1
      ipf(jmax-1) = 0
      ias(jmax) = 1
      iaf(jmax) = imax
      write(6,'(5i4)')(j,ips(j),ipf(j),ias(j),iaf(j),j=1,jmax)
      print*,'jsf ',jsf
      !####
	open (21,file=trim(path_results)//'oceans')
      write(21,*)'j   ips ipf ias iaf '
	write(21,'(5i4)')(j,ips(j),ipf(j),ias(j),iaf(j),j=1,jmax)
	close (21)
      !####
! initialize psi  - barotropic stream function

      do j=0,jmax
         do i=0,imax
           psi(i,j)=0.0
         enddo
         do i=0,imax+1
            ub(1,i,j) = 0.  !barotropic velocity components
            ub(2,i,j) = 0.
         enddo
      enddo

! seabed depth h

      do j=jmax+1,0,-1
         do i=0,imax+1
            if(k1(i,j).le.kmax)then  !kmax - max number of
               do k=k1(i,j),kmax                  !vertical levels
                  h(3,i,j) = h(3,i,j) + dz(k)
               enddo
               rh(3,i,j) = 1.0/h(3,i,j) !reciprocal of ocean depth at
			                          ! u,v,T points resp.(=1,2,3)
            endif
         enddo
! write(6,'(i4,40f5.0)')j,(dsc*h(3,i,j),i=0,imax+1)
      enddo

      do j=0,jmax+1
         do i=0,imax
            h(1,i,j) = min(h(3,i,j),h(3,i+1,j))  ! u points
            if(max(k1(i,j),k1(i+1,j)).le.kmax)rh(1,i,j) = 1.0/h(1,i,j)
         enddo
      enddo

      do j=0,jmax
         do i=0,imax+1
            h(2,i,j) = min(h(3,i,j),h(3,i,j+1))  ! v points
            if(max(k1(i,j),k1(i,j+1)).le.kmax)rh(2,i,j) = 1.0/h(2,i,j)
         enddo
      enddo

      do 120 j=1,jmax
         do 120 i=1,imax !k1-array defining the ocean domain ie geometry
	                    ! and topography wet at i,j for k=k1(i,j),...,kmax
			    !ku -  similar array derived from k1
	                    ! defining wet velocity points
            ku(1,i,j) = max(k1(i,j),k1(i+1,j))
            ku(2,i,j) = max(k1(i,j),k1(i,j+1))
  120 continue
      tv2 = 0.

! read interpolated Levitus data

! open(30,file=trim(path_source)//'tempann.silo')
! read(30,*)(((tdata(i,j,k),k=1,kmax),i=1,imax),j=1,jmax)
! close(30)
! open(30,file=trim(path_source)//'saliann.silo')
! read(30,*)(((sdata(i,j,k),k=1,kmax),i=1,imax),j=1,jmax)
! close(30)
! open(30,file=trim(path_source)//'tmp.levi')
! do j=1,jmax
! do i=1,imax
! do k=1,kmax
! if(k.ge.k1(i,j))then
! write(30,'(f11.4)')tdata(i,j,k),sdata(i,j,k)
! else
! write(30,'(f11.4)')0.0,0.0
! endif
! enddo
! enddo
! enddo
! close(30)

! set up drag and diffusion values
! drag takes the value adrag in the interior, rising twice by factor
! drgf per gridpoint close to equator and in regions of
! shallow water (k1>kmxdrg) ie land in the case kmxdrg=kmax
! jeb = 1/2 width of equatorial region of maximum drag

      adrag = 1.0/2.5/86400./fsc ! basic drag coefficient;
                                ! fsc-Earth rotation rate
! cross equator need * 4 if drag is constant ie if drgf=1
! for variable drag
! adrag = adrag * 1.0
! drgf = 3.0
! kmxdrg = kmax/2
! jeb = 1
! or
! adrag = adrag * 2.0
! drgf = 2.0
! kmxdrg = kmax
! jeb = 0
! for constant drag (using Bob's parm set)
! adrag = adrag * 8.0
! drgf = 1.0
! kmxdrg = kmax
! jeb = 0

! read drag parameters
      print*,'scl_adrag,drgf,kmxdrg,jeb'
      read(5,*)scl_adrag,drgf,kmxdrg,jeb
!###  scl_adrag=1.
!###	drgf=3.
!###	kmxdrg=4
!###	jeb=1
      print*,scl_adrag,drgf,kmxdrg,jeb

      adrag = adrag * scl_adrag
! drag(2,i,j), with basic value adrag (local variable) set in routine
!   drgset.f -  subroutine to define drag matrix drag(2,i,j)
      call drgset(adrag,drgf,kmxdrg,jeb)
      print*,'drag ds',adrag,ds,adrag/ds

! diff(1) = 3000./rsc/usc
! diff(1) = diff(1) * 4
! diff(2) = 1e-4/usc/dsc/dsc*rsc
! diff(2) = diff(2) * 2.0

! Stefan Rahmstorf & Andrew Weaver use:
      diff(1) = 2000./rsc/usc   !diffusion in phi (1) and s (2) dir'ns
      diff(2) = 1.e-4/usc/dsc/dsc*rsc

! scale ocean horiz diffusivity by factor scl_kh,
! vertical diffusivity by factor scl_kz

      diff(1) = scl_kh*diff(1)
      diff(2) = scl_kz*diff(2)

! arrays for efficiency

      do j=1,jmax !rtv and rtv3 in velc become arrays
         do i=1,imax
            rtv(i,j) = 1.0/(s(j)*s(j) + drag(1,i,j)*drag(1,i,j))
            rtv3(i,j) = 1.0/(sv(j)*sv(j) + drag(2,i,j)*drag(2,i,j))
         enddo
      enddo

      print*,'dphi ds diff(1) diff(2)'
      print*,dphi,ds,diff(1),diff(2)

! initialize some arrays to zero

      do i=0,imax
         do j=0,jmax
            do k=1,kmax
               do l=1,3
                  u(l,i,j,k) = 0.  !ocean velocity array
                  u1(l,i,j,k) = 0. !velocity at previous timestep
               enddo
            enddo
         enddo
      enddo

! initial conditions

      do i=0,imax+1
         do j=0,jmax+1
            do k=0,kmax+1
! initial uniform temperature T0 large favours thermally direct solutions
               if(j.le.jmax/2)then
                  ts(1,i,j,k) = temp0 *0.5*(1 + sign(1,k-k1(i,j)))
!ocean temp(1) and salinity(2)arr
               else
                  ts(1,i,j,k) = temp1 *0.5*(1 + sign(1,k-k1(i,j)))
               endif
! if(i.eq.1.and.k.eq.1)print*,ts(1,i,j,k)
! initial salinity     !ocean temp(1) and salinity(2)arr
               ts(2,i,j,k) =  0.0
! ts(2,i,j,k) =  saln0
               ts1(1,i,j,k) = ts(1,i,j,k)
               ts1(2,i,j,k) = ts(2,i,j,k)
            enddo
            do k=0,kmax   !density
               rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k) + ec(3)*ts(1,i,j,k)**2 + ec(4)*ts(1,i,j,k)**3
            enddo
         enddo
      enddo

! forcing fields and some more initialisation

      do 20 i=1,imax
         do 30 j=1,jmax
! th at u and v points
            tv = asin(s(j))  ! latitude
            tv1 = asin(sv(j))
! convective frequency array
!array of convection depth (in gridpoints) averaged over a run
            cost(i,j) = 0.

            rho(i,j,0) = 0. !density
! rho(i,j,0) is never referenced but this is not easy to prove in co.f

! wind stress, tau(1,i,j) is the u component at a u point
! tau(2,i,j) is the v component at a v point
! BUT NB this is confusing notation...
! dztau(l,i,j) are the deriv's of the u,v components at u points
! dztav(l,i,j) are the deriv's of the u,v components at v points

! tau(1,i,j) = - ta0*cos(2*pi*tv/th1)
! tau(1,i,j) = 0.
! if(j.eq.1)tau(1,i,j) = tau(1,i,j) + ta0
! dztau(1,i,j) = tau(1,i,j)/dzz
! this is also needed at v points
! tau(2,i,j) = 0.
! dztau(2,i,j) = 0.
! dztav(1,i,j) = - ta0*cos(2*pi*tv1/th1)
! dztav(1,i,j) = 0.
! 1                     /dzz
! dztav(2,i,j) = 0.
   30    continue !j
   20 continue  !i

! array to determine limit of easy part of double p integral in J term
! use integer wet point indicator
! (1+sign(1,kmax-k1(i,j)))/2
! mk - array used to reduce unnecessary calculation in streamfunction eq'n
! mk is largest of surrounding wet k1 values if i,j is wet, else 0
      do 130 j=1,jmax
         do 130 i=1,imax
            mk(i,j) = max(k1(i,j)*(1+sign(1,kmax-k1(i,j)))/2,&
                      k1(i+1,j)*(1+sign(1,kmax-k1(i+1,j)))/2,&
                      k1(i-1,j)*(1+sign(1,kmax-k1(i-1,j)))/2,&
                      k1(i,j+1)*(1+sign(1,kmax-k1(i,j+1)))/2,&
                      k1(i,j-1)*(1+sign(1,kmax-k1(i,j-1)))/2)
            mk(i,j) = mk(i,j) * (1+sign(1,kmax-k1(i,j)))/2
  130 continue

! periodic b.c. required for Antarctic island integral

      do j=1,jmax
         mk(imax+1,j) = mk(1,j)
      enddo


! array to avoid J term in flat regions
! For non-trivial coasts essential to avoid adding J term at non-Psi points.
! Hence first condition ensures (i,j) is a wet Psi point, 2nd that bottom
! is not flat.

      do 140 j=1,jmax
         do 140 i=1,imax
            if( (max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.kmax)&
             .and. (k1(i,j).ne.k1(i,j+1).or.k1(i,j).ne.k1(i+1,j) .or.k1(i,j).ne.k1(i+1,j+1)))then
               getj(i,j) = .true.
            else
               getj(i,j) = .false.
            endif
  140 continue

! read island geometry file or write out for manual editing
! setting wet to zero and 1 on 1st landmass, 2 on 2nd landmass (1st island)
! etc nb narrow channels may have no wet psi points and hence not show up on
! psi grid

      open(23,file=trim(path_source)//'world.psiles',status='old')

! open(23,file=trim(path_source)//'world.psiles',status='new')
! use enquire stmt if you can remember how
! do j=jmax,0,-1
! write(23,'(180i3)')((1-sign(1,kmax-
! 1      max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1))))/2,i=1,imax)
! write(6 ,'(180i3)')((1-sign(1,kmax-
! 1      max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1))))/2,i=1,imax)
! enddo
! stop 'psi file written'
      do j=jmax,0,-1 !gbold(maxi*maxj)- initial value of gb
         read(23,*)(gbold(i + j*imax),i=1,imax)
! rotate grid to check b.c.s
! read(23,*)(gbold(i+j*imax),i=19,36),(gbold(i+j*imax),i=1,18)
      enddo
      close(23)

! read island path integral data, read isles+1 paths only if want last path
! for testing

      open(24,file=trim(path_source)//'world.paths')
! read(24,*)(npi(i),i=1,isles+1)
! do i=1,isles+1   !npi(i)- array specifying length of each path
      read(24,*)(npi(i),i=1,isles)  !isles - number of islands (0 or 1)
      do i=1,isles  !(arbitrary) maximum for length of path integrals
         read(24,*) ! around islands     mpi = 2*(maxi+maxj)
         if(npi(i).gt.mpi)stop 'path integral around island too long'
         do j=1,npi(i)
            read(24,*)lpisl(j,i), ipisl(j,i), jpisl(j,i)
   !!         jpisl(j,i)=5
! rotate grid to check b.c.s
! ipisl(j,i) = 1 + mod(ipisl(j,i)-1+18,36)

            if(abs(lpisl(j,i)).ne.1.and.abs(lpisl(j,i)).ne.2)stop  !it means =1 or 2
            if(ipisl(j,i).gt.imax.or.ipisl(j,i).lt.0)stop 'bad path'
            if(jpisl(j,i).gt.jmax.or.jpisl(j,i).lt.0)stop 'bad path'
            if(k1(ipisl(j,i),jpisl(j,i)).gt.kmax)stop 'dry path' !it means must be ocean points
         enddo
      enddo
      close(24)

      print*,'horizontal diffusivity',diff(1)*rsc*usc,' m**2/s'
      print*,'vertical diffusivity',diff(2)*usc*dsc*dsc/rsc,' m**2/s'
      print*,'basic drag coefficient',adrag*fsc,' /s'
      print*,'wind stress scale',fsc*usc*dsc,' m**2/s**2'
      print*,'or',fsc*usc*dsc*rh0sc,' N/m**2'
      print*,'density variation scale',rhosc,' kg/m**3'
      print*,'vertical velocity scale',usc*dsc/rsc,' m/s'
      print*,'time scale',rsc/usc/86400./365.,' yrs'
      print*,'overturning scale',dsc*usc*rsc*1.e-6,' Sv'
      print*,'vertical heat flux scale',dsc*usc*rh0sc*cpsc/rsc,' W/m**2'
      print*,'integrated energy scale',rh0sc*fsc*usc*rsc**3*dsc,' J'
      print*,'integrated northward heat flux scale (W)'
      write(6,'(e15.5)') usc*rh0sc*cpsc*rsc*dsc

      end
