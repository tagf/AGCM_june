
! diag.f end-of-run diagnostics for gldstn

      subroutine diagend(lout,id_mnth)

      include 'var.f90'

      character name*10,lout*3,id_mnth*2 !=Jn or Jl ####
      real sum, tv2, syr
      parameter(syr = 365*86400)
                              !#####
      integer i, j, k, iwets, iice
	!######
	real plot(imax,jmax),tmin,tmax
      character*60 title,subtitle
	!######
   1  format (1x, i3,',',i3,',',e12.4)
! write out ocean depth

      do j=1,jmax
       do i=1,imax
               plot(i,j)=k1(i,j)
       enddo
      enddo

      name='Depth'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out SST

      do j=1,jmax
       do i=1,imax
               plot(i,j)=ts(1,i,j,kmax)
       enddo
      enddo

      name='SST'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out psi-barotropic streamfunction

      do j=1,jmax
       do i=1,imax
               plot(i,j)=1592.5*psi(i,j)
       enddo
      enddo

      name='BS'
      call DataForSurOc(lout,id_mnth,plot,name)


! write out SSS (g/kg)

      do j=1,jmax
       do i=1,imax
               plot(i,j)=saln0+ts(2,i,j,kmax)
               ! exclude Med
	         !if(i.ge.27.and.j.ge.28.and.j.le.30) plot(i,j)=34.9
       enddo
      enddo

      name='SSS'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out atmos air temp(deg C)

      do j=1,jmax
       do i=1,imax
         plot(i,j)=tq(1,i,j)
       enddo
      enddo
      name='tair'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out atmos specific humidity(g/kg)  - the ratio of the mass of water vapor in air to the total mass 
! of the mixture of air and water vapor. (grams of vapour per kilogram of air)


      do j=1,jmax
       do i=1,imax
         plot(i,j)=1.e3*tq(2,i,j)
       enddo
      enddo

      name='qair'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out sea ice thickness

      do j=1,jmax
       do i=1,imax
        if(varice1(1,i,j).eq.0.0) then
         plot(i,j)=0.0 !999.999
        else
         plot(i,j)=varice1(1,i,j)
        endif
       enddo
      enddo

      name='hice'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out sea ice areal fraction

      do j=1,jmax
       do i=1,imax
        if(varice1(2,i,j).eq.0.0) then
         plot(i,j)=0.0  !999.999
        else
         plot(i,j)=varice1(2,i,j)
        endif
       enddo
      enddo

      name='aice'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out sea ice temperature

      do j=1,jmax
       do i=1,imax
        if(varice1(2,i,j).eq.0.0) then
         plot(i,j)=-50.  !999.999
        else
         plot(i,j)=tice(i,j)
        endif
       enddo
      enddo

      name='tice'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out pptn - precipitation rate (cm/yr)

      do j=1,jmax
       do i=1,imax           !syr = 365*86400
         plot(i,j)=pptn(i,j)*syr*100.
       enddo
      enddo
      name='pptn'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out evap - evaporation rate (cm/yr)

      do j=1,jmax
       do i=1,imax
            tv2 = (evap(i,j)*(1-varice1(2,i,j)) + evapsic(i,j)*varice1(2,i,j))
         plot(i,j)=100.*tv2*syr
       enddo
      enddo

      name='evap'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out runoff - runoff rate (cm/yr)

      do j=1,jmax
       do i=1,imax
         plot(i,j)=100.*runoff(i,j)*syr
       enddo
      enddo

      name='runoff'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out net freshwater flux into ocean (P-E+R+freeze/melt)(cm/yr)

      do j=1,jmax
       do i=1,imax                    !syr = 365*86400
               plot(i,j)=fwfxneto(i,j)*syr*100.
       enddo
      enddo

      name='fwflux'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out net surface heat flux over ocean (W/m**2)

      do j=1,jmax
       do i=1,imax
               plot(i,j)=fx0neto(i,j)
       enddo
      enddo

      name='fluxo'
      call DataForSurOc(lout,id_mnth,plot,name)

! write out net surface heat flux into atmos (W/m**2)

      do j=1,jmax
       do i=1,imax
               plot(i,j)=fx0a(i,j)
       enddo
      enddo

      name='fluxa'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out divergence of atmos. heat transport W/m**2

        sum = 0.
        tmin=100.
        tmax=-100.
      do j=1,jmax
       do i=1,imax
               plot(i,j)=transptq(1,i,j)
       enddo
      enddo

      name='divt'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out divergence of atmos. moisture transport (*100.)

      do j=1,jmax
       do i=1,imax
               plot(i,j)=transptq(2,i,j)*100.   !*syr
       enddo
      enddo

      name='divq'
      call DataForSurAtm(lout,id_mnth,plot,name)

! write out 3-d field of temperature

      open(20,file=trim(path_results)//lout//'.temp'//id_mnth)
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               if(k.ge.k1(i,j))then
                  write(20,*)ts(1,i,j,k)
               else
                  write(20,*)0.0
               endif
            enddo
         enddo
      enddo
      close(20)

! write out 3-d field of salinity (g/kg)

      open(20,file=trim(path_results)//lout//'.saln'//id_mnth)
      do j=1,jmax
         do i=1,imax
            do k=1,kmax
               if(k.ge.k1(i,j))then
                  write(20,*)ts(2,i,j,k)
               else
                  write(20,*)0.0
               endif
            enddo
         enddo
      enddo
      close(20)

! write out albedo

      do j=1,jmax
       do i=1,imax
               plot(i,j)=albedo(i,j)
       enddo
      enddo

      name='albedo'
      call DataForSurAtm(lout,id_mnth,plot,name)

      write(6,*) 'final CO-2 concentration',co2(1,1)

      end
