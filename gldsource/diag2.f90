
! diag2.f frequent diagnostics for program gldstn variable depth
! mo - indexes the oceans 1=Pacific, 2=Atlantic, 3=Indian, 4=Southern
! ndep - indexes the depth 0=deep, 1=upper
! ipf=end of Pacific, iaf=end of Atlantic, jsf=north end of Southern
! lmax.gt.2.allowed

      subroutine diag2(sum,avn,avs)

      include 'var.f90'

      real sum(8*maxl), avn,  avs, vol(8), tv

      integer i,j,k,l,mo,ndep

      logical first

      save vol, first

      data first/.true./
      data vol/8*0.0/

      do i=1,8*lmax  !lmax = 2  1 - temper, 2 - salin
         sum(i) = 0.
      enddo
      tv = 0.
      avn = 0.
      avs = 0.
      do 10 k=1,kmax
         do 20 j=1,jmax
            do 30 i=1,imax
               if(k.ge.k1(i,j))then
                  ndep = (2*k-1)/kmax !0=deep, 1=upper
                  if(j.le.jsf)then
                     mo = 4
                  elseif(i.ge.ips(j).and.i.le.ipf(j).and.j.ne.jmax)then
                     mo = 1
                  elseif(i.ge.ias(j).and.i.le.iaf(j).and.j.ne.jmax)then
                     mo = 2
                  else
                     mo = 3
                  endif
                  do 60 l=1,lmax  !lmax = 2
                     sum(mo + 4*ndep + 8*(l-1)) = sum(mo + 4*ndep + 8*(l-1)) + ts(l,i,j,k)*dphi*ds*dz(k)
! print*,i,j,k,l,mo + 4*ndep + 8*(l-1)
! 2                   ,sum(mo + 4*ndep + 8*(l-1))
   60             continue
                  if(first)vol(mo + 4*ndep) = vol(mo + 4*ndep) + dphi*ds*dz(k)
! print*,i,j,k,mo + 4*ndep,vol(mo + 4*ndep)
                  if(k.lt.kmax) avn = avn + abs(rho(i,j,k) - rho(i,j,k+1))/dza(k)

                  avs = avs + u(1,i,j,k)**2 + u(2,i,j,k)**2

               endif
   30       continue
   20    continue
   10 continue

! crash barrier added 111099
      if(avs.lt.1.e20)then
        continue
      else
        print*,'big avs , stopping'
        stop
      endif

      do i=1,8
         do l=1,lmax
            sum(i + 8*(l-1)) = sum(i + 8*(l-1))/vol(i)
         enddo
      enddo

      avs = sqrt(avs/ntot)

      avn = avn/(intot)

      first = .false.

      end
