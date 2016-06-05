
! diag.f diagnostics for program gldstn variable depth
! lmax > 2 allowed

      subroutine diag

      include 'var.f90'

      real sum(maxl), sumsq(maxl), umax(3), cnmax, avp, avs, vol, tmax ,tmin, ubmax(2), enmax(maxk), rdmn, tv, tmax1, tmin1

      integer i, j, k, l, imm(3,3), icp, icp2, isl, icp10, ird, jrd , ien(maxk), jen(maxk)

      vol = 0.0
      do l=1,lmax  !lmax=2
         sum(l) = 0.0
         sumsq(l) = 0.0
      enddo
      do i=1,3
         umax(i) = 0
         do j=1,3
            imm(i,j) = 0
         enddo
      enddo
      cnmax = 0
      tv = 0
      icp = 0
      icp2 = 0
      icp10 = 0
      avp = 0
      avs = 0
      cnmax = 0
      ubmax(1) = 0
      ubmax(2) = 0
      do 10 k=1,kmax
! cnmax = 0
! cnvax = 0
         do 20 j=1,jmax
            do 30 i=1,imax
               if(k.ge.k1(i,j))then
                  vol = vol + dphi*ds*dz(k)
                  do l=1,lmax
                     sum(l) = sum(l) + ts(l,i,j,k)*dphi*ds*dz(k)
                     sumsq(l) = sumsq(l) + ts(l,i,j,k)**2*dphi*ds*dz(k)
                  enddo
                  do 40 l=1,3
! if(l.eq.1.and.ts(l,i,j,k).eq.0)print*,i,j,k,'zero'
                     if(abs(u(l,i,j,k)).gt.umax(l))then
                        umax(l)=abs(u(l,i,j,k))
                        imm(1,l) = i
                        imm(2,l) = j
                        imm(3,l) = k
                     endif
   40             continue
! tv = max(abs(u(1,i,j,k))/dphi,abs(u(2,i,j,k))/ds
! 1                 ,abs(u(3,i,j,k))/dz(k))*dt(k)
                  tv = max(abs(u(1,i,j,k))*rc(j)*rdphi ,abs(u(2,i,j,k))*cv(j)*rds ,abs(u(3,i,j,k))*rdz(k))*dt(k)
                  if(tv.gt.cnmax)cnmax = tv
! tv = abs(u(3,i,j,k))/dz(k)*dt(k)
! if(tv.gt.cnvax)cnvax = tv

! if(ts(1,i,j,k).lt.-0.001)then
! print*,'T<0 at',i,j,k
! stop
! endif

                  avs = avs + u(1,i,j,k)**2 + u(2,i,j,k)**2

! tv = abs(u(1,i,j,k)*dza(k)*dza(k)/c(j)/dphi/diff(2))
! tv = abs(u(2,i,j,k)*dza(k)*dza(k)*c(j)/ds/diff(2))
                  tv = abs(u(3,i,j,k)*dza(k)/diff(2))
                  avp = avp + tv
                  if(tv.gt.2)then
                     icp = icp + 1
                     if(tv.gt.4)icp2 = icp2 + 1
                     if(tv.gt.20)icp10 = icp10 + 1
                  endif
               endif
   30       continue
   20    continue
! if(k.eq.1.or.k.eq.kmax)
! 1   print*,'k dt cn ',k,dt(k),cnmax
   10 continue
! write(6,'(5e10.3)')(((ts(2,i,j,k),i=2,10,2),j=2,10,2),k=2,10,2)

      tmax1 = -1e10
      tmin1 = 1e10
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).eq.1)then
               if(ts(1,i,j,1).gt.tmax1)tmax1=ts(1,i,j,1)
               if(ts(1,i,j,1).lt.tmin1)tmin1=ts(1,i,j,1)
            endif
         enddo
      enddo

      tmax = -1e10
      tmin = 1e10
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).le.kmax)then
               if(ts(1,i,j,kmax).gt.tmax)tmax=ts(1,i,j,kmax)
               if(ts(1,i,j,kmax).lt.tmin)tmin=ts(1,i,j,kmax)
            endif
         enddo
      enddo

! find max barotropic velocity components

      do j=1,jmax
         do i=1,imax
            do l=1,2
               if(abs(ub(l,i,j)).gt.ubmax(l))ubmax(l) = abs(ub(l,i,j))
            enddo
         enddo
      enddo

! find density differences

! rdmn = rho(1,1,k1(1,1))-rho(1,1,kmax)
! do j=1,jmax
! do i=1,imax
! if(k1(i,j).le.kmax)then
! if(rho(i,j,k1(i,j))-rho(i,j,kmax).lt.rdmn)then
! rdmn = rho(i,j,k1(i,j))-rho(i,j,kmax)
! ird = i
! jrd = j
! endif
! endif
! enddo
! enddo
! do k=1,kmax-1
! enmax(k) = -100.0
! do j=1,jmax
! do i=1,imax
! en = (rho(i,j,k+1)-rho(i,j,k))/dza(k)
! if(k.ge.k1(i,j).and.en.gt.enmax(k))then
! enmax(k) = en
! ien(k) = i
! jen(k) = j
! endif
! enddo
! enddo
! enddo

      print*,'max and min T at kmax ',tmax,tmin
      print*,'max and min T at k=1  ',tmax1,tmin1
      print*,'<T>, <S> ..., <T**2>, <S**2> ...' ,(sum(l)/vol,l=1,lmax),(sumsq(l)/vol,l=1,lmax)
! print*,'ts(1,3,3,4)'
      print*,'ts(1,1,1,1)' ,ts(1,1,1,1)
      print*,'Cn Kh*dt/dx**2 Kv*dt/dz**2' ,cnmax,diff(1)/min(ds,dphi)**2 *dt(kmax), diff(2)/dz(kmax)**2 *dt(kmax)
      print*,'Dmax',dmax
      print*,'flux limited at ',limps,' points'
      print*,'max absolute velocities' ,((imm(i,l),i=1,3),umax(l),l=1,3)
      print*,'Peclet number exceeds 2  ',icp,' times'
! print*,'Peclet number exceeds 4  ',icp2,' times'
      print*,'Peclet number exceeds 20 ',icp10,' times'
      print*,'average Peclet number ',avp/ntot
      print*,'maximum horizontal Peclet nos. ', umax(1)*dphi/diff(1),umax(2)*ds/diff(1)
      print*,'r.m.s. horiz. flow speed ',sqrt(avs/ntot)
      print*,'max barotropic velocities',(ubmax(i),i=1,2)
! print*,'min top to bottom density diff',rdmn,ird,jrd
! write(6,'(i4,e12.5,2i4)')(k,enmax(k),ien(k),jen(k),k=1,12)
! isl=isles+1 path is purely for testing consistency of scheme
! do isl=1,isles+1
      do isl=1,isles
         call island(ub,tv,isl,1)
         print*,'path integral error on island ',isl,tv
      enddo

! test to find appropriate forcings

! do k=16,14,-1
! do j=1,jmax
! tmp(k,j) = 0
! do i=1,imax
! tmp(k,j) = tmp(k,j) + ts(1,i,j,k)
! enddo
! tmp(k,j) = tmp(k,j)/imax
! print*,k,j,tsa(1,1,j)-tmp(k,j),tmp(16,j)-tmp(k,j)
! enddo
! enddo

! check inversion for barotropic streamfunction

! n = imax-1
! do 50 i=1,imax-1
! do 50 j=1,jmax-1
! tv = 0
! k = i + (j-1)*(imax-1)
! if(j.gt.1)
! 1         tv = tv + gapold(k,2)*gb(k-n)
! if(i.gt.1)
! 1         tv = tv + gapold(k,n+1)*gb(k-1)

! tv = tv + gapold(k,n+2)*gb(k)

! if(i.lt.imax-1)
! 1         tv = tv + gb(k+1)*gapold(k,n+3)
! if(j.lt.jmax-1)
! 1         tv = tv + gb(k+n)*gapold(k,2*n+2)
! tv = tv - gbold(k)
! if(abs(tv).gt.1e-7)print*,i,j,k,tv,gbold(k)
! 50 continue

      end
