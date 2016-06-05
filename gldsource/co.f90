! co.f convection code simplified form for program gldstn
! suitable for arbitrary functions rho(T,S) variable depth
! cost counts occurences of mixing at each point not including the point
! at the top of each mixed region
! if cost is 2-d it counts the average number of convecting points at
! each horizontal point (divide by nsteps in mains.f)
! lmax.gt.2 allowed
! ts array passed as argument

      subroutine co(tv)

      include 'var.f90'

      real dzm(maxk), sum(maxl), tv(maxl,0:maxi+1,0:maxj+1,0:maxk+1)

      integer i, j, k(0:maxk), l, lastmix, m, n, ni

      do 10 j=1,jmax
         do 10 i=1,imax

! initialize the index array k and mixed region sizes dzm
! wet points only
           if(k1(i,j).le.kmax)then

            k(k1(i,j)-1) = 0
            do 20 m=k1(i,j),kmax
               k(m) = m
               dzm(m) = dz(m)
   20       continue

            m = kmax
            lastmix = 0

! main loop 'normally' decreasing in m

            do 30 while (k(m-1).gt.0.or.(lastmix.ne.0.and.k(m).ne.kmax))

               if(rho(i,j,k(m)).lt.rho(i,j,k(m-1)).or.k(m-1).eq.0)then
! this may need changing, unless as rho(i,j,0) dimensioned
                  if(lastmix.eq.0.or.k(m).eq.kmax)then
                     m = m-1
                  else
                     m = m+1
                  endif
                  lastmix = 0
               else
                  lastmix = 1

! look for instability before mixing

                  n = m-1
                  do 40 while (k(n-1).gt.0.and. rho(i,j,k(n)).ge.rho(i,j,k(n-1)))
                     n = n-1
   40             enddo
                  do 80 l=1,lmax
                  sum(l) = tv(l,i,j,k(m))*dzm(k(m))
   80             continue
                  do 60 ni=1,m-n
                     do 70 l=1,lmax
                     sum(l) = sum(l) + tv(l,i,j,k(m-ni))*dzm(k(m-ni))
   70                continue
                     dzm(k(m)) = dzm(k(m)) + dzm(k(m-ni))
   60             continue
                  do 90 l=1,lmax
                     tv(l,i,j,k(m)) = sum(l)/dzm(k(m))
   90             continue
                  rho(i,j,k(m)) = ec(1)*tv(1,i,j,k(m)) + ec(2)*tv(2,i,j,k(m)) + ec(3)*tv(1,i,j,k(m))**2 + ec(4)*tv(1,i,j,k(m))**3
! reindex k(m)
                  ni = m-1
                  do 50 while (k(ni+1).gt.0)
                     k(ni) = k(ni-m+n)
                     ni = ni-1
   50             enddo
               endif
   30       enddo

! fill in T,S values in mixed regions

            m = kmax-1
            do 100 n=kmax-1,k1(i,j),-1
               if(n.gt.k(m))then
                  do 110 l=1,lmax
                     tv(l,i,j,n) = tv(l,i,j,k(m+1))
  110             continue
                  rho(i,j,n) = ec(1)*tv(1,i,j,n) + ec(2)*tv(2,i,j,n) + ec(3)*tv(1,i,j,n)**2 + ec(4)*tv(1,i,j,n)**3
                  cost(i,j) = cost(i,j) + 1.0
! cost(i,j,n) = cost(i,j,n) + 1.0
               else
                  m = m-1
               endif
  100       continue
! wet points only
        endif
   10 continue
      end