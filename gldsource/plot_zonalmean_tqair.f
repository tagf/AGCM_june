	subroutine plot_zonalmeans (tair,MonthNumber)
! 
! Sb plots zonal mean T air
!

      integer i,j,imax,jmax,MonthNumber
	parameter(imax=72,jmax=72)
	real tair(imax,jmax),tair_zonalmean(jmax)
	real sv(jmax),s(jmax),alat(jmax)
      real pi,th0,th1,s0,s1,ds,sum
! get model latitudes
	pi=4.*atan(1.0)
	th0 = - 0.5*pi
	th1 = 0.5*pi
	s0 = sin(th0)
	s1 = sin(th1)
	ds = (s1-s0)/jmax
	do  j=1,jmax
	  sv(j) = s0 + j*ds
        s(j) = sv(j) - ds/2.
	  alat(j)=(180./pi)*asin(s(j))
      enddo
!	print*, 'zonal mean tair'
	do  j=1,jmax
	 sum=0.
	  do  i=1,imax
	   sum=sum+tair(i,j)
        enddo
	 tair_zonalmean(j)=sum/imax
 	write (266,'(i3,2f7.2)') MonthNumber,alat(j),tair_zonalmean(j)
      enddo
	 ! print*, alat(j),tair_zonalmean(j)

	return
	end
