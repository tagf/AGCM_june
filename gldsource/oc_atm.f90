	subroutine oc_atm(yoc,yatm,imax,jmax)
       ! interpolate yoc(0:imax+1,0:jmax+1) to yatm(1:74,1:46)
	real yoc(0:imax+1,0:jmax+1),yatm(1:74,1:46), yshift(0:imax+1,0:jmax+1),&
	 ytemp(1:74,0:jmax+1),ylon(0:imax+1),ylat(0:jmax+1), slat(1:46)
      real x, ynewgrid, pi, sinlat
	integer i,j

	yshift(1:29,:)=	yoc(8:imax,:)
	yshift(30:imax,:)=yoc(1:7,:)
	yshift(0,:)=	yoc(imax,:)
	yshift(imax+1,:)=	yoc(1,:)

      pi=4.*atan(1.0)

	do j=0,jmax+1
	  ylon(:)= yshift(:,j)
	  do i=2,73
	    x=360./72.*(i-2)-180.
	    ytemp(i,j)=ynewgrid(x,-185.,175.,ylon,imax)
       enddo
     enddo

	do j=1,46
	 slat(j)=sin(pi/45.*(j-1)-pi/2.)
      enddo

	do i=2,73
	  ylat(:)= ytemp(i,:)
	  ylat(0)= ylat(1)
	  do j=1,45
	    yatm(i,j+1)=ynewgrid(slat(j),-1.-2./jmax,1.+2./jmax,ylat, jmax+1)
        enddo
        yatm(i,1)=yatm(i,2)
      enddo

	yatm(1,:)=yatm(73,:)
	yatm(74,:)=yatm(2,:)

	return
	end subroutine

	real function ynewgrid(x,a,b,y,n)
         ! linear interpolation at x point
	real x,a,b,dx,xk,xk1
	real y(0:n+1)
	integer n,k

	if ((x.lt.a).or.(x.gt.b)) then
	  stop 'ynewgrid: x out of range'
      endif

	dx=(b-a)/n
	k=floor((x-a)/dx+1.e-7)+1
!	print *,k
	xk=a+dx*(k-1)
	xk1=xk+dx

	if (abs(x-b).le.1.e-7) then
        ynewgrid=y(n)
      else
	 if (abs(x-xk).ge.1.e-7) then
	  if ((abs(y(k)).ge.1.e-7).and.(abs(y(k+1)).ge.1.e-7)) then
	    ynewgrid=y(k)+(y(k+1)-y(k))*(x-xk)/dx
        else
	     if ((abs(y(k)).ge.1.e-7).and.(abs(y(k+1)).lt.1.e-7)) then
	       ynewgrid=y(k)
	     endif
	     if ((abs(y(k)).lt.1.e-7).and.(abs(y(k+1)).ge.1.e-7)) then
	       ynewgrid=y(k+1)
	     endif
	     if ((abs(y(k)).lt.1.e-7).and.(abs(y(k+1)).lt.1.e-7)) then
	      ynewgrid=0.
	     endif
	  endif
       else
        ynewgrid=y(k)
	 endif
	endif

	return
	end function
