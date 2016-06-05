! DataForSurAtm.f - GLDSTN Atmosphere data for surfer lat lon  /10/2008
	!map  -260 - +180 longitude    ####################
	!Точки сетки расположены в полуцелых значениях:
!i,j,=0, 0.5, 1.5, …,35.5, 36.0(или 72)- чтобы границы ячеек проходили
! через целые значения.

      subroutine DataForSurAtm(lout,id_mnth,plot,name)

      include 'var.f90'

      character (len=*) lout
      character (len=*) name
      character (len=*) id_mnth !=Jn or Jl ####
      real sum, tv2, syr, arg
      parameter(syr = 365*86400)
      integer i, j, k, iwets, iice
	real dlon
	real plot(imax,jmax),plot1(0:imax+1,0:jmax+1),tmin,tmax

        sum = 0.
        tmin=100.
        tmax=-100.
      do j=1,jmax
       do i=1,imax
	      plot1(i,j)=plot(i,j)
            if( plot(i,j).lt.tmin) tmin=plot(i,j)
            if( plot(i,j).gt.tmax) tmax=plot(i,j)
            sum = sum +plot(i,j)
       enddo
      enddo

       do i=0,imax+1
	  plot1(i,0)=plot1(i,1) !south pole
	  plot1(i,jmax+1)=plot1(i,jmax) !north pole
       enddo
       do j=0,jmax+1
	  plot1(0,j)=plot1(1,j)  !left border
	  plot1(imax+1,j)=plot1(imax,j) ! right border
       enddo

      open(20,file=trim(path_results)//lout//'.'//trim(name)//id_mnth)

      write(20,*) 'gldst: ',trim(name)
! write(6,*) 'gldst: ',trim(name)
      write(20,*)'min = ',tmin,'  max = ',tmax,'average ', sum/imax/jmax
! print*, 'min = ',tmin,'  max = ',tmax,'average ', sum/imax/jmax
!      write(20,*)'average ', sum/imax/jmax
!      write(6,*)'average ', sum/imax/jmax
      write(20,*)' sdedy = ',sdedy,' sdeyr = ',sdeyr

      dlon=360./real(imax)

      do j=0,jmax+1   !from -260 to +100 lon
	 arg=amin1(2./jmax*(j-0.5)-1.,1.)
	 if (j==0) arg=-1.
	 if (j==jmax+1) arg=1.

        write(20,1 ) dlon*0.-260.,180.*asin(arg)/pi,plot1(0,j)
       do i=1,imax
         write(20,1 ) dlon*(i-0.5)-260.,180.*asin(arg)/pi,plot1(i,j)
       enddo
        write(20,1 ) dlon*imax-260.,180.*asin(arg)/pi,plot1(imax+1,j)
      enddo
   1  format (1x, f7.2,',',f7.2,',',f9.3)

      close(20)
	return
      end
