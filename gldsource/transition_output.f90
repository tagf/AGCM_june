			subroutine TransitionDataOutForSurn(plot,imax,jmax,id_mnth)
			!перевод атм данных с повернутой оси на старую для surfer
			!запись в файл: долгота, широта, температура
			!в дир console1
		 integer lon_old(imax,jmax),lat_old(imax,jmax)
		 real plot(imax,jmax)
		  character (len=*) id_mnth !=Jn or Jl ####

		 open(62,file='../lon_map_Rot.dat', status='old') !after rotation
		 open(63,file='../lat_map_Rot.dat', status='old') !after rotation
		   do j=jmax,1,-1 !for data presentation on old map
			 !чтение  долгота, широта old map, соответсвующих
			!т. (i,j) новой карты
			! эти данные - из программы в файле data_for_map260.f в dir NewMap
			read(62,'(80f7.1)')(lon_old(i,j),i=1,imax)
			 read(63,'(80f7.1)')(lat_old(i,j),i=1,imax)
		   enddo
		   close (62)
		   close (63)

		  open(64,file='../TairRot.dat'//id_mnth, status='new') !after rotation
		 do j=1,jmax   !from -260 to +100 lon
		   do i=1,imax
			 write(64,1 ) lon_old(i,j),lat_old(i,j),plot(i,j)
		   enddo
		  enddo
		   close (64)
	   1  format (1x, f7.2,',',f7.2,',',f9.3)
		  return
		  end
