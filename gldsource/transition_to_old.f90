			subroutine TransitionData(plot,imax,jmax)
			!перевод атм данных с повернутой оси на старую для surfer
			! № точки долгота, № точки широта, температура
			! - неточное представление получается
		 integer i_map(imax,jmax),j_map(imax,jmax)
		 real plot(imax,jmax)

		 open(60,file='../i_map_Rot.dat', status='old') !after rotation
		 open(61,file='../j_map_Rot.dat', status='old') !after rotation
		   do j=jmax,1,-1 !for data presentation on old map
			!чтение № точки долгота, № точки широта old map, соответсвующих
			!т. (i,j) новой карты
			! эти данные - из программы в файле data_for_map260.f в dir NewMap
			 read(60,'(80i3)')(i_map(i,j),i=1,imax)
			 read(61,'(80i3)')(j_map(i,j),i=1,imax)
		   enddo
		   close (60)
		   close (61)
		  do i=1,imax
			do j=1,jmax
			 plot(i,j)=plot(i_map(i,j),j_map(i,j))
			enddo
		   enddo
		  return
		  end
