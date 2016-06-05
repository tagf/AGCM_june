
!// subroutine interpolate.f for program gldstn
!// interpolates temperature of ocean on land linearly along latitudes
!// to include in mains:  call interpolate(ts(1,:,:,kmax), k1) before stop 'ts'
!!// result is shifted along altituted by 1 cell t_AGCM(0,:) cells correspond to ONLY_OCEAN longtitude
!!// warning t_AGCM(0, :) cells not affected

    subroutine interpolate(toc,tAGCM) ! toc, k_norm are parameters, t_AGCM = result, shifted by 1 cell

    include 'var.f90'
     
    real toc(0:maxi+1,0:maxj+1), toc_shift(0:maxi,0:maxj+1) ! only 0:maxi dimension for i required in shifted array
	integer k_norm(0:maxi+1,0:maxj+1), k_shift(0:maxi, 0:maxj+1)
    integer i, j
    
    ! Tagir required
    integer i1, i2, normal_shift, strange_shift ! start and end of land on interpolation
    parameter (normal_shift = 18)
    parameter (strange_shift = 0) ! if everything is correct, not required
    logical boolLand
    ! for interpolation on AGCM grid
    real tAGCM(0:maxi+1, 1:46)
    real hj_coord, hj_coord_floor_dist
    real PI_180
    real pole_north, pole_south
    parameter (PI_180 = atan(1.)/45.)
    
    ! oc(1:56,:)=ts_oc_for_atm(17:72,:)             !Shift ts_oc_for_atm grid
    ! oc(57:72,:)=ts_oc_for_atm(1:16,:)
    ! k1(1:56,:)=k2(17:72,:)                        !Shift k1
    ! k1(57:72,:)=k2(1:16,:)

    k_shift(0:maxi-normal_shift,:)=k1(normal_shift:maxi,:) ! Tagir
    k_shift(maxi-normal_shift+1: maxi,:)=k1(1:normal_shift, : ) ! Tagir
  !  k_shift(maxi-normal_shift+1: maxi-1,:)=k1(0:normal_shift-2, : ) ! Tagir
  !  k_shift(maxi, : )=k1(normal_shift, : ) ! only 0..72 array actually required
    
    ! strange_shift = 0
    toc_shift(0:maxi-normal_shift-strange_shift,:)=toc(normal_shift+strange_shift: maxi, : ) ! Tagir
    toc_shift(maxi-normal_shift-strange_shift+1: maxi, : )=toc(1:normal_shift+strange_shift, : ) ! Tagir
  !  toc_shift(maxi-normal_shift-strange_shift+1: maxi-1, : )=toc(0:normal_shift+strange_shift-2, : ) ! Tagir
  !  toc_shift(maxi, : )=toc(normal_shift+strange_shift, : ) ! Tagir
    
    ! without interpolation
    !do j=72,1,-1 ! 
    !  do i=1,72 ! altitude 
    !    write (1446, *) i, j, toc_shift(i,j)
    !  enddo
    !enddo
        
    ! ocean surf type 1 - 8  (ocean to atm interp begin, Tagir)
    do j=1,maxj  
      do i=0,maxi ! altitude 
        if (k_shift(i,j).gt.8) then ! not ocean
          if (boolLand) then ! we've stepped and know bounds i1 - i2 of Land
            ! T(i) = interpolate(i; i1, i2, T1, T2)
            toc_shift(i,j)=toc_shift(i1-1,j)+float(i-i1+1)*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
          else
            i1 = i ! start of land
            boolLand = .true. ! we've stepped onLand flag, for interpolating further
            do i2=i,72 ! we search for the end of land (ocean)
              if ((k_shift(i2,j).ge.1).and.(k_shift(i2,j).le.8)) then ! ocean, found end of land
                ! first cell on land interpolation
                toc_shift(i1,j)=toc_shift(i1-1,j) + 1*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
                ! write (1448, *), i2, j, toc_shift(i2, j) ! purely diagnostics of land-altitude-line end
                exit
              endif
            enddo
          endif ! boolLand
        else ! ocean
          boolLand = .false.
          ! temperature on ocean without changes
        endif ! ocean / not ocean
      enddo
    enddo
    ! ocean to atm interp procedure end (Tagir)

    ! with interpolation    
    !do j=72,1,-1 ! 
    !  do i=1,72 ! altitude 
    !    write (1447, *) i, j, toc_shift(i,j)
    !  enddo
    !enddo

  !!!  interpolation to AGCM grid
  pole_north = sum(toc_shift(1:maxi, 72))/72
  pole_south = sum(toc_shift(1:maxi, 1))/72
  
  tAGCM(:,46) = pole_north
  tAGCM(:,1) = pole_south
  
  do j = 2,3 ! near southern pole
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5
    hj_coord_floor_dist = hj_coord - floor(hj_coord)
    do i=0,maxi
       tAGCM(i+1,j) = (1.-hj_coord_floor_dist) * pole_south + &
                     (   hj_coord_floor_dist) * toc_shift(i, 1) 
    enddo
  enddo 

  do j = 44,45 ! near northern pole
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5
    hj_coord_floor_dist = hj_coord - floor(hj_coord)
    do i=0,maxi
       tAGCM(i+1,j) = (1.-hj_coord_floor_dist) * toc_shift(i, 72) + &
                     (   hj_coord_floor_dist) * pole_north 
    enddo
  enddo
  
  do j=4,43 ! only along longtitudes, 4..43 -- inside gldstn grid
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5 ! [0., 73.]
    hj_coord_floor_dist = hj_coord - floor(hj_coord) ! [0., 1.]
    do i=0,maxi
       tAGCM(i+1,j) = (1.-hj_coord_floor_dist) * toc_shift(i, floor(hj_coord)) + &
                     (   hj_coord_floor_dist) * toc_shift(i, floor(hj_coord) + 1)  
    enddo
  enddo
  
  tAGCM(0,:) = tAGCM(72,:)
  tAGCM(73,:) = tAGCM(1,:)   ! возможно избыточно
  
  ! diagnostics with interpolation on AGCM grid    
  ! do j=46,1,-1 ! 
  !   do i=1,maxi ! altitude
  !     write (1550, *) i, j, t_AGCM(i,j)
  !   enddo
  ! enddo
        
  return
  end subroutine	



  !subroutine interpolate purely diagnostical mod
  subroutine interpolate1(toc, k_norm, t_AGCM, toc_shift) ! toc, k_norm are parameters, t_AGCM = result, shifted by 1 cell
    
    include 'var.f90'
     
    real toc(0:maxi+1,0:maxj+1), toc_shift(0:maxi,0:maxj+1) ! only 0:maxi dimension for i required in shifted array
	integer k_norm(0:maxi+1,0:maxj+1), k_shift(0:maxi, 0:maxj+1)
    integer i, j
    
    ! Tagir required
    integer i1, i2, normal_shift, strange_shift ! start and end of land on interpolation
    parameter (normal_shift = 18)
    parameter (strange_shift = 0) ! if everything is correct, not required
    logical boolLand
    ! for interpolation on AGCM grid
    real t_AGCM(0:maxi+1, 1:46)
    real hj_coord, hj_coord_floor_dist
    real PI_180
    real pole_north, pole_south
    parameter (PI_180 = atan(1.)/45.)
    
    ! oc(1:56,:)=ts_oc_for_atm(17:72,:)             !Shift ts_oc_for_atm grid
    ! oc(57:72,:)=ts_oc_for_atm(1:16,:)
    ! k1(1:56,:)=k2(17:72,:)                        !Shift k1
    ! k1(57:72,:)=k2(1:16,:)

    k_shift(0:maxi-normal_shift,:)=k_norm(normal_shift:maxi,:) ! Tagir
    k_shift(maxi-normal_shift+1: maxi-1,:)=k_norm(0:normal_shift-2, : ) ! Tagir
    k_shift(maxi, : )=k_norm(normal_shift, : ) ! only 0..72 array actually required
    
    ! strange_shift = 0
    toc_shift(0:maxi-normal_shift-strange_shift,:)=toc(normal_shift+strange_shift: maxi, : ) ! Tagir
    toc_shift(maxi-normal_shift-strange_shift+1: maxi-1, : )=toc(0:normal_shift+strange_shift-2, : ) ! Tagir
    toc_shift(maxi, : )=toc(normal_shift+strange_shift, : ) ! Tagir
    
    ! without interpolation
    !do j=72,1,-1 ! 
    !  do i=1,72 ! altitude 
    !    write (1446, *) i, j, toc_shift(i,j)
    !  enddo
    !enddo
        
    ! ocean surf type 1 - 8  (ocean to atm interp begin, Tagir)
    do j=1,maxj ! 
      do i=0,maxi ! altitude 
        if (k_shift(i,j).gt.8) then ! not ocean
          if (boolLand) then ! we've stepped and know bounds i1 - i2 of Land
            ! T(i) = interpolate(i; i1, i2, T1, T2)
            toc_shift(i,j)=toc_shift(i1-1,j)+float(i-i1+1)*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
          else
            i1 = i ! start of land
            boolLand = .true. ! we've stepped onLand flag, for interpolating further
            do i2=i,72 ! we search for the end of land (ocean)
              if ((k_shift(i2,j).ge.1).and.(k_shift(i2,j).le.8)) then ! ocean, found end of land
                ! first cell on land interpolation
                toc_shift(i1,j)=toc_shift(i1-1,j) + 1*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
                ! write (1448, *), i2, j, toc_shift(i2, j) ! purely diagnostics of land-altitude-line end
                exit
              endif
            enddo
          endif ! boolLand
        else ! ocean
          boolLand = .false.
          ! temperature on ocean without changes
        endif ! ocean / not ocean
      enddo
    enddo
    ! ocean to atm interp procedure end (Tagir)

    ! with interpolation    
    !do j=72,1,-1 ! 
    !  do i=1,72 ! altitude 
    !    write (1447, *) i, j, toc_shift(i,j)
    !  enddo
    !enddo

  !!!  interpolation to AGCM grid
  pole_north = sum(toc_shift(1:maxi, 72))/72
  pole_south = sum(toc_shift(1:maxi, 1))/72
  
  t_AGCM(:,46) = pole_north
  t_AGCM(:,1) = pole_south
  
  do j = 2,3 ! neart southern pole
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5
    hj_coord_floor_dist = hj_coord - floor(hj_coord)
    do i=0,maxi
       t_AGCM(i+1,j) = (1.-hj_coord_floor_dist) * pole_south + &
                     (   hj_coord_floor_dist) * toc_shift(i, 1) 
    enddo
  enddo 

  do j = 44,45 ! neart northern pole
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5
    hj_coord_floor_dist = hj_coord - floor(hj_coord)
    do i=0,maxi
       t_AGCM(i+1,j) = (1.-hj_coord_floor_dist) * toc_shift(i, 72) + &
                     (   hj_coord_floor_dist) * pole_north 
    enddo
  enddo
  
  do j=4,43 ! only along longitudes, 4..43 -- inside gldstn grid
    hj_coord = (sin(float(4*(j - 1) - 90) * PI_180) + 1.) * 73. * 0.5 ! [0., 73.]
    hj_coord_floor_dist = hj_coord - floor(hj_coord) ! [0., 1.]
    do i=0,maxi
       t_AGCM(i+1,j) = (1.-hj_coord_floor_dist) * toc_shift(i, floor(hj_coord)) + &
                     (   hj_coord_floor_dist) * toc_shift(i, floor(hj_coord) + 1)  
    enddo
  enddo
  
  ! diagnostics with interpolation on AGCM grid    
  ! do j=46,1,-1 ! 
  !   do i=1,maxi ! altitude
  !     write (1550, *) i, j, t_AGCM(i,j)
  !   enddo
  ! enddo
        
  return
  end subroutine	


!// subroutine interpolate_back for program gldstn
!// interpolates temperature of AGCM grid to gldstn grid
!!// warning : some border cells are not affected toc(maxi+1,:) toc(0,maxj+1,:)

    subroutine interpolate_back(t_AGCM, toc) ! t_AGCM parameter, toc = result, no shifts
    
    include 'var.f90'
     
    real toc(1:maxi,1:maxj), toc_shift(1:maxi,1:maxj) !! warning : some border cells are not affected 
    real t_AGCM(1:maxi,1:46)
    ! Tagir required
    integer normal_shift
    parameter (normal_shift = 17)
    real hj_coord, hj_coord_floor_dist
    real invPI_180
    parameter (invPI_180 = 45./atan(1.))
            
    ! interpolating back from AGCM to gldstn grid
    integer i,j
    do j=1,maxj !
      hj_coord =  1. + 0.25*(90. + (asin(-1. + 2.*float(j)/73.) * invPI_180)) ! (3, 44)
      hj_coord_floor_dist = hj_coord - floor(hj_coord)
      do i=1,maxi ! longitude 
        toc(i,j) = (1.-hj_coord_floor_dist) * t_AGCM(i, floor(hj_coord)) + &
                   (   hj_coord_floor_dist) * t_AGCM(i, floor(hj_coord) + 1)   
      enddo
    enddo
    ! end of back interpolation
    
    ! strange_shift = 0
    toc_shift(normal_shift: maxi, : )=toc(1:maxi-normal_shift+1,:) ! Tagir
    toc_shift(1:normal_shift-1, : )=toc(maxi-normal_shift+2: maxi, : ) ! Tagir
    toc=toc_shift
! back interpolation diagnostics   
  !do j=72,1,-1 ! 
  !  do i=0,72 ! altitude 
  !    write (1521, *) i, j, toc(i,j)
  !  enddo
  !enddo
    
  return
  end subroutine	

!// subroutine interpolate_diag for program gldstn
    subroutine interpolate_diag() ! t_AGCM parameter, toc = result, no shifts
        !// ts and k1 are defined in var.f90
        include 'var.f90'

    	!// warning: in my interp subroutenes some cells are not affected, and shift presented!
        real t_AGCM(0:maxi+1, 1:46) ! Tagir for interpolation on AGCM grid
        real toc(0:maxi+1,0:maxj+1) ! Tagir for interpolation back on gld grid
        
        ! temporary diagnostics
        real toc_shift(0:maxi,0:maxj+1)
        integer i,j

            
        !// Dignostics temperature and quit point 
        do j=1,maxi ! instead of imax
            do i=1,maxj  ! instead of jmax
                write(1531,*) i,j,k1(i,j),ts(1,i,j,8)
            enddo
        enddo

        !// Diagnostics interpolation 
        call interpolate1(ts(1,:,:,kmax), k1, t_AGCM, toc_shift)
        
        ! with interpolation on AGCM grid    
        do j=46,1,-1 ! 
          do i=1,maxi ! altitude 
            write (1532, *) i, j, t_AGCM(i,j)
          enddo
        enddo
        
        ! back interpolation
        call interpolate_back(t_AGCM, toc)
        ! back interpolation diagnostics   
        do j=72,1,-1 ! 
          do i=0,72 ! altitude 
            write (1533, *) i, j, toc(i,j)
          enddo
        enddo
  
        ! diagnostics of differenceA
        do j=72,1,-1 ! 
          do i=0,72 ! altitude 
            write (1534, *) i, j, toc(i,j) - ts(1,i,j,kmax)
          enddo
        enddo
        
        ! diagnostics of differenceB    
        do j=72,1,-1 ! 
          do i=0,72 ! altitude 
            write (1535, *) i, j, toc(i,j) - toc_shift(i,j)
          enddo
        enddo
    
  return
  end subroutine	