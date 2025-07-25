  implicit none

   integer, parameter :: nx=177
   integer, parameter :: ny=168
   integer, parameter :: nh=24
   integer, parameter :: recl_f=1*nx*ny  ! for intel ifort
!   integer, parameter :: recl_f=4*nx*ny  ! for gnu gfortran
 
   integer, parameter :: ip=62  ! x coord. of phx
   integer, parameter :: jp=81  ! y coord. of phx

!======
    real,dimension(nx,ny):: temp,tmin,tmax,dt,tav
    real,dimension(nh):: rt_phx
    integer :: ih,ihh,ix,iy
  character(len=3) cenar

  cenar="fdx"
!  cenar="0dx"
!===============================================================
!===============================================================
    open(44,file='temp_'//cenar//'.dat',status='old'  &
     ,form='unformatted',access='direct',recl=recl_f)
        
     tmin=99999.
     tmax=-99999.
     tav=0.

     DO ih=1,nh

        read(44,rec=ih) temp

        rt_phx(ih)=temp(ip,jp)

        do iy=1,ny
        do ix=1,nx

           if(temp(ix,iy).gt.tmax(ix,iy)) tmax(ix,iy)=temp(ix,iy)
           if(temp(ix,iy).lt.tmin(ix,iy)) tmin(ix,iy)=temp(ix,iy)
           tav(ix,iy)=tav(ix,iy)+temp(ix,iy)

        enddo
        enddo

      ENDDO

      tav=tav/24.
      close(44)

      dt=tmax-tmin
!!======================================
     open(44,file='tmin_max_dt_tav_'//cenar//'.dat',status='unknown'  &
      ,form='unformatted',access='direct',recl=recl_f)

       write(44,rec=1) tmin
       write(44,rec=2) tmax
       write(44,rec=3) dt
       write(44,rec=4) tav

       close(44)
!===============================================================
!===============================================================

       open(99,file='phx_temp_'//cenar//'.dat',status='unknown')

        do ih=1,nh
          write(99,88) float(ih-1),rt_phx(ih)
        enddo

       close(99)

 88    format( 2f10.2)
! 88    format(i12, 6f16.7)


!
      stop
      end
!============================


