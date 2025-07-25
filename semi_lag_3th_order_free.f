
  implicit none


  integer, parameter :: nx=177
  integer, parameter :: ny=168
  integer, parameter :: nh=24
   integer, parameter :: recl_f=1*nx*ny  ! for intel ifort
!   integer, parameter :: recl_f=4*nx*ny  ! for gnu gfortran


  integer, parameter :: ip=62  ! x coord. of phx
  integer, parameter :: jp=81  ! y coord. of phx
  real, parameter :: rz0=341.7279  ! this is hgt(ip,jp)

  real, parameter :: rth_0=308.   ! at z=0
  real, parameter :: rn2=1.*1.e-4

  real, parameter :: dx=2000.
  real, parameter :: dy=2000.
  integer, parameter :: idt=36
  real, parameter :: dt=float(idt)
  real, parameter :: dth=dt/3600.
  integer, parameter :: nt=3*24*3600/idt+2
  real, parameter :: rg=9.81
  real, parameter :: rgaz=287.
  real, parameter :: rcp=1004.5

  real, parameter :: acoe=0.00012
  real, parameter :: acoep=1+acoe*dt/2.
  real, parameter :: acoem=1-acoe*dt/2.
  character(len=3) cenar
!===========================================================
!===========================================================

  real,dimension (nx,ny,nh+1):: ru,rv
  real,dimension (nh+1):: rh
  real,dimension (nh+1):: th_wrf_u,th_wrf_r
  real,dimension(nx,ny):: alon,alat,rlui,hgt,frac_u,frac_r
  real,dimension(nx):: rx
  real,dimension(ny):: ry


  real,dimension(nx,ny):: rth1,rth2,rtemp
  real,dimension(nx,ny):: ruu,rvv,uindex
  real:: rth_u,rth_r,rz,dist,th_u,th_r
  integer ix,iy,it,ih,ihh
  integer irec_ar,irec,ih1,ih2
  real :: time,time_u
  real :: rpi_u,rpi_r
  real :: xf,yf,xi,yi,thf
  integer ixf,iyf
!===========================================================
!===========================================================
  open(77,file="forcing_theta_rz0.d",status="old")
!=================================================
  DO ih=1,nh
   ! here th_wrf_u is from observation at rz0 as generated in
   ! forcing_theta_rz0.d
    read(77,*)  ihh,rh(ih),th_wrf_u(ih),th_wrf_r(ih) 
  ENDDO
  close(77)


  th_wrf_u(nh+1)=th_wrf_u(1)
  th_wrf_r(nh+1)=th_wrf_r(1)
  rh(nh+1)=24.

!===========================================================
!===========================================================
!===========================================================


  irec_ar=0
    

!=================================================
  open(32,file='WRF_DATA/landuse_fdx.dat',status='old'  &
     ,form='unformatted',access='direct',recl=recl_f)
        read(32,rec=1) alon
        read(32,rec=2) alat
        read(32,rec=3) rlui
        read(32,rec=4) hgt
   close(32)

   open(32,file='WRF_DATA/landuseF_fdx.dat',status='old' &
     ,form='unformatted',access='direct',recl=recl_f)
       read(32,rec=1) frac_r
       read(32,rec=2) frac_u
   close(32)

!=================================================
    uindex=0.

    do iy=1,ny
    do ix=1,nx
      if ((rlui(ix,iy).ge.23).and.(rlui(ix,iy).le.26.)) &
       uindex(ix,iy)=1.
    enddo
    enddo

!================= scenarios================================
!     do iy=1,ny
!     do ix=1,nx
!     dist=sqrt(float(ix-ip)**2+float(iy-jp)**2)
!     if (dist.gt.2.) then 
!         uindex(ix,iy)=0.
!         frac_u(ix,iy)=0.
!         frac_r(ix,iy)=1.
!    endif
!    enddo
!    enddo
!======================== all rural =========================
!======================== all rural =========================
!    uindex=0.  ! no urb
!    frac_u=0.  ! no urb
!    frac_r=1.  ! no urb
!

!    open(55,file='rt_sim_fdx.dat',status='unknown' &    ! 

    cenar="fdx"
    open(55,file='temp_'//cenar//'.dat',status='unknown' &    !
     ,form='unformatted',access='direct',recl=recl_f)
!=================================================
!=================================================

!=================================================
    do ih=1,nh
      rh(ih)=float(ih-1)
    enddo
    rh(nh+1)=float(nh)
8  format(i8,5f12.4)

!=================================================
!=================================================


7   format(i8,3f12.4)
!=================================================

!=================================================
!=================================================
!=================================================
    open(32,file='WRF_DATA/ru_fdx_wrf.dat',status='old'  &
     ,form='unformatted',access='direct',recl=recl_f)

      irec=0
      do ih=1,nh
        irec=irec+1
        read(32,rec=irec) ru(:,:,ih)
      enddo
    close(32)

    do iy=1,ny
    do ix=1,nx
      ru(ix,iy,nh+1)=ru(ix,iy,1)
    enddo
    enddo
!=================================================
    open(32,file='WRF_DATA/rv_fdx_wrf.dat',status='old' &
     ,form='unformatted',access='direct',recl=recl_f)
      irec=0
      do ih=1,nh
        irec=irec+1
        read(32,rec=irec) rv(:,:,ih)
      enddo
    close(32)
    do iy=1,ny
    do ix=1,nx
      rv(ix,iy,nh+1)=rv(ix,iy,1)
    enddo
    enddo
!!=================================================
  print*,dt,nt

  time=0.
!
  do iy=1,ny
  do ix=1,nx
    rz=hgt(ix,iy)
    rth1(ix,iy)=rth_0*exp(rn2/rg*(rz-rz0))
  enddo
  enddo
  rth2=rth1
!
!
  do ix=1,nx
    rx(ix)=float(ix-1)*dx
  enddo
  do iy=1,ny
    ry(iy)=float(iy-1)*dy
  enddo
!
  DO it=1,nt
!
!
!
      time_u=time+dth/2.
      time=time+dth
      if (time.ge.24.) time=time-24.
      if (time_u.ge.24.) time_u=time_u-24.
!
      ih1=int(time)+1
      ih2=ih1+1
!
!
      call inter_uv(ru,ruu,rh,time_u,nx,ny,nh) 
      call inter_uv(rv,rvv,rh,time_u,nx,ny,nh)

!
      rth2=rth1
!
!
      do iy=3,ny-2
      do ix=3,nx-2

        xf=rx(ix)
        yf=ry(iy)
        ixf=ix
        iyf=iy
        thf=rth2(ix,iy)
        call departure(ruu,rvv,rx,ry,xf,yf,xi,yi,ixf,iyf,dt,dx,dy,nx,ny)
        call semi_temp_lag(rth1,thf,rx,ry,xi,yi,dx,dy,nx,ny)

        rth2(ix,iy)=thf
!
      enddo
      enddo
!

      call inter_th(th_wrf_u,rth_u,rh,time,nh) 
      call inter_th(th_wrf_r,rth_r,rh,time,nh) 


      do iy=1,ny
      do ix=1,nx

         rz=hgt(ix,iy)
         th_u=rth_u*exp(rn2/rg*(rz-rz0))
         th_r=rth_r*exp(rn2/rg*(rz-rz0))

!
        rth2(ix,iy)=  &
       (acoem*rth2(ix,iy)+dt*(acoe*(frac_u(ix,iy)*th_u  &
          +frac_r(ix,iy)*th_r)))/acoep

!        
      enddo
      enddo
!!

      if(mod(it-1,100)==0) print 3,it,ih1,rth1(ip,jp)  &
      ,rth_r,rth_u,rn2*1.e4,time-dth
!
      IF(it.ge.4801.and.mod(it-1,100)==0) THEN

        call theta_to_temp(rtemp,rth1,hgt,rg,rcp,rn2,nx,ny)

        irec_ar=irec_ar+1
!        write(55,rec=irec_ar) rth1
        write(55,rec=irec_ar) rtemp
        print *, 'ARCH  ',irec_ar,it,ih1,rth2(ip,jp),time-dth

      ENDIF

3     format(2i10,5f12.4) 
!
      rth1=rth2
!
  ENDDO
!
   close(55)
      

  stop
  end
!==============================
!==============================
  subroutine inter_th(th_wrf,rth,rh,time,nh) 

    implicit none
    integer , intent(in) :: nh
    real , intent(in) :: time
    real,dimension (nh+1), intent(in) :: rh,th_wrf
    real , intent(out) :: rth

    integer :: ih,ih1,ih2

    do ih=1,nh
  
      if ( (time.ge.rh(ih)).and.(time.lt.rh(ih+1)) ) then
!             ih1=int(time)+1
         ih1=ih
         ih2=ih1+1
         go to 40
      endif

    enddo

40    continue


    rth=th_wrf(ih1)+(th_wrf(ih2)-th_wrf(ih1)) &
    *(time-rh(ih1))/(rh(ih2)-rh(ih1))


  end subroutine inter_th
!==============================
!==============================
  subroutine inter_uv(ru,ruu,rh,time,nx,ny,nh)

    implicit none
    integer, intent(in) :: nx,ny,nh
    real,dimension (nx,ny,nh+1), intent(in):: ru
    real,dimension (nh+1), intent(in):: rh
    real,dimension(nx,ny), intent(out):: ruu
    real , intent(in):: time

    integer :: ih1,ih2,ix,iy

    ih1=int(time)+1
    ih2=ih1+1

    do iy=1,ny
    do ix=1,nx

      ruu(ix,iy)=ru(ix,iy,ih1)+(ru(ix,iy,ih2)-ru(ix,iy,ih1))  &
      *(time-rh(ih1))/(rh(ih2)-rh(ih1))

    enddo
    enddo

  end subroutine inter_uv
!==============================
!==============================

  subroutine departure(ru,rv,rx,ry,xf,yf,xi,yi,ixf,iyf,dt,dx,dy,nx,ny)

    implicit none
    integer :: ixf,iyf,nx,ny
    real,dimension (nx,ny):: ru,rv
    real,dimension(nx):: rx
    real,dimension(ny):: ry
    real :: xf,yf,xi,yi,dt
    real :: dx,dy,uf,vf,acon
    real :: ui1,ui2,ui,vi1,vi2,vi
    integer :: ir,ix1,ix2,iy1,iy2



    xi=xf
    yi=yf

    uf=ru(ixf,iyf)
    vf=rv(ixf,iyf)


    do ir=1,5

      acon=(xi-rx(1))/dx+0.000001
      ix1=int(acon)+1
      ix2=ix1+1

      acon=(yi-ry(1))/dy+0.000001
      iy1=int(acon)+1
      iy2=iy1+1

!!!!   u adv
        ui1=ru(ix1,iy1)+(ru(ix2,iy1)-ru(ix1,iy1)) &
      *(xi-rx(ix1))/(rx(ix2)-rx(ix1))
    
          ui2=ru(ix1,iy2)+(ru(ix2,iy2)-ru(ix1,iy2)) &
      *(xi-rx(ix1))/(rx(ix2)-rx(ix1))

      ui=ui1+(ui2-ui1)*(yi-ry(iy1))/(ry(iy2)-ry(iy1))
!!!!   v adv
      vi1=rv(ix1,iy1)+(rv(ix2,iy1)-rv(ix1,iy1)) &
      *(xi-rx(ix1))/(rx(ix2)-rx(ix1))

      vi2=rv(ix1,iy2)+(rv(ix2,iy2)-rv(ix1,iy2)) &
      *(xi-rx(ix1))/(rx(ix2)-rx(ix1))

      vi=vi1+(vi2-vi1)*(yi-ry(iy1))/(ry(iy2)-ry(iy1))


      xi=xf-dt/2.*(uf+ui)
      yi=yf-dt/2.*(vf+vi)

!   print *,xi/2000.,yi/2000.,ix1,iy1,ui,vi

    enddo

  end subroutine departure
!==============================
!==============================
  subroutine semi_temp_lag(rt,tf,rx,ry,xi,yi,dx,dy,nx,ny)

    implicit none
    integer ,intent(in) :: nx,ny
    real,dimension (nx,ny),intent(in) :: rt
    real,dimension(nx),intent(in) :: rx
    real,dimension(ny),intent(in) :: ry
    real ,intent(in) :: xi,yi
    real ,intent(in) :: dx,dy 
    real ,intent(out) :: tf
!=========================
    real   :: dd,adom,anum,fmin,fmax
    real   rf_t(4,4)
    real  rg_t(4)
    real  rh_t
    real  x(4),y(4),coefx(4),coefy(4)
!=========================
    integer :: ix,ii,is   ! is=1--> 4,  ii=ix-1-->ix+2
    integer :: iy,jj,js   ! js=1--> 4,  jj=iy-1-->iy+2
!=========================




!======================================
    dd=(xi-rx(1))/dx+1.0e-6
    ix=int(dd)+1
    ix =max(ix,2)
    ix =min(ix,nx-2)

    dd=(yi-ry(1))/dy+1.0e-6
    iy=int(dd)+1
    iy =max(iy,2)
    iy =min(iy,ny-2)
!======================================
!======================================

    if((ix.lt.2).or.(ix.gt.nx-2).or.(iy.lt.2).or.(iy.gt.ny-2)) then 
  
        print *, 'THERE IS A PROB'
        print *, ix,iy
        stop
    endif
!=========================================
!======================================
    DO jj=iy-1,iy+2
       js=jj-(iy-1)+1
       y(js)=ry(jj)
    ENDDO

    DO ii=ix-1,ix+2
       is=ii-(ix-1)+1
       x(is)=rx(ii)
    ENDDO
!======================================

!======================================
    DO jj=iy-1,iy+2  ! start jj loop

       js=jj-(iy-1)+1
  
       DO ii=ix-1,ix+2    ! start ii loop
         is=ii-(ix-1)+1
         rf_t(is,js)=rt(ii,jj)
       ENDDO  ! end ii loop

    ENDDO  ! end jj loop
!======================================

    adom=(y(1)-y(2))*(y(1)-y(3))*(y(1)-y(4))
    anum=(yi-y(2))*(yi-y(3))*(yi-y(4))
    coefy(1)=anum/adom

    adom=(y(2)-y(1))*(y(2)-y(3))*(y(2)-y(4))
    anum=(yi-y(1))*(yi-y(3))*(yi-y(4))
    coefy(2)=anum/adom

    adom=(y(3)-y(2))*(y(3)-y(1))*(y(3)-y(4))
    anum=(yi-y(2))*(yi-y(1))*(yi-y(4))
    coefy(3)=anum/adom
  
    adom=(y(4)-y(2))*(y(4)-y(3))*(y(4)-y(1))
    anum=(yi-y(2))*(yi-y(3))*(yi-y(1))
    coefy(4)=anum/adom
!======================================

    adom=(x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))
    anum=(xi-x(2))*(xi-x(3))*(xi-x(4))
    coefx(1)=anum/adom

    adom=(x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))
    anum=(xi-x(1))*(xi-x(3))*(xi-x(4))
    coefx(2)=anum/adom

    adom=(x(3)-x(2))*(x(3)-x(1))*(x(3)-x(4))
    anum=(xi-x(2))*(xi-x(1))*(xi-x(4))
    coefx(3)=anum/adom

    adom=(x(4)-x(2))*(x(4)-x(3))*(x(4)-x(1))
    anum=(xi-x(2))*(xi-x(3))*(xi-x(1))
    coefx(4)=anum/adom
  
!=========================================
    DO js=1,4

      rg_t(js)=0.

      DO is=1,4
          rg_t(js)=rg_t(js)+coefx(is)*rf_t(is,js)
      ENDDO

    ENDDO  ! end jj loop

    rh_t =0.
    DO js=1,4
             rh_t=rh_t+coefy(js)*rg_t(js)
    ENDDO


    fmin=minval(rf_t)
    fmax=maxval(rf_t)
    tf=max(fmin,min(fmax,rh_t))


  end subroutine semi_temp_lag
!============================================================
!============================================================
  subroutine theta_to_temp(temp,theta,hgt,rg,rcp,rn2,nx,ny)

    implicit none
    integer ,intent(in) :: nx,ny
    real ,intent(in) :: rg,rcp,rn2
    real,dimension (nx,ny),intent(in) :: theta,hgt
    real,dimension (nx,ny),intent(out) :: temp


    integer :: ix,iy
    real :: rz,fac,rpi

        do iy=1,ny
        do ix=1,nx

        rz=hgt(ix,iy)

        fac = ( 1.-exp(rn2*rz/rg) )/theta(ix,iy)

        rpi=1.+ rg*rg/(rcp*rn2) * fac

        temp(ix,iy)=(theta(ix,iy)*rpi ) - 273.15  ! convert to deg C

        enddo
        enddo

  end subroutine theta_to_temp
!============================================================
!============================================================
