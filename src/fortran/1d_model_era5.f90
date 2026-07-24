!implicit none
program main
  use NETCDF

  integer, parameter :: IK = selected_int_kind(18)
  integer(IK) :: c0, cc1, rate, cmax, ticks
  real(8) :: elapsed
  integer :: hrs, mins, secs


  real :: L2_rad, K_rad
  real :: rad1, rad2, rad3
  real :: rhs1, rhs2, rhs3
  real :: dInvTauDz
  real :: dDelUDz1, dDelUDz2, dDelUDz3
  real :: delu_kp, delu_km
  real :: delu2_kp, delu2_km
  real :: delu1_kp, delu1_km
  integer :: k40


!  **** 1D FAWA equation in the vertical ****
  !integer,parameter :: years=42,ry=42, kmax = 501,myear = 1190400, mmax = myear*ry+1,mout = 372*10, msize = 3200*ry/10, mera = myear/992 !!!!!!!!!! dt = 18
  !integer,parameter :: years=42,ry=42, kmax = 501,myear = 535680, mmax = myear*ry+1,mout = 1674, msize = 13440, mera = myear/992  !!!!!!!!!! dt = 40
  !integer,parameter :: years=42,ry=42, kmax = 501,myear = 476160 , mmax = myear*ry+1,mout = 1488, msize = 13440, mera = myear/992  !!!!!!!!!! dt = 45
  !integer,parameter :: years=42,ry=42, kmax = 501,myear = 238080 , mmax = myear*ry+1,mout = 744, msize = 13440, mera = myear/992  !!!!!!!!!! dt = 90
  integer,parameter :: years=42,ry=42, kmax = 501,myear = 178560 , mmax = myear*ry+1,mout = 558, msize = 13440, mera = myear/992  !!!!!!!!!! dt = 120
  dimension p(kmax,4),fact(kmax),flux(kmax,3)
  dimension ucos(kmax,4),ucos0(kmax),tau(kmax),damp(kmax),inv_tau(kmax),tau_a(kmax),damp_a(kmax)
  dimension uref(kmax,4),uref0(kmax),cgzz(2,kmax,4),tusu(2,kmax,4),fawa_init(years,kmax),ubar_init(years,kmax),uref_init(years,kmax)
  common /array/ a(msize,kmax),u(msize,kmax),ur(msize,kmax),tau3(kmax),fz(kmax,4),fzn(kmax,4),cgz(msize,kmax,2)
  common /brray/ epz12(992,years),fawa_interp(992,years),para(6,kmax),urad(992,years,kmax),fzo(msize,kmax),fznn(msize,kmax)

  dimension aterm1_out(msize,kmax), aterm3_out(msize,kmax), aterm4_out(msize,kmax), atot_out(msize,kmax)
  dimension uterm1_out(msize,kmax), uterm2_out(msize,kmax), uterm4_out(msize,kmax), utot_out(msize,kmax)
  dimension dadt_out(msize,kmax), dudt_out(msize,kmax)

  real :: pt11, pt12, pt13, pt14
  real :: ut1, ut2, ut3, ut4
  real :: aterm1, aterm3, aterm4, atot
  real :: uterm1, uterm2, uterm4, utot


  integer :: vid_a1, vid_a3, vid_a4, vid_atot, vid_dadt
  integer :: vid_u1, vid_u2, vid_u4, vid_utot, vid_dudt

  integer :: ncid, status,nDim,nVar,nAtt,uDimID,inq
  integer :: lonID,latID,vid2,varID
  integer :: tt,aa,yy
  real :: f,bt,qy,lwn,m1s,m2s,ssw,rr
  real :: m1,m2,c1,c2,ud,coe
  character(len=5)  :: charI
  character(len=128) :: fname

  ! ==== Command-line parameters ====
  ! Added: damp1_days, u_thres, flux_extra
  real :: beta, s_tau, damp1_days
  real :: u_thres, flux_extra
  real :: alpha, const1, ppow, aexp
  
  character(len=128) :: outtag
  character(len=64)  :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9
  character(len=256) :: arg10
  character(len=16) :: bstr, sstr, d1str, uthstr, fxstr
  character(len=16) :: astr, c1str, pstr, aestr

  integer :: arglen, i
  logical :: has_tag


  ! ==== Other variables (10 hPa index) ====
  real :: pi, dz, dt, h, shear, rkappa, const2, gamma, L, c
  real :: ontime, offtime, dtime, dp, tau2, eps, phi, omega, bvf, radius, p0, p10, z10, bt_loc
  integer :: iout, m, k, k10
  real :: xi, uref_eff, shape_fac
  real, parameter :: tiny_u = 1.0e-6


  ! --- sponge layer profile & params ---
  real, parameter :: sponge_frac=0.75, sponge_pow=2.0
  real, parameter :: r_s_u_max=1.0/(24*3600.)   ! top-level u sponge: 1-day e-fold
  real, parameter :: r_s_p_max=1.0/(6*3600.)       ! top-level p sponge: 1-day e-fold
  real, parameter :: kappa_u_amp=5.0, kappa_p_amp=1.5
  dimension w_sponge(kmax), r_s_u(kmax), r_s_p(kmax), wloc(kmax)

  real :: tau_relax
  

  real, parameter :: umin = 5.0         ! empirical min wind near tropopause in winter
  real, parameter :: tau_b = 6.*3600.   ! lower-boundary temporal smoothing 6 h
  real :: p1_new

  real, parameter :: pneg = -5.0     ! allow small negative wave activity density
  real, parameter :: fneg = -2.0     ! allow small negative flux



  call system_clock(count_rate=rate, count_max=cmax)
  call system_clock(count=c0)

  ! ==== Default parameters (used if no CLI args) ====
  beta       = 0.233  !0.238
  s_tau      = 0.5  !0.509
  damp1_days = 5 !5 !4 !5.0 !3
  u_thres    = 7 !9.0 !7
  flux_extra = 16.0 !13 
  alpha      = 0.5 !0.4
  const1     = 2.2e-3
  ppow       = 4  ! 2.5
  aexp       = 0.3
  outtag     = ''
  has_tag    = .false.
  
! Read command-line arguments:
! 1 beta, 2 s_tau, 3 damp1_days, 4 u_thres, 5 flux_extra,
! 6 alpha, 7 const1, 8 ppow, 9 aexp, 10 tag (optional)

  call get_command_argument(1, arg1, length=arglen)
  if (arglen > 0) read(arg1, *) beta
  
  call get_command_argument(2, arg2, length=arglen)
  if (arglen > 0) read(arg2, *) s_tau
  
  call get_command_argument(3, arg3, length=arglen)
  if (arglen > 0) read(arg3, *) damp1_days
  
  call get_command_argument(4, arg4, length=arglen)
  if (arglen > 0) read(arg4, *) u_thres
  
  call get_command_argument(5, arg5, length=arglen)
  if (arglen > 0) read(arg5, *) flux_extra
  
  call get_command_argument(6, arg6, length=arglen)
  if (arglen > 0) read(arg6, *) alpha
  
  call get_command_argument(7, arg7, length=arglen)
  if (arglen > 0) read(arg7, *) const1
  
  call get_command_argument(8, arg8, length=arglen)
  if (arglen > 0) read(arg8, *) ppow
  
  call get_command_argument(9, arg9, length=arglen)
  if (arglen > 0) read(arg9, *) aexp
  
  call get_command_argument(10, arg10, length=arglen)
  if (arglen > 0) then
     outtag = trim(arg10)
     has_tag = .true.
  end if
  
  ! =====  tag“”b0p25s1p4_d5u9f14 =====
  ! Requires these declarations:
  !   character(len=16) :: bstr, sstr
  !   character(len=8)  :: d1str, uthstr, fxstr
  !   integer :: i
  
  if (.not. has_tag) then
    write(bstr  ,'(F4.2)') beta
    write(sstr  ,'(F3.1)') s_tau
    write(d1str ,'(I0)') nint(damp1_days)

    write(uthstr,'(F4.1)') u_thres
    write(fxstr ,'(F4.1)') flux_extra
    

    write(astr ,'(F4.2)') alpha
    write(c1str,'(ES10.3)') const1
    write(pstr ,'(F4.1)') ppow
    write(aestr,'(F4.3)') aexp



    bstr   = adjustl(bstr)
    sstr   = adjustl(sstr)
    d1str  = adjustl(d1str)
    uthstr = adjustl(uthstr)
    fxstr  = adjustl(fxstr)
    astr  = adjustl(astr)
    c1str = adjustl(c1str)
    pstr  = adjustl(pstr)
    aestr = adjustl(aestr)
 
    do i=1,len_trim(bstr)
      if (bstr(i:i)=='.') bstr(i:i)='p'
    end do
    do i=1,len_trim(sstr)
      if (sstr(i:i)=='.') sstr(i:i)='p'
    end do
    do i=1,len_trim(uthstr)
      if (uthstr(i:i)=='.') uthstr(i:i)='p'
    end do
    do i=1,len_trim(fxstr)
      if (fxstr(i:i)=='.') fxstr(i:i)='p'
    end do
    do i=1,len_trim(astr)
      if (astr(i:i)=='.') astr(i:i)='p'
    end do
    do i=1,len_trim(c1str)
      if (c1str(i:i)=='.') c1str(i:i)='p'
      if (c1str(i:i)=='+') c1str(i:i)='_'
    end do
    do i=1,len_trim(pstr)
      if (pstr(i:i)=='.') pstr(i:i)='p'
    end do
    do i=1,len_trim(aestr)
      if (aestr(i:i)=='.') aestr(i:i)='p'
    end do

    outtag = 'b'//trim(bstr)//'s'//trim(sstr)// &
            '_d'//trim(d1str)// &
            '_u'//trim(uthstr)//'f'//trim(fxstr)// &
            '_a'//trim(astr)//'_c'//trim(c1str)// &
            '_p'//trim(pstr)//'_e'//trim(aestr)
 end if
  
  fname = '_'//trim(outtag)   !  fname


  

  ! ======  ======
  status = nf90_open('/nas/winds-home/smliu01/0.1dmodel_2025/real_forcing.nc',  nf90_nowrite, ncid)
!  status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/real_w1w2.nc',  nf90_nowrite, ncid)
  status = nf90_inq_varid(ncid,"epz_32",varID)
  status = nf90_get_var(ncid,varID,epz12)
  status = nf90_close(ncid)

  status = nf90_open('/nas/winds-home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1dmodel_2025/fawa_interp.nc',  nf90_nowrite, ncid)
!  status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/real_w1w2.nc',  nf90_nowrite, ncid)
  status = nf90_inq_varid(ncid,"fawa_interp",varID)
  status = nf90_get_var(ncid,varID,fawa_interp)
  status = nf90_close(ncid)

  ! status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/para.nc',  nf90_nowrite, ncid)
  ! status = nf90_inq_varid(ncid,"para",varID)
  ! status = nf90_get_var(ncid,varID,para)
  ! status = nf90_close(ncid)

  !status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/graybox_mixing/urad_501.nc',  nf90_nowrite, ncid)
  !status = nf90_open('/mnt/winds/home/smliu01/0.1dmodel_2025/urad/urad_hm_notaper_501_42.nc',  nf90_nowrite, ncid)
  !status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/urad.nc',  nf90_nowrite, ncid)
  !status = nf90_open('/mnt/winds/home/smliu01/0.1dmodel_2025/urad/urad_piecewise_asym_piecewiseLinear_501_42.nc',  nf90_nowrite, ncid)
  !status = nf90_open('/mnt/winds/home/smliu01/0.1dmodel_2025/urad/urad_fit110_plus_HMshear_501_42.nc',  nf90_nowrite, ncid)
  status = nf90_open('/nas/winds-home/smliu01/0.1dmodel_2025/urad/urad_below32_fit_above32_HMshear_501_42.nc',  nf90_nowrite, ncid)
  status = nf90_inq_varid(ncid,"urad",varID)
  status = nf90_get_var(ncid,varID,urad)
  status = nf90_close(ncid)

  status = nf90_open('/nas/winds-home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/graybox_mixing_2026/init_501.nc',  nf90_nowrite, ncid)
  !status = nf90_open('/mnt/winds/home/smliu01/1.reproduce_nfl20_10.2020-1.2021/1d_refine/1d_para/init.nc',  nf90_nowrite, ncid)
  status = nf90_inq_varid(ncid,"fawa",varID)
  status = nf90_get_var(ncid,varID,fawa_init)
  status = nf90_inq_varid(ncid,"ubar",varID)
  status = nf90_get_var(ncid,varID,ubar_init)
  status = nf90_inq_varid(ncid,"uref",varID)
  status = nf90_get_var(ncid,varID,uref_init)
  status = nf90_close(ncid)

  urad = urad 

  ! ======  ======
  pi = acos(-1.)
  dz = int(100000/(kmax-1))
  dt = 60*2 !120 !90 !45 !18 !18.
  h  = 7000.
  !alpha = 0.5       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  shear = 0.0013
  rkappa = 30.!*(dz/100)**2
  !const1 = 2.2e-3 !2.23e-3    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  const2 = 4.11e-3
  gamma = 3.
  L = 1.e7
  c = 5.

  ! damp1 = 5
  ! damp1 = 1./damp1
  ! damp2 = 30.
  ! damp2 = 1./damp2

  ontime = 3600.*48.
  offtime = 3600.*168.
  dtime = 3600.*6.
  dp = 3.
  tau2 = 2.*3600.*24.

  eps = 8./3./pi

  p(:,:) = 0.
  ud = 0.
  coe = 1.2
  rr = 1.

  delta = 1.
  phi = pi/3.
  omega = 7.29e-5
  f = 2*omega*sin(phi)
  bvf = 2.1e-2
  radius = 6.378e6
  bt = 2*omega*cos(phi)/radius
  lwn = 1.73e-7

  ! ===  10 hPa  SSW ===
  p0  = 1000.0
  p10 = 10.0
  z10 = h * log(p0/p10)
  k10 = nint( (z10 - 10000.) / dz ) + 1
  if (k10 < 2) k10 = 2
  if (k10 > kmax-1) k10 = kmax-1

  iout = 1
  ssw = 0.

  !! sponge layer
  k_sp = kmax - int(sponge_frac*real(kmax-1))  ! 441
  do k=1,kmax
    if (k>k_sp) then
      w = real(k-k_sp)/real(kmax-k_sp)
      w = w**sponge_pow
      r_s_u(k)   = r_s_u_max * w
      r_s_p(k)   = r_s_p_max * w
      !write(*,*) w
    end if
    !  u/p 
  end do


  sigma_k = 36
  do k = 1,kmax
    wloc(k) = exp( -0.5 * ((real(k)-real(110))/sigma_k)**2 )
  end do



! ====== precompute tau(z) and inv_tau(z) ======
  damp1 = 1.0 / damp1_days
  damp2 = 1.0 / 30.0
  k40   = int((kmax-1)*0.4) + 1
  k70   = int((kmax-1)*0.7) + 1
  
  do k = 1, kmax
    if (k <= k40) then
      damp(k) = damp2 + real(k-1) * (damp1-damp2) / real(k40-1)
    else
      damp(k) = damp1
    endif
  
    tau(k)     = 3600.0 * 24.0 / damp(k)
    inv_tau(k) = 1.0 / tau(k)

    ! damp1_a = 9/20 !1.0/5.0
    ! damp2_a = 1.0/30.0


    ! damp_a(k) = (k-1) / 1200.0 + damp2_a


    ! z = dz*float(k-1) + 10000.
    ! zkm = z / 1000.0

    ! ! 32 km 
    ! damp32 = 22/110*(kmax-1) / 1200.0 + damp2_a
    ! !  0.125 day^-1

    ! z0      = 32.0     ! km
    ! Hdamp   = 1.0      ! km 6–10
    ! damp_top = 1.0/1.0 ! day^-1 damping 3 

    ! if (zkm <= z0) then
    !     damp_a(k) = (k-1) / 1200.0 + damp2_a
    ! else
    !     damp_a(k) = damp32 !* exp((zkm - z0)/Hdamp)
    !     ! damp_a(k) = min(damp_a(k), damp_top)
    ! endif

    ! tau_a(k) = 3600.0 * 24.0 / damp_a(k)



    ! z = dz*float(k-1) + 10000.
    ! zkm = z / 1000.0
    
    ! damp_lin = (k-1)/1200.0 + damp2_a
    
    ! z0    = 32.0
    ! Hdamp = 8.0
    ! Aexp  = 0.02
    ! w0    = 2.0
    
    ! gate = 0.5 * (1.0 + tanh((zkm - z0)/w0))
    
    ! damp_a(k) = damp_lin + gate * Aexp * (exp((zkm - z0)/Hdamp) - 1.0)
    ! damp_a(k) = min(damp_a(k), 1.0/2.5)
    
    ! tau_a(k) = 3600.0 * 24.0 / damp_a(k)



    ! z = dz*float(k-1) + 10000.
    ! zkm = z / 1000.0
    
    ! damp_lin = (k-1)/1200.0 + damp2_a
    
    ! z0    = 32.0
    ! w0    = 5.0
    ! Hdamp = 12.0
    ! Aexp  = 0.18
    
    ! gate = 0.5 * (1.0 + tanh((zkm - z0)/w0))
    
    ! damp_raw = damp_lin + gate * Aexp * (1.0 - exp(-(zkm - z0)/Hdamp))
    
    ! damp_cap = 1.0/2.5
    
    ! damp_a(k) = damp_cap * tanh(damp_raw / damp_cap)
    
    ! tau_a(k) = 3600.0 * 24.0 / damp_a(k)




    z = dz*float(k-1) + 10000.0
    zkm = z / 1000.0
    damp2_a = 1.0 / 30.0
    
    zbot = 10.0
    z0   = 32.0
    
    ! lower part: convex below 32 km
    !ppow = 3.0              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! upper part: concave / saturating above 32 km
    Hdamp = 12.0
    !Aexp  = 0.25 !0.25 !0.18      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    w0    = 5.0
    ! smooth cap
    damp_cap = 1.0 / 2.5
    ! value at 32 km, matched to old linear profile
    damp32 = ((z0 - zbot) / (dz/1000.0)) / 1200.0 + damp2_a
    ! lower branch: convex
    xlow = (zkm - zbot) / (z0 - zbot)
    xlow = max(0.0, min(1.0, xlow))
    damp_low = damp2_a + (damp32 - damp2_a) * xlow**ppow
    
    ! upper branch: concave / saturating
    damp_high = damp32 + aexp * (1.0 - exp(-(zkm - z0) / Hdamp))
    ! smooth transition around 32 km
    gate = 0.5 * (1.0 + tanh((zkm - z0) / w0))
    ! blend low/high
    damp_raw = (1.0 - gate) * damp_low + gate * damp_high
    ! smooth cap
    damp_a(k) = damp_cap * tanh(damp_raw / damp_cap)
    
    tau_a(k) = 3600.0 * 24.0 / damp_a(k)




    ! z = dz*float(k-1)+10000.
    ! zkm = z / 1000.0
    
    ! if (zkm <= 32.0) then
    !     damp_a(k) = 1.0/30.0
    ! else if (zkm <= 40.0) then
    !     damp_a(k) = 1.0/30.0 + (zkm-32.0)*(1.0/12.0 - 1.0/30.0)/(40.0-32.0)
    ! else if (zkm <= 50.0) then
    !     damp_a(k) = 1.0/12.0 + (zkm-40.0)*(1.0/6.0 - 1.0/12.0)/(50.0-40.0)
    ! else
    !     damp_a(k) = 1.0/6.0
    ! endif
    
    ! tau_a(k) = 3600.0 * 24.0 / damp_a(k)

    !inv_tau_a(k) = 1.0 / tau_a(k)


  enddo
  
  inv_tau(1)    = inv_tau(2)
  inv_tau(kmax) = inv_tau(kmax-1)
  
  ! Noboru note: use L^2 = 2.1e12 m^2
  L2_rad = 2.1e12
  K_rad  = f*f*L2_rad / (bvf*bvf)


!   **** main loop starts here ****
  do m = 1,mmax
    t = dt*float(m-1)
    aa = int((m-1)/mera)+1
    tt = mod(aa,992)+1
    yy = int(aa/992)+1

    if (tt .lt. 2) ssw = 0.

    if (mod(m,myear).eq.1) then
      do k = 1,kmax
        z = dz*float(k-1)+10000.
        fact(k) = exp(-z/h)
        if(k.le. int((kmax-1)*0.4)+1) then
          !damp(k) = damp2 + float(k-1)*(damp1-damp2)/int((kmax-1)*0.4)
          ucos(k,:) =  ubar_init(yy,k)
          uref(k,:) =  uref_init(yy,k)
          ucos0(k)  =  ubar_init(yy,k)
          !tau(k) = 3600.*24./damp(k)
          p(k,:) = 0.
          !p(2,k,:) = 0.
        else
          !damp(k) = damp1
          !tau(k) = 3600.*24./damp1
          ucos(k,:) =  ubar_init(yy,k)
          uref(k,:) =  uref_init(yy,k)
          ucos0(k)  =  ubar_init(yy,k)
          p(k,:) = 0.
          !p(2,k,:) = 0.
        endif
      enddo
    end if

    flux(1,3) = rr*epz12(tt,yy)
    ! flux(2,1,3) = rr*epz12(tt,yy,2)
    p(1,3) = flux(1,3)/const1/(ucos(1,3)) !fawa_interp(tt,yy)*fact(1) !flux(1,3)/const1/(urad(tt,yy,1))  !flux(1,3)/const1/(ucos(1,3))
    ! p(2,1,3) = flux(2,1,3)/const2/(ucos(1,3))
    if (p(1,3) .lt. 0) then
      p(1,3) = 0
    endif
    if (flux(1,3) .lt. 0) then
      flux(1,3) = 0
    endif
    ! if (p(2,1,3) .lt. 0) then
    !   p(2,1,3) = 0
    ! endif
    
    do k = 1,kmax
      if (k .eq. 1) then
        c1 = const1
        c2 = const2
      else if (k .eq. kmax) then
        c1 = 0.
        c2 = 0.
      else
        qy = bt+eps*(lwn**2*uref(k,3)     &
          - 1/fact(k)*f**2/bvf*( fact(k+1)*(uref(k+1,3)-uref(k,3))/dz- fact(k)*(uref(k,3)-uref(k-1,3))/dz )/dz  &
          - lwn**2*(p(k,3))/fact(k))
        m1s = bvf**2/f**2*(qy/ucos(k,3)-((1/radius)**2+lwn**2)-f**2/4/bvf**2/h**2)
        m2s = bvf**2/f**2*(qy/ucos(k,3)-((2/radius)**2+lwn**2)-f**2/4/bvf**2/h**2)
        m1 = SQRT(max(0.,m1s))
        m2 = SQRT(max(0.,m2s))
        c1 = (2*f**2/bvf**2)*(ucos(k,3)/cos(phi)/qy)*(1./radius)*m1
        c2 = (2*f**2/bvf**2)*(ucos(k,3)/cos(phi)/qy)*(2./radius)*m2
      end if

      ! 
      c1 = const1
      c2 = const2


      ! z = dz*float(k-1)+10000.
      ! zkm = z / 1000.0
      ! gcz = 0.5*(1.0 + tanh((zkm - 40.0)/4.0))
      ! c1 = const1 * (1.0 - 0.4*gcz)

      tusu(1,k,3) = m1s
      tusu(2,k,3) = m2s
      cgzz(1,k,3) = c1
      cgzz(2,k,3) = c2

      if (k .ne. 1) then
        if ( p(k,3) .lt. 0) then
          p(k,3) = 0.
          flux(k,3) = 0.
          ! p(k,3)    = max(pneg, p(k,3))
          ! flux(k,3) = max(fneg, flux(k,3))
        else

          wT = 1


          ! z1_c = 10.0
          ! z2_c = 30.0
          ! c1_bot_fac = 5.0
          ! c1_top_fac = 1.0

          ! z   = dz*float(k-1) + 10000.
          ! zkm = z / 1000.0
          
          ! if (zkm <= z1_c) then
          !   c1_eff = c1_bot_fac * const1
          ! else if (zkm >= z2_c) then
          !   c1_eff = c1_top_fac * const1
          ! else
          !   frac = (zkm - z1_c) / (z2_c - z1_c)
          !   c1_eff = (c1_bot_fac + (c1_top_fac - c1_bot_fac)*frac) * const1
          ! end if

          flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + flux_extra * (1-(1/(1+exp(-(ucos(k,3)-u_thres)/wT)))) ) * (p(k,3) -  fact(k) * flux_extra * (1-(1/(1+exp(-(ucos(k,3) - u_thres )/wT))))  )

          if (flux(k,3) .lt. 0) then
            flux(k,3) = 0.
          end if

          ! z = dz*float(k-1)+10000.
          ! zkm = z / 1000.0
          ! gbot = 1.0 - 0.30 * exp(-((zkm - 20.0)/7.0)**2)
          ! c1   = const1 * gbot * (1.0 - 0.4*gcz)


          !flux(k,3) = c1_eff * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) ) * p(k,3)

          ! uref_eff = uref(k,3) !sign(max(abs(uref(k,3)), tiny_u), uref(k,3))
          ! xi = alpha * p(k,3) / (fact(k) * uref_eff)
      
          ! shape_fac = 1.0 - (xi - sigma_ep)**n_ep
      
          ! flux(k,3) = muC * fact(k) * (uref_eff*uref_eff/alpha) * xi * shape_fac
      

          ! if (ucos(k,3) .lt.  0 .and.  k .eq. 110) then
          !   flux(k,3) = c1 * ( abs(uref(k,3)-alpha*(p(k,3)/fact(k) )) ) * p(k,3)
          ! end if

          !flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + delta_u * (1-(1/(1+exp(-(ucos(k,3)-u_thres)/wT)))) ) * (p(k,3) )
          !flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + delta_u * (1-(1/(1+exp(-(ucos(k,3)-u_thres)/wT)))) ) * (  p(k,3) - fact(k) * delta_p * (1-(1/(1+exp(-(ucos(k,3))/wT)))) ) !* (1-(1/(1+exp(-(urad(tt,yy,k)-15)/wT)))) )!* (1 / (1+exp(-(p(k,3)/fact(k)-20)/1))) )



          ! aaa = -0.2-0.2**2 !0.2 !0.18 !
          ! bbb = 1+0.2*2  ! -1 !! !1
          ! ccc = -1.  !!-1.  !3.66    !-1
          ! ddd = 0.    !!-2.58 

          ! term0 = c1 * aaa * (uref(k,3)**2/alpha) * fact(k)
          ! term1 = c1 * bbb * uref(k,3) * p(k,3)
          ! term2 = c1 * ccc * (alpha/fact(k)) * (p(k,3)**2)
          ! !  u  0  epsilon
          ! eps  = 1.0d-8
          ! ueff = sign(max(abs(uref(k,3)), eps), uref(k,3))   ! 0
          ! term3 = c1 * ddd * (alpha*alpha/(fact(k)*fact(k))) * (p(k,3)**3) / ueff

          ! flux(k,3) = term0 + term1 + term2 + term3
          ! flux(k,3) = c1 * uref(k,3) ** 2 / alpha * fact(k) * (aaa + bbb*(alpha*p(k,3)/fact(k)/uref(k,3))+ ccc*(alpha*p(k,3)/fact(k)/uref(k,3))**2 + ddd*(alpha*p(k,3)/fact(k)/uref(k,3))**3)



          !coeff = 0.3
          
          ! flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + coeff*uref(k,3) ) * (p(k,3) - coeff*alpha*uref(k,3)*fact(k))
          ! !flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) ) * p(k,3)

          ! if (ucos(k,3) .lt.  u_thres) then
          !   flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + flux_extra ) * p(k,3)
          ! else
          !   flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) ) * p(k,3)
          ! end if

          ! u_thres = 10.0
          ! wT = 2.0
          ! flux_extra = 14.0
          ! u_here = uref(k,3) - alpha*(p(k,3)/fact(k))   ! ≈ u

          ! u_gate = ucos(k,3)                           !  u uref - α·p/fact
          ! arg = (u_thres - u_gate)/wT
          ! arg = max(-20.0, min(20.0, arg))             ! ←  [-20,20]tanh 
          ! H   = 0.5*(1.0 + tanh(arg))
          ! H   = max(0.0, min(1.0, H))                  ! ←  H  [0,1]

          ! ! write(*,*) u_here, H, flux_extra*H

          ! !  bb
          ! flux(k,3) = c1 * ( u_here + flux_extra * H ) * p(k,3)


          ! u0 = 10.0
          ! w = 3.0
          ! c_small = 0.0
          ! c_big = 20
          ! Hpos = 0.5 * ( 1.0 + tanh( ( ucos(k,3) - 0.0 ) / 2.0 ) ) 
          ! S = 1.0 - 0.5*( 1.0 + tanh( (abs(ucos(k,3)) - u0)/w ) ) 
          ! c_eff = (c_small + (c_big - c_small)*S) * Hpos
          ! flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + c_eff) * p(k,3)


          ! u_opt = 10.0
          ! w_opt = 6.0
          ! c_peak = 14.0
          ! U_sat = 12.0
          ! u_here = uref(k,3) - alpha*(p(k,3)/fact(k))   ! = ucos(k,3)
          ! Hpos  = 0.5*(1.0 + tanh( u_here/2.0 ))             ! u>0≈1, u<0≈0
          ! Gbell = exp( -0.5*((u_here - u_opt)/w_opt)**2 )    ! 
          ! c_eff = c_peak * Gbell * Hpos
          ! u_star = (u_here + c_eff)
          ! u_star = U_sat * tanh( u_star / U_sat )
          ! flux(k,3) = c1 * u_star * p(k,3)



          !flux(k,3) = c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) )) * p(k,3)
        end if

        ! if ( p(2,k,3) .lt. 0) then
        !   p(2,k,3) = 0.
        !   flux(2,k,3) = 0.
        ! else
        !   flux(2,k,3) = c2 * ( uref(k,3)-alpha*(p(2,k,3)/fact(k) ) ) * p(2,k,3)
        ! end if


        ! if (flux(2,k,3) .lt. 0) then
        !   flux(2,k,3) = 0.
        ! end if
      end if

      if(m.eq.1) then
        flux(k,2) = flux(k,3)
        flux(k,1) = flux(k,2)
        p(k,2) = p(k,3)
        p(k,1) = p(k,2)
        ucos(k,2) = ucos(k,3)
        ucos(k,1) = ucos(k,2)
        uref(k,2) = uref(k,3)
        uref(k,1) = uref(k,2)
      endif
    enddo

    flux(kmax,1) = flux(kmax-1,1)
    flux(kmax,2) = flux(kmax-1,2)
    flux(kmax,3) = flux(kmax-1,3)

    if(mod(m,mout).eq.1) then
      if(ucos(k10,3) .lt. 0.) then
        ssw = ssw + 1.
      endif
    endif

    do k = 2,kmax-1

      rkappa_k = rkappa !* kappa_fac(k)

      z = dz*float(k-1)+10000.

      ! damp1 = damp1_days
      ! damp1 = 1./damp1
      ! damp2 = 30.
      ! damp2 = 1./damp2

      ! if(k.le.int((kmax-1)*0.4)+1) then
      !   damp(k) = damp2 + float(k-1)*(damp1-damp2)/int((kmax-1)*0.4)
      !   tau(k) = 3600.*24./damp(k)
      ! else
      !   damp(k) = damp1
      !   tau(k) = 3600.*24./damp1
      ! endif

      pt11 = - dt*(23.*(flux(k+1,3)-flux(k-1,3)) - 16.*(flux(k+1,2)-flux(k-1,2)) + 5.*(flux(k+1,1)-flux(k-1,1)))/(12.*2.*dz)
      pt12 = 0.

      ! ===  beta  s_tau ===
      eps_r  = 1.0e-8
      ratio  = ucos(k,3) / max(abs(uref(k,3)), eps_r)
      tau3_nom = tau_a(k) * (1.0 - beta + beta * (ratio*ratio))
      tau3_adj = s_tau * tau3_nom
      pt13 =  - dt*(23.0* (p(k,3) -fact(k)*0  )  - 16.0*(p(k,2)-fact(k)*0 ) + 5.0*(p(k,1)-fact(k)*0 )) / (12.0*tau3_adj)
      !pt13 = - dt*(23.0*p(k,3) - 16.0*p(k,2) + 5.0*p(k,1)) / (12.0*tau(k))
      pt14 = dt*rkappa*(23.*p(k+1,3)-16.*p(k+1,2)+5.*p(k+1,1) + 23.*p(k-1,3)-16.*p(k-1,2)+5.*p(k-1,1) - 2.*(23.*p(k,3)-16.*p(k,2)+5.*p(k,1)))/(12.*dz*dz)
 
      !! 
      if ( (p(k,3))/fact(k) .gt. 170) then
        pt13 =   - dt*(23.*p(k,3)-16.*p(k,2)+5.*p(k,1))/(12.*tau3_adj)*200
      endif

      !p(k,4) =  p(k,3) + pt11 + pt12 + pt13 + pt14

      !rkappa_p = rkappa * (1.0 + kappa_p_amp * w_sponge(k))
      !pt14 = dt*rkappa_p*(23.*p(k+1,3)-16.*p(k+1,2)+5.*p(k+1,1) + 23.*p(k-1,3)-16.*p(k-1,2)+5.*p(k-1,1) - 2.*(23.*p(k,3)-16.*p(k,2)+5.*p(k,1)))/(12.*dz*dz)
      p(k,4)    = p(k,3)    +  pt11 + pt12 + pt13 + pt14 !- dt*r_s_p(k)*p(k,3)
      ! p(k,4)    = p(k,3)    +  2*pt11 + pt12 + pt13 + pt14 !- dt*r_s_p(k)*p(k,3)


      ! pt21 = - dt*(23.*(flux(2,k+1,3)-flux(2,k-1,3)) - 16.*(flux(2,k+1,2)-flux(2,k-1,2)) + 5.*(flux(2,k+1,1)-flux(2,k-1,1)))/(12.*2.*dz)
      ! pt22 = 0.

      ! ratio  = ucos(k,3) / max(abs(uref(k,3)), eps_r)
      ! tau3_nom = tau(k) * (1.0 - beta + beta * (ratio*ratio))
      ! tau3_adj = s_tau * tau3_nom
      ! pt23 = - dt*(23.0*p(2,k,3) - 16.0*p(2,k,2) + 5.0*p(2,k,1)) / (12.0*tau3_adj)
      ! pt24 = dt*rkappa*(23.*p(2,k+1,3)-16.*p(2,k+1,2)+5.*p(2,k+1,1) + 23.*p(2,k-1,3)-16.*p(2,k-1,2)+5.*p(2,k-1,1) - 2.*(23.*p(2,k,3)-16.*p(2,k,2)+5.*p(2,k,1)))/(12.*dz*dz)

      ! if ( (p(1,k,3) + p(2,k,3))/fact(k) .gt. 170) then
      !   pt23 =   - dt*(23.*p(2,k,3)-16.*p(2,k,2)+5.*p(2,k,1))/(12.*tau3_adj)*20
      ! endif



      ! p(2,k,4) = p(2,k,3) + pt21 + pt22 + pt23 + pt24 !- dt * r_s_p(k) * p(2,k,3)

      A_here = p(k,3)/fact(k)
      A_star = 150
      p_pow = 2.0
      alpha_eff = 0.4 !0.4 / (1.0D0 + (A_here/A_star)**p_pow) 

      ut1 = -(alpha_eff/fact(k)) * ( - dt*(23.*(flux(k+1,3)-flux(k-1,3)) - 16.*(flux(k+1,2)-flux(k-1,2)) + 5.*(flux(k+1,1)-flux(k-1,1)))/(12.*2.*dz) )


      !   ! 
      !  k10 σ≈30–40  ≈ 6–8 km
      gate_neg = 0.5  - 0.5*tanh( (ucos(k,3) )/2 ) !0.5*(1.0 - tanh( (ucos(k,3) + 5.0)/10.0 ))    !   ! u<0→≈1, u>0→≈0

      ! × (1 + A*wloc)
      tau_relax = tau(k) / (1.0 + 2*wloc(k)*gate_neg)    !  2.0  1–3
      !urad(tt,yy,k) = urad(tt,yy,k) - 0.005*p(k,3)/fact(k)
      !ut2 = -dt*(23.*( ucos(k,3)-urad(tt,yy,k) ) - 16.*(ucos(k,2)-urad(tt,yy,k)) + 5.*(ucos(k,1)-urad(tt,yy,k))) /  (12 * tau_relax)
      ! ! !  


      ! ! 
      ! u_z = 6.0        ! 5–7
      ! sigma_k = 30.0   ! ~6 km
      ! k10 = 110
      
      ! ! 
      ! boost_pos = -0.7   ! u>0 “”(=)-1.0-0.5
      ! boost_neg =  2.0   ! u<0 “”1.5–2.5
      
      ! wloc      = exp(-0.5*((real(k)-real(k10))/sigma_k)**2)
      ! gate_zero = exp(-0.5*(ucos(k,3)/u_z)**2)          ! |u| 
      ! Hpos      = 0.5*(1.0 + tanh(ucos(k,3)/2.0))       ! u>0
      ! Hneg      = 1.0 - Hpos                         ! u<0
      
      ! ! u>0 →  (tau × [1 + |boost_pos| * wloc*gate_zero])
      ! !         u<0 →  (tau / [1 + boost_neg  * wloc*gate_zero])
      ! tau_zero_pos = tau(k) * (1.0 + abs(boost_pos)*wloc*gate_zero)
      ! tau_zero_neg = tau(k) / (1.0 +      boost_neg *wloc*gate_zero)
      
      ! tau_relax = Hpos * tau_zero_pos + Hneg * tau_zero_neg
      
      ! ut2 = -dt*(23.*( ucos(k,3)-urad(tt,yy,k) ) - 16.*( ucos(k,2)-urad(tt,yy,k) ) + 5.*( ucos(k,1)-urad(tt,yy,k) )) / (12.*tau_relax)
      


      ! ut2 =  -dt*( 23.*(ucos(k,3)-urad(tt,yy,k)) - 16.*(ucos(k,2)-urad(tt,yy,k)) + 5.*(ucos(k,1)-urad(tt,yy,k)) ) / (12.*tau_relax)

      ! inputs: urad0, A_here
      A0   = 70
      wA   = 8
      Ahalf= 30
      kA   = 1        ! >0  urad<0  urad

      S    = 1.0 / (1.0 + exp(-(A_here - A0)/wA))
      Aex  = max(0.0, A_here - A0)
      !urad(tt,yy,k) = urad(tt,yy,k) - kA * S * ( Aex / (Ahalf + Aex) )


      ! u_center  = 4.0d0      ! “”m/s -u_center
      ! u_width   = 5.0d0      ! m/s4–7 
      ! wT_gate   = 3.0d0      ! /
      ! cpush_day = 0.35d0     ! m/s/day0.2–0.6 
      ! cpush     = cpush_day / 86400.d0   !  dt “” dt “”
      
      ! ! ====== u<0 →≈1u>0 →≈0======
      ! gate_neg_only = 0.5d0 - 0.5d0 * tanh( ucos(k,3) / wT_gate )

      ! ! ======  u≈-u_center  ======
      ! ! |u + u_center| << u_width → 
      ! gate_band = exp( - ((ucos(k,3) + u_center)/u_width)**2 )

      ! ! ====== “” ======
      ! ut2_push = -dt * cpush * gate_neg_only * gate_band *100




      ! ==========================================================
      ! New radiative relaxation term from Noboru's Eq. (3):
      ! -(1-alpha) * [ (u-urad)/tau + K_rad * d(1/tau)/dz * d(u-urad)/dz ]
      ! Apply AB3 to the whole bracket
      ! ==========================================================

      dInvTauDz = (inv_tau(k+1) - inv_tau(k-1)) / (2.0*dz)

      ! ---- time level n (slot 3) ----
      delu_kp  = ucos(k+1,3) - urad(tt,yy,k+1)
      delu_km  = ucos(k-1,3) - urad(tt,yy,k-1)
      dDelUDz3 = (delu_kp - delu_km) / (2.0*dz)

      rhs3 = (ucos(k,3) - urad(tt,yy,k)) * inv_tau(k)  &
           + K_rad * dInvTauDz * dDelUDz3

      ! ---- time level n-1 (slot 2) ----
      delu2_kp = ucos(k+1,2) - urad(tt -1,yy,k+1)
      delu2_km = ucos(k-1,2) - urad(tt -1,yy,k-1)
      dDelUDz2 = (delu2_kp - delu2_km) / (2.0*dz)

      rhs2 = (ucos(k,2) - urad(tt -1,yy,k)) * inv_tau(k)  &
           + K_rad * dInvTauDz * dDelUDz2

      ! ---- time level n-2 (slot 1) ----
      delu1_kp = ucos(k+1,1) - urad(tt -2,yy,k+1)
      delu1_km = ucos(k-1,1) - urad(tt -2,yy,k-1)
      dDelUDz1 = (delu1_kp - delu1_km) / (2.0*dz)

      rhs1 = (ucos(k,1) - urad(tt -2,yy,k)) * inv_tau(k)  &
           + K_rad * dInvTauDz * dDelUDz1

      ! ut2 = -(1.0 - alpha) * dt * (23.0*rhs3 - 16.0*rhs2 + 5.0*rhs1) / 12.0




      ut2 = -dt*(23.*(ucos(k,3)-urad(tt,yy,k)) -16.*(ucos(k,2)-urad(tt-1,yy,k)) +5.*(ucos(k,1)-urad(tt-2,yy,k))) / (12. *tau(k)) ! + ut2_push





      ut3 = 0.
      ut4 = dt*rkappa*(23.*ucos(k+1,3) -16.*ucos(k+1,2)+ 5.*ucos(k+1,1) + 23.*ucos(k-1,3) -16.*ucos(k-1,2) +5.*ucos(k-1,1) - 2.*(23.*ucos(k,3)-16.*ucos(k,2)+ 5.*ucos(k,1)))/(12.*dz*dz)
      ! ucos(k,4) = ucos(k,3) + ut1 + ut2 + ut3 + ut4

      ! if (k == kmax-1) then
      !   u_p3 = ucos(k   ,3); u_m3 = ucos(k-1,3); u0_3 = ucos(k,3)
      !   u_p2 = ucos(k   ,2); u_m2 = ucos(k-1,2); u0_2 = ucos(k,2)
      !   u_p1 = ucos(k   ,1); u_m1 = ucos(k-1,1); u0_1 = ucos(k,1)
      ! else
      !   u_p3 = ucos(k+1,3); u_m3 = ucos(k-1,3); u0_3 = ucos(k,3)
      !   u_p2 = ucos(k+1,2); u_m2 = ucos(k-1,2); u0_2 = ucos(k,2)
      !   u_p1 = ucos(k+1,1); u_m1 = ucos(k-1,1); u0_1 = ucos(k,1)
      ! endif
      
      ! lap3 = (u_p3 - 2.0*u0_3 + u_m3)/dz2
      ! lap2 = (u_p2 - 2.0*u0_2 + u_m2)/dz2
      ! lap1 = (u_p1 - 2.0*u0_1 + u_m1)/dz2
      
      ! ut4  = dt * rkappa * (23.0*lap3 - 16.0*lap2 + 5.0*lap1)/12.0
      

      !rkappa_u = rkappa * (1.0 + kappa_u_amp * w_sponge(k))
      !ut4       = dt*rkappa_u*(23.*ucos(k+1,3) -16.*ucos(k+1,2)+ 5.*ucos(k+1,1) + 23.*ucos(k-1,3) -16.*ucos(k-1,2) +5.*ucos(k-1,1) - 2.*(23.*ucos(k,3)-16.*ucos(k,2)+ 5.*ucos(k,1)))/(12.*dz*dz)
      ucos(k,4) = ucos(k,3) + ut1 + ut2 + ut3 + ut4    !- dt*r_s_u(k)*ucos(k,3)

      ! tau_lp = 36.0*3600.0       ! 36 h e-fold; try 24–48 h
      ! lambda = min( 0.3, dt/tau_lp )   ! stability guard

      ! ucos(k,4) = (1.0 - lambda)*ucos(k,4) + lambda*ucos(k,3)
      ! uref(k,4) = (1.0 - lambda)*uref(k,4) + lambda*uref(k,3)

      
      uref(k,4) = ucos(k,4)+alpha*(p(k,4))/fact(k) + (ud)

      fz(k,3) =  flux(k,3) !c1 * ( uref(k,3)-alpha*(p(k,3)/fact(k) ) + flux_extra * (1-(1/(1+exp(-(ucos(k,3)-u_thres)/wT)))) ) * (p(k,3) -  fact(k) * flux_extra * (1-(1/(1+exp(-(ucos(k,3) - u_thres )/wT))))  ) !cgzz(1,k,3)*p(k,3)* ucos(k,3) !( uref(k,3)-alpha*(p(k,3)/fact(k) ) ) 

      fzn(k,3) = (cgzz(1,k,3)*p(k,3)*( ucos0(k)-alpha*(p(k,3)/fact(k) ) ) ) /fact(k)*alpha/uref(k,3)/uref(k,3)





      aterm1 = (pt11 / fact(k)) / dt
      aterm3 = (pt13 / fact(k)) / dt
      aterm4 = (pt14 / fact(k)) / dt
      atot   = ((pt11 + pt13 + pt14) / fact(k)) / dt

      uterm1 = ut1 / dt
      uterm2 = ut2 / dt
      uterm4 = ut4 / dt
      utot   = (ut1 + ut2 + ut4) / dt

      if (mod(m,mout).eq.0) then
        aterm1_out(iout,k) = aterm1
        aterm3_out(iout,k) = aterm3
        aterm4_out(iout,k) = aterm4
        atot_out(iout,k)   = atot

        uterm1_out(iout,k) = uterm1
        uterm2_out(iout,k) = uterm2
        uterm4_out(iout,k) = uterm4
        utot_out(iout,k)   = utot

        dadt_out(iout,k) = ((p(k,4)-p(k,3))/fact(k)) / dt
        dudt_out(iout,k) = (ucos(k,4)-ucos(k,3)) / dt
      endif


      
    enddo

    !write(*,*) m, ucos(2,4), p(1,2,4), p(2,2,4)


    aa = int(m/mera)+1
    tt = mod(aa,992)+1
    yy = int(aa/992)+1
    flux(1,3) = rr*epz12(tt,yy)
    ! flux(2,1,4) = rr*epz12(tt,yy,2)


    ! p1_new = flux(1,3) / (const1 * max(abs(ucos(1,4)), umin))
    ! if (p1_new < 0.) p1_new = 0.
    ! p(1,4) = p(1,4) + (dt/tau_b) * (p1_new - p(1,4))   ! 
    ! ucos(1,4) = ucos0(1) - alpha * p(1,4)/fact(1)


    p(1,4) = flux(1,3)/const1/ucos(1,4)!*2
    if (p(1,4) .lt. 0) then
      p(1,4) = 0
    end if
    ucos(1,4) = ucos0(1)-alpha*(p(1,4))/fact(1)



    ! ucos(kmax,4) = 2.*ucos(kmax-1,3)-ucos(kmax-2,2)
    ! p(:,kmax,4) = 2.*p(:,kmax-1,3)-p(:,kmax-2,2)
    

    ucos(kmax,4) = ucos(kmax-1,4)
    p(kmax,4)  = p(kmax-1,4)
    ! p(2,kmax,4)  = p(2,kmax-1,4)    

    uref(kmax,4) = ucos(kmax,4)+alpha*(p(kmax,4))/fact(kmax) + (ud)

    if(mod(m,mout).eq.0) then
      a(iout,:)     = (p(:,3))/fact(:)
      u(iout,:)     = ucos(:,3)
      ur(iout,:)    = uref(:,3)
      fzo(iout,:)   = fz(:,3)
      fznn(iout,:)  = fzn(:,3)
      cgz(iout,:,1) = cgzz(1,:,3)
      cgz(iout,:,2) = cgzz(2,:,3)
      iout = iout+1
    endif

    p(:,1) = p(:,2)
    p(:,2) = p(:,3)
    p(:,3) = p(:,4)
    ucos(:,1) = ucos(:,2)
    ucos(:,2) = ucos(:,3)
    ucos(:,3) = ucos(:,4)
    uref(:,1) = uref(:,2)
    uref(:,2) = uref(:,3)
    uref(:,3) = uref(:,4)
    flux(:,1) = flux(:,2)
    flux(:,2) = flux(:,3)

    if (mod(m,myear).eq.0) then
      write(6,*) '*********************************************************************'
      write(6,*) yy,p(1,3)/fact(1),ucos(110,3),uref(110,3)
      write(6,*) '*********************************************************************'
    endif
  enddo

  ! ======  tag ======
  status = nf90_create('aout'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  status = nf90_def_dim(ncid2,"time",msize,it)
  status = nf90_def_dim(ncid2,"height",kmax,iz)
  status = nf90_def_var(ncid2,"fawa",nf90_float,(/it,iz/), vid2)
  status = nf90_put_att(ncid2,vid2,"title",'aout_real_urad'//trim(fname)//'.nc')
  status = nf90_enddef(ncid2)
  status = nf90_put_var(ncid2,vid2,a)
  status = nf90_close(ncid2)

  status = nf90_open('aout'//trim(fname)//'.nc',nf90_nowrite,ncid)
  status = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
  status = nf90_inq_varid(ncid,"fawa",varID)
  status = nf90_get_var(ncid,varID,v)
  status = nf90_close(ncid)
  !write(6,*) a(160,200)

  status = nf90_create('uout'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  status = nf90_def_dim(ncid2,"time",msize,it)
  status = nf90_def_dim(ncid2,"height",kmax,iz)
  status = nf90_def_var(ncid2,"ubar",nf90_float,(/it,iz/), vid2)
  status = nf90_put_att(ncid2,vid2,"title",'uout_real_urad'//trim(fname)//'.nc')
  status = nf90_enddef(ncid2)
  status = nf90_put_var(ncid2,vid2,u)
  status = nf90_close(ncid2)

  status = nf90_open('uout'//trim(fname)//'.nc',nf90_nowrite,ncid)
  status = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
  status = nf90_inq_varid(ncid,"ubar",varID)
  status = nf90_get_var(ncid,varID,v)
  status = nf90_close(ncid)
  !write(6,*) u(160,200)

  status = nf90_create('urout'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  status = nf90_def_dim(ncid2,"time",msize,it)
  status = nf90_def_dim(ncid2,"height",kmax,iz)
  status = nf90_def_var(ncid2,"uref",nf90_float,(/it,iz/), vid2)
  status = nf90_put_att(ncid2,vid2,"title",'urout_real_urad'//trim(fname)//'.nc')
  status = nf90_enddef(ncid2)
  status = nf90_put_var(ncid2,vid2,ur)
  status = nf90_close(ncid2)

  status = nf90_open('urout'//trim(fname)//'.nc',nf90_nowrite,ncid)
  status = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
  status = nf90_inq_varid(ncid,"uref",varID)
  status = nf90_get_var(ncid,varID,v)
  status = nf90_close(ncid)
  !write(6,*) ur(160,200)

  status = nf90_create('fzout'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  status = nf90_def_dim(ncid2,"time",msize,it)
  status = nf90_def_dim(ncid2,"height",kmax,iz)
  status = nf90_def_var(ncid2,"fz",nf90_float,(/it,iz/), vid2)
  status = nf90_put_att(ncid2,vid2,"title",'fzout_real_urad'//trim(fname)//'.nc')
  status = nf90_enddef(ncid2)
  status = nf90_put_var(ncid2,vid2,fzo)
  status = nf90_close(ncid2)

  status = nf90_open('fzout'//trim(fname)//'.nc',nf90_nowrite,ncid)
  status = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
  status = nf90_inq_varid(ncid,"fz",varID)
  status = nf90_get_var(ncid,varID,v)
  status = nf90_close(ncid)

  status = nf90_create('fznout'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  status = nf90_def_dim(ncid2,"time",msize,it)
  status = nf90_def_dim(ncid2,"height",kmax,iz)
  status = nf90_def_var(ncid2,"fzn",nf90_float,(/it,iz/), vid2)
  status = nf90_put_att(ncid2,vid2,"title",'fznout_real_urad'//trim(fname)//'.nc')
  status = nf90_enddef(ncid2)
  status = nf90_put_var(ncid2,vid2,fznn)
  status = nf90_close(ncid2)

  status = nf90_open('fznout'//trim(fname)//'.nc',nf90_nowrite,ncid)
  status = nf90_inquire(ncid,nDim,nVar,nAtt,uDimID)
  status = nf90_inq_varid(ncid,"fzn",varID)
  status = nf90_get_var(ncid,varID,v)
  status = nf90_close(ncid)



  ! status = nf90_create('abudget'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  ! status = nf90_def_dim(ncid2,"time",msize,it)
  ! status = nf90_def_dim(ncid2,"height",kmax,iz)

  ! status = nf90_def_var(ncid2,"aterm1_fluxdiv",nf90_float,(/it,iz/),vid_a1)
  ! status = nf90_def_var(ncid2,"aterm3_damp",   nf90_float,(/it,iz/),vid_a3)
  ! status = nf90_def_var(ncid2,"aterm4_mix",    nf90_float,(/it,iz/),vid_a4)
  ! status = nf90_def_var(ncid2,"atot",          nf90_float,(/it,iz/),vid_atot)
  ! status = nf90_def_var(ncid2,"dadt",          nf90_float,(/it,iz/),vid_dadt)

  ! status = nf90_put_att(ncid2,vid_a1,"units","A s-1")
  ! status = nf90_put_att(ncid2,vid_a3,"units","A s-1")
  ! status = nf90_put_att(ncid2,vid_a4,"units","A s-1")
  ! status = nf90_put_att(ncid2,vid_atot,"units","A s-1")
  ! status = nf90_put_att(ncid2,vid_dadt,"units","A s-1")

  ! status = nf90_enddef(ncid2)

  ! status = nf90_put_var(ncid2,vid_a1,  aterm1_out)
  ! status = nf90_put_var(ncid2,vid_a3,  aterm3_out)
  ! status = nf90_put_var(ncid2,vid_a4,  aterm4_out)
  ! status = nf90_put_var(ncid2,vid_atot,atot_out)
  ! status = nf90_put_var(ncid2,vid_dadt,dadt_out)

  ! status = nf90_close(ncid2)



  ! status = nf90_create('ubudget'//trim(fname)//'.nc',nf90_noclobber,ncid2)
  ! status = nf90_def_dim(ncid2,"time",msize,it)
  ! status = nf90_def_dim(ncid2,"height",kmax,iz)

  ! status = nf90_def_var(ncid2,"uterm1_epfd",nf90_float,(/it,iz/),vid_u1)
  ! status = nf90_def_var(ncid2,"uterm2_rad", nf90_float,(/it,iz/),vid_u2)
  ! status = nf90_def_var(ncid2,"uterm4_mix", nf90_float,(/it,iz/),vid_u4)
  ! status = nf90_def_var(ncid2,"utot",       nf90_float,(/it,iz/),vid_utot)
  ! status = nf90_def_var(ncid2,"dudt",       nf90_float,(/it,iz/),vid_dudt)

  ! status = nf90_put_att(ncid2,vid_u1,"units","m s-2")
  ! status = nf90_put_att(ncid2,vid_u2,"units","m s-2")
  ! status = nf90_put_att(ncid2,vid_u4,"units","m s-2")
  ! status = nf90_put_att(ncid2,vid_utot,"units","m s-2")
  ! status = nf90_put_att(ncid2,vid_dudt,"units","m s-2")

  ! status = nf90_enddef(ncid2)

  ! status = nf90_put_var(ncid2,vid_u1,  uterm1_out)
  ! status = nf90_put_var(ncid2,vid_u2,  uterm2_out)
  ! status = nf90_put_var(ncid2,vid_u4,  uterm4_out)
  ! status = nf90_put_var(ncid2,vid_utot,utot_out)
  ! status = nf90_put_var(ncid2,vid_dudt,dudt_out)

  ! status = nf90_close(ncid2)


  ! call system_clock(count_end)
  ! runtime = real(count_end - count_start) / real(count_rate)

  ! secs = int(runtime)
  ! hrs = secs / 3600
  ! mins = mod(secs, 3600) / 60
  ! secs = mod(secs, 60)
  ! write(6,*) 'Total run time: ', hrs, 'h ', mins, 'm ', secs, 's'

  call system_clock(count=cc1)

  if (cc1 >= c0) then
    ticks = cc1 - c0
  else
    !  c0  cmax 0  c1
    ticks = (cmax - c0) + cc1 + 1_IK
  end if
  
  elapsed = real(ticks, kind=8) / real(rate, kind=8)
  
  hrs  = int(elapsed) / 3600
  mins = mod(int(elapsed), 3600) / 60
  secs = mod(int(elapsed), 60)
  
  write(6,'(A,I0,A,1X,I0,A,1X,I0,A)') 'Total run time: ', hrs, ' h', mins, ' m', secs, ' s'

  stop

  ! contains

  ! subroutine ck(stat, where)
  !   use netcdf
  !   implicit none
  !   integer, intent(in) :: stat
  !   character(len=*), intent(in) :: where
  !   if (stat /= nf90_noerr) then
  !     write(*,*) 'NETCDF ERROR at: ', trim(where)
  !     write(*,*) trim(nf90_strerror(stat))
  !     stop 1
  !   end if
  ! end subroutine ck
  






end
