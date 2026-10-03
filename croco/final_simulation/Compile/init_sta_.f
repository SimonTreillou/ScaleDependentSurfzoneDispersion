       subroutine init_sta
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=300, MMm0=1000, N=10)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=6,  NP_ETA=20,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      real D_wetdry
      parameter (D_wetdry=0.05D0)
       integer*4 NS
       parameter (NS=5)
       integer*4 Msta
       parameter (Msta=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=35)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      integer*4   ntrc_mld
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=0)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=1)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , itpas
      parameter (itpas=itemp+ntrc_salt+ntrc_mld+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=15)
      parameter (ntrc_surf=0)
      real h(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real hinv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real f(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real fomn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real xp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real xr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real pmon_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real grdscl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask_wet(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmask_wet(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real umask_wet(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real vmask_wet(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rmask_wet_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real Dcrit(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real wetdry(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mask_r_wet/rmask_wet /mask_p_wet/pmask_wet
      common /mask_u_wet/umask_wet /mask_v_wet/vmask_wet
      common /mask_r_wet_avg/rmask_wet_avg
      common /Dcrit_wet/Dcrit
      common /wetdry_wet/wetdry
      real zob(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /Z0B_VAR/zob
      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real wz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,3)
      common /ocean_wz/wz
      real Hzr(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hzr/Hzr
      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      integer*4 stafield
      parameter(stafield=7)
      integer*4 indxstaGrd, indxstaTemp, indxstaSalt,
     &        indxstaRho, indxstaVel, indxstaVrt
     &      , indxstaPtr
      parameter (     indxstaGrd=1, indxstaTemp=2,
     &                indxstaSalt=3, indxstaRho=4,  indxstaVel=5,
     &                indxstaVrt=6, indxstaPtr=7)
      integer*4 ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
     &      , staVrt
     &      , staX,       staY
     &      , staDepth,   staDen
     &      , staTemp
     &      , staPtr(NT)
      logical wrtsta(stafield)
      common/incscrum_sta/
     &        ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
     &      , staVrt
     &      , staX,       staY
     &      , staDepth,   staDen
     &      ,   staTemp
     &      , staPtr
     &      , wrtsta
      character*80  staname,   staposname
      common /cncscrum_sta/ staname,   staposname
      integer*4 NSTAVARS,
     &        istagrd,           istatstr,
     &        istaxgrd,          istaygrd,        istazgrd,
     &        istalon,           istalat,         istadpt,
     &        istatem,           istasal,         istaden,
     &        istau,             istav,           istaz,
     &        istavrt
     &      , istaptr, nptr_sta
      parameter (nptr_sta=NT-1)
      parameter (NSTAVARS=12+nptr_sta,
     &        istagrd=-1,        istatstr=0,
     &        istaxgrd=1,        istaygrd=2,      istazgrd=3,
     &        istalon=4,         istalat=5,       istadpt=6,
     &        istatem=7,         istasal=8,       istaden=9,
     &        istau=10,          istav=11,        istaz=12,
     &        istavrt=13
     &      , istaptr=14
     &        )
      logical diagsta
      integer*4 nstas0,nstas, stagrd(Msta)
      common /stan/ nstas0,nstas, diagsta
      real stainfo(istagrd:istazgrd,Msta)
      common /sta_info/ stainfo
      real staspval, stadeltap2c
      common /sta_scalars/ staspval, stadeltap2c
      real stadata(1:NSTAVARS,Msta), staSigm(1:NSTAVARS,Msta,N)
      common /sta_data/ stadata, staSigm
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real R0,T0,S0, Tcoef, Scoef
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      integer*4 nwrtdiags_eddy
      integer*4 ntsdiags_eddy_avg, nwrtdiags_eddy_avg
      integer*4 nsta, nrpfsta
      logical ldefhis
      logical got_tini(NT)
      logical ldefdiags_eddy
      logical ldefdiags_eddy_avg
      logical ldefsta
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , nsta, nrpfsta
     &                      , got_tini
     &                      , ldefdiags_eddy, nwrtdiags_eddy
     &                      , ldefdiags_eddy_avg
     &                      , nwrtdiags_eddy_avg
     &                      , ntsdiags_eddy_avg
     &                      , ldefsta
     &                      , ldefhis
      real Akv_bak
      common /scalars_akv/ Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akt_bak
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      real stat0(Msta), stax0(Msta), stay0(Msta), staz0(Msta)
      common /ncrealsta/ stat0, stax0, stay0, staz0
      integer*4  stacoor(Msta), STgrd(Msta)
      common /ncintsta/ stacoor,STgrd
      character*80 statitle
      common /nccharsta/ statitle
      integer*4 Ncount, i, icard, iunit, j, k, lstr, mc, nc
      integer*4 index(Msta)
      integer*4 lenstr
      real xsta, ysta, zsta
      real Ip0(Msta), Jp0(Msta), lat(Msta), lon(Msta)
      character*35 frmt
      parameter (iunit=50)
      k=0
      lstr=lenstr(staposname)
      open(iunit,file=staposname(1:lstr),form='formatted',
     &  status='old', err=195)
      icard=0
      do while (icard.lt.99)
        if (icard.eq.1) then
          read(iunit,'(a)',err=70) statitle
          lstr=lenstr(statitle)
        elseif (icard.eq.2) then
          i=0
          nstas=0
          do while (.true.)
            i=i+1
            read(iunit,*,err=40) stat0(i), stax0(i), stay0(i), staz0(i),
     &                           STgrd(i),stacoor(i)
            nstas=i
          enddo
  40      if (stat0(i).ne.99.0D0) then
            write(stdout,50) icard, i, staposname
  50        format(/,
     &            ' INIT_STATIONS - error while reading input card: ',
     &            i2, ', STATIONS location entry: ',i3,/,15x,
     &            'in input script: ',a)
            stop
          else
            Ncount=i-1
            goto 90
          endif
          if (Msta.lt.nstas) then
            write(stdout,60) Msta, nstas
  60        format(/,' INIT_STATIONS - too small dimension parameter,',
     &             ' Msta',2i6,/,15x,'change file  param.h  and',
     &             ' recompile.')
            stop
          endif
        endif
        read(iunit,*,err=70) icard
      enddo
      goto 90
  70  write(stdout,80) icard, staposname
  80  format(/,' INIT_STATIONS - error while reading input card: ',
     &       i2,15x,'in input script: ',a)
  90  close(iunit)
      if (mynode.eq.0) write(stdout,100) nstas
 100  format(2x,' INIT_STATIONS -- Number of stations (nstas) : ',i6,/)
      mc=0
      nc=0
      do i=1,Ncount
          nc=nc+1
          stainfo(istatstr,nc)=(tdays+stat0(i))*day2sec
          stainfo(istazgrd,nc)=staz0(i)
          stainfo(istagrd,i)=FLOAT(STgrd(i))
          if (stacoor(i).eq.0) then
            stax0(i)=stax0(i)-iminmpi+1
            stay0(i)=stay0(i)-jminmpi+1
            if (stax0(i).ge.1 .and. stax0(i).lt.Lmmpi+1 .and.
     &          stay0(i).ge.1 .and. stay0(i).lt.Mmmpi+1) then
              stainfo(istaxgrd,nc)=MIN(MAX(0.5D0,stax0(i)),
     &                                     FLOAT(Lmmpi)+0.5D0)
              stainfo(istaygrd,nc)=MIN(MAX(0.5D0,stay0(i)),
     &                                     FLOAT(Mmmpi)+0.5D0)
            else
              nc=nc-1
            endif
          else
          endif
      enddo
      nstas0=nstas
      nstas=nc
      if (mc.gt.0) then
          call hindices (Ip0,Jp0,lon,lat,mc,xr,yr)
        do i=1,mc
          nc=index(i)
          stainfo(istaxgrd,nc)=MIN(MAX(0.5D0,Ip0(i)),FLOAT(Lmmpi)+0.5D0)
          stainfo(istaygrd,nc)=MIN(MAX(0.5D0,Jp0(i)),FLOAT(Mmmpi)+0.5D0)
        enddo
      endif
      return
 195  write(stdout,205) staposname
 205  format(/,'station file ',A,/,
     &   ' not found => station initialization using restart file')
      return
      end
      subroutine hindices (Ipos,Jpos,Xpos,Ypos,Npos,Xgrd,Ygrd)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=300, MMm0=1000, N=10)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=6,  NP_ETA=20,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      real D_wetdry
      parameter (D_wetdry=0.05D0)
       integer*4 NS
       parameter (NS=5)
       integer*4 Msta
       parameter (Msta=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=35)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      integer*4   ntrc_mld
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=0)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=1)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , itpas
      parameter (itpas=itemp+ntrc_salt+ntrc_mld+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=15)
      parameter (ntrc_surf=0)
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real R0,T0,S0, Tcoef, Scoef
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      integer*4 nwrtdiags_eddy
      integer*4 ntsdiags_eddy_avg, nwrtdiags_eddy_avg
      integer*4 nsta, nrpfsta
      logical ldefhis
      logical got_tini(NT)
      logical ldefdiags_eddy
      logical ldefdiags_eddy_avg
      logical ldefsta
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , nsta, nrpfsta
     &                      , got_tini
     &                      , ldefdiags_eddy, nwrtdiags_eddy
     &                      , ldefdiags_eddy_avg
     &                      , nwrtdiags_eddy_avg
     &                      , ntsdiags_eddy_avg
     &                      , ldefsta
     &                      , ldefhis
      real Akv_bak
      common /scalars_akv/ Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akt_bak
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      logical found
      logical Try_Range
      integer*4 Imax, Imin, Jmax, Jmin, Npos, i0, j0, k
      real Ipos(Npos), Jpos(Npos), Xpos(Npos), Ypos(Npos),
     &  deltax, deltay, dx1, dy1, dx2, dy2, c1, c2
      real Xgrd(-2:Lm+3+padd_X,-2:Mm+3+padd_E),
     &     Ygrd(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      do k=1,Npos
        Ipos(k)=0.0D0
        Jpos(k)=0.0D0
      enddo
      do k=1,Npos
        Imin=0
        Imax=Lmmpi+1
        Jmin=0
        Jmax=Mmmpi+1
        found=Try_Range(Imin,Imax,Jmin,Jmax,Xpos(k),Ypos(k),Xgrd,Ygrd)
        if (found) then
          do while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
            if ((Imax-Imin).gt.1) then
              i0=(Imin+Imax)/2
              found=Try_Range(Imin,i0,Jmin,Jmax,Xpos(k),Ypos(k),
     &                        Xgrd,Ygrd)
              if (found) then
                Imax=i0
              else
                Imin=i0
              endif
            endif
            if ((Jmax-Jmin).gt.1) then
              j0=(Jmin+Jmax)/2
              found=Try_Range(Imin,Imax,Jmin,j0,Xpos(k),Ypos(k),
     &                        Xgrd,Ygrd)
              if (found) then
                Jmax=j0
              else
                Jmin=j0
              endif
            endif
          enddo
          dy1=Ygrd(Imin,Jmin+1)-Ygrd(Imin,Jmin)
          dx1=Xgrd(Imin,Jmin+1)-Xgrd(Imin,Jmin)
          dy2=Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin)
          dx2=Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin)
          c1=Xpos(k)        *dy1-Ypos(k)        *dx1
          c2=Xgrd(Imin,Jmin)*dy2-Ygrd(Imin,Jmin)*dx2
          deltax=(c1*dx2-c2*dx1)/(dx2*dy1-dy2*dx1)
          deltax=(deltax-Xgrd(Imin,Jmin))/dx2
          Ipos(k)=FLOAT(Imin)+MIN(MAX(0.0D0,deltax),1.0D0)
          c1=Xgrd(Imin,Jmin)*dy1-Ygrd(Imin,Jmin)*dx1
          c2=Xpos(k)        *dy2-Ypos(k)        *dx2
          deltay=(c1*dy2-c2*dy1)/(dx2*dy1-dy2*dx1)
          deltay=(deltay-Ygrd(Imin,Jmin))/dy1
          Jpos(k)=FLOAT(Jmin)+MIN(MAX(0.0D0,deltay),1.0D0)
        endif
      enddo
      return
      end
      function Try_Range (Imin,Imax,Jmin,Jmax,Xo,Yo,Xgrd,Ygrd)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=300, MMm0=1000, N=10)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=6,  NP_ETA=20,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      real D_wetdry
      parameter (D_wetdry=0.05D0)
       integer*4 NS
       parameter (NS=5)
       integer*4 Msta
       parameter (Msta=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=35)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      integer*4   ntrc_mld
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=0)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=1)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , itpas
      parameter (itpas=itemp+ntrc_salt+ntrc_mld+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=15)
      parameter (ntrc_surf=0)
      logical Try_Range, inside
      integer*4 Imax, Imin, Jmax, Jmin, Nb, NX, i, j, shft
      parameter(NX=2*(LLm0+2)+2*(MMm0+2)+1)
      real Xb(NX), Yb(NX), Xo, Yo
      real Xgrd(-2:Lm+3+padd_X,-2:Mm+3+padd_E),
     &     Ygrd(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      do i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      enddo
      shft=1-Jmin+Imax-Imin
      do j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      enddo
      shft=1+Jmax-Jmin+2*Imax-Imin
      do i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      enddo
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      do j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      enddo
      Try_Range=inside(Xo,Yo,Xb,Yb,Nb)
      return
      end
      function inside (Xo,Yo,Xb,Yb,Nb)
      implicit none
      logical inside
      integer*4 Nb, Nstep, crossings, i, inc, k, kk, nc
      parameter (Nstep=128)
      integer*4 index(Nstep)
      real Xb(Nb+1), Yb(Nb+1), Xo, Yo, dx1, dx2, dxy
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
      do kk=0,Nb-1,Nstep
        nc=0
        do k=kk+1,MIN(kk+Nstep,Nb)
          if (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0D0).and.
     &        (Xb(k).ne.Xb(k+1))) then
            nc=nc+1
            index(nc)=k
          endif
        enddo
        do i=1,nc
          k=index(i)
          if (Xb(k).ne.Xb(k+1)) then
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            if ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) then
              crossings=1
              goto 10
            elseif (((dx1.eq.0.0D0).and.(Yo.ge.Yb(k  ))).or.
     &              ((dx2.eq.0.0D0).and.(Yo.ge.Yb(k+1)))) then
              inc=1
            elseif ((dx1*dx2.gt.0.0D0).and.
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0D0)) then
              inc=2
            endif
            if (Xb(k+1).gt.Xb(k)) then
              crossings=crossings+inc
            else
              crossings=crossings-inc
            endif
          endif
        enddo
      enddo
  10  if (crossings.eq.0) then
        inside=.false.
      else
        inside=.true.
      endif
      return
      end
