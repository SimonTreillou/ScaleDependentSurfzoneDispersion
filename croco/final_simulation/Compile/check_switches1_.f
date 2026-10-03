      subroutine check_switches1 (ierr)
      implicit none
      integer*4 ierr, is,ie, iexample
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
      parameter (itemp=0)
      parameter (ntrc_temp=0)
      parameter (ntrc_salt=0)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=0)
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
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=15)
      parameter (ntrc_surf=0)
      integer*4 max_opt_size
      parameter (max_opt_size=5000)
      character*5000 Coptions,srcs
      common /strings/ Coptions,srcs
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
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , nsta, nrpfsta
     &                      , ldefdiags_eddy, nwrtdiags_eddy
     &                      , ldefdiags_eddy_avg
     &                      , nwrtdiags_eddy_avg
     &                      , ntsdiags_eddy_avg
     &                      , ldefsta
     &                      , ldefhis
      real Akv_bak
      common /scalars_akv/ Akv_bak
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
      if (mynode.eq.0) write(stdout,'(/1x,A/)')
     &      'Activated C-preprocessing Options:'
      do is=1,max_opt_size
        Coptions(is:is)=' '
      enddo
      iexample=0
      is=1
      iexample=iexample+1
      if (mynode.eq.0) write(stdout,'(10x,A)') 'FLASH_RIP'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FLASH_RIP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MPI'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MPI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_WEST'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'STATIONS'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='STATIONS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NC4PAR'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NC4PAR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WET_DRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WET_DRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NEW_S_COORD'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NEW_S_COORD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SOLVE3D'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SOLVE3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_ADV'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_ADV'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'FRC_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'Z_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='Z_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M2_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M3_FRC_BRY'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'T_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='T_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_HADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_VIS2'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_VADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_HADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'TS_VADV_WENO5'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='TS_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPONGE'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'LIMIT_BSTRESS'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='LIMIT_BSTRESS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_MIXING'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BTFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BTFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_M2CHARACT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M2CHARACT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_M3ORLANSKI'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_M3ORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_TORLANSKI'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_TORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AVERAGES'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'AVERAGES_K'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='AVERAGES_K'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'DIAGNOSTICS_EDDY'
      ie=is +15
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DIAGNOSTICS_EDDY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ALL_SIGMA'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ALL_SIGMA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'W_HADV_WENO5'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_HADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'W_VADV_WENO5'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_VADV_WENO5'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_INITIAL'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_INITIAL'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SSFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SSFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_STFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_STFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_BRY'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M2FILTER_NONE'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M2FILTER_NONE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_GRID'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_GRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SMFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SMFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SRFLUX'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SRFLUX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ANA_SST'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ANA_SST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NS_PERIODIC'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NS_PERIODIC'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_KOMEGA'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_KOMEGA'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_PRECISE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_PRECISE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_MIXING_3D'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_MIXING_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER_SPECTRUM'
      ie=is +18
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_SPECTRUM'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER_DSPREAD'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_DSPREAD'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_SPECIFIED_WEST'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_SPECIFIED_WEST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'MR'
      ie=is + 1
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='MR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPINUP'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPINUP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NO_TRACER'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NO_TRACER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NO_TEMPERATURE'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NO_TEMPERATURE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'BATHY_PERTURBATION'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BATHY_PERTURBATION'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER_FILE'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_FILE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER_OBLIQUE'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_OBLIQUE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'GLS_KEPSILON'
      ie=is +11
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='GLS_KEPSILON'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'CLARK_R3'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CLARK_R3'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'M3FAST'
      ie=is + 5
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='M3FAST'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_IMP'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_IMP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'HZR'
      ie=is + 2
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='HZR'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_MASS'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_MASS'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_HZCORRECT'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_HZCORRECT'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_NBQ'
      ie=is + 6
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NBQ'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_NBQORLANSKI'
      ie=is +14
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_NBQORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'OBC_WORLANSKI'
      ie=is +12
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='OBC_WORLANSKI'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_NUDGING'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_NUDGING'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQCLIMATOLOGY'
      ie=is +13
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQCLIMATOLOGY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NBQ_FRC_BRY'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NBQ_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'W_FRC_BRY'
      ie=is + 8
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='W_FRC_BRY'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VAR_RHO_2D'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VAR_RHO_2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'RESET_RHO0'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='RESET_RHO0'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'PGF_BASIC_JACOBIAN'
      ie=is +17
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='PGF_BASIC_JACOBIAN'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WJ_GRADP'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WJ_GRADP'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'UV_MIX_S'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='UV_MIX_S'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'VIS_COEF_3D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='VIS_COEF_3D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NTRA_T3DMIX'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NTRA_T3DMIX'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPONGE_GRID'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_GRID'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPONGE_DIF2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_DIF2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'SPONGE_VIS2'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='SPONGE_VIS2'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'CANUTO_A'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='CANUTO_A'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'WAVE_MAKER_DSPREAD_PER'
      ie=is +21
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='WAVE_MAKER_DSPREAD_PER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'NF_CLOBBER'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='NF_CLOBBER'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'BD_EXCHANGE'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='BD_EXCHANGE'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'ENDDOLOOP2D'
      ie=is +10
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='ENDDOLOOP2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'DOLOOP2D_R'
      ie=is + 9
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DOLOOP2D_R'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(10x,A)') 'DOLOOP2D'
      ie=is + 7
      if (ie.ge.max_opt_size) goto 99
      Coptions(is:ie)='DOLOOP2D'
      Coptions(ie+1:ie+1)=' '
      is=ie+2
      if (mynode.eq.0) write(stdout,'(/)')
      if (iexample.eq.0) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR in "cppdefs.h": no configuration is specified.'
        ierr=ierr+1
      elseif (iexample.gt.1) then
        if (mynode.eq.0) write(stdout,'(1x,A)')
     & 'ERROR: more than one configuration in "cppdefs.h".'
        ierr=ierr+1
      endif
      return
  99  if (mynode.eq.0) write(stdout,'(/1x,A,A/14x,A)')
     &  'CHECKDEFS -- ERROR: Unsufficient size of string Coptions',
     &  'in file "strings.h".', 'Increase the size it and recompile.'
      ierr=ierr+1
      return
      end
