      subroutine vnbq_bc_tile (Istr,Iend,Jstr,Jend,grad)
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
      real zeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real ubar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real vbar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
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
      logical M2bc_nbq_flag
      common /nbq_M2bc/ M2bc_nbq_flag
      integer*4 iteration_nbq_max
      common /nbq_var1/ iteration_nbq_max
      integer*4 iteration_nbq
      common /nbq_var2/ iteration_nbq
      integer*4 ifl_nbq
      common /nbq_var3/ ifl_nbq
      integer*4 slip_nbq
      common /nbq_var4/ slip_nbq
      real soundspeed_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_param1/ soundspeed_nbq
      real soundspeed2_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_param2/ soundspeed2_nbq
      double precision time_nbq
      common /nbq_param3/ time_nbq
      double precision csvisc1_nbq
      common /nbq_param4/ csvisc1_nbq
      double precision csvisc2_nbq
      common /nbq_param5/ csvisc2_nbq
      double precision cw_int_nbq
      common /nbq_param6/ cw_int_nbq
      double precision ifl_imp_nbq
      common /nbq_param7/ ifl_imp_nbq
      integer*4 ndtnbq
      common /time_nbq1/ ndtnbq
      real dtnbq
      common /time_nbq2/ dtnbq
      real csound_nbq
      common /nbq_csound/ csound_nbq
      real visc2_nbq
      common /nbq_visc2/ visc2_nbq
      real dtgrid_nbq
      common /nbq_dtgrid/ dtgrid_nbq
      real qdmu_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_qdmu_nbq/ qdmu_nbq
      real qdmv_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_qdmv_nbq/ qdmv_nbq
      real qdmw_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /nbq_qdmw_nbq/ qdmw_nbq
      real thetadiv_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_thetadiv_nbq/ thetadiv_nbq
      real thetadiv2_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_thetadiv2_nbq/ thetadiv2_nbq
      real thetadiv3_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_thetadiv3_nbq/ thetadiv3_nbq
      real ru_int_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_ruint/ ru_int_nbq
      real rv_int_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_rvint/ rv_int_nbq
      real ru_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_ru/ ru_nbq
      real rv_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /nbq_rv/ rv_nbq
      real ru_nbq_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /avg2_runbq/ ru_nbq_avg2
      real rv_nbq_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /avg2_rvnbq/ rv_nbq_avg2
      real rw_int_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /nbq_rwint/ rw_int_nbq
      real rw_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /nbq_rw/ rw_nbq
      real rw_nbq_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /avg2_rwnbq/ rw_nbq_avg2
      real rho_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common/nbq_rho_nbq/rho_nbq
      integer*4 inc_faststep
      common/nbq_inc_faststep/inc_faststep
      integer*4 nb_faststep
      common/nbq_nb_faststep/nb_faststep
      real DU_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_DU_nbq/ DU_nbq
      real DV_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_DV_nbq/ DV_nbq
      real ru_int_nbq_2d (-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_ruint_2d/ru_int_nbq_2d
      real rv_int_nbq_2d (-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_rvint_2d/rv_int_nbq_2d
      real rho_grd(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common/nbq_rho_grd/rho_grd
      real rho_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common/nbq_rho_bak/rho_bak
      real rho_nbq_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /avg1_rhonbq/ rho_nbq_avg1
      real rhobar_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      common /nbq_rhobar/ rhobar_nbq
      real rhobar_nbq_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_rhobar_AVG1/ rhobar_nbq_avg1
      real Hzw_half_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hzw_half_nbq/ Hzw_half_nbq
      real zw_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,4)
      common /nbq_zw/ zw_nbq
       real Hz_correct(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
       common /grid_Hz_correct/ Hz_correct
      real dthetadiv_nbqdz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /nbq_nods3/ dthetadiv_nbqdz
      real dZdxq_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N+1)
      common /nbq_nods5/ dZdxq_w
      real dZdyq_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N+1)
      common /nbq_nods7/ dZdyq_w
      real wsurf_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_wsurf/ wsurf_nbq
      real usurf_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_usurf/ usurf_nbq
      real vsurf_nbq(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_vsurf/ vsurf_nbq
      real qdmu_nbq_west(-2:Mm+3+padd_E,N,2)
      common /bry_unbq_west/ qdmu_nbq_west
      real qdmv_nbq_west(-2:Mm+3+padd_E,N,2)
      common /bry_vnbq_west/ qdmv_nbq_west
      real qdmw_nbq_west(-2:Mm+3+padd_E,0:N,2)
      common /bry_wnbq_west/ qdmw_nbq_west
      real  rho_nbq_west(-2:Mm+3+padd_E,N,2)
      common /bry_rnbq_west/ rho_nbq_west
      real NBQnudgcof(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_nudg/ NBQnudgcof
      real unbqclm(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real vnbqclm(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /climat_unbqclm/unbqclm
      common /climat_vnbqclm/vnbqclm
      real wnbqclm(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real rnbqclm(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /climat_wnbqclm/wnbqclm
      common /climat_rnbqclm/rnbqclm
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
      logical got_tbry(NT)
      common /bry_logical/ got_tbry
      real zetabry_west(-2:Mm+3+padd_E),
     &    zetabry_west_dt(-2:Mm+3+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(-2:Mm+3+padd_E),
     &    ubarbry_west_dt(-2:Mm+3+padd_E,2)
     &    ,vbarbry_west(-2:Mm+3+padd_E),
     &    vbarbry_west_dt(-2:Mm+3+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real ubry_west(-2:Mm+3+padd_E,N),
     &    ubry_west_dt(-2:Mm+3+padd_E,N,2)
     &    ,vbry_west(-2:Mm+3+padd_E,N),
     &    vbry_west_dt(-2:Mm+3+padd_E,N,2)
      common /bry_u_west/ ubry_west, ubry_west_dt,
     &                    vbry_west, vbry_west_dt
      real tbry_west(-2:Mm+3+padd_E,N,NT),
     &     tbry_west_dt(-2:Mm+3+padd_E,N,2,NT)
      common /bry_t_west/ tbry_west, tbry_west_dt
      real unbqbry_west(-2:Mm+3+padd_E,N),
     &     vnbqbry_west(-2:Mm+3+padd_E,N)
      common /bry_nbq_west1/ unbqbry_west,
     &                       vnbqbry_west
      real wnbqbry_west(-2:Mm+3+padd_E,0:N),
     &     rnbqbry_west(-2:Mm+3+padd_E,N)
      common /bry_nbq_west2/ wnbqbry_west,
     &                       rnbqbry_west
      real wbry_west(-2:Mm+3+padd_E,0:N)
      common /bry_w_west/ wbry_west
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2), cff,eps,cff1,cff2,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.D-20)
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        JstrR=Jstr-2
      else
        JstrR=Jstr
      endif
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        JendR=Jend+2
      else
        JendR=Jend
      endif
      tau_in =dtfast*tauM_in
      tau_out=dtfast*tauM_out
!$acc kernels if(compute_on_device) default(present)
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        do k=1,N
          do j=Jstr,Jend
            qdmv_nbq(Istr-1,j,k)=vnbqbry_west(j,k)
     &                          *vmask_wet(Istr-1,j)
          enddo
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        do k=1,N
          do j=Jstr,Jend
            qdmv_nbq(Iend+1,j,k)=gamma2*qdmv_nbq(Iend,j,k)
          enddo
        enddo
      endif
!$acc end kernels
      return
      end
