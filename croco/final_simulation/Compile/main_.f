      program main
      use buffer, only: init_buffer
      integer*4 tile, subs, trd, ierr
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
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,9,0:NPP-1)
     &    ,A3dHz(N3dHz,4,0:NPP-1)
      common/private_scratch/ A2d,A3d,A3dHz
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
      real NBQnudgcof(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nbq_nudg/ NBQnudgcof
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
     &       ,filetype_abl, filetype_abl_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16,
     &           filetype_abl=18, filetype_abl_avg=19)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4, dimension(ntrc_pas) :: indxTPAS
     & =(/(iloop,iloop=indxV+ntrc_temp+ntrc_salt+ntrc_mld+1,
     &  indxV+ntrc_temp+ntrc_salt+ntrc_mld+ntrc_pas)/)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_mld+ntrc_pas+
     &           ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxeddyzz,
     &        indxeddyuu,indxeddyvv,indxeddyuv,indxeddyub,
     &        indxeddyvb,indxeddywb,indxeddyuw,indxeddyvw,
     &        indxeddyubu,indxeddyvbv,
     &        indxeddyusu,indxeddyvsv,
     &        indxeddyugsu,indxeddyvgsv
      parameter (indxeddyzz=indxV+ntrc_temp+ntrc_salt+ntrc_mld
     &                           +ntrc_pas+ntrc_bio+ntrc_sed
     &                           +ntrc_diats+ntrc_diauv+ntrc_diavrt
     &                           +ntrc_diaek+ntrc_diapv+400,
     &           indxeddyuu=indxeddyzz+1,
     &           indxeddyvv=indxeddyuu+1,
     &           indxeddyuv=indxeddyvv+1,
     &           indxeddyub=indxeddyuv+1,
     &           indxeddyvb=indxeddyub+1,
     &           indxeddywb=indxeddyvb+1,
     &           indxeddyuw=indxeddywb+1,
     &           indxeddyvw=indxeddyuw+1,
     &           indxeddyubu=indxeddyvw+1,
     &           indxeddyvbv=indxeddyubu+1,
     &           indxeddyusu=indxeddyvbv+1,
     &           indxeddyvsv=indxeddyusu+1,
     &           indxeddyugsu=indxeddyvsv+1,
     &           indxeddyvgsv=indxeddyugsu+1)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt
     &                      +ntrc_mld+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxTke
      parameter (indxTke=indxAkv+ntrc_temp+7)
      integer*4 indxGls
      parameter (indxGls=indxAkv+ntrc_temp+8)
      integer*4 indxLsc
      parameter (indxLsc=indxAkv+ntrc_temp+9)
      integer*4 indxAkk
      parameter (indxAkk=indxAkv+ntrc_temp+10)
      integer*4 indxAkp
      parameter (indxAkp=indxAkv+ntrc_temp+11)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      parameter (indxShflx_rsw=indxShflx+1)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstAkv,rstAkt
      integer*4 rstTke,rstGls
      integer*4 rstBustr, rstBvstr
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
      integer*4 hisT(NT)
      integer*4 nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy
     &      , diags_eddyTime, diags_eddyTime2, diags_eddyTstep
     &      , diags_eddyzz(2)
     &      , diags_eddyuu(2), diags_eddyvv(2), diags_eddyuv(2)
     &      , diags_eddyub(2), diags_eddyvb(2), diags_eddywb(2)
     &      , diags_eddyuw(2), diags_eddyvw(2)
     &      , diags_eddyubu(2), diags_eddyvbv(2)
     &      , diags_eddyusu(2), diags_eddyvsv(2)
     &      , diags_eddyugsu(2), diags_eddyvgsv(2)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
      integer*4 avgT(NT)
       integer*4 nciddiags_eddy_avg, nrecdiags_eddy_avg
     &      , nrpfdiags_eddy_avg
     &      , diags_eddyTime_avg, diags_eddyTime2_avg
     &      , diags_eddyTstep_avg
     &      , diags_eddyzz_avg(2)
     &      , diags_eddyuu_avg(2), diags_eddyvv_avg(2)
     &      , diags_eddyuv_avg(2)
     &      , diags_eddyub_avg(2), diags_eddyvb_avg(2)
     &      , diags_eddywb_avg(2)
     &      , diags_eddyuw_avg(2), diags_eddyvw_avg(2)
     &      , diags_eddyubu_avg(2), diags_eddyvbv_avg(2)
     &      , diags_eddyusu_avg(2), diags_eddyvsv_avg(2)
     &      , diags_eddyugsu_avg(2), diags_eddyvgsv_avg(2)
      logical wrthis(1000+NT)
     &      , wrtavg(1000+NT)
     &      , wrtdiags_eddy(3)
     &      , wrtdiags_eddy_avg(3)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstAkv,rstAkt
     &      , rstTke,rstGls
     &      , rstBustr,rstBvstr
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
     &      , nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy
     &      , diags_eddyTime, diags_eddyTstep
     &      , diags_eddyzz
     &      , diags_eddyuu, diags_eddyvv, diags_eddyuv, diags_eddyub
     &      , diags_eddyvb, diags_eddywb, diags_eddyuw, diags_eddyvw
     &      , diags_eddyubu, diags_eddyvbv
     &      , diags_eddyusu, diags_eddyvsv
     &      , diags_eddyugsu, diags_eddyvgsv
     &      , nciddiags_eddy_avg, nrecdiags_eddy_avg
     &      , nrpfdiags_eddy_avg
     &      , diags_eddyTime_avg, diags_eddyTime2_avg
     &      , diags_eddyTstep_avg
     &      , diags_eddyzz_avg
     &      , diags_eddyuu_avg, diags_eddyvv_avg, diags_eddyuv_avg
     &      , diags_eddyub_avg, diags_eddyvb_avg, diags_eddywb_avg
     &      , diags_eddyuw_avg, diags_eddyvw_avg
     &      , diags_eddyubu_avg, diags_eddyvbv_avg
     &      , diags_eddyusu_avg, diags_eddyvsv_avg
     &      , diags_eddyugsu_avg, diags_eddyvgsv_avg
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
     &      , wrthis
     &      , wrtavg
     &      , wrtdiags_eddy
     &      , wrtdiags_eddy_avg
      character*80 title
      character*80 origin_date, start_date_run, xios_origin_date
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8) :: origin_date_in_sec, xios_origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,  diags_eddyname
     &                                ,  diags_eddyname_avg
     &                                ,   bry_file
      character*75  vname(20, 1000)
      common /cncscrum/   title
     &         ,   origin_date, start_date_run
     &         ,   xios_origin_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   xios_origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,  diags_eddyname
     &                                ,  diags_eddyname_avg
     &                                ,   bry_file
     &                                ,   vname
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
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
      include 'mpif.h'
      integer*4 :: iifroot, iicroot
      call MPI_Init (ierr)
      call MPI_Setup (ierr)
      if (ierr.ne.0) goto 100
      call init_buffer
      call read_inp (ierr)
      if (ierr.ne.0) goto 100
      call init_scalars (ierr)
      if (ierr.ne.0) goto 100
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call init_arrays (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call ana_grid (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call setup_grid1 (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call setup_grid2 (tile)
         enddo
       enddo
      call set_scoord
      call set_weights
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call set_depth (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call grid_stiffness (tile)
         enddo
       enddo
      if (nrrec.ne.0) then
        call get_initial (nrrec, 1)
      else
      endif
      time=start_time
      tdays=time*sec2day
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call ana_initial (tile)
         enddo
       enddo
      if (may_day_flag.ne.0) goto 99
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call set_depth (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call set_HUV (tile)
        call reset_rho0 (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call omega (tile)
        call rho_eos (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call set_nudgcof (tile)
         enddo
       enddo
      call get_vbc
      if (may_day_flag.ne.0) goto 99
      nrecsta=0
      ncidsta=-1
      staspval=1.D15
      stadeltap2c=2.5D0
      call init_arrays_sta
      call init_sta
      call interp_r2d_sta_ini (  xr(-2,-2), istalon)
      call interp_r2d_sta_ini (  yr(-2,-2), istalat)
      call fill_sta_ini
      if (ldefsta) call wrt_sta
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call set_depth (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*trd,subs*(trd+1)-1,+1
        call set_HUV (tile)
         enddo
       enddo
C$OMP PARALLEL DO PRIVATE(trd,subs, tile), SHARED(numthreads)
       do trd=0,numthreads-1
        subs=NSUB_X*NSUB_E/numthreads
         do tile=subs*(trd+1)-1,subs*trd,-1
        call initial_nbq(tile)
         enddo
       enddo
      if (ldefhis .and. wrthis(indxTime)) call wrt_his
      if (may_day_flag.ne.0) goto 99
      if (mynode.eq.0) write(stdout,'(/1x,A27/)')
     &                'MAIN: started time-stepping.'
      next_kstp=kstp
      time_start=time
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      iif = -1
      nbstep3d = 0
      iic = ntstart
      do iicroot=ntstart,ntimes+1
        do iifroot = 0,nfast+2
          call step()
        enddo
        if (may_day_flag.ne.0) goto 99
      enddo
  99  continue
      call closecdf
 100  continue
      if (ierr.ne.0) call mpi_abort (MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Finalize (ierr)
      stop
      end
