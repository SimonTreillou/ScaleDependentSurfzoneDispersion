      subroutine wrt_avg
      integer*4 ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite, type
     &            , itrc
      include 'mpif.h'
      integer*4 status(MPI_STATUS_SIZE), blank
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
      real zeta_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real ubar_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real vbar_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_zeta/zeta_avg
     &       /avg_ubar/ubar_avg
     &       /avg_vbar/vbar_avg
      real bostr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_bostr/bostr_avg
      real bustr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_bustr/bustr_avg
      real bvstr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_bvstr/bvstr_avg
      real wstr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_wstr/wstr_avg
      real sustr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_sustr/sustr_avg
      real svstr_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_svstr/svstr_avg
      real srflx_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_srflx/srflx_avg
      real u_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real v_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real t_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,NT)
      real rho_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real bvf_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real omega_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real w_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &       /avg_rho/rho_avg /avg_omega/omega_avg
     &       /avg_bvf/bvf_avg
     &       /avg_w/w_avg
      real stflx_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /avg_stflx/stflx_avg
      real btflx_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /avg_btflx/btflx_avg
      real hbl_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /avg_hbl/hbl_avg
      real tke_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real gls_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Lscale_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /avg_tke/tke_avg
      common /avg_gls/gls_avg
      common /avg_Lscale/Lscale_avg
      real visc3d_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /avg_visc3d/visc3d_avg
      real Akv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /avg_Akv/Akv_avg /avg_Akt/Akt_avg
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      real work(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /work3d/ work
      real workr(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      common /work3d_r/ workr
      real work2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /work2d/ work2d
      real work2d2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /work2d2/ work2d2
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_64bit_data
      integer*4 nf_cdf5
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_format_64bit_offset
      integer*4 nf_format_64bit_data
      integer*4 nf_format_cdf5
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_int64
      external        nf_put_att_int64
      integer*4         nf_get_att_int64
      external        nf_get_att_int64
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_compact
      parameter (nf_compact = 2)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_szip
      external nf_def_var_szip
      integer*4 nf_inq_var_szip
      external nf_inq_var_szip
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_def_var_filter
      external nf_def_var_filter
      integer*4 nf_inq_var_filter
      external nf_inq_var_filter
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      integer*4 nf_set_log_level
      external nf_set_log_level
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
      call def_avg (ncidavg, nrecavg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(avgname)
      nrecavg=max(nrecavg,1)
      if (nrpfavg.eq.0) then
        record=nrecavg
      else
        record=1+mod(nrecavg-1, nrpfavg)
      endif
      type=filetype_avg
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecavg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncidavg, avgTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
         write(stdout,1) 'time_step', record,ierr ,' mynode =', mynode
         goto 99
      endif
      ierr=nf_put_var1_double (ncidavg, avgTime, record, time_avg
     &  )
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxTime))
         write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &        ,' mynode =', mynode
         goto 99
      endif
      ierr=nf_put_var1_double (ncidavg, avgTime2, record, time_avg
     &  )
      if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxTime2))
         write(stdout,1) vname(1,indxTime2)(1:lvar), record, ierr
     &        ,' mynode =', mynode
         goto 99
      endif
      if (wrtavg(indxZ)) then
!$acc update if(compute_on_device) host(zeta_avg)
         work2d=zeta_avg
         do j=0,Mm
            do i=0,Lm
               if (h(i,j) .le. Dcrit(i,j)) then
                 work2d(i,j)=work2d(i,j)+h(i,j)
               endif
            enddo
         enddo
         call fillvalue2d(work2d,ncidavg,avgZ,vname(1,indxZ),
     &                   record,r2dvar,type)
      endif
      if (wrtavg(indxUb)) then
!$acc update if(compute_on_device) host(ubar_avg)
         work2d=ubar_avg
         call fillvalue2d(work2d,ncidavg,avgUb,vname(1,indxUb),
     &                      record,u2dvar,type)
      endif
      if (wrtavg(indxVb)) then
!$acc update if(compute_on_device) host(vbar_avg)
         work2d=vbar_avg
         call fillvalue2d(work2d,ncidavg,avgVb,vname(1,indxVb),
     &                   record,v2dvar,type)
      endif
      if (wrtavg(indxBostr)) then
!$acc update if(compute_on_device) host(bostr_avg)
         work2d=bostr_avg
         call fillvalue2d(work2d,ncidavg,avgBostr,vname(1,indxBostr),
     &        record,r2dvar,type)
      endif
      if (wrtavg(indxBustr)) then
!$acc update if(compute_on_device) host(bustr_avg)
         work2d=bustr_avg
         call fillvalue2d(work2d,ncidavg,avgBustr,vname(1,indxBustr),
     &        record,u2dvar,type)
      endif
      if (wrtavg(indxBvstr)) then
!$acc update if(compute_on_device) host(bvstr_avg)
         work2d=bvstr_avg
         call fillvalue2d(work2d,ncidavg,avgBvstr,vname(1,indxBvstr),
     &        record,v2dvar,type)
      endif
      if (wrtavg(indxWstr)) then
!$acc update if(compute_on_device) host(wstr_avg)
         ierr=nf_fwrite(wstr_avg, ncidavg, avgWstr, record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxWstr))
            write(stdout,1) vname(1,indxWstr)(1:lvar), record, ierr
     &           ,' mynode =', mynode
            goto 99
         endif
      endif
      if (wrtavg(indxUWstr)) then
!$acc update if(compute_on_device) host(sustr_avg)
         ierr=nf_fwrite(sustr_avg, ncidavg, avgUWstr, record, u2dvar)
        if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxUWstr))
            write(stdout,1) vname(1,indxUWstr)(1:lvar), record, ierr
     &           ,' mynode =', mynode
            goto 99
         endif
      endif
      if (wrtavg(indxVWstr)) then
!$acc update if(compute_on_device) host(svstr_avg)
         ierr=nf_fwrite(svstr_avg, ncidavg, avgVWstr, record, v2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxVWstr))
            write(stdout,1) vname(1,indxVWstr)(1:lvar), record, ierr
     &           ,' mynode =', mynode
            goto 99
         endif
      endif
      if (wrtavg(indxU)) then
!$acc update if(compute_on_device) host(u_avg)
         workr=u_avg
         call fillvalue3d(workr,ncidavg,avgU,
     &        vname(1,indxU),record,u3dvar,type)
      endif
      if (wrtavg(indxV)) then
!$acc update if(compute_on_device) host(v_avg)
         workr=v_avg
         call fillvalue3d(workr,ncidavg,avgV,
     &        vname(1,indxV),record,v3dvar,type)
      endif
      do itrc=1,NT
         if (wrtavg(indxV+itrc)) then
!$acc update if(compute_on_device) host(t_avg(:,:,:,itrc))
            workr=t_avg(:,:,:,itrc)
            call fillvalue3d(workr,ncidavg,avgT(itrc),
     &        vname(1,indxV+itrc), record,r3dvar,type)
         endif
      enddo
      if (wrtavg(indxR)) then
!$acc update if(compute_on_device) host(rho_avg)
         workr=rho_avg+rho0-1000.D0
         call fillvalue3d(workr,ncidavg,avgR,
     &        vname(1,indxR), record,r3dvar,type)
      endif
      if (wrtavg(indxbvf)) then
!$acc update if(compute_on_device) host(bvf_avg)
         work=bvf_avg
         call fillvalue3d_w(work,ncidavg,avgbvf,
     &        vname(1,indxbvf), record,w3dvar,type)
      endif
      if (wrtavg(indxO)) then
!$acc update if(compute_on_device) host(omega_avg)
         work=omega_avg
         call fillvalue3d_w(work,ncidavg,avgO,
     &        vname(1,indxO), record,w3dvar,type)
      endif
      if (wrtavg(indxW)) then
!$acc update if(compute_on_device) host(w_avg)
         work=w_avg
         call fillvalue3d(work ,ncidavg,avgW,
     &        vname(1,indxW), record,w3dvar,type)
      endif
      if (wrtavg(indxVisc)) then
!$acc update if(compute_on_device) host(visc3d_avg)
         workr=visc3d_avg
         call fillvalue3d(workr,ncidavg,avgVisc,
     &        vname(1,indxVisc), record,r3dvar,type)
      endif
      if (wrtavg(indxAkv)) then
!$acc update if(compute_on_device) host(Akv_avg)
         work=Akv_avg
         call fillvalue3d_w(work,ncidavg,avgAkv,
     &        vname(1,indxAkv), record,w3dvar,type)
      endif
      if (wrtavg(indxAkt)) then
!$acc update if(compute_on_device) host(Akt_avg(:,:,:,itemp))
         work=Akt_avg(:,:,:,itemp)
         call fillvalue3d_w(work,ncidavg,avgAkt,
     &        vname(1,indxAkt), record,w3dvar,type)
      endif
      if (wrtavg(indxHbl)) then
!$acc update if(compute_on_device) host(hbl_avg)
         work2d=hbl_avg
         call fillvalue2d(work2d,ncidavg,avgHbl,
     &        vname(1,indxHbl), record,r2dvar,type)
      endif
      if (wrtavg(indxTke)) then
!$acc update if(compute_on_device) host(tke_avg)
         work=tke_avg
         call fillvalue3d_w(work,ncidavg,avgTke,
     &        vname(1,indxTke), record,w3dvar,type)
      endif
      if (wrtavg(indxGls)) then
!$acc update if(compute_on_device) host(gls_avg)
         work=gls_avg
         call fillvalue3d_w(work,ncidavg,avgGls,
     &        vname(1,indxGls), record,w3dvar,type)
      endif
      if (wrtavg(indxLsc)) then
!$acc update if(compute_on_device) host(Lscale_avg)
         work=Lscale_avg
         call fillvalue3d_w(work,ncidavg,avgLsc,
     &        vname(1,indxLsc), record,w3dvar,type)
      endif
      if (wrtavg(indxShflx)) then
!$acc update if(compute_on_device) host(stflx_avg(:,:,itemp))
         work2d=stflx_avg(:,:,itemp)
         ierr=nf_fwrite(work2d, ncidavg, avgShflx, record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxShflx))
            write(stdout,1) vname(1,indxShflx)(1:lvar), record, ierr
     &           ,' mynode =', mynode
            goto 99
         endif
       endif
      if (wrtavg(indxShflx_rsw)) then
!$acc update if(compute_on_device) host(srflx_avg)
         ierr=nf_fwrite(srflx_avg, ncidavg, avgShflx_rsw,
     &                  record, r2dvar)
         if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxShflx_rsw))
            write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), record, ierr
     &           ,' mynode =', mynode
            goto 99
         endif
      endif
 1    format(/' WRT_AVG - ERROR while writing variable(',1x,a,1x,
     &     ')into averages file.',/,11x,'Time record:',
     &     i6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
 99   may_day_flag=3
 100  continue
      if (nrpfavg.gt.0 .and. record.ge.nrpfavg) then
         ierr=nf_close(ncidavg)
         ncidavg=-1
      else
         ierr=nf_sync(ncidavg)
      endif
      if (ierr .eq. nf_noerr) then
      if (mynode.eq.0) write(stdout,'(6x,A,2(A,I4,1x),A,I3)')
     &        'WRT_AVG -- wrote ',
     &        'averaged fields into time record =', record, '/',
     &        nrecavg  ,' mynode =', mynode
      else
       if (mynode.eq.0) write(stdout,'(/1x,2A/)')
     &        'WRT_AVG ERROR: Cannot ',
     &        'synchronize/close averages netCDF file.'
         may_day_flag=3
      endif
      return
      end
