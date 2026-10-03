      subroutine init_scalars (ierr)
      implicit none
      integer*4 ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3
      integer*4 omp_get_num_threads
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
      real visc2_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real diff2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real visc3d_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_visc3d_r/visc3d_r
      real visc3d_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_visc3d_p/visc3d_p
      real Akv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real Akv_old(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt_old(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_Akvold/Akv_old /mixing_Aktold/Akt_old
      real bvf(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_bvf/ bvf
      real trb(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2,NGLS)
      common /gls_trb/trb
      real Lscale(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /gls_lsc/Lscale
      real Eps_gls(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /gls_eps/Eps_gls
      integer*4 kbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /gls_kbl/ kbl
      real hbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E  )
      common /gls_hbl/ hbl
      real cm0
      common /gls_cm0/ cm0
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
      real timediags_eddy_avg
      real eddyzz_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyuu_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyvv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyuv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyub_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyvb_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddywb_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyuw_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyvw_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real eddyubu_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyvbv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyusu_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyvsv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyugsu_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real eddyvgsv_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /diag_timediags_eddy_avg/timediags_eddy_avg
      common /diag_eddyzz_avg/eddyzz_avg
     &       /diag_eddyuu_avg/eddyuu_avg
     &       /diag_eddyvv_avg/eddyvv_avg
     &       /diag_eddyuv_avg/eddyuv_avg
     &       /diag_eddyub_avg/eddyub_avg
     &       /diag_eddyvb_avg/eddyvb_avg
     &       /diag_eddywb_avg/eddywb_avg
     &       /diag_eddyuw_avg/eddyuw_avg
     &       /diag_eddyvw_avg/eddyvw_avg
     &       /diag_eddyubu_avg/eddyubu_avg
     &       /diag_eddyvbv_avg/eddyvbv_avg
     &       /diag_eddyusu_avg/eddyusu_avg
     &       /diag_eddyvsv_avg/eddyvsv_avg
     &       /diag_eddyugsu_avg/eddyugsu_avg
     &       /diag_eddyvgsv_avg/eddyvgsv_avg
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
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
      if (mynode.eq.0) write(stdout,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
      if (numthreads.gt.NPP) then
        if (mynode.eq.0) write(stdout,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds NPP =', NPP
        ierr=ierr+1
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
        if (mynode.eq.0) write(stdout,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1
      endif
      time=0.D0
      tdays=0.D0
      PREDICTOR_2D_STEP=.FALSE.
      iic=0
      kstp=1
      krhs=1
      knew=1
      ntstart=1
      nstp=1
      nrhs=1
      nnew=1
      nfast=1
      synchro_flag=.true.
      first_time=0
      may_day_flag=0
      nrecrst=0
      nrechis=0
      nrecavg=0
      tile_count=0
      bc_count=0
      avgke=0.D0
      avgpe=0.D0
      avgkp=0.D0
      volume=0.D0
      hmin=+1.D+20
      hmax=-1.D+20
      grdmin=+1.D+20
      grdmax=-1.D+20
      Cu_min=+1.D+20
      Cu_max=-1.D+20
      bc_crss=0.D0
      rx0=-1.D+20
      rx1=-1.D+20
      ncidrst=-1
      ncidhis=-1
      ncidavg=-1
      ncidfrc=-1
      ncidbulk=-1
      ncidclm=-1
      ncidqbar=-1
      ncidbtf=-1
      nciddiags_eddy=-1
      nrecdiags_eddy=0
      nciddiags_eddy_avg=-1
      nrecdiags_eddy_avg=0
      vname(:,:)='  '
      vname(1,indxTime)='scrum_time'
      vname(2,indxTime)='time since initialization'
      vname(3,indxTime)='second'
      vname(4,indxTime)='time, scalar, series'
      vname(5,indxTime)='time'
      vname(6,indxTime)='    '
      vname(7,indxTime)='T'
      vname(8,indxTime)='.NO.'
      vname(9,indxTime)=' '
      vname(1,indxTime2)='time'
      vname(2,indxTime2)='time since initialization'
      vname(3,indxTime2)='second'
      vname(4,indxTime2)='time, scalar, series'
      vname(5,indxTime2)='time'
      vname(6,indxTime2)='  '
      vname(7,indxTime2)='T '
      vname(8,indxTime2)='.NO.'
      vname(9,indxTime2)=' '
      vname(1,indxZ)='zeta                                        '
      vname(2,indxZ)='free-surface                                '
      vname(3,indxZ)='meter                                       '
      vname(4,indxZ)='free-surface, scalar, series                '
      vname(5,indxZ)='sea_surface_height                          '
      vname(6,indxZ)='x_rho y_rho                                 '
      vname(7,indxZ)='                                            '
      vname(1,indxUb)='ubar                                       '
      vname(2,indxUb)='vertically integrated u-momentum component '
      vname(3,indxUb)='meter second-1                             '
      vname(4,indxUb)='ubar-velocity, scalar, series              '
      vname(5,indxUb)='barotropic_sea_water_x_'
     &                                  // 'velocity_at_u_location'
      vname(7,indxUb)='                                           '
      vname(1,indxVb)='vbar                                       '
      vname(2,indxVb)='vertically integrated v-momentum component '
      vname(3,indxVb)='meter second-1                             '
      vname(4,indxVb)='vbar-velocity, scalar, series              '
      vname(5,indxVb)='barotropic_sea_water_y_'
     &                                  // 'velocity_at_v_location'
      vname(7,indxVb)='                                           '
      vname(11,indxVb)=' '
      vname(1,indxBostr)='bostr                                   '
      vname(2,indxBostr)='Kinematic bottom stress                 '
      vname(3,indxBostr)='N/m2                                    '
      vname(4,indxBostr)='                                        '
      vname(5,indxBostr)='                                        '
      vname(6,indxBostr)='x_rho y_rho                             '
      vname(7,indxBostr)='                                        '
      vname(1,indxBustr)='bustr                                   '
      vname(2,indxBustr)='Kinematic u bottom stress component     '
      vname(3,indxBustr)='N/m2                                    '
      vname(4,indxBustr)='                                        '
      vname(5,indxBustr)='                                        '
      vname(7,indxBustr)='                                        '
      vname(1,indxBvstr)='bvstr                                   '
      vname(2,indxBvstr)='Kinematic v bottom stress component     '
      vname(3,indxBvstr)='N/m2                                    '
      vname(4,indxBvstr)='                                        '
      vname(5,indxBvstr)='                                        '
      vname(7,indxBvstr)='                                        '
      vname(1,indxWstr)='wstr                                     '
      vname(2,indxWstr)='Kinematic wind stress                    '
      vname(3,indxWstr)='N/m2                                     '
      vname(4,indxWstr)='                                         '
      vname(5,indxWstr)='magnitude_of_surface_downward_stress     '
      vname(6,indxWstr)='x_rho y_rho                              '
      vname(7,indxWstr)='                                         '
      vname(1,indxUWstr)='sustr                                   '
      vname(2,indxUWstr)='Kinematic u wind stress component       '
      vname(3,indxUWstr)='N/m2                                    '
      vname(4,indxUWstr)='                                        '
      vname(5,indxUWstr)='surface_downward_eastward_stress        '
      vname(7,indxUWstr)='                                        '
      vname(8,indxUWstr)='                                        '
      vname(1,indxVWstr)='svstr                                   '
      vname(2,indxVWstr)='Kinematic v wind stress component       '
      vname(3,indxVWstr)='N/m2                                    '
      vname(4,indxVWstr)='                                        '
      vname(5,indxVWstr)='surface_downward_northward_stress       '
      vname(7,indxVWstr)='                                        '
      vname(1,indxU)='u                                           '
      vname(2,indxU)='u-momentum component                        '
      vname(3,indxU)='meter second-1                              '
      vname(4,indxU)='u-velocity, scalar, series                  '
      vname(5,indxU)='sea_water_x_velocity_at_u_location          '
      vname(7,indxU)='                                            '
      vname(1,indxV)='v                                           '
      vname(2,indxV)='v-momentum component                        '
      vname(3,indxV)='meter second-1                              '
      vname(4,indxV)='v-velocity, scalar, series                  '
      vname(5,indxV)='sea_water_y_velocity_at_v_location          '
      vname(7,indxV)='                                            '
      vname(1,indxT)='temp                                        '
      vname(2,indxT)='potential temperature                       '
      vname(3,indxT)='Celsius                                     '
      vname(4,indxT)='temperature, scalar, series                 '
      vname(5,indxT)='sea_water_potential_temperature             '
      vname(6,indxT)='x_rho y_rho                                 '
      vname(7,indxT)='                                            '
      do i=1,ntrc_pas
        if( i < 10) then
           write(nametrc,'(A,I1)') 'tpas0',i
        else
           write(nametrc,'(A,I2)') 'tpas',i
        endif
        vname(1,indxTPAS(i))=TRIM(nametrc)
        vname(2,indxTPAS(i))=TRIM(nametrc)//" passive tracer"
        vname(3,indxTPAS(i))='no unit                             '
        vname(4,indxTPAS(i))='passive tracer, scalar, series      '
        vname(5,indxTPAS(i))='                                   '
        vname(6,indxTPAS(i))='x_rho y_rho                         '
        vname(7,indxTPAS(i))='                                    '
      enddo
      vname(1,indxShflx)='shflux                                  '
      vname(2,indxShflx)='surface net heat flux                   '
      vname(3,indxShflx)='Watts meter-2                           '
      vname(4,indxShflx)='surface heat flux, scalar, series       '
      vname(6,indxShflx)='x_rho y_rho                             '
      vname(7,indxShflx)='                                        '
      vname(1,indxShflx_rsw)='swrad                               '
      vname(2,indxShflx_rsw)='Short-wave surface radiation        '
      vname(3,indxShflx_rsw)='Watts meter-2                       '
      vname(4,indxShflx_rsw)='                                    '
      vname(5,indxShflx_rsw)='                                    '
      vname(6,indxShflx_rsw)='lat_rho lon_rho                     '
      vname(7,indxShflx_rsw)='                                    '
      vname(1,indxO)='omega                                       '
      vname(2,indxO)='S-coordinate vertical momentum component    '
      vname(3,indxO)='meter second-1                              '
      vname(4,indxO)='omega, scalar, series                       '
      vname(5,indxO)='                                            '
      vname(6,indxO)='lat_rho lon_rho                             '
      vname(7,indxO)='                                            '
      vname(1,indxW)='w                                           '
      vname(2,indxW)='vertical momentum component                 '
      vname(3,indxW)='meter second-1                              '
      vname(4,indxW)='w-velocity, scalar, series                  '
      vname(5,indxW)='upward_sea_water_velocity                   '
      vname(6,indxW)='lat_rho lon_rho                             '
      vname(1,indxR)='rho                                         '
      vname(2,indxR)='density anomaly                             '
      vname(3,indxR)='kilogram meter-3                            '
      vname(4,indxR)='density, scalar, series                     '
      vname(5,indxR)='sea_water_sigma_t                           '
      vname(7,indxR)='                                            '
      vname(1,indxbvf)='bvf                                       '
      vname(2,indxbvf)='Brunt Vaisala Frequency                   '
      vname(3,indxbvf)='second-1                                  '
      vname(4,indxbvf)='bvf, scalar, series                       '
      vname(5,indxbvf)='brunt_vaisala_frequency                   '
      vname(7,indxbvf)='                                          '
      vname(1,indxVisc)='visc3d                                   '
      vname(2,indxVisc)='horizontal viscosity coefficient         '
      vname(3,indxVisc)='meter2 second-1                          '
      vname(5,indxVisc)='ocean_momentum_xy_laplacian              '
      vname(4,indxVisc)='visc3d, scalar, series                   '
      vname(6,indxVisc)='lat_rho lon_rho                          '
      vname(7,indxVisc)='                                         '
      vname(1,indxAkv)='AKv                                       '
      vname(2,indxAkv)='vertical viscosity coefficient            '
      vname(3,indxAkv)='meter2 second-1                           '
      vname(4,indxAkv)='AKv, scalar, series                       '
      vname(5,indxAkv)='ocean_vertical_momentum_diffusivity_'
     &                                          // 'at_w_location '
      vname(6,indxAkv)='lat_rho lon_rho                           '
      vname(7,indxAkv)='                                          '
      vname(1,indxAkt)='AKt                                       '
      vname(2,indxAkt)='temperature vertical diffusion coefficient'
      vname(3,indxAkt)='meter2 second-1                           '
      vname(4,indxAkt)='AKt, scalar, series                       '
      vname(5,indxAkt)='ocean_vertical_heat_diffusivity_'
     &                                         //  'at_w_location '
      vname(6,indxAkt)='lat_rho lon_rho                           '
      vname(7,indxAkt)='                                          '
      vname(1,indxHbl)='hbl                                       '
      vname(2,indxHbl)='depth of planetary boundary layer         '
      vname(3,indxHbl)='meter                                     '
      vname(4,indxHbl)='hbl, scalar, series                       '
      vname(5,indxHbl)='ocean_mixed_layer_thickness_defined_'
     &                                       // 'by_mixing_scheme '
      vname(6,indxHbl)='lat_rho lon_rho                           '
      vname(7,indxHbl)='                                          '
      vname(1,indxTke)='tke                                       '
      vname(2,indxTke)='turbulent kinetic energy                  '
      vname(3,indxTke)='meter2 second-2                           '
      vname(4,indxTke)='tke, scalar, series                       '
      vname(5,indxTke)='ocean_turbulent_kinetic_energy'
     &                                         // '_at_w_location '
      vname(6,indxTke)='lat_rho lon_rho                           '
      vname(7,indxTke)='                                          '
      vname(1,indxGls)='gls                                       '
      vname(2,indxGls)='turbulent generic length scale            '
      vname(3,indxGls)='meter3 second-2                           '
      vname(4,indxGls)='gls, scalar, series                       '
      vname(5,indxGls)='ocean_generic_length_scale'
     &                                         // '_at_w_location '
      vname(6,indxGls)='lat_rho lon_rho                           '
      vname(7,indxGls)='                                          '
      vname(1,indxLsc)='Lscale                                    '
      vname(2,indxLsc)='vertical mixing turbulent length scale    '
      vname(3,indxLsc)='meter                                     '
      vname(4,indxLsc)='Lscale, scalar, series                    '
      vname(5,indxLsc)='ocean_turbulent_length_scale'
     &                                         // '_at_w_location '
      vname(6,indxLsc)='lat_rho lon_rho                           '
      vname(7,indxLsc)='                                          '
      vname(1,indxSSH)='SSH                                       '
      vname(2,indxSSH)='sea surface height                        '
      vname(3,indxSSH)='meter                                     '
      vname(4,indxSSH)='SSH, scalar, series                       '
      vname(5,indxSSH)='sea_surface_height_above_sea_level        '
      vname(6,indxSSH)='lat_rho lon_rho                           '
      vname(7,indxSSH)='                                          '
      vname(1,indxSUSTR)='sustr                                   '
      vname(2,indxSUSTR)='surface u-momentum stress               '
      vname(3,indxSUSTR)='Newton meter-2                          '
      vname(4,indxSUSTR)='surface u-mom. stress, scalar, series   '
      vname(5,indxSUSTR)='surface_downward_x_stress               '
      vname(6,indxSUSTR)='lat_u lon_u                             '
      vname(7,indxSUSTR)='                                        '
      vname(1,indxSVSTR)='svstr                                   '
      vname(2,indxSVSTR)='surface v-momentum stress               '
      vname(3,indxSVSTR)='Newton meter-2                          '
      vname(4,indxSVSTR)='surface v-mom. stress, scalar, series   '
      vname(5,indxSVSTR)='surface_downward_y_stress               '
      vname(6,indxSVSTR)='lat_v lon_v                             '
      vname(7,indxSVSTR)='                                        '
       write(vname1,*) 'zz                                  '
       vname(1,indxeddyzz)=vname1
       vname(2,indxeddyzz)='zz                              '
       write(vname3,*) 'meter2                              '
       vname(3,indxeddyzz)=vname3
       vname(4,indxeddyzz)='                                '
       vname(5,indxeddyzz)='                                '
       vname(7,indxeddyzz)='                                '
       write(vname1,*) 'uu                                        '
       vname(1,indxeddyuu)=vname1
       vname(2,indxeddyuu)='uu                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyuu)=vname3
       vname(4,indxeddyuu)='                                      '
       vname(5,indxeddyuu)='                                      '
       vname(7,indxeddyuu)='                                      '
       write(vname1,*) 'vv                                        '
       vname(1,indxeddyvv)=vname1
       vname(2,indxeddyvv)='vv                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyvv)=vname3
       vname(4,indxeddyvv)='                                      '
       vname(5,indxeddyvv)='                                      '
       vname(7,indxeddyvv)='                                      '
       write(vname1,*) 'uv                                        '
       vname(1,indxeddyuv)=vname1
       vname(2,indxeddyuv)='uv                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyuv)=vname3
       vname(4,indxeddyuv)='                                      '
       vname(5,indxeddyuv)='                                      '
       vname(7,indxeddyuv)='                                      '
       write(vname1,*) 'uw                                        '
       vname(1,indxeddyuw)=vname1
       vname(2,indxeddyuw)='uw                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyuw)=vname3
       vname(4,indxeddyuw)='                                      '
       vname(5,indxeddyuw)='                                      '
       vname(7,indxeddyuw)='                                      '
       write(vname1,*) 'vw                                        '
       vname(1,indxeddyvw)=vname1
       vname(2,indxeddyvw)='vw                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyvw)=vname3
       vname(4,indxeddyvw)='                                      '
       vname(5,indxeddyvw)='                                      '
       vname(7,indxeddyvw)='                                      '
       write(vname1,*) 'wb                                        '
       vname(1,indxeddywb)=vname1
       vname(2,indxeddywb)='wb                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddywb)=vname3
       vname(4,indxeddywb)='                                      '
       vname(5,indxeddywb)='                                      '
       vname(7,indxeddywb)='                                      '
       write(vname1,*) 'ub                                        '
       vname(1,indxeddyub)=vname1
       vname(2,indxeddyub)='ub                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyub)=vname3
       vname(4,indxeddyub)='                                      '
       vname(5,indxeddyub)='                                      '
       vname(7,indxeddyub)='                                      '
       write(vname1,*) 'vb                                        '
       vname(1,indxeddyvb)=vname1
       vname(2,indxeddyvb)='vb                                    '
       write(vname3,*) 'meter2 second-2                           '
       vname(3,indxeddyvb)=vname3
       vname(4,indxeddyvb)='                                '
       vname(5,indxeddyvb)='                                '
       vname(7,indxeddyvb)='                                '
       write(vname1,*) 'ubustr                              '
       vname(1,indxeddyubu)=vname1
       vname(2,indxeddyubu)='ubustr                         '
       write(vname3,*) 'meter3 second-3                     '
       vname(3,indxeddyubu)=vname3
       vname(4,indxeddyubu)='                               '
       vname(5,indxeddyubu)='                               '
       vname(7,indxeddyubu)='                               '
       write(vname1,*) 'vbvstr                              '
       vname(1,indxeddyvbv)=vname1
       vname(2,indxeddyvbv)='vbvstr                         '
       write(vname3,*) 'meter3 second-3                     '
       vname(3,indxeddyvbv)=vname3
       vname(4,indxeddyvbv)='                               '
       vname(5,indxeddyvbv)='                               '
       vname(7,indxeddyvbv)='                               '
       write(vname1,*) 'usustr                              '
       vname(1,indxeddyusu)=vname1
       vname(2,indxeddyusu)='usustr                         '
       write(vname3,*) 'meter3 second-3                     '
       vname(3,indxeddyusu)=vname3
       vname(4,indxeddyusu)='                               '
       vname(5,indxeddyusu)='                               '
       vname(7,indxeddyusu)='                               '
       write(vname1,*) 'vsvstr                              '
       vname(1,indxeddyvsv)=vname1
       vname(2,indxeddyvsv)='vsvstr                         '
       write(vname3,*) 'meter3 second-3                     '
       vname(3,indxeddyvsv)=vname3
       vname(4,indxeddyvsv)='                               '
       vname(5,indxeddyvsv)='                               '
       vname(7,indxeddyvsv)='                               '
       write(vname1,*) 'ugsustr                              '
       vname(1,indxeddyugsu)=vname1
       vname(2,indxeddyugsu)='ugsustr                        '
       write(vname3,*) 'meter3 second-3                      '
       vname(3,indxeddyugsu)=vname3
       vname(4,indxeddyugsu)='                               '
       vname(5,indxeddyugsu)='                               '
       vname(7,indxeddyugsu)='                               '
       write(vname1,*) 'vgsvstr                              '
       vname(1,indxeddyvgsv)=vname1
       vname(2,indxeddyvgsv)='vgsvstr                        '
       write(vname3,*) 'meter3 second-3                      '
       vname(3,indxeddyvgsv)=vname3
       vname(4,indxeddyvgsv)='                               '
       vname(5,indxeddyvgsv)='                               '
       vname(7,indxeddyvgsv)='                               '
      return
      end
