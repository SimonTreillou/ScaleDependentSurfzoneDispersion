      subroutine read_inp (ierr)
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
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
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
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
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
      integer*4 stafield
      parameter(stafield=6)
      integer*4 indxstaGrd, indxstaTemp, indxstaSalt,
     &        indxstaRho, indxstaVel, indxstaVrt
      parameter (     indxstaGrd=1, indxstaTemp=2,
     &                indxstaSalt=3, indxstaRho=4,  indxstaVel=5,
     &                indxstaVrt=6)
      integer*4 ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
     &      , staVrt
     &      , staX,       staY
     &      , staDepth,   staDen
      logical wrtsta(stafield)
      common/incscrum_sta/
     &        ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
     &      , staVrt
     &      , staX,       staY
     &      , staDepth,   staDen
     &      , wrtsta
      character*80  staname,   staposname
      common /cncscrum_sta/ staname,   staposname
      real sustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real svstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real bvstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real bvstrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_stflx/stflx
      real srflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_srflx/srflx
      integer*4 Nfrq, Ndir
      integer*4 Nfrq0
      parameter (Nfrq0=200, Ndir=31, Nfrq=Nfrq0*Ndir)
      real wf_bry(Nfrq), wk_bry(Nfrq), wa_bry(Nfrq)
      real wd_bry(Nfrq), wa_bry_d(Nfrq), wa_bry_f(Nfrq)
      real wkx_bry(Nfrq), wky_bry(Nfrq)
      real wpha_bry(Nfrq)
      common /wave_maker/ wf_bry, wk_bry, wa_bry
      common /wave_maker_d/ wd_bry, wa_bry_d
      common /wave_maker_f/ wa_bry_f
      common /wave_maker_k/ wkx_bry, wky_bry
      common /wave_maker_pha/ wpha_bry
      real wmaker_amp, wmaker_prd, wmaker_dir
      real wmaker_dsp, wmaker_fsp
      common /wave_maker_par/ wmaker_amp, wmaker_prd, wmaker_dir
      common /wave_maker_par/ wmaker_dsp, wmaker_fsp
      real coswd,sinwd,coswds,sinwds
      common /wave_maker_cos/ coswd,sinwd,coswds,sinwds
      include 'mpif.h'
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      integer*4 kwsize, testunit, input
      parameter (kwsize=32, testunit=40, input=15)
      character end_signal*3, keyword*32, fname*180
      parameter (end_signal='end')
      integer*4 ierr, iargc, is,ie, kwlen, lstr, lenstr
     &                                       , itrc
      logical dumboolean
      fname='TEST_CASES/croco.in.Flashrip'
      if (mynode.eq.0 .and. iargc().GT.0) call getarg(1,fname)
      call MPI_Bcast(fname,64,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
      wrthis(indxTime)=.false.
      wrtavg(indxTime)=.false.
      wrtdiags_eddy(3)=.false.
      wrtdiags_eddy_avg(3)=.false.
      ierr=0
      call setup_kwds (ierr)
      open (input,file=fname,status='old',form='formatted',err=97)
   1  keyword='                                '
      read(input,'(A)',err=1,end=99) keyword
      if (ichar(keyword(1:1)).eq.33) goto 1
      is=1
   2  if (is.eq.kwsize) then
        goto 1
      elseif (keyword(is:is).eq.' ') then
        is=is+1
        goto 2
      endif
      ie=is
   3  if (keyword(ie:ie).eq.':') then
        keyword(ie:ie)=' '
        goto 4
      elseif (keyword(ie:ie).ne.' ' .and. ie.lt.kwsize) then
        ie=ie+1
        goto 3
      endif
      goto 1
   4  kwlen=ie-is
      if (is.gt.1) keyword(1:kwlen)=keyword(is:is+kwlen-1)
      if (keyword(1:kwlen).eq.'title') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,'(A)',err=95) title
      elseif (keyword(1:kwlen).eq.'time_stepping') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntimes,dt,ndtfast, ninfo
        dtfast=dt/float(ndtfast)
        if (NWEIGHT.lt.(2*ndtfast-1)) then
          write(stdout,'(a,i3)')
     &    ' Error - Number of 2D timesteps (2*ndtfast-1): ',
     &    2*ndtfast-1
          write(stdout,'(a,i3)')
     &    '           exceeds barotopic weight dimension: ',NWEIGHT
          goto 95
        endif
      elseif (keyword(1:kwlen).eq.'time_stepping_nbq') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ndtnbq, csound_nbq, visc2_nbq
        dtnbq=dtfast
        ndtnbq=1
      elseif (keyword(1:kwlen).eq.'S-coord') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) theta_s, theta_b, Tcline
      elseif (keyword(1:kwlen).eq.'initial') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrrec
        if (nrrec.gt.0) then
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          open (testunit, file=fname(1:lstr), status='old', err=97)
          close(testunit)
          ininame=fname(1:lstr)
        endif
      elseif (keyword(1:kwlen).eq.'restart') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) nrst, nrpfrst
        read(input,'(A)',err=95)  fname
        lstr=lenstr(fname)
        rstname=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'history') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefhis, nwrt, nrpfhis
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        hisname=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ntsavg, navg, nrpfavg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        avgname=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'diags_eddy') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_eddy, nwrtdiags_eddy, 
     &                                                    nrpfdiags_eddy
        if (nwrtdiags_eddy.eq.0) nwrtdiags_eddy = nwrt
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        diags_eddyname=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'diags_eddy_avg') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) ldefdiags_eddy_avg, ntsdiags_eddy_avg,
     &                      nwrtdiags_eddy_avg,  nrpfdiags_eddy_avg
        if (nwrtdiags_eddy_avg.eq.0) nwrtdiags_eddy_avg = navg
        read(input,'(A)',err=95) fname
        lstr=lenstr(fname)
        diags_eddyname_avg=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'stations') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) ldefsta, nsta, nrpfsta
          read(input,'(A)',err=95) staposname
          read(input,'(A)',err=95) fname
          lstr=lenstr(fname)
          staname=fname(1:lstr)
      elseif (keyword(1:kwlen).eq.'primary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxZ),  wrthis(indxUb)
     &                                       ,  wrthis(indxVb)
     &                    ,  wrthis(indxU),  wrthis(indxV)
        if ( wrthis(indxZ) .or. wrthis(indxUb) .or. wrthis(indxVb)
     &                        .or. wrthis(indxU) .or. wrthis(indxV)
     &     ) wrthis(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'auxiliary_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)
     &                                             wrthis(indxR)
     &                                          ,  wrthis(indxO)
     &                                          ,  wrthis(indxW)
     &                                          ,  wrthis(indxAkv)
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  wrthis(indxbvf)
     &                                          ,  wrthis(indxVisc)
     &                                          ,  dumboolean
     &                                          ,  wrthis(indxHbl)
     &                                          ,  dumboolean
     &                                          ,  wrthis(indxBostr)
     &                                          ,  wrthis(indxBustr)
     &                                          ,  wrthis(indxBvstr)
     &                                          ,  wrthis(indxWstr)
     &                                          ,  wrthis(indxUWstr)
     &                                          ,  wrthis(indxVWstr)
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
        if ( wrthis(indxR)
     &                                        .or. wrthis(indxO)
     &                                        .or. wrthis(indxW)
     &                                        .or. wrthis(indxAkv)
     &                                        .or. wrthis(indxbvf)
     &                                        .or. wrthis(indxVisc)
     &                                        .or. wrthis(indxHbl)
     &                                        .or. wrthis(indxBostr)
     &                                        .or. wrthis(indxBustr)
     &                                        .or. wrthis(indxBvstr)
     &                                        .or. wrthis(indxWstr)
     &                                        .or. wrthis(indxUWstr)
     &                                        .or. wrthis(indxVWstr)
     &     ) wrthis(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'gls_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrthis(indxTke),  wrthis(indxGls)
     &                    ,  wrthis(indxLsc)
        if (wrthis(indxTke) .or. wrthis(indxGls) .or. wrthis(indxLsc)
     &     ) wrthis(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'primary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxZ),  wrtavg(indxUb)
     &                                    ,  wrtavg(indxVb)
     &                    ,  wrtavg(indxU),  wrtavg(indxV)
        if ( wrtavg(indxZ) .or. wrtavg(indxUb) .or. wrtavg(indxVb)
     &                     .or. wrtavg(indxU)  .or. wrtavg(indxV)
     &     ) wrtavg(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'auxiliary_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxR), wrtavg(indxO)
     &        ,  wrtavg(indxW),  wrtavg(indxAkv)
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  wrtavg(indxbvf)
     &                                          ,  wrtavg(indxVisc)
     &                                          ,  dumboolean
     &                                          ,  wrtavg(indxHbl)
     &                                          ,  dumboolean
     &                                          ,  wrtavg(indxBostr)
     &                                          ,  wrtavg(indxBustr)
     &                                          ,  wrtavg(indxBvstr)
     &                                          ,  wrtavg(indxWstr)
     &                                          ,  wrtavg(indxUWstr)
     &                                          ,  wrtavg(indxVWstr)
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          ,  dumboolean
     &                                          , dumboolean
     &                                          , dumboolean
     &                                          ,  dumboolean
        if ( wrtavg(indxR) .or. wrtavg(indxO) .or. wrtavg(indxW)
     &                   .or. wrtavg(indxAkv)
     &                                        .or. wrtavg(indxbvf)
     &                                        .or. wrtavg(indxVisc)
     &                                        .or. wrtavg(indxHbl)
     &                                        .or. wrtavg(indxBostr)
     &                                        .or. wrtavg(indxBustr)
     &                                        .or. wrtavg(indxBvstr)
     &                                        .or. wrtavg(indxWstr)
     &                                        .or. wrtavg(indxUWstr)
     &                                        .or. wrtavg(indxVWstr)
     &     ) wrtavg(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'gls_averages') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxTke),  wrtavg(indxGls)
     &                    ,  wrtavg(indxLsc)
        if ( wrtavg(indxAkk) .or. wrtavg(indxAkp) .or. wrtavg(indxTke)
     &                       .or. wrtavg(indxGls) .or. wrtavg(indxLsc)
     &     ) wrtavg(indxTime)=.true.
      elseif (keyword(1:kwlen).eq.'diags_eddy_history_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95)  wrtdiags_eddy(1)
         do itrc=1,1
           if (wrtdiags_eddy(itrc)) wrtdiags_eddy(3)=.true.
        enddo
      elseif (keyword(1:kwlen).eq.'diags_eddy_average_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtdiags_eddy_avg(1)
        do itrc=1,1
          if (wrtdiags_eddy_avg(itrc)) wrtdiags_eddy_avg(3)=.true.
       enddo
      elseif (keyword(1:kwlen).eq.'station_fields') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) wrtsta(indxstaGrd),
     &       wrtsta(indxstaTemp), wrtsta(indxstaSalt),
     &       wrtsta(indxstaRho), wrtsta(indxstaVel),
     &       wrtsta(indxstaVrt),  wrtsta(indxstaPtr)
      elseif (keyword(1:kwlen).eq.'rho0') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rho0
      elseif (keyword(1:kwlen).eq.'lateral_visc') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) visc2, visc4
      elseif (keyword(1:kwlen).eq.'bottom_drag') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) rdrg, rdrg2, Zobt, Cdb_min, Cdb_max
      elseif (keyword(1:kwlen).eq.'gamma2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) gamma2
      elseif (keyword(1:kwlen).eq.'tracer_diff2') then
        call cancel_kwd (keyword(1:kwlen), ierr)
      elseif (keyword(1:kwlen).eq.'sponge') then
         call cancel_kwd (keyword(1:kwlen), ierr)
      elseif (keyword(1:kwlen).eq.'nudg_cof') then
        call cancel_kwd (keyword(1:kwlen), ierr)
          read(input,*,err=95) tauT_in,tauT_out,tauM_in,tauM_out
          tauT_in =1.D0/(tauT_in *86400.D0)
          tauT_out=1.D0/(tauT_out*86400.D0)
          tauM_in =1.D0/(tauM_in *86400.D0)
          tauM_out=1.D0/(tauM_out*86400.D0)
      elseif (keyword(1:kwlen).eq.'lin_EOS_cff') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) R0, T0, S0, Tcoef, Scoef
      elseif (keyword(1:kwlen).eq.'wave_maker') then
        call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wmaker_amp,wmaker_prd,wmaker_dir,
     &                       wmaker_dsp,wmaker_fsp
      else
        if (mynode.eq.0) write(stdout,'(/3(1x,A)/)')
     &                  'WARNING: Unrecognized keyword:',
     &                   keyword(1:kwlen),' --> DISREGARDED.'
      endif
      if (keyword(1:kwlen) .eq. end_signal) goto 99
      goto 1
  95  write(stdout,'(/1x,4A/)') 'READ_INP ERROR while reading block',
     &                    ' with keyword ''', keyword(1:kwlen), '''.'
      ierr=ierr+1
      goto 99
  97  write(stdout,'(/1x,4A/)') 'READ_INP ERROR: Cannot find input ',
     &                                'file ''', fname(1:lstr), '''.'
      ierr=ierr+1
  99  close (input)
      if (ierr.eq.0) then
        call check_kwds (ierr)
        call check_srcs
        call check_switches1 (ierr)
        call check_switches2 (ierr)
      endif
      if (ierr.eq.0) then
          lstr=lenstr(title)
          if (mynode.eq.0) write(stdout,'(/1x,A)') title(1:lstr)
          if (mynode.eq.0) write(stdout,
     & '(I10,2x,A,1x,A /F10.2,2x,A,2(/I10,2x,A,1x,A)/F10.4,2x,A)')
     &    ntimes,  'ntimes   Total number of timesteps for',
     &                                            '3D equations.',
     &    dt,      'dt       Timestep [sec] for 3D equations',
     &    ndtfast, 'ndtfast  Number of 2D timesteps within each',
     &                                                 '3D step.',
     &    ninfo,   'ninfo    Number of timesteps between',
     &                                     'runtime diagnostics.'
          if (mynode.eq.0) write(stdout,
     &    '(I10,2x,A,1x,A,/F10.2,2x,A,/1pe10.3,2x,A,1x,A/)')
     &   ndtnbq,    'ndtnbq      Number of NBQ timesteps within each',
     &                                                 '2D step.',
     &   csound_nbq,'csound_nbq  Sound wave celerity.',
     &   visc2_nbq, 'visc2_nbq   Second viscosity coefficient for',
     &                                      'compressible fluids.'
          if (mynode.eq.0) write(stdout,
     &                        '(3(1pe10.3,2x,A,1x,A/),32x,A)')
     &    theta_s, 'theta_s  S-coordinate surface control',
     &                                               'parameter.',
     &    theta_b, 'theta_b  S-coordinate bottom control',
     &                                               'parameter.',
     &    Tcline,  'Tcline   S-coordinate surface/bottom layer',
     &  'width used in', 'vertical coordinate stretching, meters.'
          if (nrrec.gt.0) then
          lstr=lenstr(ininame)
          if (mynode.eq.0) write(stdout,'(1x,A,2x,A,4x,A,I3)')
     &     'Initial State File:', ininame(1:lstr), 'Record:',nrrec
          endif
          lstr=lenstr(rstname)
          if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,4x,A,I6,4x,A,I4)')
     &             'Restart File:', rstname(1:lstr),
     &             'nrst =', nrst, 'rec/file: ', nrpfrst
          lstr=lenstr(hisname)
          if (mynode.eq.0) write(stdout,
     &             '(7x,A,2x,A,2x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &       'History File:', hisname(1:lstr),  'Create new:',
     &       ldefhis, 'nwrt =', nwrt, 'rec/file =', nrpfhis
          lstr=lenstr(avgname)
          if (mynode.eq.0) write(stdout,
     &         '(2(I10,2x,A,1x,A/32x,A/),6x,A,2x,A,1x,A,I3)')
     &      ntsavg, 'ntsavg      Starting timestep for the',
     &         'accumulation of output', 'time-averaged data.',
     &      navg,   'navg        Number of timesteps between',
     &     'writing of time-averaged','data into averages file.',
     &     'Averages File:', avgname(1:lstr),
     &     'rec/file =', nrpfavg
          lstr=lenstr(diags_eddyname)
          if (mynode.eq.0) write(stdout,
     &    '(5x,A,2x,A/,32x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &    'eddy Diag File:', diags_eddyname(1:lstr),
     &    'Create new:', ldefdiags_eddy, 'nwrt =', nwrtdiags_eddy,
     &    'rec/file =', nrpfdiags_eddy
          lstr=lenstr(diags_eddyname_avg)
          if (mynode.eq.0) write(stdout,
     &    '(5x,A,2x,A/,32x,A,1x,L1,2x,A,I4,2x,A,I3,/32x,A,I10)')
     &    'eddy AVG Diag File:',diags_eddyname_avg(1:lstr),
     &    'Create new:', ldefdiags_eddy_avg,'nwrt =',nwrtdiags_eddy_avg,
     &    'rec/file =',nrpfdiags_eddy_avg,
     &    'Starting timestep = ',ntsdiags_eddy_avg
          lstr=lenstr(staname)
          if (mynode.eq.0) write(stdout,
     &              '(9x,A,2x,A,2x,A,1x,L1,2x,A,I4,2x,A,I3)')
     &        'Station File:',staname(1:lstr),  'Create new:',
     &        ldefsta, 'nsta =', nsta, 'rec/file =', nrpfsta
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &    , wrthis(indxZ),  'write zeta ', 'free-surface.'
     &    , wrthis(indxUb), 'write UBAR ', '2D U-momentum component.'
     &    , wrthis(indxVb), 'write VBAR ', '2D V-momentum component.'
     &    , wrthis(indxU),  'write U    ', '3D U-momentum component.'
     &    , wrthis(indxV),  'write V    ', '3D V-momentum component.'
          if (mynode.eq.0) write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrthis(indxR),    'write RHO  ', 'Density anomaly.'
     &  , wrthis(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrthis(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrthis(indxAkv),  'write Akv  ', 'Vertical viscosity.'
     &  , wrthis(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
     &  , wrthis(indxVisc),  'write Visc3d', 'Horizontal viscosity.'
     &  , wrthis(indxHbl),  'write Hbl  ',
     &                      'Depth of model boundary layer.'
     &  , wrthis(indxBostr), 'write Bostr', 'Bottom Stress.'
     &  , wrthis(indxBustr), 'write Bustr', 'U-Bottom Stress.'
     &  , wrthis(indxBvstr), 'write Bvstr', 'V-Bottom Stress.'
     &  , wrthis(indxWstr),  'write Wstress', 'Wind Stress.'
     &  , wrthis(indxUWstr), 'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrthis(indxVWstr), 'write V-Wstress comp.', 'V-Wind Stress.'
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in history file: (T/F)'
     &   , wrthis(indxTke), 'write TKE ', 'turbulent kinetic energy.  '
     &   , wrthis(indxGls), 'write GLS ', 'generic length scale.'
     &   , wrthis(indxLsc), 'write Lscale ',
     &                                  'vertical mixing length scale.'
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &  'Fields to be saved in averages file: (T/F)'
     &  , wrtavg(indxZ),  'write zeta ', 'free-surface.'
     &  , wrtavg(indxUb), 'write UBAR ', '2D U-momentum component.'
     &  , wrtavg(indxVb), 'write VBAR ', '2D V-momentum component.'
     &  , wrtavg(indxU),  'write U    ', '3D U-momentum component.'
     &  , wrtavg(indxV),  'write V    ', '3D V-momentum component.'
          if (mynode.eq.0) write(stdout,'(8(/6x,l1,2x,A,1x,A))')
     &    wrtavg(indxR),    'write RHO  ', 'Density anomaly'
     &  , wrtavg(indxO),    'write Omega', 'Omega vertical velocity.'
     &  , wrtavg(indxW),    'write W    ', 'True vertical velocity.'
     &  , wrtavg(indxAkv),  'write Akv  ', 'Vertical viscosity'
     &  , wrtavg(indxbvf),  'write bvf  ',
     &                         'Brunt Vaisala Frequency.'
     &  , wrtavg(indxVisc),'write visc3d', 'Horizontal viscosity'
     &  , wrtavg(indxHbl),  'write Hbl  ',
     &                          'Depth of model boundary layer'
     &  , wrtavg(indxBostr),'write Bostr', 'Bottom Stress.'
     &  , wrtavg(indxBustr),'write Bustr', 'U-Bottom Stress.'
     &  , wrtavg(indxBvstr),'write Bvstr', 'V-Bottom Stress.'
     &  , wrtavg(indxWstr), 'write Wstr', 'Wind Stress.'
     &  , wrtavg(indxUWstr),'write U-Wstress comp.', 'U-Wind Stress.'
     &  , wrtavg(indxVWstr),'write V-Wstress comp.', 'V-Wind Stress.'
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A,1x,A))')
     &    'Fields to be saved in average file: (T/F)'
     &   , wrtavg(indxTke), 'write TKE ', 'turbulent kinetic energy.  '
     &   , wrtavg(indxGls), 'write GLS ', 'generic length scale.'
     &   , wrtavg(indxLsc), 'write Lscale ',
     &                                  'vertical mixing length scale.'
          do itrc=1,1
          if (mynode.eq.0) write(stdout, '(6x,L1,2x,A,2x,A,I2)')
     &         wrtdiags_eddy(itrc),
     &         'write Reynolds stress terms '
          enddo
          do itrc=1,1
          if (mynode.eq.0) write(stdout, '(6x,L1,2x,A,2x,A,I2)')
     &        wrtdiags_eddy_avg(itrc),
     &         'write Reynolds stress terms '
          enddo
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/6x,l1,2x,A))')
     &      'Fields to be saved in stations output  (T/F)'
     &     , wrtsta(indxstaGrd),   'write Grid location variables'
     &     , wrtsta(indxstaTemp),  'write temperature.'
     &     , wrtsta(indxstaSalt),  'write salinity.'
     &     , wrtsta(indxstaRho),   'write density.'
     &     , wrtsta(indxstaVel),   'write mean station velocity'
     &     , wrtsta(indxstaVrt),   'write vorticity'
     &     , wrtsta(indxstaPtr),   'write tracer concentrations'
          if (mynode.eq.0) write(stdout,'(F10.4,2x,A,1x,A)')
     &        rho0, 'rho0     Boussinesq approximation',
     &                           'mean density, kg/m3.'
          if (mynode.eq.0) write(stdout,9) visc2
 9    format(1pe10.3,2x,'visc2    Horizontal Laplacian ',
     &       'mixing coefficient [m2/s]',/,32x,'for momentum.')
          if (mynode.eq.0) write(stdout,'(5(1pe10.3,2x,A/))')
     &     rdrg, 'rdrg     Linear bottom drag coefficient (m/si).',
     &    rdrg2, 'rdrg2    Quadratic bottom drag coefficient.',
     &     Zobt, 'Zobt     Bottom roughness for logarithmic law (m).',
     &  Cdb_min, 'Cdb_min  Minimum bottom drag coefficient.',
     &  Cdb_max, 'Cdb_max  Maximum bottom drag coefficient.'
          if (mynode.eq.0) write(stdout,'(f10.2,2x,A,1x,A)')
     &     gamma2, 'gamma2   Slipperiness parameter:',
     &                     'free-slip +1, or no-slip -1.'
          if (mynode.eq.0) write(stdout,'(/,1x,A,/,25x,A/)')
     &   'SPONGE_GRID is defined: x_sponge parameter in sponge/nudging',
     &   'layer is set generically in set_nudgcof.F routine'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauT_in,'tauT_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &       tauT_out,'tauT_out Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A)')
     &        tauM_in,'tauM_in  Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(1pe10.3,2x,A/)')
     &       tauM_out,'tauM_out Nudging coefficients [sec^-1]'
          if (mynode.eq.0) write(stdout,'(5(f10.4,2x,A,1x,A/))')
     &       T0, 'T0       Background value for potential',
     &                                     'temperature (Celsius).',
     &       S0, 'S0       Background salinity (PSU),', 'constant.',
     &       R0, 'R0       Background density (kg/m3) used in',
     &                                                'linear EOS.',
     &    Tcoef, 'Tcoef    Thermal expansion coefficient',
     &                                           '(kg/m3/Celsius).',
     &    Scoef, 'Scoef    Saline contraction coefficient',
     &                                                '(kg/m3/PSU).'
          if (mynode.eq.0) write(stdout,'(/1x,A,5(/3x,f10.4,2x,A))')
     &    'Offshore wavemaker parameters for wave modeling',
     &   wmaker_amp, 'wmaker_amp  RMS wave amplitude [m]',
     &   wmaker_prd, 'wmaker_prd  peak wave period [s]',
     &   wmaker_dir, 'wmaker_dir  mean wave angle [deg]',
     &   wmaker_dsp, 'wmaker_dsp  directional spread [deg]',
     &   wmaker_fsp, 'wmaker_fsp  freq. spread (gamma in JONSWAP)'
      endif
      if (ierr.ne.0) then
        write(stdout,'(/1x,2A,I3,1x,A/)') 'READ_INP ERROR: ',
     & 'A total of', ierr, 'configuration errors discovered.'
        return
      endif
      return
      end
