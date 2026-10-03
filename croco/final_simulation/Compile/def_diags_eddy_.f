      subroutine def_diags_eddy (ncid, total_rec, ierr)
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
      integer*4 max_opt_size
      parameter (max_opt_size=5000)
      character*5000 Coptions,srcs
      common /strings/ Coptions,srcs
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
      logical create_new_file, res
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
     &      , p2dgrd(3), p3dgrd(4), pw3dgrd(4)
     &      , csize,cmode
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      character*60 text
      if (may_day_flag.ne.0) return
      ierr=0
      lstr=lenstr(diags_eddyname)
      if (nrpfdiags_eddy.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfdiags_eddy))
        call insert_time_index (diags_eddyname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefdiags_eddy
      if (ncid.ne.-1) create_new_file=.false.
 10   if (create_new_file) then
        lstr=lenstr(diags_eddyname)
        cmode = ior(nf_netcdf4, nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        if (mynode.eq.0) write(stdout,*)
     &       'CREATE DIAGS_EDDY HIS/AVG NC4 PARALLEL FILE'
        ierr=nf_create_par(diags_eddyname(1:lstr), cmode,
     &                     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) diags_eddyname(1:lstr)
          may_day_flag=3
          return
        endif
        call put_global_atts (ncid, ierr)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) diags_eddyname(1:lstr)
          may_day_flag=3
          return
        endif
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        p2dgrd(1)=u2dgrd(1)
        p2dgrd(2)=v2dgrd(2)
        p2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        p3dgrd(1)=u2dgrd(1)
        p3dgrd(2)=v2dgrd(2)
        p3dgrd(3)=r2dgrd(3)
        p3dgrd(4)=timedim
        pw3dgrd(1)=u2dgrd(1)
        pw3dgrd(2)=v2dgrd(2)
        pw3dgrd(3)=w3dgrd(3)
        pw3dgrd(4)=timedim
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 diags_eddyTstep)
        ierr=nf_var_par_access(ncid, diags_eddyTstep, nf_collective)
        ierr=nf_put_att_text (ncid, diags_eddyTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, diags_eddyTime)
        ierr=nf_var_par_access(ncid, diags_eddyTime, nf_collective)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime, 'long_name',
     &       lvar, vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, diags_eddyTime, vname(:,indxTime), 
     &                                                                5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, 
     &                                                  diags_eddyTime2)
        ierr=nf_var_par_access(ncid, diags_eddyTime2, nf_collective)
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2, 'long_name',
     &       lvar, vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, diags_eddyTime2, 
     &                                            vname(:,indxTime2), 5,
     &       NF_REAL, ierr)
        itrc=1
        if (wrtdiags_eddy(itrc)) then
           lvar=lenstr(vname(1,indxeddyzz+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                      NF_REAL, 3, r2dgrd, diags_eddyzz(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyzz(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyzz(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyzz+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyzz(itrc), 'units',
     &                          lvar,vname(3,indxeddyzz+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyzz(itrc), 'field',
     &                       lvar, vname(4,indxeddyzz+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyzz(itrc)
     &            ,vname(:,indxeddyzz+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuu(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuu(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuu+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuu(itrc), 'units',
     &                          lvar,vname(3,indxeddyuu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuu(itrc), 'field',
     &                       lvar, vname(4,indxeddyuu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuu(itrc)
     &             ,vname(:,indxeddyuu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvv(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvv(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvv(itrc), 'units',
     &                          lvar,vname(3,indxeddyvv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvv(itrc)
     &              ,vname(:,indxeddyvv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuv(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuv(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuv(itrc), 'units',
     &                          lvar,vname(3,indxeddyuv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuv(itrc), 'field',
     &                       lvar, vname(4,indxeddyuv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuv(itrc)
     &               ,vname(:,indxeddyuv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuw(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuw(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuw(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyuw+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuw(itrc), 'units',
     &                          lvar,vname(3,indxeddyuw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuw(itrc), 'field',
     &                       lvar, vname(4,indxeddyuw+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuw(itrc)
     &          ,vname(:,indxeddyuw+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvw(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvw(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvw(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvw+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvw(itrc), 'units',
     &                          lvar,vname(3,indxeddyvw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvw(itrc), 'field',
     &                       lvar, vname(4,indxeddyvw+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvw(itrc)
     &            ,vname(:,indxeddyvw+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyub+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyub(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyub(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyub(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyub+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyub(itrc), 'units',
     &                          lvar,vname(3,indxeddyub+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyub(itrc), 'field',
     &                       lvar, vname(4,indxeddyub+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyub(itrc)
     &          ,vname(:,indxeddyub+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvb(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvb(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvb(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvb+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvb(itrc), 'units',
     &                          lvar,vname(3,indxeddyvb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvb(itrc), 'field',
     &                       lvar, vname(4,indxeddyvb+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvb(itrc)
     &          ,vname(:,indxeddyvb+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddywb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddywb(itrc))
           ierr=nf_var_par_access(ncid,diags_eddywb(itrc),nf_collective)
           lvar=lenstr(vname(2,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddywb(itrc), 'long_name',
     &                       lvar, vname(2,indxeddywb+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddywb(itrc), 'units',
     &                          lvar,vname(3,indxeddywb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddywb(itrc), 'field',
     &                       lvar, vname(4,indxeddywb+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddywb(itrc)
     &          ,vname(:,indxeddywb+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyubu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, diags_eddyubu(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyubu(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyubu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyubu+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyubu(itrc), 'units',
     &                         lvar,vname(3,indxeddyubu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyubu(itrc), 'field',
     &                       lvar, vname(4,indxeddyubu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyubu(itrc)
     &           ,vname(:,indxeddyubu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, diags_eddyvbv(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvbv(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvbv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvbv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvbv(itrc), 'units',
     &                         lvar,vname(3,indxeddyvbv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvbv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvbv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvbv(itrc)
     &          ,vname(:,indxeddyvbv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyusu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, diags_eddyusu(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyusu(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyusu(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyusu+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyusu(itrc), 'units',
     &                         lvar,vname(3,indxeddyusu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyusu(itrc), 'field',
     &                       lvar, vname(4,indxeddyusu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyusu(itrc)
     &           ,vname(:,indxeddyusu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, diags_eddyvsv(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvsv(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvsv(itrc), 'long_name',
     &                       lvar, vname(2,indxeddyvsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvsv(itrc), 'units',
     &                         lvar,vname(3,indxeddyvsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvsv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvsv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvsv(itrc)
     &        ,vname(:,indxeddyvsv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, diags_eddyugsu(itrc))
          ierr=nf_var_par_access(ncid,diags_eddyugsu(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyugsu(itrc), 
     &                                                      'long_name',
     &                       lvar, vname(2,indxeddyugsu+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyugsu(itrc), 'units',
     &            lvar,vname(3,indxeddyugsu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyugsu(itrc), 'field',
     &                       lvar, vname(4,indxeddyugsu+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diags_eddyugsu(itrc)
     &          ,vname(:,indxeddyugsu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, diags_eddyvgsv(itrc))
          ierr=nf_var_par_access(ncid,diags_eddyvgsv(itrc),
     &                                                    nf_collective)
           lvar=lenstr(vname(2,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvgsv(itrc), 
     &                                                      'long_name',
     &                       lvar, vname(2,indxeddyvgsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvgsv(itrc), 'units',
     &               lvar,vname(3,indxeddyvgsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvgsv(itrc), 'field',
     &                       lvar, vname(4,indxeddyvgsv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diags_eddyvgsv(itrc)
     &          ,vname(:,indxeddyvgsv+itrc-1),5,NF_REAL, ierr)
        endif
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &        'DEF_DIAGS_EDDY - Created ',
     &                'new netCDF file ''',
     &                 diags_eddyname(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open_par (diags_eddyname(1:lstr), IOR(nf_write, 
     &                                                        nf_mpiio),
     &                    MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
        if (ierr .eq. nf_noerr) then
          ierr=checkdims (ncid, diags_eddyname, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfdiags_eddy.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfdiags_eddy))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &           ) 'WARNING: def_diags_eddy: Actual number of records',
     &               rec,  'in netCDF file',  '''',  
     &                                           diags_eddyname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfdiags_eddy.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          create_new_file=.true.
          goto 10
        endif
        ierr=nf_inq_varid (ncid, 'time_step', diags_eddyTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', diags_eddyname(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTstep, nf_collective)
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),
     &                                                   diags_eddyTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), 
     &                                            diags_eddyname(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTime, nf_collective)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),
     &                                                  diags_eddyTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), 
     &                                            diags_eddyname(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTime2, nf_collective)
        itrc=1
        if (wrtdiags_eddy(itrc)) then
          lvar=lenstr(vname(1,indxeddyzz+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                            diags_eddyzz(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyzz(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                            diags_eddyuu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuu(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                            diags_eddyvv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvv(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                            diags_eddyuv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuv(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                            diags_eddyuw(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuw(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                            diags_eddyvw(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvw(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyub+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                            diags_eddyub(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyub(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                            diags_eddyvb(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvb(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddywb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                            diags_eddywb(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddywb(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyubu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                            diags_eddyubu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                                diags_eddyubu(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                            diags_eddyvbv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                                diags_eddyvbv(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyusu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                            diags_eddyusu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                                diags_eddyusu(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                            diags_eddyvsv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                                diags_eddyvsv(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                            diags_eddyugsu(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                               diags_eddyugsu(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                            diags_eddyvgsv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                         diags_eddyname(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                               diags_eddyvgsv(itrc),nf_collective)
        endif
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &    'def_diags_eddy ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'def_diags_eddy ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
  11  format(/' def_diags_eddy - unable to create diag file: ',a)
  20  format(/' def_diags_eddy - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
  99  return
      end
      subroutine def_diags_eddy_avg(ncid, total_rec, ierr)
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
      integer*4 max_opt_size
      parameter (max_opt_size=5000)
      character*5000 Coptions,srcs
      common /strings/ Coptions,srcs
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
      logical create_new_file, res
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
     &      , p2dgrd(3), p3dgrd(4), pw3dgrd(4)
     &      , csize,cmode
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      character*60 text
      if (may_day_flag.ne.0) return
      ierr=0
      lstr=lenstr(diags_eddyname_avg)
      if (nrpfdiags_eddy_avg.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfdiags_eddy_avg))
        call insert_time_index (diags_eddyname_avg, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefdiags_eddy_avg
      if (ncid.ne.-1) create_new_file=.false.
 10   if (create_new_file) then
        lstr=lenstr(diags_eddyname_avg)
        cmode = ior(nf_netcdf4, nf_classic_model)
        cmode = ior(cmode, nf_mpiio)
        csize=xi_rho*eta_rho/NNODES
        if (mynode.eq.0) write(stdout,*)
     &       'CREATE DIAGS_EDDY HIS/AVG NC4 PARALLEL FILE'
        ierr=nf_create_par(diags_eddyname_avg(1:lstr), cmode,
     &                     MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) diags_eddyname_avg(1:lstr)
          may_day_flag=3
          return
        endif
        call put_global_atts (ncid, ierr)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) diags_eddyname_avg(1:lstr)
          may_day_flag=3
          return
        endif
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        p2dgrd(1)=u2dgrd(1)
        p2dgrd(2)=v2dgrd(2)
        p2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
        w3dgrd(1)=r2dgrd(1)
        w3dgrd(2)=r2dgrd(2)
        w3dgrd(4)=timedim
        p3dgrd(1)=u2dgrd(1)
        p3dgrd(2)=v2dgrd(2)
        p3dgrd(3)=r2dgrd(3)
        p3dgrd(4)=timedim
        pw3dgrd(1)=u2dgrd(1)
        pw3dgrd(2)=v2dgrd(2)
        pw3dgrd(3)=w3dgrd(3)
        pw3dgrd(4)=timedim
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 
     &                                              diags_eddyTstep_avg)
        ierr=nf_var_par_access(ncid, diags_eddyTstep_avg, nf_collective)
        ierr=nf_put_att_text (ncid, diags_eddyTstep_avg, 'long_name', 
     &                                                               48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, 
     &                                               diags_eddyTime_avg)
        ierr=nf_var_par_access(ncid, diags_eddyTime_avg, nf_collective)
        text='avg'/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, diags_eddyTime_avg, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime_avg, 'long_name', 
     &                                                             lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime_avg, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, diags_eddyTime_avg, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, diags_eddyTime_avg, 
     &                                             vname(:,indxTime), 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, 
     &                                              diags_eddyTime2_avg)
        ierr=nf_var_par_access(ncid, diags_eddyTime2_avg, nf_collective)
        text='avg'/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, diags_eddyTime2_avg, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2_avg, 'long_name', 
     &                                                             lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2_avg, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, diags_eddyTime2_avg, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, diags_eddyTime2_avg, 
     &                                            vname(:,indxTime2), 5,
     &       NF_REAL, ierr)
        itrc=1
        if (wrtdiags_eddy_avg(itrc)) then
           lvar=lenstr(vname(1,indxeddyzz+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                      NF_REAL, 3, r2dgrd, diags_eddyzz_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyzz_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyzz+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyzz_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyzz_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyzz+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyzz+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyzz_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyzz+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyzz_avg(itrc)
     &            ,vname(:,indxeddyzz+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuu_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuu_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyuu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyuu_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuu_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyuu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuu_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyuu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuu_avg(itrc)
     &             ,vname(:,indxeddyuu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvv_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvv_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvv_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvv_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyvv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvv_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvv_avg(itrc)
     &              ,vname(:,indxeddyvv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuv_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuv_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyuv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyuv_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuv_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyuv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuv_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyuv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuv_avg(itrc)
     &               ,vname(:,indxeddyuv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyuw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyuw_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyuw_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyuw+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyuw_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyuw_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyuw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyuw+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyuw_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyuw+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyuw_avg(itrc)
     &          ,vname(:,indxeddyuw+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvw+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvw_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvw_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvw+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvw_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvw_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyvw+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvw+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvw_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvw+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvw_avg(itrc)
     &            ,vname(:,indxeddyvw+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyub+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyub_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyub_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyub+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyub_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyub_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyub+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyub+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyub_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyub+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyub_avg(itrc)
     &          ,vname(:,indxeddyub+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddyvb_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvb_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvb+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvb_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvb_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddyvb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvb+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvb_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvb+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvb_avg(itrc)
     &          ,vname(:,indxeddyvb+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddywb+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                      NF_REAL, 4, r3dgrd, diags_eddywb_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddywb_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddywb+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddywb_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddywb_avg(itrc), 'units',
     &                          lvar,vname(3,indxeddywb+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddywb+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddywb_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddywb+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddywb_avg(itrc)
     &          ,vname(:,indxeddywb+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyubu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, diags_eddyubu_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyubu_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyubu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyubu_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyubu_avg(itrc), 'units',
     &                         lvar,vname(3,indxeddyubu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyubu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyubu_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyubu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyubu_avg(itrc)
     &           ,vname(:,indxeddyubu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, diags_eddyvbv_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvbv_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvbv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvbv_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvbv_avg(itrc), 'units',
     &                         lvar,vname(3,indxeddyvbv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvbv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvbv_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvbv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvbv_avg(itrc)
     &          ,vname(:,indxeddyvbv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyusu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, diags_eddyusu_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyusu_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyusu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyusu_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyusu_avg(itrc), 'units',
     &                         lvar,vname(3,indxeddyusu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyusu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyusu_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyusu+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyusu_avg(itrc)
     &           ,vname(:,indxeddyusu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, diags_eddyvsv_avg(itrc))
           ierr=nf_var_par_access(ncid,diags_eddyvsv_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvsv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvsv_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvsv_avg(itrc), 'units',
     &                         lvar,vname(3,indxeddyvsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvsv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvsv_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvsv+itrc-1)(1:lvar))
           call nf_add_attribute(ncid,diags_eddyvsv_avg(itrc)
     &        ,vname(:,indxeddyvsv+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                      NF_REAL, 3, u2dgrd, 
     &                                         diags_eddyugsu_avg(itrc))
          ierr=nf_var_par_access(ncid,diags_eddyugsu_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyugsu+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyugsu_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyugsu_avg(itrc), 
     &                                                          'units',
     &            lvar,vname(3,indxeddyugsu+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyugsu+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyugsu_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyugsu+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diags_eddyugsu_avg(itrc)
     &          ,vname(:,indxeddyugsu+itrc-1),5,NF_REAL, ierr)
           lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                      NF_REAL, 3, v2dgrd, 
     &                                         diags_eddyvgsv_avg(itrc))
          ierr=nf_var_par_access(ncid,diags_eddyvgsv_avg(itrc),
     &                                                    nf_collective)
           text='averaged '/ /vname(2,indxeddyvgsv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diags_eddyvgsv_avg(itrc), 
     &                                                      'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid, diags_eddyvgsv_avg(itrc), 
     &                                                          'units',
     &               lvar,vname(3,indxeddyvgsv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxeddyvgsv+itrc-1))
           ierr=nf_put_att_text (ncid,diags_eddyvgsv_avg(itrc), 'field',
     &                       lvar, vname(4,indxeddyvgsv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diags_eddyvgsv_avg(itrc)
     &          ,vname(:,indxeddyvgsv+itrc-1),5,NF_REAL, ierr)
        endif
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &        'DEF_DIAGS_EDDY_AVG - Created ',
     &                'new netCDF file ''',
     &                 diags_eddyname_avg(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open_par (diags_eddyname_avg(1:lstr), IOR(nf_write, 
     &                                                        nf_mpiio),
     &                    MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
        if (ierr .eq. nf_noerr) then
          ierr=checkdims (ncid, diags_eddyname_avg, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfdiags_eddy_avg.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfdiags_eddy_avg))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &           ) 'WARNING: def_diags_eddy: Actual number of records',
     &               rec,  'in netCDF file',  '''',  
     &                                       diags_eddyname_avg(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfdiags_eddy_avg.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          create_new_file=.true.
          goto 10
        endif
        ierr=nf_inq_varid (ncid, 'time_step', diags_eddyTstep_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', diags_eddyname_avg(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTstep_avg, nf_collective)
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),
     &                                               diags_eddyTime_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), 
     &                                        diags_eddyname_avg(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTime_avg, nf_collective)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),
     &                                              diags_eddyTime2_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), 
     &                                        diags_eddyname_avg(1:lstr)
          goto 99
        endif
        ierr=nf_var_par_access(ncid, diags_eddyTime2_avg, nf_collective)
        itrc=1
        if (wrtdiags_eddy_avg(itrc)) then
          lvar=lenstr(vname(1,indxeddyzz+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyzz_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyzz+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyzz_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyuu_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuu+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuu_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyvv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvv+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvv_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyuv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuv+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuv_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyuw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyuw_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyuw+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyuw_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvw+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyvw_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvw+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvw_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyub+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyub_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyub+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyub_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyvb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddyvb_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvb+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddyvb_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddywb+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                            
     &                                           diags_eddywb_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddywb+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, diags_eddywb_avg(itrc), 
     &                                                    nf_collective)
          lvar=lenstr(vname(1,indxeddyubu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                            
     &                                          diags_eddyubu_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyubu+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                            diags_eddyubu_avg(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvbv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                            
     &                                          diags_eddyvbv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvbv+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                            diags_eddyvbv_avg(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyusu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                            
     &                                          diags_eddyusu_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyusu+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                            diags_eddyusu_avg(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                            
     &                                          diags_eddyvsv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvsv+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                            diags_eddyvsv_avg(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyugsu+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                            
     &                                         diags_eddyugsu_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyugsu+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                           diags_eddyugsu_avg(itrc),nf_collective)
          lvar=lenstr(vname(1,indxeddyvgsv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                            
     &                                         diags_eddyvgsv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxeddyvgsv+itrc-1)(1:lvar),
     &                                         
     &                                        diags_eddyname_avg(1:lstr)
            goto 99
          endif
          ierr=nf_var_par_access(ncid, 
     &                           diags_eddyvgsv_avg(itrc),nf_collective)
        endif
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'def_diags_eddy: -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &    'def_diags_eddy ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'def_diags_eddy ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
  11  format(/' def_diags_eddy - unable to create diag file: ',a)
  20  format(/' def_diags_eddy - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
  99  return
      end
