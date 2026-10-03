      module ncscrum
      use param
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
      end module ncscrum
