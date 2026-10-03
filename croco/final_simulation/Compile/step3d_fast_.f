      subroutine step3d_fast (tile)
      implicit none
      integer*4 tile, trd
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
C$    integer*4 omp_get_thread_num
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=0
C$    trd=omp_get_thread_num()
      call step3d_fast_tile ( Istr,Iend,Jstr,Jend
     &                   , A2d(1, 1,trd), A2d(1, 2,trd), A2d(1, 3,trd)
     &                   , A2d(1, 4,trd), A2d(1, 5,trd), A2d(1, 6,trd)
     &                   , A2d(1, 7,trd), A2d(1, 8,trd)
     &                   , A2d(1, 9,trd), A2d(1,10,trd)
     &                   , A2d(1,11,trd), A2d(1,12,trd)
     &                   , A2d(1,13,trd), A2d(1,14,trd)
     &                   , A2d(1,15,trd), A2d(1,16,trd), A2d(1,17,trd)
     &                   , A2d(1,18,trd), A2d(1,20,trd), A2d(1,22,trd)
     &                   , A2d(1,24,trd), A2d(1,25,trd), A2d(1,26,trd)
     &                   , A2d(1,27,trd), A2d(1,28,trd), A2d(1,29,trd)
     &                   , A3d(1, 1,trd), A3d(1, 2,trd)
     &                   , A3d(1, 3,trd), A3d(1, 4,trd)
     &                   , A3d(1, 5,trd), A3d(1, 6,trd)
     &                   , A3d(1, 7,trd), A3d(1, 8,trd)
     &                   )
      return
      end
      subroutine step3d_fast_tile (Istr,Iend,Jstr,Jend
     &                          ,Dnew,rubar,rvbar
     &                          ,Drhs, UFx,UFe
     &                          ,VFx,VFe
     &                          ,urhs,vrhs
     &                          ,DUon,DVom
     &                          ,ru_ext_nbq_sum, rv_ext_nbq_sum
     &                          ,ru_ext_nbq_old, rv_ext_nbq_old,work
     &                          ,dthetadiv_nbqdz_w
     &                          ,dZdxq_u,dZdyq_v
     &                          ,FX,FY,FC,DC,CF,BC
     &                          ,Hzw_half_nbq_inv,  Hzr_half_nbq_inv
     &                          ,Hzw_half_nbq_inv_u,Hzw_half_nbq_inv_v
     &                          ,Hzu_half_qdmu,     Hzv_half_qdmv
     &                          ,FX_3D, FY_3D
     &                          )
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
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, kbak,kold,
     &        imin,imax,jmin,jmax,
     &        k1, k2, kp1
     &       ,err
      integer*4,dimension(2) :: ijmax
      real    mybeta,myalpha,myepsilon,mygamma,
     &        VMAX,VMAXL, cff,cff0,cff1,cff2,cff3,
     &        DUnew,DVnew, dum_s
      real,parameter :: thetaimp_nbq = 1.D0
      real,parameter :: gammau=1.0D0
      real,parameter :: gammau_2=(1.D0/3.D0)*(1.D0-gammau)
      real
     &                Dnew(Istr-2:Iend+2,Jstr-2:Jend+2),
     &               rubar(Istr-2:Iend+2,Jstr-2:Jend+2),
     &               rvbar(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                Drhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                 UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                 UFe(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                 VFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                 VFe(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                urhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                vrhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                DUon(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                DVom(Istr-2:Iend+2,Jstr-2:Jend+2)
      real
     &       ru_ext_nbq_sum(Istr-2:Iend+2,Jstr-2:Jend+2),
     &       rv_ext_nbq_sum(Istr-2:Iend+2,Jstr-2:Jend+2),
     &       ru_ext_nbq_old(Istr-2:Iend+2,Jstr-2:Jend+2),
     &       rv_ext_nbq_old(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                 work(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    Hzu_half_qdmu(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &        Hzv_half_qdmv(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     Hzw_half_nbq_inv(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     Hzr_half_nbq_inv(Istr-2:Iend+2,Jstr-2:Jend+2,  N),
     &   Hzw_half_nbq_inv_u(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &   Hzw_half_nbq_inv_v(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
      real
     &    dthetadiv_nbqdz_w(Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &              dZdxq_u(Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &              dZdyq_v(Istr-2:Iend+2,Jstr-2:Jend+2,2),
     &                   FX(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                   FY(Istr-2:Iend+2,Jstr-2:Jend+2),
     &                   FC(Istr-2:Iend+2,0:N),
     &                   DC(Istr-2:Iend+2,0:N),
     &                   CF(Istr-2:Iend+2,0:N),
     &                   BC(Istr-2:Iend+2,0:N),
     &         FX_3D(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &         FY_3D(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
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
      real rhoA(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rhoS(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rvfrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rufrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real rvfrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DU_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DV_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DU_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DV_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
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
      real btflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_btflx/btflx
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
      real sum_nbq
      real cff1_cst
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
      real cff1_WD,cff2_WD
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
      include 'mpif.h'
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
!$acc kernels if(compute_on_device) default(present)
      if (iic.eq.ntstart.and.iif.eq.1) then
         do k=1,N
           do j=Jstr-2,Jend+1
             do i=IstrU-2,Iend+1
               rho_bak(i,j,k)=rho(i,j,k)
             enddo
            enddo
          enddo
      endif
      do j=JstrR,JendR
          do k=1,N
            do i=IstrR,IendR
              rho_grd(i,j,k)=(rho(i,j,k)*1.5D0-rho_bak(i,j,k)*0.5D0)
     &                                                             /rho0
            enddo
          enddo
      enddo
!$acc end kernels
      mybeta=0.281105D0
       if (iif.eq.1) then
        kbak=kstp
        kold=kstp
        cff1= 1.0D0
        cff2= 0.0D0
        cff3= 0.0D0
       elseif (iif.eq.1+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1= 1.5D0
        cff2=-0.5D0
        cff3= 0.0D0
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.5D0+mybeta
        cff2=-2.0D0*mybeta-0.5D0
        cff3= mybeta
      endif
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr-2,Jend+1
        do i=IstrU-2,Iend+1
          Drhs(i,j)=cff1*(zeta(i,j,kstp)+h(i,j))*rhobar_nbq(i,j,kstp)
     &             +cff2*(zeta(i,j,kbak)+h(i,j))*rhobar_nbq(i,j,kbak)
     &             +cff3*(zeta(i,j,kold)+h(i,j))*rhobar_nbq(i,j,kold)
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=IstrU-1,Iend+1
          urhs(i,j)=cff1*ubar(i,j,kstp)
     &             +cff2*ubar(i,j,kbak)
     &             +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5D0*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
     &                                                              )
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp)
     &             +cff2*vbar(i,j,kbak)
     &             +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5D0*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
     &                                                              )
        enddo
      enddo
!$acc end kernels
      if (iif.eq.1) then
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv)
      endif
!$acc kernels if(compute_on_device) default(present)
      if (IstrU.le.Iend) then
        do j=Jstr,Jend
          do i=IstrU-1,Iend+1
            usurf_nbq(i,j)=(qdmu_nbq(i,j,N)
     &                     *2.D0/(Hz(i,j,N)+Hz(i-1,j,N))
     &                     )
          enddo
        enddo
      endif
      if (Jstr.le.Jend) then
        do j=Jstr-1,Jend+1
          do i=Istr,Iend
            vsurf_nbq(i,j)=(qdmv_nbq(i,j,N)
     &                     *2.D0/(Hz(i,j,N)+Hz(i,j-1,N))
     &                      )
          enddo
        enddo
      endif
      do j=Jstr-1,Jend
        do i=IstrU-1,Iend
          wsurf_nbq(i,j)=(qdmw_nbq(i,j,N)
     &                       /(1.D0+rho_grd(i,j,N))
     &                       *Hzw_half_nbq_inv(i,j,N)
     &                    )
        enddo
      enddo
      if (iif.eq.1) then
        cff1 = 1.0D0
        cff2 = 0.0D0
        cff3 = 0.0D0
      elseif (iif.eq.1+1) then
        cff1 = 1.5D0
        cff2 =-0.5D0
        cff3 = 0.0D0
      else
        cff1 = 1.5D0+mybeta
        cff2 =-2.0D0*mybeta-0.5D0
        cff3 = mybeta
      endif
      do j=Jstr-2,Jend+1
        do i=IstrU-2,Iend+1
          UFx(i,j) =   cff1 * zeta(i,j,kstp)
     &                + cff2 * zeta(i,j,kbak)
     &                + cff3 * zeta(i,j,kold)
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=IstrU-1,Iend
          zeta(i,j,knew)=(zeta(i,j,kstp) + dtfast*( wsurf_nbq(i,j)
     &                     -0.5D0*(usurf_nbq(i  ,j)
     &                                     *(UFx(i  ,j)
     &                                      -UFx(i-1,j))*pm_u(i,j)
     &                          +usurf_nbq(i+1,j)
     &                                     *(UFx(i+1,j)
     &                                      -UFx(i  ,j))*pm_u(i+1,j)
     &                          )
     &                     -0.5D0*(vsurf_nbq(i  ,j)
     &                                     *(UFx(i,j  )
     &                                      -UFx(i,j-1))*pn_v(i,j)
     &                          +vsurf_nbq(i,j+1)
     &                                     *(UFx(i,j+1)
     &                                      -UFx(i,j  ))*pn_v(i,j+1)
     &                          )
     &                                                ) )
        enddo
      enddo
!$acc end kernels
      call wetdry_tile (Istr,Iend,Jstr,Jend)
      call zetabc_tile (Istr,Iend,Jstr,Jend)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(-2,-2,knew))
      if ( (iic.eq.ntstart .and. iif.eq.1) .or. mod(iif,1).eq.0 ) then
        call set_depth_tile(Istr,Iend,Jstr,Jend)
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv )
      endif
!$acc kernels if(compute_on_device) default(present)
      myalpha   = 0.25D0
      myepsilon = 0.00976186D0 - 0.13451357D0*myalpha
      mygamma   = 0.08344500D0 - 0.51358400D0*myalpha
      if (iif.eq.1) then
        cff0=0.D0
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iif.eq.1+1) then
        cff0= 1.0833333333333D0
        cff1=-0.1666666666666D0
        cff2= 0.0833333333333D0
        cff3= 0.D0
      else
        cff0=0.5D0+2.D0*myepsilon+mygamma+2.D0*myalpha
        cff1=1.D0-cff0-mygamma-myepsilon
        cff2=mygamma
        cff3=myepsilon
      endif
      do j=Jstr-1,Jend
        do i=IstrU-1,Iend
          UFx(i,j)=(cff0*zeta(i,j,knew)+cff1*zeta(i,j,kstp)
     &              +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold))
     &                                *rhobar_nbq(i,j,kstp)
          UFe(i,j)=(1.D0+rhoS(i,j))*UFx(i,j)
          VFe(i,j)=UFe(i,j)*UFx(i,j)
          VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
        enddo
      enddo
      cff=0.5D0*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*( (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(UFx(i-1,j)-UFx(i,j)))
     &                                                              )
          rvbar(i,j)=cff*om_v(i,j)*( (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(UFx(i,j-1)-UFx(i,j)))
     &                                                              )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25D0*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))
          VFx(i+1,j)=0.25D0*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25D0*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))
          UFe(i,j+1)=0.25D0*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)
          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo
      if (rdrg2.gt.0.D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25D0*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            cff=0.25D0*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
      if (iif.eq.1) then
        if (iic.eq.ntstart) then
          cff3=0.D0
          cff2=0.D0
          cff1=1.D0
        elseif (iic.eq.ntstart+1) then
          cff3= 0.D0
          cff2=-0.5D0
          cff1= 1.5D0
        else
          cff3= 0.281105D0
          cff2=-0.5D0-2.D0*cff3
          cff1= 1.5D0+cff3
        endif
        do j=Jstr,Jend
         do i=IstrU,Iend
           cff=rufrc(i,j)-rubar(i,j)
           rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                         + cff3*rufrc_bak(i,j,nstp)
           rufrc_bak(i,j,nstp)=cff
         enddo
        enddo
        do j=Jstr,Jend
         do i=Istr,Iend
           cff=rvfrc(i,j)-rvbar(i,j)
           rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                         + cff3*rvfrc_bak(i,j,nstp)
           rvfrc_bak(i,j,nstp)=cff
         enddo
        enddo
      endif
      if (iif.eq.1) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_ext_nbq_sum(i,j)=0.D0
            ru_ext_nbq_old(i,j)=0.D0
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            rv_ext_nbq_sum(i,j)=0.D0
            rv_ext_nbq_old(i,j)=0.D0
          enddo
        enddo
      endif
      do j=Jstr,Jend
        do i=IstrU,Iend
          UFx(i,j)=(rufrc(i,j)+rubar(i,j))
     &                                   *pm_u(i,j)*pn_u(i,j)
     &                                  /(Drhs(i,j)+Drhs(i-1,j))
          ru_ext_nbq_old(i,j)=UFx(i,j)-ru_ext_nbq_old(i,j)
          ru_ext_nbq_sum(i,j)=ru_ext_nbq_sum(i,j)+UFx(i,j)
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru_int_nbq(i,j,k)=ru_int_nbq(i,j,k)+
     &                        ru_ext_nbq_old(i,j)
     &                        *(Hz(i-1,j,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          ru_ext_nbq_old(i,j)=UFx(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          UFx(i,j)=(rvfrc(i,j)+rvbar(i,j))
     &                                   *pm_v(i,j)*pn_v(i,j)
     &                                   /(Drhs(i,j)+Drhs(i,j-1))
          rv_ext_nbq_old(i,j)=UFx(i,j)-rv_ext_nbq_old(i,j)
          rv_ext_nbq_sum(i,j)=rv_ext_nbq_sum(i,j)+UFx(i,j)
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
            rv_int_nbq(i,j,k)=rv_int_nbq(i,j,k)+
     &                        rv_ext_nbq_old(i,j)
     &                        *(Hz(i,j-1,k)+Hz(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rv_ext_nbq_old(i,j)=UFx(i,j)
        enddo
      enddo
      if (iif.eq.1) then
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_nbq_avg2(i,j,k)=qdmu_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rv_nbq_avg2(i,j,k)=qdmv_nbq(i,j,k)
            enddo
          enddo
        enddo
        do k=0,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rw_nbq_avg2(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (iic.eq.ntstart .and. iif.eq.1) then
        do j=JstrR,JendR
          do k=0,N
            do i=IstrR,IendR
              zw_nbq (i,j,k,:)=z_w(i,j,k)
            enddo
          enddo
        enddo
      endif
      do j=Jstr,Jend
        do i=Istr,Iend
          work(i,j)=pm(i,j)*pn(i,j)
          DU_nbq(i,j)=0.D0
          DV_nbq(i,j)=0.D0
        enddo
      enddo
      if (iif.eq.nfast) then
        do k=0,N
          do j=Jstr,Jend
            do i=Istr,Iend
               rw_nbq(i,j,k)=qdmw_nbq(i,j,k)
            enddo
          enddo
        enddo
      endif
      if (iic.gt.ntstart .and. iif.eq.1) then
        do k=1,N
          do j=Jstr-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz(i,j,k) = Hz_correct(i,j,k)
            enddo
          enddo
        enddo
      endif
      call nbq_bry_store_tile (Istr,Iend,Jstr,Jend)
      do k=1,N
        do j=Jstr-2,Jend+1
          do i=IstrU-2,Iend+1
            thetadiv_nbq(i,j,k)=( -visc2_nbq*( thetadiv_nbq(i,j,k)
     &                                       +thetadiv3_nbq(i,j,k))
     &                         +soundspeed2_nbq(i,j)*rho_nbq(i,j,k) )
     &                                     *Hzr_half_nbq_inv(i,j,k)
          enddo
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr-1,Iend
          dthetadiv_nbqdz(i,j,N)    = - thetadiv_nbq(i,j,N)
        enddo
      enddo
      do k=1,N-1
        do j=Jstr-1,Jend
          do i=Istr-1,Iend
                  dthetadiv_nbqdz(i,j,k)   = thetadiv_nbq(i,j,k+1)
     &                                   - thetadiv_nbq(i  ,j,k)
             enddo
          enddo
      enddo
      do j=Jstr,Jend
        do i=Istr-1,Iend
          FY_3D(i,j,0)=0.D0
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend
                FY_3D(i,j,k)=Hzw_half_nbq_inv_u(i,j,k)*(
     &                                    dthetadiv_nbqdz(i,j,k)
     &                                   +dthetadiv_nbqdz(i-1,j,k))
          enddo
        enddo
      enddo
      do k=1,N
        do j=Jstr,Jend
          do i=IstrU,Iend
            if (k.gt.1.and.k.lt.N) then
              dum_s=(z_r(i,j,k)-z_r(i-1,j,k))
     &              *(FY_3D(i,j,k)+
     &                FY_3D(i,j,k-1))
          elseif (k.gt.1) then
            dum_s=(z_r(i,j,k)-z_r(i-1,j,k))
     &                *FY_3D(i,j,k-1)
     &                +(z_w(i,j,N)-z_w(i-1,j,N))
     &                *FY_3D(i,j,k)
          else
            dum_s=(z_r(i,j,k)-z_r(i-1,j,k))
     &              *2.D0*FY_3D(i,j,k)
          endif
            dum_s=dum_s
     &                -(gammau  *thetadiv_nbq(i  ,j,k)+
     &                  gammau_2*thetadiv_nbq(i+1,j,k)-
     &                  gammau  *thetadiv_nbq(i-1,j,k)-
     &                  gammau_2*thetadiv_nbq(i-2,j,k))
          dum_s=dum_s*Hzu_half_qdmu(i,j,k)
          qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k) + dtnbq * (
     &                        dum_s + ru_int_nbq(i,j,k) )
          ru_nbq(i,j,k)=dum_s/work(i,j)
          qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*(1.D0-NBQnudgcof(i,j))
     &                        +u(i,j,k,nrhs)*Hzu_half_qdmu(i,j,k)
     &                                           *NBQnudgcof(i,j)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          sum_nbq=0.D0
!$acc loop reduction(+:sum_nbq)
          do k=1,N
             sum_nbq=sum_nbq+qdmu_nbq(i,j,k)
          enddo
          DU_nbq(i,j)=DU_nbq(i,j)+sum_nbq
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          FY_3D(i,j,0)=0.D0
        enddo
      enddo
      do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
                 FY_3D(i,j,k)=Hzw_half_nbq_inv_v(i,j,k)*(
     &                                     dthetadiv_nbqdz(i,j,k)
     &                                    +dthetadiv_nbqdz(i,j-1,k))
            enddo
          enddo
      enddo
        do k=1,N
          do j=Jstr,Jend
            if (j.ge.Jstr) then
            do i=Istr,Iend
              if (k.gt.1.and.k.lt.N) then
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))
     &                *(FY_3D(i,j,k)+
     &                  FY_3D(i,j,k-1))
              elseif (k.gt.1) then
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))
     &                *FY_3D(i,j,k-1)
     &                +(z_w(i,j,N)-z_w(i,j-1,N))
     &                *FY_3D(i,j,k)
              else
                dum_s=(z_r(i,j,k)-z_r(i,j-1,k))
     &                *2.D0*FY_3D(i,j,k)
              endif
              dum_s=dum_s
     &             -(gammau  *thetadiv_nbq(i,j  ,k)+
     &               gammau_2*thetadiv_nbq(i,j+1,k)-
     &               gammau  *thetadiv_nbq(i,j-1,k)-
     &               gammau_2*thetadiv_nbq(i,j-2,k))
              dum_s=dum_s*Hzv_half_qdmv(i,j,k)
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k) + dtnbq * (
     &                        dum_s + rv_int_nbq(i,j,k) )
              rv_nbq(i,j,k)=dum_s/work(i,j)
              qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*(1.D0-NBQnudgcof(i,j))
     &                        +v(i,j,k,nrhs)*Hzv_half_qdmv(i,j,k)
     &                                           *NBQnudgcof(i,j)
            enddo
           endif
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            sum_nbq=0.D0
!$acc loop reduction(+:sum_nbq)
            do k=1,N
               sum_nbq=sum_nbq+qdmv_nbq(i,j,k)
            enddo
            DV_nbq(i,j)=DV_nbq(i,j)+sum_nbq
          enddo
        enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1_WD=ABS(ABS(umask_wet(i,j))-1.D0)
          cff2_WD=0.5D0+SIGN(0.5D0,DU_nbq(i,j))*umask_wet(i,j)
          umask_wet(i,j)=0.5D0*umask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.D0-cff1_WD)
          DU_nbq(i,j)=DU_nbq(i,j)*umask_wet(i,j)
          do k=1,N
            qdmu_nbq(i,j,k)=qdmu_nbq(i,j,k)*umask_wet(i,j)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          cff1_WD=ABS(ABS(vmask_wet(i,j))-1.D0)
          cff2_WD=0.5D0+SIGN(0.5D0,DV_nbq(i,j))*vmask_wet(i,j)
          vmask_wet(i,j)=0.5D0*vmask_wet(i,j)*cff1_WD
     &                         +cff2_WD*(1.D0-cff1_WD)
          DV_nbq(i,j)=DV_nbq(i,j)*vmask_wet(i,j)
          do k=1,N
            qdmv_nbq(i,j,k)=qdmv_nbq(i,j,k)*vmask_wet(i,j)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          ubar(i,j,knew)=urhs(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          vbar(i,j,knew)=vrhs(i,j)
        enddo
      enddo
!$acc end kernels
      M2bc_nbq_flag=.true.
      call u2dbc_tile   (Istr,Iend,Jstr,Jend, work)
      call v2dbc_tile   (Istr,Iend,Jstr,Jend, work)
      call unbq_bc_tile (Istr,Iend,Jstr,Jend, work)
      call vnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmu_nbq(-2,-2,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        qdmv_nbq(-2,-2,1))
!$acc kernels if(compute_on_device) default(present)
      if (IstrU.gt.Iend) then
        do j=Jstr,Jend
          do i=Istr,Iend+1
            FX(i,j)=0.D0
          enddo
        enddo
      endif
      if (Jstr.gt.Jend) then
        do j=Jstr,Jend+1
          do i=Istr,Iend
            FY(i,j)=0.D0
          enddo
        enddo
      endif
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend+1
            dthetadiv_nbqdz(i,j,k)=(z_r(i,j,k)-z_r(i-1,j,k))
     &                           *qdmu_nbq(i,j,k)
          enddo
        enddo
      enddo
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_w(i,j,0)=0.D0
                qdmw_nbq(i,j,0)=0.D0
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_w(i,j,N )= (z_w(i,j,N)-z_w(i-1,j,N))
     &                            *qdmu_nbq(i,j,N)
     &                            /(Hzr(i,j,N)+Hzr(i-1,j,N))
              enddo
            enddo
      do k=1,N-1
            do j=Jstr,Jend
              do i=Istr,Iend+1
                dZdxq_w(i,j,k )=Hzw_half_nbq_inv_u(i,j,k)
     &                                 *(dthetadiv_nbqdz(i,j,k)+
     &                                   dthetadiv_nbqdz(i,j,k+1))
              enddo
            enddo
      enddo
      do k=1,N
          if (IstrU.le.Iend) then
            do j=Jstr,Jend
              do i=Istr,Iend+1
                FY_3D(i,j,k)=-pm_u(i,j)*
     &       (dZdxq_w(i,j,k )-dZdxq_w(i,j,k-1))
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=
     &  FY_3D(i,j,k)+FY_3D(i+1,j,k)
              enddo
            enddo
          else
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=0.D0
              enddo
            enddo
          endif
      enddo
      do k=1,N
        do j=Jstr,Jend+1
          do i=Istr,Iend
            dthetadiv_nbqdz(i,j,k)=(z_r(i,j,k)-z_r(i,j-1,k))
     &                       *qdmv_nbq(i,j,k)
          enddo
        enddo
      enddo
            do j=Jstr,Jend +1
              do i=Istr,Iend
                dZdyq_w(i,j,0)=0.D0
                qdmw_nbq(i,j,0)=0.D0
              enddo
            enddo
        do j=Jstr,Jend+1
          do i=Istr,Iend
                dZdyq_w(i,j,N )= (z_w(i,j,N)-z_w(i,j-1,N))
     &                           *qdmv_nbq(i,j,N)
     &                           /(Hzr(i,j,N)+Hzr(i,j-1,N))
          enddo
        enddo
      do k=1,N-1
        do j=Jstr,Jend+1
          do i=Istr,Iend
                dZdyq_w(i,j,k )=Hzw_half_nbq_inv_v(i,j,k)
     &                           *(dthetadiv_nbqdz(i,j,k)+
     &                             dthetadiv_nbqdz(i,j,k+1))
          enddo
        enddo
      enddo
          if (Jstr.le.Jend) then
          do k=1,N
            do j=Jstr,Jend+1
              do i=Istr,Iend
                FY_3D(i,j,k)=-pn_v(i,j)*
     &   (dZdyq_w(i,j,k)-dZdyq_w(i,j,k-1))
              enddo
            enddo
            do j=Jstr,Jend
              do i=Istr,Iend
                thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &    +FY_3D(i,j,k)+FY_3D(i,j+1,k)
              enddo
            enddo
          enddo
          endif
      do k=1,N
        if (IstrU.le.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend+1
              FX(i,j)=on_u(i,j)*qdmu_nbq(i,j,k)
            enddo
          enddo
        endif
        if (Jstr.le.Jend) then
          do j=Jstr,Jend+1
            do i=Istr,Iend
              FY(i,j)=om_v(i,j)*qdmv_nbq(i,j,k)
            enddo
          enddo
        endif
        if (IstrU.gt.Iend) then
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                           +pm(i,j)*pn(i,j)*(FY(i,j+1)-FY(i,j))
            enddo
          enddo
        elseif (Jstr.gt.Jend) then
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j))
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                            +pm(i,j)*pn(i,j)*(FX(i+1,j)-FX(i,j)+
     &                                              FY(i,j+1)-FY(i,j))
            enddo
          enddo
        endif
      enddo
        do k=0,N
          do j=Jstr-2,Jend+1
            do i=IstrU-2,Iend+1
              zw_nbq(i,j,k,knew)=z_w(i,j,k)
            enddo
          enddo
        enddo
        do j=Jstr-2,Jend+1
          do i=IstrU-2,Iend+1
            FY_3D(i,j,0)=0.D0
          enddo
        enddo
        do k=1,N-1
          do j=Jstr-2,Jend+1
            do i=IstrU-2,Iend+1
              FX_3D(i,j,k) =
     &          -(zw_nbq(i,j,k,knew)-zw_nbq(i,j,k,kstp))/dtgrid_nbq
              FY_3D(i,j,k)=
     &          FX_3D(i,j,k)
     &          *0.5D0*(    rho_grd(i,j,k  )
     &                   +rho_nbq(i,j,k  )*Hzr_half_nbq_inv(i,j,k  )
     &                   +rho_grd(i,j,k+1)
     &                   +rho_nbq(i,j,k+1)*Hzr_half_nbq_inv(i,j,k+1) )
            enddo
          enddo
        enddo
        do j=Jstr-2,Jend+1
          do i=IstrU-2,Iend+1
            FY_3D(i,j,N)=
     &        -(zw_nbq(i,j,N,knew)-zw_nbq(i,j,N,kstp))/dtgrid_nbq
     &         *( rho_grd(i,j,N)
     &          )
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              thetadiv3_nbq(i,j,k)=
     &    FY_3D(i,j,k)-FY_3D(i,j,k-1)
            enddo
          enddo
        enddo
      do k=0,N
        do j=Jstr,Jend
          do i=Istr,Iend
            FX_3D(i,j,k)=qdmw_nbq(i,j,k)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do k=1,N
          do i=Istr,Iend
            FC(i,k)= soundspeed2_nbq(i,j)*rho_nbq(i,j,k)
     &              -(thetaimp_nbq*soundspeed2_nbq(i,j)*dtnbq+visc2_nbq)
     &                       *(thetadiv_nbq(i,j,k)+thetadiv3_nbq(i,j,k))
            FC(i,k)=FC(i,k)*Hzr_half_nbq_inv(i,j,k)
          enddo
        enddo
        do k=1,N-1
          do i=Istr,Iend
            dum_s = FC(i,k) - FC(i,k+1)
            qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)
     &                        + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
          enddo
        enddo
        k=N
        do i=Istr,Iend
          dum_s =   FC(i,k)
          qdmw_nbq(i,j,k) = qdmw_nbq(i,j,k)
     &                      + dtnbq * ( dum_s + rw_int_nbq(i,j,k) )
        enddo
        enddo
      cff1_cst = 1.D0/(dtnbq*(thetaimp_nbq**2*csound_nbq**2
     &            *dtnbq+visc2_nbq))
      do j=Jstr,Jend
      do k=1,N
        do i=Istr,Iend
          FC(i,k) =
     &      Hzw_half_nbq_inv(i,j,k)*Hzr_half_nbq_inv(i,j,k)
          BC(i,k) =
     &      Hzw_half_nbq_inv(i,j,k-1)*Hzr_half_nbq_inv(i,j,k)
        enddo
      enddo
        k=1
        do i=Istr,Iend
          cff1 = cff1_cst
          cff=(cff1+FC(i,1)+BC(i,2))
          CF(i,1)=(-FC(i,2))
     &                                                          /cff
          DC(i,1)=qdmw_nbq(i,j,1)*cff1/cff
     &           +qdmw_nbq(i,j,0)/cff*BC(i,1)
        enddo
        do k=2,N-1
          do i=Istr,Iend
          cff1 = cff1_cst
            cff=(cff1+
     &               FC(i,k)+BC(i,k+1)
     &              +BC(i,k)*CF(i,k-1))
            CF(i,k)=(-FC(i,k+1))/cff
            DC(i,k)=(qdmw_nbq(i,j,k)*cff1+
     &               BC(i,k)*DC(i,k-1)) /cff
          enddo
        enddo
        k=N
        do i=Istr,Iend
          cff1 = cff1_cst
          cff=(cff1+FC(i,N)
     &             +BC(i,N)*CF(i,N-1))
          CF(i,N)=0.D0
          DC(i,k)=(qdmw_nbq(i,j,N)*cff1+
     &             BC(i,N)*DC(i,N-1))/cff
        enddo
        do i=Istr,Iend
          qdmw_nbq(i,j,N)=DC(i,k)
        enddo
        do k=N-1,1,-1
          do i=Istr,Iend
            qdmw_nbq(i,j,k)=DC(i,k)-CF(i,k)*qdmw_nbq(i,j,k+1)
          enddo
        enddo
      enddo
      call wnbq_bc_tile (Istr,Iend,Jstr,Jend, work)
        do j=Jstr,Jend
        do i=Istr,Iend
          FC(i,0)=0.D0
        enddo
        do k=1,N-1
          do i=Istr,Iend
            FC(i,k)=Hzw_half_nbq_inv(i,j,k) *
     &              (thetaimp_nbq*qdmw_nbq(i,j,k)
     &              +(1.D0-thetaimp_nbq)*FX_3D(i,j,k))
            thetadiv_nbq(i,j,k)=thetadiv_nbq(i,j,k)
     &                         +FC(i,k)-FC(i,k-1)
          enddo
        enddo
        do i=Istr,Iend
          FC(i,N)=Hzw_half_nbq_inv(i,j,N) *
     &            (thetaimp_nbq*qdmw_nbq(i,j,N)
     &           +(1.D0-thetaimp_nbq)*FX_3D(i,j,N))
          thetadiv_nbq(i,j,N)=thetadiv_nbq(i,j,N)
     &                        +FC(i,N)-FC(i,N-1)
        enddo
        enddo
!$acc end kernels
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv_nbq(-2,-2,1))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        thetadiv3_nbq(-2,-2,1))
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
        do j=Jstr-2,Jend+1
          do i=IstrU-2,Iend+1
            rho_nbq(i,j,k) = rho_nbq(i,j,k)
     &                       - dtfast*(thetadiv_nbq(i,j,k)+
     &                                 thetadiv3_nbq(i,j,k))
          enddo
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew)=0.D0
        enddo
      enddo
      do k=1,N
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            rhobar_nbq(i,j,knew)= rhobar_nbq(i,j,knew)
     &                             +rho_nbq(i,j,k)
     &                             +rho_grd(i,j,k)*Hzr(i,j,k)
          enddo
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=Istr-1,Iend+1
          rhobar_nbq(i,j,knew) = 1.D0+(rhobar_nbq(i,j,knew))
     &                            / (z_w(i,j,N)-z_w(i,j,0))
        enddo
      enddo
         call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                           rhobar_nbq(-2,-2,knew))
      do j=Jstr-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
     &               *rhobar_nbq(i,j,knew)
        enddo
      enddo
      if (iif.eq.nfast) then
        do k=0,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rw_nbq(i,j,k)=((qdmw_nbq(i,j,k)-rw_nbq(i,j,k))
     &                              /dtnbq-ndtnbq*rw_int_nbq(i,j,k))
     &                                            /(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo
      endif
      if (iif.eq.nfast) then
        do j=Jstr,Jend
          do i=Istr,Iend
            rhobar_nbq_avg1(i,j)=rhobar_nbq(i,j,knew)
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rho_nbq_avg1(i,j,k)=1.D0
     &                            + ( rho_nbq(i,j,k)/Hzr(i,j,k)
     &                               +rho_grd(i,j,k) )
            enddo
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=IstrU,Iend
              ru_int_nbq(i,j,k) = ru_int_nbq(i,j,k)
     &                      -ru_ext_nbq_old(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
              ru_nbq_avg2(i,j,k)=
     &                          ((qdmu_nbq(i,j,k)-ru_nbq_avg2(i,j,k))/dt
     &                   -ru_int_nbq(i,j,k)-(ru_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i-1,j,k)))*on_u(i,j)*om_u(i,j)
            enddo
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rv_int_nbq(i,j,k) = rv_int_nbq(i,j,k)
     &                      -rv_ext_nbq_old(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
              rv_nbq_avg2(i,j,k)=
     &                          ((qdmv_nbq(i,j,k)-rv_nbq_avg2(i,j,k))/dt
     &                   -rv_int_nbq(i,j,k)-(rv_ext_nbq_sum(i,j)/nfast)*
     &                      (Hz(i,j,k)+Hz(i,j-1,k)))*on_v(i,j)*om_v(i,j)
            enddo
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              rw_nbq_avg2(i,j,k)=
     &                          ((qdmw_nbq(i,j,k)-rw_nbq_avg2(i,j,k))/dt
     &                           -rw_int_nbq(i,j,k))*on_r(i,j)*om_r(i,j)
            enddo
          enddo
        enddo
      endif
!$acc end kernels
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        rhobar_nbq_avg1(-2,-2))
      call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                        rho_nbq_avg1(-2,-2,1))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                        rw_nbq_avg2(-2,-2,0))
!$acc kernels if(compute_on_device) default(present)
      do j=Jstr-1,Jend
        do i=IstrU-1,Iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
!$acc end kernels
      cff=2*dtfast
      cff1=0.5D0*weight(1,iif)
      cff2=0.5D0*weight(2,iif)
!$acc kernels if(compute_on_device) default(present)
      if (iif.eq.1) then
         DU_avg1=0.D0
         DV_avg1=0.D0
         DU_avg2=0.D0
         DV_avg2=0.D0
      endif
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew=DU_nbq(i,j) *2.D0
     &                                     *umask_wet(i,j)
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew)
     &                             +cff1*on_u(i,j)*( DUnew
     &                                                   )
          DU_avg2(i,j)=DU_avg2(i,j)+cff2*on_u(i,j)*( DUnew
     &                                                   )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          DVnew=DV_nbq(i,j) *2.D0
     &                                     *vmask_wet(i,j)
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew)
     &                              +cff1*om_v(i,j)*(DVnew
     &                                                   )
          DV_avg2(i,j)=DV_avg2(i,j)+cff2*om_v(i,j)*( DVnew
     &                                                   )
        enddo
      enddo
      M2bc_nbq_flag=.false.
!$acc end kernels
      call u2dbc_tile (Istr,Iend,Jstr,Jend, work)
      call v2dbc_tile (Istr,Iend,Jstr,Jend, work)
!$acc kernels if(compute_on_device) default(present)
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        DO j=Jstr,Jend
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask_wet(Istr,j)
        END DO
        DO j=Jstr,Jend
          vbar(Istr-1,j,knew)=vbar(Istr-1,j,knew)*vmask_wet(Istr-1,j)
        END DO
      END IF
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        DO j=Jstr,Jend
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask_wet(Iend+1,j)
        END DO
        DO j=Jstr,Jend
          vbar(Iend+1,j,knew)=vbar(Iend+1,j,knew)*vmask_wet(Iend+1,j)
        END DO
      END IF
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        DO i=IstrU,Iend
          ubar(i,Jstr-1,knew)=ubar(i,Jstr-1,knew)*umask_wet(i,Jstr-1)
        END DO
        DO i=IstrU,Iend
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask_wet(i,Jstr)
        END DO
      END IF
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        DO i=Istr,Iend
          ubar(i,Jend+1,knew)=ubar(i,Jend+1,knew)*umask_wet(i,Jend+1)
        END DO
        DO i=Istr,Iend
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask_wet(i,Jend+1)
        END DO
      END IF
!$acc end kernels
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        ubar(-2,-2,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        vbar(-2,-2,knew))
!$acc kernels if(compute_on_device) default(present)
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=(h(Istr-1,j)+zeta(Istr-1,j,knew))
     &                          *rhobar_nbq(Istr-1,j,knew)
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=(h(Iend+1,j)+zeta(Iend+1,j,knew))
     &                          *rhobar_nbq(Iend+1,j,knew)
        enddo
      endif
!$acc end kernels
      cff1=0.5D0*weight(1,iif)
      cff2=0.5D0*weight(2,iif)
!$acc kernels if(compute_on_device) default(present)
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        do j=JstrR,JendR
          cff=(Dnew(IstrU-1,j)+Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
     &                                              )*on_u(IstrU-1,j)
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)+cff1*cff
          DU_avg2(IstrU-1,j)=DU_avg2(IstrU-1,j)+cff2*cff
        enddo
        do j=Jstr,Jend
          cff=(Dnew(Istr-1,j)+Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
     &                                               )*om_v(Istr-1,j)
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)+cff1*cff
          DV_avg2(Istr-1,j)=DV_avg2(Istr-1,j)+cff2*cff
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        do j=JstrR,JendR
          cff=(Dnew(Iend+1,j)+Dnew(Iend,j))*(ubar(Iend+1,j,knew)
     &                                          )*on_u(Iend+1,j)
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)+cff1*cff
          DU_avg2(Iend+1,j)=DU_avg2(Iend+1,j)+cff2*cff
        enddo
        do j=Jstr,Jend
          cff=(Dnew(Iend+1,j)+Dnew(Iend+1,j-1))*(vbar(Iend+1,j,knew)
     &                                              )*om_v(Iend+1,j)
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)+cff1*cff
          DV_avg2(Iend+1,j)=DV_avg2(Iend+1,j)+cff2*cff
        enddo
      endif
        do j=Jstr-1,Jend+1
          do i=Istr-1,Iend+1
            Dnew(i,j)=(zeta(i,j,knew)+h(i,j))
     &                                          *rhobar_nbq(i,j,knew)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            zeta(i,j,knew)=( (h(i,j)+zeta(i,j,kstp))
     &                                          *rhobar_nbq(i,j,kstp)
     &                                + (dtfast*pm(i,j)*pn(i,j)*0.5D0*(
     &      (Dnew(i  ,j)+Dnew(i-1,j))*(ubar(i  ,j,knew)
     &                                                 )*on_u(i  ,j)
     &     -(Dnew(i+1,j)+Dnew(i  ,j))*(ubar(i+1,j,knew)
     &                                                 )*on_u(i+1,j)
     &     +(Dnew(i,j  )+Dnew(i,j-1))*(vbar(i,j  ,knew)
     &                                                 )*om_v(i,j  )
     &     -(Dnew(i,j+1)+Dnew(i,j  ))*(vbar(i,j+1,knew)
     &                                                 )*om_v(i,j+1))) )
     &                                             /rhobar_nbq(i,j,knew)
     &                                                          - h(i,j)
          enddo
        enddo
!$acc end kernels
      if (iif.eq.nfast) then
!$acc parallel loop if(compute_on_device) default(present)
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=zeta(i,j,knew)
          enddo
        enddo
      endif
        call set_depth_tile(Istr,Iend,Jstr,Jend)
        call grid_nbq_tile(Istr,Iend,Jstr,Jend,
     &                     Hzw_half_nbq_inv,   Hzr_half_nbq_inv,
     &                     Hzw_half_nbq_inv_u, Hzw_half_nbq_inv_v,
     &                     Hzu_half_qdmu,      Hzv_half_qdmv)
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                        zeta(-2,-2,knew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg1(-2,-2,nnew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg1(-2,-2,nnew))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                        DU_avg2(-2,-2))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                        DV_avg2(-2,-2))
      if (iif.eq.nfast) then
!$acc parallel loop if(compute_on_device) default(present)
        do k=1,N
          do j=Jstr-2,Jend+1
            do i=IstrU-2,Iend+1
              rho_bak(i,j,k)=rho(i,j,k)
            enddo
           enddo
         enddo
      endif
      if (iif.eq.nfast) then
!$acc parallel loop if(compute_on_device) default(present)
        do k=1,N
          do j=Jstr-2,Jend+1
            do i=IstrU-2,Iend+1
              Hz_correct(i,j,k)=Hz(i,j,k)
            enddo
          enddo
        enddo
        call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                          DU_avg2(-2,-2))
        call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                          DV_avg2(-2,-2))
!$acc kernels if(compute_on_device) default(present)
        do j=Jstr,Jend
          do i=Istr,Iend
            dum_s=0.D0
            do k=1,N
              dum_s=dum_s+Hz(i,j,k)-Hz_bak(i,j,k)
            enddo
            do k=1,N
              Hz(i,j,k)=Hz(i,j,k)
     &                -(dum_s
     &             + (DU_avg2(i+1,j)-DU_avg2(i,j)
     &               +DV_avg2(i,j+1)-DV_avg2(i,j)
     &               )*pm(i,j)*pn(i,j)
     &               *dt)
     &               /(z_w(i,j,N)-z_w(i,j,0))
     &               *(z_w(i,j,k)-z_w(i,j,k-1))
            enddo
          enddo
        enddo
!$acc end kernels
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          Hz(-2,-2,1))
      endif
!$acc kernels if(compute_on_device) default(present)
      VMAXL=100.D0
      VMAX=0.D0
      do j=Jstr,Jend
        do i=Istr,Iend
          VMAX=max(VMAX,abs(ubar(i,j,knew)))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          VMAX=max(VMAX,abs(vbar(i,j,knew)))
        enddo
      enddo
!$acc end kernels
      IF (VMAX.GT.VMAXL) THEN
!$acc update if(compute_on_device) host(ubar(Istr:Iend,Jstr:Jend,knew))
!$acc update if(compute_on_device) host(vbar(Istr:Iend,Jstr:Jend,knew))
        cff1=maxval(abs(ubar(Istr:Iend,Jstr:Jend,knew)))
        cff2=maxval(abs(vbar(Istr:Iend,Jstr:Jend,knew)))
        IF (cff1>cff2) THEN
          ijmax= maxloc(abs(ubar(Istr:Iend,Jstr:Jend,knew)))
        ELSE
          ijmax= maxloc(abs(vbar(Istr:Iend,Jstr:Jend,knew)))
        ENDIF
        imax=ijmax(1)
        jmax=ijmax(2)
        imax=imax+iminmpi-1
        jmax=jmax+jminmpi-1
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =  STEP3D_FAST:   ABNORMAL JOB END    = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '
        write(stdout,'(A,I4)') '      mynode =',mynode
        if (VMAX.eq.666.D0) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')
     &                                       '  IMAX JMAX   =',imax,jmax
        write(stdout,'(A,2I6/)')
     &                                       '  IINT IFAST  =',iic,iif
        may_day_flag=1
        call mpi_abort (MPI_COMM_WORLD, err)
      ENDIF
      return
      end
