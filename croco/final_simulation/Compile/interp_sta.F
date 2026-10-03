!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#ifdef STATIONS
      subroutine interp_r3d_sta (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "sta.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2, cff1, cff2

      do id=1,nfltmax
        iflt=indx(id)

        k1=int(stainfo(istazgrd,iflt))
        r2=stainfo(istazgrd,iflt) - float(k1)
        r1=1.-r2
        k2=min(k1+1, N)

        i1=int(stainfo(istaxgrd,iflt))
        i2=i1+1
        p2=stainfo(istaxgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(stainfo(istaygrd,iflt))
        j2=j1+1
        q2=stainfo(istaygrd,iflt)-float(j1)
        q1=1.0-q2

# ifdef MASKING
        cff1=p1*q1*rmask(i1,j1)*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &      +p2*q1*rmask(i2,j1)*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &      +p1*q2*rmask(i1,j2)*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &      +p2*q2*rmask(i2,j2)*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2.gt.0.) then
          stadata(ifld,iflt)=cff1/cff2
        else
          stadata(ifld,iflt)=0.0
        endif
# else
        stadata(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                    +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                    +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                    +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
# endif
      enddo
      return
      end
!
!======================================================================
!
      subroutine interp_w3d_sta (A, ifld, nfltmax,indx)
      implicit none
# include "param.h"
# include "sta.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,0:N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2

      do id=1,nfltmax
        iflt=indx(id)

        k1=int(stainfo(istazgrd,iflt))
        r2=stainfo(istazgrd,iflt) - float(k1)
        r1=1.-r2
        k1=max(k1,   0)
        k2=min(k1+1, N)

        i1=int(stainfo(istaxgrd,iflt))
        i2=i1+1
        p2=stainfo(istaxgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(stainfo(istaygrd,iflt))
        j2=j1+1
        q2=stainfo(istaygrd,iflt)-float(j1)
        q1=1.0-q2

        stadata(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                    +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                    +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                    +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
      enddo
      return
      end
!
!======================================================================
!
      subroutine interp_r2d_sta (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "sta.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY)
      integer id,iflt, i1,i2, j1,j2
      real p1,p2, q1,q2

      do id=1,nfltmax
        iflt=indx(id)

        i1=int(stainfo(istaxgrd,iflt))
        i2=i1+1
        p2=stainfo(istaxgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(stainfo(istaygrd,iflt))
        j2=j1+1
        q2=stainfo(istaygrd,iflt)-float(j1)
        q1=1.0-q2

        stadata(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                    +q2*(p1*A(i1,j2) + p2*A(i2,j2))
      enddo
      return
      end
!
!======================================================================
!
      subroutine fill_sta_ini
      implicit none
# include "param.h"
# include "sta.h"
# include "scalars.h"
# include "nc_sta.h"
      integer iflt
      do iflt=1,nstas
         stadata(istadpt,iflt)=stainfo(istazgrd,iflt)    ! Either depth or
         if (wrtsta(indxstaGrd)) then
           stadata(istazgrd,iflt)=stainfo(istazgrd,iflt) ! sigma level position
                                                         ! is unavailaible
           stadata(istaxgrd,iflt)=stainfo(istaxgrd,iflt) ! for initial storage
           stadata(istaygrd,iflt)=stainfo(istaygrd,iflt)
         endif
      enddo
      return
      end
!
!======================================================================
!
      subroutine interp_r2d_sta_ini (A, ifld)
      implicit none
# include "param.h"
# include "sta.h"
# include "scalars.h"
      integer ifld
      real A(GLOBAL_2D_ARRAY)
      integer iflt, i1,i2, j1,j2
      real p1,p2, q1,q2

      do iflt=1,nstas

        i1=int(stainfo(istaxgrd,iflt))
        i2=i1+1
        p2=stainfo(istaxgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(stainfo(istaygrd,iflt))
        j2=j1+1
        q2=stainfo(istaygrd,iflt)-float(j1)
        q1=1.0-q2

        stadata(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                    +q2*(p1*A(i1,j2) + p2*A(i2,j2))
      enddo
      return
      end


!======================================================================
! Simon Treillou, 02/2026 - add vorticity to station
# ifdef SOLVE3D
      subroutine interp_r3d_sta_vrt (ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "sta.h"
# include "scalars.h"
# include "ocean3d.h"
      integer ifld, nfltmax, indx(nfltmax)
      integer id, iflt
      integer i1,i2,j1,j2,k1,k2
      real p1,p2,q1,q2,r1,r2
      real z11k1,z21k1,z12k1,z22k1
      real z11k2,z21k2,z12k2,z22k2
      real cff1, cff2
      real zeta_rho

      do id=1,nfltmax
        iflt=indx(id)

        !---- same weights as interp_r3d_sta
        k1=int(stainfo(istazgrd,iflt))
        r2=stainfo(istazgrd,iflt) - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(stainfo(istaxgrd,iflt))
        i2=i1+1
        p2=stainfo(istaxgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(stainfo(istaygrd,iflt))
        j2=j1+1
        q2=stainfo(istaygrd,iflt)-float(j1)
        q1=1.0-q2

        !---- boundary safety: vorticity stencil needs (i-1) and (j-1)
        !     If stations can sit on edges, avoid illegal access.
        if ( (i1 .lt. 2) .or. (j1 .lt. 2) ) then
          stadata(ifld,iflt)=0.0
        else

          !---- compute zeta at the 4 corners for k1 and k2
          z11k1 = zeta_rho(i1,j1,k1)
          z21k1 = zeta_rho(i2,j1,k1)
          z12k1 = zeta_rho(i1,j2,k1)
          z22k1 = zeta_rho(i2,j2,k1)

          z11k2 = zeta_rho(i1,j1,k2)
          z21k2 = zeta_rho(i2,j1,k2)
          z12k2 = zeta_rho(i1,j2,k2)
          z22k2 = zeta_rho(i2,j2,k2)

# ifdef MASKING
          ! Apply rmask weighting consistent with interp_r3d_sta
          cff1 = p1*q1*rmask(i1,j1)*(r1*z11k1 + r2*z11k2)
     &         +p2*q1*rmask(i2,j1)*(r1*z21k1 + r2*z21k2)
     &         +p1*q2*rmask(i1,j2)*(r1*z12k1 + r2*z12k2)
     &         +p2*q2*rmask(i2,j2)*(r1*z22k1 + r2*z22k2)

          cff2 = q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &         +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

          if (cff2.gt.0.) then
            stadata(ifld,iflt)=cff1/cff2
          else
            stadata(ifld,iflt)=0.0
          endif
# else
          stadata(ifld,iflt)=
     &        p1*q1*(r1*z11k1 + r2*z11k2)
     &       +p2*q1*(r1*z21k1 + r2*z21k2)
     &       +p1*q2*(r1*z12k1 + r2*z12k2)
     &       +p2*q2*(r1*z22k1 + r2*z22k2)
# endif

        endif  ! bounds check

      enddo
      return
      end
# endif /* SOLVE3D */


      real function zeta_rho(i,j,k)
      implicit none
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean3d.h"
      integer i,j,k
      real cff1, cff2

      ! pm and pn are inverse grid spacing; average to RHO like you did
      cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
      cff2 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

      zeta_rho =
     &      ( v(i,j,k,nnew) - v(i-1,j,k,nnew) ) * cff1
     &    - ( u(i,j,k,nnew) - u(i,j-1,k,nnew) ) * cff2

      return
      end
!======================================================================

#else
      subroutine interp_sta_empty
      end
#endif
