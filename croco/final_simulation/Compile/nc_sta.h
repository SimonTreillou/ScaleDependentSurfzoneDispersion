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
! This is include file "nc_sta.h".
! ==== == ======= ==== ============
!
! stafield     Number of station fields for output
! wrtsta       Logical vector with flags for output
! indxsta[...] Index of logical flag to output several fields
!       Grd  - grid level
!       Temp - temp
!       Salt - Salt
!       Rho  - Density
!       Vel  - u and v components
! ncidsta      id of station output file
! nrecsta      step to output station data
! sta[...]     several reference names of netcdf output
! staname      station output filename
! staposname   station input data filename

      integer stafield
!======================================================================
! Simon Treillou, 02/2026 - add tracer concentration and vorticity to station
#ifdef PASSIVE_TRACER
      parameter(stafield=7)
#else
      parameter(stafield=6)
#endif

      integer indxstaGrd, indxstaTemp, indxstaSalt,
     &        indxstaRho, indxstaVel, indxstaVrt
#ifdef PASSIVE_TRACER
     &      , indxstaPtr
      parameter (     indxstaGrd=1, indxstaTemp=2,
     &                indxstaSalt=3, indxstaRho=4,  indxstaVel=5,
     &                indxstaVrt=6, indxstaPtr=7)
#else
      parameter (     indxstaGrd=1, indxstaTemp=2,
     &                indxstaSalt=3, indxstaRho=4,  indxstaVel=5,
     &                indxstaVrt=6)
#endif
!======================================================================


      integer ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
!======================================================================
! Simon Treillou, 02/2026 - add vorticity to station
     &      , staVrt
!======================================================================
#ifdef SPHERICAL
     &      , staLon,     staLat
#else
     &      , staX,       staY
#endif
#ifdef SOLVE3D
     &      , staDepth,   staDen
# ifdef TEMPERATURE
     &      , staTemp
# endif
# ifdef SALINITY
     &      , staSal
# endif
!======================================================================
! Simon Treillou, 02/2026 - add tracer concentration to station
# ifdef PASSIVE_TRACER
     &      , staPtr(NT)
# endif
!======================================================================
# ifdef MUSTANG
     &      , staMUS(NT-2)
# endif
#endif
      logical wrtsta(stafield)

      common/incscrum_sta/
     &        ncidsta,    nrecsta,    staGlevel
     &      , staTstep,   staTime,    staXgrd,   staYgrd
     &      , staZgrd,    staZeta,    staU,      staV
!======================================================================
! Simon Treillou, 02/2026 - add vorticity to station
     &      , staVrt
!======================================================================
#ifdef SPHERICAL
     &      , staLon,     staLat
#else
     &      , staX,       staY
#endif
#ifdef SOLVE3D
     &      , staDepth,   staDen
# ifdef TEMPERATURE
     &      ,   staTemp
# endif
# ifdef SALINITY
     &      , staSal
# endif
!======================================================================
! Simon Treillou, 02/2026 - add tracer concentration to station
# ifdef PASSIVE_TRACER
     &      , staPtr
# endif
!======================================================================
# ifdef MUSTANG
     &      , staMUS
# endif

#endif
     &      , wrtsta


      character*80  staname,   staposname
      common /cncscrum_sta/ staname,   staposname
