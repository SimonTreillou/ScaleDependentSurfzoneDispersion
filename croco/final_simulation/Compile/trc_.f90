










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







! PROCESSES:

! Vertical shear (VS) only
! Short-Crested waves (SC) only
!  and SC
!  and Longshore Drift (LD)
! , LD and SC

! Tracer is active or not

! Point source release if active, longshore uniform else


! When only , need for a perturbation 
! to trigger the instability

! VSSC

! SC

! VSLD

! VSLDSC

!# define GLS_KOMEGA
! Common bathymetry:
!# define CLARK 
! Case specific bathymetry


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





















!







!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh









MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :   OPA  !  1996-01  (M. Levy)  Original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   NEMO      1.0  !  2004-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------

   !!======================================================================
END MODULE trc
