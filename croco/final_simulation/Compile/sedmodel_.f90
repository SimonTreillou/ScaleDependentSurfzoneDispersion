










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










MODULE sedmodel
   !!======================================================================
   !!                       ***  MODULE sedmodel   ***
   !!   Sediment model : Main routine of sediment model 
   !!======================================================================
CONTAINS

   SUBROUTINE sed_model ( kt, Kbb, Kmm, Krhs )
      INTEGER, INTENT(in) ::   kt               ! number of iteration
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level indices
   END SUBROUTINE sed_model

END MODULE sedmodel
