# data
``data`` contains data displayed in Clark et al. (2010) and Feddersen et al. (2011). Data was taken from figures.

# croco
## final simulation
``croco/final_simulation`` contains all the necessary Fortran and C routines to run CROCO and reproduce the simulations displayed in the paper. 

- ``input_...txt`` are input spectrum (for energy, wave angle and directional spread), read in wave_maker.h
- ``stations.in`` defines where are located stations with higher frequency output.
- ``cppdefs.h`` defines the keys to activate part of the code:
  
  - ``FLASH_RIP`` is the configuration used here.
  - ``VS``, ``SC``, ``VSSC``, ``VSLD`` or ``VSLDSC`` determine which forcings are active.
  - ``SPINUP`` determines if tracer equation is active or not

- ``croco.in`` sets the parameters:

  - ``time_steppping``: 600000 NTIMES for spinup, 1000000 or 1400000 for tracer simulations.
  - ``wave_maker``: mean wave angle and directional spread to change for VS, VS+SC and SC runs (other parameters are not important as we use observed spectrum)
  - ``initial``: spinup simulation
 
- ``jobcomp`` is used to compile the routines (``./jobcomp`` after changing with correct paths)

- ``jobsub`` to use on a cluster

**Note**: If you want to be sure to use the same CROCO code as the one used for the paper, copy all ``*.F`` and ``*.h`` files from ``Compile`` in your simulation directory. These are all the files used to produce the final code and they are the ones that will be compiled. 

# analysis
``analysis`` contains all the code used to generate the figures in the manuscript.
