# hadstruct-progs
### Hadron structure main programs, run scripts and related work flows

Main programs and run scripts for running hadron structure measurements:
* Main programs call [QuaHoG](https://github.com/g-koutsou/QuaHoG) library functions for source preparation, contractions and I/O
* Includes interfaces to [tmLQCD](https://github.com/etmc/tmLQCD) and [DDalphaAMG](https://github.com/sbacchio/DDalphaAMG) solvers
* Includes post processing bash scripts which call [QuaHoG](https://github.com/g-koutsou/QuaHoG)/utils python scripts to handle hdf5 outputs
