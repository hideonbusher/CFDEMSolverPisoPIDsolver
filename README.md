This repository provides a modified cfdemSolverPisoPIDsolver based on the open-source CFDEMcoupling framework.

The solver introduces a PID-controlled streamwise body force (acceleration [m/s^2] in the solver) to maintain a target flow discharge in periodic CFD-DEM simulations. 

The repository includes the modified solver files, the corresponding configuration examples, and the configuration files for the proposed CG semi-resolved algorithm used in the study.
We also provide tutorials for the particle-laden channel flow based on the proposed framework.

The CG-files should be compiled on the standard version of CFDEM.

Standard version CFDEM framework:  
https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC

Contact email for any related questions: yuxiang.liu1998@gmail.com
