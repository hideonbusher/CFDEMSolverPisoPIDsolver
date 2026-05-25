This repository provides a modified cfdemSolverPisoPIDsolver based on the open-source CFDEMcoupling framework.

The solver introduces a PID-controlled streamwise body force (actually the acceleration [m/s^2] in the solver) to maintain a target flow discharge in periodic CFD-DEM simulations. 

The repository includes the modified solver files, the corresponding configuration examples, and the configuration files for the proposed point-cloud-based coarse-graining semi-resolved algorithm used in the study.
We also provide the configuration files of the improved semi-resolved algorithm based on the point-cloud coarse-graining technique.

Base framework:  
https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC
