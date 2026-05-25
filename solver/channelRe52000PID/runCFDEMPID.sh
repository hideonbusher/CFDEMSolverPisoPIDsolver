#!/bin/bash

. ~/.bashrc

source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cd $casePath/CFD
blockMesh

logpath=$casePath
headerText="CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPisoPIDsolver"
nrProcs="4"
machineFileName="none"
debugMode="off"
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"

parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode "true"
