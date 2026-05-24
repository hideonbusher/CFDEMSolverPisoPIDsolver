#!/bin/bash

block_file=$(find . -type f -name "blockMeshDict" | head -n 1)

if [[ ! -f "$block_file" ]]; then
    exit 1
fi

read a b c < <(grep -oP 'domain\s*\(\K[^\)]+' "$block_file")

decomp_file=$(find . -type f -name "decomposeParDict" | head -n 1)

if [[ ! -f "$decomp_file" ]]; then
    exit 1
fi

read d e f < <(awk '/simpleCoeffs/,/}/ {if ($1=="n") {gsub(/[()]/,""); print $2, $3, $4}}' "$decomp_file")

cp_file=$(find . -type f -name "couplingProperties" | head -n 1)

if [[ ! -f "$cp_file" ]]; then
    exit 1
fi

cp "$cp_file" "$cp_file.bak"

sed -i -E "s/(meshx\s*)[0-9.eE+-]+;/\1$a;/" "$cp_file"
sed -i -E "s/(meshy\s*)[0-9.eE+-]+;/\1$b;/" "$cp_file"
sed -i -E "s/(meshz\s*)[0-9.eE+-]+;/\1$c;/" "$cp_file"

sed -i -E "s/(decx\s*)[0-9.eE+-]+;/\1$d;/" "$cp_file"
sed -i -E "s/(decy\s*)[0-9.eE+-]+;/\1$e;/" "$cp_file"
sed -i -E "s/(decz\s*)[0-9.eE+-]+;/\1$f;/" "$cp_file"

control_file=$(find . -type f -name "controlDict" | head -n 1)

if [[ ! -f "$control_file" ]]; then
    exit 1
fi

max_x=$(awk '/vertices/,/^\)/ {
    gsub(/[()]/,"");
    if ($0 ~ /^[[:space:]]*[0-9.+-]/) {
        count++;
        if (count==3) print $1;
    }
}' "$block_file")

convert_to_meters=$(awk '/convertToMeters/ {print $2}' "$block_file" | tr -d ';')

if [[ -z "$max_x" || -z "$a" || -z "$convert_to_meters" ]]; then
    exit 1
fi

length_ratio=$(echo "$max_x / $a * $convert_to_meters" | bc -l)

cp "$control_file" "$control_file.bak"

sed -i -E "s/^\s*(x_plane\s+)[0-9.eE+-]+;/\1$length_ratio;/" "$control_file"


#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cd $casePath/CFD
blockMesh

#--------------------------------------------------------------------------------#
#- define variables
logpath=$casePath
headerText="CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPisoPIDsolver"
nrProcs="6"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode "true"

cd $casePath/CFD/postflow ; python2 $HOME/LIGGGHTS/LPP/src/lpp.py dump.flowflow
cd $casePath/CFD
foamToVTK



