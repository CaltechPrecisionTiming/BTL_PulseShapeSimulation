#!/bin/tcsh

hostname

echo "Arguments: $*"
set config=$1
set NEvents=$2
set outputFile=$3
set outputDirectory=$4

echo " "; echo "Initialize CMSSW"; echo " "
#setenv KRB5CCNAME /home/sixie/.krb5/ticket
set workDir=`pwd`

#setenv SCRAM_ARCH slc6_amd64_gcc472
cd    /afs/cern.ch/work/s/sixie/public/releases/run2/Timing/CMSSW_9_0_2/src/
eval `scramv1 runtime -sh`
cd -

pwd
# env

cp $CMSSW_BASE/src/BTL_PulseShapeSimulation/SimulatePulseShape ./

echo " "; echo "Show where we are"; echo " "
hostname
pwd
## env

klist

#setenv STAGE_SVCCLASS cmsprod

# Get ready to run in your home directory
echo " "; echo "Starting job now"; echo " ";
echo ./SimulatePulseShape --config_file=${config} --output_file=${outputFile} --n_experiments=${NEvents}
./SimulatePulseShape --config_file=${config} --output_file=${outputFile} --n_experiments=${NEvents}

ls -ltr 

echo $outputFile 
echo $outputDirectory

mkdir -p $outputDirectory
cp -v $outputFile $outputDirectory

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
