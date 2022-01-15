#!/bin/bash

#set filelist=${1}
#set outfile=${2}
#set data=${3}

#Check if you need any of such similar sourcing as next two lines
#source /uscmst1/prod/sw/cms/cshrc prod
#setenv SCRAM_ARCH slc5_amd64_gcc462

#I think you need to set your root environment. Best is to checkout a CMSSW package, cmsenv will give you access to root
cd /home/work/spandey/public/HGCAL/CERN_TB/octTB/CMSSW_9_3_0/src/pion_analysis/newTree_showerShapes_inclusive_SS_faster
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

make

echo $PWD
./analyzeHGCOctTB  $1 $2 $3 $4
