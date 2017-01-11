#!/bin/bash

# requires 4 argument inputs:   
# 1: UNIQUE_ID - any unique string identifier  
# 2: CONDOR_PROCESS - condor process number  
# RUN_DIR - running directory (CMSSW_X_Y_Z/subdir)   
# mode.  should be 0 for background or 1 for signal

#
# header 
#

UNIQUE_ID=$1
CONDOR_PROCESS=$2
RUN_DIR=$3
OPT_MODE1=$4
OPT_MODE2=$5
RUN_MODE=$6

echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"


#
# setup CMSSW software environment at UMD
#
export VO_CMS_SW_DIR=/sharesoft/cmssw
. $VO_CMS_SW_DIR/cmsset_default.sh
cd $RUN_DIR
eval `scramv1 runtime -sh`

FINAL_PREFIX_NAME=`echo ${UNIQUE_ID}_${CONDOR_PROCESS}`
FINAL_LOG=`echo $FINAL_PREFIX_NAME.log`

#
# run c
#
./main $OPT_MODE1 $OPT_MODE2 $RUN_MODE >> $FINAL_LOG 2>&1



#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
