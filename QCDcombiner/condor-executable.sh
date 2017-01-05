#!/bin/bash

# requires 2 argument inputs:   
# RUN_DIR - running directory (CMSSW_X_Y_Z/subdir)   
# mode.  should be 1 or 2

#
# header 
#


RUN_DIR=$1
RUN_MODE=$2

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



#
# run c
#
./main $RUN_MODE >& main.log



#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
