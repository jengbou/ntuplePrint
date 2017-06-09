#!/bin/sh
## NOTE: if running 74 signals, remember to change to HLT_PFHT900 and
## don't run signals and QCD

PRODTAG="analysis_20170523_v0_p20170609_UMD_test1p1"

backupcmd=`tar -czvf ./backup/${PRODTAG}.tgz ./*.cc ./*.h ./*.py ./*.sh ./Plotter/*.C ./Plotter/*.h ./Plotter/*.py`
#echo $backupcmd

python submitJobs_QCD80all_noMerge_UMD.py -t $PRODTAG

#sleep 1800
sleep 600

COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
until [  $COUNTER -lt 7  ] && [  $COUNTER1 -lt 1  ]; do
    sleep 10
    echo COUNTER $COUNTER
    echo Jobs $COUNTER1
    COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
    COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
done
echo "Processing done! Submitting jobs for merging no norm now!"


## merging
## QCDHT100to1000
python submitJobs_QCD80all_MergeNoNorm_UMD.py -t $PRODTAG
#sleep 60
## model A and B
python submitJobs_QCD80all_Merge_UMD.py -t $PRODTAG

#sleep 1800
sleep 120

COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
until [  $COUNTER -lt 7  ] && [  $COUNTER1 -lt 1  ]; do
    sleep 10
    echo COUNTER $COUNTER
    echo Jobs $COUNTER1
    COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
    COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
done
echo "Processing done! Submitting jobs for merging norm now!"

python submitJobs_QCD80all_MergeNorm_UMD.py -t $PRODTAG
