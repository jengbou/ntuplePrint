#!/bin/sh
## 
python submitJobs_QCD76all_noMerge.py
sleep 1800

COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
until [  $COUNTER -lt 7  ] && [  $COUNTER1 -lt 1  ]; do
    sleep 10
    echo COUNTER $COUNTER
    echo Jobs $COUNTER1
    COUNTER=`condor_q -submitter jengbou |wc |awk '$1>0 {print $1}'`
    COUNTER1=`condor_q -submitter jengbou |grep jobs|awk '{print $1}'`
done
echo "Processing done! Submitting jobs for merging now!"


## merging
## QCDHT100to1000
python submitJobs_QCD76all_MergeNoNorm.py
sleep 1
## model A and B and QCDHT1000+
python submitJobs_QCD76all_Merge.py

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
echo "Processing done! Submitting jobs for merging norm (HT100to200 only) now!"

python submitJobs_QCD76all_MergeNorm.py
