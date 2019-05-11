#!/bin/sh

PRODTAG="analysis_20180126_v0_p20180427_UMD_cuts1to8_unblinded_r1"
#PRODTAG="analysis_20180126_v0_p20180423_UMD_cut2"

#CutsToRun=("2")
CutsToRun=("1" "2" "3" "4" "5" "6" "7" "8")
#CutsToRun=("9" "10")

backupcmd=`tar -czvf ./backup/${PRODTAG}.tgz ./Makefile ./*.cc ./*.h ./*.py ./*.sh ./*.cpp ./Plotter/*.C ./Plotter/*.h ./Plotter/*.py ./FakeRateHistograms/*.root ./BtagFiles/*.csv`
#echo $backupcmd
tempdir="./temp/$PRODTAG"
mkdir -p $tempdir
cp ./backup/${PRODTAG}.tgz $tempdir/
cd $tempdir
unzipcmd=`tar -xzvf ./${PRODTAG}.tgz`
status=`make`
#echo $status

starttime=`date "+%H:%M:%S"`
echo "starttime:"$starttime

sleep 10
echo "Done compiling in temp working dir! Submitting jobs in 60 sec!"
sleep 60

python submitJobs_QCD80all_noMerge_UMD.py -t $PRODTAG

#sleep 1800
sleep 60

jobMin=`condor_q -submitter jengbou |grep jengbou|grep condor- |awk -v stime=$starttime '$4>=stime { line[NR] = $1 }END {print line[1]}'`
jobMax=`condor_q -submitter jengbou |grep jengbou|grep condor- |awk -v stime=$starttime '$4>=stime { line[NR] = $1 }END {print line[NR]}'`

lastjobidfile="/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/lastcondorjobsID.txt"

lastJobID=0
while read -r line
  do
  lastJobID="$line"
  echo "Last condor JobID in previous runs: $lastJobID"
done < "$lastjobidfile"

echo "jobMin = "$jobMin
if [[ -z "$jobMin" ]]; then
    echo "jobMin = "$jobMin
    echo "jobMin doesn't exists. Use lastJobID from previous run."
    jobMin=$lastJobID
fi

if [ ${jobMin%.*} -lt ${lastJobID%.*} ]; then
    jobMin=$lastJobID
fi

if [[ -z "$jobMax" ]]; then
    jobMax=$( bc <<< "$lastJobID+1.0")
fi

echo "Job numbers: ["$jobMin" to "$jobMax"]"
echo "$jobMax" > "$lastjobidfile"

COUNTER0=`condor_q -submitter jengbou |awk -v jobmin="$jobMin" -v jobmax="$jobMax" '$1>jobmin && $1<=jobmax {print $1}'|wc|awk '{print $1}'`
until [  $COUNTER0 -lt 1  ]; do
    sleep 10
    echo COUNTER $COUNTER0
    COUNTER0=`condor_q -submitter jengbou |awk -v jobmin="$jobMin" -v jobmax="$jobMax" '$1>jobmin && $1<=jobmax {print $1}'|wc|awk '{print $1}'`
done
echo "Processing done! Submitting jobs for merging no norm now!"


starttime=`date "+%H:%M:%S"`
echo "starttime:"$starttime
sleep 60

## merging
## Data and QCDHT100toInf

python submitJobs_QCD80all_MergeNoNorm_UMD.py -t $PRODTAG

#sleep 60
## model A and B
#python submitJobs_QCD80all_Merge_UMD.py -t $PRODTAG

#sleep 1800
#sleep 10

jobMin=`condor_q -submitter jengbou |grep jengbou|grep condor- |awk -v stime=$starttime '$4>=stime { line[NR] = $1 }END {print line[1]}'`
jobMax=`condor_q -submitter jengbou |grep jengbou|grep condor- |awk -v stime=$starttime '$4>=stime { line[NR] = $1 }END {print line[NR]}'`

while read -r line
  do
  lastJobID="$line"
  echo "Last condor JobID: $lastJobID"
done < "$lastjobidfile"

echo "jobMin = "$jobMin
if [[ -z "$jobMin" ]]; then
    echo "jobMin = "$jobMin
    echo "jobMin doesn't exists. Use lastJobID from previous run."
    jobMin=$lastJobID
fi

if [ ${jobMin%.*} -lt ${lastJobID%.*} ]; then
    jobMin=$lastJobID
fi

if [[ -z "$jobMax" ]]; then
    jobMax=$( bc <<< "$lastJobID+1.0")
fi

echo "Job numbers: ["$jobMin" to "$jobMax"]"
echo "Last condor JobID in this run: $jobMax"
echo "$jobMax" > "$lastjobidfile"

COUNTER0=`condor_q -submitter jengbou |awk -v jobmin="$jobMin" -v jobmax="$jobMax" '$1>jobmin && $1<=jobmax {print $1}'|wc|awk '{print $1}'`
until [  $COUNTER0 -lt 1  ]; do
    sleep 10
    echo COUNTER $COUNTER0
    COUNTER0=`condor_q -submitter jengbou |awk -v jobmin="$jobMin" -v jobmax="$jobMax" '$1>jobmin && $1<=jobmax {print $1}'|wc|awk '{print $1}'`
done
echo "Processing done! Submitting jobs for merging norm now!"

for element in ${CutsToRun[@]}
  do
  echo "MergeNorm for Cutset["$element"]"
  python submitJobs_QCD80all_MergeNorm_UMD.py -t $PRODTAG -c $element
done
