to compile:

make

to run on background
./main 0

to run on signal
./main 1

(right now it doesn't run on the full sample.  edit main.cc and it should be obvious now to change this)

to run on condor (since thie job takes a while)
first edit condor-executable.sh and condor_jobs.jdl and change all the directory names
Also, for the arguments, set the first to SIGNAL and the last to 1 to run on signal
set the first to BACKGROUND and the last to 0 to run on background


then do

condor_submit condor_jobs.jdl


to see your job running do

condor_q -submitter your_user_name
