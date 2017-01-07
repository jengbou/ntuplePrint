to compile and run:

make
./main


to run on condor (since thie job takes a while)
first edit condor-executable.sh and condor_jobs.jdl and change all the directory names

then do

condor_submit condor_jobs.jdl


to see your job running do

condor_q -submitter your_user_name
