name="QCD76_HT500to700"

hostarea = "/data/users/eno/outputQCD/tmp"+name+"/"
exearea = "/data/users/eno/em6/CMSSW_7_6_3/src/EmergingJetAnalysis/"
thisarea = "/data/users/eno/DARKJOBS/"

npfile=input("number of files per file ")


filename = "local_"+name+'.txt'
f = open(filename)



a =1
line = f.readline()
while line:
    txtfile = open (name+"_"+str(a)+".txt", 'w') # here I save each line as a file 
    print a
    b=0
    while b<npfile and line:
       print line
       line2='file:'+line
       txtfile.write(line2) 
       line = f.readline()
       b=b+1
    txtfile.close()
    name2=name+"_"+str(a)
    jdlfile = open("condor-jobs-"+name+"_"+str(a)+".jdl","w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable = condor-executable.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Output = "+hostarea+name2+"_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name2+"_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name2+"_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = "+name2+" $(process) "+hostarea+" "+exearea+" Configuration/test/condor_cfg.py 1234567 100000 "+thisarea+name2+".txt data=0"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    a += 1 # This variable will change the file number

f.close()
