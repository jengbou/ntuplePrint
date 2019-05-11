#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date
import time
#sys.path.append(os.path.abspath(os.path.curdir))
from collections import OrderedDict
from Plotter import parseInputArgs
options = parseInputArgs()

Pmode = 1
ProdTag = options.outtag
OutDir  = "/data/users/jengbou/histos"
WorkDir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/temp/"+ProdTag

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")

sTags = {}
sTags["Data"]=["10",6767,60]# G, H only
##sTags["modelA"]=["1",30,15]
##sTags["modelB"]=["2",30,15]
## YH's ntuple 8/23
##sTags["QCD80HT500to700"]=["14",1075,100]
##sTags["QCD80HT700to1000"]=["15",1904,30]
##sTags["QCD80HT1000to1500"]=["16",1788,60]
##sTags["QCD80HT1500to2000"]=["17",1708,60]
##sTags["QCD80HT2000toInf"]=["18",751,40]
## YH's ntuple 11/03
##sTags["QCD80HT500to700"]=["14",1092,100]
##sTags["QCD80HT700to1000"]=["15",2666,30]
##sTags["QCD80HT1000to1500"]=["16",2218,30]
##sTags["QCD80HT1500to2000"]=["17",1942,50]
##sTags["QCD80HT2000toInf"]=["18",864,40]

#test
##sTags["Data"]=["10",30,10]
##sTags["QCD80HT700to1000"]=["15",30,10]

jobTags = OrderedDict(sorted(sTags.items(), key=lambda x: x[1]))

##backupcmd = "tar -czvf ./backup/%s.tgz ./*.cc ./*.h ./*.py ./*.sh ./Plotter/*.C ./Plotter/*.h ./Plotter/*.py"%(ProdTag)
##print 'backupcmd: ', backupcmd
##print 'Executing backupcmd'
##p=subprocess.Popen(backupcmd, shell=True)
##p.wait()
##exit()
#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = condor-executableNoMerge.sh
+IsLocalJob = true
+IsHighPriorityJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet" && machine !="compute-0-5.privnet" && machine !="compute-0-6.privnet" && machine !="compute-0-8.privnet" && machine !="r510-0-6.privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 %(IMODE)s 0 0 %(OUTDIR)s %(PMODE)s %(RNGL)s %(RNGH)s
Queue 1
"""
#########################################


for k,v in jobTags.items():
    print "Submitting jobs for [%-20s]"%k
    sTag_ = k
    Imode_ = v[0]
    numFiles_ = v[1]
    filesPerJob_ = v[2]

    dirname = "jobs/%s_noMerge_%s"%(sTag_,fTag)

    try:
        os.makedirs(dirname)
    except:
        pass


    try:
        os.makedirs(OutDir)
    except:
        pass

    try:
        os.makedirs(OutDir+"/"+ProdTag)
    except:
        pass

    try:
        os.makedirs(OutDir+"/"+ProdTag+"/logs")
    except:
        pass

    kw = {}

    kw["MYPREFIX"]  = ProdTag
    kw["WORKDIR"]   = WorkDir
    kw["OUTDIR"]    = OutDir
    kw["PMODE"]     = Pmode
    kw["IMODE"]     = Imode_

    count = 0 
    for fn in xrange(1,numFiles_+1):
        if fn%filesPerJob_ == 1:
            rUp = fn+filesPerJob_ - 1
            if rUp >= numFiles_: rUp = numFiles_

            kw["RNGL"] = "%i"%(fn)
            kw["RNGH"] = "%i"%(rUp)
            kw["SAMPLENAME"]    = "%s_noMerge_%ito%i"%(sTag_,fn,rUp)
            #print kw["RNGL"], kw["RNGH"]

            script_str = condor_script_template % kw
            f = open("%s/condor_jobs_%s_noMerge_%ito%i.jdl"%(dirname,sTag_,fn,rUp), 'w')
            f.write(script_str)
            f.close()

            condorcmd = "condor_submit %s/condor_jobs_%s_noMerge_%ito%i.jdl"%(dirname,sTag_,fn,rUp)
            print 'condorcmd: ', condorcmd
            print 'Executing condorcmd'
            p=subprocess.Popen(condorcmd, shell=True)
            p.wait()
            time.sleep(0.2)
            #if count == 20: time.sleep(300)
            count += 1
            print "Histos output dir: %s/%s"%(OutDir,ProdTag)

