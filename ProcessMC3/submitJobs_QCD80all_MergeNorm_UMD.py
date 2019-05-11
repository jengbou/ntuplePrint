#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date
import time
#sys.path.append(os.path.abspath(os.path.curdir))
from collections import OrderedDict
from Plotter import parseInputArgs
options = parseInputArgs()

Pmode = 9
ProdTag = options.outtag
CutSet  = "Cutset"+options.cutset
OutDir  = "/data/users/jengbou/histos"
WorkDir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/temp/"+ProdTag

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")

sTags = {}
sTags["Data"]=["10",6767,120] # G, H only
## YH's ntuple 8/23
##sTags["QCD80HT500to700"]=["14",1075,100]
##sTags["QCD80HT700to1000"]=["15",1904,50]
##sTags["QCD80HT1000to1500"]=["16",1788,60]
##sTags["QCD80HT1500to2000"]=["17",1708,60]
##sTags["QCD80HT2000toInf"]=["18",751,40]
## YH's ntuple 11/03
##sTags["QCD80HT500to700"]=["14",1092,100]
##sTags["QCD80HT700to1000"]=["15",2666,60]
##sTags["QCD80HT1000to1500"]=["16",2218,60]
##sTags["QCD80HT1500to2000"]=["17",1942,50]
##sTags["QCD80HT2000toInf"]=["18",864,40]
## Note: the num of files per job should be the same as in NoNorm.py
#test
##sTags["Data"]=["10",30,20]
##sTags["QCD80HT700to1000"]=["15",30,20]

jobTags = OrderedDict(sorted(sTags.items(), key=lambda x: x[1]))

#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = condor-executableMergeNorm_New.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet" && machine !="compute-0-7.privnet" && machine !="compute-0-10.privnet" && machine !="compute-0-8.privnet" && machine !="r510-0-6.privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/%(CUTIDX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/%(CUTIDX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/%(CUTIDX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 %(IMODE)s 0 0 %(OUTDIR)s %(PMODE)s %(RNGL)s %(RNGH)s %(FIDX)s %(CUTIDX)s
Queue 1
"""
#########################################


for k,v in jobTags.items():
    print "Submitting jobs for [%-20s]"%k
    sTag_ = k
    Imode_ = v[0]
    numFiles_ = v[1]
    filesPerJob_ = v[2]
    numMerg_ = int(numFiles_/filesPerJob_) if numFiles_%filesPerJob_==0 else int(numFiles_/filesPerJob_)+1

    dirname = "jobs/%s_MergeNorm_%s_%s"%(sTag_,CutSet,fTag)

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
        os.makedirs(OutDir+"/"+ProdTag+"/"+CutSet)
    except:
        pass

    try:
        os.makedirs(OutDir+"/"+ProdTag+"/"+CutSet+"/logs")
    except:
        pass

    kw = {}

    kw["MYPREFIX"]   = ProdTag
    kw["WORKDIR"]    = WorkDir
    kw["OUTDIR"]     = OutDir
    kw["PMODE"]      = Pmode
    kw["IMODE"]      = Imode_
    kw["FIDX"]       = "%i"%(numMerg_)
    kw["RNGL"]       = "1"
    kw["RNGH"]       = "%i"%(numMerg_)
    kw["CUTIDX"]     = CutSet

    kw["SAMPLENAME"] = "%s_MergeNorm"%(sTag_)
    script_str = condor_script_template % kw
    f = open("%s/condor_jobs_%s_MergeNorm.jdl"%(dirname,sTag_), 'w')
    f.write(script_str)
    f.close()

    condorcmd = "condor_submit %s/condor_jobs_%s_MergeNorm.jdl"%(dirname,sTag_)
    print 'condorcmd: ', condorcmd
    print 'Executing condorcmd'
    p=subprocess.Popen(condorcmd, shell=True)
    p.wait()
    time.sleep(2)

    print "Histos output dir: %s/%s/%s"%(OutDir,ProdTag,CutSet)

