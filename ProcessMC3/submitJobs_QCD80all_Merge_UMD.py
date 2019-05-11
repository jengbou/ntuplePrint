#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date
import time
#sys.path.append(os.path.abspath(os.path.curdir))
from collections import OrderedDict
from Plotter import parseInputArgs
options = parseInputArgs()

Pmode = 2
ProdTag = options.outtag
OutDir  = "/data/users/jengbou/histos"
WorkDir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/temp/"+ProdTag

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")

sTags = {}
## Remember to change corresponding nfiles in main.cc
sTags["modelA"]=["1"]
sTags["modelB"]=["2"]


jobTags = OrderedDict(sorted(sTags.items(), key=lambda x: x[1]))

#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = condor-executableMerge.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet" && machine != "r510-0-1.privnet" && machine !="r510-0-9.privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 %(IMODE)s 0 0 %(OUTDIR)s %(PMODE)s
Queue 1
"""
#########################################


for k,v in jobTags.items():
    print "Submitting jobs for [%-20s]"%k
    sTag_ = k
    Imode_ = v[0]

    dirname = "jobs/%s_Merge_%s"%(sTag_,fTag)

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

    kw["MYPREFIX"]   = ProdTag
    kw["WORKDIR"]    = WorkDir
    kw["OUTDIR"]     = OutDir
    kw["PMODE"]      = Pmode
    kw["IMODE"]      = Imode_
    kw["SAMPLENAME"] = "%s_Merge"%(sTag_)

    script_str = condor_script_template % kw
    f = open("%s/condor_jobs_%s_Merge.jdl"%(dirname,sTag_), 'w')
    f.write(script_str)
    f.close()

    condorcmd = "condor_submit %s/condor_jobs_%s_Merge.jdl"%(dirname,sTag_)
    print 'condorcmd: ', condorcmd
    print 'Executing condorcmd'
    p=subprocess.Popen(condorcmd, shell=True)
    p.wait()
    time.sleep(2)

    print "Histos output dir: %s/%s"%(OutDir,ProdTag)

