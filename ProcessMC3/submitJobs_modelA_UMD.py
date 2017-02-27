#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date, time
#sys.path.append(os.path.abspath(os.path.curdir))

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")
dirname = "jobs/ModelA_%s"%fTag

try:
    os.makedirs(dirname)
except:
    pass

ProdTag = "analysis_20170223_v0_p20170227_UMD"
OutDir  = "/data/users/jengbou/histos"
WorkDir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD"

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


#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = condor-executableNew.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 1 0 1 %(OUTDIR)s
Queue 1
"""
#########################################

kw = {}

kw["MYPREFIX"]  = ProdTag
kw["WORKDIR"]   = WorkDir
kw["OUTDIR"]    = OutDir
kw["SAMPLENAME"] = "modelA"

script_str = condor_script_template % kw
f = open("%s/condor_jobs_%s.jdl"%(dirname,kw["SAMPLENAME"]), 'w')
f.write(script_str)
f.close()

condorcmd = "condor_submit %s/condor_jobs_%s.jdl"%(dirname,kw["SAMPLENAME"])
print 'condorcmd: ', condorcmd
print 'Executing condorcmd'
p=subprocess.Popen(condorcmd, shell=True)
p.wait()

print "Histos output dir: %s/%s"%(OutDir,ProdTag)
