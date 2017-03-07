#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date, time
#sys.path.append(os.path.abspath(os.path.curdir))

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")
sTag = "QCD76HT100to200"
dirname = "jobs/%s_MergeNorm_%s"%(sTag,fTag)

try:
    os.makedirs(dirname)
except:
    pass

Pmode = 9
NumFiles = 10 # 3951%400+1

#ProdTag = "analysis_20170223_v0_test"
ProdTag = "analysis_20170301_v0_p20170306"
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
Executable = condor-executableMergeNorm.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 11 0 1 %(OUTDIR)s %(PMODE)s %(RNGL)s %(RNGH)s %(FIDX)s
Queue 1
"""
#########################################

kw = {}

kw["MYPREFIX"]  = ProdTag
kw["WORKDIR"]   = WorkDir
kw["OUTDIR"]    = OutDir
kw["PMODE"]     = Pmode

kw["FIDX"] = "%i"%(NumFiles)
kw["RNGL"] = "1"
kw["RNGH"] = "%i"%(NumFiles)

kw["SAMPLENAME"]    = "%s_MergeNorm"%sTag
script_str = condor_script_template % kw
f = open("%s/condor_jobs_%s_MergeNorm.jdl"%(dirname,sTag), 'w')
f.write(script_str)
f.close()

condorcmd = "condor_submit %s/condor_jobs_%s_MergeNorm.jdl"%(dirname,sTag)
print 'condorcmd: ', condorcmd
print 'Executing condorcmd'
p=subprocess.Popen(condorcmd, shell=True)
p.wait()

print "Histos output dir: %s/%s"%(OutDir,ProdTag)

