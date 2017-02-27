#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date, time
#sys.path.append(os.path.abspath(os.path.curdir))

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")
sTag = "QCD74HT2000toInf"
dirname = "jobs/%s_%s"%(sTag,fTag)

try:
    os.makedirs(dirname)
except:
    pass

Pmode = 1
NumFiles = 261
FilesPerJob = 20

ProdTag = "analysis_20170223_v0_p20170227"
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
Executable = condor-executableNoMerge.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 7 0 1 %(OUTDIR)s %(PMODE)s %(RNGL)s %(RNGH)s
Queue 1
"""
#########################################

kw = {}

kw["MYPREFIX"]  = ProdTag
kw["WORKDIR"]   = WorkDir
kw["OUTDIR"]    = OutDir
kw["PMODE"]     = Pmode

for fn in xrange(1,NumFiles+1):
    if fn%FilesPerJob == 1:
        rUp = fn+FilesPerJob - 1
        if rUp >= NumFiles: rUp = NumFiles

        kw["RNGL"] = "%i"%(fn)
        kw["RNGH"] = "%i"%(rUp)
        kw["SAMPLENAME"]    = "%s_%ito%i"%(sTag,fn,rUp)
        #print kw["RNGL"], kw["RNGH"]

        script_str = condor_script_template % kw
        f = open("%s/condor_jobs_%s_noMerge_%ito%i.jdl"%(dirname,sTag,fn,rUp), 'w')
        f.write(script_str)
        f.close()

        condorcmd = "condor_submit %s/condor_jobs_%s_noMerge_%ito%i.jdl"%(dirname,sTag,fn,rUp)
        print 'condorcmd: ', condorcmd
        print 'Executing condorcmd'
        p=subprocess.Popen(condorcmd, shell=True)
        p.wait()

        print "Histos output dir: %s/%s"%(OutDir,ProdTag)

