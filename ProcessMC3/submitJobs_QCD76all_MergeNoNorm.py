#!/usr/bin/python
import os, sys
import shlex, subprocess
from datetime import datetime, date, time
#sys.path.append(os.path.abspath(os.path.curdir))
from collections import OrderedDict

Pmode = 8
ProdTag = "analysis_20170301_v0_p20170306_r1"
OutDir  = "/data/users/jengbou/histos"
WorkDir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD"

JobTime = datetime.now()
fTag = JobTime.strftime("%Y%m%d_%H%M%S")

sTags = {}
sTags["QCD76HT100to200"]=["11",3951,400]
sTags["QCD76HT200to300"]=["12",1262,400]
sTags["QCD76HT300to500"]=["13",1012,600]
sTags["QCD76HT500to700"]=["14",1375,300]
sTags["QCD76HT700to1000"]=["15",1325,300]


jobTags = OrderedDict(sorted(sTags.items(), key=lambda x: x[1]))

#########################################
# make sure OutDir is the same in main.cc
#########################################
condor_script_template = """
universe = vanilla
Executable = condor-executableMergeNoNorm.sh
+IsLocalJob = true
Should_transfer_files = NO
Requirements = TARGET.FileSystemDomain == "privnet"
Output = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stdout
Error  = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).stderr
Log    = %(OUTDIR)s/%(MYPREFIX)s/logs/%(SAMPLENAME)s_sce_$(cluster)_$(process).condor
Arguments = %(MYPREFIX)s %(SAMPLENAME)s_$(process) %(WORKDIR)s 0 0 %(IMODE)s 0 1 %(OUTDIR)s %(PMODE)s %(RNGL)s %(RNGH)s %(FIDX)s
Queue 1
"""
#########################################


for k,v in jobTags.items():
    print "Submitting jobs for [%-20s]"%k
    sTag_ = k
    Imode_ = v[0]
    numFiles_ = v[1]
    filesPerJob_ = v[2]

    dirname = "jobs/%s_MergeNoNorm_%s"%(sTag_,fTag)

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

    fidx = 0
    for fn in xrange(1,numFiles_+1):
        if fn%filesPerJob_ == 1:
            rUp = fn+filesPerJob_ - 1
            if rUp >= numFiles_: rUp = numFiles_

            kw["FIDX"] = "%i"%(fidx)
            kw["RNGL"] = "%i"%(fn)
            kw["RNGH"] = "%i"%(rUp)
            kw["SAMPLENAME"]    = "%s_MergeNoNorm_%ito%i"%(sTag_,fn,rUp)
            #print kw["RNGL"], kw["RNGH"]

            script_str = condor_script_template % kw
            f = open("%s/condor_jobs_%s_MergeNoNorm_%ito%i.jdl"%(dirname,sTag_,fn,rUp), 'w')
            f.write(script_str)
            f.close()

            condorcmd = "condor_submit %s/condor_jobs_%s_MergeNoNorm_%ito%i.jdl"%(dirname,sTag_,fn,rUp)
            print 'condorcmd: ', condorcmd
            print 'Executing condorcmd'
            #p=subprocess.Popen(condorcmd, shell=True)
            #p.wait()

            print "Histos output dir: %s/%s"%(OutDir,ProdTag)
            fidx+=1
