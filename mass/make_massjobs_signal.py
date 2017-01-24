# for some reason now you have to run it twice.  
# it seems the .txt1 files are not quite there when the cat commands comes along
medmass = [400, 600, 800, 1000, 1500, 2000]
dpmass = [1, 2, 5, 10]
dplt = ["0p001", "0p1", "1", "5", "150", "300"]

import subprocess

#medmass = [800]
#dpmass = [1]
#dplt = ["0p001"]


exearea = "/data/users/eno/em6/CMSSW_7_6_3/src/EmergingJetAnalysis/"
thisarea = "/data/users/eno/DARKJOBS/"

f = open("massjobs.sh",'w')


for i in medmass:
  for j in dpmass:
    for k in dplt:
#      print "HAHAHAHAHA"
      bname = "EmergingJets_mass_X_d_"+str(i)+"_mass_pi_d_"+str(j)+"_tau_pi_d_"+k
      scname = "Xd"+str(i)+"masspid"+str(j)+"taupid"+k
      aname = "/mnt/hadoop/cms//store/user/abelloni/EmJetMC/AODSIM-v1/"+bname+"_TuneCUETP8M1_13TeV_pythia8Mod/AODSIM-v1/"
#      print aname
#      command = "find "+aname+" -name aodsim_1.root -print | grep -FzZ '000/aodsim_1.root'"
      command = "find "+aname+" -name aodsim_1.root -print | grep '000/aodsim_1.root'"
#      command = "find "+aname+" -name aodsim_1.root -print"
#      print command
      p = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
      (namek, err) = p.communicate();
#      print "namek is "+namek
      namev = namek.replace("aodsim_1.root"+'\n',"")
#      print "namev is "+namev
      command2 = "mkdir /data/users/eno/outputQCD/"+scname
      q = subprocess.Popen(command2,stdout=subprocess.PIPE, shell=True)
      command3= "ls "+namev+"*.root > /data/users/eno/DARKJOBS/"+scname+".txt1"
#      print command3
      q = subprocess.Popen(command3,stdout=subprocess.PIPE, shell=True)
      command4="cat /data/users/eno/DARKJOBS/"+scname+".txt1"+" | sed -e s~/mnt~file:/mnt~ > /data/users/eno/DARKJOBS/"+scname+".txt"
#      print command4
      q = subprocess.Popen(command4,stdout=subprocess.PIPE, shell=True)

      name = scname
      hostarea = "/data/users/eno/outputQCD/"+name+"/"
      jdlfile = open("condor-jobs-"+name+".jdl","w")
      jdlfile.write("universe = vanilla"+'\n')
      jdlfile.write("Executable = condor-executable.sh"+'\n')
      jdlfile.write("should_transfer_files = NO"+'\n')
      jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
      jdlfile.write("Output = "+hostarea+name+"_sce_$(cluster)_$(process).stdout"+'\n')
      jdlfile.write("Error = "+hostarea+name+"_sce_$(cluster)_$(process).stderr"+'\n')
      jdlfile.write("Log = "+hostarea+name+"_sce_$(cluster)_$(process).condor"+'\n')
      jdlfile.write("Arguments = "+name+" $(process) "+hostarea+" "+exearea+" Configuration/test/condor_cfg.py \
1234567 100000 "+thisarea+name+".txt data=0"+'\n')
      jdlfile.write("Queue 1"+'\n')
      jdlfile.close()
      f.write("condor_submit condor-jobs-"+name+'.jdl'+'\n')




    


f.close()
