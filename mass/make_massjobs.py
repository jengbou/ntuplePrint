name="QCD_HT500to700"

nstart=input("enter starting file ")
ndo=input("enter number to do ")


f = open("massjobs.sh",'w')

i=0
while (i<ndo):
    j=i+nstart
    f.write("condor_submit condor-jobs-"+name+"_"+str(j)+'.jdl'+'\n')
    i=i+1


f.close()
