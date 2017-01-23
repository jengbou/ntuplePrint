name="QCD74_HT500to700"


nstart=input("enter starting file ")
ndo=input("enter number to do ")

filename = r'/mnt/hadoop/cms/store/DARKQCD/'+name+'.txt'
print "reading from file "+filename
f = open(filename)

outname = open("masscopy.sh",'w')


outdir='/mnt/hadoop/cms/store/DARKQCD/'+name+'/'

a =1
line = f.readline()
while line:
    if( a>=nstart ):
        if(a<(nstart+ndo)):
            print a
            line=line.rstrip('\r\n')
            print line
            comd='xrdcp root://cmsxrootd.fnal.gov/'+line+' '+outdir+"."+'\n'
            outname.write(comd)
    line=f.readline()
    a += 1 # This variable will change the file number

f.close()
outname.close()
