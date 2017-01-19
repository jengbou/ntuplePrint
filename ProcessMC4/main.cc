
#include <iostream>
#include <string>
#include <map>

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta) ;



int main(int argc, char *argv[])
{ 
  int dooptk =*(argv[1])-'0';
  int doopta =*(argv[2])-'0';
  int imode=*(argv[3])-'0';

  if(dooptk==0) {
    std::cout<<"not doing kinematic optimization"<<std::endl;
  } else if(dooptk==1) {
    std::cout<<"doing kinematic optimization"<<std::endl;
  } else {
    std::cout<<"invalid choice"<<std::endl;
  }

  if(doopta==0) {
    std::cout<<"not doing alphamax optimization"<<std::endl;
  } else if(doopta==1) {
    std::cout<<"doing alphamax optimization"<<std::endl;
  } else {
    std::cout<<"invalid choice"<<std::endl;
  }

  if(imode==0) {
    std::cout<<"doing background"<<std::endl;
  } else if(imode==1) {
    std::cout<<"doing modelA"<<std::endl;
  } else if(imode==2) {
    std::cout<<"doing modelB"<<std::endl;
  } else if(imode==3) {
    std::cout<<"doing debug sample"<<std::endl;
  } else {
    std::cout<<"invalid choice"<<std::endl;
  }


float goalintlum=20; // fb-1                                                                                        
 std::string aaname = "/data/users/eno/outputQCD/";  // area containing subdirectors with YHS's ntuples

// for background 
//const int nbin=2; // 500-700,700-1000,1000-1500,1500-2000,200toInf
//float xsec[nbin]={29370000,6524000}; // fb 
//int nfiles[nbin]={1,1};                                                                                
//std::string binnames[nbin]={"QCD_HT1500to2000","QCD_HT2000toInf"};

const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb 
int nfiles[nbin]={138,133,50,40,23};
std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};


// for signal models A.  mediat mass is 1000
const int anbin=1; 
float axsec[nbin]={18.45}; // fb 
int anfiles[nbin]={50}; 
//int anfiles[nbin]={5}; 
std::string abinnames[nbin]={"modelA"};

// for signal models B.  mediat mass is 1000
const int bnbin=1; 
float bxsec[nbin]={18.45}; // fb 
int bnfiles[nbin]={50}; 
//int bnfiles[nbin]={5}; 
std::string bbinnames[nbin]={"modelB"};


// for debugging
const int dnbin=1; 
float dxsec[nbin]={18.45}; // fb 
 int dnfiles[nbin]={1}; 
 if(imode==3) aaname = "/home/eno/em5/EmergingJetAnalysis/";
std::string dbinnames[nbin]={"tmpStore"};



 if(imode==0) {
   QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root",dooptk,doopta);
 } else if (imode==1) {
   QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,aaname,"SumHistsModelA.root",dooptk,doopta);
 } else if (imode==2) {
   QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname,"SumHistsModelB.root",dooptk,doopta);
 } else if (imode==3) {
   QCDhists(goalintlum,dnbin,dxsec,dnfiles,dbinnames,aaname,"SumHistsDebug.root",0,0);
 }
}