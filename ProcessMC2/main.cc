
#include <iostream>
#include <string>
#include <map>

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname) ;



int main(int argc, char *argv[])
{ 
  int imode=*(argv[1])-'0';
  if(imode==0) {
    std::cout<<"doing background"<<std::endl;
  } else if(imode==1) {
    std::cout<<"doing signal"<<std::endl;
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


// for signal models A and  B.  mediat mass is 1000
const int snbin=1; 
float sxsec[nbin]={18.45}; // fb 
//int snfiles[nbin]={50}; 
int snfiles[nbin]={30}; 
std::string sbinnames[nbin]={"modelB"};



 if(imode==0) {
   QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root");
 } else if (imode==1) {
   QCDhists(goalintlum,snbin,sxsec,snfiles,sbinnames,aaname,"SumHistsSignal.root");
 }
}
