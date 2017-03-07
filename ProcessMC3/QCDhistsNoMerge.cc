#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"


void QCDhistsNoMerge(int nrange[2], std::string binname,std::string aaname, bool hasPre, bool blind, bool b16003, std::string bbname, bool crabformat)
{
    std::string inputfile;
    std::string outputfile;
  
    // opt
    float DHTcut=1000;
    float Dpt1cut=400;
    float Dpt2cut=200;
    float Dpt3cut=200;
    float Dpt4cut=100;
    float Dalphacut=0.04;
    float DmaxIPcut=0.4;
    float Djetacut = 2.;
    // dont forget there is a hidden cut nalmostemergin<4!!!!!!!!!!!!!!!!!
    int Dnemcut=2;
    int Dntrk1=0;
    // for alpha max scan

    // first make histograms for each file in each bin for the qcd sample
    std::cout<<"making histograms for each file in each bin"<<std::endl;
    mkdir((bbname+binname).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    for(int j=nrange[0]-1;j<nrange[1];j++) { //for each file for that bin
        std::cout << "Processing file # " << j << std::endl;
        inputfile=aaname+binname+"/"+binname+"_"+std::to_string(j+1)+"_0.histo.root";//condor format
        if (crabformat) inputfile=aaname+binname+"/ntuple_"+std::to_string(j+1)+".root";//crab format
        std::cout<<"input file is "<<inputfile<<std::endl;
        outputfile=bbname+binname+"/histos"+binname+"_"+std::to_string(j)+".root";
        std::cout<<"output file is "<<outputfile<<std::endl;
        int itmp;
        if(!b16003) {
            itmp = EMJselect(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,DmaxIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
        } else {
            itmp = EMJ16003(true,hasPre,inputfile.c_str(),outputfile.c_str());
        }
    }

    return;
}
