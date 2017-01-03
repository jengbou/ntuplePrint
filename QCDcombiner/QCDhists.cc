#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
using std::vector;


#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void EMJselect(const char* inputfilename,const char* outputfilename,
	       float HTcut, float alphaMaxcut, float NemfracCut,float CemfracCut,int NemergingCut\
	       );
void  HistNorm(double* norm);
TH1F* HistMan(std::string thisHIST,double* histnorm);


// need to update this section below

float goalintlum=20; // fb-1
const int nbin=2;
float xsec[nbin]={29370000,6524000};
const int nfiles[nbin]={1,1};
const std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000"};
std::string aaname = "/data/users/eno/outputQCDdummy/";


//const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
//float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb



void QCDhists() {

  std::string inputfile;
  std::string outputfile;

  // first make histograms for each file in each bin for the qcd sample
  std::cout<<"making histograms for each file in each bin"<<std::endl;
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin

      inputfile=aaname+binnames[i]+"/"+"ntuple_"+std::to_string(j)+".root";
      std::cout<<"input file is "<<inputfile<<std::endl;
      outputfile="histos"+binnames[i]+std::to_string(j)+".root";
      std::cout<<"output file is "<<outputfile<<std::endl;
      EMJselect(inputfile.c_str(),outputfile.c_str(),1000., 0.2,0.9,0.9,1);
    }
  }

  // get normalization
  double norm[nbin];
  HistNorm(norm);
  for(int i=0;i<nbin;i++) {
    std::cout<<"total number events in bin "<<i<<" is "<<norm[i]<<std::endl;
  }
  


  //make and  output summed and renormalized histograms
  std::cout<<"doing the stuff"<<std::endl;
  TH1F* SUMH_T = HistMan("H_T",norm);
  TH1F* SUMhpt = HistMan("hpt",norm);
  TH1F* SUMheta = HistMan("heta",norm);


  std::cout<<"outputting histograms"<<std::endl;
  outputfile="SumHistos.root";
  TFile out(outputfile.c_str(),"RECREATE");

  SUMH_T->Write();
  SUMhpt->Write();
  SUMheta->Write();

  return;
}



TH1F* HistMan(std::string thisHIST,double* norm) {

  std::string inputfile;
  std::string outputfile;


  // now add up all the files for one bin
  std::cout<<" adding up histos within a bin"<<std::endl;
  vector<TH1F*> sum(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+std::to_string(j)+".root";
      TFile *in = new TFile(inputfile.c_str());
      if(j==0) {
	sum[i] = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
	sum[i]->Add(tmp);
      }
    }
  }

  // reweight to int lum
  std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
  for(int i=0;i<nbin;i++) {
    // get total number of events before filter
    float ntotal = norm[i];
    std::cout<<" for bin "<<i<<" number of pretrigger events is "<<ntotal<<std::endl;
    float fileLum= ntotal/xsec[i];
    std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
    float norm = goalintlum/fileLum;
    std::cout<<" scaling by a factor of "<<norm<<std::endl;
    sum[i]->Scale(norm);
  }


  //add the bins
  std::cout<<" adding bins"<<std::endl;
  TH1F* SUM=static_cast<TH1F*>((sum[0])->Clone());
  for(int i=1;i<nbin;i++) {
    SUM->Add(sum[i]);
  }

  return SUM;
}

void  HistNorm(double* norm) {

  std::string inputfile;


  // now add up all the files for one bin
  vector<TH1F*> sum(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+std::to_string(j)+".root";
      TFile *in = new TFile(inputfile.c_str());
      if(j==0) {
	sum[i] = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
	sum[i]->Add(tmp);
      }
    }
  }

  // reweight to int lum

  for(int i=0;i<nbin;i++) {
    // get total number of events before filter
    norm[i] = sum[i]->GetBinContent(2);
  }

  return;
}
