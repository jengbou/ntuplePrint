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


float goalintlum=20; // fb-1


//const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
//float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb

const int nbin=2;
float xsec[nbin]={29370000,6524000};
const int nfiles[nbin]={1,1};
const std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000"};
std::string aaname = "/data/users/eno/outputQCDdummy/";




void QCDhists() {


  std::string inputfile;
  std::string outputfile;
  //const char* aa1=inputfile.c_str();
  //const char* bb1=outputfile.c_str();

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

  // now add up all the files for one bin
  std::cout<<" adding up histos within a bin"<<std::endl;
  vector<TH1F*> sumeventCountPreTrigger(nbin);
  vector<TH1F*> sumacount(nbin);
  vector<TH1F*> sumH_T(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+std::to_string(j)+".root";
      TFile *in = new TFile(inputfile.c_str());
      if(j==0) {
	sumeventCountPreTrigger[i] = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
	sumacount[i] = static_cast<TH1F*>(in->Get("acount")->Clone());
	sumH_T[i] = static_cast<TH1F*>(in->Get("H_T")->Clone());
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
	sumeventCountPreTrigger[i]->Add(tmp);
	tmp = static_cast<TH1F*>(in->Get("acount")->Clone());
	sumacount[i]->Add(tmp);
	tmp = static_cast<TH1F*>(in->Get("H_T")->Clone());
	sumH_T[i]->Add(tmp);
      }
    }
  }

  // reweight to int lum
  std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
  for(int i=0;i<nbin;i++) {
    // get total number of events before filter
    float ntotal = sumeventCountPreTrigger[i]->GetBinContent(2);
    std::cout<<" for bin "<<i<<" number of pretrigger events is "<<ntotal<<std::endl;
    float fileLum= ntotal/xsec[i];
    float norm = goalintlum/fileLum;
    sumeventCountPreTrigger[i]->Scale(norm);
    sumacount[i]->Scale(norm);
    sumH_T[i]->Scale(norm);
  }


  //add the bins
  std::cout<<" adding bins"<<std::endl;
  TH1F* SUMeventCountPreTrigger=static_cast<TH1F*>((sumeventCountPreTrigger[0])->Clone());
  TH1F* SUMacount=static_cast<TH1F*>((sumacount[0])->Clone());
  TH1F* SUMH_T=static_cast<TH1F*>((sumH_T[0])->Clone());
  for(int i=1;i<nbin;i++) {
    SUMeventCountPreTrigger->Add(sumeventCountPreTrigger[i]);
    SUMacount->Add(sumacount[i]);
    SUMH_T->Add(sumH_T[i]);
  }

  // output histograms
  std::cout<<"outputting histograms"<<std::endl;
  outputfile="SumHistos.root";
  TFile out(outputfile.c_str(),"RECREATE");

  SUMeventCountPreTrigger->Write();
  SUMacount->Write();
  SUMH_T->Write();

  return;
}
