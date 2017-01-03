#include <iostream>
#include <iomanip>
#include <locale>

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
  cout<<"making histograms for each file in each bin"<<endl;
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin

      inputfile=aaname+binnames[i]+"/"+"ntuple_"+std::to_string(j)+".root";
      std::cout<<"input file is "<<inputfile<<endl;
      outputfile="histos"+binnames[i]+std::to_string(j)+".root";
      std::cout<<"output file is "<<outputfile<<endl;
      EMJselect(inputfile.c_str(),outputfile.c_str(),1000., 0.2,0.9,0.9,1);
    }
  }

  // now add up all the files for one bin
  cout<<" adding up histos within a bin"<<endl;
  vector<TH1F*> sumacount(nbin);
  vector<TH1F*> sumH_T(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+std::to_string(j)+".root";
      TFile *in = new TFile(inputfile.c_str());
      if(j==0) {
	sumacount[i] = static_cast<TH1F*>(in->Get("acount")->Clone());
	sumH_T[i] = static_cast<TH1F*>(in->Get("H_T")->Clone());
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get("acount")->Clone());
	sumacount[i]->Add(tmp);
	tmp = static_cast<TH1F*>(in->Get("H_T")->Clone());
	sumH_T[i]->Add(tmp);
      }
    }
  }
  // reweight to int lum
  cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<endl;

  //add the bins
  cout<<" adding bins"<<endl;
  TH1F* SUMacount=static_cast<TH1F*>((sumacount[0])->Clone());
  TH1F* SUMH_T=static_cast<TH1F*>((sumH_T[0])->Clone());
  for(int i=1;i<nbin;i++) {
    SUMacount->Add(sumacount[i]);
    SUMH_T->Add(sumH_T[i]);
  }

  // output histograms
cout<<"outputting histograms"<<endl;
  outputfile="SumHistos.root";
  TFile out(outputfile.c_str(),"RECREATE");

  SUMacount->Write();
  SUMH_T->Write();

  return;
}
