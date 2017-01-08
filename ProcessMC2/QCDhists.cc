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

vector<int> EMJscan(const char* inputfilename,
		    float HTcutmin,int NHTcut, float HTcutSS,
		    float pt1cutmin, int Npt1cut, float pt1cutSS,
		    float pt2cutmin,  int Npt2cut,float pt2cutSS,
		    float pt3cutmin,  int Npt3cut,float pt3cutSS,
		    float pt4cutmin, int Npt4cut,float pt4cutSS,
		    int NemergingCutmin, int NNemergingCut,
                    float alphaMaxcut, float NemfracCut,float CemfracCut,int ntrk1cut);
int EMJselect(bool otfile, const char* inputfilename,const char* outputfilename,
float HTcut, float pt1cut, float pt2cut,float pt3cut, float pt4cut,
float alphaMaxcut, float NemfracCut,float CemfracCut,
	       int ntrk1cut, int NemergingCut	\
	       );
void  HistNorm(float goalintlum,vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames);
TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames);


// need to update this section below
/*
float goalintlum=20; // fb-1
const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb
const int nfiles[nbin]={2,2,2,2,2};
//const int nfiles[nbin]={138,133,50,40,23};
const std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};
std::string aaname = "/data/users/eno/outputQCD/";
*/



std::string bbname = "./";




//void QCDhists() 
void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname) 
{

  std::string inputfile;
  std::string outputfile;


  // first make histograms for each file in each bin for the qcd sample

  std::cout<<"making histograms for each file in each bin"<<std::endl;
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;
      outputfile=bbname+"histos"+binnames[i]+"_"+std::to_string(j)+".root";
      std::cout<<"output file is "<<outputfile<<std::endl;
      int itmp = EMJselect(true,inputfile.c_str(),outputfile.c_str(),1000., 400.,300.,170.,100.,0.2,0.9,0.9,0,2);
    }
  }


  // do some cut optimization on cuts not related to choosing the emerging jets
  const int nkincut=6;
  int nstep[nkincut]={3,3,3,3,3,3};
  int iicut =nstep[0];
  for (int hh=1;hh<nkincut;hh++) iicut*=nstep[hh];
  vector < vector <int> > nnpass(iicut,vector<int>(nbin,0));
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;

      vector<int> npass = EMJscan(inputfile.c_str(),
			      1000.,nstep[0],100,
                              400,nstep[0],50,
                              200,nstep[0],50,
                              120,nstep[0],50, 
                              50,nstep[0],50,
                              0,nstep[0],
                              0.04,0.9,0.9,0);
      for(int tt=0;tt<iicut;tt++) {
	nnpass[tt][i]=nnpass[tt][i]+npass[tt];
      }
    }
  }

  //do some cut optimization on alpha max

  const int ncutscan=5;
  float acut=0.2;
  //  int ipass[ncutscan][nbin];
  vector < vector <int> > ipass(ncutscan, vector<int>(nbin,0));
  for(int k=0;k<ncutscan;k++) {
    float acut2=(acut/(ncutscan))*(k+1);
    std::cout<<" cut value is "<<acut2<<std::endl;
    for(int i=0;i<nbin;i++) ipass[k][i]=0;
    for(int i=0;i<nbin;i++) {  // for each bin
      for(int j=0;j<nfiles[i];j++) { //for each file for that bin
	std::cout<<"k i j="<<k<<" "<<i<<" "<<j<<std::endl;
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;
	int iii = EMJselect(false,inputfile.c_str(),outputfile.c_str(),1000., 400.,300.,170.,100.,acut2,0.9,0.9,0,2);
	ipass[k][i]+=iii;
	std::cout<<" iii ipass  is "<<iii<<" "<<ipass[k][i]<<std::endl;
      }
    }
  }

  // get normalization
  //  double norm[nbin];
  vector<double> norm(nbin);
  HistNorm(goalintlum,norm,nbin,xsec,nfiles,binnames);
  for(int i=0;i<nbin;i++) {
    std::cout<<"total number events in bin "<<i<<" is "<<norm[i]<<std::endl;
  }
  TH1F* normhst = new TH1F("normhst","counts pretrigger by bin",nbin,0.,nbin);
  for(int i=0;i<nbin;i++){
    normhst->AddBinContent(i+1,norm[i]);
  }


  //make and  output summed and renormalized histograms
  std::cout<<"doing the stuff"<<std::endl;
  const int nhist=19;
  std::vector<TH1F*> vv(nhist);
  std::string histnames[nhist]={"acount","hjetcut","hjetchf","h_nemg","hnjet","hpt","heta","heta2","halpha","H_T","H_T2","hbcut_ntrkpt1","hacut_ntrkpt1","hbcut_nef","hacut_nef","hbcut_cef","hacut_cef","hbcut_alphamax","hacut_alphamax"};
  vector<double> outnorm(nbin);
  for(int i=0;i<nhist;i++) {
    vv[i]=HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,binnames);
  }


  // normalize cut scan and sum bins
  double ffpass[iicut];
  for(int i=0;i<iicut;i++) ffpass[i]=0;
  for(int k=0;k<iicut;k++) {
    for(int i=0;i<nbin;i++) {
      ffpass[k]+=nnpass[k][i]*outnorm[i];
    }
    std::cout<<" ffpass"<<k<<" "<<ffpass[k]<<std::endl;
  }
  TH1F* kcutscan = new TH1F("kcutscan","n pass versus cut kin cuts",iicut,0.,iicut);
  for(int i=0;i<iicut;i++){
    kcutscan->AddBinContent(i+1,ffpass[i]);
  }


  // normalize cut scan and sum bins
  double fpass[ncutscan];
  for(int i=0;i<ncutscan;i++) fpass[i]=0;
  for(int k=0;k<ncutscan;k++) {
    for(int i=0;i<nbin;i++) {
      fpass[k]+=ipass[k][i]*outnorm[i];
    }
    std::cout<<" fpass"<<k<<" "<<fpass[k]<<std::endl;
  }
  TH1F* cutscan = new TH1F("cutscan","n pass versus cut",ncutscan,0.,ncutscan);
  for(int i=0;i<ncutscan;i++){
    cutscan->AddBinContent(i+1,fpass[i]);
  }



  std::cout<<"outputting histograms"<<std::endl;
  outputfile=bbname+ohname;
  TFile out(outputfile.c_str(),"RECREATE");
  normhst->Write();
  cutscan->Write();
  kcutscan->Write();
  for(int i=0;i<nhist;i++) {
    vv[i]->Write();
  }


  return;
}



TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames) {

  std::string inputfile;


  // now add up all the files for one bin
  std::cout<<" adding up histos within a bin"<<std::endl;
  vector<TH1F> sum(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+"_"+std::to_string(j)+".root";
      TFile* in = new TFile(inputfile.c_str());
      if(j==0) {
	sum[i] = *(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone()));
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
	sum[i].Add(tmp);
      }
      in->Close();
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
    outnorm[i] = goalintlum/fileLum;
    std::cout<<" scaling by a factor of "<<outnorm[i]<<std::endl;
    sum[i].Scale(outnorm[i]);
  }


  //add the bins
  std::cout<<" adding bins"<<std::endl;
  TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone());
  for(int i=1;i<nbin;i++) {
    SUM->Add(&sum[i]);
  }


  return SUM;
}

void  HistNorm(float goalintlum,vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames) {

  std::cout<<"entering HistNorm"<<std::endl; 

  std::string inputfile;
  TFile * in;

  // now add up all the files for one bin
  vector<TH1F> sum(nbin);
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile="histos"+binnames[i]+"_"+std::to_string(j)+".root";
      std::cout<<i<<" "<<j<<" "<<inputfile<<std::endl;
      in = new TFile(inputfile.c_str());
      if(j==0) {
	sum[i] = *(static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone()));
      } else {
	TH1F* tmp = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
	sum[i].Add(tmp);
      }
      in->Close();
    }
  }

  // reweight to int lum

  for(int i=0;i<nbin;i++) {
    // get total number of events before filter
    norm[i] = sum[i].GetBinContent(2);
    std::cout<<"norm "<<i<<" "<<norm[i]<<std::endl;
  }


  return;
}
