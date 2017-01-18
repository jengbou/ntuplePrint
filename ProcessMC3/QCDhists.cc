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


vector<float> Decode(int cutindex, int ncut,vector<int> nstep, vector<float> stepsize);

vector<int> EMJscan(const char* inputfilename,
		    float HTcutmin,int NHTcut, float HTcutSS,
		    float pt1cutmin, int Npt1cut, float pt1cutSS,
		    float pt2cutmin,  int Npt2cut,float pt2cutSS,
		    float pt3cutmin,  int Npt3cut,float pt3cutSS,
		    float pt4cutmin, int Npt4cut,float pt4cutSS,
		    int NemergingCutmin, int NNemergingCut,int NNemergingCutSS,
		    float jetacut,
                    float alphaMaxcut, float meanIPcut,
                    float NemfracCut,float CemfracCut,int ntrk1cut);
int EMJselect(bool otfile, const char* inputfilename,const char* outputfilename,
	      float HTcut, float pt1cut, float pt2cut,float pt3cut, float pt4cut, float jetacut,float alphaMaxcut, float meanIPcut, float NemfracCut,float CemfracCut,
	       int ntrk1cut, int NemergingCut	\
	       );
void  HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames);
TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames);


std::string bbname = "./";




//void QCDhists() 
void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta) 
{

  std::string inputfile;
  std::string outputfile;

  //for kine scan

  /*
  // YH default
  float DHTcut=1000;
  float Dpt1cut=400;
  float Dpt2cut=200;
  float Dpt3cut=125;
  float Dpt4cut=50;
  float Dalphacut=0.2;
  float DmeanIPcut=-1.;
  int Dnemcut=1;
  float Djetacut = 2.;
  // for alpha max scan
  const int ncutscan=5;
  */

    
  // opt
  float DHTcut=1000;
  float Dpt1cut=400;
  float Dpt2cut=200;
  float Dpt3cut=200;
  float Dpt4cut=100;
  float Dalphacut=0.04;
  float DmeanIPcut=-1;
  float Djetacut = 2.;
  int Dnemcut=2;
  // for alpha max scan
  const int ncutscan=5;
  //const int ncutscan=1;
  

  // first make histograms for each file in each bin for the qcd sample

  std::cout<<"making histograms for each file in each bin"<<std::endl;
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;
      outputfile=bbname+"histos"+binnames[i]+"_"+std::to_string(j)+".root";
      std::cout<<"output file is "<<outputfile<<std::endl;
      int itmp = EMJselect(true,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,DmeanIPcut,0.9,0.9,0,Dnemcut);
    }
  }


  // do some cut optimization on cuts not related to choosing the emerging jets
  const int nkincut=6;
  vector<int> nstep {10,10,10,10,10,3};
  //vector<int> nstep {2,2,2,2,2,2};
  // ht pt1 pt2 pt3 pt4 nemerging
  vector<float> cutmin {1000.,400.,200.,200.,100.,0};
  vector<float> ss {50,10,10,10,10,1};


  int iicut =nstep[0];
  for (int hh=1;hh<nkincut;hh++) iicut*=nstep[hh];
  vector < vector <int> > nnpass(iicut,vector<int>(nbin,0));
  if(dooptk==1) {
  for(int i=0;i<nbin;i++) {  // for each bin
    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;


      vector<int> npass = EMJscan(inputfile.c_str(),
			      cutmin[0],nstep[0],ss[0],
                              cutmin[1],nstep[1],ss[1],
                              cutmin[2],nstep[2],ss[2],
                              cutmin[3],nstep[3],ss[3], 
                              cutmin[4],nstep[4],ss[4],
			      cutmin[5],nstep[5],ss[5],
				  Djetacut,
				  0.04,-1.,0.9,0.9,0);
      for(int tt=0;tt<iicut;tt++) {
	nnpass[tt][i]=nnpass[tt][i]+npass[tt];
      }
    }
  }}


  //vector<float> decode(6);
  //decode = Decode(icut,6,nstep,stepsize);

  //do some cut optimization on alpha max


  //const int ncutscan=2;
  float acut=0.3;
  //  int ipass[ncutscan][nbin];
  vector < vector <int> > ipass(ncutscan, vector<int>(nbin,0));
  if(doopta==1) {
  for(int k=0;k<ncutscan;k++) {
    float acut2=(acut/(ncutscan))*(k);
    std::cout<<" cut value is "<<acut2<<std::endl;
    for(int i=0;i<nbin;i++) ipass[k][i]=0;
    for(int i=0;i<nbin;i++) {  // for each bin
      for(int j=0;j<nfiles[i];j++) { //for each file for that bin
	std::cout<<"k i j="<<k<<" "<<i<<" "<<j<<std::endl;
      inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";
      std::cout<<"input file is "<<inputfile<<std::endl;
      int iii = EMJselect(false,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,acut2,0.9,0.9,0,Dnemcut);
	ipass[k][i]+=iii;
	std::cout<<" iii ipass  is "<<iii<<" "<<ipass[k][i]<<std::endl;
      }
    }
  }}

  // get normalization
  //  double norm[nbin];
  vector<double> norm(nbin);
  HistNorm(norm,nbin,xsec,nfiles,binnames);
  for(int i=0;i<nbin;i++) {
    std::cout<<"total number events in bin "<<i<<" is "<<norm[i]<<std::endl;
  }
  TH1F* normhst = new TH1F("normhst","counts pretrigger by bin",nbin,0.,nbin);
  for(int i=0;i<nbin;i++){
    normhst->AddBinContent(i+1,norm[i]);
  }


  //make and  output summed and renormalized histograms
  std::cout<<"normalizing histograms"<<std::endl;
  const int nhist=39;
  std::vector<TH1F*> vv(nhist);
  std::string histnames[nhist]={
"count","acount","hjetcut","hjetchf","h_nemg",
"hnjet","hpt","heta","heta2","halpha",
"H_T","H_T2","hpt1","hpt2","hpt3",
"hpt4","hbcut_ntrkpt1","hacut_ntrkpt1","hbcut_nef","hacut_nef",
"hbcut_cef","hacut_cef","hbcut_alphamax","hacut_alphamax","hHTnm1",
"hpt1nm1","hpt2nm1","hpt3nm1","hpt4nm1","halphanm1",
"hnemnm1","hipXYEJ","hipXYnEJ","htvwEJ","htvw",
"hipXYSigEJ","hipXYSignEJ","hmeanipXYEJ","hmeanipXYnEJ"
};
  vector<double> outnorm(nbin);
  for(int i=0;i<nhist;i++) {
    std::cout<<" enering Histman with i = "<<i<<std::endl;
    vv[i]=HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,binnames);
  }
  // output total event count
  std::cout<<" initial event count before and after norm is"<<std::endl;
  double ttotal=0;
  for(int i=0;i<nbin;i++) {
    std::cout<<" bin "<<i<<" norm "<<norm[i]<<" times outnorm is "<<norm[i]*outnorm[i]<<std::endl;
    ttotal = ttotal + norm[i]*outnorm[i];
  }
  std::cout<<"total is "<<ttotal<<std::endl;;


  // normalize cut scan and sum bins
  std::cout<<"normalizing kinematic cut stuff"<<std::endl;
  std::cout<<"iicut "<<iicut<<std::endl;
  double ffpass[iicut];
  for(int i=0;i<iicut;i++) ffpass[i]=0;
  for(int k=0;k<iicut;k++) {
    for(int i=0;i<nbin;i++) {
      ffpass[k]+=nnpass[k][i]*outnorm[i];
    }
    std::cout<<" output kinematic cut scan "<<k<<" "<<ffpass[k]<<std::endl;
  }
  TH1F* kcutscan = new TH1F("kcutscan","n pass versus cut kin cuts",iicut,0.,iicut);
  for(int i=0;i<iicut;i++){
    kcutscan->AddBinContent(i+1,ffpass[i]);
  }


  // normalize cut scan and sum bins
  std::cout<<"normalizing alpha scan stuff"<<std::endl;
  double fpass[ncutscan];
  for(int i=0;i<ncutscan;i++) fpass[i]=0;
  for(int k=0;k<ncutscan;k++) {
    for(int i=0;i<nbin;i++) {
      fpass[k]+=ipass[k][i]*outnorm[i];
    }
    std::cout<<" output alphamax scan "<<k<<" "<<fpass[k]<<std::endl;
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

void  HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames) {

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
