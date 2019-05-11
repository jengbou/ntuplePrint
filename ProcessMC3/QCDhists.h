#ifndef QCDhists_h
#define QCDhists_h

#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMatrixD.h>

#include "vector"
using std::vector;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sys/stat.h>


class QCDhists {
 public:

  TFile* frFile;
  TH1F* hfrB;
  TH1F* hfrL;

  void initialize();
  ~QCDhists(){}
  QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* samplenames,std::string indir,std::string ohname, int dooptk, int doopta,bool hasPre,bool norm, bool blind, bool b16003,std::string,bool,bool);

  Int_t GetNtrkBin( const int jet_nTrack );
  vector<float> Decode(int cutindex, int ncut,vector<int> nstep, vector<float> stepsize);

  void  HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* samplenames,std::string dirname);
  TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname);

  TH1F* HistMerge(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname);

  TH2F* HistMan2(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname);

  double fakerate(double,double,int,int varType=1);//pt,eta,ntrk,type
  double fakerateTP(double,double,int,int varType=1);//pt,eta,ntrk,type
  double fakerateF(double,double,int,int varType=1,int flav=1,std::string cutset="0");//pt,eta,ntrk,type,flav(5:b; 8: gbb; else:g,u,d,c,s)
  double fakerateFD(double,double,int,int varType=1,int flav=1,std::string cutset="0");
  double ntrkrewgt(int nTrk);
  double nGJrewgt(int nGoodJet);
  double UnfoldWgt(int bUnfType, int nbtrue, int nbtagged, int ptType=1);
  double UnfoldWgtD(int bUnfType, int nbtrue, int nbtagged, int ptType=1);
  //double UnfoldWgtPtDep(int bUnfType, int nbtrue, int nbtagged, vector<float> *jetpt, bool isData=false);
  bool UnfoldWgtPtDep(TMatrixD&, int bUnfType, vector<float> *jetpt, bool isData);
  double effDeepCSV(int bUnfType, double pt, int flav, bool isData=false);
  double effCSVv2(int bUnfType, double pt, int flav, bool isData=false);
  double effBTag(int bUnfType, double pt, int flav, bool isData=false);
  double effBTagPara(int bUnfType, double pt, int flav, bool isData);
  double frWeight(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<vector<float> >*track_pt, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &track_pt, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightF1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flav, int njetscut=-1, double jptcut=100.0, int varType=1, std::string cutset="0");
  double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightFT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flav, int njetscut=-1, double jptcut=100.0, int varType=1, std::string cutset="0");
  double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightFT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flav, int njetscut=-1, double jptcut=100.0, int varType=1, std::string cutset="0");
  double frWeightFT12(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, double jptcut, int bUnfType, int nbTagged, int ptType=1, int varType=1, bool isData=false, std::string cutset="0");
  void frWeightUFT1(double (&frwgts)[5],vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int* flav, int varType=1, bool isData=false, std::string cutset="0");
  double frWeightT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightT21(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightFT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flav, int njetscut=-1, double jptcut=100.0, int varType=1, std::string cutset="0");
  double frWeightFT22(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, double jptcut, int bUnfType, int nbTagged, int ptType=1, int varType=1, bool isData=false, std::string cutset="0");
  double frWeightUFT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &jetIdx, vector<int> &ntrack, vector<int> &flav, int varType=1, bool isData=false, std::string cutset="0");
  void frWeightUFT22(double (&frwgts)[25],vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int* flav, int varType=1, bool isData=false, std::string cutset="0");
  void frWeightUFT23(double (&frwgts)[7],vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int* flav, int varType=1, bool isData=false, std::string cutset="0");
  double frWeightT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
  double frWeightFT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flav, int njetscut=-1, double jptcut=100.0, int varType=1, std::string cutset="0");

  double frWeight4(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, double jptcut=100.0, int varType=1);

  double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight);
  double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight, vector<float> &track_ref_zs, float pv_z, float pilecut=5000.0);
  double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_iPXYSigs);
  double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_iPXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut=5000.0);
  double GetAlpha3Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_iPXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut=5000.0, float sigzcut=4.0);
  void fillFRPlots(std::string cutname, char * hnames[9],
                   vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                   vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                   vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                   vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                   double *varwgt, double scale,double ncomb=1);

  void fillFRPlots2(std::string cutname, char * hnames[9],
                    vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                    vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                    vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                    vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                    double *varwgt, double scale,double ncomb=1);

  void fillFRPlots22(std::string cutname, char * hnames[9],
                     vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                     vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                     vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                     vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                     double *varwgt, double scale,double ncomb=1);

 private:

};

#endif
