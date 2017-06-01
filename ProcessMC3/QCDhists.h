#ifndef QCDhists_h
#define QCDhists_h

#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
using std::vector;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sys/stat.h>

vector<float> Decode(int cutindex, int ncut,vector<int> nstep, vector<float> stepsize);

void  HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames,std::string bbname);
TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname);

TH1F* HistMerge(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname);

TH2F* HistMan2(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname);

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta,bool hasPre,bool norm, bool blind, bool b16003,std::string,bool);

double fakerate(double,double,int,int varType=1);//pt,eta,ntrk,type
double fakerateTP(double,double,int,int varType=1);//pt,eta,ntrk,type
double ntrkrewgt(int nTrk);
double nGJrewgt(int nGoodJet);
double frWeight(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<vector<float> >*track_pt, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &track_pt, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT21(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeightT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut=-1, double jptcut=100.0, int varType=1);
double frWeight4(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, double jptcut=100.0, int varType=1);

double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight);
double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight, vector<float> &track_ref_zs, float pv_z, float pilecut=5000.0);
double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_iPXYSigs);
double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_iPXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut=5000.0);

#endif
