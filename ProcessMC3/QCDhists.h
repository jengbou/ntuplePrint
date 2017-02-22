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

TH2F* HistMan2(float goalintlum,std::string thisHIST,vector<double>& histnorm, vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname);

void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta,bool hasPre,bool norm, bool blind, bool b16003,std::string,bool);


#endif
