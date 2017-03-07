#ifndef QCDhistsNoMerge_h
#define QCDhistsNoMerge_h

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
#include "QCDhists.h"


void QCDhistsNoMerge(int nrange[2], std::string binnames, std::string aaname, bool hasPre, bool blind, bool b16003, std::string bbname, bool crabformat);

#endif
