#ifndef MergeHists_h
#define MergeHists_h

#include <iostream>
#include <iomanip>
#include <locale>

#include "vector"
using std::vector;

#include <sys/stat.h>


void MergeHists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* samplenames,std::string indir,std::string ohname,bool donorm, std::string outdir="./");


#endif
