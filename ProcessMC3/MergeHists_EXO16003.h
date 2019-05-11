#ifndef MergeHists_EXO16003_h
#define MergeHists_EXO16003_h

#include <iostream>
#include <iomanip>
#include <locale>

#include "vector"
using std::vector;

#include <sys/stat.h>


void MergeHists_EXO16003(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname,bool donorm, std::string bbname="./");


#endif
