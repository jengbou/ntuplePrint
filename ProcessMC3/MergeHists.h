#ifndef MergeHists_h
#define MergeHists_h

#include <iostream>
#include <iomanip>
#include <locale>

#include "vector"
using std::vector;

#include <sys/stat.h>


void MergeHists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname,bool donorm, std::string bbname="./");


#endif
