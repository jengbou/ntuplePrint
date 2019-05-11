#ifndef MergeHistsN_h
#define MergeHistsN_h


void MergeHistsN(float goalintlum, float xsec, int nfiles, std::string samplename, std::string ohname, bool donorm, std::string dirname);
void  HistNormN(double &norm, float xsec, int nfiles, std::string samplename, std::string dirname);
TH1F* HistManN(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname);
TH2F* HistMan2N(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname);


#endif
