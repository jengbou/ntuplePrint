#ifndef MergeHistsNorm_h
#define MergeHistsNorm_h


void MergeHistsNorm(float goalintlum, float xsec, int nfiles, std::string samplename, std::string ohname, bool donorm, std::string dirname);
void  HistNormNew(double &norm, float xsec, int nfiles, std::string samplename, std::string dirname);
TH1F* HistManNorm(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname);
TH2F* HistMan2Norm(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname);


#endif
