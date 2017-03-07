#ifndef MergeHistsN_h
#define MergeHistsN_h


void MergeHistsN(float goalintlum, float xsec, int nfiles, std::string binname, std::string aaname, std::string ohname, bool donorm, std::string bbname);
void  HistNormN(double &norm, float xsec, int nfiles, std::string binname, std::string bbname);
TH1F* HistManN(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname);
TH2F* HistMan2N(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname);


#endif
