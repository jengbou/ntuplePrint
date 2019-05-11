#ifndef MergeHistsN_EXO16003_h
#define MergeHistsN_EXO16003_h


void MergeHistsN_EXO16003(float goalintlum, float xsec, int nfiles, std::string binname, std::string ohname, bool donorm, std::string bbname);
void  HistNormN_EXO16003(double &norm, float xsec, int nfiles, std::string binname, std::string bbname);
TH1F* HistManN_EXO16003(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname);
TH2F* HistMan2N_EXO16003(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname);


#endif
