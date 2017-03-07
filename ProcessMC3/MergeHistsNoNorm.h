#ifndef MergeHistsNoNorm_h
#define MergeHistsNoNorm_h

void MergeHistsNoNorm(int fidx, int nrange[2], std::string binname,std::string aaname, std::string bbname);
TH1F* HistManNoNorm(std::string thisHIST, int nrange[2], std::string binname, std::string bbname);
TH2F* HistMan2NoNorm(std::string thisHIST, int nrange[2], std::string binname, std::string bbname);


#endif
