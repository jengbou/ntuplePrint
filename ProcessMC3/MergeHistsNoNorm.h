#ifndef MergeHistsNoNorm_h
#define MergeHistsNoNorm_h

void MergeHistsNoNorm(int fidx, int nrange[2], std::string samplename, std::string dirname);
TH1F* HistManNoNorm(std::string thisHIST, int nrange[2], std::string samplename, std::string dirname);
TH2F* HistMan2NoNorm(std::string thisHIST, int nrange[2], std::string samplename, std::string dirname);


#endif
