#ifndef MergeHistsNoNormNew_h
#define MergeHistsNoNormNew_h

void MergeHistsNoNormNew(int fidx, int nrange[2], std::string samplename, std::string dirname);
TH1F* HistManNoNormNew(std::string thisHIST, std::string thisCUT, int nrange[2], std::string samplename, std::string dirname);
TH2F* HistMan2NoNormNew(std::string thisHIST, std::string thisCUT,  int nrange[2], std::string samplename, std::string dirname);


#endif
