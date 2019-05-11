#ifndef EMJbkg_h
#define EMJbkg_h

int EMJbkg(bool otfile, bool hasPre, std::string cutsetname, const char* inputfilename, const char* outputfilename,
           float HTcut, float pt1cut, float pt2cut, float pt3cut, float pt4cut, float jetacut,
           float alphaMaxcut, float maxIPcut, float NemfracCut, float CemfracCut, int ntrk1cut, 
           int ntagType, float PUdzCut, float sigzCut, float METcut, bool blind, bool isData);

float DeltaR(float eta1, float phi1, float eta2, float phi2);

#endif
