#ifndef EMJbkg_h
#define EMJbkg_h

int EMJbkg(bool otfile, bool hasPre, const char* inputfilename,const char* outputfilename,
           float HTcut, float pt1cut, float pt2cut, float pt3cut, float pt4cut, float jetacut,
           float alphaMaxcut, float maxIPcut, float NemfracCut,float CemfracCut,int ntrk1cut, 
           int NemergingCut,bool blind);

#endif
