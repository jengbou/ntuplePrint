#ifndef EMJscan_h
#define EMJscan_h

#include <iostream>
#include <iomanip>
#include <locale>

#include "vector"
using std::vector;
#include "algorithm"



vector<float> Decode(int cutindex, int ncut,vector<int> nstep, vector<float> stepsize);



vector<int> EMJscan(const char* inputfilename,
                    float HTcutmin,int NHTcut, float HTcutSS,
                    float pt1cutmin, int Npt1cut, float pt1cutSS,
                    float pt2cutmin,  int Npt2cut,float pt2cutSS,
                    float pt3cutmin,  int Npt3cut,float pt3cutSS,
                    float pt4cutmin, int Npt4cut,float pt4cutSS,
		    int NemergingCutmin, int NNemergingCut, int NNemergingCutSS,
		    float jetacut,
		    float alphaMaxcut, float maxIPcut,
		    float NemfracCut,float CemfracCut,int ntrk1cut,bool blind);

#endif
