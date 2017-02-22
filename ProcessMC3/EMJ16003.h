#ifndef EMJ16003_h
#define EMJ16003_h

#include "vector"
using std::vector;
#include "algorithm"


float CalcMedian(std::vector<float> vec);
int EMJ16003(bool otfile, bool hasPre, const char* inputfilename,const char* outputfilename);

#endif
