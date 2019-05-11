#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"
#include "EMJbkgNew.h"
#include "BtagEff.h"


void QCDhistsNoMergeNew(int nrange[2], std::string samplename,std::string indir, bool hasPre, bool blind, bool b16003, std::string outdir, bool crabformat, bool isData, std::string runyr)
{
    std::string inputfile;
    std::string outputfile;

    QCDhists *qcdtools_=0;
    // first make histograms for each file in each bin for the qcd sample
    std::cout<<"making histograms for each file in each bin"<<std::endl;
    struct stat info;
    if( stat( (outdir+samplename).c_str(), &info ) != 0 ) {
        //printf( "%s does not exist.\n", (outdir+samplename).c_str() );
        printf( "Creating directory: [%s]\n", (outdir+samplename).c_str() );
        mkdir((outdir+samplename).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    for(int j=nrange[0]-1;j<nrange[1];j++) { //for each file for that bin
        bool printCuts=false;
        if (j==nrange[0]-1){
            qcdtools_->initialize();
            printCuts=true;
        }

        std::cout << "Processing file # " << j << std::endl;
        inputfile=indir+samplename+"/"+samplename+"_"+std::to_string(j+1)+"_0.histo.root";//condor format
        if (crabformat) inputfile=indir+samplename+"/ntuple_"+std::to_string(j+1)+".root";//crab format
        std::cout<<"input file is "<<inputfile<<std::endl;
        outputfile=outdir+samplename+"/histos"+samplename+"_"+std::to_string(j)+".root";
        std::cout<<"output file is "<<outputfile<<std::endl;
        int itmp;

        itmp = EMJbkgNew(true,hasPre,inputfile.c_str(),outputfile.c_str(),blind,isData,printCuts, runyr);
        std::cout<<"total number of events passing cuts is "<< itmp <<std::endl;
    }

    return;
}
