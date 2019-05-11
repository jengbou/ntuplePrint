#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"
#include "EMJbkg.h"
#include "BtagEff.h"

//Remember to change runJobsUMD80New_allcuts.sh accordingly
//std::string cutstorun[] = {"1","2","3","4","5","6","7","8"};
//std::string cutstorun[] = {"1","2"};
std::string cutstorun[] = {"9","10"};

void QCDhistsNoMerge(int nrange[2], std::string samplename,std::string indir, bool hasPre, bool blind, bool b16003, std::string outdir, bool crabformat, bool isData)
{

    int numcuts = sizeof(cutstorun)/sizeof(cutstorun[0]);
    std::cout << "Number of cuts to run = " << numcuts << std::endl;
    for (int i=0;i<numcuts;i++) {
        std::cout << "cut: " << cutstorun[i] << std::endl;
    }

    QCDhists *qcdtools_=0;
    for (int icut=0;icut<numcuts;icut++) {
        std::string CutSet=cutstorun[icut];
        std::cout << "Running cut: " << CutSet << std::endl;
        std::string cutname="Cutset"+CutSet+"/";//Note: if this change, need to also change "CUTIDX" in submitJobs_*py

        // create directory under "ProdTag" level split by "samplename"
        struct stat info;
        if( stat( (outdir+cutname).c_str(), &info ) != 0 ) {
            //printf( "%s does not exist.\n", (outdir+cutname).c_str() );
            printf( "Creating directory: [%s]\n", (outdir+cutname).c_str() );
            mkdir((outdir+cutname).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
//         else if( info.st_mode & S_IFDIR )
//             printf( "%s is a directory\n", (outdir+cutname).c_str() );
//         else {
//             printf( "%s is no directory\n", (outdir+cutname).c_str() );
//         }

        if (icut==0){
            qcdtools_->initialize();
        }

        std::string inputfile;
        std::string outputfile;

        // opt CutSet: 1-4,9
        float DHTcut=900;
        float Dpt1cut=225;
        float Dpt2cut=100;
        float Dpt3cut=100;
        float Dpt4cut=100;
        float DMETcut=0;
        int DntagType=22;
        int Dnemcut=2;//This maybe superseded by ntagType in EMJbkg.cc
        float Djetacut = 2.0;
        int Dntrk1=0;


        //emcut1
        float DPUdzCut=2.5;
        float DsigzCut=4.0;
        float DmedIPcut=0.05;

        float DalphaCut=0.25;

        if (CutSet == "9" || CutSet == "10") DalphaCut = 0.4;

        if (CutSet == "2") {
            DPUdzCut = 4.0;
            DmedIPcut = 0.10;
        }
        else if (CutSet == "3") {
            // emcut3
            DPUdzCut = 15.0;
            DsigzCut = 30.0;
            DmedIPcut = 0.15;

            DntagType=12;
            Dnemcut=1;
            DMETcut=150.0;
        }
        else if (CutSet == "4") {
            // emcut4
            DPUdzCut = 4.0;
            DsigzCut = 20.0;
            DmedIPcut = 0.25;

            DntagType=12;
            Dnemcut=1;
            DMETcut=200.0;
        }
        else if (CutSet == "5" || CutSet == "10"){
            DHTcut=1100;
            Dpt1cut=275;
            Dpt2cut=250;
            Dpt3cut=150;
            Dpt4cut=150;

            //emcut1
        }
        else if (CutSet == "6"){
            DHTcut=1000;
            Dpt1cut=250;
            Dpt2cut=150;
            Dpt3cut=100;
            Dpt4cut=100;

            // emcut5
            DsigzCut = 4.0;
            DmedIPcut = 0.10;
        }
        else if (CutSet == "7"){
            DHTcut=1000;
            Dpt1cut=250;
            Dpt2cut=150;
            Dpt3cut=100;
            Dpt4cut=100;

            // emcut6
            DsigzCut = 20.0;
        }
        else if (CutSet == "8"){
            DHTcut=1200;
            Dpt1cut=300;
            Dpt2cut=250;
            Dpt3cut=200;
            Dpt4cut=150;

            // emcut7
            DsigzCut = 10.0;
        }
        else if (CutSet == "X") {
            DHTcut=900;
            Dpt1cut=100;
            Dpt2cut=100;
            Dpt3cut=100;
            Dpt4cut=100;

            // emcutX
            DPUdzCut = 2.5;
            DsigzCut = 8.0;
            DmedIPcut = 0.08;
            DalphaCut = 0.30;
        }

        std::cout << "Cutset[" << std::setw(2) << CutSet << "]: HT pt1 pt2 pt3 pt4 MET nemg(ntagType) PUdz 3Dsig medIP alpha3Dcut"
                  << std::endl
                  << std::setw(12) << " "
                  << DHTcut << " "
                  << Dpt1cut << " "
                  << Dpt2cut << " "
                  << Dpt3cut << " "
                  << Dpt4cut << " "
                  << DMETcut << " "
                  << Dnemcut << "(" << DntagType << ") "
                  << DPUdzCut << " "
                  << DsigzCut << " "
                  << DmedIPcut << " "
                  << DalphaCut << std::endl;

        // first make histograms for each file in each bin for the qcd sample
        std::cout<<"making histograms for each file in each bin"<<std::endl;

        // create directory under "CutSet" level split by "samplename"
        struct stat info1;
        if( stat( (outdir+cutname+samplename).c_str(), &info1 ) != 0 ) {
            //printf( "%s does not exist.\n", (outdir+cutname+samplename).c_str() );
            printf( "Creating directory: [%s]\n", (outdir+cutname+samplename).c_str() );
            mkdir((outdir+cutname+samplename).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        for(int j=nrange[0]-1;j<nrange[1];j++) { //for each file for that bin

            std::cout << "Processing file # " << j << std::endl;
            inputfile=indir+samplename+"/"+samplename+"_"+std::to_string(j+1)+"_0.histo.root";//condor format
            if (crabformat) inputfile=indir+samplename+"/ntuple_"+std::to_string(j+1)+".root";//crab format
            std::cout<<"input file is "<<inputfile<<std::endl;
            outputfile=outdir+cutname+samplename+"/histos"+samplename+"_"+std::to_string(j)+".root";
            std::cout<<"output file is "<<outputfile<<std::endl;
            int itmp;
            if(!b16003) {
                itmp = EMJbkg(true,hasPre,CutSet.c_str(),inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,DmedIPcut,0.9,0.9,Dntrk1,DntagType,DPUdzCut,DsigzCut,DMETcut,blind,isData);
                //itmp = BtagEff(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,DmedIPcut,0.9,0.9,Dntrk1,Dnemcut,blind,isData);
            } else {
                itmp = EMJ16003(true,hasPre,inputfile.c_str(),outputfile.c_str());
            }
            std::cout<<"total number of events passing cuts is "<< itmp <<std::endl;
        }
    }// End of icut loops
    return;
}
