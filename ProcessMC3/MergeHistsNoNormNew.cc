#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include "QCDhists.h"
#include "MergeHistsNoNormNew.h"

void MergeHistsNoNormNew(int fidx, int nrange[2], std::string samplename, std::string dirname)
{

    std::string Cutstorun[] = {"1","21"};
    //std::string Cutstorun[] = {"1","2","3","4","5","6","7","8"};
    //std::string Cutstorun[] = {"9","10","11"};
    //std::string Cutstorun[] = {"11","11a"};
    //std::string Cutstorun[] = {"1","2","5","6","7","8","9"};
    const int nCuts = sizeof(Cutstorun)/sizeof(Cutstorun[0]);

    std::string inputfile;
    std::string outputfile;

    // first make histograms for each file in each bin for the qcd sample
    // get normalization


    std::cout<<"normalizing histograms"<<std::endl;
    // merge 1D histograms
    const int nhist=57;
    std::string histnames[nhist]={
        "eventCountPreTrigger",
        "count","acount",
        "h_nemg","h_nalemg",
        "H_TSR","H_TFR",
        "halpha","halphaZero",
        "hnjet",
        "hntrkSR","hntrkFR",
        "hmassSR","hmassFR",
        "hjptaSR","hjptaFR",
        "hetaaSR","hetaaFR",
        "h_nemgSR","h_nemgFR",
        "hnjetSR","hnjetFR",
        "hmedipXYSigSR","hmedtheta2DSR",
        "hlogmedipXYSigSR","hlogmedtheta2DSR",
        "hmedipXYSigFR","hmedtheta2DFR",
        "hlogmedipXYSigFR","hlogmedtheta2DFR",
        "hjpt1SR","hjpt1FR",
        "hjpt2SR","hjpt2FR",
        "hjpt3SR","hjpt3FR",
        "hjpt4SR","hjpt4FR",
        "heta1SR","heta1FR",
        "heta2SR","heta2FR",
        "heta3SR","heta3FR",
        "heta4SR","heta4FR",
        "METSR","METFR",
        "dRjtrk0SR","dRjtrk0FR",
        "dRjtrkSR","dRjtrkFR",
        "djtrkSR","djtrkFR",
        "chiAll","dzAll",
        "medipXY",
    };

    // merge 2D histograms
    const int nhist2=2;
    std::string histnames2[nhist2]={
        "htheta2DvipXYSigSR",
        "htheta2DvipXYSigFR",
    };

    for (int iCut=0;iCut<nCuts;iCut++) {
        std::string histcutname="Cutset"+Cutstorun[iCut];
        std::vector<TH1F*> vv(nhist);
        for(int i=0;i<nhist;i++) {
            std::cout<<" entering HistmanNoNormNew with i = "<<i<<": "<<histnames[i]<<std::endl;
            vv[i]=HistManNoNormNew(histnames[i],histcutname,nrange,samplename,dirname);
        }

        std::vector<TH2F*> vv2(nhist2);
        for(int i=0;i<nhist2;i++) {
            std::cout<<" entering Histman2NoNormNew with i = "<<i<<": "<<histnames2[i]<<std::endl;
            vv2[i]=HistMan2NoNormNew(histnames2[i],histcutname,nrange,samplename,dirname);
        }

        std::cout<<"outputting histograms:"<<std::endl;

        struct stat info;
        if( stat( (dirname+"/"+histcutname).c_str(), &info ) != 0 ) {
            //printf( "%s does not exist.\n", (dirname+"/"+histcutname).c_str() );
            printf( "Creating directory: [%s]\n", (dirname+"/"+histcutname).c_str() );
            mkdir((dirname+"/"+histcutname).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
        struct stat info1;
        if( stat( (dirname+"/"+histcutname+"/"+samplename).c_str(), &info1 ) != 0 ) {
            //printf( "%s does not exist.\n", (dirname+"/"+histcutname+"/"+samplename).c_str() );
            printf( "Creating directory: [%s]\n", (dirname+"/"+histcutname+"/"+samplename).c_str() );
            mkdir((dirname+"/"+histcutname+"/"+samplename).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        outputfile=dirname+"/"+histcutname+"/"+samplename+"/SumHistsNoNorm"+samplename+"_"+std::to_string(fidx)+".root";
        std::cout<<outputfile<<std::endl;
        TFile out(outputfile.c_str(),"RECREATE");

        for(int i=0;i<nhist;i++) {
            vv[i]->Write();
        }
        for(int i=0;i<nhist2;i++) {
            vv2[i]->Write();
        }
    }

    return;
}


TH1F* HistManNoNormNew(std::string thisHIST, std::string thisCUT, int nrange[2], std::string samplename, std::string dirname) {

    std::string inputfile;

    // now add up all the files for one bin
    vector<TH1F> sum(1);
    int idx=0;
    for(int j=nrange[0]-1;j<nrange[1];j++) {
        inputfile=dirname+samplename+"/histos"+samplename+"_"+std::to_string(j)+".root";
        std::cout << inputfile << std::endl;
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if (!in->GetListOfKeys()->Contains((thisHIST+"_"+thisCUT).c_str())) return (new TH1F(thisHIST.c_str(),"dummy empty hist",10,0.,10.));

        if(idx==0) {
            std::cout<<" adding up histos within a bin"<<std::endl;
            sum[0] = *(static_cast<TH1F*>(in->Get((thisHIST+"_"+thisCUT).c_str())->Clone(thisHIST.c_str())));
        } else {
            TH1F* tmp = static_cast<TH1F*>(in->Get((thisHIST+"_"+thisCUT).c_str())->Clone(thisHIST.c_str()));
            sum[0].Add(tmp);
        }
        in->Close();
        idx++;
    }
    TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone(thisHIST.c_str()));
    return SUM;
}

TH2F* HistMan2NoNormNew(std::string thisHIST, std::string thisCUT, int nrange[2], std::string samplename, std::string dirname) {

    std::string inputfile;

    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(1);
    int idx=0;
    for(int j=nrange[0]-1;j<nrange[1];j++) {
        inputfile=dirname+samplename+"/histos"+samplename+"_"+std::to_string(j)+".root";
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if(idx==0) {
            sum[0] = *(static_cast<TH2F*>(in->Get((thisHIST+"_"+thisCUT).c_str())->Clone(thisHIST.c_str())));
        } else {
            TH2F* tmp = static_cast<TH2F*>(in->Get((thisHIST+"_"+thisCUT).c_str())->Clone(thisHIST.c_str()));
            sum[0].Add(tmp);
        }
        in->Close();
        idx++;
    }

    TH2F* SUM=static_cast<TH2F*>((sum[0]).Clone(thisHIST.c_str()));
    return SUM;
}
