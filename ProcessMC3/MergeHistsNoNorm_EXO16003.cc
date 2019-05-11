#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include "QCDhists.h"
#include "MergeHistsNoNorm_EXO16003.h"

void MergeHistsNoNorm_EXO16003(int fidx, int nrange[2], std::string binname, std::string bbname)
{

    std::string inputfile;
    std::string outputfile;

    // first make histograms for each file in each bin for the qcd sample
    // get normalization


    std::cout<<"normalizing histograms"<<std::endl;
    // merge 1D histograms
    const int nhist=101;
    std::vector<TH1F*> vv(nhist);
    std::string histnames[nhist]={
        "eventCountPreTrigger",
        "count","acount","hfrwgt",
        "hjetcut","hjetchf","h_nemg","h_nalemg",
        "hnjet","hpt","heta","heta2",
        "H_T","H_T0","H_T1","H_T2","H_T3","H_T4","H_TFR",
        "hpt1","hpt2","hpt3",
        "hpt4","hbcut_ntrkpt1","hacut_ntrkpt1","hbcut_nef","hacut_nef",
        "hbcut_cef","hacut_cef","hbcut_alphamax","hacut_alphamax",
//         "hHTnm1","hpt1nm1","hpt2nm1","hpt3nm1","hpt4nm1","halphanm1",
//         "hmaxipnm1","hnHitsnm1","hntrk1nm1","hnemnm1",
        "hipXYEJ","hipXYnEJ","htvwEJ","htvw","hipXYSigEJ","hipXYSignEJ",
        "hmaxipXYEJ","hmaxipXYnEJ","hmeanipXYEJ","hmeanipXYnEJ","hnmaxipnm1",
//         "hn2maxipnm1","hjptfrb","hjptfra1",
//         "hjptfra2","hjptfrbc","hjptfra1c","hjptfra2c","hjptb",
        "hjpta","haMgj","hHTko","hpt1ko","hpt2ko",
        "hpt3ko","hpt4ko","hmass","hmassFR","hlogmedipXYSigEJ","hlogmedipXYSignEJ","hlogmeanipXYSigEJ","hlogmeanipXYSignEJ",
        "hmedipXYSigEJ","hmedipXYSignEJ","hmeanipXYSigEJ","hmeanipXYSignEJ","hmedipXYEJ","hmedipXYnEJ",
        "hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n",
        "h_ntag",
        "h_nloosetag",
        "halpha","halphaPS",
        "halphaZero","halphaZeroPS",
        "hmedtheta2DEJ","hmedtheta2DnEJ","hlogmedtheta2DEJ","hlogmedtheta2DnEJ",
        "hmedtheta2DPS","hlogmedtheta2DPS","hmedipXYSigPS","hlogmedipXYSigPS",
        "hmedtheta2DSR","hlogmedtheta2DSR","hmedipXYSigSR","hlogmedipXYSigSR",
        "hmedtheta2DFR","hlogmedtheta2DFR","hmedipXYSigFR","hlogmedipXYSigFR",
        "hfr_ntrkpt1d","hfr_ntrkpt1n","hntrkSR","hntrkFR",
        "hjptaSR","hjptaFR","hetaaSR","hetaaFR","h_nemgSR","h_nemgFR",
        "hnjetSR","hnjetFR",
    };

    for(int i=0;i<nhist;i++) {
        std::cout<<" entering HistmanNoNorm with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=HistManNoNorm_EXO16003(histnames[i],nrange,binname,bbname);
    }

    // merge 2D histograms
    const int nhist2=11;
    std::vector<TH2F*> vv2(nhist2);
    std::string histnames2[nhist2]={
        "aMip","haMvjpt","haMvHT","haMvnvtx",
        "halphavtheta2D","halphavipXYSig","htheta2DvipXYSig",
        "halphavtheta2DPS","halphavipXYSigPS","htheta2DvipXYSigPS",
        "htheta2DvipXYSigSR",
        //"htheta2DvipXYSigFR",
    };

    for(int i=0;i<nhist2;i++) {
        std::cout<<" entering Histman2NoNorm with i = "<<i<<": "<<histnames2[i]<<std::endl;
        vv2[i]=HistMan2NoNorm_EXO16003(histnames2[i],nrange,binname,bbname);
    }

    std::cout<<"outputting histograms:"<<std::endl;
    outputfile=bbname+binname+"/SumHistsNoNorm"+binname+"_"+std::to_string(fidx)+".root";
    std::cout<<outputfile<<std::endl;
    TFile out(outputfile.c_str(),"RECREATE");

    for(int i=0;i<nhist;i++) {
        vv[i]->Write();
    }
    for(int i=0;i<nhist2;i++) {
        vv2[i]->Write();
    }

    return;
}


TH1F* HistManNoNorm_EXO16003(std::string thisHIST, int nrange[2], std::string binname, std::string bbname) {

    std::string inputfile;

    // now add up all the files for one bin
    vector<TH1F> sum(1);
    int idx=0;
    for(int j=nrange[0]-1;j<nrange[1];j++) {
        inputfile=bbname+binname+"/histos"+binname+"_"+std::to_string(j)+".root";
        std::cout << inputfile << std::endl;
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if (!in->GetListOfKeys()->Contains(thisHIST.c_str())) return (new TH1F(thisHIST.c_str(),"dummy empty hist",10,0.,10.));

        if(idx==0) {
            std::cout<<" adding up histos within a bin"<<std::endl;
            sum[0] = *(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone()));
        } else {
            TH1F* tmp = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
            sum[0].Add(tmp);
        }
        in->Close();
        idx++;
    }
    TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone());
    return SUM;
}

TH2F* HistMan2NoNorm_EXO16003(std::string thisHIST, int nrange[2], std::string binname, std::string bbname) {

    std::string inputfile;

    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(1);
    int idx=0;
    for(int j=nrange[0]-1;j<nrange[1];j++) {
        inputfile=bbname+binname+"/histos"+binname+"_"+std::to_string(j)+".root";
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if(idx==0) {
            sum[0] = *(static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone()));
        } else {
            TH2F* tmp = static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone());
            sum[0].Add(tmp);
        }
        in->Close();
        idx++;
    }

    TH2F* SUM=static_cast<TH2F*>((sum[0]).Clone());
    return SUM;
}
