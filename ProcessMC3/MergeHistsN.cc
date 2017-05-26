#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include "QCDhists.h"
#include "MergeHistsN.h"

void MergeHistsN(float goalintlum, float xsec, int nfiles, std::string binname, std::string ohname, bool donorm, std::string bbname)
{

    std::string inputfile;
    std::string outputfile;

    // first make histograms for each file in each bin for the qcd sample
    // get normalization
    double norm;
    if(donorm) {
        HistNormN(norm,xsec,nfiles,binname,bbname);  // this gives the total number of events in each bin before all selections using the eventCountPreTrigger histogram
    } else{
        norm=1.;
    }
    std::cout<<"total number events is "<<norm<<std::endl;

    std::cout<<"normalizing histograms"<<std::endl;
    // merge 1D histograms
    const int nhist=100;
    std::vector<TH1F*> vv(nhist);
    std::string histnames[nhist]={
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
    double outnorm;
    for(int i=0;i<nhist;i++) {
        std::cout<<" entering Histman with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=HistManN(goalintlum,histnames[i],norm,outnorm,xsec,nfiles,binname,donorm,bbname);
    }

    // merge 2D histograms
    const int nhist2=12;
    std::vector<TH2F*> vv2(nhist2);
    std::string histnames2[nhist2]={
        "aMip","haMvjpt","haMvHT","haMvnvtx",
        "halphavtheta2D","halphavipXYSig","htheta2DvipXYSig",
        "halphavtheta2DPS","halphavipXYSigPS","htheta2DvipXYSigPS",
        "htheta2DvipXYSigSR",
        "htheta2DvipXYSigFR",
    };
    double outnorm2;
    for(int i=0;i<nhist2;i++) {
        std::cout<<" entering Histman2 with i = "<<i<<": "<<histnames2[i]<<std::endl;
        vv2[i]=HistMan2N(goalintlum,histnames2[i],norm,outnorm2,xsec,nfiles,binname,donorm,bbname);
    }

    // output total event count
    std::cout<<" initial event count before and after norm is"<<std::endl;
    double ttotal=0;
    std::cout<<" norm "<<norm<<" times outnorm is "<<norm*outnorm<<std::endl;
    ttotal += norm*outnorm;

    std::cout<<"total is "<<ttotal<<std::endl;;

    std::cout<<"outputting histograms"<<std::endl;
    outputfile=bbname+ohname;
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


void  HistNormN(double &norm, float xsec, int nfiles, std::string binname, std::string bbname) {

    std::cout<<"entering HistNorm"<<std::endl; 

    std::string inputfile;
    TFile * in;

    // now add up all the files for one bin
    vector<TH1F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=bbname+binname+"/SumHistsNoNorm"+binname+"_"+std::to_string(j)+".root";
        std::cout<<j<<" "<<inputfile<<std::endl;
        in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if(j==0) {
            sum[0] = *(static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone()));
        } else {
            TH1F* tmp = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
            sum[0].Add(tmp);
        }
        in->Close();
    }

    // reweight to int lum

    // get total number of events before filter
    norm = sum[0].GetBinContent(2);
    std::cout<<"norm "<<" "<<norm<<std::endl;


    return;
}



TH1F* HistManN(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname) {

    std::string inputfile;


    // now add up all the files for one bin
    vector<TH1F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=bbname+binname+"/SumHistsNoNorm"+binname+"_"+std::to_string(j)+".root";
        std::cout << inputfile << std::endl;
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if (!in->GetListOfKeys()->Contains(thisHIST.c_str())) return (new TH1F(thisHIST.c_str(),"dummy empty hist",10,0.,10.));

        if(j==0) {
            std::cout<<" adding up histos within a bin"<<std::endl;
            sum[0] = *(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone()));
        } else {
            TH1F* tmp = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
            sum[0].Add(tmp);
        }
        in->Close();
    }

    if(donorm) {
        // reweight to int lum
        std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
        // get total number of events before filter
        float ntotal = norm;
        std::cout<<"number of pretrigger events is "<<ntotal<<std::endl;
        float fileLum= ntotal/xsec;
        std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
        outnorm = goalintlum/fileLum;
        std::cout<<" scaling by a factor of "<<outnorm<<std::endl;
        sum[0].Scale(outnorm);
    }

    TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone());
    return SUM;
}

TH2F* HistMan2N(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string binname,bool donorm,std::string bbname) {

    std::string inputfile;


    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=bbname+binname+"/SumHistsNoNorm"+binname+"_"+std::to_string(j)+".root";
        TFile* in = new TFile(inputfile.c_str());
        if (in->IsZombie()) {in->Close(); continue;}
        if(j==0) {
            sum[0] = *(static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone()));
        } else {
            TH2F* tmp = static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone());
            sum[0].Add(tmp);
        }
        in->Close();
    }

    if(donorm) {
        // reweight to int lum
        std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
        // get total number of events before filter
        float ntotal = norm;
        std::cout<<"number of pretrigger events is "<<ntotal<<std::endl;
        float fileLum= ntotal/xsec;
        std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
        outnorm = goalintlum/fileLum;
        std::cout<<" scaling by a factor of "<<outnorm<<std::endl;
        sum[0].Scale(outnorm);
    }

    TH2F* SUM=static_cast<TH2F*>((sum[0]).Clone());
    return SUM;

}
