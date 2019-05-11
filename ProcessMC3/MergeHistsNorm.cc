#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include "QCDhists.h"
#include "MergeHistsNorm.h"

void MergeHistsNorm(float goalintlum, float xsec, int nfiles, std::string samplename, std::string ohname, bool donorm, std::string dirname)
{

    std::string inputfile;
    std::string outputfile;

    // first make histograms for each file in each bin for the qcd sample
    // get normalization
    double norm;
    if(donorm) {
        HistNormNew(norm,xsec,nfiles,samplename,dirname);  // this gives the total number of events in each bin before all selections using the eventCountPreTrigger histogram
    } else{
        norm=1.;
    }
    std::cout<<"total number events is "<<norm<<std::endl;

    std::cout<<"normalizing histograms"<<std::endl;
    // merge 1D histograms
    const int nhist=56;
    std::vector<TH1F*> vv(nhist);
    std::string histnames[nhist]={
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

    // runtype==1
//     const int nhist=2;
//     std::vector<TH1F*> vv(nhist);
//     std::string histnames[nhist]={"count","acount"};

    double outnorm;
    for(int i=0;i<nhist;i++) {
        std::cout<<" entering Histman with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=HistManNorm(goalintlum,histnames[i],norm,outnorm,xsec,nfiles,samplename,donorm,dirname);
    }

    // merge 2D histograms
    const int nhist2=2;
    std::vector<TH2F*> vv2(nhist2);
    std::string histnames2[nhist2]={
        "htheta2DvipXYSigSR",
        "htheta2DvipXYSigFR",
    };

    //runtype==1
//     const int nhist2=8;
//     std::vector<TH2F*> vv2(nhist2);
//     std::string histnames2[nhist2]={
//         "btagEff_Den_b",
//         "btagEff_Den_udsg",
//         "btagEff_Num_b",
//         "btagEff_Num_udsg",
//         "btagEff_Den_b_all",
//         "btagEff_Den_udsg_all",
//         "btagEff_Num_b_all",
//         "btagEff_Num_udsg_all",
//     };

    double outnorm2;
    for(int i=0;i<nhist2;i++) {
        std::cout<<" entering Histman2 with i = "<<i<<": "<<histnames2[i]<<std::endl;
        vv2[i]=HistMan2Norm(goalintlum,histnames2[i],norm,outnorm2,xsec,nfiles,samplename,donorm,dirname);
    }

    // output total event count
    std::cout<<" initial event count before and after norm is"<<std::endl;
    double ttotal=0;
    std::cout<<" norm "<<norm<<" times outnorm is "<<norm*outnorm<<std::endl;
    ttotal += norm*outnorm;

    std::cout<<"total is "<<ttotal<<std::endl;;

    std::cout<<"outputting histograms"<<std::endl;
    outputfile=dirname+ohname;
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


void  HistNormNew(double &norm, float xsec, int nfiles, std::string samplename, std::string dirname) {

    std::cout<<"entering HistNorm"<<std::endl; 

    std::string inputfile;
    TFile * in;

    // now add up all the files for one bin
    vector<TH1F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=dirname+samplename+"/SumHistsNoNorm"+samplename+"_"+std::to_string(j)+".root";
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



TH1F* HistManNorm(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname) {

    std::string inputfile;


    // now add up all the files for one bin
    vector<TH1F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=dirname+samplename+"/SumHistsNoNorm"+samplename+"_"+std::to_string(j)+".root";
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

TH2F* HistMan2Norm(float goalintlum,std::string thisHIST,double& norm,double& outnorm,float xsec, int nfiles, std::string samplename,bool donorm,std::string dirname) {

    std::string inputfile;


    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(1);
    for(int j=0;j<nfiles;j++) { //for each file for that bin
        inputfile=dirname+samplename+"/SumHistsNoNorm"+samplename+"_"+std::to_string(j)+".root";
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
