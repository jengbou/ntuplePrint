#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include "QCDhists.h"
#include "MergeHists.h"


void MergeHists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* samplenames,std::string indir,std::string ohname,bool donorm, std::string outdir)
{
    QCDhists* qcdtools=0;
    std::string inputfile;
    std::string outputfile;

    // first make histograms for each file in each bin for the qcd sample
    // get normalization
    vector<double> norm(nbin);
    if(donorm) {
        qcdtools->HistNorm(norm,nbin,xsec,nfiles,samplenames,outdir);  // this gives the total number of events in each bin before all selections using the eventCountPreTrigger histogram
    } else{
        for(int i=0;i<nbin;i++) norm[i]=1.;
    }
    for(int i=0;i<nbin;i++) {
        std::cout<<"total number events in bin "<<i<<" is "<<norm[i]<<std::endl;
    }
    TH1F* countclone = new TH1F("countclone","unnormalized count",20,0,20);
    for(int i=0;i<nbin;i++){
        countclone->AddBinContent(i+1,norm[i]);
    }
  
    TH1F* normhst = new TH1F("normhst","counts pretrigger by bin",nbin,0.,nbin);
    for(int i=0;i<nbin;i++){
        normhst->AddBinContent(i+1,norm[i]);
    }


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
        //"hHTnm1","hpt1nm1","hpt2nm1","hpt3nm1","hpt4nm1","halphanm1",
        //"hmaxipnm1","hnHitsnm1","hntrk1nm1","hnemnm1",
        "hipXYEJ","hipXYnEJ","htvwEJ","htvw","hipXYSigEJ","hipXYSignEJ",
        "hmaxipXYEJ","hmaxipXYnEJ","hmeanipXYEJ","hmeanipXYnEJ","hnmaxipnm1",
        //"hn2maxipnm1","hjptfrb","hjptfra1",
        //"hjptfra2","hjptfrbc","hjptfra1c","hjptfra2c","hjptb",
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

    // runtype==1
//     const int nhist=2;
//     std::vector<TH1F*> vv(nhist);
//     std::string histnames[nhist]={"count","acount"};

    vector<double> outnorm(nbin);
    for(int i=0;i<nhist;i++) {
        std::cout<<" entering Histman with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=qcdtools->HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,samplenames,donorm,outdir);
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

    vector<double> outnorm2(nbin);
    for(int i=0;i<nhist2;i++) {
        std::cout<<" entering Histman2 with i = "<<i<<": "<<histnames2[i]<<std::endl;
        vv2[i]=qcdtools->HistMan2(goalintlum,histnames2[i],norm,outnorm2,nbin,xsec,nfiles,samplenames,donorm,outdir);
    }

    // output total event count
    std::cout<<" initial event count before and after norm is"<<std::endl;
    double ttotal=0;
    for(int i=0;i<nbin;i++) {
        std::cout<<" bin "<<i<<" norm "<<norm[i]<<" times outnorm is "<<norm[i]*outnorm[i]<<std::endl;
        ttotal = ttotal + norm[i]*outnorm[i];
    }
    std::cout<<"total is "<<ttotal<<std::endl;;

    std::cout<<"outputting histograms"<<std::endl;
    outputfile=outdir+ohname;
    TFile out(outputfile.c_str(),"RECREATE");
    normhst->Write();

    countclone->Write();
    for(int i=0;i<nhist;i++) {
        vv[i]->Write();
    }
    for(int i=0;i<nhist2;i++) {
        vv2[i]->Write();
    }

    return;
}
