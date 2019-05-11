#include <TList.h>
#include "QCDhists.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"
#include "EMJbkg.h"
#include "BtagEff.h"
#include "TH1.h"
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "BTagCalibrationStandalone.cpp"
#include "TRandom.h"

//const int runtype=1;//0: EMJbkg; 1: BtagEff

const bool mergeOnly = false;
bool scanCuts = false;
bool verbose = false;

//void QCDhists()
QCDhists::QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* samplenames,std::string indir,std::string ohname, int dooptk, int doopta, bool hasPre,bool donorm, bool blind, bool b16003, std::string outdir="./", bool crabformat=true,bool isData=false)
{
    std::string CutSet="1";

    if (b16003) scanCuts = false;
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

    if (CutSet == "9") DalphaCut = 0.4;
    if (CutSet == "10" || CutSet == "11" || CutSet == "11a") DalphaCut = 0.5;

    if (CutSet == "2") {
        DPUdzCut = 4.0;
        DmedIPcut = 0.10;
    }
    else if (CutSet == "3" || CutSet == "10") {
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
    else if (CutSet == "5"){
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
    else if (CutSet == "11") {
        // emcut11
        DPUdzCut = 4.0;
        DsigzCut = 20.0;
        DmedIPcut = 0.10;

        DntagType=12;
        Dnemcut=1;
        DMETcut=200.0;
    }
    else if (CutSet == "11a") {
        // emcut11a
        DPUdzCut = 4.0;
        DsigzCut = 20.0;
        DmedIPcut = 0.10;

        DntagType=1;
        Dnemcut=1;
        DMETcut=200.0;
    }

    const int ncutscan=3;
    //const int ncutscan=1;

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

    // scan cut related
    const int nkincut=6;
    vector<int> nstep {4,4,4,4,4,3};
    //vector<int> nstep {2,2,2,2,2,2};
    // ht pt1 pt2 pt3 pt4 nemerging
    vector<float> cutmin {1000.,400.,200.,200.,100.,0};
    vector<float> ss {50,10,10,10,10,1};
    int iicut =nstep[0];
    for (int hh=1;hh<nkincut;hh++) iicut*=nstep[hh];
    vector < vector <int> > nnpass(iicut,vector<int>(nbin,0));  
    vector < vector <int> > ipass(ncutscan, vector<int>(nbin,0));
    TH1F* cutscan = new TH1F("cutscan","n pass versus cut",ncutscan,0.,ncutscan);
    TH1F* kcutscan = new TH1F("kcutscan","n pass versus cut kin cuts",iicut,0.,iicut);

    // first make histograms for each file in each bin for the qcd sample
    if (!mergeOnly){
        std::cout<<"making histograms for each file in each bin"<<std::endl;
        for(int i=0;i<nbin;i++) {  // for each bin
            mkdir((outdir+samplenames[i]).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            for(int j=0;j<nfiles[i];j++) { //for each file for that bin
                if (j==0){
                    initialize();
                }
                //inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                //inputfile=indir+samplenames[i]+"/ntuple_"+samplenames[i]+"_"+std::to_string(j+1)+"_hlt1p1.root";

                inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.histo.root";//condor format
                if (crabformat) inputfile=indir+samplenames[i]+"/ntuple_"+std::to_string(j+1)+".root";//crab format
                std::cout<<"input file is "<<inputfile<<std::endl;
                outputfile=outdir+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
                std::cout<<"output file is "<<outputfile<<std::endl;
                int itmp;
                if(!b16003) {
                    //itmp = EMJselect(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,DmedIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
                    itmp = EMJbkg(true,hasPre,CutSet.c_str(),inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,DmedIPcut,0.9,0.9,Dntrk1,DntagType,DPUdzCut,DsigzCut,DMETcut,blind,isData);
                    //itmp = BtagEff(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,DmedIPcut,0.9,0.9,Dntrk1,Dnemcut,blind,isData);
                } else {
                    itmp = EMJ16003(true,hasPre,inputfile.c_str(),outputfile.c_str());
                }
            }
        }
    }
    if (scanCuts){
        // do some cut optimization on cuts not related to choosing the emerging jets

        if(dooptk==1) {
            for(int i=0;i<nbin;i++) {  // for each bin
                for(int j=0;j<nfiles[i];j++) { //for each file for that bin
                    //inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                    inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.histo.root";// condor format
                    if (crabformat) inputfile=indir+samplenames[i]+"/ntuple_"+std::to_string(j+1)+".root";// crab format
                    std::cout<<"input file is "<<inputfile<<std::endl;


                    vector<int> npass = EMJscan(inputfile.c_str(),
                                                cutmin[0],nstep[0],ss[0],
                                                cutmin[1],nstep[1],ss[1],
                                                cutmin[2],nstep[2],ss[2],
                                                cutmin[3],nstep[3],ss[3], 
                                                cutmin[4],nstep[4],ss[4],
                                                cutmin[5],nstep[5],ss[5],
                                                Djetacut,
                                                0.04,-1.,0.9,0.9,0,blind);
                    for(int tt=0;tt<iicut;tt++) {
                        nnpass[tt][i]=nnpass[tt][i]+npass[tt];
                    }
                }
            }
        }


        //vector<float> decode(6);
        //decode = Decode(icut,6,nstep,stepsize);

        //do some cut optimization on alpha max

        float acut = 1.0;
        if(doopta==2) acut=1.2;
        //  int ipass[ncutscan][nbin];

        if(doopta>0) {
            for(int k=0;k<ncutscan;k++) {
                float acut2=(acut/(ncutscan))*(k+1);
                std::cout<<" cut value is "<<acut2<<std::endl;
                for(int i=0;i<nbin;i++) ipass[k][i]=0;
                for(int i=0;i<nbin;i++) {  // for each bin
                    for(int j=0;j<nfiles[i];j++) { //for each file for that bin
                        std::cout<<"k i j="<<k<<" "<<i<<" "<<j<<std::endl;
                        //inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                        //inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+".root";
                        inputfile=indir+samplenames[i]+"/"+samplenames[i]+"_"+std::to_string(j+1)+"_0.histo.root";// condor format
                        if (crabformat) inputfile=indir+samplenames[i]+"/ntuple_"+std::to_string(j+1)+".root";// crab format
                        std::cout<<"input file is "<<inputfile<<std::endl;
                        int iii=0;
                        outputfile=outdir+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
                        std::cout<<"output file is "<<outputfile<<std::endl;

                        if(doopta==1) {
                            iii = EMJselect(false,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,acut2,DmedIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
                        } else {
                            iii = EMJselect(false,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,DalphaCut,acut2,0.9,0.9,Dntrk1,Dnemcut,blind);
                        }
                        ipass[k][i]+=iii;
                        std::cout<<" iii ipass  is "<<iii<<" "<<ipass[k][i]<<std::endl;
                    }
                }
            }
        }
    }
    // get normalization
    //  double norm[nbin];
    vector<double> norm(nbin);
    if(donorm) {
        HistNorm(norm,nbin,xsec,nfiles,samplenames,outdir);  // this gives the total number of events in each bin before all selections using the eventCountPreTrigger histogram
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


    //make and  output summed and renormalized histograms
    std::cout<<"normalizing histograms"<<std::endl;

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
        "hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n","h_ntag","h_nloosetag",
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
        vv[i]=HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,samplenames,donorm,outdir);
    }


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
        vv2[i]=HistMan2(goalintlum,histnames2[i],norm,outnorm2,nbin,xsec,nfiles,samplenames,donorm,outdir);
    }

    // output total event count
    std::cout<<" initial event count before and after norm is"<<std::endl;
    double ttotal=0;
    for(int i=0;i<nbin;i++) {
        std::cout<<" bin "<<i<<" norm "<<norm[i]<<" times outnorm is "<<norm[i]*outnorm[i]<<std::endl;
        ttotal = ttotal + norm[i]*outnorm[i];
    }
    std::cout<<"total is "<<ttotal<<std::endl;;

    if (scanCuts){
        // normalize cut scan and sum bins
        std::cout<<"normalizing kinematic cut stuff"<<std::endl;
        std::cout<<"iicut "<<iicut<<std::endl;
        double ffpass[iicut];
        for(int i=0;i<iicut;i++) ffpass[i]=0;
        for(int k=0;k<iicut;k++) {
            for(int i=0;i<nbin;i++) {
                ffpass[k]+=nnpass[k][i]*outnorm[i];
            }
            std::cout<<" output kinematic cut scan "<<k<<" "<<ffpass[k]<<std::endl;
        }

        for(int i=0;i<iicut;i++){
            kcutscan->AddBinContent(i+1,ffpass[i]);
        }


        // normalize cut scan and sum bins
        std::cout<<"normalizing alpha scan stuff"<<std::endl;
        double fpass[ncutscan];
        for(int i=0;i<ncutscan;i++) fpass[i]=0;
        for(int k=0;k<ncutscan;k++) {
            for(int i=0;i<nbin;i++) {
                //      std::cout<<"k i "<<k<<" "<<i<<" "<<ipass[k][i]<<" "<<outnorm[i]<<std::endl;
                fpass[k]+=ipass[k][i]*outnorm[i];
            }
            std::cout<<" output alphamax scan "<<k<<" "<<fpass[k]<<std::endl;
        }

        for(int i=0;i<ncutscan;i++){
            cutscan->AddBinContent(i+1,fpass[i]);
        }
    }


    std::cout<<"outputting histograms"<<std::endl;
    outputfile=outdir+ohname;
    TFile out(outputfile.c_str(),"RECREATE");
    normhst->Write();
    if (scanCuts){
        cutscan->Write();
        kcutscan->Write();
    }
    countclone->Write();
    for(int i=0;i<nhist;i++) {
        vv[i]->Write();
    }
    for(int i=0;i<nhist2;i++) {
        vv2[i]->Write();
    }

    return;
}

BTagCalibrationStandalone calib_("CSVv2", "BtagFiles/CSVv2_Moriond17_mistag_G_H.csv");
BTagCalibrationStandaloneReader readerB_(BTagEntryStandalone::OP_LOOSE, "central");  // nominal
BTagCalibrationStandaloneReader readerL_(BTagEntryStandalone::OP_LOOSE, "central");  // nominal
// BTagCalibrationStandaloneReader readerB_(BTagEntryStandalone::OP_MEDIUM, "central");  // nominal
// BTagCalibrationStandaloneReader readerL_(BTagEntryStandalone::OP_MEDIUM, "central");  // nominal
// BTagCalibrationStandaloneReader readerB_(BTagEntryStandalone::OP_TIGHT, "central");  // nominal
// BTagCalibrationStandaloneReader readerL_(BTagEntryStandalone::OP_TIGHT, "central");  // nominal
TRandom *r0 = new TRandom();

void QCDhists::initialize() {
    readerB_.load(calib_, BTagEntryStandalone::FLAV_B,"comb");
    readerL_.load(calib_, BTagEntryStandalone::FLAV_UDSG,"incl");
}

TH1F* QCDhists::HistMan(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname) {

    std::string inputfile;


    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=dirname+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
            std::cout << inputfile << std::endl;
            TFile* in = new TFile(inputfile.c_str());
            if (in->IsZombie()) {in->Close(); continue;}
            if (!in->GetListOfKeys()->Contains(thisHIST.c_str())) return (new TH1F(thisHIST.c_str(),"dummy empty hist",10,0.,10.));

            if(j==0) {
                std::cout<<" adding up histos within a bin"<<std::endl;
                sum[i] = *(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone()));
            } else {
                TH1F* tmp = static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone());
                sum[i].Add(tmp);
            }
            in->Close();
        }
    }

    if(donorm) {
        // reweight to int lum
        std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
        for(int i=0;i<nbin;i++) {
            // get total number of events before filter
            float ntotal = norm[i];
            std::cout<<" for bin "<<i<<" number of pretrigger events is "<<ntotal<<std::endl;
            float fileLum= ntotal/xsec[i];
            std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
            outnorm[i] = goalintlum/fileLum;
            std::cout<<" scaling by a factor of "<<outnorm[i]<<std::endl;
            sum[i].Scale(outnorm[i]);
        }
    }

    //add the bins
    std::cout<<" adding bins"<<std::endl;
    TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone());
    for(int i=1;i<nbin;i++) {
        SUM->Add(&sum[i]);
    }


    return SUM;
}

TH2F* QCDhists::HistMan2(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname) {

    std::string inputfile;


    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=dirname+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
            TFile* in = new TFile(inputfile.c_str());
            if (in->IsZombie()) {in->Close(); continue;}
            if(j==0) {
                sum[i] = *(static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone()));
            } else {
                TH2F* tmp = static_cast<TH2F*>(in->Get(thisHIST.c_str())->Clone());
                sum[i].Add(tmp);
            }
            in->Close();
        }
    }

    if(donorm) {
        // reweight to int lum
        std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
        for(int i=0;i<nbin;i++) {
            // get total number of events before filter
            float ntotal = norm[i];
            std::cout<<" for bin "<<i<<" number of pretrigger events is "<<ntotal<<std::endl;
            float fileLum= ntotal/xsec[i];
            std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
            outnorm[i] = goalintlum/fileLum;
            std::cout<<" scaling by a factor of "<<outnorm[i]<<std::endl;
            sum[i].Scale(outnorm[i]);
        }
    }


    //add the bins
    std::cout<<" adding bins"<<std::endl;
    TH2F* SUM=static_cast<TH2F*>((sum[0]).Clone());
    for(int i=1;i<nbin;i++) {
        SUM->Add(&sum[i]);
    }


    return SUM;
}

void QCDhists::HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* samplenames,std::string dirname) {

    std::cout<<"entering HistNorm"<<std::endl; 

    std::string inputfile;
    TFile * in;

    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=dirname+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
            std::cout<<i<<" "<<j<<" "<<inputfile<<std::endl;
            in = new TFile(inputfile.c_str());
            if (in->IsZombie()) {in->Close(); continue;}
            if(j==0) {
                sum[i] = *(static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone()));
            } else {
                TH1F* tmp = static_cast<TH1F*>(in->Get("eventCountPreTrigger")->Clone());
                sum[i].Add(tmp);
            }
            in->Close();
        }
    }

    // reweight to int lum

    for(int i=0;i<nbin;i++) {
        // get total number of events before filter
        norm[i] = sum[i].GetBinContent(2);
        std::cout<<"norm "<<i<<" "<<norm[i]<<std::endl;
    }


    return;
}


TH1F* QCDhists::HistMerge(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* samplenames,bool donorm,std::string dirname) {

    std::string inputfile;

    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    TList *listO = new TList;
    for(int i=0;i<nbin;i++) {  // for each bin
        TList *listI = new TList;
        int ij=0;
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=dirname+samplenames[i]+"/histos"+samplenames[i]+"_"+std::to_string(j)+".root";
            TFile* in = new TFile(inputfile.c_str());
            if (in->IsZombie()) {in->Close(); continue;}
            if (!in->GetListOfKeys()->Contains(thisHIST.c_str())) return (new TH1F(thisHIST.c_str(),"dummy empty hist",10,0.,10.));
            if (ij==0) sum[i] = *(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone("h")));
            listI->Add(static_cast<TH1F*>(in->Get(thisHIST.c_str())->Clone()));
            in->Close();
            ij++;
        }
        sum[i].Reset();
        std::cout << "HERE0" << std::endl;
        sum[i].Merge(listI);
        std::cout << "HERE1" << std::endl;
    }

    if(donorm) {
        // reweight to int lum
        std::cout<<" reweighting to inst lum of "<<goalintlum<<" for each bin"<<std::endl;
        for(int i=0;i<nbin;i++) {
            // get total number of events before filter
            float ntotal = norm[i];
            std::cout<<" for bin "<<i<<" number of pretrigger events is "<<ntotal<<std::endl;
            float fileLum= ntotal/xsec[i];
            std::cout<<" equ lum for bin is "<<fileLum<<" fb-1"<<std::endl;
            outnorm[i] = goalintlum/fileLum;
            std::cout<<" scaling by a factor of "<<outnorm[i]<<std::endl;
            sum[i].Scale(outnorm[i]);
        }
    }

    //add the bins
    std::cout<<" adding bins"<<std::endl;

    for(int i=1;i<nbin;i++) {
        listO->Add(&sum[i]);
    }
    TH1F* SUM=static_cast<TH1F*>((sum[0]).Clone("SUM"));
    SUM->Reset();
    SUM->Merge(listO);

    return SUM;
}

double QCDhists::fakerate(double jet_pt, double jet_eta, int jet_nTrack, int varType){
    double fakerate = 0.0;

    if (varType == 1) {//alpha (tracksource = 0, trackquality highPurity, pvWeight>0)
        if( jet_pt>=100 && jet_pt<150 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.411328136921;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.110956951976;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.040011562407;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0200386364013;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0114627024159;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00724920025095;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00589098641649;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00409679533914;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00336773553863;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00276778638363;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00236453558318;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00189406401478;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00185049965512;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00158481043763;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00145675137173;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00138681638055;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00134717230685;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00123968662228;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00144343671855;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00160576647613;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00139793439303;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00170455907937;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00135290774051;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0013701407006;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00147025857586;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00149507762399;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00060332933208;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.002577773761;
            else if( jet_nTrack>=38 ) fakerate = 0.00587211269885;
        }
        else if( jet_pt>=150 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.45132869482;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.101133942604;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0410790778697;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0210691113025;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.012582000345;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00855491217226;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00611928943545;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00495016062632;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00387488841079;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0029952081386;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00244221882895;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00210818881169;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00169900327455;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0016031388659;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0016477368772;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00125536310952;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00129489053506;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00118776073214;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00122364691924;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00110992544796;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0010274410015;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00095766002778;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00128242862411;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00115970673505;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00111562979873;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00131704146042;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00149643281475;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.002198030008;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00244665308855;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00127782102209;
            else if( jet_nTrack>=46 ) fakerate = 0.00909058563411;
        }
        else if( jet_pt>=200 && jet_pt<250 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.404711425304;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0976014062762;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0454151332378;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0210127774626;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0136418100446;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00967902503908;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00740078045055;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00637179519981;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00436860881746;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00387860834599;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00280092400499;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0026380897034;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00235352898017;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00188532669563;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0015695883194;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00157465902157;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00150334683713;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00157660653349;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00127295986749;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00111181696411;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00129996147007;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00121080689132;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0013805431081;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00109967123717;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.000975408998784;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.000885513203684;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.000963868689723;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00105696090031;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0019351686351;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.000475604232633;
            else if( jet_nTrack>=46 ) fakerate = 0.00167973421048;
        }
        else if( jet_pt>=250 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.400169342756;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0994481369853;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0399819202721;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0219794176519;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0142903169617;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0102302450687;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00747192790732;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00598222529516;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00522055942565;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00407708110288;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00329083902761;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00272695790045;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00235966220498;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00184452848043;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00185582612175;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00167160248384;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00154473213479;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00164583197329;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00151590362657;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00133423495572;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.000874941179063;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00113894836977;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0013736380497;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00099275528919;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00128106120974;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00117876229342;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.000909059250262;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.000724673038349;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.000521186098922;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00221360684372;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00117480754852;
            else if( jet_nTrack>=50 ) fakerate = 0.00163900339976;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.387598574162;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.103158503771;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0389530956745;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.025333231315;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0141781233251;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00992115493864;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00793446693569;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00655185151845;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00532865710557;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00469963299111;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00419264798984;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00319076376036;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00282761850394;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00240248208866;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00224672211334;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00189148029312;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00193505862262;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00167585350573;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00168165867217;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00135798368137;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00169014523271;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00139151338954;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00116392422933;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00122381700203;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00147913338151;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00128917978145;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00141418108251;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00109150633216;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00108981470112;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00106288318057;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0019850989338;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.000513938721269;
            else if( jet_nTrack>=60 ) fakerate = 0.0102659985423;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.216925904155;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0888387784362;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0363200306892;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0218701343983;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0138519853354;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0104828095064;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00837497971952;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00713236350566;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00606932863593;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00481897918507;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00486017158255;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00402672030032;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00332744116895;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00308391521685;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00271835410967;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00240200012922;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00232656858861;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00221279705875;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00208962056786;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0016031927662;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00154155795462;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00144156778697;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00147982721683;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00151919363998;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00151698605623;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00136860064231;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00134961248841;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00143579591531;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.000827852229122;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00134823273402;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00107333355118;
            else if( jet_nTrack>=50 ) fakerate = 0.00203973311;
        }
        else if( jet_pt>=500 && jet_pt<700 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.376151800156;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0868468135595;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.035598449409;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0229212380946;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0142381386831;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0109982527792;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00950829312205;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00778950750828;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0068740863353;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00605094013736;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00545303616673;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00465875724331;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00427696295083;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0036800140515;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00331812747754;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00304454704747;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00283369026147;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00236481241882;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00217493646778;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00192248390522;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00197803368792;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00191504415125;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00184300413821;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00169333920348;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00151543680113;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00154388754163;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00144024461042;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00128579442389;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00135043927003;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00127844105009;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00116159627214;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00127593008801;
            else if( jet_nTrack>=60 ) fakerate = 0.000503823743202;
        }
        else if( jet_pt>=700 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.267393410206;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0866660252213;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0433187969029;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0251880865544;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0164148863405;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.013462588191;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0108908601105;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00977971404791;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00873505417258;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00782956369221;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00703095179051;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00587085355073;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00545358890668;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00522325839847;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00449662003666;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00438663177192;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00392517028376;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00353423459455;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00334319029935;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00275669223629;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00288099888712;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00265783490613;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00237119430676;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00235787499696;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0022879501339;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00201869383454;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00186536891852;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00173735222779;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00171718758065;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00148752564564;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00148196378723;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00130013783928;
            else if( jet_nTrack>=60 ) fakerate = 0.0013450642582;
        }
    }
    else {//alpha2Dsig (tracksource=0, trackQuality HighPurity, track ipXYsig<4)
        if( jet_pt>=100 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0677589178085;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0584921687841;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0434360653162;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.030286507681;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0215533468872;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0175257064402;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0123068392277;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0110569149256;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00865289010108;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00668580783531;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00518961623311;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00365586904809;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00420937594026;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00243917782791;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00183067307808;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.00197910983115;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00123751547653;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.000416403345298;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.000500752590597;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.000224228337174;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 1.97723202291e-05;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 0.0;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=200 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.113108091056;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0587297044694;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0477516464889;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0354919433594;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0295128300786;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0271997079253;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0223652496934;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0186233781278;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0156621076167;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0123292161152;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00974194146693;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00718888407573;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00714933639392;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00626720581204;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00505017722026;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.00287377391942;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00193187827244;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.00110103562474;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.00144552567508;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.000448211299954;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 0.000859102408867;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 3.95192801079e-05;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.110977262259;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0547714754939;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0439901910722;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0369902215898;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0313876010478;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.026527710259;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0258763767779;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0221502240747;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.018551832065;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0164150875062;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0132394433022;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0104881953448;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00969817768782;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00730892596766;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00694784848019;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.0053078760393;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00314273918048;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.00417654309422;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.00156928133219;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.000763362564612;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 0.000623483967502;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 0.00017194612883;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0676606819034;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0514590106905;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0371716059744;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0300994291902;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.026263512671;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0317170806229;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.022155592218;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0210395660251;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0188111457974;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0169373173267;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.012868299149;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0149185825139;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0123065654188;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00943822599947;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00862506590784;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.00719480169937;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00507904915139;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.00342770968564;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.00332551216707;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.00164829485584;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 0.000904305721633;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 0.000225199313718;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=500 && jet_pt<700 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0462760552764;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0484475977719;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0411542616785;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0307706091553;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0288863796741;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0241593420506;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0240439511836;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0216179843992;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0184720549732;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0166816432029;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0167829915881;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0148193519562;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0122329294682;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0111497864127;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00999608263373;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.0079288603738;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00645595556125;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.00497179245576;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.00356202735566;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.00232572574168;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 0.00190258270595;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 0.000663979910314;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.000107081490569;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.000327471614582;
        }
        else if( jet_pt>=700 && jet_pt<1000 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0634981170297;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0450021922588;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0383410677314;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0310143921524;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0303523428738;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0263274647295;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0253953244537;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0239105354995;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.022131152451;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0203762669116;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0189634338021;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0182483512908;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0167915876955;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0144937383011;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0138714052737;
            else if( jet_nTrack>=21 && jet_nTrack<23 ) fakerate = 0.0121630979702;
            else if( jet_nTrack>=23 && jet_nTrack<25 ) fakerate = 0.00995420850813;
            else if( jet_nTrack>=25 && jet_nTrack<27 ) fakerate = 0.00799527484924;
            else if( jet_nTrack>=27 && jet_nTrack<29 ) fakerate = 0.00635117152706;
            else if( jet_nTrack>=29 && jet_nTrack<32 ) fakerate = 0.00470072636381;
            else if( jet_nTrack>=32 && jet_nTrack<38 ) fakerate = 0.00316254841164;
            else if( jet_nTrack>=38 && jet_nTrack<50 ) fakerate = 0.0015876092948;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.000807894044556;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.000572000455577;
        }
    }

    return fakerate;

}

double QCDhists::fakerateTP(double jet_pt, double jet_eta, int jet_nTrack, int varType){
    double fakerate = 0.0;

    if (varType == 1) {//alpha (tracksource = 0, trackquality highPurity, pvWeight>0)
        if( jet_pt>=100 && jet_pt<150 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.411328136921;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.110956951976;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.040011562407;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0200386364013;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0114627024159;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00724920025095;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00589098641649;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00409679533914;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00336773553863;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00276778638363;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00236453558318;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00189406401478;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00185049965512;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00158481043763;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00145675137173;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00138681638055;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00134717230685;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00123968662228;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00144343671855;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00160576647613;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00139793439303;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00170455907937;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00135290774051;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0013701407006;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00147025857586;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00149507762399;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00060332933208;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.002577773761;
            else if( jet_nTrack>=38 ) fakerate = 0.00587211269885;
        }
        else if( jet_pt>=150 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.45132869482;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.101133942604;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0410790778697;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0210691113025;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.012582000345;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00855491217226;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00611928943545;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00495016062632;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00387488841079;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0029952081386;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00244221882895;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00210818881169;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00169900327455;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0016031388659;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0016477368772;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00125536310952;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00129489053506;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00118776073214;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00122364691924;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00110992544796;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0010274410015;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00095766002778;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00128242862411;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00115970673505;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00111562979873;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00131704146042;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00149643281475;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.002198030008;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00244665308855;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00127782102209;
            else if( jet_nTrack>=46 ) fakerate = 0.00909058563411;
        }
        else if( jet_pt>=200 && jet_pt<250 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.404711425304;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0976014062762;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0454151332378;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0210127774626;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0136418100446;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00967902503908;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00740078045055;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00637179519981;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00436860881746;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00387860834599;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00280092400499;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0026380897034;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00235352898017;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00188532669563;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0015695883194;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00157465902157;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00150334683713;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00157660653349;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00127295986749;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00111181696411;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00129996147007;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00121080689132;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0013805431081;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00109967123717;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.000975408998784;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.000885513203684;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.000963868689723;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00105696090031;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0019351686351;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.000475604232633;
            else if( jet_nTrack>=46 ) fakerate = 0.00167973421048;
        }
        else if( jet_pt>=250 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.400169342756;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0994481369853;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0399819202721;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0219794176519;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0142903169617;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0102302450687;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00747192790732;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00598222529516;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00522055942565;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00407708110288;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00329083902761;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00272695790045;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00235966220498;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00184452848043;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00185582612175;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00167160248384;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00154473213479;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00164583197329;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00151590362657;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00133423495572;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.000874941179063;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00113894836977;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0013736380497;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00099275528919;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00128106120974;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00117876229342;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.000909059250262;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.000724673038349;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.000521186098922;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00221360684372;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00117480754852;
            else if( jet_nTrack>=50 ) fakerate = 0.00163900339976;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.387598574162;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.103158503771;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0389530956745;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.025333231315;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0141781233251;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00992115493864;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00793446693569;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00655185151845;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00532865710557;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00469963299111;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00419264798984;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00319076376036;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00282761850394;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00240248208866;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00224672211334;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00189148029312;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00193505862262;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00167585350573;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00168165867217;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00135798368137;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00169014523271;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00139151338954;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00116392422933;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00122381700203;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00147913338151;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00128917978145;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00141418108251;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00109150633216;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00108981470112;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00106288318057;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0019850989338;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.000513938721269;
            else if( jet_nTrack>=60 ) fakerate = 0.0102659985423;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.216925904155;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0888387784362;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0363200306892;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0218701343983;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0138519853354;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0104828095064;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00837497971952;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00713236350566;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00606932863593;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00481897918507;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00486017158255;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00402672030032;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00332744116895;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00308391521685;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00271835410967;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00240200012922;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00232656858861;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00221279705875;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00208962056786;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0016031927662;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00154155795462;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00144156778697;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00147982721683;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00151919363998;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00151698605623;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00136860064231;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00134961248841;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00143579591531;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.000827852229122;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00134823273402;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00107333355118;
            else if( jet_nTrack>=50 ) fakerate = 0.00203973311;
        }
        else if( jet_pt>=500 && jet_pt<700 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.376151800156;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0868468135595;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.035598449409;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0229212380946;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0142381386831;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0109982527792;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00950829312205;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00778950750828;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0068740863353;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00605094013736;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00545303616673;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00465875724331;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00427696295083;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0036800140515;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00331812747754;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00304454704747;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00283369026147;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00236481241882;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00217493646778;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00192248390522;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00197803368792;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00191504415125;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00184300413821;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00169333920348;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00151543680113;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00154388754163;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00144024461042;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00128579442389;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00135043927003;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00127844105009;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00116159627214;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00127593008801;
            else if( jet_nTrack>=60 ) fakerate = 0.000503823743202;
        }
        else if( jet_pt>=700 ){
            if( jet_nTrack>=0 && jet_nTrack<2 ) fakerate = 0.267393410206;
            else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0866660252213;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0433187969029;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0251880865544;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0164148863405;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.013462588191;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0108908601105;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00977971404791;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00873505417258;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00782956369221;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00703095179051;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00587085355073;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00545358890668;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00522325839847;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00449662003666;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00438663177192;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00392517028376;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00353423459455;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00334319029935;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00275669223629;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00288099888712;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00265783490613;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00237119430676;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00235787499696;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0022879501339;
            else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00201869383454;
            else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00186536891852;
            else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00173735222779;
            else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00171718758065;
            else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00148752564564;
            else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00148196378723;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00130013783928;
            else if( jet_nTrack>=60 ) fakerate = 0.0013450642582;
        }
    }
    else {//alpha2Dsig (tracksource=0, trackQuality HighPurity, track ipXYsig<4)
        if( jet_pt>=100 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.121065996587;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0790413767099;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0997192487121;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0591917671263;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0293076150119;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0264087319374;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0224130842835;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0183669030666;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0213598310947;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0174206923693;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00726037798449;
            else if( jet_nTrack>=17 && jet_nTrack<19 ) fakerate = 0.010630573146;
            else if( jet_nTrack>=19 && jet_nTrack<21 ) fakerate = 0.00367591320537;
            else if( jet_nTrack>=21 && jet_nTrack<24 ) fakerate = 0.00923292059451;
            else if( jet_nTrack>=24 && jet_nTrack<27 ) fakerate = 0.00130184134468;
            else if( jet_nTrack>=27 && jet_nTrack<32 ) fakerate = 0.00247625564225;
            else if( jet_nTrack>=32 && jet_nTrack<40 ) fakerate = 0.0;
            else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 0.0;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=200 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.31605091691;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.131103932858;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.128241091967;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0749500766397;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0626299083233;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0524722002447;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.042688947171;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0482788868248;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0379444025457;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0459026992321;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0169816799462;
            else if( jet_nTrack>=17 && jet_nTrack<19 ) fakerate = 0.0132858417928;
            else if( jet_nTrack>=19 && jet_nTrack<21 ) fakerate = 0.00981686450541;
            else if( jet_nTrack>=21 && jet_nTrack<24 ) fakerate = 0.00289142271504;
            else if( jet_nTrack>=24 && jet_nTrack<27 ) fakerate = 0.00541092362255;
            else if( jet_nTrack>=27 && jet_nTrack<32 ) fakerate = 0.000751811719965;
            else if( jet_nTrack>=32 && jet_nTrack<40 ) fakerate = 0.0;
            else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 0.0;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.181784406304;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.167348787189;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0949432253838;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0872481763363;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0701444298029;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0772162899375;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0575019456446;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0488184243441;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0337308980525;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0308888703585;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0398111864924;
            else if( jet_nTrack>=17 && jet_nTrack<19 ) fakerate = 0.0151519188657;
            else if( jet_nTrack>=19 && jet_nTrack<21 ) fakerate = 0.014994725585;
            else if( jet_nTrack>=21 && jet_nTrack<24 ) fakerate = 0.00596855627373;
            else if( jet_nTrack>=24 && jet_nTrack<27 ) fakerate = 0.003117055865;
            else if( jet_nTrack>=27 && jet_nTrack<32 ) fakerate = 0.00155009666923;
            else if( jet_nTrack>=32 && jet_nTrack<40 ) fakerate = 0.00164627935737;
            else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 0.0;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.00536637520418;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0684101954103;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0841860026121;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0590422004461;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0542525202036;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0604502037168;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0472973957658;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0419289618731;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0192966535687;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.023137755692;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.028977362439;
            else if( jet_nTrack>=17 && jet_nTrack<19 ) fakerate = 0.0309908036143;
            else if( jet_nTrack>=19 && jet_nTrack<21 ) fakerate = 0.0147506343201;
            else if( jet_nTrack>=21 && jet_nTrack<24 ) fakerate = 0.00558003596961;
            else if( jet_nTrack>=24 && jet_nTrack<27 ) fakerate = 0.00555703230202;
            else if( jet_nTrack>=27 && jet_nTrack<32 ) fakerate = 0.00172243244015;
            else if( jet_nTrack>=32 && jet_nTrack<40 ) fakerate = 0.000458588823676;
            else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 0.000292478565825;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
        else if( jet_pt>=500 && jet_pt<1000 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0300417393446;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0578078515828;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0600556693971;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0448701195419;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0499757453799;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0342401452363;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0464037358761;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.030430605635;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0289333984256;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0264333747327;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0178076997399;
            else if( jet_nTrack>=17 && jet_nTrack<19 ) fakerate = 0.0200277492404;
            else if( jet_nTrack>=19 && jet_nTrack<21 ) fakerate = 0.0173008013517;
            else if( jet_nTrack>=21 && jet_nTrack<24 ) fakerate = 0.0131498826668;
            else if( jet_nTrack>=24 && jet_nTrack<27 ) fakerate = 0.00763271423057;
            else if( jet_nTrack>=27 && jet_nTrack<32 ) fakerate = 0.00373955839314;
            else if( jet_nTrack>=32 && jet_nTrack<40 ) fakerate = 0.00151756347623;
            else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 0.000578296836466;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.000300995947327;
            else if( jet_nTrack>=60 && jet_nTrack<70 ) fakerate = 0.0;
        }
    }

    return fakerate;

}

//[MC] FakeRate
double QCDhists::fakerateF(double jet_pt, double jet_eta, int jet_nTrack, int varType, int flav, std::string cutset_){
    double fakerate = 0.0;
    double error = 0.0;

    if (cutset_=="X"){
        //~/Dropbox/UMD\ Analysis/FakeRate/20171115/fakerate_QCDMC_3DSig_8_Med_8.txt
        //MC fake rate bin in nTrk, flavor; alpha3dsig only; D1 cuts (5mm); a3Dsig cut=0.25
        if (fabs(flav)==5 || fabs(flav)==8) {// default: b and g->bb

            if( jet_nTrack>=0.0 && jet_nTrack<4.0 ) {
                fakerate = 0.0764036476612;
            }
            if( jet_nTrack>=4.0 && jet_nTrack<8.0 ) {
                fakerate = 0.0240180138499;
            }
            if( jet_nTrack>=8.0 && jet_nTrack<12.0 ) {
                fakerate = 0.00735678197816;
            }
            if( jet_nTrack>=12.0 && jet_nTrack<16.0 ) {
                fakerate = 0.00186479813419;
            }
            if( jet_nTrack>=16.0 && jet_nTrack<20.0 ) {
                fakerate = 0.000393063208321;
            }
            if( jet_nTrack>=20.0 && jet_nTrack<24.0 ) {
                fakerate = 0.000105553546746;
            }
            if( jet_nTrack>=24.0 && jet_nTrack<40.0 ) {
                fakerate = 1.39190433401e-06;
            }
            if( jet_nTrack>=40.0 && jet_nTrack<80.0 ) {
                fakerate = 0.0;
            }
        }
        else {//g,u,d,c,s
            if( jet_nTrack>=0.0 && jet_nTrack<4.0 ) {
                fakerate = 0.0272194761783;
            }
            if( jet_nTrack>=4.0 && jet_nTrack<8.0 ) {
                fakerate = 0.0024541572202;
            }
            if( jet_nTrack>=8.0 && jet_nTrack<12.0 ) {
                fakerate = 0.000505852687638;
            }
            if( jet_nTrack>=12.0 && jet_nTrack<16.0 ) {
                fakerate = 0.00015330662427;
            }
            if( jet_nTrack>=16.0 && jet_nTrack<20.0 ) {
                fakerate = 5.15911160619e-05;
            }
            if( jet_nTrack>=20.0 && jet_nTrack<24.0 ) {
                fakerate = 2.135487739e-05;
            }
            if( jet_nTrack>=24.0 && jet_nTrack<40.0 ) {
                fakerate = 6.16831357547e-06;
            }
            if( jet_nTrack>=40.0 && jet_nTrack<80.0 ) {
                fakerate = 0.0;
            }

        }
    }//end of cutsetX
    else if (cutset_=="1"){
        //std::cout<<"[MC] CutSet = 1"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset1.txt
        //cutset1
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.091007791;
                error = 0.056625968;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.034747239;
                error = 0.010749298;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.035219494;
                error = 0.025553014;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000241418;
                error = 0.003399285;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.006634058;
                error = 0.005712825;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.011468501;
                error = 0.001658203;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002515812;
                error = 0.000401359;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000047768;
                error = 0.000790234;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000293478;
                error = 0.000174225;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset1
    else if (cutset_=="2"){
        //std::cout<<"[MC] CutSet = 2"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset2.txt
        //cutset2
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.017776422;
                error = 0.059319302;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.002533249;
//                 error = 0.000382620;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.002772487;
                error = 0.005547458;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.000318333;
//                 error = 0.000133305;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.009511134;
                error = 0.009995830;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.000078965;
//                 error = 0.000080139;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.009646106;
                error = 0.001920828;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002533249;
                error = 0.000382620;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000509282;
                error = 0.000212996;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000318333;
                error = 0.000133305;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.000078965;
//                 error = 0.000080139;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000078965;
                error = 0.000080139;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset2
    else if (cutset_=="3"){
        //std::cout<<"[MC] CutSet = 3"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset3.txt
        //cutset3
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.005435236;
                error = 0.006502964;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000474873;
                error = 0.002526085;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.005688290;
                error = 0.001619985;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000778228;
                error = 0.000277043;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000138929;
                error = 0.000096240;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000066653;
                error = 0.000052848;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000068944;
                error = 0.000069848;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset3
    else if (cutset_=="4"){
        //std::cout<<"[MC] CutSet = 4"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset4.txt
        //cutset4
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.001413448;
                error = 0.048206738;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.004730882;
                error = 0.004399605;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.004409454;
                error = 0.004343413;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.004691387;
                error = 0.001487461;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000344542;
                error = 0.000177416;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset4
    else if (cutset_=="5"){
        //std::cout<<"[MC] CutSet = 5"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset5.txt
        //cutset5
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.106563225;
                error = 0.068887879;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.034559544;
                error = 0.010682442;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.035705920;
                error = 0.025923635;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000246726;
                error = 0.003036555;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.009337957;
                error = 0.008216448;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.011129123;
                error = 0.001893431;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002525200;
                error = 0.000398519;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000015756;
                error = 0.000814237;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000293273;
                error = 0.000161504;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset5
    else if (cutset_=="6"){
        //std::cout<<"[MC] CutSet = 6"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset6.txt
        //cutset6
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.005767023;
                error = 0.005732841;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.001942574;
                error = 0.003015861;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012686262;
                error = 0.004026754;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002139245;
                error = 0.000328302;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000242342;
                error = 0.000208354;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset6
    else if (cutset_=="7"){
        //std::cout<<"[MC] CutSet = 7"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset7.txt
        //cutset7
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.085752450;
                error = 0.092576154;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.005402870;
                error = 0.005250144;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.003395952;
                error = 0.004667172;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.005398412;
                error = 0.003012322;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001008174;
                error = 0.000246674;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000141152;
                error = 0.000179217;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000054247;
                error = 0.000039226;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset7
    else if (cutset_=="8"){
        //std::cout<<"[MC] CutSet = 8"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset8.txt
        //cutset8
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.082552657;
                error = 0.070909274;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.017333694;
                error = 0.010369442;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.004979199;
                error = 0.006896357;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.001863467;
                error = 0.003561052;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.008585087;
                error = 0.002023294;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001776625;
                error = 0.000400376;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000364279;
                error = 0.000253261;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000100283;
                error = 0.000154937;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset8
    else if (cutset_=="9"){
        //std::cout<<"[MC] CutSet = 9"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180201r1/FakeRate_cutset9.txt
        //cutset9
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.140636355;
                error = 0.095017506;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.029160412;
                error = 0.014766987;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.018177882;
                error = 0.014256610;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.014818208;
                error = 0.012781939;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.014810834;
                error = 0.002667873;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.005120015;
                error = 0.000621572;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.001360814;
                error = 0.000501770;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000159184;
                error = 0.000502774;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000195149;
                error = 0.000141333;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000250109;
                error = 0.000253132;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset9
    else if (cutset_=="10"){
        //std::cout<<"[MC] CutSet = 10"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180222_cut10/FakeRate_cutset10.txt
        //cutset10
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.061317001;
                error = 0.090330640;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002033191;
                error = 0.008962274;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000432861;
                error = 0.011457711;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000432861;
                error = 0.011457711;
//                 fakerate = 0.014920635;
//                 error = 0.012054389;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.007742672;
                error = 0.002782491;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002660405;
                error = 0.000440500;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.001222754;
                error = 0.000256548;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000535164;
                error = 0.000256586;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000592041;
                error = 0.000576590;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset10 (=cutset 3 except a3d cut =0.5)
    else if (cutset_=="11" || cutset_=="11a"){
        //std::cout<<"[MC] CutSet = 11"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180227_cut11/FakeRate_cutset11.txt
        //cutset11
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.063960016;
                error = 0.075843471;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.001469939;
                error = 0.007620897;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.007549996;
                error = 0.008793359;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.010295879;
                error = 0.002199902;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002575297;
                error = 0.000380299;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000903147;
                error = 0.000281890;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000573658;
                error = 0.000199477;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000074672;
                error = 0.000075485;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset11 (=cutset 4 except a3d cut =0.5 & medID cut=0.10)
    else if (cutset_=="21"){
        //std::cout<<"[MC] CutSet = 21"<<std::endl;
        //FakeRateHistograms_20181029/fakerate_p20181030/FakeRate_cutset1.txt
        //cutset21 == cutset1 2017
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.085364841;
                error = 0.052403670;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.032302987;
                error = 0.009883264;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.036249205;
                error = 0.026352946;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000249382;
                error = 0.002855125;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.008695275;
                error = 0.007620107;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.023867251;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.671910291;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.011556502;
                error = 0.001602595;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002607194;
                error = 0.000373927;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000003050;
                error = 0.000824383;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000293200;
                error = 0.000157078;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000999705;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.028620991;
            }
        }
    }//end of cutset21 (=cutset 1 except 2017 FR)
    else{
        std::cout<<"CutSet["<<cutset_<<"] not defined!"<<std::endl;
    }

//     if (fakerate>0.0) {
//         r0->SetSeed(ULong_t(jet_pt));
//         double tmpfr=fakerate-r0->Gaus(0,error);
//         if (tmpfr>0.0) return tmpfr;
//         else return 0.0;
//     }
//     else return fakerate;

    return fakerate;
}//End of fakerateF

Int_t QCDhists::GetNtrkBin( const int jet_nTrack ){
    Int_t bin=-1;
    if( jet_nTrack>=0 && jet_nTrack<4  ) bin = 1;
    if( jet_nTrack>=4 && jet_nTrack<8  ) bin = 2;
    if( jet_nTrack>=8 && jet_nTrack<12  ) bin = 3;
    if( jet_nTrack>=12 && jet_nTrack<16  ) bin = 4;
    if( jet_nTrack>=16 && jet_nTrack<20  ) bin = 5;
    if( jet_nTrack>=20 && jet_nTrack<24  ) bin = 6;
    if( jet_nTrack>=24 && jet_nTrack<40  ) bin = 7;
    if( jet_nTrack>=40 ) bin = 8;
    return bin;
}

// gamma+jets data
double QCDhists::fakerateFD(double jet_pt, double jet_eta, int jet_nTrack, int varType, int flav, std::string cutSet){
    double fakerate = 0.0;
    double error = 0.0;

    if (cutSet=="X"){
        //~/Dropbox/UMD\ Analysis/FakeRate/20171115/fakerate_QCDMC_3DSig_8_Med_8.txt
        //MC fake rate bin in nTrk, flavor; alpha3dsig only; D1 cuts (5mm); a3Dsig cut=0.25
        if (fabs(flav)==5 || fabs(flav)==8) {// default: b and g->bb

            if( jet_nTrack>=0.0 && jet_nTrack<4.0 ) {
                fakerate = 0.0764036476612;
            }
            if( jet_nTrack>=4.0 && jet_nTrack<8.0 ) {
                fakerate = 0.0240180138499;
            }
            if( jet_nTrack>=8.0 && jet_nTrack<12.0 ) {
                fakerate = 0.00735678197816;
            }
            if( jet_nTrack>=12.0 && jet_nTrack<16.0 ) {
                fakerate = 0.00186479813419;
            }
            if( jet_nTrack>=16.0 && jet_nTrack<20.0 ) {
                fakerate = 0.000393063208321;
            }
            if( jet_nTrack>=20.0 && jet_nTrack<24.0 ) {
                fakerate = 0.000105553546746;
            }
            if( jet_nTrack>=24.0 && jet_nTrack<40.0 ) {
                fakerate = 1.39190433401e-06;
            }
            if( jet_nTrack>=40.0 && jet_nTrack<80.0 ) {
                fakerate = 0.0;
            }
        }
        else {//g,u,d,c,s
            if( jet_nTrack>=0.0 && jet_nTrack<4.0 ) {
                fakerate = 0.0272194761783;
            }
            if( jet_nTrack>=4.0 && jet_nTrack<8.0 ) {
                fakerate = 0.0024541572202;
            }
            if( jet_nTrack>=8.0 && jet_nTrack<12.0 ) {
                fakerate = 0.000505852687638;
            }
            if( jet_nTrack>=12.0 && jet_nTrack<16.0 ) {
                fakerate = 0.00015330662427;
            }
            if( jet_nTrack>=16.0 && jet_nTrack<20.0 ) {
                fakerate = 5.15911160619e-05;
            }
            if( jet_nTrack>=20.0 && jet_nTrack<24.0 ) {
                fakerate = 2.135487739e-05;
            }
            if( jet_nTrack>=24.0 && jet_nTrack<40.0 ) {
                fakerate = 6.16831357547e-06;
            }
            if( jet_nTrack>=40.0 && jet_nTrack<80.0 ) {
                fakerate = 0.0;
            }

        }
    }//end of cutsetX
    else if (cutSet=="1"){
        //std::cout<<"[Data] CutSet = 1"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset1.txt
        //cutset1
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.044195049;
                error = 0.062926506;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.066300757;
                error = 0.011223864;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.025045995;
                error = 0.007403220;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.005860697;
                error = 0.007175669;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.015725104;
                error = 0.001647199;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002351884;
                error = 0.000400395;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000344560;
                error = 0.000323935;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000211086;
                error = 0.000396216;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000417405;
                error = 0.000301608;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }

    }//end of cutset1
    else if (cutSet=="2"){
        //std::cout<<"[Data] CutSet = 2"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset2.txt
        //cutset2
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.013937147;
//                 error = 0.003622609;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002611315;
                error = 0.008099194;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.012256928;
                error = 0.007141428;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
//                 fakerate = 0.000467706;
//                 error = 0.000182418;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.013937147;
                error = 0.003622609;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002585000;
                error = 0.000351653;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000236379;
                error = 0.000306443;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000467706;
                error = 0.000182418;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset2
    else if (cutSet=="3"){
        //std::cout<<"[Data] CutSet = 3"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset3.txt
        //cutset3
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.002884209;
                error = 0.003948002;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.024049988;
                error = 0.085467658;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001426918;
                error = 0.000248509;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000280365;
                error = 0.000176596;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000101748;
                error = 0.000072359;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000344517;
                error = 0.000363641;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset3
    else if (cutSet=="4"){
        //std::cout<<"[Data] CutSet = 4"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset4.txt
        //cutset4
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.003019468;
                error = 0.003124473;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.005227647;
                error = 0.002085540;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000747301;
                error = 0.000169026;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset4
    else if (cutSet=="5"){
        //std::cout<<"[Data] CutSet = 5"<<std::endl;
        //FakeRateHistograms_v2/fakerate_p20171211r1/FakeRate_cutset5.txt
        //cutset5
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.044195049;
                error = 0.062926506;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.066300757;
                error = 0.011223864;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.025045995;
                error = 0.007403220;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.005860697;
                error = 0.007175669;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.015725104;
                error = 0.001647199;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002351884;
                error = 0.000400395;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000344560;
                error = 0.000323935;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000211086;
                error = 0.000396216;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000417405;
                error = 0.000301608;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset5
    else if (cutSet=="6"){
        //std::cout<<"[Data] CutSet = 6"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset6.txt
        //cutset6
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.004874171;
                error = 0.007225038;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.008500207;
                error = 0.005179613;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012395947;
                error = 0.003749585;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001960682;
                error = 0.000306980;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000238004;
                error = 0.000140176;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset6
    else if (cutSet=="7"){
        //std::cout<<"[Data] CutSet = 7"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset7.txt
        //cutset7
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.021617584;
                error = 0.008291794;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000469282;
                error = 0.003321747;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.009204083;
                error = 0.003265770;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000780240;
                error = 0.000294516;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000275936;
                error = 0.000165483;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000238004;
                error = 0.000140176;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset7
    else if (cutSet=="8"){
        //std::cout<<"[Data] CutSet = 8"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180109r1/FakeRate_cutset8.txt
        //cutset8
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.016383125;
                error = 0.144820405;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.041638434;
                error = 0.012032606;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.004972943;
                error = 0.006165526;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.001672277;
                error = 0.005924211;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012647440;
                error = 0.003912600;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001743738;
                error = 0.000431009;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000595314;
                error = 0.000292252;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000191787;
                error = 0.000344736;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset8
    else if (cutSet=="9"){
        //std::cout<<"[Data] CutSet = 9"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180201r1/FakeRate_cutset9.txt
        //cutset9
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.046485364;
                error = 0.182139700;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.063739181;
                error = 0.016713895;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.037514411;
                error = 0.013565677;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.015383565;
                error = 0.014730966;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.019020380;
                error = 0.004923453;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.004988618;
                error = 0.000626857;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000913211;
                error = 0.000604533;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000516047;
                error = 0.000832754;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.001137016;
                error = 0.000932514;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset9 (=cutset 1 except a3d cut =0.4)
    else if (cutSet=="10"){
        //std::cout<<"[Data] CutSet = 10"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180222_cut10/FakeRate_cutset10.txt
        //cutset10
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.014767224;
                error = 0.008692570;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.008138913;
                error = 0.029221028;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012542081;
                error = 0.027295911;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.003518312;
                error = 0.000431324;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.001209196;
                error = 0.000384633;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.001321519;
                error = 0.000324626;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000197349;
                error = 0.002066180;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset10 (=cutset 3 except a3d cut =0.5)
    else if (cutSet=="11" || cutSet=="11a"){
        //std::cout<<"[Data] CutSet = 11"<<std::endl;
        //FakeRateHistograms_Cut0105/fakerate_p20180227_cut11/FakeRate_cutset11.txt
        //cutset11
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.016085995;
                error = 0.010242065;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.012783732;
                error = 0.008152422;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000260225;
                error = 0.004937175;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012391529;
                error = 0.003500236;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.002809802;
                error = 0.000409830;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000652494;
                error = 0.000361470;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000313428;
                error = 0.000290809;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000923301;
                error = 0.000731900;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
        }
    }//end of cutset11 (=cutset 4 except a3d cut =0.5 & medID cut=0.10)
    else if (cutSet=="21"){
        //std::cout<<"[Data] CutSet = 21"<<std::endl;
        //FakeRateHistograms_20181029/fakerate_p20181030/FakeRate_cutset1.txt
        //cutset21 == cutset1 2017
        if (fabs(flav)==5 || fabs(flav)==8) {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.042640701;
                error = 0.014196423;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.014077252;
                error = 0.009043920;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.003157260;
                error = 0.005000407;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.068417704;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.034407973;
                error = 0.044103475;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 2.333466706;
            }
        }
        else {
            if( jet_nTrack>=1.00 && jet_nTrack<6.00 ) {
                fakerate = 0.012757472;
                error = 0.002188627;
            }
            else if( jet_nTrack>=6.00 && jet_nTrack<11.00 ) {
                fakerate = 0.001609106;
                error = 0.000529965;
            }
            else if( jet_nTrack>=11.00 && jet_nTrack<16.00 ) {
                fakerate = 0.000573322;
                error = 0.000469009;
            }
            else if( jet_nTrack>=16.00 && jet_nTrack<21.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=21.00 && jet_nTrack<24.00 ) {
                fakerate = 0.000000000;
                error = 0.007098993;
            }
            else if( jet_nTrack>=24.00 && jet_nTrack<40.00 ) {
                fakerate = 0.000000000;
                error = 0.000000000;
            }
            else if( jet_nTrack>=40.00 && jet_nTrack<80.00 ) {
                fakerate = 0.000000000;
                error = 0.244031412;
            }
        }
    }//end of cutset21 (=cutset 1 except 2017 FR)
    else {
        std::cout<<"CutSet["<<cutSet<<"] not defined!"<<std::endl;
    }

//     if (fakerate>0.0) {
//         r0->SetSeed(ULong_t(jet_pt));
//         double tmpfr=fakerate-r0->Gaus(0,error);
//         if (tmpfr>0.0) return tmpfr;
//         else return 0.0;
//     }
//     else return fakerate;

    return fakerate;
}//end of fakerateFD


// MC fake rate bin in eta, nTrk; alpha3dsig only; D1 cuts (5mm); a3Dsig cut=0.25
// // MC fake rate bin in nTrk only; alpha3dsig only
// double fakerateF(double jet_pt, double jet_eta, int jet_nTrack, int varType, int flav){
//     double fakerate = 0.0;
//     double error = 0.0;
//     if (fabs(flav)==5 || fabs(flav)==8) {// default: b and g->bb

//         if( fabs(jet_eta)<1.4 ){

//         }
//         else {

//         }


//     }
//     else {//g,u,d,c,s
//         if( fabs(jet_eta)<1.4 ){

//         }
//         else {

//         }


//     }
//     return fakerate;
// }

double ntrkrewgt(int nTrk){
    double wgt = 1.0;
//     if (nTrk<5) wgt = 0.02550;
//     else if (nTrk>=5  && nTrk<10) wgt = 0.10118;
//     else if (nTrk>=10 && nTrk<15) wgt = 0.30487;
//     else if (nTrk>=15 && nTrk<20) wgt = 0.25302;
//     else if (nTrk>=20 && nTrk<25) wgt = 0.22333;
//     else if (nTrk>=25 && nTrk<30) wgt = 0.25405;
//     else if (nTrk>=30 && nTrk<35) wgt = 0.50399;
//     else if (nTrk>=35 && nTrk<40) wgt = 0.34021;
//     else if (nTrk>=40 && nTrk<45) wgt = 1.12203;
//     else if (nTrk>=45 && nTrk<50) wgt = 0.36343;
//     else if (nTrk>=50) wgt = 0.35776;

    if (nTrk<5) wgt = 0.06009;
    else if (nTrk>=5  && nTrk<10) wgt = 0.30426;
    else if (nTrk>=10 && nTrk<15) wgt = 0.47897;
    else if (nTrk>=15 && nTrk<20) wgt = 0.37181;
    else if (nTrk>=20 && nTrk<25) wgt = 0.47475;
    else if (nTrk>=25 && nTrk<30) wgt = 0.74978;
    else if (nTrk>=30 && nTrk<35) wgt = 0.40890;
    else if (nTrk>=35 && nTrk<40) wgt = 0.87315;
    else if (nTrk>=40 && nTrk<45) wgt = 0.34007;
    else if (nTrk>=45 && nTrk<50) wgt = 0.18420;
    else if (nTrk>=50) wgt = 1.0;

    return wgt;
}


void QCDhists::fillFRPlots(std::string cutname, char * hnames[9],
                           vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                           vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                           vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                           vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                           double *frwgts, double scale, double ncomb)
{
    for (int i=0;i<4;i++) {
        double varWgt = scale*frwgts[i];
        int jdx0 = goodjetIdx[i];

        ((TH1F*)gDirectory->Get(TString(std::string(hnames[0])+"_"+cutname)))->Fill(jet_pt->at(jdx0),varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[1])+"_"+cutname)))->Fill(jet_eta->at(jdx0),varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[2])+"_"+cutname)))->Fill(jntrack[jdx0],varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[3])+"_"+cutname)))->Fill(jet_medipsig[jdx0],varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[4])+"_"+cutname)))->Fill(jet_logmedipsig[jdx0],varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[5])+"_"+cutname)))->Fill(jet_medtheta2D[jdx0],varWgt);
        ((TH1F*)gDirectory->Get(TString(std::string(hnames[6])+"_"+cutname)))->Fill(jet_logmedtheta2D[jdx0],varWgt);

        // 2D
        ((TH2F*)gDirectory->Get(TString(std::string(hnames[8])+"_"+cutname)))->Fill(jet_medtheta2D[jdx0],jet_medipsig[jdx0],varWgt);

        for (int j=0;j<4;j++) {
            if (i==j) continue;
            int jdx1 = goodjetIdx[j];
            double mass = sqrt(
                               pow((jet_e[jdx0]+jet_e[jdx1]),2) -
                               pow((jet_px[jdx0]+jet_px[jdx1]),2) -
                               pow((jet_py[jdx0]+jet_py[jdx1]),2) -
                               pow((jet_pz[jdx0]+jet_pz[jdx1]),2)
                               );
            //(TH1*)R__H(hnames[7])->Fill(mass,varWgt);//only work for TH1
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[7])+"_"+cutname)))->Fill(mass,varWgt/ncomb);
        }
    }
}

void QCDhists::fillFRPlots2(std::string cutname, char * hnames[9],
                            vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                            vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                            vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                            vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                            double *frwgts, double scale, double ncomb)
{
    int nthComb=0;
    for(int i0=0;i0<4;i0++) {
        int idx0 = goodjetIdx[i0];
        for(Int_t i1=0; i1<4; i1++) {
            if (i1==i0) continue;
            int idx1 = goodjetIdx[i1];
            for(Int_t i2=0; i2<4; i2++) {
                if (i2==i0 || i2==i1) continue;
//                 int idx2 = goodjetIdx[i2];
//                 for(Int_t i3=0; i3<4; i3++) {
//                     if (i3==i0 || i3==i1 || i3==i2) continue;
//                     int idx3 = goodjetIdx[i3];
//                 }
                double varWgt = scale*frwgts[nthComb];
                nthComb++;

                ((TH1F*)gDirectory->Get(TString(std::string(hnames[0])+"_"+cutname)))->Fill(jet_pt->at(idx0),varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[1])+"_"+cutname)))->Fill(jet_eta->at(idx0),varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[2])+"_"+cutname)))->Fill(jntrack[idx0],varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[3])+"_"+cutname)))->Fill(jet_medipsig[idx0],varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[4])+"_"+cutname)))->Fill(jet_logmedipsig[idx0],varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[5])+"_"+cutname)))->Fill(jet_medtheta2D[idx0],varWgt);
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[6])+"_"+cutname)))->Fill(jet_logmedtheta2D[idx0],varWgt);

                // 2D
                ((TH2F*)gDirectory->Get(TString(std::string(hnames[8])+"_"+cutname)))->Fill(jet_medtheta2D[idx0],jet_medipsig[idx0],varWgt);

                double mass = sqrt(
                                   pow((jet_e[idx0]+jet_e[idx1]),2) -
                                   pow((jet_px[idx0]+jet_px[idx1]),2) -
                                   pow((jet_py[idx0]+jet_py[idx1]),2) -
                                   pow((jet_pz[idx0]+jet_pz[idx1]),2)
                                   );
                //(TH1*)R__H(hnames[7])->Fill(mass,varWgt);//only work for TH1
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[7])+"_"+cutname)))->Fill(mass,scale*frwgts[nthComb]/ncomb);
            }
        }
    }
}

void QCDhists::fillFRPlots22(std::string cutname, char * hnames[9],
                             vector<float> *jet_pt, vector<float> *jet_eta, vector<int> &goodjetIdx,
                             vector<int> &jntrack, vector<float> &jet_medipsig, vector<float> &jet_logmedipsig,
                             vector<float> &jet_medtheta2D, vector<float> &jet_logmedtheta2D,
                             vector<float> &jet_e,vector<float> &jet_px,vector<float> &jet_py, vector<float> &jet_pz,
                             double *frwgts, double scale, double ncomb)
{
    int nthComb=0;
    for(int i0=0;i0<4;i0++) {
        int idx0 = goodjetIdx[i0];
        for(Int_t i1=i0+1; i1<4; i1++) {
            int idx1 = goodjetIdx[i1];
            double varWgt = scale*frwgts[nthComb];
            //std::cout<<"varWgt["<<nthComb<<"]="<<varWgt<<std::endl;
            //TString teststr = TString(std::string(hnames[0])+"_"+cutname);
            //std::cout << teststr << std::endl;
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[0])+"_"+cutname)))->Fill(jet_pt->at(idx0),varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[1])+"_"+cutname)))->Fill(jet_eta->at(idx0),varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[2])+"_"+cutname)))->Fill(jntrack[idx0],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[3])+"_"+cutname)))->Fill(jet_medipsig[idx0],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[4])+"_"+cutname)))->Fill(jet_logmedipsig[idx0],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[5])+"_"+cutname)))->Fill(jet_medtheta2D[idx0],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[6])+"_"+cutname)))->Fill(jet_logmedtheta2D[idx0],varWgt);

            ((TH1F*)gDirectory->Get(TString(std::string(hnames[0])+"_"+cutname)))->Fill(jet_pt->at(idx1),varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[1])+"_"+cutname)))->Fill(jet_eta->at(idx1),varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[2])+"_"+cutname)))->Fill(jntrack[idx1],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[3])+"_"+cutname)))->Fill(jet_medipsig[idx1],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[4])+"_"+cutname)))->Fill(jet_logmedipsig[idx1],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[5])+"_"+cutname)))->Fill(jet_medtheta2D[idx1],varWgt);
            ((TH1F*)gDirectory->Get(TString(std::string(hnames[6])+"_"+cutname)))->Fill(jet_logmedtheta2D[idx1],varWgt);

            // 2D
            ((TH2F*)gDirectory->Get(TString(std::string(hnames[8])+"_"+cutname)))->Fill(jet_medtheta2D[idx0],jet_medipsig[idx0],varWgt);
            ((TH2F*)gDirectory->Get(TString(std::string(hnames[8])+"_"+cutname)))->Fill(jet_medtheta2D[idx1],jet_medipsig[idx1],varWgt);

            for(Int_t i2=0; i2<4; i2++) {
                if (i2==i0 || i2==i1) continue;
                int idx2 = goodjetIdx[i2];
                double mass = sqrt(
                                   pow((jet_e[idx0]+jet_e[idx2]),2) -
                                   pow((jet_px[idx0]+jet_px[idx2]),2) -
                                   pow((jet_py[idx0]+jet_py[idx2]),2) -
                                   pow((jet_pz[idx0]+jet_pz[idx2]),2)
                                   );
                //(TH1*)R__H(hnames[7])->Fill(mass,varWgt);//only work for TH1
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[7])+"_"+cutname)))->Fill(mass,varWgt/ncomb);
                double mass1 = sqrt(
                                   pow((jet_e[idx1]+jet_e[idx2]),2) -
                                   pow((jet_px[idx1]+jet_px[idx2]),2) -
                                   pow((jet_py[idx1]+jet_py[idx2]),2) -
                                   pow((jet_pz[idx1]+jet_pz[idx2]),2)
                                   );
                //(TH1*)R__H(hnames[7])->Fill(mass1,varWgt);//only work for TH1
                ((TH1F*)gDirectory->Get(TString(std::string(hnames[7])+"_"+cutname)))->Fill(mass1,varWgt/ncomb);
            }
            //std::cout<<"End of "<<nthComb<<"-th comb"<<std::endl;
            nthComb++;
        }
    }
}

double QCDhists::frWeight(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<vector<float> > *track_pt, int njetscut, double jptcut, int varType){
    double p0 = 1.;
    double p1 = 0.;
    int njets = (njetscut == -1 ? jetpt->size() : std::min(njetscut,(int)jetpt->size()));
    if (verbose) {
        std::cout << "[frWeight] varType = " << varType << std::endl;
        std::cout << "[frWeight] njets = " << njets << std::endl;
        std::cout << "[frWeight] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        if (!basicjet->at(j) || jetpt->at(j)<jptcut) {
            if (verbose) std::cout << "[frWeight] Jet[" << j << "] pT = " << jetpt->at(j) << " doesn't pass basic jet!" << std::endl;
            continue;
        }
        const vector<float> *track_pts = &(track_pt->at(j));
        int ntrks1 = track_pts->size();
        if (verbose) std::cout << "[frWeight] Jet pT = " << jetpt->at(j) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(j),jeteta->at(j),ntrks1,varType));
        double p11 = fakerate(jetpt->at(j),jeteta->at(j),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            const vector<float> *track_pts1 = &(track_pt->at(k));
            int ntrks2 = track_pts1->size();
            if (!basicjet->at(k) || jetpt->at(k)<jptcut) continue;
            if (k!=j) p11 *= (1.0 - fakerate(jetpt->at(k),jeteta->at(k),ntrks2,varType));
        }
        p1 += p11;
    }
    return (1.0 - p0 - p1);
}

double QCDhists::frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p0 = 1.;
    double p1 = 0.;
    int njets = (njetscut == -1 ? jetpt->size() : std::min(njetscut,(int)jetpt->size()));
    if (verbose) {
        std::cout << "[frWeight1] varType = " << varType << std::endl;
        std::cout << "[frWeight1] njets = " << njets << std::endl;
        std::cout << "[frWeight1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        if (!basicjet->at(j) || jetpt->at(j)<jptcut) {
            if (verbose) std::cout << "[frWeight1] Jet[" << j << "] pT = " << jetpt->at(j) << " doesn't pass basic jet!" << std::endl;
            continue;
        }
        int ntrks1 = ntrack[j];
        if (verbose) std::cout << "[frWeight1] Jet pT = " << jetpt->at(j) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(j),jeteta->at(j),ntrks1,varType));
        double p11 = fakerate(jetpt->at(j),jeteta->at(j),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            int ntrks2 = ntrack[k];
            if (!basicjet->at(k) || jetpt->at(k)<jptcut) continue;
            if (k!=j) p11 *= (1.0 - fakerate(jetpt->at(k),jeteta->at(k),ntrks2,varType));
        }
        p1 += p11;
    }
    return (1.0 - p0 - p1);
}

double QCDhists::frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p0 = 1.;
    double p1 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeight1] varType = " << varType << std::endl;
        std::cout << "[frWeight1] njets = " << njets << std::endl;
        std::cout << "[frWeight1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        if (verbose) std::cout << "[frWeight1] Jet pT = " << jetpt->at(jdx) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType));
        double p11 = fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            if (kdx!=jdx) p11 *= (1.0 - fakerate(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType));
        }
        p1 += p11;
    }
    return (1.0 - p0 - p1);
}

double QCDhists::frWeightF1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType, std::string cutset){
    double p0 = 1.;
    double p1 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightF1] varType = " << varType << std::endl;
        std::cout << "[frWeightF1] njets = " << njets << std::endl;
        std::cout << "[frWeightF1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        int flav1 = flavor[jdx];
        if (verbose) std::cout << "[frWeightF1] Jet pT = " << jetpt->at(jdx) << std::endl;
        p0 *= (1.0 - fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1,cutset));
        double p11 = fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1,cutset);
        for(Int_t k=0; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            if (kdx!=jdx) p11 *= (1.0 - fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2,cutset));
        }
        p1 += p11;
    }
    return (1.0 - p0 - p1);
}

double QCDhists::frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p0 = 1.;
    int njets = (njetscut == -1 ? jetpt->size() : std::min(njetscut,(int)jetpt->size()));
    if (verbose) {
        std::cout << "[frWeightT0] varType = " << varType << std::endl;
        std::cout << "[frWeightT0] njets = " << njets << std::endl;
        std::cout << "[frWeightT0] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        if (!basicjet->at(j) || jetpt->at(j)<jptcut) {
            if (verbose) std::cout << "[frWeightT0] Jet[" << j << "] pT = " << jetpt->at(j) << " doesn't pass basic jet!" << std::endl;
            continue;
        }
        int ntrks = ntrack[j];
        if (verbose) std::cout << "[frWeightT0] Jet pT = " << jetpt->at(j) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(j),jeteta->at(j),ntrks,varType));
    }
    return p0;
}

double QCDhists::frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p0 = 1.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightT0] varType = " << varType << std::endl;
        std::cout << "[frWeightT0] njets = " << njets << std::endl;
        std::cout << "[frWeightT0] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks = ntrack[jdx];
        if (verbose) std::cout << "[frWeightT0] Jet pT = " << jetpt->at(jdx) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks,varType));
    }
    return p0;
}

double QCDhists::frWeightFT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType, std::string cutset){
    double p0 = 1.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightFT0] varType = " << varType << std::endl;
        std::cout << "[frWeightFT0] njets = " << njets << std::endl;
        std::cout << "[frWeightFT0] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks = ntrack[jdx];
        int flav = flavor[jdx];
        if (verbose) std::cout << "[frWeightFT0] Jet pT = " << jetpt->at(jdx) << std::endl;
        p0 *= (1.0 - fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks,varType,flav,cutset));
    }
    return p0;
}

double QCDhists::frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p1 = 0.;
    int njets = (njetscut == -1 ? jetpt->size() : std::min(njetscut,(int)jetpt->size()));
    if (verbose) {
        std::cout << "[frWeightT1] varType = " << varType << std::endl;
        std::cout << "[frWeightT1] njets = " << njets << std::endl;
        std::cout << "[frWeightT1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        if (!basicjet->at(j) || jetpt->at(j)<jptcut) {
            if (verbose) std::cout << "[frWeightT1] Jet[" << j << "] pT = " << jetpt->at(j) << " doesn't pass basic jet!" << std::endl;
            continue;
        }
        int ntrks1 = ntrack[j];
        if (verbose) std::cout << "[frWeightT1] Jet pT = " << jetpt->at(j) << std::endl;
        double p11 = fakerate(jetpt->at(j),jeteta->at(j),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            int ntrks2 = ntrack[k];
            if (!basicjet->at(k) || jetpt->at(k)<jptcut) continue;
            if (k!=j) p11 *= (1.0 - fakerate(jetpt->at(k),jeteta->at(k),ntrks2,varType));
        }
        p1 += p11;
    }
    return p1;
}

double QCDhists::frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p1 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightT1] varType = " << varType << std::endl;
        std::cout << "[frWeightT1] njets = " << njets << std::endl;
        std::cout << "[frWeightT1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        if (verbose) std::cout << "[frWeightT1] Jet pT = " << jetpt->at(jdx) << std::endl;
        double p11 = fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            if (kdx!=jdx) p11 *= (1.0 - fakerate(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType));
        }
        p1 += p11;
    }
    return p1;
}

double QCDhists::frWeightFT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType, std::string cutset){
    double p1 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightFT1] varType = " << varType << std::endl;
        std::cout << "[frWeightFT1] njets = " << njets << std::endl;
        std::cout << "[frWeightFT1] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        int flav1 = flavor[jdx];
        if (verbose) std::cout << "[frWeightFT1] Jet pT = " << jetpt->at(jdx) << std::endl;
        double p11 = fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1,cutset);
        for(Int_t k=0; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            if (kdx!=jdx) p11 *= (1.0 - fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2,cutset));
        }
        p1 += p11;
    }
    return p1;
}

//Ntag==1; unfolded b flavor dependent; only for 4-jet events
double QCDhists::frWeightFT12(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, double jptcut, int bUnfType, int nbTagged, int ptType, int varType, bool isData, std::string cutset){
    double frwgttmp=0.;
    if (verbose) {
        std::cout << "[frWeightFT12] varType = " << varType << std::endl;
        std::cout << "[frWeightFT12] jet pt cut = " << jptcut << std::endl;
    }
    // loop # unfolded "true" b jets
    TMatrixD Mwgt(5,5);
    bool goodWgt = UnfoldWgtPtDep(Mwgt, bUnfType, jetpt, isData);
    if (!goodWgt) {
        std::cout << "Bad unfolding matrix inversion!" << std::endl;
    }
    for (int nbT=0;nbT<5;nbT++){
        double evtWgt = 0.0;
//         if (isData) evtWgt=UnfoldWgtD(bUnfType, nbT, nbTagged, ptType);
//         else evtWgt=UnfoldWgt(bUnfType, nbT, nbTagged, ptType);
//         evtWgt=UnfoldWgtPtDep(bUnfType, nbT, nbTagged, jetpt, isData);
        evtWgt=Mwgt(nbT, nbTagged);

        double frwgts[5];
        int ncomb=1;
        switch (nbT) {
        case 1: ncomb=4; break;
        case 2: ncomb=6; break;
        case 3: ncomb=4; break;
        default: break;
        }
        evtWgt/=ncomb;
        switch (nbT) {
        case 0: {
            int flavors[] = {0,0,0,0};
            frWeightUFT1(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
            frwgttmp+=(evtWgt*frwgts[4]);
            break;}
        case 1: {
            for (int fidx=0;fidx<4;fidx++) {
                int flavors[] = {0,0,0,0};
                flavors[fidx]=5;
                frWeightUFT1(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                frwgttmp+=(evtWgt*frwgts[4]);
            }
            break;}
        case 2: {
            //flavor double counted due to looping
            for (int fidx0=0;fidx0<4;fidx0++) {
                for (int fidx1=fidx0+1;fidx1<4;fidx1++) {
                    int flavors[] = {0,0,0,0};
                    flavors[fidx0]=5;
                    flavors[fidx1]=5;
                    frWeightUFT1(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                    frwgttmp+=(evtWgt*frwgts[4]);
                }
            }
            break;}
        case 3: {
            for (int fidx=0;fidx<4;fidx++) {
                int flavors[] = {5,5,5,5};
                flavors[fidx]=0;
                frWeightUFT1(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                frwgttmp+=(evtWgt*frwgts[4]);
            }
            break;}
        case 4: {
            int flavors[] = {5,5,5,5};
            frWeightUFT1(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
            frwgttmp+=(evtWgt*frwgts[4]);
            break;}
        }
    }   
    return frwgttmp;
}

//Ntag==1; unfolded b flavor dependent; only for 4-jet events
void QCDhists::frWeightUFT1(double (&frwgts)[5], vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int* flavor, int varType, bool isData, std::string cutset){
    double frwgttmp=0.0;
    for(int i0=0;i0<4;i0++) {
        int idx0 = goodjetIdx[i0];
        double jfr = (isData ? fakerateFD(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset) :
                      fakerateF(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset));
        for(Int_t i1=0; i1<4; i1++) {
            int idx1 = goodjetIdx[i1];
            if (i0!=i1) jfr *= (isData ? (1.0 - fakerateFD(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset)):
                                (1.0 - fakerateF(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset)));
        }
        frwgts[i0]=jfr;
        frwgttmp+=jfr;
    }
    frwgts[4]=frwgttmp;
}


double QCDhists::frWeightT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p2 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightT2] varType = " << varType << std::endl;
        std::cout << "[frWeightT2] njets = " << njets << std::endl;
        std::cout << "[frWeightT2] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        double p21 = fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            if (k==j) continue;
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            double p22 = p21 * fakerate(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType);
            for(Int_t l=0; l<njets; l++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                p22 *= (1.0-fakerate(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType));
            }
            p2 += p22/2.0;
        }
    }
    return p2;
}

//Ntag==21; tag-and-probe
double QCDhists::frWeightT21(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p2 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightT21] varType = " << varType << std::endl;
        std::cout << "[frWeightT21] njets = " << njets << std::endl;
        std::cout << "[frWeightT21] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        double p11=fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType);
        for(Int_t k=0; k<njets; k++) {
            if (k==j) continue;
            double p12 = p11;
//             for(Int_t j1=0; j1<njets; j1++) {
//                 if (j1==j || j1==k) continue;
//                 int jdx1 = goodjetIdx[j1];
//                 int ntrks11 = ntrack[jdx1];
//                 p12 *= (1.0-fakerate(jetpt->at(jdx1),jeteta->at(jdx1),ntrks11,varType));
//             }

            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            double p22 = p12*fakerateTP(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType);
            for(Int_t l=0; l<njets; l++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                p22 *= (1.0-fakerateTP(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType));
            }
            p2 += p22/2.0;
        }
    }
    return p2;
}


//Ntag==2
double QCDhists::frWeightFT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType, std::string cutset){
    double p2 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightFT2] varType = " << varType << std::endl;
        std::cout << "[frWeightFT2] njets = " << njets << std::endl;
        std::cout << "[frWeightFT2] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        int flav1 = flavor[jdx];
        double p11=fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1,cutset);
        for(Int_t k=0; k<njets; k++) {
            if (k==j) continue;
            double p12 = p11;
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            double p22 = p12*fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2,cutset);
            int nComb = 0;
            for(Int_t l=0; l<njets; l++, nComb++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                int flav3 = flavor[ldx];
                p22 *= (1.0-fakerateF(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType,flav3,cutset));
            }
            //p2 += p22/nComb;
            p2 += p22/(njets-2);
        }
    }
    return p2;
}


//Ntag==2; unfolded b flavor dependent; only for 4-jet events
double QCDhists::frWeightFT22(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, double jptcut, int bUnfType, int nbTagged, int ptType, int varType, bool isData, std::string cutset){
    double frwgttmp=0.;
    if (verbose) {
        std::cout << "[frWeightFT22] varType = " << varType << std::endl;
        std::cout << "[frWeightFT22] jet pt cut = " << jptcut << std::endl;
    }

    TMatrixD Mwgt(5,5);
    bool goodWgt = UnfoldWgtPtDep(Mwgt, bUnfType, jetpt, isData);
    if (!goodWgt) {
        std::cout << "Bad unfolding matrix inversion!" << std::endl;
    }
    // loop # unfolded "true" b jets
    for (int nbT=0;nbT<5;nbT++){
        double evtWgt = 0.0;
//         if (isData) evtWgt=UnfoldWgtD(bUnfType, nbT, nbTagged, ptType);
//         else evtWgt=UnfoldWgt(bUnfType, nbT, nbTagged, ptType);
//         evtWgt=UnfoldWgtPtDep(bUnfType, nbT, nbTagged, jetpt, isData);
        evtWgt=Mwgt(nbT, nbTagged);

        double frwgts[7];
        int ncomb=1;
        switch (nbT) {
        case 1: ncomb=4; break;
        case 2: ncomb=6; break;
        case 3: ncomb=4; break;
        default: break;
        }
        evtWgt/=ncomb;
        switch (nbT) {
        case 0: {
            int flavors[] = {0,0,0,0};
            frWeightUFT23(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
            frwgttmp+=(evtWgt*frwgts[6]);
            break;}
        case 1: {
            for (int fidx=0;fidx<4;fidx++) {
                int flavors[] = {0,0,0,0};
                flavors[fidx]=5;
                frWeightUFT23(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                frwgttmp+=(evtWgt*frwgts[6]);
            }
            break;}
        case 2: {
            //flavor double counted due to looping
            for (int fidx0=0;fidx0<4;fidx0++) {
                for (int fidx1=fidx0+1;fidx1<4;fidx1++) {
                    int flavors[] = {0,0,0,0};
                    flavors[fidx0]=5;
                    flavors[fidx1]=5;
                    frWeightUFT23(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                    frwgttmp+=(evtWgt*frwgts[6]);
                }
            }
            break;}
        case 3: {
            for (int fidx=0;fidx<4;fidx++) {
                int flavors[] = {5,5,5,5};
                flavors[fidx]=0;
                frWeightUFT23(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
                frwgttmp+=(evtWgt*frwgts[6]);
            }
            break;}
        case 4: {
            int flavors[] = {5,5,5,5};
            frWeightUFT23(frwgts,jetpt,jeteta,goodjetIdx,ntrack,flavors,varType,isData,cutset);
            frwgttmp+=(evtWgt*frwgts[6]);
            break;}
        }
    }   
    return frwgttmp;
}


//Ntag==2; unfolded b flavor dependent; only for 4-jet events
double QCDhists::frWeightUFT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &jetIdx, vector<int> &ntrack, vector<int> &flavor, int varType, bool isData, std::string cutset){
    if (verbose) {
        for (int i=0;i<4;i++) {
            std::cout << "flavor [" << i << "]= " << flavor[i] << std::endl;
        }
    }
    double jfr=0.0;
    if (isData) {
        jfr = fakerateFD(jetpt->at(jetIdx[0]),jeteta->at(jetIdx[0]),ntrack[jetIdx[0]],varType,flavor[0],cutset);
        jfr *= (1.0-fakerateFD(jetpt->at(jetIdx[1]),jeteta->at(jetIdx[1]),ntrack[jetIdx[1]],varType,flavor[1],cutset));
        jfr *= fakerateFD(jetpt->at(jetIdx[2]),jeteta->at(jetIdx[2]),ntrack[jetIdx[2]],varType,flavor[2],cutset);
        jfr *= (1.0-fakerateFD(jetpt->at(jetIdx[3]),jeteta->at(jetIdx[3]),ntrack[jetIdx[3]],varType,flavor[3],cutset));
    }
    else {
        jfr = fakerateF(jetpt->at(jetIdx[0]),jeteta->at(jetIdx[0]),ntrack[jetIdx[0]],varType,flavor[0],cutset);
        jfr *= (1.0-fakerateF(jetpt->at(jetIdx[1]),jeteta->at(jetIdx[1]),ntrack[jetIdx[1]],varType,flavor[1],cutset));
        jfr *= fakerateF(jetpt->at(jetIdx[2]),jeteta->at(jetIdx[2]),ntrack[jetIdx[2]],varType,flavor[2],cutset);
        jfr *= (1.0-fakerateF(jetpt->at(jetIdx[3]),jeteta->at(jetIdx[3]),ntrack[jetIdx[3]],varType,flavor[3],cutset));
    }
    return jfr;
}

//Ntag==2; unfolded b flavor dependent; only for 4-jet events
void QCDhists::frWeightUFT22(double (&frwgts)[25], vector<float> *jetpt, vector<float> *jeteta,
                   vector<int> &goodjetIdx, vector<int> &ntrack, int* flavor, int varType, bool isData, std::string cutset){
    double frwgttmp=0.0;
    int nthComb=0;
    for(int i0=0;i0<4;i0++) {
        int idx0 = goodjetIdx[i0];
        double jfr = (isData? fakerateFD(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset):
                      fakerateF(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset));

        for(Int_t i1=0; i1<4; i1++) {
            if (i1==i0) continue;
            int idx1 = goodjetIdx[i1];
            double kfr = jfr;
            if (isData) kfr*=(1.0 - fakerateFD(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset));
            else kfr*=(1.0 - fakerateF(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset));

            for(Int_t i2=0; i2<4; i2++) {
                if (i2==i0 || i2==i1) continue;
                int idx2 = goodjetIdx[i2];
                double lfr = kfr;
                if (isData) lfr*=fakerateFD(jetpt->at(idx2),jeteta->at(idx2),ntrack[idx2],varType,flavor[i2],cutset);
                else lfr*=fakerateF(jetpt->at(idx2),jeteta->at(idx2),ntrack[idx2],varType,flavor[i2],cutset);

                for(Int_t i3=0; i3<4; i3++) {
                    if (i3==i0 || i3==i1 || i3==i2) continue;
                    int idx3 = goodjetIdx[i3];
                    if (isData) lfr *= (1.0 - fakerateFD(jetpt->at(idx3),jeteta->at(idx3),ntrack[idx3],varType,flavor[i3],cutset));
                    else lfr *= (1.0 - fakerateF(jetpt->at(idx3),jeteta->at(idx3),ntrack[idx3],varType,flavor[i3],cutset));
                }
                frwgts[nthComb]=lfr;
                frwgttmp+=lfr/2;
                nthComb++;
            }
        }
    }
    //std::cout<<"nTotalComb = "<<nthComb<<std::endl;
    frwgts[24]=frwgttmp;
}

//Ntag==2; unfolded b flavor dependent; only for 4-jet events
void QCDhists::frWeightUFT23(double (&frwgts)[7], vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int* flavor, int varType, bool isData, std::string cutset){
    double frwgttmp=0.0;
    int nthComb=0;
    for(int i0=0;i0<4;i0++) {
        int comb[] = {0,0,0,0};
        comb[0]=i0;
        int idx0 = goodjetIdx[i0];
        double jfr = (isData ? fakerateFD(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset):
                      fakerateF(jetpt->at(idx0),jeteta->at(idx0),ntrack[idx0],varType,flavor[i0],cutset));

        for(Int_t i1=i0+1; i1<4; i1++) {
            comb[1]=i1;
            int idx1 = goodjetIdx[i1];
            double kfr = jfr;
            if (isData) kfr*=fakerateFD(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset);
            else kfr*=fakerateF(jetpt->at(idx1),jeteta->at(idx1),ntrack[idx1],varType,flavor[i1],cutset);
            int nTemp=0;
            for(Int_t i2=0; i2<4; i2++) {
                if (i2==i0 || i2==i1) continue;
                comb[nTemp+2]=i2;
                int idx2 = goodjetIdx[i2];
                if (isData) kfr*=(1.0-fakerateFD(jetpt->at(idx2),jeteta->at(idx2),ntrack[idx2],varType,flavor[i2],cutset));
                else kfr*=(1.0-fakerateF(jetpt->at(idx2),jeteta->at(idx2),ntrack[idx2],varType,flavor[i2],cutset));
                nTemp++;
            }
            frwgts[nthComb]=kfr;
            frwgttmp+=kfr;
            //std::cout<<"["<<nthComb<<"] Comb = ["<<comb[0]+1<<","<<comb[1]+1<<","<<comb[2]+1<<","<<comb[3]+1<<"]" <<std::endl;
            nthComb++;
        }
    }
    //std::cout<<"nTotalComb = "<<nthComb<<std::endl;
    frwgts[6]=frwgttmp;
}


double QCDhists::frWeightT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
    double p3 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightT3] varType = " << varType << std::endl;
        std::cout << "[frWeightT3] njets = " << njets << std::endl;
        std::cout << "[frWeightT3] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        double p31 = fakerate(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType);
        for(Int_t k=j+1; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            double p32 = p31 * fakerateTP(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType);
            for(Int_t l=k+1; l<njets; l++) {
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                double p33 = p32 * fakerateTP(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType);
                for(Int_t m=0; m<njets; m++) {
                    if (m==j || m==k || m ==l) continue;
                    int mdx = goodjetIdx[m];
                    int ntrks4 = ntrack[mdx];
                    p33 *= (1.0-fakerateTP(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType));
                    //p33 *= (1.0-fakerate(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType));
                }
                p3 += p33;
            }
        }
    }
    return p3;
}

double QCDhists::frWeightFT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType, std::string cutset){
    double p3 = 0.;
    int njets = (njetscut == -1 ? goodjetIdx.size() : std::min(njetscut,(int)goodjetIdx.size()));
    if (verbose) {
        std::cout << "[frWeightFT3] varType = " << varType << std::endl;
        std::cout << "[frWeightFT3] njets = " << njets << std::endl;
        std::cout << "[frWeightFT3] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<njets; j++) {
        int jdx = goodjetIdx[j];
        int ntrks1 = ntrack[jdx];
        int flav1 = flavor[jdx];
        double p31 = fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1,cutset);
        for(Int_t k=0; k<njets; k++) {
            if (k==j) continue;
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            double p32 = p31 * fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2,cutset);
            for(Int_t l=0; l<njets; l++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                int flav3 = flavor[ldx];
                double p33 = p32 * fakerateF(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType,flav3,cutset);
                for(Int_t m=0; m<njets; m++) {
                    if (m==j || m==k || m ==l) continue;
                    int mdx = goodjetIdx[m];
                    int ntrks4 = ntrack[mdx];
                    int flav4 = flavor[mdx];
                    p33 *= (1.0-fakerateF(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType,flav4,cutset));
                }
                p3 += p33/6.0;
            }
        }
    }
    return p3;
}

// Only 4 leading jets
double QCDhists::frWeight4(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, double jptcut, int varType){
    double p0 = 1.;
    double p1 = 0.;
    double p4 = 1.;
    if (verbose) {
        std::cout << "[frWeight4] varType = " << varType << std::endl;
        std::cout << "[frWeight4] jet pt cut = " << jptcut << std::endl;
    }
    for(Int_t j=0; j<4; j++) {
        if (!basicjet->at(j) || jetpt->at(j)<jptcut) {
            if (verbose) std::cout << "[frWeight4] Jet[" << j << "] pT = " << jetpt->at(j) << " doesn't pass basic jet!" << std::endl;
            continue;
        }
        int ntrks = ntrack[j];
        if (verbose) std::cout << "[frWeight4] Jet pT = " << jetpt->at(j) << std::endl;
        p0 *= (1.0 - fakerate(jetpt->at(j),jeteta->at(j),ntrks,varType));
        p4 *= fakerate(jetpt->at(j),jeteta->at(j),ntrks,varType);
        double p11 = 1.;
        for(Int_t k=0; k<4; k++) {
            int ntrks1 = ntrack[k];
            if (!basicjet->at(k) || jetpt->at(k)<jptcut) continue;
            if (k!=j) p11 *= (1.0 - fakerate(jetpt->at(k),jeteta->at(k),ntrks1,varType));
        }
        p1 += (fakerate(jetpt->at(j),jeteta->at(j),ntrks,varType)*p11);
    }
    double fr1 = 1.0 - p0 - p1;
    double fr2 = 1.0 - p0 - p1 - p4;
    if (verbose) {
        std::cout << std::fixed << std::setprecision(9)
                  << "fr1 = " << fr1 << std::endl
                  << "fr2 = " << fr2 << std::endl
                  << std::resetiosflags(std::ios::fixed);
    }
    return (1.0 - p0 - p1 - p4);
}


double QCDhists::GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        ptsum_total += track_pt.at(itk);
        if ( track_pvWeight.at(itk) > 0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha] alpha = " << alpha << std::endl;
    return alpha;
}

double QCDhists::GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
                vector<float> &track_pvWeight, vector<float> &track_ref_zs, float pv_z, float pilecut)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        if (fabs(pv_z-track_ref_zs.at(itk))>pilecut) continue;// remove tracks with exceedingly large z
        ptsum_total += track_pt.at(itk);
        if ( track_pvWeight.at(itk) > 0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha] alpha = " << alpha << std::endl;
    return alpha;
}

double QCDhists::GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_ipXYSigs)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        ptsum_total += track_pt.at(itk);
        if ( fabs(track_ipXYSigs.at(itk)) < 4.0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha2Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

// default for analysis_20170523_v0
double QCDhists::GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
                     vector<float> &track_ipXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        if (fabs(pv_z-track_ref_zs.at(itk))>pilecut) continue;// remove tracks with exceedingly large z
        ptsum_total += track_pt.at(itk);
        if ( fabs(track_ipXYSigs.at(itk)) < 4.0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha2Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

double QCDhists::GetAlpha3Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
                               vector<float> &track_ipXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut, float sigzcut)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        if (fabs(pv_z-track_ref_zs.at(itk))>pilecut) continue;// remove tracks with exceedingly large z
        ptsum_total += track_pt.at(itk);
        double ip3Dsimp = sqrt(pow((pv_z-track_ref_zs.at(itk))/0.01,2)+pow(track_ipXYSigs.at(itk),2));
        if ( ip3Dsimp < sigzcut ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha3Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

double QCDhists::nGJrewgt(int nGoodJet){
    double wgt = 1.0;

    if (nGoodJet==4)       wgt = 0.99161;
    else if (nGoodJet==5)  wgt = 1.35307;
    else if (nGoodJet==6)  wgt = 1.49794;
    else if (nGoodJet==7)  wgt = 1.74202;
    else if (nGoodJet==8)  wgt = 1.68261;
    else if (nGoodJet==9)  wgt = 1.20239;
    else if (nGoodJet==10) wgt = 1.60758;

    return wgt;
}

//~/Dropbox/UMD\ Analysis/FakeRate/20170926/responsematrixandinvert_MC_lowpT.txt
//~/Dropbox/UMD\ Analysis/FakeRate/20170926/responsematrixandinvert_MC_highpT.txt
//Fixed with correct MC mis-tag rate
double QCDhists::UnfoldWgt(int bUnfType, int nbtrue, int nbtagged, int ptType){
    // ptType = 1: low pT mistag; 2: high pT mistag
    if (ptType==1) {
        if (bUnfType==1) {//CSVv2L
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.4772270061;
                case 1: return -0.1995346353;
                case 2: return 0.0269518974;
                case 3: return -0.0036404947;
                case 4: return 0.0004917354;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.5491470731;
                case 1: return 1.5951064013;
                case 2: return -0.4208949539;
                case 3: return 0.0846011464;
                case 4: return -0.0151755960;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0765528521;
                case 1: return -0.4343855485;
                case 2: return 1.6818795409;
                case 3: return -0.6579461414;
                case 4: return 0.1756270180;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0047429829;
                case 1: return 0.0400494935;
                case 2: return -0.3017932177;
                case 3: return 1.7323689470;
                case 4: return -0.9033465320;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0001101978;
                case 1: return -0.0012357111;
                case 2: return 0.0138567333;
                case 3: return -0.1553834574;
                case 4: return 1.7424033746;
                default: break;
                }
            }
        }
        else if (bUnfType==2){//CSVv2M
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0489790138;
                case 1: return -0.4119945430;
                case 2: return 0.1618140127;
                case 3: return -0.0635536930;
                case 4: return 0.0249611998;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0498607489;
                case 1: return 1.4631957629;
                case 2: return -1.1416712401;
                case 3: return 0.6710897498;
                case 4: return -0.3510386959;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0008887550;
                case 1: return -0.0518131102;
                case 2: return 2.0272730472;
                case 3: return -2.3647570517;
                case 4: return 1.8512957186;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0000070408;
                case 1: return 0.0006143209;
                case 2: return -0.0476982683;
                case 3: return 2.7900427223;
                case 4: return -4.3392441127;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.000000020917;
                case 1: return -0.0000024306;
                case 2: return 0.0002824485;
                case 3: return -0.0328217274;
                case 4: return 3.8140258903;
                default: break;
                }
            }

        }
        else if (bUnfType==3){//CSVv2T
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0095747547;
                case 1: return -1.0382442203;
                case 2: return 1.0677278290;
                case 3: return -1.0980487004;
                case 4: return 1.1292306108;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0096089967;
                case 1: return 2.0528281275;
                case 2: return -4.2120844092;
                case 3: return 6.4923204909;
                case 4: return -8.8986661650;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0000342964;
                case 1: return -0.0146186361;
                case 2: return 4.1640781981;
                case 3: return -12.8006388748;
                case 4: return 26.2965306067;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.000000054405;
                case 1: return 0.0000347566;
                case 2: return -0.0197450810;
                case 3: return 8.4263449007;
                case 4: return -34.5373869211;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 3.23636E-11;
                case 1: return -2.75563E-08;
                case 2: return 0.0000234630;
                case 3: return -0.0199778164;
                case 4: return 17.0102918686;
                default: break;
                }
            }
        }
        else{
            std::cout << "b unfolding type type not defined!" << std::endl;
        }
    }
    else if (ptType==2) {//High pT mistag
        if (bUnfType==1) {//CSVv2L
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.6686791061;
                case 1: return -0.2253947941;
                case 2: return 0.0304449268;
                case 3: return -0.0041123113;
                case 4: return 0.0005554654;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.8019921327;
                case 1: return 1.7748219486;
                case 2: return -0.4648315310;
                case 3: return 0.0931916086;
                case 4: return -0.0166946688;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.1445435296;
                case 1: return -0.6202304367;
                case 2: return 1.8292279883;
                case 3: return -0.7076502549;
                case 4: return 0.1881610991;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0115782970;
                case 1: return 0.0737410436;
                case 2: return -0.4196562002;
                case 3: return 1.8281778706;
                case 4: return -0.9425390614;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0003477939;
                case 1: return -0.0029377614;
                case 2: return 0.0248148162;
                case 3: return -0.2096069130;
                case 4: return 1.7705171657;
                default: break;
                }
            }
        }
        else if (bUnfType==2){//CSVv2M
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0741125093;
                case 1: return -0.4218659159;
                case 2: return 0.1656910700;
                case 3: return -0.0650764370;
                case 4: return 0.0255592691;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0761111830;
                case 1: return 1.4993705671;
                case 2: return -1.1660377483;
                case 3: return 0.6846497414;
                case 4: return -0.3579315341;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0020224530;
                case 1: return -0.0788892654;
                case 2: return 2.0720267640;
                case 3: return -2.4050276457;
                case 4: return 1.8796749795;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0000238850;
                case 1: return 0.0013928249;
                case 2: return -0.0723174088;
                case 3: return 2.8349237925;
                case 4: return -4.3871556331;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.000000105780;
                case 1: return -0.0000082107;
                case 2: return 0.0006373230;
                case 3: return -0.0494694512;
                case 4: return 3.8398529186;
                default: break;
                }
            }

        }
        else if (bUnfType==3){//CSVv2T
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0130127955;
                case 1: return -1.0417798931;
                case 2: return 1.0713639063;
                case 3: return -1.1017880334;
                case 4: return 1.1330761318;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0130759535;
                case 1: return 2.0616091592;
                case 2: return -4.2264785337;
                case 3: return 6.5126394234;
                case 4: return -8.9252347979;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0000632941;
                case 1: return -0.0198933228;
                case 2: return 4.1819296106;
                case 3: return -12.8390099106;
                case 4: return 26.3640105347;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -1.36166E-07;
                case 1: return 6.41257E-05;
                case 2: return -0.0268582480;
                case 3: return 8.4553102939;
                case 4: return -34.6115289818;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 1.09852E-10;
                case 1: return -6.89401E-08;
                case 2: return 0.0000432648;
                case 3: return -0.0271517732;
                case 4: return 17.0396771132;
                default: break;
                }
            }
        }
        else{
            std::cout << "b unfolding type type not defined!" << std::endl;
        }
    }
    else{
        std::cout << "b mis-tag rate pT type type not defined!" << std::endl;
    }

    return 0;
}

// Unfold for data
double QCDhists::UnfoldWgtD(int bUnfType, int nbtrue, int nbtagged, int ptType){
    // ptType = 1: low pT mistag; 2: high pT mistag
    if (ptType==1) {
        if (bUnfType==1) {//CSVv2L
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.5578131183;
                case 1: return -0.2599384783;
                case 2: return 0.0433736317;
                case 3: return -0.0072373738;
                case 4: return 0.0012076365;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.6536634742;
                case 1: return 1.7361390042;
                case 2: return -0.5611885615;
                case 3: return 0.1389424348;
                case 4: return -0.0307432153;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.1028547485;
                case 1: return -0.5292047379;
                case 2: return 1.8731453134;
                case 3: return -0.8939399822;
                case 4: return 0.2934902128;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0071930314;
                case 1: return 0.0549138905;
                case 2: return -0.3746629576;
                case 3: return 1.9579476642;
                case 4: return -1.2452468188;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0001886387;
                case 1: return -0.0019096785;
                case 2: return 0.0193325740;
                case 3: return -0.1957127430;
                case 4: return 1.9812921849;
                default: break;
                }
            }
        }
        else if (bUnfType==2){//CSVv2M
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0575929199;
                case 1: return -0.5232634626;
                case 2: return 0.2588941796;
                case 3: return -0.1280926359;
                case 4: return 0.0633761771;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0588078574;
                case 1: return 1.5879766164;
                case 2: return -1.5569651070;
                case 3: return 1.1519441743;
                case 4: return -0.7587526123;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0012262625;
                case 1: return -0.0656183180;
                case 2: return 2.3624040681;
                case 3: return -3.4584872903;
                case 4: return 3.4064783068;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0000113645;
                case 1: return 0.0009093715;
                case 2: return -0.0647813329;
                case 3: return 3.4823802686;
                case 4: return -6.7971763515;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0000000395;
                case 1: return -0.0000042073;
                case 2: return 0.0004481922;
                case 3: return -0.0477445167;
                case 4: return 5.0860744799;
                default: break;
                }
            }

        }
//         else if (bUnfType==3){//CSVv2T
//             if (nbtrue==0){
//                 switch (nbtagged) {
//                 case 0: return 1.0095747547;
//                 case 1: return -1.0382442203;
//                 case 2: return 1.0677278290;
//                 case 3: return -1.0980487004;
//                 case 4: return 1.1292306108;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==1){
//                 switch (nbtagged) {
//                 case 0: return -0.0096089967;
//                 case 1: return 2.0528281275;
//                 case 2: return -4.2120844092;
//                 case 3: return 6.4923204909;
//                 case 4: return -8.8986661650;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==2){
//                 switch (nbtagged) {
//                 case 0: return 0.0000342964;
//                 case 1: return -0.0146186361;
//                 case 2: return 4.1640781981;
//                 case 3: return -12.8006388748;
//                 case 4: return 26.2965306067;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==3){
//                 switch (nbtagged) {
//                 case 0: return -0.000000054405;
//                 case 1: return 0.0000347566;
//                 case 2: return -0.0197450810;
//                 case 3: return 8.4263449007;
//                 case 4: return -34.5373869211;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==4){
//                 switch (nbtagged) {
//                 case 0: return 3.23636E-11;
//                 case 1: return -2.75563E-08;
//                 case 2: return 0.0000234630;
//                 case 3: return -0.0199778164;
//                 case 4: return 17.0102918686;
//                 default: break;
//                 }
//             }
//         }
        else{
            std::cout << "b unfolding type type not defined!" << std::endl;
        }
    }
    else if (ptType==2) {//High pT mistag
        if (bUnfType==1) {//CSVv2L
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.7862707518;
                case 1: return -0.2980591803;
                case 2: return 0.0497344957;
                case 3: return -0.0082987548;
                case 4: return 0.0013847397;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.9646278926;
                case 1: return 1.9638921439;
                case 2: return -0.6285367895;
                case 3: return 0.1550767905;
                case 4: return -0.0342524569;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.1953455902;
                case 1: return -0.7628151036;
                case 2: return 2.0688806023;
                case 3: return -0.9728385349;
                case 4: return 0.3177214785;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0175818629;
                case 1: return 0.1015176246;
                case 2: return -0.5247432097;
                case 3: return 2.0910058076;
                case 4: return -1.3098419741;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0005934135;
                case 1: return -0.0045354846;
                case 2: return 0.0346649013;
                case 3: return -0.2649453084;
                case 4: return 2.0249882127;
                default: break;
                }
            }
        }
        else if (bUnfType==2){//CSVv2M
            if (nbtrue==0){
                switch (nbtagged) {
                case 0: return 1.0916058829;
                case 1: return -0.5400919989;
                case 2: return 0.2672204060;
                case 3: return -0.1322121889;
                case 4: return 0.0654144014;
                default: break;
                }
            }
            else if (nbtrue==1){
                switch (nbtagged) {
                case 0: return -0.0946384771;
                case 1: return 1.6431563521;
                case 2: return -1.6027962953;
                case 3: return 1.1837880525;
                case 4: return -0.7790440191;
                default: break;
                }
            }
            else if (nbtrue==2){
                switch (nbtagged) {
                case 0: return 0.0030768115;
                case 1: return -0.1053196961;
                case 2: return 2.4378936595;
                case 3: return -3.5416048185;
                case 4: return 3.4792200036;
                default: break;
                }
            }
            else if (nbtrue==3){
                switch (nbtagged) {
                case 0: return -0.0000444582;
                case 1: return 0.0022717156;
                case 2: return -0.1034305536;
                case 3: return 3.5656598577;
                case 4: return -6.9058838138;
                default: break;
                }
            }
            else if (nbtrue==4){
                switch (nbtagged) {
                case 0: return 0.0000002409;
                case 1: return -0.0000163728;
                case 2: return 0.0011127835;
                case 3: return -0.0756309028;
                case 4: return 5.1402934279;
                default: break;
                }
            }

        }
//         else if (bUnfType==3){//CSVv2T
//             if (nbtrue==0){
//                 switch (nbtagged) {
//                 case 0: return 1.0130127955;
//                 case 1: return -1.0417798931;
//                 case 2: return 1.0713639063;
//                 case 3: return -1.1017880334;
//                 case 4: return 1.1330761318;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==1){
//                 switch (nbtagged) {
//                 case 0: return -0.0130759535;
//                 case 1: return 2.0616091592;
//                 case 2: return -4.2264785337;
//                 case 3: return 6.5126394234;
//                 case 4: return -8.9252347979;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==2){
//                 switch (nbtagged) {
//                 case 0: return 0.0000632941;
//                 case 1: return -0.0198933228;
//                 case 2: return 4.1819296106;
//                 case 3: return -12.8390099106;
//                 case 4: return 26.3640105347;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==3){
//                 switch (nbtagged) {
//                 case 0: return -1.36166E-07;
//                 case 1: return 6.41257E-05;
//                 case 2: return -0.0268582480;
//                 case 3: return 8.4553102939;
//                 case 4: return -34.6115289818;
//                 default: break;
//                 }
//             }
//             else if (nbtrue==4){
//                 switch (nbtagged) {
//                 case 0: return 1.09852E-10;
//                 case 1: return -6.89401E-08;
//                 case 2: return 0.0000432648;
//                 case 3: return -0.0271517732;
//                 case 4: return 17.0396771132;
//                 default: break;
//                 }
//             }
//         }
        else{
            std::cout << "b unfolding type type not defined!" << std::endl;
        }
    }
    else{
        std::cout << "b mis-tag rate pT type type not defined!" << std::endl;
    }

    return 0;
}


// Ref: AN-17-018 (table not shown in the note but in text:
// /Users/jengbou/Dropbox/MyResearch/MyAnalysis/CMS/Notes/tdr2/notes/AN-17-018/trunk/MCefficiencyQCD.tex
double QCDhists::effCSVv2(int bUnfType, double pt, int flav, bool isData) {
    double eff_ = 0;
    double sf_ = 1.0;
    if (pt>=1000) pt = 1000.;
    if (isData) {
        if (flav==0) {//b and g->bb
            sf_ =  readerB_.eval(BTagEntryStandalone::FLAV_B, 0.0, pt);
            //std::cout<< "scalefactor b jet = " << sf_ << std::endl;
        }
        else if (flav==1) {//udsg
            sf_ =  readerL_.eval(BTagEntryStandalone::FLAV_UDSG, 0.0, pt);
            //std::cout<< "scalefactor l jet = " << sf_ << std::endl;
        }
        else{//c
        }
    }
    if (bUnfType==1) {
        //CSVv2L
        if (flav==0) {//b and g->bb
            if (pt>=30 && pt<165){
                eff_ = 0.625404+0.00629137*pt-7.44931*pow(10,-5)*pow(pt,2)+3.79928*pow(10,-7)*pow(pt,3)-7.23654*pow(10,-10)*pow(pt,4);
            }
            else if (pt>=165){
                eff_ = 0.845173-0.000226895*pt-2.12743*pow(10,-8)*pow(pt,2);
            }
        }
        else if (flav==1) {//udsg
            if (pt>=30 && pt<195){
                eff_ = 0.239697-0.0060077*pt+7.44621*pow(10,-5)*pow(pt,2)-3.7838*pow(10,-7)*pow(pt,3)+6.96297*pow(10,-10)*pow(pt,4);
            }
            else if (pt>=195){
                eff_ = 0.0652218+0.000191387*pt-5.82485*pow(10,-8)*pow(pt,2);
            }
        }
    }
    else if (bUnfType==2) {//CSVv2M
    }
    else if (bUnfType==3) {//CSVv2T
    }
    else {//Undefined
    }

    return eff_*sf_;
}

// Ref: BTV-16-002 p.96
double QCDhists::effDeepCSV(int bUnfType, double pt, int flav, bool isData) {
    double eff_ = 0;
    if (bUnfType==1) {
        //DeepCSVL
        if (isData) {
            if (flav==0) {//b and g->bb
                if (pt>=20 && pt<160){
                    eff_=0.4344+0.02069*pt-0.0004429*pow(pt,2)+5.137*pow(10,-6)*pow(pt,3)-3.406*pow(10,-8)*pow(pt,4)+1.285*pow(10,-10)*pow(pt,5)-2.559*pow(10,-13)*pow(pt,6)+2.084*pow(10,-16)*pow(pt,7);
                }
                else if (pt>=160 && pt<300){
                    eff_=0.714+0.002617*pt-1.656*pow(10,-5)*pow(pt,2)+4.767*pow(10,-8)*pow(pt,3)-6.431*pow(10,-11)*pow(pt,4)+3.287*pow(10,-14)*pow(pt,5);
                }
                else if (pt>=300){// removed upper limit (1000) to avoid sigular matrix due to zero efficiency
                    eff_=0.872-6.885*pow(10,-5)*pt+4.34*pow(10,-8)*pow(pt,2);
                }
            }
            else if (flav==1) {//udsg
                if (pt>=20 && pt<150){
                    eff_=0.245-0.0054*pt+6.92*pow(10,-5)*pow(pt,2)-3.89*pow(10,-7)*pow(pt,3)+1.021*pow(10,-9)*pow(pt,4)-1.007*pow(10,-12)*pow(pt,5);
                }
                else if (pt>=150){// removed upper limit (1000) to avoid sigular matrix due to zero efficiency
                    eff_=0.0558+0.000428*pt-1.0*pow(10,-7)*pow(pt,2);
                }
            }
            else{//c
            }
        }
        else {//ttbar MC
            if (flav==0) {//b and g->bb
                if (pt>=20 && pt<100){
                    eff_=0.491+0.0191*pt-0.0004172*pow(pt,2)+4.893*pow(10,-6)*pow(pt,3)-3.266*pow(10,-8)*pow(pt,4)+1.238*pow(10,-10)*pow(pt,5)-2.474*pow(10,-13)*pow(pt,6)+2.021*pow(10,-16)*pow(pt,7);
                }
                else if (pt>=100 && pt<300){
                    eff_=0.912-0.001846*pt+2.479*pow(10,-5)*pow(pt,2)-1.417*pow(10,-7)*pow(pt,3)+3.617*pow(10,-10)*pow(pt,4)-3.433*pow(10,-13)*pow(pt,5);
                }
                else if (pt>=300){// removed upper limit (1000) to avoid sigular matrix due to zero efficiency
                    eff_=0.892-0.00014*pt+1.01*pow(10,-7)*pow(pt,2);
                }
            }
            else if (flav==1) {//udsg
                if (pt>=20 && pt<250){
                    eff_=0.2407-0.00593*pt+8.5*pow(10,-5)*pow(pt,2)-5.658*pow(10,-7)*pow(pt,3)+1.828*pow(10,-9)*pow(pt,4)-2.287*pow(10,-12)*pow(pt,5);
                }
                else if (pt>=250){// removed upper limit (1000) to avoid sigular matrix due to zero efficiency
                    eff_=0.0541+0.00036*pt-7.392*pow(10,-8)*pow(pt,2);
                }
            }
            else{//c
            }
        }
    }
    else if (bUnfType==2) {//DeepCSVM
        if (isData) {}
        else {}
    }
    else if (bUnfType==3) {//DeepCSVT
        if (isData) {}
        else {}
    }
    else {//Undefined
    }

    return eff_;
}


double QCDhists::effBTag(int bUnfType, double pt, int flav, bool isData) {
    double eff_ = 0;
    double sf_ = 1.0;
    if (isData) {
        if (pt>=1000) pt = 999.;
        if (flav==0) {//b and g->bb
            sf_ =  readerB_.eval(BTagEntryStandalone::FLAV_B, 0.0, pt);
            //std::cout<< "scalefactor b jet = " << sf_ << std::endl;
        }
        else if (flav==1) {//udsg
            sf_ =  readerL_.eval(BTagEntryStandalone::FLAV_UDSG, 0.0, pt);
            //std::cout<< "scalefactor l jet = " << sf_ << std::endl;
        }
        else{//c
        }
    }
    if (bUnfType==1) {
        //CSVv2L
        if (flav==0) {//b and g->bb
            // systematic 500toInf all selected jets
//             if (pt>=20 && pt<30){
//                 eff_ = 0.0;
//             }
//             else if (pt>=30 && pt<30){
//                 eff_ = 0.801967;
//             }
//             else if (pt>=50 && pt<70){
//                 eff_ = 0.799573;
//             }
//             else if (pt>=70 && pt<100){
//                 eff_ = 0.80511;
//             }
//             else if (pt>=100 && pt<140){
//                 eff_ = 0.804147;
//             }
//             else if (pt>=140 && pt<200){
//                 eff_ = 0.809298;
//             }
//             else if (pt>=200 && pt<300){
//                 eff_ = 0.809617;
//             }
//             else if (pt>=300 && pt<600){
//                 eff_ = 0.816883;
//             }
//             else if (pt>=600){
//                 eff_ = 0.803702;
//             }

            // Default 500toInf passing pre-selection
            if (pt>=100 && pt<140){
                eff_ = 0.807549;//0.802547;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.804981;//0.801777;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.803957;//0.802073;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.8076;//0.811238;
            }
            else if (pt>=600){
                eff_ = 0.78868;//0.803135;
            }

        }
        else if (flav==1) {//udsgc
            // systematic 500toInf all selected jets
//             if (pt>=20 && pt<30){
//                 eff_ = 0.0;
//             }
//             else if (pt>=30 && pt<30){
//                 eff_ = 0.122245;
//             }
//             else if (pt>=50 && pt<70){
//                 eff_ = 0.124183;
//             }
//             else if (pt>=70 && pt<100){
//                 eff_ = 0.128752;
//             }
//             else if (pt>=100 && pt<140){
//                 eff_ = 0.133395;
//             }
//             else if (pt>=140 && pt<200){
//                 eff_ = 0.145496;
//             }
//             else if (pt>=200 && pt<300){
//                 eff_ = 0.156708;
//             }
//             else if (pt>=300 && pt<600){
//                 eff_ = 0.179475;
//             }
//             else if (pt>=600){
//                 eff_ = 0.222216;
//             }

            // 500toInf passing pre-selection
            if (pt>=100 && pt<140){
                eff_ = 0.127032;//0.131396;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.138622;//0.142532;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.151775;//0.154968;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.182386;//0.18302;
            }
            else if (pt>=600){
                eff_ = 0.22567;//0.223272;
            }

        }
        else{//c
        }
    }//End CSVv2L
    else if (bUnfType==2) {//CSVv2M
        if (flav==0) {//b and g->bb
            if (pt>=100 && pt<140){
                eff_ = 0.640644;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.629018;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.611541;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.58964;
            }
            else if (pt>=600){
                eff_ = 0.519469;
            }
        }
        else if (flav==1) {//udsgc
            if (pt>=100 && pt<140){
                eff_ = 0.0226542;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.0264172;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.0305033;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.038221;
            }
            else if (pt>=600){
                eff_ = 0.045846;
            }
        }
    }//End CSVv2M
    else if (bUnfType==3) {//CSVv2T
        if (flav==0) {//b and g->bb
            if (pt>=100 && pt<140){
                eff_ = 0.442714;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.427096;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.401281;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.36345;
            }
            else if (pt>=600){
                eff_ = 0.318183;
            }
        }
        else if (flav==1) {//udsgc
            if (pt>=100 && pt<140){
                eff_ = 0.00395776;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.00502544;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.00597841;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.00768997;
            }
            else if (pt>=600){
                eff_ = 0.0123835;
            }
        }

    }//End CSVv2T
    else if (bUnfType==4) {//CSVv2 User
        if (flav==0) {//b and g->bb
            if (pt>=100 && pt<140){
                eff_ = 0.77394;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.7736;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.77264;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.769626;
            }
            else if (pt>=600){
                eff_ = 0.729124;
            }
        }
        else if (flav==1) {//udsgc
            if (pt>=100 && pt<140){
                eff_ = 0.0880548;
            }
            else if (pt>=140 && pt<200){
                eff_ = 0.101232;
            }
            else if (pt>=200 && pt<300){
                eff_ = 0.115253;
            }
            else if (pt>=300 && pt<600){
                eff_ = 0.13782;
            }
            else if (pt>=600){
                eff_ = 0.159787;
            }
        }
    }//End CSVv2 user
    else {//Undefined
    }

    return eff_*sf_;
}


// UMD Emerging Jets group estimation
double QCDhists::effBTagPara(int bUnfType, double pt, int flav, bool isData) {
    double eff_ = 0;
    double sf_ = 1.0;
    if (pt>=1000) pt = 1000.;
    if (isData) {
        if (flav==0) {//b and g->bb
            sf_ =  readerB_.eval(BTagEntryStandalone::FLAV_B, 0.0, pt);
            //std::cout<< "scalefactor b jet = " << sf_ << std::endl;
        }
        else if (flav==1) {//udsg
            sf_ =  readerL_.eval(BTagEntryStandalone::FLAV_UDSG, 0.0, pt);
            //std::cout<< "scalefactor l jet = " << sf_ << std::endl;
        }
        else{//c
        }
    }
    if (bUnfType==1) {
        //CSVv2L
        if (flav==0) {//b and g->bb
            //if (pt>=100 && pt<170){
            if (pt>=100 && pt<175){//analysis_20180126_v0_p20180219_UMD_BtagEff_r1
                //eff_ = 1.13252-8.51292*pow(10,-3)*pt+7.06146*pow(10,-5)*pow(pt,2)-1.88610*pow(10,-7)*pow(pt,3);
                eff_ = 1.00020-5.43908*pow(10,-3)*pt+4.81536*pow(10,-5)*pow(pt,2)-1.35117*pow(10,-7)*pow(pt,3);
            }
            //else if (pt>=170 && pt<290){
            else if (pt>=175 && pt<290){//analysis_20180126_v0_p20180219_UMD_BtagEff_r1
                //eff_ = 1.17946-5.28402*pow(10,-3)*pt+2.38701*pow(10,-5)*pow(pt,2)-3.49534*pow(10,-8)*pow(pt,3);
                eff_ = 0.910730-1.71551*pow(10,-3)*pt+8.38286*pow(10,-6)*pow(pt,2)-1.28359*pow(10,-8)*pow(pt,3);
            }
            else if (pt>=290 && pt<575){
                //eff_ = -0.298764+1.01243*pow(10,-2)*pt-3.42044*pow(10,-5)*pow(pt,2)+5.08782*pow(10,-8)*pow(pt,3)-2.81822*pow(10,-11)*pow(pt,4);
                eff_ = 0.211333+5.30634*pow(10,-3)*pt-1.72964*pow(10,-5)*pow(pt,2)+2.47381*pow(10,-8)*pow(pt,3)-1.32797*pow(10,-11)*pow(pt,4);
            }
            else if (pt>=575){
                //eff_ = 0.224663+2.34533*pow(10,-3)*pt-3.04098*pow(10,-6)*pow(pt,2)+1.25291*pow(10,-9)*pow(pt,3);
                eff_ = 0.267642+2.14461*pow(10,-3)*pt-2.77930*pow(10,-6)*pow(pt,2)+1.12336*pow(10,-9)*pow(pt,3);
            }
        }
        else if (flav==1) {//udsg+c
            //if (pt>=100 && pt<170){
            if (pt>=100 && pt<165){//analysis_20180126_v0_p20180219_UMD_BtagEff_r1
                //eff_ = -0.790549+3.15976*pow(10,-2)*pt-4.13510*pow(10,-4)*pow(pt,2)+2.40094*pow(10,-6)*pow(pt,3)-5.13356*pow(10,-9)*pow(pt,4);
                eff_ = -1.85241+6.54939*pow(10,-2)*pt-8.18323*pow(10,-4)*pow(pt,2)+4.53507*pow(10,-6)*pow(pt,3)-9.32263*pow(10,-9)*pow(pt,4);
            }
            //else if (pt>=170 && pt<270){
            else if (pt>=165 && pt<270){//analysis_20180126_v0_p20180219_UMD_BtagEff_r1
                //eff_ = 1.08587-1.31000*pow(10,-2)*pt+4.27287*pow(10,-5)*pow(pt,2)+1.30308*pow(10,-7)*pow(pt,3)-9.42234*pow(10,-10)*pow(pt,4)+1.35103*pow(10,-12)*pow(pt,5);
                eff_ = 2.09104-2.76730*pow(10,-2)*pt+9.39869*pow(10,-5)*pow(pt,2)+2.77580*pow(10,-7)*pow(pt,3)-2.13824*pow(10,-9)*pow(pt,4)+3.19607*pow(10,-12)*pow(pt,5);
            }
            else if (pt>=270 && pt<475){
                //eff_ = 0.197115-6.87138*pow(10,-4)*pt+2.86535*pow(10,-6)*pow(pt,2)-3.07054*pow(10,-9)*pow(pt,3);
                eff_ = 0.202451-7.84820*pow(10,-4)*pt+3.19867*pow(10,-6)*pow(pt,2)-3.38650*pow(10,-9)*pow(pt,3);
            }
            else if (pt>=475){
                //eff_ = 0.395237-1.48947*pow(10,-3)*pt+3.47386*pow(10,-6)*pow(pt,2)-3.10196*pow(10,-9)*pow(pt,3)+9.63860*pow(10,-13)*pow(pt,4);
                eff_ = 0.386765-1.46702*pow(10,-3)*pt+3.47448*pow(10,-6)*pow(pt,2)-3.12835*pow(10,-9)*pow(pt,3)+9.79634*pow(10,-13)*pow(pt,4);
            }
        }
        else{//c
//             if (pt>=100 && pt<170){
//                 eff_ = *pow(10,)*pt  *pow(10,)*pow(pt,2) *pow(10,)*pow(pt,3) *pow(10,)*pow(pt,4) *pow(10,)*pow(pt,5) *pow(10,)*pow(pt,6);
//             }
//             else if (pt>=170 && pt<290){
//                 eff_ = *pow(10,)*pt  *pow(10,)*pow(pt,2) *pow(10,)*pow(pt,3) *pow(10,)*pow(pt,4) *pow(10,)*pow(pt,5) *pow(10,)*pow(pt,6);
//             }
//             else if (pt>=290 && pt<575){
//                 eff_ = *pow(10,)*pt  *pow(10,)*pow(pt,2) *pow(10,)*pow(pt,3) *pow(10,)*pow(pt,4) *pow(10,)*pow(pt,5) *pow(10,)*pow(pt,6);
//             }
//             else if (pt>=575){
//                 eff_ = *pow(10,)*pt  *pow(10,)*pow(pt,2) *pow(10,)*pow(pt,3) *pow(10,)*pow(pt,4) *pow(10,)*pow(pt,5) *pow(10,)*pow(pt,6);
//             }
        }
    }
    else if (bUnfType==2) {//CSVv2M
    }
    else if (bUnfType==3) {//CSVv2T
    }
    else {//Undefined
    }

    return eff_*sf_;
}


//double QCDhists::UnfoldWgtPtDep(int bUnfType, int nbtrue, int nbtagged, vector<float> *jetpt, bool isData){
bool QCDhists::UnfoldWgtPtDep(TMatrixD& MB, int bUnfType, vector<float> *jetpt, bool isData){

    bool status = true;
    double effs_[2][4];
    for (int i=0;i<4;i++){
        effs_[0][i]=effBTag(bUnfType,jetpt->at(i),0,isData);//b
        effs_[1][i]=effBTag(bUnfType,jetpt->at(i),1,isData);//light
//         effs_[0][i]=effBTagPara(bUnfType,jetpt->at(i),0,isData);//b
//         effs_[1][i]=effBTagPara(bUnfType,jetpt->at(i),1,isData);//light
//         effs_[0][i]=effCSVv2(bUnfType,jetpt->at(i),0,isData);//b
//         effs_[1][i]=effCSVv2(bUnfType,jetpt->at(i),1,isData);//light
//         effs_[0][i]=effDeepCSV(bUnfType,jetpt->at(i),0,isData);//b
//         effs_[1][i]=effDeepCSV(bUnfType,jetpt->at(i),1,isData);//light
//         effs_[0][i]=0.857;//data
//         effs_[1][i]=0.0899;//data
    }

    // Start of Matrix
    TMatrixD MA(5,5);

    //MA(0,0) = pow(1.0-effs_[1][0],4);
    MA(0,0) = (1.0-effs_[1][0])*(1.0-effs_[1][1])*(1.0-effs_[1][2])*(1.0-effs_[1][3]);

    //MA(0,1) = pow(1.0-effs_[0][0],1)*pow(1.0-effs_[1][0],3);
    MA(0,1)=0.0;
    for (int i=0;i<4;i++){
        double tmp=1.0-effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(1.0-effs_[1][j]);
            }
        }
        MA(0,1)+=(tmp/4.0);
    }

    //MA(0,2) = pow(1.0-effs_[0][0],2)*pow(1.0-effs_[1][0],2);
    MA(0,2)=0.0;
    for (int i=0;i<4;i++){
        double tmp=1.0-effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp1=tmp*(1.0-effs_[0][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp1*=(1.0-effs_[1][k]);
                    }
                }
                MA(0,2)+=(tmp1/12.0);
            }
        }
    }

    //MA(0,3) = pow(1.0-effs_[0][0],3)*pow(1.0-effs_[1][0],1);
    MA(0,3)=0.0;
    for (int i=0;i<4;i++){
        double tmp=1.0-effs_[1][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(1.0-effs_[0][j]);
            }
        }
        MA(0,3)+=(tmp/4.0);
    }

    //MA(0,4) = pow(1.0-effs_[0][0],4);
    MA(0,4) = (1.0-effs_[0][0])*(1.0-effs_[0][1])*(1.0-effs_[0][2])*(1.0-effs_[0][3]);

    //MA(1,0) = 4.0*pow(effs_[1][0],1)*pow(1.0-effs_[1][0],3);
    MA(1,0)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[1][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(1.0-effs_[1][j]);
            }
        }
        MA(1,0)+=tmp;
    }

    //MA(1,1) = effs_[0][0]*pow(1.0-effs_[1][0],3)+3.0*(1.0-effs_[0][0])*effs_[1][0]*pow(1.0-effs_[1][0],2);
    MA(1,1)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[0][i];
        double tmp10=(1-effs_[0][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp00*=(1.0-effs_[1][j]);
                double tmp11=tmp10*effs_[1][j];
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp11*=(1.0-effs_[1][k]);
                    }
                }
                MA(1,1)+=(tmp11/4.0);
            }
        }
        MA(1,1)+=(tmp00/4.0);
    }

    //MA(1,2) = 2.0*effs_[0][0]*(1.0-effs_[0][0])*pow(1.0-effs_[1][0],2)+2.0*pow(1.0-effs_[0][0],2)*effs_[1][0]*(1.0-effs_[1][0]);
    MA(1,2)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[0][i];
        double tmp10=effs_[1][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp01=tmp00*(1.0-effs_[0][j]);
                double tmp11=tmp10*(1.0-effs_[1][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp01*=(1.0-effs_[1][k]);
                        tmp11*=(1.0-effs_[0][k]);
                    }
                }
                MA(1,2)+=(tmp01/6.0);
                MA(1,2)+=(tmp11/6.0);
            }
        }
    }

    //MA(1,3) = effs_[1][0]*pow(1.0-effs_[0][0],3)+3.0*(1.0-effs_[1][0])*effs_[0][0]*pow(1.0-effs_[0][0],2);
    MA(1,3)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[1][i];
        double tmp10=(1-effs_[1][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp00*=(1.0-effs_[0][j]);
                double tmp11=tmp10*effs_[0][j];
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp11*=(1.0-effs_[0][k]);
                    }
                }
                MA(1,3)+=(tmp11/4.0);
            }
        }
        MA(1,3)+=(tmp00/4.0);
    }

    //MA(1,4) = 4.0*pow(effs_[0][0],1)*pow(1.0-effs_[0][0],3);
    MA(1,4)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(1.0-effs_[0][j]);
            }
        }
        MA(1,4)+=tmp;
    }

    //MA(2,0) = 6.0*pow(effs_[1][0],2)*pow(1.0-effs_[1][0],2);
    MA(2,0)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[1][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp1=tmp*effs_[1][j];
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp1*=(1.0-effs_[1][k]);
                    }
                }
                MA(2,0)+=(tmp1/2.0);
            }
        }
    }

    //MA(2,1) = 3.0*effs_[0][0]*effs_[1][0]*pow(1.0-effs_[1][0],2)    +3.0*pow(effs_[1][0],2)*(1.0-effs_[0][0])*(1.0-effs_[1][0]);
    MA(2,1)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[0][i];
        double tmp10=(1-effs_[0][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp01=tmp00*effs_[1][j];
                double tmp11=tmp10*(1.0-effs_[1][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp01*=(1.0-effs_[1][k]);
                        tmp11*=(effs_[1][k]);
                    }
                }
                MA(2,1)+=(tmp01/4.0);
                MA(2,1)+=(tmp11/4.0);
            }
        }
    }

    //MA(2,2) = pow(effs_[0][0],2)*pow(1.0-effs_[1][0],2)+4.0*effs_[0][0]*(1.0-effs_[0][0])*effs_[1][0]*(1.0-effs_[1][0])+pow(1.0-effs_[0][0],2)*pow(effs_[1][0],2);
    MA(2,2)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[0][i];
        double tmp10=(1-effs_[0][i]);
        double tmp20=effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp01=tmp00*effs_[0][j];
                double tmp11=tmp10*(1.0-effs_[0][j]);
                double tmp21=tmp20*(1.0-effs_[0][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp01*=(1.0-effs_[1][k]);
                        tmp11*=(effs_[1][k]);
                        double tmp22 = tmp21*(effs_[1][k]);
                        for (int l=0;l<4;l++){
                            if (l!=j && l!=i && l!=k){
                                tmp22*=(1-effs_[1][l]);
                            }
                        }
                        MA(2,2)+=(tmp22/6.0);
                    }
                }
                MA(2,2)+=(tmp01/12.0);
                MA(2,2)+=(tmp11/12.0);
            }
        }
    }

    //MA(2,3) = 3.0*effs_[1][0]*effs_[0][0]*pow(1.0-effs_[0][0],2)+3.0*(1.0-effs_[1][0])*(1.0-effs_[0][0])*pow(effs_[0][0],2);
    MA(2,3)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=effs_[1][i];
        double tmp10=(1-effs_[1][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp01=tmp00*effs_[0][j];
                double tmp11=tmp10*(1.0-effs_[0][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp01*=(1.0-effs_[0][k]);
                        tmp11*=(effs_[0][k]);
                    }
                }
                MA(2,3)+=(tmp01/4.0);
                MA(2,3)+=(tmp11/4.0);
            }
        }
    }

    //MA(2,4) = 6.0*pow(effs_[0][0],2)*pow(1.0-effs_[0][0],2);
    MA(2,4)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp1=tmp*effs_[0][j];
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp1*=(1.0-effs_[0][k]);
                    }
                }
                MA(2,4)+=(tmp1/2.0);
            }
        }
    }

    //MA(3,0) = 4.0*pow(effs_[1][0],3)*(1.0-effs_[1][0]);
    MA(3,0)=0.0;
    for (int i=0;i<4;i++){
        double tmp=(1-effs_[1][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(effs_[1][j]);
            }
        }
        MA(3,0)+=tmp;
    }

    //MA(3,1) = 3.0*effs_[0][0]*(1.0-effs_[1][0])*pow(effs_[1][0],2)+(1.0-effs_[0][0])*pow(effs_[1][0],3);
    MA(3,1)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=(1-effs_[0][i]);
        double tmp10=(effs_[0][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp00*=(effs_[1][j]);
                double tmp11=tmp10*(1-effs_[1][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp11*=(effs_[1][k]);
                    }
                }
                MA(3,1)+=(tmp11/4.0);
            }
        }
        MA(3,1)+=(tmp00/4.0);
    }


    //MA(3,2) = 2.0*pow(effs_[0][0],2)*(1.0-effs_[1][0])*effs_[1][0]+2.0*effs_[0][0]*(1.0-effs_[0][0])*pow(effs_[1][0],2);
    MA(3,2)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=(1-effs_[0][i]);
        double tmp10=(1-effs_[1][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp01=tmp00*(effs_[0][j]);
                double tmp11=tmp10*(effs_[1][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp01*=(effs_[1][k]);
                        tmp11*=(effs_[0][k]);
                    }
                }
                MA(3,2)+=(tmp01/6.0);
                MA(3,2)+=(tmp11/6.0);
            }
        }
    }


    //MA(3,3) = pow(effs_[0][0],3)*(1.0-effs_[1][0])+3.0*pow(effs_[0][0],2)*(1.0-effs_[0][0])*effs_[1][0];
    MA(3,3)=0.0;
    for (int i=0;i<4;i++){
        double tmp00=(1-effs_[1][i]);
        double tmp10=(effs_[1][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp00*=(effs_[0][j]);
                double tmp11=tmp10*(1-effs_[0][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp11*=(effs_[0][k]);
                    }
                }
                MA(3,3)+=(tmp11/4.0);
            }
        }
        MA(3,3)+=(tmp00/4.0);
    }

    //MA(3,4) = 4.0*pow(effs_[0][0],3)*(1.0-effs_[0][0]);
    MA(3,4)=0.0;
    for (int i=0;i<4;i++){
        double tmp=(1-effs_[0][i]);
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(effs_[0][j]);
            }
        }
        MA(3,4)+=tmp;
    }

    //MA(4,0) = pow(effs_[1][0],4);
    MA(4,0) = (effs_[1][0])*(effs_[1][1])*(effs_[1][2])*(effs_[1][3]);

    //MA(4,1) = pow(effs_[1][0],3)*pow(effs_[0][0],1);
    MA(4,1)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(effs_[1][j]);
            }
        }
        MA(4,1)+=(tmp/4.0);
    }

    //MA(4,2) = pow(effs_[1][0],2)*pow(effs_[0][0],2);
    MA(4,2)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[0][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                double tmp1=tmp*(effs_[0][j]);
                for (int k=0;k<4;k++){
                    if (k!=j && k!=i){
                        tmp1*=(effs_[1][k]);
                    }
                }
                MA(4,2)+=(tmp1/12.0);
            }
        }
    }

    //MA(4,3) = pow(effs_[1][0],1)*pow(effs_[0][0],3);
    MA(4,3)=0.0;
    for (int i=0;i<4;i++){
        double tmp=effs_[1][i];
        for (int j=0;j<4;j++){
            if (j!=i){
                tmp*=(effs_[0][j]);
            }
        }
        MA(4,3)+=(tmp/4.0);
    }

    //MA(4,4) = pow(effs_[0][0],4);
    MA(4,4) = (effs_[0][0])*(effs_[0][1])*(effs_[0][2])*(effs_[0][3]);
    // End of Matrix

    MA.SetTol(1.e-23);

    // Note: SVD can manipulate sigular matrix but potentially can also hide real problem.
    // Therefore, use TDecompLU as default matrix inversion check.
    TDecompLU lu(MA);
    //TDecompSVD svd(MA);
    //TMatrixD MB(5,5);
    if (!lu.Decompose()){
    //if (!svd.Decompose()){
        std::cout << "[Decompose ERROR]: Decomposition failed, matrix singular ?" << std::endl;
        status = false;
        for (int i=0;i<2;i++){
            for (int j=0;j<4;j++){
                std::cout << "[Decompose ERROR] Effs(" << i << "," << j <<")=" << effs_[i][j] << std::endl;
            }
        }
        for (int i=0;i<5;i++){
            for (int j=0;j<5;j++){
                std::cout << "[Decompose ERROR] MA(" << i << "," << j <<")=" << MA(i,j) << std::endl;
            }
        }
        for (int i=0;i<4;i++){
            std::cout << "[Decompose ERROR] Jet[" << i << "] pt = " << jetpt->at(i) << std::endl;
        }
        //if (isData) return UnfoldWgtD(bUnfType, nbtrue, nbtagged, 1);
        //else return UnfoldWgt(bUnfType, nbtrue, nbtagged, 1);
        if (isData) {
            for (int nbT=0;nbT<5;nbT++){
                for (int nbTagged=0;nbTagged<5;nbTagged++){
                    MB(nbT,nbTagged) = UnfoldWgtD(bUnfType, nbT, nbTagged, 1);
                }
            }
        } else {
            for (int nbT=0;nbT<5;nbT++){
                for (int nbTagged=0;nbTagged<5;nbTagged++){
                    MB(nbT,nbTagged) = UnfoldWgt(bUnfType, nbT, nbTagged, 1);
                }
            }
        }
    }
    else {
        lu.Invert(MB);
        //svd.Invert(MB);
//         for (int i=0;i<5;i++){
//             for (int j=0;j<5;j++){
//                 std::cout << "MB(" << i << "," << j <<")=" << MB(i,j) << std::endl;
//             }
//         }
    }

    //if (nbtagged>4) nbtagged=4;
    //return MB(nbtrue,nbtagged);
    return status;
}
