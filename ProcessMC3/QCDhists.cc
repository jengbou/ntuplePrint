#include <TList.h>
#include "QCDhists.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"
#include "EMJbkg.h"


const bool mergeOnly = false;
bool scanCuts = false;
bool verbose = false;

//void QCDhists()
void QCDhists(float goalintlum,int nbin, float* xsec, int* nfiles, std::string* binnames,std::string aaname,std::string ohname, int dooptk, int doopta, bool hasPre,bool donorm, bool blind, bool b16003, std::string bbname="./", bool crabformat=true)
{

    if (b16003) scanCuts = false;
    std::string inputfile;
    std::string outputfile;
  
    // opt
    float DHTcut=1000;
    float Dpt1cut=400;
    float Dpt2cut=200;
    float Dpt3cut=200;
    float Dpt4cut=100;
    float Dalphacut=0.04;
    float DmaxIPcut=0.4;
    float Djetacut = 2.;
    // dont forget there is a hidden cut nalmostemergin<4!!!!!!!!!!!!!!!!!
    int Dnemcut=2;
    int Dntrk1=0;
    // for alpha max scan
    const int ncutscan=3;
    //const int ncutscan=1;

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
            mkdir((bbname+binnames[i]).c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            for(int j=0;j<nfiles[i];j++) { //for each file for that bin
                //inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                //inputfile=aaname+binnames[i]+"/ntuple_"+binnames[i]+"_"+std::to_string(j+1)+"_hlt1p1.root";

                inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";//condor format
                if (crabformat) inputfile=aaname+binnames[i]+"/ntuple_"+std::to_string(j+1)+".root";//crab format
                std::cout<<"input file is "<<inputfile<<std::endl;
                outputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
                std::cout<<"output file is "<<outputfile<<std::endl;
                int itmp;
                if(!b16003) {
                    //itmp = EMJselect(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,DmaxIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
                    itmp = EMJbkg(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,DmaxIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
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
                    //inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                    inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";// condor format
                    if (crabformat) inputfile=aaname+binnames[i]+"/ntuple_"+std::to_string(j+1)+".root";// crab format
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
                        //inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.ntpl.root";
                        //inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+".root";
                        inputfile=aaname+binnames[i]+"/"+binnames[i]+"_"+std::to_string(j+1)+"_0.histo.root";// condor format
                        if (crabformat) inputfile=aaname+binnames[i]+"/ntuple_"+std::to_string(j+1)+".root";// crab format
                        std::cout<<"input file is "<<inputfile<<std::endl;
                        int iii=0;
                        outputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
                        std::cout<<"output file is "<<outputfile<<std::endl;

                        if(doopta==1) {
                            iii = EMJselect(false,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,acut2,DmaxIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
                        } else {
                            iii = EMJselect(false,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,acut2,0.9,0.9,Dntrk1,Dnemcut,blind);
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
        HistNorm(norm,nbin,xsec,nfiles,binnames,bbname);  // this gives the total number of events in each bin before all selections using the eventCountPreTrigger histogram
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
    vector<double> outnorm(nbin);
    for(int i=0;i<nhist;i++) {
        std::cout<<" entering Histman with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,binnames,donorm,bbname);
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
    vector<double> outnorm2(nbin);
    for(int i=0;i<nhist2;i++) {
        std::cout<<" entering Histman2 with i = "<<i<<": "<<histnames2[i]<<std::endl;
        vv2[i]=HistMan2(goalintlum,histnames2[i],norm,outnorm2,nbin,xsec,nfiles,binnames,donorm,bbname);
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
    outputfile=bbname+ohname;
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



TH1F* HistMan(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname) {

    std::string inputfile;


    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
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

TH2F* HistMan2(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname) {

    std::string inputfile;


    // now add up all the files for one bin
    std::cout<<" adding up histos within a bin"<<std::endl;
    vector<TH2F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
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

void  HistNorm(vector<double>& norm,int nbin,float* xsec, int* nfiles, std::string* binnames,std::string bbname) {

    std::cout<<"entering HistNorm"<<std::endl; 

    std::string inputfile;
    TFile * in;

    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    for(int i=0;i<nbin;i++) {  // for each bin
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
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


TH1F* HistMerge(float goalintlum,std::string thisHIST,vector<double>& norm,vector<double>& outnorm,int nbin,float* xsec, int* nfiles, std::string* binnames,bool donorm,std::string bbname) {

    std::string inputfile;

    // now add up all the files for one bin
    vector<TH1F> sum(nbin);
    TList *listO = new TList;
    for(int i=0;i<nbin;i++) {  // for each bin
        TList *listI = new TList;
        int ij=0;
        for(int j=0;j<nfiles[i];j++) { //for each file for that bin
            inputfile=bbname+binnames[i]+"/histos"+binnames[i]+"_"+std::to_string(j)+".root";
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

double fakerate(double jet_pt, double jet_eta, int jet_nTrack, int varType){
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
    else {//alpha2Dsig (tracksource=0, trackQuality HighPurity, track ipXYsig<3)
        if( fabs(jet_eta)<0.5 ){
            if( jet_pt>=100 && jet_pt<150 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.151857331395;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0989927500486;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.071523591876;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0583955533803;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0425348058343;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0295034591109;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0250466018915;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0202143471688;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0157502554357;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0120992455631;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0105512058362;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00902389641851;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00738134095445;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00649250904098;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00651684915647;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00583861209452;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00436769472435;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00473926169798;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00483116414398;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00298961787485;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00403584120795;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00540086813271;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00356606277637;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.002758448245;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00229456066154;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00256172404625;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00197994941846;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00971114449203;
                else if( jet_nTrack>=38 ) fakerate = 0.00296219531447;
            }
            else if( jet_pt>=150 && jet_pt<200 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.245900705457;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0934824049473;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0756023824215;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0618086606264;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0510251745582;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0403037928045;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0316348560154;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0263390913606;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0213516950607;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0182356443256;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0150372749195;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0120538091287;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0104210646823;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00909145642072;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00835578795522;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00565933668986;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00608934508637;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00485468888655;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00419934466481;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00377890979871;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00330415018834;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00306431250647;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00365268625319;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00324747688137;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00248104054481;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00133285927586;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00337457493879;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00430448306724;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00881730392575;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00483456626534;
                else if( jet_nTrack>=46 ) fakerate = 0.0163383018225;
            }
            else if( jet_pt>=200 && jet_pt<250 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.0603950731456;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0844824910164;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0948974788189;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0710218325257;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0576264187694;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0449043884873;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0417783968151;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0344189368188;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0275405570865;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0241658147424;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0207691341639;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0163192208856;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0136830052361;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0106671229005;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0100081181154;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00874896999449;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00668188603595;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00703618163243;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00593793205917;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00421985564753;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00547888828442;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00415584398434;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00480628293008;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00284560234286;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00375892105512;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00371308531612;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00295971496962;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00459376303479;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00030099612195;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00168187043164;
                else if( jet_nTrack>=46 ) fakerate = 0.00481773074716;
            }
            else if( jet_pt>=250 && jet_pt<300 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.00975985918194;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0889163687825;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.077287659049;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0619161836803;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0559390969574;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0472660847008;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0423096604645;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0382587201893;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0318798087537;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0257359798998;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0246854331344;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0197968482971;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0162832122296;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0146863460541;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0115514984354;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0101168779656;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00805226154625;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00857106037438;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0061517721042;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00660997955129;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00533122057095;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00489646662027;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00543434964493;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00495341140777;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00285317189991;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00411512888968;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00317823514342;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00194222375285;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00184705748688;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00366246490739;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0141051160172;
                else if( jet_nTrack>=50 ) fakerate = 0.00587621284649;
            }
            else if( jet_pt>=300 && jet_pt<400 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.189693003893;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0999341383576;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0727653950453;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0670783743262;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0549201965332;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0496036447585;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0437100119889;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0381193235517;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0323593355715;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0283249281347;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0240285694599;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0217919554561;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0186498723924;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.01554155536;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0138020562008;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0123977260664;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0103293675929;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00921337492764;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00868061278015;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00757455453277;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00661264546216;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00545204803348;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00504176737741;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0036203712225;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00445430912077;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00400694692507;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00368160195649;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00297693093307;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00169692130294;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00299716647714;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00236323359422;
                else if( jet_nTrack>=50 ) fakerate = 0.00256323674694;
            }
            else if( jet_pt>=400 && jet_pt<500 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.0175273567438;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.104379989207;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0692531764507;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.060599733144;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0490125864744;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0403060913086;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0348743721843;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0362033769488;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0316405333579;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0278960652649;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0242136958987;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0206768065691;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0193699002266;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0159613061696;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0156798008829;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0140598500147;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0120726674795;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0104231303558;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00981436111033;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00796204246581;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00635648239404;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00686525320634;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00606501242146;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00546380365267;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00564965326339;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00441046291962;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00390057801269;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00367193785496;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0024451338686;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00355843361467;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0027879839763;
                else if( jet_nTrack>=50 ) fakerate = 0.00263872276992;
            }
            else if( jet_pt>=500 && jet_pt<700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.175584703684;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0815467834473;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0669844672084;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0555565804243;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0450505875051;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0407588295639;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0362579748034;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0347638614476;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0303578712046;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0283873155713;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0251670498401;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0229269545525;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0204853601754;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0185879711062;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0162506587803;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0143141672015;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0139370402321;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0119601730257;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0112093212083;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0097376331687;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00934199616313;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00766484672204;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00648667570204;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.007010191679;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00661429250613;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00550322467461;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00518739642575;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00392627762631;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00348046771251;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00359245366417;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00278067705221;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00315450155176;
                else if( jet_nTrack>=60 ) fakerate = 0.00293253967538;
            }
            else if( jet_pt>=700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.114222630858;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0949421077967;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0694894194603;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0600205697119;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0495452024043;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0466119199991;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0428516678512;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0388750955462;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0374865233898;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0342103391886;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0330101698637;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0285010822117;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0266151353717;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0254352483898;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0215899460018;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0203823503107;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0191044155508;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0173943415284;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0160061642528;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0136652048677;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0131142884493;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0119175380096;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0113620702177;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0110817234963;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0103604663163;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00892030633986;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00747158238664;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00666214060038;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00569134531543;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00506417639554;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00466517545283;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00453108316287;
                else if( jet_nTrack>=60 ) fakerate = 0.00567397009581;
            }
        }
        else if( fabs(jet_eta)>=0.5 && fabs(jet_eta)<1.0 ){
            if( jet_pt>=100 && jet_pt<150 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.170513823628;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0976020768285;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0770787000656;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.057204246521;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0431796424091;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0331309251487;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0258806552738;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0206073224545;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0177502073348;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0148945013061;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0121969422325;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00943796988577;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00878859311342;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00725841475651;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00600283732638;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00531557388604;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00528359645978;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00483246147633;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00429958151653;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00634636124596;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00380366970785;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0048054610379;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00155170378275;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00244646891952;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00327797373757;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00168662902433;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00366014172323;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0067491075024;
                else if( jet_nTrack>=38 ) fakerate = 0.0026011634618;
            }
            else if( jet_pt>=150 && jet_pt<200 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.193289726973;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.101173572242;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0788148120046;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0692893788218;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.052091460675;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0425400100648;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0343787111342;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0282762311399;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0241883434355;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0182577539235;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0158818569034;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0141678936779;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.011212346144;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00922064296901;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00868694856763;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00755132874474;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00618250481784;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00589641043916;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00478685693815;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00460787676275;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0044647725299;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00517804035917;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00306901498698;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00398483220488;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00332207675092;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00575256440789;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00278268195689;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0035510500893;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00698330486193;
                else if( jet_nTrack>=42 ) fakerate = 0.0239578485489;
            }
            else if( jet_pt>=200 && jet_pt<250 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.143404647708;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0862248986959;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0786747038364;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0644471570849;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0570120960474;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0512683168054;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0424475371838;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0364868082106;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0299212373793;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0259703379124;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0197230428457;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0183239504695;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0159998033196;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0123596023768;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0108837746084;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0100736077875;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00872416049242;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00725273089483;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00591913377866;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00663360022008;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00580902863294;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00373453716747;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00520370062441;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00309303961694;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00333341071382;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0030821559485;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00398039352149;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00316792936064;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00333702308126;
                else if( jet_nTrack>=42 ) fakerate = 0.00241131777875;
            }
            else if( jet_pt>=250 && jet_pt<300 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.075891725719;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0698479264975;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0712507665157;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0714282989502;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0541564971209;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0515945032239;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0399643741548;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0362113416195;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0342532284558;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0248468089849;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.02255159989;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0205303970724;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0187212917954;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.015768866986;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0124216713011;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0117093548179;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00974681973457;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00831395294517;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00759842665866;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00644537620246;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00605357484892;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00601290073246;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0058848913759;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00485557597131;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00425033317879;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.003253323026;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00400642305613;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0049508754164;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00432567112148;
                else if( jet_nTrack>=42 ) fakerate = 0.0029760915786;
            }
            else if( jet_pt>=300 && jet_pt<400 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.15862031281;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0868663415313;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0702392011881;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0555565021932;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0510096848011;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0470600053668;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0398038849235;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.035717073828;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0344335250556;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0317068770528;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0260328073055;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0230134148151;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0194245949388;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0177654530853;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0153127908707;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0140515724197;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0109336785972;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0105117568746;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00937966350466;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00881176069379;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00753300078213;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00707871979102;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00491958623752;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00604562880471;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00566201796755;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00477601448074;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00444400589913;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00350999482907;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00357440277003;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0028207777068;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00100049667526;
                else if( jet_nTrack>=50 ) fakerate = 0.00834739021957;
            }
            else if( jet_pt>=400 && jet_pt<500 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.189870551229;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0981664136052;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0608892850578;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0536837279797;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0459804013371;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0416405163705;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0385137908161;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0336772501469;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0299430340528;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0290754213929;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0263321977109;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0237992573529;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0214492175728;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0195790473372;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0171814430505;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0159290563315;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0122665930539;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0126533117145;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0114986775443;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0100250476971;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00867976713926;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00832493044436;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00765977380797;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00644849054515;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00755128730088;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00551326526329;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00409991573542;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00456806458533;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0032770796679;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00403824215755;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00544906873256;
                else if( jet_nTrack>=50 ) fakerate = 0.00482672406361;
            }
            else if( jet_pt>=500 && jet_pt<700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.0224321167916;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0800021067262;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0585300773382;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0507373847067;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0452921055257;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0403737574816;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0382602475584;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0345542766154;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0323972478509;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0316301882267;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0280776135623;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.024643311277;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0231151971966;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0206130538136;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0189243871719;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.017852185294;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.015711504966;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0141834244132;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0131395412609;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0121249342337;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0114516075701;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0118307014927;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0100426701829;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00848437286913;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00682390108705;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00746431527659;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00603669742122;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00504211196676;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00516580091789;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0040119541809;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00456259539351;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00465831626207;
                else if( jet_nTrack>=60 ) fakerate = 0.00161913910415;
            }
            else if( jet_pt>=700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.079474106431;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.104029551148;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0704998746514;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0593690797687;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.050481069833;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0518011376262;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0463353469968;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0437184050679;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0418292284012;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0377959311008;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0357478633523;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0335546061397;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0317607596517;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0295002646744;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0265270397067;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0246847197413;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0237637981772;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0222049783915;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0185312107205;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0179794281721;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.016377305612;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0160262174904;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0147328441963;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0128872776404;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0116054145619;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0112696290016;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00962768029422;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00889393594116;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00760719878599;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00676980242133;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00640687486157;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00602737301961;
                else if( jet_nTrack>=60 ) fakerate = 0.00693205185235;
            }
        }
        else if( fabs(jet_eta)>=1.0 && fabs(jet_eta)<1.5 ){
            if( jet_pt>=100 && jet_pt<150 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.222389340401;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.102774843574;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0729825049639;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0541573613882;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0394843369722;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0311686042696;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0232803560793;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0190794952214;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0168103724718;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0136360125616;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0113313077018;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00978357437998;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0093952966854;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00796641316265;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00880806613714;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00751326512545;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00712850922719;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00688026752323;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00514424219728;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00762578239664;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00358663499355;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00658151227981;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00834454316646;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00514760380611;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00686286576092;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00422340724617;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00245882291347;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00291538354941;
                else if( jet_nTrack>=38 ) fakerate = 0.00604299549013;
            }
            else if( jet_pt>=150 && jet_pt<200 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.216000571847;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.092906832695;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0784682407975;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0631284490228;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.047132730484;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0405228622258;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0319824032485;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0283982586116;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0231099892408;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0219556875527;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0168165564537;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0155014405027;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0122931124642;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0115858409554;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.01004816778;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00848408322781;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00826018303633;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00774392066523;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0098865441978;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00443980377167;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00650695478544;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00457735685632;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.006362025626;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00535798678175;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00475359428674;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00579366600141;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0035765664652;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00669197738171;
                else if( jet_nTrack>=38 ) fakerate = 0.00604566326365;
            }
            else if( jet_pt>=200 && jet_pt<250 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.13898165524;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0907359942794;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0835166424513;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0614678375423;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0514428056777;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0435580871999;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0394417271018;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0338076464832;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0306124053895;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0263771042228;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0218295268714;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.020446959883;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0159538891166;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0163309667259;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0142868449911;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0118168992922;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0108088608831;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.010772690177;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00861599948257;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0103352107108;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00851218122989;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00826506596059;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00677710957825;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00686643319204;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00796253699809;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00571914436296;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00657793600112;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00504098599777;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00451089628041;
                else if( jet_nTrack>=42 ) fakerate = 0.00639342423528;
            }
            else if( jet_pt>=250 && jet_pt<300 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.311084747314;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0864741951227;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0698335170746;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0587332397699;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0555995553732;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0447252243757;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0423851050436;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0380578897893;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0338206775486;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0313702747226;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0262907408178;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0225215982646;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0219041313976;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0190496481955;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0174232274294;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0158133674413;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0135601898655;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.012844350189;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0112554160878;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0119822239503;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00840862281621;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00960650853813;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00811697728932;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00707910768688;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00946616567671;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0083574084565;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00680798431858;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0102671552449;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00448811519891;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0142447762191;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00423789443448;
                else if( jet_nTrack>=50 ) fakerate = 0.0187719147652;
            }
            else if( jet_pt>=300 && jet_pt<400 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.151099383831;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.113407030702;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0652764737606;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0610016696155;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0563226491213;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0517815984786;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.045717574656;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0408288240433;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0376935303211;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0342921279371;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.033896394074;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0310007613152;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0261331815273;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0244656354189;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0217231437564;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0211399272084;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0177857484668;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0163103658706;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0157126467675;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0147317247465;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0130281485617;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.012510092929;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0133888004348;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0101588591933;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0122444722801;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00998898688704;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00888221431524;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00991842709482;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00910107884556;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00726334564388;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0153241390362;
                else if( jet_nTrack>=50 ) fakerate = 0.000721593736671;
            }
            else if( jet_pt>=400 && jet_pt<500 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.201494991779;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0997964143753;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0723950490355;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0684645697474;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0575515441597;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0518182143569;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0470187179744;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0482326038182;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0432221665978;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0401954390109;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0385078899562;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0341330841184;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0318082757294;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0319540128112;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0293744541705;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0257902275771;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0249920673668;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0219820644706;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0205920916051;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0197712946683;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0169254224747;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.018394337967;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0181673131883;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0176311060786;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0150367449969;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0149011537433;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0131340082735;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0124378288165;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0111462585628;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.00981917977333;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0137569811195;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0170360393822;
                else if( jet_nTrack>=60 ) fakerate = 0.00993147026747;
            }
            else if( jet_pt>=500 && jet_pt<700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.293515563011;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.109293736517;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0831679552794;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0695674419403;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0629605799913;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0595296211541;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0585890635848;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0518032684922;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0502379797399;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0484657734632;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0470886975527;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0433531031013;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0428042151034;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0414185971022;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0380354225636;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0356326773763;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0316658429801;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0317451730371;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0286677796394;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0267685577273;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0276807099581;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0242429710925;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0231446418911;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.021554633975;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0209287572652;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0195045061409;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0185017306358;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.01716238074;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0173185002059;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0169545002282;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0210290849209;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0171519648284;
                else if( jet_nTrack>=60 ) fakerate = 0.00762183777988;
            }
            else if( jet_pt>=700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.140940442681;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.114723570645;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.107684403658;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0811650902033;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0740848481655;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.076260946691;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0686476528645;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0689130201936;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0680393353105;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0665683820844;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0643426701427;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0637082308531;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0625579729676;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0601998493075;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0585605315864;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0579168386757;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0551403313875;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0502668581903;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0457830242813;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0462974533439;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0457593090832;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.045802809298;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0387262776494;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0384150668979;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0376512035728;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0344886891544;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0350197218359;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0301985722035;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0274092946202;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0270156636834;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0261100456119;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0278373602778;
                else if( jet_nTrack>=60 ) fakerate = 0.0275077018887;
            }
        }
        else if( fabs(jet_eta)>=1.5 && fabs(jet_eta)<2.0 ){
            if( jet_pt>=100 && jet_pt<150 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.123499549925;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.11568620801;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0869610905647;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0739129707217;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.059928715229;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0476401783526;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0414832830429;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0356137678027;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0300451815128;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0258584488183;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0199110563844;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0185232404619;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0155907198787;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0136446589604;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0122634544969;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0110163362697;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0106465369463;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00958545226604;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00967862829566;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00791541300714;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00867203529924;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00701516726986;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00763392588124;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00529406731948;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00272674416192;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00818560365587;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00207397108898;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0100735463202;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00146801641677;
                else if( jet_nTrack>=42 ) fakerate = 0.0377774909139;
            }
            else if( jet_pt>=150 && jet_pt<200 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.21923814714;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.123739510775;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.084182061255;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0744320452213;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0615046061575;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0542178861797;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.05254271999;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0426488891244;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0384915806353;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0369292907417;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0305441897362;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0255582090467;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0231840927154;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0233526621014;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0204962901771;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0152162667364;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0137805938721;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.012135906145;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.011738512665;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0115396054462;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00953314453363;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0116667980328;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00593753857538;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00420593144372;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00971141643822;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00451125763357;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00656423019245;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00443779537454;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0103693818673;
                else if( jet_nTrack>=42 ) fakerate = 0.0458075776696;
            }
            else if( jet_pt>=200 && jet_pt<250 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.264200180769;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.15654283762;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0910488814116;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0742888525128;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0624251142144;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0568093135953;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0563772730529;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0474453791976;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0446477718651;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.043346952647;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0376808047295;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0344673544168;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0334555022418;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0282824914902;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0262836068869;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.019550235942;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0190424602479;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0198956653476;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0156062142923;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0140284448862;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0140475817025;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0111948428676;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0108061162755;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.011536010541;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0125310616568;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.00768739031628;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.00990822911263;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.00581849506125;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0163407344371;
                else if( jet_nTrack>=42 ) fakerate = 0.0014725306537;
            }
            else if( jet_pt>=250 && jet_pt<300 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.105601385236;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.125779345632;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0784260630608;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0641886740923;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0596063695848;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0563928224146;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.055113542825;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0519326552749;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0486539378762;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0453642345965;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0418092831969;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0376684777439;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0348563492298;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0325388200581;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.029140220955;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0253713838756;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0256521385163;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0210781376809;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0232215560973;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0200850050896;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0169396121055;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0170668326318;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.015671633184;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0145198553801;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0148998592049;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0128045333549;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.010105419904;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0120368897915;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00550971133634;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0126441279426;
                else if( jet_nTrack>=46 ) fakerate = 0.0201187320054;
            }
            else if( jet_pt>=300 && jet_pt<400 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.073428593576;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.165598109365;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.085615478456;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0742686837912;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.058540482074;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0584048703313;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0593656860292;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0541160739958;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0544582866132;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0480604469776;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0458510406315;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0437284857035;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0384068787098;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0393913984299;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0336638465524;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0335330180824;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0314447954297;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0286570414901;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0255417190492;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0236482191831;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0243055168539;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.023555861786;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0180069599301;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0182526800781;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0159663949162;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.016789957881;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0141255129129;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0095970723778;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.00915038585663;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0100325206295;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0085044587031;
                else if( jet_nTrack>=50 ) fakerate = 0.00298135099001;
            }
            else if( jet_pt>=400 && jet_pt<500 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.0300003960729;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.130307257175;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.082386687398;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0723518654704;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0663996860385;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0645938441157;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0536122322083;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0543504543602;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0556614510715;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0510120950639;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0516452677548;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0456601604819;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0483054742217;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0424055084586;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0436770804226;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.040022470057;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0395106561482;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0351606719196;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0361736305058;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0287842154503;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0285102576017;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0260046906769;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0268041938543;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0262372978032;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0246757231653;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0197156481445;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0174688734114;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0159055367112;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0110204778612;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0137908607721;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.00727520510554;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0159557051957;
                else if( jet_nTrack>=60 ) fakerate = 0.00560536561534;
            }
            else if( jet_pt>=500 && jet_pt<700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.227905154228;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.121943034232;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0818708539009;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0638139173388;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0691694915295;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0667762607336;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0673333033919;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0639580786228;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0676628798246;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0640733242035;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0608435422182;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0612069144845;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0599221922457;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0536738485098;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0532978065312;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0543997958302;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0480499975383;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0441527664661;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0433960147202;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0416267812252;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0397669412196;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0380749218166;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0399314947426;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0338355302811;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0271122027189;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0271069444716;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0278930552304;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.019773805514;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0186952482909;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.015772011131;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.0155677152798;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0116945821792;
                else if( jet_nTrack>=60 ) fakerate = 0.0146919060498;
            }
            else if( jet_pt>=700 ){
                if( jet_nTrack>0 && jet_nTrack<2 ) fakerate = 0.0;
                else if( jet_nTrack>=2 && jet_nTrack<4 ) fakerate = 0.0728413388133;
                else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0911043807864;
                else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0964451804757;
                else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0842428132892;
                else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0872284695506;
                else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0843918398023;
                else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.083222463727;
                else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0860547497869;
                else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0784920081496;
                else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0818098932505;
                else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0795335397124;
                else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0788054019213;
                else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0782911255956;
                else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0730749890208;
                else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0718063637614;
                else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0706673264503;
                else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0696748867631;
                else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0629977211356;
                else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0630459487438;
                else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0629019141197;
                else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0564953461289;
                else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0545506551862;
                else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0524924844503;
                else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0501857474446;
                else if( jet_nTrack>=30 && jet_nTrack<32 ) fakerate = 0.0477095879614;
                else if( jet_nTrack>=32 && jet_nTrack<35 ) fakerate = 0.0449750721455;
                else if( jet_nTrack>=35 && jet_nTrack<38 ) fakerate = 0.0384937077761;
                else if( jet_nTrack>=38 && jet_nTrack<42 ) fakerate = 0.0351780094206;
                else if( jet_nTrack>=42 && jet_nTrack<46 ) fakerate = 0.0269066989422;
                else if( jet_nTrack>=46 && jet_nTrack<50 ) fakerate = 0.023744482547;
                else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.020927419886;
                else if( jet_nTrack>=60 ) fakerate = 0.0228255633265;
            }
        }
    }

    return fakerate;

}

double fakerateTP(double jet_pt, double jet_eta, int jet_nTrack, int varType){
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
    else {//alpha2Dsig (tracksource=0, trackQuality HighPurity, track ipXYsig<3)
        if( jet_pt>=100 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.146489977837;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.124317288399;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.103396132588;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0840845927596;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0670311823487;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0543308518827;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0459702052176;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0382093004882;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0299763418734;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0242423731834;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0204917676747;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0173804014921;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0152069283649;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0141417263076;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0124218799174;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0114593682811;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0124008366838;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0116524109617;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00932095758617;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00963907688856;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0120211420581;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0112675875425;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00953529030085;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0134123973548;
            else if( jet_nTrack>=30 && jet_nTrack<35 ) fakerate = 0.010936117731;
            else if( jet_nTrack>=35 && jet_nTrack<42 ) fakerate = 0.0220378506929;
            else if( jet_nTrack>=42 ) fakerate = 0.0937880277634;
        }
        else if( jet_pt>=200 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.130365163088;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.132416114211;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.114442609251;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0946234315634;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0883239805698;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0799155980349;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0667744278908;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0589411668479;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0517020821571;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0400869362056;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0329704545438;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0291685108095;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0287834238261;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.022603482008;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0186357628554;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0181780718267;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0145326927304;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.010152939707;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.013935171999;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0098467124626;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0108909485862;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0119141079485;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0104302698746;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0100729605183;
            else if( jet_nTrack>=30 && jet_nTrack<35 ) fakerate = 0.00885492935777;
            else if( jet_nTrack>=35 && jet_nTrack<42 ) fakerate = 0.0122162709013;
            else if( jet_nTrack>=42 && jet_nTrack<50 ) fakerate = 0.0432729423046;
            else if( jet_nTrack>=50 ) fakerate = 0.0196891464293;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.123541176319;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.121213160455;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.114801801741;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.09071765095;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0896960049868;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.082162424922;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0651457458735;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0622268691659;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0543217472732;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0465599000454;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.040523018688;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.03957490623;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0329116247594;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0251397248358;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0245438702404;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0211482346058;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0184361729771;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.018400894478;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0138794146478;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.015748500824;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.013408601284;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0113888029009;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0096627715975;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0133152222261;
            else if( jet_nTrack>=30 && jet_nTrack<35 ) fakerate = 0.0119673470035;
            else if( jet_nTrack>=35 && jet_nTrack<42 ) fakerate = 0.010244095698;
            else if( jet_nTrack>=42 && jet_nTrack<50 ) fakerate = 0.00985997635871;
            else if( jet_nTrack>=50 ) fakerate = 0.0117689948529;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.113848514855;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0781108811498;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.092665605247;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.07217040658;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0651525780559;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0624618008733;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0606939457357;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0517875105143;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0477401874959;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0435847081244;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0359434783459;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0345578230917;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0319542661309;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0298171639442;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0249006766826;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0221967771649;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0177872478962;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0185889396816;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0133231794462;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0154745690525;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0120643116534;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0126705411822;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0145116290078;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0110222976655;
            else if( jet_nTrack>=30 && jet_nTrack<35 ) fakerate = 0.00861669424921;
            else if( jet_nTrack>=35 && jet_nTrack<42 ) fakerate = 0.00974339526147;
            else if( jet_nTrack>=42 && jet_nTrack<50 ) fakerate = 0.00983470957726;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.00937225110829;
            else if( jet_nTrack>=60 ) fakerate = 0.0143514527008;
        }
        else if( jet_pt>=500 ){
            if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.127319455147;
            else if( jet_nTrack>=4 && jet_nTrack<6 ) fakerate = 0.0954762846231;
            else if( jet_nTrack>=6 && jet_nTrack<8 ) fakerate = 0.0602521598339;
            else if( jet_nTrack>=8 && jet_nTrack<10 ) fakerate = 0.0661811977625;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0658358260989;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0596187226474;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0538003332913;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0530491471291;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0521784685552;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0447088964283;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0415754318237;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0405207909644;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0344916060567;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0315645970404;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0302998442203;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0261623077095;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0266867168248;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0240785628557;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0200217217207;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0212977267802;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0177422259003;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0185490883887;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0163516066968;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0139737566933;
            else if( jet_nTrack>=30 && jet_nTrack<35 ) fakerate = 0.0123072881252;
            else if( jet_nTrack>=35 && jet_nTrack<42 ) fakerate = 0.0116289332509;
            else if( jet_nTrack>=42 && jet_nTrack<50 ) fakerate = 0.0109352841973;
            else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0107339816168;
            else if( jet_nTrack>=60 ) fakerate = 0.0145338568836;
        }
    }
    return fakerate;

}

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

double frWeight(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<vector<float> > *track_pt, int njetscut, double jptcut, int varType){
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

double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeight1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeightT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeightT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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

double frWeightT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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
        for(Int_t k=j+1; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            double p22 = p21 * fakerate(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType);
            for(Int_t l=0; l<njets; l++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                p22 *= (1.0-fakerate(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType));
            }
            p2 += p22;
        }
    }
    return p2;
}

//Ntag==2; tag-and-probe
double frWeightT21(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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
            for(Int_t j1=0; j1<njets; j1++) {
                if (j1==j || j1==k) continue;
                int jdx1 = goodjetIdx[j1];
                int ntrks11 = ntrack[jdx1];
                p12 *= (1.0-fakerate(jetpt->at(jdx1),jeteta->at(jdx1),ntrks11,varType));
            }

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

double frWeightT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, int njetscut, double jptcut, int varType){
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
                    p33 *= (1.0-fakerate(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType));
                }
                p3 += p33;
            }
        }
    }
    return p3;
}

// Only 4 leading jets
double frWeight4(vector<float> *jetpt, vector<float> *jeteta, vector<bool> *basicjet, vector<int> &ntrack, double jptcut, int varType){
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


double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_pvWeight)
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

double GetAlpha(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
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

double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality, vector<float> &track_ipXYSigs)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        ptsum_total += track_pt.at(itk);
        if ( track_ipXYSigs.at(itk) < 3.0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha2Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
                     vector<float> &track_ipXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        if (fabs(pv_z-track_ref_zs.at(itk))>pilecut) continue;// remove tracks with exceedingly large z
        ptsum_total += track_pt.at(itk);
        if ( track_ipXYSigs.at(itk) < 3.0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha2Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

double nGJrewgt(int nGoodJet){
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
