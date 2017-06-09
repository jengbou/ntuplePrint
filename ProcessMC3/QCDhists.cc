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


double fakerateF(double jet_pt, double jet_eta, int jet_nTrack, int varType, int flav){
    double fakerate = 0.0;
    if (fabs(flav)==5 || fabs(flav)==8) {// default: b and g->bb
        if( jet_pt>=100 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.715059876442;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.629313707352;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.528343737125;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.507715523243;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.469951748848;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.439262956381;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.401671767235;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.363309115171;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.326585620642;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.296461105347;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.258569836617;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.222324073315;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.195010304451;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.163752108812;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.140868544579;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.119272492826;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.100279137492;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0826535075903;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0638188943267;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0559257231653;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0424351431429;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0365692898631;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0293366611004;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.022599350661;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0171766523272;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0168513003737;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0122290300205;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00733736064285;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00623644841835;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.00793581176549;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.00518085062504;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.00479727284983;
            else if( jet_nTrack>=33 && jet_nTrack<36 ) fakerate = 0.00976078398526;
        }
        else if( jet_pt>=200 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.766644239426;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.661502659321;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.497469037771;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.494727581739;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.495208054781;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.44789364934;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.425969421864;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.400438308716;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.362398028374;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.337267398834;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.306725084782;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.275644272566;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.246207788587;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.21221858263;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.187123671174;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.162651032209;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.14050988853;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.117521502078;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0994837731123;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.0853370204568;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0699598491192;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0596641488373;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0470349006355;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0376072451472;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0317801348865;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0238691847771;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0205276757479;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.016407025978;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0100565757602;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0121419224888;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.00575947761536;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.0057715815492;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.00471845502034;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.00305754481815;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00598808191717;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.00464069005102;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00156007974874;
            else if( jet_nTrack>=38 && jet_nTrack<40 ) fakerate = 9.03845648281e-06;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.591901659966;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.463596314192;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.585190534592;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.491322517395;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.461325109005;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.433495372534;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.42331725359;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.384311169386;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.343639194965;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.32081386447;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.294203609228;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.27275159955;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.245075434446;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.218041613698;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.194739356637;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.171788409352;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.158191621304;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.13571421802;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.119220420718;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.103109136224;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0881746262312;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0768402591348;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0597792975605;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0515190511942;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0410923399031;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0334663018584;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0284008923918;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0203177854419;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0163359735161;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0146741019562;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.0126445731148;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.00880257040262;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.00698458123952;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.00431843055412;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00560474162921;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.00341084669344;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00200200802647;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.00312445149757;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.00133896083571;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.00259156827815;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.00330586731434;
            else if( jet_nTrack>=42 && jet_nTrack<44 ) fakerate = 0.00202568154782;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.499583095312;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.450487554073;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.500702619553;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.491060495377;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.486063927412;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.439745545387;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.378054827452;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.344892859459;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.342964202166;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.309404373169;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.276457190514;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.262597471476;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.231651797891;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.224189996719;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.197939813137;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.181883990765;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.166576698422;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.141870379448;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.127526342869;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.113013081253;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.092493519187;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0846745967865;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0732051581144;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0595168806612;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0533206872642;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0408093817532;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0308560002595;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0300334058702;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0240820683539;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0201674029231;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.0141323199496;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.0135890766978;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.0101735536009;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.00839622691274;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00733977090567;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.00661996938288;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00415988545865;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.00545014487579;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.00192761630751;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.00531906681135;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.00155104300939;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 0.000190330873011;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.000231184167205;
            else if( jet_nTrack>=44 && jet_nTrack<46 ) fakerate = 0.000426233513281;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.00641486980021;
        }
        else if( jet_pt>=500 && jet_pt<700 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.920715034008;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.581622302532;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.324101954699;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.541637420654;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.455724567175;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.378196239471;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.401509702206;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.34158629179;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.326691269875;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.305788725615;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.283505499363;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.251563817263;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.252281159163;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.192850738764;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.195830658078;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.178407281637;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.170073494315;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.147526562214;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.136231586337;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.114961661398;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0962586924434;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0835715830326;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0720031931996;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0673136338592;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0601850897074;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0456812269986;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0424221903086;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0280471350998;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0296144168824;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0237792115659;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.0169200859964;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.0151808932424;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.0160291325301;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.0088355075568;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00864937249571;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.00590739864856;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00558930682018;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.0037547533866;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.00624633487314;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.00242132181302;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.00371199380606;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 0.00232599675655;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.00207602279261;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 0.00138432916719;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 0.000106653409603;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 0.000885473156814;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.00187107559759;
            else if( jet_nTrack>=48 && jet_nTrack<52 ) fakerate = 0.00132475944702;
            else if( jet_nTrack>=52 && jet_nTrack<53 ) fakerate = 0.000343284773408;
            else if( jet_nTrack>=53 && jet_nTrack<54 ) fakerate = 0.000471311243018;
            else if( jet_nTrack>=54 && jet_nTrack<55 ) fakerate = 0.000542446679901;
            else if( jet_nTrack>=55 && jet_nTrack<56 ) fakerate = 0.00513691082597;
        }
        else if( jet_pt>=700 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.0;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.272499978542;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.351407468319;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.384297341108;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.32820519805;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.333652466536;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.30006146431;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.229936897755;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.223351731896;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.219748437405;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.217988535762;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.195175051689;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.185304924846;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.179233416915;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.157421916723;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.151232481003;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.126067653298;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.117333039641;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.111748322845;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.104287564754;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.0923201143742;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0840925723314;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.0726553052664;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0625237897038;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.0589427165687;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.0544106028974;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.0398165099323;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.0375618040562;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.0290282871574;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0292799528688;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.0226878318936;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.0195199176669;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.0179017689079;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.014044621028;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00960959307849;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.0100202020258;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00780270900577;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.00687822094187;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.00579924928024;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.00496913213283;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.00317132775672;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 0.00394399557263;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.0030550006777;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 0.00284199137241;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 0.00232639769092;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 0.00255473703146;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.000801635440439;
            else if( jet_nTrack>=48 && jet_nTrack<49 ) fakerate = 0.00107172341086;
            else if( jet_nTrack>=49 && jet_nTrack<50 ) fakerate = 0.00165503879543;
            else if( jet_nTrack>=50 && jet_nTrack<51 ) fakerate = 0.00096572691109;
            else if( jet_nTrack>=51 && jet_nTrack<52 ) fakerate = 0.00229308777489;
            else if( jet_nTrack>=52 && jet_nTrack<53 ) fakerate = 0.00216134591028;
            else if( jet_nTrack>=53 && jet_nTrack<54 ) fakerate = 0.00127281574532;
            else if( jet_nTrack>=54 && jet_nTrack<55 ) fakerate = 0.00114713376388;
            else if( jet_nTrack>=55 && jet_nTrack<56 ) fakerate = 0.00040523026837;
            else if( jet_nTrack>=56 && jet_nTrack<57 ) fakerate = 0.00117483024951;
            else if( jet_nTrack>=57 && jet_nTrack<58 ) fakerate = 0.00140666321386;
            else if( jet_nTrack>=58 && jet_nTrack<59 ) fakerate = 0.00161025836132;
            else if( jet_nTrack>=59 && jet_nTrack<60 ) fakerate = 0.00293966941535;
            else if( jet_nTrack>=60 && jet_nTrack<62 ) fakerate = 0.00260410434566;
            else if( jet_nTrack>=62 && jet_nTrack<63 ) fakerate = 0.0;
            else if( jet_nTrack>=63 && jet_nTrack<64 ) fakerate = 0.0;
            else if( jet_nTrack>=64 && jet_nTrack<65 ) fakerate = 0.0;
            else if( jet_nTrack>=65 && jet_nTrack<66 ) fakerate = 0.0;
            else if( jet_nTrack>=66 && jet_nTrack<67 ) fakerate = 0.0045856917277;
            else if( jet_nTrack>=67 && jet_nTrack<68 ) fakerate = 0.0;
            else if( jet_nTrack>=68 && jet_nTrack<69 ) fakerate = 0.0;
        }
    }
    else {//g,u,d,c,s
        if( jet_pt>=100 && jet_pt<200 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.116929680109;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.0776538178325;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.0538891553879;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0415382906795;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0316988527775;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.0250606592745;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0193945094943;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0146419722587;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0112940631807;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.00823308341205;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.00598462251946;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00450083101168;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00336895207874;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00252591259778;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00190512055997;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00143598543946;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0010693782242;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.000833429046907;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.000605412409641;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.000423416495323;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.000423510908149;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.000318155100103;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00023320899345;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.000230387493502;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.000127563762362;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.000122862751596;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.000105378574517;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.000107260457298;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 4.58786780655e-05;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 3.37183860211e-06;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 4.77731509818e-06;
        }
        else if( jet_pt>=200 && jet_pt<300 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.108254872262;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.074730053544;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.0513213984668;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0410384088755;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0327457636595;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.027534108609;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0222020410001;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0185786578804;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0156719554216;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0125106573105;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0104531822726;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00823249574751;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00651287473738;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00500586302951;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00411078240722;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00320230936632;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00246857269667;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00184428191278;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00141112843994;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00110363855492;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00084279559087;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.000741443480365;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.000489442143589;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.000401988858357;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.000336873577908;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.000284847978037;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.000250481680268;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.000155415211339;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.000132717294036;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.000159239934874;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 4.78518049931e-05;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 4.54061628261e-05;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 8.17726831883e-05;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.00017799444322;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 5.16239015269e-05;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.0;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 6.78041487845e-07;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.0;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 1.45537467233e-06;
        }
        else if( jet_pt>=300 && jet_pt<400 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.0907843187451;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.0560795851052;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.049036514014;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0390446595848;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0304229017347;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.0255571380258;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0214924663305;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0184146333486;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0162382181734;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.013980647549;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0114672658965;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.00995473284274;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.00833496358246;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00704440334812;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00598440179601;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00464624864981;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00395054928958;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00314948242158;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00249475007877;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00206430442631;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00155460077804;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00131264445372;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00106300867628;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.000933959672693;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.000744890305214;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.000631943810731;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00044018618064;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.000346210144926;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.000320376100717;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.000289181567496;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.000208661906072;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.000182903779205;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.000143267461681;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.000136656686664;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 6.0026657593e-05;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 5.18191009178e-05;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 6.88793224981e-05;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.000145528509165;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 7.97821194283e-05;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 5.24160786881e-05;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 6.79302829667e-05;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 4.81825736642e-06;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.0;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 0.0;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 0.0;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 2.19233443204e-05;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.0;
            else if( jet_nTrack>=48 && jet_nTrack<49 ) fakerate = 2.24663085646e-06;
        }
        else if( jet_pt>=400 && jet_pt<500 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.0994534566998;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.0532199628651;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.0447547473013;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0352230519056;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0299788024276;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.0244853626937;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0217329263687;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0187269691378;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0166731029749;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0149432225153;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0131156900898;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0118324719369;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0102959554642;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.00897935125977;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.00799846835434;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.00683767907321;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0058213295415;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00467379018664;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00422027660534;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00333816837519;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00273886718787;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00230916379951;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00178867008071;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00155322323553;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00129058444872;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.000987972831354;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.000830858596601;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00069657311542;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.000539968721569;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.000432714936323;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.000362792372471;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.000262513203779;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.000260592787527;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.000281449873;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.000257843319559;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.000146337828483;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.000147041253513;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 5.30376455572e-05;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 9.99950643745e-05;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 7.82359347795e-05;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 1.43478227983e-05;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 9.39989986364e-05;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 2.33539194596e-05;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 2.46599865932e-05;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 0.0;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 0.0;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.0;
            else if( jet_nTrack>=48 && jet_nTrack<49 ) fakerate = 2.6267618523e-05;
            else if( jet_nTrack>=49 && jet_nTrack<50 ) fakerate = 0.0014248683583;
            else if( jet_nTrack>=50 && jet_nTrack<51 ) fakerate = 0.0;
            else if( jet_nTrack>=51 && jet_nTrack<52 ) fakerate = 7.67048986745e-05;
        }
        else if( jet_pt>=500 && jet_pt<700 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.0858066082001;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.0586804300547;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.0443582907319;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0356691963971;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0337994284928;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.0284771174192;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0246626529843;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0231386981905;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0207361020148;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0192219540477;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.0170571003109;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.015814114362;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0143355354667;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0131427785382;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0116070406511;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0100796399638;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.00938135385513;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.00794247630984;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.00697147566825;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00603941828012;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00516418647021;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.00450647808611;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00315080280416;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.0031693787314;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00259743374772;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00213916227221;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00197701016441;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00135822803713;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00134331779554;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.00108657393139;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.000713593210094;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.00102262641303;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.000524442468304;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.000803839589935;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.000555855804123;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.000357956392691;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.000303968816297;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.000543461123016;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.000191067476408;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.000214508268982;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.000183827185538;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 0.000195658096345;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.000156530193635;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 2.94724686682e-05;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 9.50700050453e-05;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 9.79013784672e-05;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.000185092809261;
            else if( jet_nTrack>=48 && jet_nTrack<49 ) fakerate = 3.53252253262e-05;
            else if( jet_nTrack>=49 && jet_nTrack<50 ) fakerate = 0.000297106569633;
            else if( jet_nTrack>=50 && jet_nTrack<51 ) fakerate = 2.82178989437e-05;
            else if( jet_nTrack>=51 && jet_nTrack<52 ) fakerate = 0.000144546211231;
            else if( jet_nTrack>=52 && jet_nTrack<53 ) fakerate = 0.0;
            else if( jet_nTrack>=53 && jet_nTrack<54 ) fakerate = 0.000250924611464;
            else if( jet_nTrack>=54 && jet_nTrack<55 ) fakerate = 4.26125334343e-05;
            else if( jet_nTrack>=55 && jet_nTrack<56 ) fakerate = 5.47255185666e-05;
            else if( jet_nTrack>=56 && jet_nTrack<57 ) fakerate = 7.0546499046e-05;
            else if( jet_nTrack>=57 && jet_nTrack<58 ) fakerate = 0.000768716272432;
            else if( jet_nTrack>=58 && jet_nTrack<59 ) fakerate = 0.0;
            else if( jet_nTrack>=59 && jet_nTrack<60 ) fakerate = 0.000177867972525;
        }
        else if( jet_pt>=700 ){
            if( jet_nTrack>=0 && jet_nTrack<1 ) fakerate = 0.0;
            else if( jet_nTrack>=1 && jet_nTrack<2 ) fakerate = 0.0447863936424;
            else if( jet_nTrack>=2 && jet_nTrack<3 ) fakerate = 0.0588819645345;
            else if( jet_nTrack>=3 && jet_nTrack<4 ) fakerate = 0.0513378120959;
            else if( jet_nTrack>=4 && jet_nTrack<5 ) fakerate = 0.0438645370305;
            else if( jet_nTrack>=5 && jet_nTrack<6 ) fakerate = 0.0386700183153;
            else if( jet_nTrack>=6 && jet_nTrack<7 ) fakerate = 0.0364851355553;
            else if( jet_nTrack>=7 && jet_nTrack<8 ) fakerate = 0.0314921587706;
            else if( jet_nTrack>=8 && jet_nTrack<9 ) fakerate = 0.0285205300897;
            else if( jet_nTrack>=9 && jet_nTrack<10 ) fakerate = 0.0260573718697;
            else if( jet_nTrack>=10 && jet_nTrack<11 ) fakerate = 0.0235692001879;
            else if( jet_nTrack>=11 && jet_nTrack<12 ) fakerate = 0.02212530002;
            else if( jet_nTrack>=12 && jet_nTrack<13 ) fakerate = 0.0207715183496;
            else if( jet_nTrack>=13 && jet_nTrack<14 ) fakerate = 0.0190450642258;
            else if( jet_nTrack>=14 && jet_nTrack<15 ) fakerate = 0.0174615960568;
            else if( jet_nTrack>=15 && jet_nTrack<16 ) fakerate = 0.0161285307258;
            else if( jet_nTrack>=16 && jet_nTrack<17 ) fakerate = 0.0139742586762;
            else if( jet_nTrack>=17 && jet_nTrack<18 ) fakerate = 0.0134326573461;
            else if( jet_nTrack>=18 && jet_nTrack<19 ) fakerate = 0.0118568539619;
            else if( jet_nTrack>=19 && jet_nTrack<20 ) fakerate = 0.0109351547435;
            else if( jet_nTrack>=20 && jet_nTrack<21 ) fakerate = 0.00968848727643;
            else if( jet_nTrack>=21 && jet_nTrack<22 ) fakerate = 0.00895960722119;
            else if( jet_nTrack>=22 && jet_nTrack<23 ) fakerate = 0.0073858378455;
            else if( jet_nTrack>=23 && jet_nTrack<24 ) fakerate = 0.00651353038847;
            else if( jet_nTrack>=24 && jet_nTrack<25 ) fakerate = 0.00586835388094;
            else if( jet_nTrack>=25 && jet_nTrack<26 ) fakerate = 0.00529484637082;
            else if( jet_nTrack>=26 && jet_nTrack<27 ) fakerate = 0.00424272380769;
            else if( jet_nTrack>=27 && jet_nTrack<28 ) fakerate = 0.00410476280376;
            else if( jet_nTrack>=28 && jet_nTrack<29 ) fakerate = 0.00343705480918;
            else if( jet_nTrack>=29 && jet_nTrack<30 ) fakerate = 0.00304065016098;
            else if( jet_nTrack>=30 && jet_nTrack<31 ) fakerate = 0.0027418883983;
            else if( jet_nTrack>=31 && jet_nTrack<32 ) fakerate = 0.00238865474239;
            else if( jet_nTrack>=32 && jet_nTrack<33 ) fakerate = 0.00205157999881;
            else if( jet_nTrack>=33 && jet_nTrack<34 ) fakerate = 0.00169514503796;
            else if( jet_nTrack>=34 && jet_nTrack<35 ) fakerate = 0.00164767552633;
            else if( jet_nTrack>=35 && jet_nTrack<36 ) fakerate = 0.00149726239033;
            else if( jet_nTrack>=36 && jet_nTrack<37 ) fakerate = 0.00120591977611;
            else if( jet_nTrack>=37 && jet_nTrack<38 ) fakerate = 0.00144696678035;
            else if( jet_nTrack>=38 && jet_nTrack<39 ) fakerate = 0.00106767576654;
            else if( jet_nTrack>=39 && jet_nTrack<40 ) fakerate = 0.000990787753835;
            else if( jet_nTrack>=40 && jet_nTrack<41 ) fakerate = 0.000736114860047;
            else if( jet_nTrack>=41 && jet_nTrack<42 ) fakerate = 0.00078924826812;
            else if( jet_nTrack>=42 && jet_nTrack<43 ) fakerate = 0.000829569704365;
            else if( jet_nTrack>=43 && jet_nTrack<44 ) fakerate = 0.000609877693933;
            else if( jet_nTrack>=44 && jet_nTrack<45 ) fakerate = 0.000497099827044;
            else if( jet_nTrack>=45 && jet_nTrack<46 ) fakerate = 0.000568752584513;
            else if( jet_nTrack>=46 && jet_nTrack<47 ) fakerate = 0.000466152909212;
            else if( jet_nTrack>=47 && jet_nTrack<48 ) fakerate = 0.000465505669126;
            else if( jet_nTrack>=48 && jet_nTrack<49 ) fakerate = 0.000449804530945;
            else if( jet_nTrack>=49 && jet_nTrack<50 ) fakerate = 0.000436621834524;
            else if( jet_nTrack>=50 && jet_nTrack<51 ) fakerate = 0.000426853308454;
            else if( jet_nTrack>=51 && jet_nTrack<52 ) fakerate = 0.000510422629304;
            else if( jet_nTrack>=52 && jet_nTrack<53 ) fakerate = 0.000362708873581;
            else if( jet_nTrack>=53 && jet_nTrack<54 ) fakerate = 0.000419475400122;
            else if( jet_nTrack>=54 && jet_nTrack<55 ) fakerate = 0.000385161780287;
            else if( jet_nTrack>=55 && jet_nTrack<56 ) fakerate = 0.000453403539723;
            else if( jet_nTrack>=56 && jet_nTrack<57 ) fakerate = 0.000172137253685;
            else if( jet_nTrack>=57 && jet_nTrack<58 ) fakerate = 4.65822886326e-05;
            else if( jet_nTrack>=58 && jet_nTrack<59 ) fakerate = 0.000253831822192;
            else if( jet_nTrack>=59 && jet_nTrack<60 ) fakerate = 0.000317007652484;
            else if( jet_nTrack>=60 && jet_nTrack<61 ) fakerate = 0.000221586102271;
            else if( jet_nTrack>=61 && jet_nTrack<62 ) fakerate = 0.000322647509165;
            else if( jet_nTrack>=62 && jet_nTrack<63 ) fakerate = 0.0;
            else if( jet_nTrack>=63 && jet_nTrack<64 ) fakerate = 0.0;
            else if( jet_nTrack>=64 && jet_nTrack<65 ) fakerate = 0.0;
            else if( jet_nTrack>=65 && jet_nTrack<66 ) fakerate = 0.0;
            else if( jet_nTrack>=66 && jet_nTrack<67 ) fakerate = 0.000709799816832;
            else if( jet_nTrack>=67 && jet_nTrack<68 ) fakerate = 0.00042257248424;
            else if( jet_nTrack>=68 && jet_nTrack<69 ) fakerate = 0.00153069198132;
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

double frWeightFT0(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType){
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
        p0 *= (1.0 - fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks,varType,flav));
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

double frWeightFT1(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType){
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
        double p11 = fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1);
        for(Int_t k=0; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            if (kdx!=jdx) p11 *= (1.0 - fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2));
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


//Ntag==22
double frWeightFT2(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType){
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
        //std::cout << "[frWeightFT2] jet flavor = " << flav1 << std::endl;
        double p11=fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1);
        for(Int_t k=0; k<njets; k++) {
            if (k==j) continue;
            double p12 = p11;
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            double p22 = p12*fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2);
            for(Int_t l=0; l<njets; l++) {
                if (l==j || l==k) continue;
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                int flav3 = flavor[ldx];
                p22 *= (1.0-fakerateF(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType,flav3));
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
                    //p33 *= (1.0-fakerate(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType));
                }
                p3 += p33;
            }
        }
    }
    return p3;
}

double frWeightFT3(vector<float> *jetpt, vector<float> *jeteta, vector<int> &goodjetIdx, vector<int> &ntrack, vector<int> &flavor, int njetscut, double jptcut, int varType){
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
        double p31 = fakerateF(jetpt->at(jdx),jeteta->at(jdx),ntrks1,varType,flav1);
        for(Int_t k=j+1; k<njets; k++) {
            int kdx = goodjetIdx[k];
            int ntrks2 = ntrack[kdx];
            int flav2 = flavor[kdx];
            double p32 = p31 * fakerateF(jetpt->at(kdx),jeteta->at(kdx),ntrks2,varType,flav2);
            for(Int_t l=k+1; l<njets; l++) {
                int ldx = goodjetIdx[l];
                int ntrks3 = ntrack[ldx];
                int flav3 = flavor[ldx];
                double p33 = p32 * fakerateF(jetpt->at(ldx),jeteta->at(ldx),ntrks3,varType,flav3);
                for(Int_t m=0; m<njets; m++) {
                    if (m==j || m==k || m ==l) continue;
                    int mdx = goodjetIdx[m];
                    int ntrks4 = ntrack[mdx];
                    int flav4 = flavor[mdx];
                    p33 *= (1.0-fakerateF(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType,flav4));
                    //p33 *= (1.0-fakerate(jetpt->at(mdx),jeteta->at(mdx),ntrks4,varType));
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
        if ( track_ipXYSigs.at(itk) < 4.0 ) ptsum += track_pt.at(itk);
    }

    double alpha = (ptsum_total > 0 ? ptsum/ptsum_total : 0.);
    if (verbose) std::cout << std::fixed << std::setprecision(6) << "[GetAlpha2Dsig] alpha = " << alpha << std::endl;
    return alpha;
}

// default for analysis_20170523_v0
double GetAlpha2Dsig(vector<float> &track_pt, vector<int> &track_source, vector<int> &track_quality,
                     vector<float> &track_ipXYSigs, vector<float> &track_ref_zs, float pv_z, float pilecut)
{
    double ptsum_total=0, ptsum=0;
    for (unsigned itk=0; itk < track_pt.size(); itk++) {
        if ( track_source.at(itk) != 0 ) continue; // Only process tracks with source=0
        if ( (track_quality.at(itk) & 4 ) == 0 ) continue; // Only process tracks with "highPurity" quality
        if (fabs(pv_z-track_ref_zs.at(itk))>pilecut) continue;// remove tracks with exceedingly large z
        ptsum_total += track_pt.at(itk);
        if ( track_ipXYSigs.at(itk) < 4.0 ) ptsum += track_pt.at(itk);
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
