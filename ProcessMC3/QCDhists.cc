#include <TList.h>
#include "QCDhists.h"
#include "EMJselect.h"
#include "EMJscan.h"
#include "EMJ16003.h"


const bool mergeOnly = false;
bool scanCuts = false;

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
                    itmp = EMJselect(true,hasPre,inputfile.c_str(),outputfile.c_str(),DHTcut, Dpt1cut,Dpt2cut,Dpt3cut,Dpt4cut,Djetacut,Dalphacut,DmaxIPcut,0.9,0.9,Dntrk1,Dnemcut,blind);
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
    const int nhist=98;
    std::vector<TH1F*> vv(nhist);
    std::string histnames[nhist]={
        "count","acount",
        "hjetcut","hjetchf","h_nemg","h_nalemg",
        "hnjet","hpt","heta","heta2",
        "H_T","H_T0","H_T1","H_T2","H_T3","H_T4",
        "hpt1","hpt2","hpt3",
        "hpt4","hbcut_ntrkpt1","hacut_ntrkpt1","hbcut_nef","hacut_nef",
        "hbcut_cef","hacut_cef","hbcut_alphamax","hacut_alphamax","hHTnm1",
        "hpt1nm1","hpt2nm1","hpt3nm1","hpt4nm1","halphanm1",
        "hmaxipnm1","hnHitsnm1","hntrk1nm1","hnemnm1","hipXYEJ",
        "hipXYnEJ","htvwEJ","htvw","hipXYSigEJ","hipXYSignEJ",
        "hmaxipXYEJ","hmaxipXYnEJ","hmeanipXYEJ","hmeanipXYnEJ","hnmaxipnm1",
        "hn2maxipnm1","hjptfrb","hjptfra1",
        "hjptfra2","hjptfrbc","hjptfra1c","hjptfra2c","hjptb",
        "hjpta","haMgj","hHTko","hpt1ko","hpt2ko",
        "hpt3ko","hpt4ko","hmass","hlogmedipXYSigEJ","hlogmedipXYSignEJ","hlogmeanipXYSigEJ","hlogmeanipXYSignEJ",
        "hmedipXYSigEJ","hmedipXYSignEJ","hmeanipXYSigEJ","hmeanipXYSignEJ","hmedipXYEJ","hmedipXYnEJ",
        "hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n","h_ntag",
        "halpha","halphaPS",
        "hmedtheta2DEJ","hmedtheta2DnEJ","hlogmedtheta2DEJ","hlogmedtheta2DnEJ",
        "hmedtheta2DPS","hlogmedtheta2DPS","hmedipXYSigPS","hlogmedipXYSigPS",
        "hmedtheta2DSR","hlogmedtheta2DSR","hmedipXYSigSR","hlogmedipXYSigSR",
        "hfr_ntrkpt1d","hfr_ntrkpt1n",
    };
    vector<double> outnorm(nbin);
    for(int i=0;i<nhist;i++) {
        std::cout<<" entering Histman with i = "<<i<<": "<<histnames[i]<<std::endl;
        vv[i]=HistMan(goalintlum,histnames[i],norm,outnorm,nbin,xsec,nfiles,binnames,donorm,bbname);
    }

    const int nhist2=11;
    std::vector<TH2F*> vv2(nhist2);
    std::string histnames2[nhist2]={
        "aMip","haMvjpt","haMvHT","haMvnvtx",
        "halphavtheta2D","halphavipXYSig","htheta2DvipXYSig",
        "halphavtheta2DPS","halphavipXYSigPS","htheta2DvipXYSigPS",
        "htheta2DvipXYSigSR",
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
