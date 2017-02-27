#include <iostream>
#include <string>
#include <map>
#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "MergeHists.h"


int main(int argc, char *argv[])
{
    int dooptk =*(argv[1])-'0';
    int doopta =*(argv[2])-'0';
    int imode=*(argv[3])-'0';
    int iblind=*(argv[4])-'0';
    int i16003=*(argv[5])-'0';

    std::string outdir;
    std::cout << "# of argvs = " << argc << std::endl;
    if (argv[6]!=NULL){
        outdir = std::string(argv[6])+"/";
    }
    else outdir = "histos/test/";

    int pmode = 0;
    if (argv[7]!=NULL){
        std::cout << "argv[7] = " << argv[7] << std::endl;
        pmode = (int)atoi(argv[7]);
    }
    std::cout << "pmode = " << pmode << std::endl;

    int nrange[2];
    if (argv[8]!=NULL && argv[9]!=NULL){
        nrange[0] = (int)atoi(argv[8]);
        nrange[1] = (int)atoi(argv[9]);
    }

    std::cout << "Results will be saved to: " << outdir << std::endl;

    std::cout << "imode = " << imode << std::endl;
    bool b16003 = false;
    if(i16003>0) b16003=true;

    bool blind=false;
    if(iblind!=0) blind=true;

    bool hasPre=true;

    if(dooptk==0) {
        std::cout<<"not doing kinematic optimization"<<std::endl;
    } else if(dooptk==1) {
        std::cout<<"doing kinematic optimization"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }

    if(doopta==0) {
        std::cout<<"not doing optimization on an emerging jet variable"<<std::endl;
    } else if(doopta==1) {
        std::cout<<"doing alphamax optimization"<<std::endl;
    } else if(doopta==2) {
        std::cout<<"doing maxIP optimization"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }

    if(imode==0) {
        std::cout<<"doing background"<<std::endl;
    } else if(imode==1) {
        std::cout<<"doing modelA"<<std::endl;
    } else if(imode==2) {
        std::cout<<"doing modelB"<<std::endl;
    } else if(imode==3) {
        std::cout<<"doing QCD74HT500to700"<<std::endl;
    } else if(imode==4) {
        std::cout<<"doing QCD74HT700to1000"<<std::endl;
    } else if(imode==5) {
        std::cout<<"doing QCD74HT1000to1500"<<std::endl;
    } else if(imode==6) {
        std::cout<<"doing QCD74HT1500toInf"<<std::endl;
    } else if(imode==7) {
        std::cout<<"doing DATA"<<std::endl;
        hasPre=true;
        blind=true;
    } else if(imode==8) {
        std::cout<<"doing QCD74"<<std::endl;
    } else if(imode==9) {
        std::cout<<"doing RelVal"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }


    float goalintlum=20; // fb-1                                                                                        
    //float goalintlum=0.07956; // fb-1                                                                                        

    std::string aaname;
    if (imode==9) aaname = "/data/users/jengbou/EmJetMC/";  // RelVal samples
    else if (imode<=7) aaname = "/data/users/jengbou/crab_output/ntuple_20170223_v0/";
    else aaname = "/data/users/eno/outputQCD/";

    std::cout << "Input directory = " << aaname << std::endl;
    const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int nfiles[nbin]={138,133,50,40,23};
    std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};

    // quick background
    const int qnbin=1; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float qxsec[qnbin]={25420}; // fb 
    int qnfiles[qnbin]={10};
    std::string qbinnames[qnbin]={"QCD_HT2000toInf"};


    // for signal models A.  mediat mass is 1000
    const int anbin=1; 
    float axsec[anbin]={18.45}; // fb 
    //int anfiles[anbin]={716};//74Full
    int anfiles[anbin]={71}; 
    std::string abinnames[anbin]={"modelA"};

    // for signal models B.  mediat mass is 1000
    const int bnbin=1; 
    float bxsec[bnbin]={18.45}; // fb 
    //int bnfiles[bnbin]={498}; //74Full
    int bnfiles[bnbin]={50};
    std::string bbinnames[bnbin]={"modelB"};


    // DATA
    const int datanbin=1; 
    float dataxsec[datanbin]={11811000}; // fb 
    int datanfiles[datanbin]={19};
    std::string databinnames[datanbin]={"DATA"};


    //QCD74
    const int q74nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int q74nfiles[q74nbin]={183,153,65,44,26};
    std::string q74binnames[q74nbin]={"QCD74_HT500to700","QCD74_HT700to1000","QCD74_HT1000to1500","QCD74_HT1500to2000","QCD74_HT2000toInf"};


    // RelVal ModelA
    // for signal models A.  mediat mass is 1000
    const int vnbin=1; 
    float vxsec[vnbin]={18.45}; // fb
    int vnfiles[vnbin]={10};//76
    //std::string vbinnames[anbin]={"Tests"};
    std::string vbinnames[anbin]={"ModelA74X_20170215_v0"};

    // mode = 3
    //QCD74
    const int q74nbin1=1; // 500-700
    float q74xsec1[q74nbin1]={29370000}; //fb
    int q74nfiles1[q74nbin1]={1843};//1843
    int q74nfrng1[q74nbin1][2]={1,1843};
    std::string q74binnames1[q74nbin1]={"QCD74_HT500to700"};

    // mode = 4
    const int q74nbin2=1; // 700-1000
    float q74xsec2[q74nbin2]={6524000}; // fb
    int q74nfiles2[q74nbin2]={1536};
    int q74nfrng2[q74nbin2][2]={1,1536};
    std::string q74binnames2[q74nbin2]={"QCD74_HT700to1000"};

    // mode = 5
    const int q74nbin3=1; // 1000-1500
    float q74xsec3[q74nbin3]={6524000}; // fb
    int q74nfiles3[q74nbin3]={651};
    int q74nfrng3[q74nbin3][2]={1,651};
    std::string q74binnames3[q74nbin3]={"QCD74_HT1000to1500"};

    // mode = 6
    const int q74nbin4=1; // 1500-2000
    float q74xsec4[q74nbin4]={121500}; // fb
    int q74nfiles4[q74nbin4]={452};
    int q74nfrng4[q74nbin4][2]={1,1452};
    std::string q74binnames4[q74nbin4]={"QCD74_HT1500to2000"};

    // mode = 7
    const int q74nbin5=1; // 200toInf
    float q74xsec5[q74nbin5]={25420}; // fb
    int q74nfiles5[q74nbin5]={261};
    int q74nfrng5[q74nbin5][2]={1,261};
    std::string q74binnames5[q74nbin5]={"QCD74_HT2000toInf"};


    if (pmode==0) {
        std::cout << "All in one go." << std::endl;
        if(imode==0) {
            QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==1) {
            QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,aaname,"SumHistsModelA.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==2) {
            QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname,"SumHistsModelB.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==3) {//QCD74 HT500to700
            QCDhists(goalintlum,1,q74xsec1,q74nfiles1,q74binnames1,aaname,"SumHistsQCD74_HT500to700.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==4) {// QCD74 HT700to1000
            QCDhists(goalintlum,1,q74xsec2,q74nfiles2,q74binnames2,aaname,"SumHistsQCD74_HT700to1000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==5) {// QCD74 HT1000to1500
            QCDhists(goalintlum,1,q74xsec3,q74nfiles3,q74binnames3,aaname,"SumHistsQCD74_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==6) {// QCD74 HT1500to2000
            QCDhists(goalintlum,1,q74xsec4,q74nfiles4,q74binnames4,aaname,"SumHistsQCD74_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==7) {// QCD74 HT2000toInf
            QCDhists(goalintlum,1,q74xsec5,q74nfiles5,q74binnames5,aaname,"SumHistsQCD74_HT200toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==8) {
            QCDhists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,aaname,"SumHistsDATA.root",0,0,hasPre,false,blind,b16003,outdir,true);
        } else if (imode==9) {
            QCDhists(goalintlum,vnbin,vxsec,vnfiles,vbinnames,aaname+"RelVal/","SumHistsModelARelValEXO16003.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        }
    }
    else if (pmode==1) {
        std::cout << "No merging." << std::endl;
        if (imode==3) {//QCD74 HT500to700
            QCDhistsNoMerge(nrange,"QCD74_HT500to700",aaname,"PartialSumHistsQCD74_HT500to700.root",hasPre,blind,b16003,outdir,true);
        } else if (imode==4) {// QCD74 HT700to1000
            QCDhistsNoMerge(nrange,"QCD74_HT700to1000",aaname,"PartialSumHistsQCD74_HT700to1000.root",hasPre,blind,b16003,outdir,true);
        } else if (imode==5) {// QCD74 HT1000to1500
            QCDhistsNoMerge(nrange,"QCD74_HT1000to1500",aaname,"PartialSumHistsQCD74_HT1000to1500.root",hasPre,blind,b16003,outdir,true);
        } else if (imode==6) {// QCD74 HT1500to2000
            QCDhistsNoMerge(nrange,"QCD74_HT1500to2000",aaname,"PartialSumHistsQCD74_HT1500to2000.root",hasPre,blind,b16003,outdir,true);
        } else if (imode==7) {// QCD74 HT2000toInf
            QCDhistsNoMerge(nrange,"QCD74_HT2000toInf",aaname,"PartialSumHistsQCD74_HT2000toInf.root",hasPre,blind,b16003,outdir,true);
        }
    }
    else{
        if (imode==3) {//QCD74 HT500to700
            MergeHists(goalintlum,q74nbin1,q74xsec1,q74nfiles1,q74binnames1,aaname,"SumHistsQCD74_HT500to700.root",true,outdir);
        } else if (imode==4) {// QCD74 HT700to1000
            MergeHists(goalintlum,q74nbin2,q74xsec2,q74nfiles2,q74binnames2,aaname,"SumHistsQCD74_HT700to1000.root",true,outdir);
        } else if (imode==5) {// QCD74 HT1000to1500
            MergeHists(goalintlum,q74nbin3,q74xsec3,q74nfiles3,q74binnames3,aaname,"SumHistsQCD74_HT1000to1500.root",true,outdir);
        } else if (imode==6) {// QCD74 HT1500to2000
            MergeHists(goalintlum,q74nbin4,q74xsec4,q74nfiles4,q74binnames4,aaname,"SumHistsQCD74_HT1500to2000.root",true,outdir);
        } else if (imode==7) {// QCD74 HT2000toInf
            MergeHists(goalintlum,q74nbin5,q74xsec5,q74nfiles5,q74binnames5,aaname,"SumHistsQCD74_HT2000toInf.root",true,outdir);
        }
    }

}

