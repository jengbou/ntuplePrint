#include <iostream>
#include <string>
#include <map>
#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "QCDhistsNoMergeNew.h"
#include "MergeHists.h"
#include "MergeHistsNoNorm.h"
#include "MergeHistsNoNormNew.h"
#include "MergeHistsNorm.h"
#include "MergeHistsN.h"
#include "MergeHists_EXO16003.h"
#include "MergeHistsNoNorm_EXO16003.h"
#include "MergeHistsN_EXO16003.h"

int main(int argc, char *argv[])
{
    int dooptk =*(argv[1])-'0';
    int doopta =*(argv[2])-'0';
    int imode=(int)atoi((argv[3]));
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
    int fidx;
    if (argv[10]!=NULL){
        fidx = (int)atoi(argv[10]);
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
        std::cout<<"doing QCD74HT1500to2000"<<std::endl;
    } else if(imode==7) {
        std::cout<<"doing QCD74HT2000toInf"<<std::endl;
    } else if(imode==8) {
        std::cout<<"doing QCD74"<<std::endl;
    } else if(imode==9) {
        std::cout<<"doing RelVal"<<std::endl;
    } else if(imode==10 || imode==100) {
        std::cout<<"doing DATA; blinding for the time being."<<std::endl;
        hasPre=true;
        blind=true;
    } else if(imode==11) {
        std::cout<<"doing QCD80HT100to200"<<std::endl;
    } else if(imode==12) {
        std::cout<<"doing QCD80HT200to300"<<std::endl;
    } else if(imode==13) {
        std::cout<<"doing QCD80HT300to500"<<std::endl;
    } else if(imode==14) {
        std::cout<<"doing QCD80HT500to700"<<std::endl;
    } else if(imode==15) {
        std::cout<<"doing QCD80HT700to1000"<<std::endl;
    } else if(imode==16) {
        std::cout<<"doing QCD80HT1000to1500"<<std::endl;
    } else if(imode==17) {
        std::cout<<"doing QCD80HT1500to2000"<<std::endl;
    } else if(imode==18) {
        std::cout<<"doing QCD80HT2000toInf"<<std::endl;
    } else if(imode==28) {
        std::cout<<"doing QCD94HT2000toInf"<<std::endl;
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }


    float goalintlum=16.132397; // fb-1 G-H 2016
    float goalintlum17=14.487;// fb-1 B+C ; 41.86 for 2017 A-F
    //float goalintlum=20.0; // fb-1
    //float goalintlum=0.07956; // fb-1

    std::string indir;
    if (imode==9) indir = "/data/users/jengbou/EmJetMC/";  // RelVal samples
    //else if (imode<=7) indir = "/data/users/jengbou/crab_output/ntuple_20170223_v0/";
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170223_v0/";//QCD76
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170301_v0/";//QCD76
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170420_v0/";
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170427_v0/";//QCD76 w/o Sarah's alphaMax calculation
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170427_v1/";//QCD76 w/ Sarah's alphaMax calculation
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170426_v2/";//QCD74&80 YH's (ntuple_x.root; use pv_indexInColl)
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170517_v0/";//QCD80 HT700+
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170523_v0/";//80 YH's
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170823_v0/";//80 YH's
    //else indir = "/data/users/jengbou/crab_output/ntuple_20171206_v0/";//QCD: 0823_v0; Signals: 1002_v2: Data: 1103_v0
    //else indir = "/data/users/jengbou/crab_output/ntuple_20180126_v0/";//QCD: 1103_v0 and 1103_v1(700to1000x); Signals: 1103_v1(b): Data: 1103_v0
    //else indir = "/data/users/jengbou/crab_output/ntuple_20180208_v0/";//same as 20180126 except original files are in /data rather than in /store
    else indir = "/data/users/jengbou/crab_output/ntuple_20181019_v1/";//2017 QCD
    //else indir = "/home/jengbou/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/temp/test_20180427/";

    // Note: QCD74vNpM are in condor format, e.g., QCD74_HT1000to1500_X_0.histo.root so remember to change QCDhistsNoMerge arguments
    // It uses pv_index so remember to change EMJbkg.cc
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170502_v0_QCD74v3p2/";//QCD74 w/ Sarah's alphaMax calculation
    //else indir = "/data/users/jengbou/crab_output/ntuple_20170505_v0_QCD74v3p3/";//QCD74

    //else indir = "/data/users/eno/outputQCD/";

    std::cout << "Input directory = " << indir << std::endl;
    const int nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float xsec[nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int nfiles[nbin]={138,133,50,40,23};
    std::string binnames[nbin]={"QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf"};

    // quick background
    const int qnbin=1; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float qxsec[qnbin]={25420}; // fb 
    int qnfiles[qnbin]={10};
    std::string qbinnames[qnbin]={"QCD_HT2000toInf"};


    // mode = 10
    // DATA
    const int datanbin=1; 
    float dataxsec[datanbin]={11811000}; // fb 
    int datanfiles[datanbin]={6748};
    std::string databinnames[datanbin]={"Data"};


    //QCD74
    const int q74nbin=5; // 500-700,700-1000,1000-1500,1500-2000,200toInf
    float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
    int q74nfiles[q74nbin]={183,153,65,44,26};
    std::string q74binnames[q74nbin]={"QCD74_HT500to700","QCD74_HT700to1000","QCD74_HT1000to1500","QCD74_HT1500to2000","QCD74_HT2000toInf"};

    // mode = 9
    // RelVal ModelA
    // for signal models A.  mediat mass is 1000
    const int vnbin=1; 
    float vxsec[vnbin]={18.45}; // fb
    int vnfiles[vnbin]={10};//76
    //std::string vbinnames[vnbin]={"Tests"};
    std::string vbinnames[vnbin]={"ModelA74X_20170215_v0"};

    // mode = 1
    // for signal models A.  mediat mass is 1000
    const int anbin=1; 
    float axsec[anbin]={18.45}; // fb 
    int anfiles[anbin]={30};//74Full 716
    //int anfiles[anbin]={71}; 
    int anfrng[anbin][2]={1,30};//716
    std::string abinnames[anbin]={"modelA"};

    // mode = 2
    // for signal models B.  mediat mass is 1000
    const int bnbin=1; 
    float bxsec[bnbin]={18.45}; // fb 
    int bnfiles[bnbin]={30}; //74Full 498
    //int bnfiles[bnbin]={50};
    int bnfrng[bnbin][2]={1,30};//498
    std::string bbinnames[bnbin]={"modelB"};

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
    float q74xsec3[q74nbin3]={1064000}; // fb
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

    // mode = 11
    //QCD80
    const int q76nbin1=1; // 100-200
    float q76xsec1[q76nbin1]={27990000000.0}; //fb
    int q76nfiles1[q76nbin1]={3951};//1843
    int q76nfrng1[q76nbin1][2]={1,3951};
    std::string q76binnames1[q76nbin1]={"QCD80_HT100to200"};

    // mode = 12
    //QCD80
    const int q76nbin2=1; // 200-300
    float q76xsec2[q76nbin2]={1712000000}; //fb
    int q76nfiles2[q76nbin2]={1262};//1843
    int q76nfrng2[q76nbin2][2]={1,1262};
    std::string q76binnames2[q76nbin2]={"QCD80_HT200to300"};

    // mode = 13
    //QCD80
    const int q76nbin3=1; // 300-500
    float q76xsec3[q76nbin3]={347700000}; //fb
    int q76nfiles3[q76nbin3]={1012};//1843
    int q76nfrng3[q76nbin3][2]={1,1012};
    std::string q76binnames3[q76nbin3]={"QCD80_HT300to500"};

    // mode = 14
    //QCD80
    const int q76nbin4=1; // 500-700
    float q76xsec4[q76nbin4]={32100000}; //fb
    int q76nfiles4[q76nbin4]={1375};//1843
    int q76nfrng4[q76nbin4][2]={1,1375};
    std::string q76binnames4[q76nbin4]={"QCD80_HT500to700"};

    // mode = 15
    //QCD80
    const int q76nbin5=1; // 700-1000
    float q76xsec5[q76nbin5]={6831000}; //fb
    int q76nfiles5[q76nbin5]={1325};//1843
    int q76nfrng5[q76nbin5][2]={1,1325};
    std::string q76binnames5[q76nbin5]={"QCD80_HT700to1000"};

    // mode = 16
    //QCD80
    const int q76nbin6=1; // 1000-1500
    float q76xsec6[q76nbin6]={1207000}; //fb
    int q76nfiles6[q76nbin6]={498};//1843
    int q76nfrng6[q76nbin6][2]={1,498};
    std::string q76binnames6[q76nbin6]={"QCD80_HT1000to1500"};

    // mode = 17
    //QCD80
    const int q76nbin7=1; // 1500-2000
    float q76xsec7[q76nbin7]={119900}; //fb
    int q76nfiles7[q76nbin7]={397};//1843
    int q76nfrng7[q76nbin7][2]={1,397};
    std::string q76binnames7[q76nbin7]={"QCD80_HT1500to2000"};

    // mode = 18
    //QCD80
    const int q76nbin8=1; // 2000-Inf
    float q76xsec8[q76nbin8]={25240}; //fb
    int q76nfiles8[q76nbin8]={223};//1843
    int q76nfrng8[q76nbin8][2]={1,223};
    std::string q76binnames8[q76nbin8]={"QCD80_HT2000toInf"};


    // 2017 MC
    // mode = 26
    //QCD94
    const int q94nbin6=1; // 1000-1500
    float q94xsec6[q94nbin6]={1088000}; //fb
    int q94nfiles6[q94nbin6]={1398};
    int q94nfrng6[q94nbin6][2]={1,1398};
    std::string q94binnames6[q94nbin6]={"QCD94_HT1000to1500"};

    // mode = 27
    //QCD94
    const int q94nbin7=1; // 1500-2000
    float q94xsec7[q94nbin7]={99110}; //fb
    int q94nfiles7[q94nbin7]={980};
    int q94nfrng7[q94nbin7][2]={1,980};
    std::string q94binnames7[q94nbin7]={"QCD94_HT1500to2000"};

    // mode = 28
    //QCD94
    const int q94nbin8=1; // 2000-Inf
    float q94xsec8[q94nbin8]={20230}; //fb
    int q94nfiles8[q94nbin8]={473};
    int q94nfrng8[q94nbin8][2]={1,473};
    std::string q94binnames8[q94nbin8]={"QCD94_HT2000toInf"};

    if (pmode==0) {
        std::cout << "All in one go." << std::endl;
        if(imode==0) {
            QCDhists(goalintlum,nbin,xsec,nfiles,binnames,indir,"SumHistsQCD.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==1) {
            QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,indir,"SumHistsModelA.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==2) {
            QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,indir,"SumHistsModelB.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==3) {// QCD74 HT500to700
            QCDhists(goalintlum,1,q74xsec1,q74nfiles1,q74binnames1,indir,"SumHistsQCD74_HT500to700.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==4) {// QCD74 HT700to1000
            QCDhists(goalintlum,1,q74xsec2,q74nfiles2,q74binnames2,indir,"SumHistsQCD74_HT700to1000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==5) {// QCD74 HT1000to1500
            QCDhists(goalintlum,1,q74xsec3,q74nfiles3,q74binnames3,indir,"SumHistsQCD74_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==6) {// QCD74 HT1500to2000
            QCDhists(goalintlum,1,q74xsec4,q74nfiles4,q74binnames4,indir,"SumHistsQCD74_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==7) {// QCD74 HT2000toInf
            QCDhists(goalintlum,1,q74xsec5,q74nfiles5,q74binnames5,indir,"SumHistsQCD74_HT200toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==8) {// QCD74
            QCDhists(goalintlum,1,q74xsec,q74nfiles,q74binnames,indir,"SumHistsQCD74_HT500toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==9) {// RelVal
            QCDhists(goalintlum,vnbin,vxsec,vnfiles,vbinnames,indir+"RelVal/","SumHistsModelARelValEXO16003.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==10) {//Data
            QCDhists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",0,0,hasPre,false,blind,b16003,outdir,true,true);
        } else if (imode==100) {//Data 2017; goalintlum is dummy as donorm = false
            QCDhists(goalintlum17,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",0,0,hasPre,false,blind,b16003,outdir,true,true);
        } else if (imode==11) {// QCD80 HT100to200
            QCDhists(goalintlum,1,q76xsec1,q76nfiles1,q76binnames1,indir,"SumHistsQCD80_HT100to200.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==12) {// QCD80 HT200to300
            QCDhists(goalintlum,1,q76xsec2,q76nfiles2,q76binnames2,indir,"SumHistsQCD80_HT200to300.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==13) {// QCD80 HT300to500
            QCDhists(goalintlum,1,q76xsec3,q76nfiles3,q76binnames3,indir,"SumHistsQCD80_HT300to500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==14) {// QCD80 HT500to700
            QCDhists(goalintlum,1,q76xsec1,q76nfiles4,q76binnames4,indir,"SumHistsQCD80_HT500to700.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==15) {// QCD80 HT700to1000
            QCDhists(goalintlum,1,q76xsec2,q76nfiles5,q76binnames5,indir,"SumHistsQCD80_HT700to1000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==16) {// QCD80 HT1000to1500
            QCDhists(goalintlum,1,q76xsec3,q76nfiles6,q76binnames6,indir,"SumHistsQCD80_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==17) {// QCD80 HT1500to2000
            QCDhists(goalintlum,1,q76xsec4,q76nfiles7,q76binnames7,indir,"SumHistsQCD80_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==18) {// QCD80 HT2000toInf
            QCDhists(goalintlum,1,q76xsec5,q76nfiles8,q76binnames8,indir,"SumHistsQCD80_HT2000toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==26) {// QCD94 HT1000to1500
            QCDhists(goalintlum17,1,q94xsec6,q94nfiles6,q94binnames6,indir,"SumHistsQCD94_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==27) {// QCD94 HT1500to2000
            QCDhists(goalintlum17,1,q94xsec7,q94nfiles7,q94binnames7,indir,"SumHistsQCD94_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        } else if (imode==28) {// QCD94 HT2000toInf
            QCDhists(goalintlum17,1,q94xsec8,q94nfiles8,q94binnames8,indir,"SumHistsQCD94_HT2000toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true,false);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==1) {
        std::cout << "No merging." << std::endl;
        if (imode==1) {// modelA
            QCDhistsNoMergeNew(nrange,"modelA",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==2) {// modelB 
            QCDhistsNoMergeNew(nrange,"modelB",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==3) {// QCD74 HT500to700
            QCDhistsNoMergeNew(nrange,"QCD74_HT500to700",indir,hasPre,blind,b16003,outdir,false,false);
        } else if (imode==4) {// QCD74 HT700to1000
            QCDhistsNoMergeNew(nrange,"QCD74_HT700to1000",indir,hasPre,blind,b16003,outdir,false,false);
        } else if (imode==5) {// QCD74 HT1000to1500
            QCDhistsNoMergeNew(nrange,"QCD74_HT1000to1500",indir,hasPre,blind,b16003,outdir,false,false);
        } else if (imode==6) {// QCD74 HT1500to2000
            QCDhistsNoMergeNew(nrange,"QCD74_HT1500to2000",indir,hasPre,blind,b16003,outdir,false,false);
        } else if (imode==7) {// QCD74 HT2000toInf
            QCDhistsNoMergeNew(nrange,"QCD74_HT2000toInf",indir,hasPre,blind,b16003,outdir,false,false);
//         } else if (imode==3) {// QCD74 HT500to700
//             QCDhistsNoMergeNew(nrange,"QCD74_HT500to700",indir,hasPre,blind,b16003,outdir,true,false);
//         } else if (imode==4) {// QCD74 HT700to1000
//             QCDhistsNoMergeNew(nrange,"QCD74_HT700to1000",indir,hasPre,blind,b16003,outdir,true,false);
//         } else if (imode==5) {// QCD74 HT1000to1500
//             QCDhistsNoMergeNew(nrange,"QCD74_HT1000to1500",indir,hasPre,blind,b16003,outdir,true,false);
//         } else if (imode==6) {// QCD74 HT1500to2000
//             QCDhistsNoMergeNew(nrange,"QCD74_HT1500to2000",indir,hasPre,blind,b16003,outdir,true,false);
//         } else if (imode==7) {// QCD74 HT2000toInf
//             QCDhistsNoMergeNew(nrange,"QCD74_HT2000toInf",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==10) {//2016
            QCDhistsNoMergeNew(nrange,"Data",indir,hasPre,blind,b16003,outdir,true,true);
        } else if (imode==100) {//2017
            QCDhistsNoMergeNew(nrange,"Data",indir,hasPre,blind,b16003,outdir,true,true,"2017");
        } else if (imode==11) {// QCD80 HT100to200
            QCDhistsNoMergeNew(nrange,"QCD80_HT100to200",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==12) {// QCD80 HT200to300
            QCDhistsNoMergeNew(nrange,"QCD80_HT200to300",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==13) {// QCD80 HT300to500
            QCDhistsNoMergeNew(nrange,"QCD80_HT300to500",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==14) {// QCD80 HT500to700
            QCDhistsNoMergeNew(nrange,"QCD80_HT500to700",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==15) {// QCD80 HT700to1000
            QCDhistsNoMergeNew(nrange,"QCD80_HT700to1000",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==16) {// QCD80 HT1000to1500
            QCDhistsNoMergeNew(nrange,"QCD80_HT1000to1500",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==17) {// QCD80 HT1500to2000
            QCDhistsNoMergeNew(nrange,"QCD80_HT1500to2000",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==18) {// QCD80 HT2000toInf
            QCDhistsNoMergeNew(nrange,"QCD80_HT2000toInf",indir,hasPre,blind,b16003,outdir,true,false);
        } else if (imode==26) {// QCD94 HT1000to1500
            QCDhistsNoMergeNew(nrange,"QCD94_HT1000to1500",indir,hasPre,blind,b16003,outdir,true,false,"2017");
        } else if (imode==27) {// QCD94 HT1500to2000
            QCDhistsNoMergeNew(nrange,"QCD94_HT1500to2000",indir,hasPre,blind,b16003,outdir,true,false,"2017");
        } else if (imode==28) {// QCD94 HT2000toInf
            QCDhistsNoMergeNew(nrange,"QCD94_HT2000toInf",indir,hasPre,blind,b16003,outdir,true,false,"2017");
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==8) {
        if (imode==11) {// QCD80 HT100to200
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT100to200",outdir);
        }
        else if (imode==12) {// QCD80 HT200to300
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT200to300",outdir);
        }
        else if (imode==13) {// QCD80 HT300to500
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT300to500",outdir);
        }
        else if (imode==14) {// QCD80 HT500to700
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT500to700",outdir);
        }
        else if (imode==15) {// QCD80 HT700to1000
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT700to1000",outdir);
        }
        else if (imode==16) {// QCD80 HT1000to1500
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT1000to1500",outdir);
        }
        else if (imode==17) {// QCD80 HT1500to2000
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT1500to2000",outdir);
        }
        else if (imode==18) {// QCD80 HT2000toInf
            MergeHistsNoNormNew(fidx,nrange,"QCD80_HT2000toInf",outdir);
        }
        else if (imode==26) {// QCD94 HT1000to1500
            MergeHistsNoNormNew(fidx,nrange,"QCD94_HT1000to1500",outdir);
        }
        else if (imode==27) {// QCD94 HT1500to2000
            MergeHistsNoNormNew(fidx,nrange,"QCD94_HT1500to2000",outdir);
        }
        else if (imode==28) {// QCD94 HT2000toInf
            MergeHistsNoNormNew(fidx,nrange,"QCD94_HT2000toInf",outdir);
        }
        else if (imode==3) {// QCD74 HT500to700
            MergeHistsNoNormNew(fidx,nrange,"QCD74_HT500to700",outdir);
        }
        else if (imode==4) {// QCD74 HT700to1000
            MergeHistsNoNormNew(fidx,nrange,"QCD74_HT700to1000",outdir);
        }
        else if (imode==5) {// QCD74 HT1000to1500
            MergeHistsNoNormNew(fidx,nrange,"QCD74_HT1000to1500",outdir);
        }
        else if (imode==6) {// QCD74 HT1500to2000
            MergeHistsNoNormNew(fidx,nrange,"QCD74_HT1500to2000",outdir);
        }
        else if (imode==7) {// QCD74 HT2000toInf
            MergeHistsNoNormNew(fidx,nrange,"QCD74_HT2000toInf",outdir);
        }
        else if (imode==10) {// Data
            MergeHistsNoNormNew(fidx,nrange,"Data",outdir);
        }
        else if (imode==100) {// Data 2017
            MergeHistsNoNormNew(fidx,nrange,"Data",outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==9) {    //float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
        if (imode==11) {// QCD80 HT100to200
            MergeHistsNorm(goalintlum,27990000000,fidx,"QCD80_HT100to200","SumHistsQCD80_HT100to200.root",true,outdir);
        }
        else if (imode==12) {// QCD80 HT200to300
            MergeHistsNorm(goalintlum,1712000000,fidx,"QCD80_HT200to300","SumHistsQCD80_HT200to300.root",true,outdir);
        }
        else if (imode==13) {// QCD80 HT300to500
            MergeHistsNorm(goalintlum,347700000,fidx,"QCD80_HT300to500","SumHistsQCD80_HT300to500.root",true,outdir);
        }
        else if (imode==14) {// QCD80 HT500to700
            //MergeHistsNorm(goalintlum,29370000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,32100000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700.root",true,outdir);
        }
        else if (imode==15) {// QCD80 HT700to1000
            //MergeHistsNorm(goalintlum,6524000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,6831000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000.root",true,outdir);
        }
        else if (imode==16) {// QCD80 HT1000to1500
            //MergeHistsNorm(goalintlum,1064000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,1207000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500.root",true,outdir);
        }
        else if (imode==17) {// QCD80 HT1500to2000
            //MergeHistsNorm(goalintlum,121500,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,119900,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000.root",true,outdir);
        }
        else if (imode==18) {// QCD80 HT2000toInf
            //MergeHistsNorm(goalintlum,25420,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,25240,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf.root",true,outdir);
        }
        else if (imode==26) {// QCD94 HT1000to1500
            MergeHistsNorm(goalintlum17,1088000,fidx,"QCD94_HT1000to1500","SumHistsQCD94_HT1000to1500.root",true,outdir);
        }
        else if (imode==27) {// QCD94 HT1500to2000
            MergeHistsNorm(goalintlum17,99110,fidx,"QCD94_HT1500to2000","SumHistsQCD94_HT1500to2000.root",true,outdir);
        }
        else if (imode==28) {// QCD94 HT2000toInf
            MergeHistsNorm(goalintlum17,20230,fidx,"QCD94_HT2000toInf","SumHistsQCD94_HT2000toInf.root",true,outdir);
        }
        else if (imode==3) {// QCD74 HT500to700
            //MergeHistsNorm(goalintlum,29370000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,32100000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700.root",true,outdir);
        }
        else if (imode==4) {// QCD74 HT700to1000
            //MergeHistsNorm(goalintlum,6524000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,6831000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000.root",true,outdir);
        }
        else if (imode==5) {// QCD74 HT1000to1500
            //MergeHistsNorm(goalintlum,1064000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,1207000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500.root",true,outdir);
        }
        else if (imode==6) {// QCD74 HT1500to2000
            //MergeHistsNorm(goalintlum,121500,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,119900,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000.root",true,outdir);
        }
        else if (imode==7) {// QCD74 HT2000toInf
            //MergeHistsNorm(goalintlum,25420,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsNorm(goalintlum,25240,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf.root",true,outdir);
        }
        else if (imode==10) {// Data; goalintlum is dummy as donorm = false
            MergeHistsNorm(goalintlum,dataxsec[0],fidx,"Data","SumHistsData.root",false,outdir);
        }
        else if (imode==100) {// Data 2017; goalintlum17 is dummy as donorm = false
            MergeHistsNorm(goalintlum17,dataxsec[0],fidx,"Data","SumHistsData.root",false,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==10) {
        if (imode==11) {// QCD80 HT100to200
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT100to200",outdir);
        }
        else if (imode==12) {// QCD80 HT200to300
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT200to300",outdir);
        }
        else if (imode==13) {// QCD80 HT300to500
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT300to500",outdir);
        }
        else if (imode==14) {// QCD80 HT500to700
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT500to700",outdir);
        }
        else if (imode==15) {// QCD80 HT700to1000
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT700to1000",outdir);
        }
        else if (imode==16) {// QCD80 HT1000to1500
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT1000to1500",outdir);
        }
        else if (imode==17) {// QCD80 HT1500to2000
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT1500to2000",outdir);
        }
        else if (imode==18) {// QCD80 HT2000toInf
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD80_HT2000toInf",outdir);
        }
        else if (imode==26) {// QCD94 HT1000to1500
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD94_HT1000to1500",outdir);
        }
        else if (imode==27) {// QCD94 HT1500to2000
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD94_HT1500to2000",outdir);
        }
        else if (imode==28) {// QCD94 HT2000toInf
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD94_HT2000toInf",outdir);
        }
        else if (imode==3) {// QCD74 HT500to700
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD74_HT500to700",outdir);
        }
        else if (imode==4) {// QCD74 HT700to1000
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD74_HT700to1000",outdir);
        }
        else if (imode==5) {// QCD74 HT1000to1500
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD74_HT1000to1500",outdir);
        }
        else if (imode==6) {// QCD74 HT1500to2000
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD74_HT1500to2000",outdir);
        }
        else if (imode==7) {// QCD74 HT2000toInf
            MergeHistsNoNorm_EXO16003(fidx,nrange,"QCD74_HT2000toInf",outdir);
        }
        else if (imode==10) {// Data
            MergeHistsNoNorm_EXO16003(fidx,nrange,"Data",outdir);
        }
        else if (imode==100) {// Data 2017
            MergeHistsNoNorm_EXO16003(fidx,nrange,"Data",outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==11) {    //float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
        if (imode==11) {// QCD80 HT100to200
            MergeHistsN_EXO16003(goalintlum,27990000000,fidx,"QCD80_HT100to200","SumHistsQCD80_HT100to200.root",true,outdir);
        }
        else if (imode==12) {// QCD80 HT200to300
            MergeHistsN_EXO16003(goalintlum,1712000000,fidx,"QCD80_HT200to300","SumHistsQCD80_HT200to300.root",true,outdir);
        }
        else if (imode==13) {// QCD80 HT300to500
            MergeHistsN_EXO16003(goalintlum,347700000,fidx,"QCD80_HT300to500","SumHistsQCD80_HT300to500.root",true,outdir);
        }
        else if (imode==14) {// QCD80 HT500to700
            //MergeHistsN_EXO16003(goalintlum,29370000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,32100000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700.root",true,outdir);
        }
        else if (imode==15) {// QCD80 HT700to1000
            //MergeHistsN_EXO16003(goalintlum,6524000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,6831000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000.root",true,outdir);
        }
        else if (imode==16) {// QCD80 HT1000to1500
            //MergeHistsN_EXO16003(goalintlum,1064000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,1207000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500.root",true,outdir);
        }
        else if (imode==17) {// QCD80 HT1500to2000
            //MergeHistsN_EXO16003(goalintlum,121500,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,119900,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000.root",true,outdir);
        }
        else if (imode==18) {// QCD80 HT2000toInf
            //MergeHistsN_EXO16003(goalintlum,25420,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,25240,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf.root",true,outdir);
        }
        else if (imode==26) {// QCD94 HT1000to1500
            MergeHistsN_EXO16003(goalintlum17,1088000,fidx,"QCD94_HT1000to1500","SumHistsQCD94_HT1000to1500.root",true,outdir);
        }
        else if (imode==27) {// QCD94 HT1500to2000
            MergeHistsN_EXO16003(goalintlum17,99110,fidx,"QCD94_HT1500to2000","SumHistsQCD94_HT1500to2000.root",true,outdir);
        }
        else if (imode==28) {// QCD94 HT2000toInf
            MergeHistsN_EXO16003(goalintlum17,20230,fidx,"QCD94_HT2000toInf","SumHistsQCD94_HT2000toInf.root",true,outdir);
        }
        else if (imode==3) {// QCD74 HT500to700
            //MergeHistsN_EXO16003(goalintlum,29370000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,32100000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700.root",true,outdir);
        }
        else if (imode==4) {// QCD74 HT700to1000
            //MergeHistsN_EXO16003(goalintlum,6524000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,6831000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000.root",true,outdir);
        }
        else if (imode==5) {// QCD74 HT1000to1500
            //MergeHistsN_EXO16003(goalintlum,1064000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,1207000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500.root",true,outdir);
        }
        else if (imode==6) {// QCD74 HT1500to2000
            //MergeHistsN_EXO16003(goalintlum,121500,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,119900,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000.root",true,outdir);
        }
        else if (imode==7) {// QCD74 HT2000toInf
            //MergeHistsN_EXO16003(goalintlum,25420,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsN_EXO16003(goalintlum,25240,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf.root",true,outdir);
        }
        else if (imode==10) {// Data; goalintlum is dummy as donorm = false
            MergeHistsN_EXO16003(goalintlum,dataxsec[0],fidx,"Data","SumHistsData.root",false,outdir);
        }
        else if (imode==100) {// Data 2017; goalintlum is dummy as donorm = false
            MergeHistsN_EXO16003(goalintlum17,dataxsec[0],fidx,"Data","SumHistsData.root",false,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==3) {//Merging only
        if (imode==1) {// modelA
            MergeHists_EXO16003(goalintlum,anbin,axsec,anfiles,abinnames,indir,"SumHistsModelA.root",true,outdir);
        } else if (imode==2) {// modelB 
            MergeHists_EXO16003(goalintlum,bnbin,bxsec,bnfiles,bbinnames,indir,"SumHistsModelB.root",true,outdir);
        } else if (imode==3) {// QCD74 HT500to700
            MergeHists_EXO16003(goalintlum,q74nbin1,q74xsec1,q74nfiles1,q74binnames1,indir,"SumHistsQCD74_HT500to700.root",true,outdir);
        } else if (imode==4) {// QCD74 HT700to1000
            MergeHists_EXO16003(goalintlum,q74nbin2,q74xsec2,q74nfiles2,q74binnames2,indir,"SumHistsQCD74_HT700to1000.root",true,outdir);
        } else if (imode==5) {// QCD74 HT1000to1500
            MergeHists_EXO16003(goalintlum,q74nbin3,q74xsec3,q74nfiles3,q74binnames3,indir,"SumHistsQCD74_HT1000to1500.root",true,outdir);
        } else if (imode==6) {// QCD74 HT1500to2000
            MergeHists_EXO16003(goalintlum,q74nbin4,q74xsec4,q74nfiles4,q74binnames4,indir,"SumHistsQCD74_HT1500to2000.root",true,outdir);
        } else if (imode==7) {// QCD74 HT2000toInf
            MergeHists_EXO16003(goalintlum,q74nbin5,q74xsec5,q74nfiles5,q74binnames5,indir,"SumHistsQCD74_HT2000toInf.root",true,outdir);
        } else if (imode==10) {// Data; goalintlum is dummy as donorm = false
            MergeHists_EXO16003(goalintlum,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",false,outdir);
        } else if (imode==100) {// Data 2017; goalintlum is dummy as donorm = false
            MergeHists_EXO16003(goalintlum17,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",false,outdir);
        } else if (imode==11) {// QCD80 HT100to200
            MergeHists_EXO16003(goalintlum,q76nbin1,q76xsec1,q76nfiles1,q76binnames1,indir,"SumHistsQCD80_HT100to200.root",true,outdir);
        } else if (imode==12) {// QCD80 HT200to300
            MergeHists_EXO16003(goalintlum,q76nbin2,q76xsec2,q76nfiles2,q76binnames2,indir,"SumHistsQCD80_HT200to300.root",true,outdir);
        } else if (imode==13) {// QCD80 HT300to500
            MergeHists_EXO16003(goalintlum,q76nbin3,q76xsec3,q76nfiles3,q76binnames3,indir,"SumHistsQCD80_HT300to500.root",true,outdir);
        } else if (imode==14) {// QCD80 HT500to700
            MergeHists_EXO16003(goalintlum,q76nbin4,q76xsec4,q76nfiles4,q76binnames4,indir,"SumHistsQCD80_HT500to700.root",true,outdir);
        } else if (imode==15) {// QCD80 HT700to1000
            MergeHists_EXO16003(goalintlum,q76nbin5,q76xsec5,q76nfiles5,q76binnames5,indir,"SumHistsQCD80_HT700to1000.root",true,outdir);
        } else if (imode==16) {// QCD80 HT1000to1500
            MergeHists_EXO16003(goalintlum,q76nbin6,q76xsec6,q76nfiles6,q76binnames6,indir,"SumHistsQCD80_HT1000to1500.root",true,outdir);
        } else if (imode==17) {// QCD80 HT1500to2000
            MergeHists_EXO16003(goalintlum,q76nbin7,q76xsec7,q76nfiles7,q76binnames7,indir,"SumHistsQCD80_HT1500to2000.root",true,outdir);
        } else if (imode==18) {// QCD80 HT2000toInf
            MergeHists_EXO16003(goalintlum,q76nbin8,q76xsec8,q76nfiles8,q76binnames8,indir,"SumHistsQCD80_HT2000toInf.root",true,outdir);
        } else if (imode==26) {// QCD94 HT1000to1500
            MergeHists_EXO16003(goalintlum17,q94nbin6,q94xsec6,q94nfiles6,q94binnames6,indir,"SumHistsQCD94_HT1000to1500.root",true,outdir);
        } else if (imode==27) {// QCD94 HT2000toInf
            MergeHists_EXO16003(goalintlum17,q94nbin7,q94xsec7,q94nfiles7,q94binnames7,indir,"SumHistsQCD94_HT1500to2000.root",true,outdir);
        } else if (imode==28) {// QCD94 HT2000toInf
            MergeHists_EXO16003(goalintlum17,q94nbin8,q94xsec8,q94nfiles8,q94binnames8,indir,"SumHistsQCD94_HT2000toInf.root",true,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else{//Merging only
        if (imode==1) {// modelA
            MergeHists(goalintlum,anbin,axsec,anfiles,abinnames,indir,"SumHistsModelA.root",true,outdir);
        } else if (imode==2) {// modelB 
            MergeHists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,indir,"SumHistsModelB.root",true,outdir);
        } else if (imode==3) {// QCD74 HT500to700
            MergeHists(goalintlum,q74nbin1,q74xsec1,q74nfiles1,q74binnames1,indir,"SumHistsQCD74_HT500to700.root",true,outdir);
        } else if (imode==4) {// QCD74 HT700to1000
            MergeHists(goalintlum,q74nbin2,q74xsec2,q74nfiles2,q74binnames2,indir,"SumHistsQCD74_HT700to1000.root",true,outdir);
        } else if (imode==5) {// QCD74 HT1000to1500
            MergeHists(goalintlum,q74nbin3,q74xsec3,q74nfiles3,q74binnames3,indir,"SumHistsQCD74_HT1000to1500.root",true,outdir);
        } else if (imode==6) {// QCD74 HT1500to2000
            MergeHists(goalintlum,q74nbin4,q74xsec4,q74nfiles4,q74binnames4,indir,"SumHistsQCD74_HT1500to2000.root",true,outdir);
        } else if (imode==7) {// QCD74 HT2000toInf
            MergeHists(goalintlum,q74nbin5,q74xsec5,q74nfiles5,q74binnames5,indir,"SumHistsQCD74_HT2000toInf.root",true,outdir);
        } else if (imode==10) {// Data; goalintlum is dummy as donorm = false
            MergeHists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",false,outdir);
        } else if (imode==100) {// Data 2017; goalintlum is dummy as donorm = false
            MergeHists(goalintlum17,datanbin,dataxsec,datanfiles,databinnames,indir,"SumHistsData.root",false,outdir);
        } else if (imode==11) {// QCD80 HT100to200
            MergeHists(goalintlum,q76nbin1,q76xsec1,q76nfiles1,q76binnames1,indir,"SumHistsQCD80_HT100to200.root",true,outdir);
        } else if (imode==12) {// QCD80 HT200to300
            MergeHists(goalintlum,q76nbin2,q76xsec2,q76nfiles2,q76binnames2,indir,"SumHistsQCD80_HT200to300.root",true,outdir);
        } else if (imode==13) {// QCD80 HT300to500
            MergeHists(goalintlum,q76nbin3,q76xsec3,q76nfiles3,q76binnames3,indir,"SumHistsQCD80_HT300to500.root",true,outdir);
        } else if (imode==14) {// QCD80 HT500to700
            MergeHists(goalintlum,q76nbin4,q76xsec4,q76nfiles4,q76binnames4,indir,"SumHistsQCD80_HT500to700.root",true,outdir);
        } else if (imode==15) {// QCD80 HT700to1000
            MergeHists(goalintlum,q76nbin5,q76xsec5,q76nfiles5,q76binnames5,indir,"SumHistsQCD80_HT700to1000.root",true,outdir);
        } else if (imode==16) {// QCD80 HT1000to1500
            MergeHists(goalintlum,q76nbin6,q76xsec6,q76nfiles6,q76binnames6,indir,"SumHistsQCD80_HT1000to1500.root",true,outdir);
        } else if (imode==17) {// QCD80 HT1500to2000
            MergeHists(goalintlum,q76nbin7,q76xsec7,q76nfiles7,q76binnames7,indir,"SumHistsQCD80_HT1500to2000.root",true,outdir);
        } else if (imode==18) {// QCD80 HT2000toInf
            MergeHists(goalintlum,q76nbin8,q76xsec8,q76nfiles8,q76binnames8,indir,"SumHistsQCD80_HT2000toInf.root",true,outdir);
        } else if (imode==26) {// QCD94 HT1000to1500
            MergeHists(goalintlum17,q94nbin6,q94xsec6,q94nfiles6,q94binnames6,indir,"SumHistsQCD94_HT1000to1500.root",true,outdir);
        } else if (imode==27) {// QCD94 HT1500to2000
            MergeHists(goalintlum17,q94nbin7,q94xsec7,q94nfiles7,q94binnames7,indir,"SumHistsQCD94_HT1500to2000.root",true,outdir);
        } else if (imode==28) {// QCD94 HT2000toInf
            MergeHists(goalintlum17,q94nbin8,q94xsec8,q94nfiles8,q94binnames8,indir,"SumHistsQCD94_HT2000toInf.root",true,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }

}

