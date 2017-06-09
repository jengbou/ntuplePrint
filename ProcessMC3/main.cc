#include <iostream>
#include <string>
#include <map>
#include "QCDhists.h"
#include "QCDhistsNoMerge.h"
#include "MergeHists.h"
#include "MergeHistsNoNorm.h"
#include "MergeHistsN.h"

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
    } else if(imode==10) {
        std::cout<<"doing DATA"<<std::endl;
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
    } else {
        std::cout<<"invalid choice"<<std::endl;
    }


    float goalintlum=20.0; // fb-1                                                                                        
    //float goalintlum=0.07956; // fb-1                                                                                        

    std::string aaname;
    if (imode==9) aaname = "/data/users/jengbou/EmJetMC/";  // RelVal samples
    else aaname = "/data/users/jengbou/crab_output/ntuple_20170523_v0/";//QCD80 HT1000+

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


    // mode = 10
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
    int anfiles[anbin]={31};//74Full 716
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

    if (pmode==0) {
        std::cout << "All in one go." << std::endl;
        if(imode==0) {
            QCDhists(goalintlum,nbin,xsec,nfiles,binnames,aaname,"SumHistsQCD.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==1) {
            QCDhists(goalintlum,anbin,axsec,anfiles,abinnames,aaname,"SumHistsModelA.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==2) {
            QCDhists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname,"SumHistsModelB.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==3) {//QCD80 HT500to700
            QCDhists(goalintlum,1,q74xsec1,q74nfiles1,q74binnames1,aaname,"SumHistsQCD80_HT500to700.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==4) {// QCD80 HT700to1000
            QCDhists(goalintlum,1,q74xsec2,q74nfiles2,q74binnames2,aaname,"SumHistsQCD80_HT700to1000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==5) {// QCD80 HT1000to1500
            QCDhists(goalintlum,1,q74xsec3,q74nfiles3,q74binnames3,aaname,"SumHistsQCD80_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==6) {// QCD80 HT1500to2000
            QCDhists(goalintlum,1,q74xsec4,q74nfiles4,q74binnames4,aaname,"SumHistsQCD80_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==7) {// QCD80 HT2000toInf
            QCDhists(goalintlum,1,q74xsec5,q74nfiles5,q74binnames5,aaname,"SumHistsQCD80_HT200toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==8) {
            QCDhists(goalintlum,datanbin,dataxsec,datanfiles,databinnames,aaname,"SumHistsDATA.root",0,0,hasPre,false,blind,b16003,outdir,true);
        } else if (imode==9) {
            QCDhists(goalintlum,vnbin,vxsec,vnfiles,vbinnames,aaname+"RelVal/","SumHistsModelARelValEXO16003.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==11) {//QCD80 HT100to200
            QCDhists(goalintlum,1,q76xsec1,q76nfiles1,q76binnames1,aaname,"SumHistsQCD80_HT100to200.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==12) {// QCD80 HT200to300
            QCDhists(goalintlum,1,q76xsec2,q76nfiles2,q76binnames2,aaname,"SumHistsQCD80_HT200to300.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==13) {// QCD80 HT300to500
            QCDhists(goalintlum,1,q76xsec3,q76nfiles3,q76binnames3,aaname,"SumHistsQCD80_HT300to500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==14) {//QCD80 HT500to700
            QCDhists(goalintlum,1,q76xsec1,q76nfiles1,q76binnames4,aaname,"SumHistsQCD80_HT500to700.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==15) {// QCD80 HT700to1000
            QCDhists(goalintlum,1,q76xsec2,q76nfiles2,q76binnames5,aaname,"SumHistsQCD80_HT700to1000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==16) {// QCD80 HT1000to1500
            QCDhists(goalintlum,1,q76xsec3,q76nfiles3,q76binnames6,aaname,"SumHistsQCD80_HT1000to1500.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==17) {// QCD80 HT1500to2000
            QCDhists(goalintlum,1,q76xsec4,q76nfiles4,q76binnames7,aaname,"SumHistsQCD80_HT1500to2000.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        } else if (imode==18) {// QCD80 HT2000toInf
            QCDhists(goalintlum,1,q76xsec5,q76nfiles5,q76binnames8,aaname,"SumHistsQCD80_HT200toInf.root",dooptk,doopta,hasPre,true,blind,b16003,outdir,true);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==1) {
        std::cout << "No merging." << std::endl;
        if (imode==1) {// modelA
            QCDhistsNoMerge(nrange,"modelA",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==2) {// modelB 
            QCDhistsNoMerge(nrange,"modelB",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==3) {//QCD74 HT500to700
            QCDhistsNoMerge(nrange,"QCD74_HT500to700",aaname,hasPre,blind,b16003,outdir,false);
        } else if (imode==4) {// QCD74 HT700to1000
            QCDhistsNoMerge(nrange,"QCD74_HT700to1000",aaname,hasPre,blind,b16003,outdir,false);
        } else if (imode==5) {// QCD74 HT1000to1500
            QCDhistsNoMerge(nrange,"QCD74_HT1000to1500",aaname,hasPre,blind,b16003,outdir,false);
        } else if (imode==6) {// QCD74 HT1500to2000
            QCDhistsNoMerge(nrange,"QCD74_HT1500to2000",aaname,hasPre,blind,b16003,outdir,false);
        } else if (imode==7) {// QCD74 HT2000toInf
            QCDhistsNoMerge(nrange,"QCD74_HT2000toInf",aaname,hasPre,blind,b16003,outdir,false);
//         } else if (imode==3) {//QCD74 HT500to700
//             QCDhistsNoMerge(nrange,"QCD74_HT500to700",aaname,hasPre,blind,b16003,outdir,true);
//         } else if (imode==4) {// QCD74 HT700to1000
//             QCDhistsNoMerge(nrange,"QCD74_HT700to1000",aaname,hasPre,blind,b16003,outdir,true);
//         } else if (imode==5) {// QCD74 HT1000to1500
//             QCDhistsNoMerge(nrange,"QCD74_HT1000to1500",aaname,hasPre,blind,b16003,outdir,true);
//         } else if (imode==6) {// QCD74 HT1500to2000
//             QCDhistsNoMerge(nrange,"QCD74_HT1500to2000",aaname,hasPre,blind,b16003,outdir,true);
//         } else if (imode==7) {// QCD74 HT2000toInf
//             QCDhistsNoMerge(nrange,"QCD74_HT2000toInf",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==11) {//QCD80 HT100to200
            QCDhistsNoMerge(nrange,"QCD80_HT100to200",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==12) {// QCD80 HT200to300
            QCDhistsNoMerge(nrange,"QCD80_HT200to300",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==13) {// QCD80 HT300to500
            QCDhistsNoMerge(nrange,"QCD80_HT300to500",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==14) {//QCD80 HT500to700
            QCDhistsNoMerge(nrange,"QCD80_HT500to700",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==15) {// QCD80 HT700to1000
            QCDhistsNoMerge(nrange,"QCD80_HT700to1000",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==16) {// QCD80 HT1000to1500
            QCDhistsNoMerge(nrange,"QCD80_HT1000to1500",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==17) {// QCD80 HT1500to2000
            QCDhistsNoMerge(nrange,"QCD80_HT1500to2000",aaname,hasPre,blind,b16003,outdir,true);
        } else if (imode==18) {// QCD80 HT2000toInf
            QCDhistsNoMerge(nrange,"QCD80_HT2000toInf",aaname,hasPre,blind,b16003,outdir,true);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==8) {
        if (imode==11) {//QCD80 HT100to200
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT100to200",aaname,outdir);
        }
        else if (imode==12) {//QCD80 HT200to300
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT200to300",aaname,outdir);
        }
        else if (imode==13) {//QCD80 HT300to500
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT300to500",aaname,outdir);
        }
        else if (imode==14) {//QCD80 HT500to700
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT500to700",aaname,outdir);
        }
        else if (imode==15) {//QCD80 HT700to1000
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT700to1000",aaname,outdir);
        }
        else if (imode==16) {//QCD80 HT1000to1500
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT1000to1500",aaname,outdir);
        }
        else if (imode==17) {//QCD80 HT1500to2000
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT1500to2000",aaname,outdir);
        }
        else if (imode==18) {//QCD80 HT2000toInf
            MergeHistsNoNorm(fidx,nrange,"QCD80_HT2000toInf",aaname,outdir);
        }
        else if (imode==3) {//QCD74 HT500to700
            MergeHistsNoNorm(fidx,nrange,"QCD74_HT500to700",aaname,outdir);
        }
        else if (imode==4) {//QCD74 HT700to1000
            MergeHistsNoNorm(fidx,nrange,"QCD74_HT700to1000",aaname,outdir);
        }
        else if (imode==5) {//QCD74 HT1000to1500
            MergeHistsNoNorm(fidx,nrange,"QCD74_HT1000to1500",aaname,outdir);
        }
        else if (imode==6) {//QCD74 HT1500to2000
            MergeHistsNoNorm(fidx,nrange,"QCD74_HT1500to2000",aaname,outdir);
        }
        else if (imode==7) {//QCD74 HT2000toInf
            MergeHistsNoNorm(fidx,nrange,"QCD74_HT2000toInf",aaname,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else if (pmode==9) {    //float q74xsec[q74nbin]={29370000,6524000,1064000,121500,25420}; // fb 
        if (imode==11) {//QCD80 HT100to200
            MergeHistsN(goalintlum,27990000000,fidx,"QCD80_HT100to200","SumHistsQCD80_HT100to200.root",true,outdir);
        }
        else if (imode==12) {//QCD80 HT200to300
            MergeHistsN(goalintlum,1712000000,fidx,"QCD80_HT200to300","SumHistsQCD80_HT200to300.root",true,outdir);
        }
        else if (imode==13) {//QCD80 HT300to500
            MergeHistsN(goalintlum,347700000,fidx,"QCD80_HT300to500","SumHistsQCD80_HT300to500.root",true,outdir);
        }
        else if (imode==14) {//QCD80 HT500to700
            //MergeHistsN(goalintlum,29370000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,32100000,fidx,"QCD80_HT500to700","SumHistsQCD80_HT500to700.root",true,outdir);
        }
        else if (imode==15) {//QCD80 HT700to1000
            //MergeHistsN(goalintlum,6524000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,6831000,fidx,"QCD80_HT700to1000","SumHistsQCD80_HT700to1000.root",true,outdir);
        }
        else if (imode==16) {//QCD80 HT1000to1500
            //MergeHistsN(goalintlum,1064000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,1207000,fidx,"QCD80_HT1000to1500","SumHistsQCD80_HT1000to1500.root",true,outdir);
        }
        else if (imode==17) {//QCD80 HT1500to2000
            //MergeHistsN(goalintlum,121500,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,119900,fidx,"QCD80_HT1500to2000","SumHistsQCD80_HT1500to2000.root",true,outdir);
        }
        else if (imode==18) {//QCD80 HT2000toInf
            //MergeHistsN(goalintlum,25420,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,25240,fidx,"QCD80_HT2000toInf","SumHistsQCD80_HT2000toInf.root",true,outdir);
        }
        else if (imode==3) {//QCD74 HT500to700
            //MergeHistsN(goalintlum,29370000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,32100000,fidx,"QCD74_HT500to700","SumHistsQCD74_HT500to700.root",true,outdir);
        }
        else if (imode==4) {//QCD74 HT700to1000
            //MergeHistsN(goalintlum,6524000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,6831000,fidx,"QCD74_HT700to1000","SumHistsQCD74_HT700to1000.root",true,outdir);
        }
        else if (imode==5) {//QCD74 HT1000to1500
            //MergeHistsN(goalintlum,1064000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,1207000,fidx,"QCD74_HT1000to1500","SumHistsQCD74_HT1000to1500.root",true,outdir);
        }
        else if (imode==6) {//QCD74 HT1500to2000
            //MergeHistsN(goalintlum,121500,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,119900,fidx,"QCD74_HT1500to2000","SumHistsQCD74_HT1500to2000.root",true,outdir);
        }
        else if (imode==7) {//QCD74 HT2000toInf
            //MergeHistsN(goalintlum,25420,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf_oldXsec.root",true,outdir);
            MergeHistsN(goalintlum,25240,fidx,"QCD74_HT2000toInf","SumHistsQCD74_HT2000toInf.root",true,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }
    else{//Merging only
        if (imode==1) {// modelA
            MergeHists(goalintlum,anbin,axsec,anfiles,abinnames,aaname,"SumHistsModelA.root",true,outdir);
        } else if (imode==2) {// modelB 
            MergeHists(goalintlum,bnbin,bxsec,bnfiles,bbinnames,aaname,"SumHistsModelB.root",true,outdir);
        } else if (imode==3) {//QCD74 HT500to700
            MergeHists(goalintlum,q74nbin1,q74xsec1,q74nfiles1,q74binnames1,aaname,"SumHistsQCD74_HT500to700.root",true,outdir);
        } else if (imode==4) {// QCD74 HT700to1000
            MergeHists(goalintlum,q74nbin2,q74xsec2,q74nfiles2,q74binnames2,aaname,"SumHistsQCD74_HT700to1000.root",true,outdir);
        } else if (imode==5) {// QCD74 HT1000to1500
            MergeHists(goalintlum,q74nbin3,q74xsec3,q74nfiles3,q74binnames3,aaname,"SumHistsQCD74_HT1000to1500.root",true,outdir);
        } else if (imode==6) {// QCD74 HT1500to2000
            MergeHists(goalintlum,q74nbin4,q74xsec4,q74nfiles4,q74binnames4,aaname,"SumHistsQCD74_HT1500to2000.root",true,outdir);
        } else if (imode==7) {// QCD74 HT2000toInf
            MergeHists(goalintlum,q74nbin5,q74xsec5,q74nfiles5,q74binnames5,aaname,"SumHistsQCD74_HT2000toInf.root",true,outdir);
        } else if (imode==11) {//QCD80 HT100to200
            MergeHists(goalintlum,q76nbin1,q76xsec1,q76nfiles1,q76binnames1,aaname,"SumHistsQCD80_HT100to200.root",true,outdir);
        } else if (imode==12) {// QCD80 HT200to300
            MergeHists(goalintlum,q76nbin2,q76xsec2,q76nfiles2,q76binnames2,aaname,"SumHistsQCD80_HT200to300.root",true,outdir);
        } else if (imode==13) {// QCD80 HT300to500
            MergeHists(goalintlum,q76nbin3,q76xsec3,q76nfiles3,q76binnames3,aaname,"SumHistsQCD80_HT300to500.root",true,outdir);
        } else if (imode==14) {//QCD80 HT500to700
            MergeHists(goalintlum,q76nbin4,q76xsec4,q76nfiles4,q76binnames4,aaname,"SumHistsQCD80_HT500to700.root",true,outdir);
        } else if (imode==15) {// QCD80 HT700to1000
            MergeHists(goalintlum,q76nbin5,q76xsec5,q76nfiles5,q76binnames5,aaname,"SumHistsQCD80_HT700to1000.root",true,outdir);
        } else if (imode==16) {// QCD80 HT1000to1500
            MergeHists(goalintlum,q76nbin6,q76xsec6,q76nfiles6,q76binnames6,aaname,"SumHistsQCD80_HT1000to1500.root",true,outdir);
        } else if (imode==17) {// QCD80 HT1500to2000
            MergeHists(goalintlum,q76nbin7,q76xsec7,q76nfiles7,q76binnames7,aaname,"SumHistsQCD80_HT1500to2000.root",true,outdir);
        } else if (imode==18) {// QCD80 HT2000toInf
            MergeHists(goalintlum,q76nbin8,q76xsec8,q76nfiles8,q76binnames8,aaname,"SumHistsQCD80_HT2000toInf.root",true,outdir);
        }
        else {
            std::cout<<"invalid choice"<<std::endl;
        }
    }

}

