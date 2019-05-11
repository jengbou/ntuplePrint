#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TFile.h>
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"

#include "vector"
#include <map>
using std::vector;

#include <TStyle.h>
#include <TCanvas.h>
#include "Style.hh"

double effDeepCSV(int bUnfType, double pt, int flav) {
    double f = 0;
    if (bUnfType==1) {
        //DeepCSVL
        if (flav==0) {//b and g->bb
            if (pt>=20 && pt<160){
                f=0.4344+0.02069*pt-0.0004429*pow(pt,2)+5.137*pow(10,-6)*pow(pt,3)-3.406*pow(10,-8)*pow(pt,4)+1.285*pow(10,-10)*pow(pt,5)-2.559*pow(10,-13)*pow(pt,6)+2.084*pow(10,-16)*pow(pt,7);
            }
            else if (pt>=160 && pt<300){
                f=0.714+0.002617*pt-1.656*pow(10,-5)*pow(pt,2)+4.767*pow(10,-8)*pow(pt,3)-6.431*pow(10,-11)*pow(pt,4)+3.287*pow(10,-14)*pow(pt,5);
            }
            else if (pt>=300 && pt<1000){
                f=0.872-6.885*pow(10,-5)*pt+4.34*pow(10,-8)*pow(pt,2);
            }
        }
        else if (flav==1) {//udsg
            if (pt>=20 && pt<150){
                f=0.245-0.0054*pt+6.92*pow(10,-5)*pow(pt,2)-3.89*pow(10,-7)*pow(pt,3)+1.021*pow(10,-9)*pow(pt,4)-1.007*pow(10,-12)*pow(pt,5);
            }
            else if (pt>=150 && pt<1000){
                f=0.0558+0.000428*pt-1.0*pow(10,-7)*pow(pt,2);
            }
        }
        else{//c
        }
    }
    else if (bUnfType==2) {//DeepCSVM
    }
    else if (bUnfType==3) {//DeepCSVT
    }
    else {//Undefined
    }

    return f;
}


void UnfoldWgt(){
    std::vector<float> jetpt(4);
    jetpt[0]= 225.;
    jetpt[1]= 150.;
    jetpt[2]= 100.;
    jetpt[3]= 100.;

    double effs_[2][4];
    for (int i=0;i<4;i++){
//         effs_[0][i]=effDeepCSV(1,jetpt[i],0);//b
//         effs_[1][i]=effDeepCSV(1,jetpt[i],1);//light
        effs_[0][i]=0.857;
        effs_[1][i]=0.0899;
    }

//     effs_[0][0]=0.857;
//     effs_[1][0]=0.0899;

    for (int i=0;i<2;i++){
        for (int j=0;j<4;j++){
            std::cout << "Effs(" << i << "," << j <<")=" << effs_[i][j] << std::endl;
        }
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                                //std::cout<<"ijkl=" << i << j << k << l <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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
                        //std::cout<<"ijk=" << i << j << k <<std::endl;
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

    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++){
            std::cout << "MA(" << i << "," << j <<")=" << MA(i,j) << std::endl;
        }
    }

    //MA.SetTol(1.e-23);
    //TDecompSVD svd(MA);
    TDecompLU lu(MA);
    TMatrixD MB(5,5);
    //if (!svd.Decompose()){
    if (!lu.Decompose()){
        for (int i=0;i<5;i++){
            for (int j=0;j<5;j++){
                std::cout << "MA(" << i << "," << j <<")=" << MA(i,j) << std::endl;
            }
        }
    }
    else {
        //MB=svd.Invert();
        MB=lu.Invert();
    }


    //TMatrixD MB=MA.Invert();
    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++){
            std::cout << "MB(" << i << "," << j <<")=" << MB(i,j) << std::endl;
        }
    }
}
