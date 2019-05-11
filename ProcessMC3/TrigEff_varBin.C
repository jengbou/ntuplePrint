#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TFile.h>

#include "vector"
using std::vector;

#include <TStyle.h>
#include <TCanvas.h>


void FillVariantBinHist (TH1F& hin, TH1F& hout){
    for (int i=0;i<hin.GetNbinsX();i++){
        //std::cout << hin.GetBinContent(i+1) << std::endl;
        for (int j=0;j<hout.GetNbinsX();j++){
            if (j<hout.GetNbinsX()-1){
                if (hin.GetXaxis()->GetBinLowEdge(i+1)>=hout.GetXaxis()->GetBinLowEdge(j+1) &&
                    hin.GetXaxis()->GetBinLowEdge(i+2)<hout.GetXaxis()->GetBinLowEdge(j+2)){
                    double tmpVal = hout.GetBinContent(j+1)+hin.GetBinContent(i+1);
                    hout.SetBinContent(j+1,tmpVal);
                    //hout.Fill(j+1,hin.GetBinContent(i+1);
                }
            }
            else{
                if (hin.GetXaxis()->GetBinLowEdge(i+1)>=hout.GetXaxis()->GetBinLowEdge(j+1)){
                    double tmpVal = hout.GetBinContent(j+1)+hin.GetBinContent(i+1);
                    hout.SetBinContent(j+1,tmpVal);
                    //hout.Fill(j+1,hin.GetBinContent(i+1);
                }
            }
        }
    }
}


void TrigEff_varBin() {

    gStyle->SetOptStat(111111);

    TString canvName = "Fig_";
    canvName += "HLTeff";

    //  if( writeExtraText ) canvName += "-prelim";
    //if( iPos%10==0 ) canvName += "-out";                                        
    //else if( iPos%10==1 ) canvName += "-left";                                  
    //else if( iPos%10==2 )  canvName += "-center";                               
    //else if( iPos%10==3 )  canvName += "-right";                                
    int W = 800;
    int H = 600;
    TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
    // references for T, B, L, R                                                  
    float T = 0.08*H;
    float B = 0.12*H;
    float L = 0.12*W;
    float R = 0.04*W;

    canv->SetFillColor(0);
    canv->SetBorderMode(0);
    canv->SetFrameFillStyle(0);
    canv->SetFrameBorderMode(0);
    canv->SetLeftMargin( L/W );
    canv->SetRightMargin( R/W );
    canv->SetTopMargin( T/H );
    canv->SetBottomMargin( B/H );
    canv->SetTickx(0);
    canv->SetTicky(0);


    int n_ = 2;

    float x1_l = 0.9;
    //  float x1_l = 0.75;                                                                  
    float y1_l = 0.80;

    float dx_l = 0.60;
    float dy_l = 0.1;
    float x0_l = x1_l-dx_l;
    float y0_l = y1_l-dy_l;

    TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l);
    lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);

    TFile *f = new TFile("histos/analysis_20170223_v0_UMD/SumHistsQCD74.root");
    float hltBins[26] = {0.,500.,600.,650.,700.,720.,740.,760.,780.,800.,820.,840.,860.,880.,900.,925.,950.,975.,1000.,1050.,1100.,1200.,1300.,1400.,1600.,2000.};
    TH1F* hTrig3num = new TH1F("hTrig3num","HLT_PFHT900",25,hltBins);
    TH1F* hTrig3den = new TH1F("hTrig3den","HLT_PFHT600",25,hltBins);
    hTrig3num->Sumw2();
    hTrig3den->Sumw2();

    TH1F* nhist = static_cast<TH1F*>(f->Get("hTrig3n")->Clone("num_HLT"));
    TH1F* dhist = static_cast<TH1F*>(f->Get("hTrig3d")->Clone("den_HLT"));

    FillVariantBinHist(*nhist,*hTrig3num);
    FillVariantBinHist(*dhist,*hTrig3den);
    TH1F* rhTrig  = static_cast<TH1F*>(hTrig3num->Clone("effHLT"));
    //rhTrig->Sumw2();

    bool divideOK = rhTrig->Divide(rhTrig, hTrig3den, 1., .1, "cl=0.683 cp");
    //bool divideOK = rhTrig->Divide(hTrig3den);
//     for (int i=0;i<rhTrig->GetNbinsX();i++){
//         std::cout << rhTrig->GetBinContent(i+1) << std::endl;
//     }
    if (!divideOK) return;
    rhTrig->SetMaximum(1.2);
    rhTrig->SetMinimum(0.0);
    rhTrig->GetXaxis()->SetTitle("H_{T}");
    rhTrig->GetYaxis()->SetTitle("#epsilon_{HLT}");

    rhTrig->SetMarkerStyle(20);
    rhTrig->SetMarkerColor(2);
    //rhTrig->SetMarkerSize(1);
    rhTrig->SetLineColor(2);
    rhTrig->SetLineWidth(2);
    rhTrig->SetStats(0);
    rhTrig->Draw("p0 e0");

    lgd->AddEntry(rhTrig, "CR: N_{tagged jets} < 2", "l");

    lgd->Draw();


    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();


    return;


}
