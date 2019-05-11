#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TFile.h>
#include <TH2F.h>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"
#include "TObject.h"
#include "TFitResult.h"
#include "TF1.h"

#include "vector"
#include <map>
using std::vector;

#include <TStyle.h>
#include <TCanvas.h>
#include "Style.hh"
#include <TCutG.h>

#include "CMS_lumi.C"

// Taken from TriggerTutorial/SimpleHLTAnalyzer
TGraphAsymmErrors* makeEffGraph(TString gtitle, TH1D* pass, TH1D* total, bool debug=false) {
    // make sure <pass> and <total> have the same binning!
    int npoints = total->GetNbinsX();

    float x[npoints], y[npoints], errx[npoints], erryl[npoints], erryh[npoints];

    float npass = 0.0;
    float ntotal = 0.0;

    for(int ibin = 1; ibin < npoints+1; ibin++) {
        x[ibin-1] = total->GetBinCenter(ibin);
        npass = pass->GetBinContent(ibin);
        ntotal = total->GetBinContent(ibin);
        y[ibin-1] = ntotal < 1.0 ? 0.0 : npass/ntotal;
        std::cout << "eff [" << ibin << "] = " << y[ibin-1] << std::endl;
        errx[ibin-1] = (total->GetXaxis()->GetBinUpEdge(ibin)-total->GetXaxis()->GetBinLowEdge(ibin))/2.0;
        if(y[ibin-1]==0.0) {
            erryl[ibin-1] = 0.0; erryh[ibin-1] = 0.0;
        } else {
            if(debug) printf("npass = %3.1f, ntotal = %3.1f, eff = %4.2f", npass, ntotal, y[ibin-1]);
            erryl[ibin-1] = y[ibin-1] - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, false);
            erryh[ibin-1] = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, true) - y[ibin-1];
        }
    }

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(npoints, x, y, errx, errx, erryl, erryh);
    gr->SetTitle(gtitle);

    return gr;

}


void GetBtagEff_CSVv2All_L() {
    //gStyle->SetOptStat(111111);
    std::map<TString,TFile*> files_;
    std::map<TString,TH1D*> hists_;
    std::map<TString,TH2F*> hists2_;
    std::map<TString,TGraphAsymmErrors*> geff_;
    std::map<TString,int> colors_;

    // original binning
//     Int_t  ptbinnum = 31;
//     Float_t ptbins[32] = {20, 30, 40, 50, 60, 70, 80, 90, 100,
//                           110, 120, 130, 140, 150, 160, 180, 200,
//                           220, 240, 260, 280, 300, 350, 400, 450,
//                           500, 550, 600, 700, 800, 900, 1000};

//     Int_t  ptbinnum = 33;
//     Float_t ptbins[34] = {20, 30, 40, 50, 60, 70, 80, 90, 100,
//                           110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
//                           220, 240, 260, 280, 300, 350, 400, 450,
//                           500, 550, 600, 700, 800, 900, 1000};

    Int_t  etabinnum = 4;
    Float_t etabins[5] = {-2.0,-1.4442,0.0,1.4442,2.0};
    // end original binning

    Int_t  ptbinnum = 9;
    Double_t ptbins[10] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
//     Int_t  etabinnum = 6;
//     Float_t etabins[7] = {-3.0,-2.0,-1.4442,0.0,1.4442,2.0,3.0};

    colors_["CSVv2L"] = 3;
    colors_["CSVv2M"] = 4;
    colors_["CSVv2T"] = 2;

    files_["CSVv2L"] = new TFile("histos/analysis_20180126_v0_p20180219_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//CSVv2L
    files_["CSVv2M"] = new TFile("histos/analysis_20180126_v0_p20180225_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//CSVv2M
    files_["CSVv2T"] = new TFile("histos/analysis_20180126_v0_p20180302_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//CSVv2T

    // ploting
    SetStyle();
    //TLegend *lgd = new TLegend(0.55,0.65,0.74,0.89);
    TLegend *lgd = new TLegend(0.22,0.2,0.52,0.45);
    //lgd->SetHeader("CSVv2L");
    lgd->SetHeader("light jets; UMD incl. c jets");
    lgd->SetBorderSize(0); lgd->SetTextSize(0.038); lgd->SetFillColor(0); // lgd->SetTextFont(62);

    TCanvas *cvs = MakeCanvas("effcanvas","",800,600);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";p_{T} [GeV];Relative Efficiency");

    const double ptmin_=20.0;
    const double ptmax_=1050.0;

    for(auto const &ent1 : files_) {
        TString fName = ent1.first;
        TString gName = "QCD "+ent1.first;
        //TString gName = "#splitline{QCD "+ent1.first+"}{(UMD kine. cuts)}";
        hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_udsg")->Clone("num_"+fName));
        hists2_["num"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c")->Clone("num_"+fName+"_c")));
        hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_udsg")->Clone("den_"+fName));
        hists2_["den"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c")->Clone("den_"+fName+"_c")));
        //hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c")->Clone("num_"+fName));
        //hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_c")->Clone("den_"+fName));
//         hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_udsg_all")->Clone("num_"+fName));
//         hists2_["num"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c_all")->Clone("num_"+fName+"_c")));
//         hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_udsg_all")->Clone("den_"+fName));
//         hists2_["den"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c_all")->Clone("den_"+fName+"_c")));

        hists_["num"+fName] = new TH1D("num"+fName, ";p_{T} [GeV];#epsilon", ptbinnum, ptbins);
        hists_["den"+fName] = new TH1D("den"+fName, ";p_{T} [GeV];#epsilon", ptbinnum, ptbins);

        hists2_["num"+fName]->ProjectionX("num"+fName);
        hists2_["den"+fName]->ProjectionX("den"+fName);

        geff_[fName] = makeEffGraph("eff_"+fName,hists_["num"+fName],hists_["den"+fName]);
        InitGraph(geff_[fName],"eff_"+fName,"H_{T} [GeV]","Relative Efficiency",ptmin_,ptmax_,0.0,1.1,colors_[fName]);
        geff_[fName]->SetMarkerStyle(kFullCircle);
        geff_[fName]->SetMarkerSize(0.8);
        mg->Add(geff_[fName]);
        lgd->AddEntry(geff_[fName],gName, "l");
    }
    mg->Draw("ap");
    //lgd->Draw();
    mg->GetXaxis()->SetLimits(ptmin_,ptmax_);
//     mg->SetMaximum(0.25);
//     mg->SetMinimum(0.10);
    //CSVv2L
    mg->SetMaximum(6e-1);
    mg->SetMinimum(1.e-4);

    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridy();
    gPad->Update();
    gPad->RedrawAxis();
    lgd->Draw();


    int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
    int iPos  = 11;
    CMS_lumi( cvs, iPeriod, iPos );
    gPad->Update();


    // coarse binning
    cvs->SaveAs("eff_btag_l_CSVv2All_vs_pt_BtagEff_r1_p20180312_500toInf_fine.pdf");
    cvs->SaveAs("eff_btag_l_CSVv2All_vs_pt_BtagEff_r1_p20180312_500toInf_fine.png");

    return;


}
