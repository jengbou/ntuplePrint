#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TFile.h>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"

#include "vector"
#include <map>
using std::vector;

#include <TStyle.h>
#include <TCanvas.h>
#include "Style.hh"

#include "CMS_lumi.C"

// Taken from TriggerTutorial/SimpleHLTAnalyzer
TGraphAsymmErrors* makeEffGraph(TString gtitle, TH1F* pass, TH1F* total, bool debug=false) {
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
        errx[ibin-1] = 0.0;
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


void TrigEff_EXO() {
    //gStyle->SetOptStat(111111);
    std::map<TString,TFile*> files_;
    std::map<TString,TH1F*> hists_;
    std::map<TString,TGraphAsymmErrors*> geff_;
    std::map<TString,int> colors_;

    colors_["QCD"] = 1;
    colors_["ModelA"] = 2;
    colors_["ModelB"] = 4;

    files_["QCD"] = new TFile("histos/analysis_20170426_v2_p20170517_UMD_test7/SumHistsQCD80_HT1000toInf.root");
    files_["ModelA"] = new TFile("histos/analysis_20170426_v2_p20170517_UMD_test7/SumHistsModelA.root");
    files_["ModelB"] = new TFile("histos/analysis_20170426_v2_p20170517_UMD_test7/SumHistsModelB.root");

    // ploting
    SetStyle();
    TLegend *lgd = new TLegend(0.65,0.2,0.9,0.4);
    //lgd->SetHeader("HLT_PFHT900");//74
    lgd->SetHeader("HLT_PFHT800");//76&80
    lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);

    TCanvas *cvs = MakeCanvas("effcanvas","",800,600);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";H_{T} [GeV];Relative Efficiency");

    const double htmin_=300.0;
    const double htmax_=1800.0;

    for(auto const &ent1 : files_) {
        TString fName = ent1.first;
        hists_["num"+fName] = static_cast<TH1F*>(ent1.second->Get("hTrig3n")->Clone("num_"+fName));
        hists_["den"+fName] = static_cast<TH1F*>(ent1.second->Get("hTrig3d")->Clone("den_"+fName));
        hists_["num"+fName]->Rebin(2);
        hists_["den"+fName]->Rebin(2);

        geff_[fName] = makeEffGraph("eff_"+fName,hists_["num"+fName],hists_["den"+fName]);
        InitGraph(geff_[fName],"eff_"+fName,"H_{T} [GeV]","Relative Efficiency",htmin_,htmax_,0.0,1.1,colors_[fName]);
        geff_[fName]->SetMarkerStyle(kFullCircle);
        geff_[fName]->SetMarkerSize(0.8);
        mg->Add(geff_[fName]);
        lgd->AddEntry(geff_[fName],fName, "l");
    }
    mg->Draw("ap");
    lgd->Draw();
    mg->GetXaxis()->SetLimits(htmin_,htmax_);
    gPad->Update();

    gPad->RedrawAxis();

    int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
    int iPos  = 11;
    CMS_lumi( cvs, iPeriod, iPos );
    gPad->Update();

    cvs->SaveAs("eff_trig_vs_ht_20170426_v0_p20170517_test7.pdf");
    cvs->SaveAs("eff_trig_vs_ht_20170426_v0_p20170517_test7.png");

    return;


}
