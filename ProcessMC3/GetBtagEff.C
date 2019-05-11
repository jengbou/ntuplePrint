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


void GetBtagEff() {
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

    Int_t  ptbinnum = 33;
    Float_t ptbins[34] = {20, 30, 40, 50, 60, 70, 80, 90, 100,
                          110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                          220, 240, 260, 280, 300, 350, 400, 450,
                          500, 550, 600, 700, 800, 900, 1000};

    Int_t  etabinnum = 4;
    Float_t etabins[5] = {-2.0,-1.4442,0.0,1.4442,2.0};
    // end original binning

//     Int_t  ptbinnum = 9;
//     Float_t ptbins[10] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};
//     Int_t  etabinnum = 6;
//     Float_t etabins[7] = {-3.0,-2.0,-1.4442,0.0,1.4442,2.0,3.0};

    colors_["QCD"] = 1;

    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180201_UMD_BtagEff_r2/SumHistsQCD80_HT500toInf.root");
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180201_UMD_BtagEff_r2/SumHistsQCD80_HT700toInf.root");
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180206_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//finer binning
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180206_UMD_BtagEff_r1/SumHistsQCD80_HT700toInf.root");//finer binning
    files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180219_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//new finer binning
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180219_UMD_BtagEff_r1/SumHistsQCD80_HT700toInf.root");//new finer binning
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180221_UMD_BtagEff_r1/SumHistsQCD80_HT700toInf.root");//new finer binning
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180221_UMD_BtagEff_r1/SumHistsQCD80_HT700toInf.root");//new finer binning
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180225_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//CSVv2M
    //files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180225_UMD_BtagEff_r2/SumHistsQCD80_HT500toInf.root");//CSVv2 user

    // ploting
    SetStyle();
    TLegend *lgd = new TLegend(0.55,0.65,0.74,0.89);
    //lgd->SetHeader("CSVv2L");
    lgd->SetHeader("b jets");
    lgd->SetBorderSize(0); lgd->SetTextSize(0.038); lgd->SetFillColor(0); // lgd->SetTextFont(62);

    TCanvas *cvs = MakeCanvas("effcanvas","",800,600);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";p_{T} [GeV];Relative Efficiency");

    const double ptmin_=50.0;
    const double ptmax_=1050.0;

    for(auto const &ent1 : files_) {
        TString fName = ent1.first;
        TString gName = "#splitline{"+ent1.first+" CSVv2L}{(UMD kine. cuts)}";
        hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_b")->Clone("num_"+fName));
        hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_b")->Clone("den_"+fName));
//         hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_b_all")->Clone("num_"+fName));
//         hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_b_all")->Clone("den_"+fName));

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
    mg->SetMaximum(1.00);
    mg->SetMinimum(0.58);
//     mg->SetMaximum(0.85);
//     mg->SetMinimum(0.65);
    gPad->Update();

    //gPad->RedrawAxis();


    TF1 *func1 = new TF1("func1","pol3",100.0,175.0);
    TF1 *func2 = new TF1("func2","pol3",175.0,290.0);
    TF1 *func3 = new TF1("func3","pol4",290.0,575.0);
    TF1 *func4 = new TF1("func4","pol3",575.0,1000.0);
    //TF1 *func5 = new TF1("func5","pol3",525.0,1000.0);

    func1->SetLineColor(colors_["QCD"]);
    func2->SetLineColor(colors_["QCD"]);
    func3->SetLineColor(colors_["QCD"]);
    func4->SetLineColor(colors_["QCD"]);

    TFitResultPtr resfit1 = mg->Fit("func1","SEMR","");
    TFitResultPtr resfit2 = mg->Fit("func2","SEMR+","");
    TFitResultPtr resfit3 = mg->Fit("func3","SEMR+","");
    TFitResultPtr resfit4 = mg->Fit("func4","SEMR+","");
    //TFitResultPtr resfit5 = mg->Fit("func5","SEMR+","");

    cout << "[FIT 1]: chi2/ndof = " << resfit1->Chi2()/(resfit1->Ndf()+1.) << endl;
    cout << "[FIT 2]: chi2/ndof = " << resfit2->Chi2()/(resfit2->Ndf()+1.) << endl;
    cout << "[FIT 3]: chi2/ndof = " << resfit3->Chi2()/(resfit3->Ndf()+1.) << endl;
    cout << "[FIT 4]: chi2/ndof = " << resfit4->Chi2()/(resfit4->Ndf()+1.) << endl;
    //cout << "[FIT 5]: chi2/ndof = " << resfit5->Chi2()/(resfit5->Ndf()+1.) << endl;

    //Plot BTagWG eff for demonstration
    // CSVv2L eff QCD MC
    TF1 *fCSVv2L1 = new TF1("fCSVv2L_30to165","0.625404+0.00629137*x-7.44931*pow(10,-5)*pow(x,2)+3.79928*pow(10,-7)*pow(x,3)-7.23654*pow(10,-10)*pow(x,4)",30,165);
    fCSVv2L1->SetLineColor(4);
    fCSVv2L1->Draw("same");

    TF1 *fCSVv2L2 = new TF1("fCSVv2L_165to1000","0.845173-0.000226895*x-2.12743*pow(10,-8)*pow(x,2)",165,1000);
    fCSVv2L2->SetLineColor(4);
    fCSVv2L2->Draw("same");

    lgd->AddEntry(fCSVv2L2,"QCD CSVv2L (BTagWG)", "l");

    //gPad->SetLogx();
    // End CSVv2L eff QCD MC


    // DeepCSVL eff ttbar MC
    TF1 *fDeepCSVL1 = new TF1("fDeepCSVL_20to100",
                              "0.491+0.0191*x-0.0004172*pow(x,2)+4.893*pow(10,-6)*pow(x,3)-3.266*pow(10,-8)*pow(x,4)+1.238*pow(10,-10)*pow(x,5)-2.474*pow(10,-13)*pow(x,6)+2.021*pow(10,-16)*pow(x,7)",20,100);
    fDeepCSVL1->SetLineColor(2);
    fDeepCSVL1->Draw("same");
    TF1 *fDeepCSVL2 = new TF1("fDeepCSVL_100to300","0.912-0.001846*x+2.479*pow(10,-5)*pow(x,2)-1.417*pow(10,-7)*pow(x,3)+3.617*pow(10,-10)*pow(x,4)-3.433*pow(10,-13)*pow(x,5)",100,300);
    fDeepCSVL2->SetLineColor(2);
    fDeepCSVL2->Draw("same");

    TF1 *fDeepCSVL3 = new TF1("fDeepCSVL_300to1000","0.892-0.00014*x+1.01*pow(10,-7)*pow(x,2)",300,1000);
    fDeepCSVL3->SetLineColor(2);
    fDeepCSVL3->Draw("same");

    lgd->AddEntry(fDeepCSVL2,"t#bar{t} DeepCSVL (BTagWG)", "l");
    // End DeepCSVL eff ttbar MC


    gPad->RedrawAxis();
    lgd->Draw();

//     TFitResultPtr resfit = mg->Fit("pol2","SEM0","");

//     TString sfunc;
//     sfunc.Form("#epsilon (p_{T}) = %7.6g + %7.6g #times p_{T} + %7.6g #times p_{T}^{2}",resfit->Parameter(0),resfit->Parameter(1),resfit->Parameter(2));
//     sfunc.ReplaceAll("e-08","x10^{-8}");

//     TF1 *func_ = new TF1("func_","[0]+[1]*x+[2]*pow(x,2)",ptmin_,ptmax_);
//     func_->SetParameters(resfit->Parameter(0),resfit->Parameter(1),resfit->Parameter(2));
//     func_->SetLineColor(2);
//     func_->Draw("same");

//     TLatex *tfunc = new TLatex();
//     tfunc->SetNDC(true);
//     tfunc->SetTextAlign(12);
//     tfunc->SetTextSize(0.03);
//     tfunc->SetTextColor(2);
//     tfunc->DrawLatex(0.4,0.42,sfunc);

//     cvs->Update();

//     TLegend *lgdf = new TLegend(0.45,0.4,0.9,0.5);
//     lgdf->SetBorderSize(0); lgdf->SetTextSize(0.025); lgdf->SetTextFont(82); lgdf->SetFillColor(0);
//     lgdf->AddEntry(func_,sfunc.Data(),"l");
//     lgdf->Draw();

    int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
    int iPos  = 11;
    CMS_lumi( cvs, iPeriod, iPos );
    gPad->Update();


//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180206_700toInf.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180206_700toInf.png");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180213_700toInf_fine.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180213_700toInf_fine.png");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180213_700toInf_fine_all.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_20180213_700toInf_fine_all.png");


    // coarse binning
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf_all.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf_all.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_700toInf.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_700toInf.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_700toInf_all.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_700toInf_all.png");

    // fine binning
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf_fine_UMDonly.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_p20180220_500toInf_fine_UMDonly.png");

    cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_500toInf_fine.pdf");
    cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_500toInf_fine.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_500toInf_fine_all.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_500toInf_fine_all.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_700toInf_fine.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_700toInf_fine.png");

//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_700toInf_fine_all.pdf");
//     cvs->SaveAs("eff_btag_b_CSVv2L_vs_pt_analysis_20180126_v0_p20180219_UMD_BtagEff_r1_p20180301_700toInf_fine_all.png");

    return;


}
