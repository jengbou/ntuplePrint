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


void GetBtagEff_CSVv2M_L_CoarseBin() {
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

    colors_["QCD"] = 1;

    files_["QCD"] = new TFile("histos/analysis_20180126_v0_p20180225_UMD_BtagEff_r1/SumHistsQCD80_HT500toInf.root");//CSVv2M

    // ploting
    SetStyle();
    //TLegend *lgd = new TLegend(0.55,0.2,0.8,0.44);
    TLegend *lgd = new TLegend(0.55,0.6,0.74,0.89);
    //lgd->SetHeader("CSVv2M");
    lgd->SetHeader("light jets; UMD incl. c jets");
    lgd->SetBorderSize(0); lgd->SetTextSize(0.03); lgd->SetFillColor(0);// lgd->SetTextFont(62);

    TCanvas *cvs = MakeCanvas("effcanvas","",800,600);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";p_{T} [GeV];Relative Efficiency");

    const double ptmin_=50.0;
    const double ptmax_=1050.0;

    for(auto const &ent1 : files_) {
        TString fName = ent1.first;
        TString gName = "#splitline{"+ent1.first+" CSVv2M}{(UMD kine. cuts)}";
        hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_udsg")->Clone("num_"+fName));
        hists2_["num"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c")->Clone("num_"+fName+"_c")));
        hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_udsg")->Clone("den_"+fName));
        hists2_["den"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c")->Clone("den_"+fName+"_c")));
//         hists2_["num"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Num_udsg_all")->Clone("num_"+fName));
//         hists2_["num"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c_all")->Clone("num_"+fName+"_c")));
//         hists2_["den"+fName] = static_cast<TH2F*>(ent1.second->Get("btagEff_Den_udsg_all")->Clone("den_"+fName));
//         hists2_["den"+fName]->Add(static_cast<TH2F*>(ent1.second->Get("btagEff_Num_c_all")->Clone("den_"+fName+"_c")));

        hists_["num"+fName] = new TH1D("num"+fName, ";p_{T} [GeV];#epsilon", ptbinnum, ptbins);
        hists_["den"+fName] = new TH1D("den"+fName, ";p_{T} [GeV];#epsilon", ptbinnum, ptbins);

        hists2_["num"+fName]->ProjectionX("num"+fName);
        hists2_["den"+fName]->ProjectionX("den"+fName);

        hists_["new_num"+fName] = (TH1D*) hists_["num"+fName]->Rebin(ptbinnum, TString("new_num"+fName).Data(), ptbins);
        hists_["new_den"+fName] = (TH1D*) hists_["den"+fName]->Rebin(ptbinnum, TString("new_num"+fName).Data(), ptbins);

        geff_[fName] = makeEffGraph("eff_"+fName,hists_["new_num"+fName],hists_["new_den"+fName]);
        InitGraph(geff_[fName],"eff_"+fName,"H_{T} [GeV]","Relative Efficiency",ptmin_,ptmax_,0.0,1.1,colors_[fName]);
        geff_[fName]->SetMarkerStyle(kFullCircle);
        geff_[fName]->SetMarkerSize(0.8);
        mg->Add(geff_[fName]);
        lgd->AddEntry(geff_[fName],gName, "l");
    }
    mg->Draw("ap");
    lgd->Draw();
    mg->GetXaxis()->SetLimits(ptmin_,ptmax_);

    // CSVv2M
    mg->SetMaximum(0.20);
    mg->SetMinimum(0.001);


    gPad->Update();

    gPad->RedrawAxis();

    //TF1 *func1 = new TF1("func1","[0]+[1]*x+[2]*pow(x,2)",ptmin_,160.0);
    //TF1 *func2 = new TF1("func2","[0]+[1]*x+[2]*pow(x,2)",160.0,450.0);
    //TF1 *func3 = new TF1("func3","[0]+[1]*x+[2]*pow(x,2)",550.0,1000.0);
    //TF1 *func2 = new TF1("func2","pol7",170.0,450.0);
    //TF1 *func3 = new TF1("func3","pol3",270.0,500.0);//250 475
    //TF1 *functot = new TF1("functot","pol2(0)+pol2(3)+pol2(6)",100.0,1000.0);
    //TF1 *functot = new TF1("functot","pol8",100.0,1000.0);

    TF1 *func1 = new TF1("func1","pol4",100.0,165.0);
    TF1 *func2 = new TF1("func2","pol5",165.0,270.0);
    TF1 *func3 = new TF1("func3","pol3",270.0,475.0);
    TF1 *func4 = new TF1("func4","pol4",475.0,1000.0);

    func1->SetLineColor(colors_["QCD"]);
    func2->SetLineColor(colors_["QCD"]);
    func3->SetLineColor(colors_["QCD"]);
    func4->SetLineColor(colors_["QCD"]);

//     TFitResultPtr resfit1 = mg->Fit("func1","SEMR","");
//     TFitResultPtr resfit2 = mg->Fit("func2","SEMR+","");
//     TFitResultPtr resfit3 = mg->Fit("func3","SEMR+","");
//     TFitResultPtr resfit4 = mg->Fit("func4","SEMR+","");

//     cout << "[FIT 1]: chi2/ndof = " << resfit1->Chi2()/(resfit1->Ndf()+1.) << endl;
//     cout << "[FIT 2]: chi2/ndof = " << resfit2->Chi2()/(resfit2->Ndf()+1.) << endl;
//     cout << "[FIT 3]: chi2/ndof = " << resfit3->Chi2()/(resfit3->Ndf()+1.) << endl;
//     cout << "[FIT 4]: chi2/ndof = " << resfit4->Chi2()/(resfit4->Ndf()+1.) << endl;


    //Plot BTagWG eff for demonstration
    // CSVv2M eff QCD MC
    TF1 *fCSVv2M1 = new TF1("fCSVv2M_30to65","0.0301219-0.00098619*x+1.61105*pow(10,-5)*pow(x,2)-1.13902*pow(10,-7)*pow(x,3)+2.993*pow(10,-10)*pow(x,4)",30,65);
    fCSVv2M1->SetLineColor(4);
//     fCSVv2M1->Draw("same");

    TF1 *fCSVv2M2 = new TF1("fCSVv2M_65to1000","0.00677158+2.04218*pow(10,-5)*x+1.50644*pow(10,-9)*pow(x,2)",65,1000);
    fCSVv2M2->SetLineColor(4);
    fCSVv2M2->Draw("same");

    lgd->AddEntry(fCSVv2M2,"QCD CSVv2M udsg jet (BTagWG)", "l");

    //gPad->SetLogx();
    // End CSVv2M eff QCD MC light jets


    // CSVv2M eff QCD MC c jets
    TF1 *fCSVv2Mc1 = new TF1("fCSVv2Mc_30to1000","0.157019-0.000175927*x+7.95035*pow(10,-8)*pow(x,2)",30,1000);
    fCSVv2Mc1->SetLineColor(3);
    fCSVv2Mc1->Draw("same");

    lgd->AddEntry(fCSVv2Mc1,"QCD CSVv2M c jet (BTagWG)", "l");

    //gPad->SetLogx();
    // End CSVv2M eff QCD MC c jets



    // DeepCSVM eff ttbar MC X SF
//     TF1 *fDeepCSVM1 = new TF1("fDeepCSVM_20to225",
//                               "0.019-0.00031*x+3.39*pow(10,-6)*pow(x,2)-1.47*pow(10,-8)*pow(x,3)+2.92*pow(10,-11)*pow(x,4)-2.12*pow(10,-14)*pow(x,5)",20,225);
    // DeepCSVM eff ttbar MC
    TF1 *fDeepCSVM1 = new TF1("fDeepCSVM_20to240",
                              "0.0186-0.000321*x+3.614*pow(10,-6)*pow(x,2)-1.683*pow(10,-8)*pow(x,3)+3.625*pow(10,-11)*pow(x,4)-2.875*pow(10,-14)*pow(x,5)",20,240);
    fDeepCSVM1->SetLineColor(2);
    fDeepCSVM1->Draw("same");

    //TF1 *fDeepCSVM2 = new TF1("fDeepCSVM_225to1000","0.00328+5.7*pow(10,-5)*x+4.7*pow(10,-9)*pow(x,2)",225,1000);//MC X SF
    TF1 *fDeepCSVM2 = new TF1("fDeepCSVM_240to1000","-0.0003747+6.858*pow(10,-5)*x-2.7*pow(10,-8)*pow(x,2)",240,1000);//MC
    fDeepCSVM2->SetLineColor(2);
    fDeepCSVM2->Draw("same");

    lgd->AddEntry(fDeepCSVM2,"t#bar{t} DeepCSVM udsg jet (BTagWG)", "l");
    // End DeepCSVM eff ttbar MC

    gPad->RedrawAxis();
    lgd->Draw();


//     TF1 *functot = new TF1("functot","pol6(0)+pol5(6)+pol3(11)+pol4(14)",100.0,1000.0);
//     functot->SetLineColor(2);

//     functot->SetParameters(func1->GetParameter(0),
//                            func1->GetParameter(1),
//                            func1->GetParameter(2),
//                            func1->GetParameter(3),
//                            func1->GetParameter(4),
//                            func1->GetParameter(5),
//                            func2->GetParameter(0),
//                            func2->GetParameter(1),
//                            func2->GetParameter(2),
//                            func2->GetParameter(3),
//                            func2->GetParameter(4)
//                            );
//     functot->SetParameter(12, func3->GetParameter(0));
//     functot->SetParameter(13, func3->GetParameter(1));
//     functot->SetParameter(14, func3->GetParameter(2));
//     functot->SetParameter(15, func4->GetParameter(0));
//     functot->SetParameter(16, func4->GetParameter(1));
//     functot->SetParameter(17, func4->GetParameter(2));
//     functot->SetParameter(18, func4->GetParameter(3));

//     TFitResultPtr resfittot = mg->Fit("functot","SEMR+","");

//     functot->Draw("same");

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

    // coarse bin
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_500toInf.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_500toInf.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_500toInf_all.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_500toInf_all.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_700toInf.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_700toInf.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_700toInf_all.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_p20180220_700toInf_all.png");

    // fine bin
    cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_500toInf.pdf");
    cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_500toInf.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_500toInf_all.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_500toInf_all.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_700toInf.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_700toInf.png");

//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_700toInf_all.pdf");
//     cvs->SaveAs("eff_btag_l_CSVv2M_vs_pt_analysis_20180126_v0_p20180225_UMD_BtagEff_r1_p20180306_700toInf_all.png");

    return;


}
