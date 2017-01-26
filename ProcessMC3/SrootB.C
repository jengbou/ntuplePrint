  
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

#include "vector"


bool pairCompare(const std::pair<int,double>& firstElem, const std::pair<int,double>& secondElem) {
  return firstElem.second < secondElem.second;

}


int dolog=0;
void SrootB() 
{ 
 
  gStyle->SetOptStat(0);
 
  TString canvName = "Fig_";
  canvName += "hptdp_A_B";
  
  if( writeExtraText ) canvName += "-prelim";
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
  
  //canv = new TCanvas(canvName,canvName,50,50,W,H);
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
  
  if (dolog) canv->SetLogy();

  //TH1* h_pt = new TH1F("h_pt"," ",100,0,500);
  //h_pt->GetXaxis()->SetNdivisions(6,5,0);
  //h_pt->GetXaxis()->SetTitle("Dark Pion p_{T} (GeV)");  
  //h_pt->GetXaxis()->SetTitleSize(0.05);  
  //h_pt->GetYaxis()->SetNdivisions(6,5,0);
  //h_pt->GetYaxis()->SetTitleOffset(1);
  //h_pt->GetYaxis()->SetTitle("Events / 5 GeV");  
  //h_pt->GetYaxis()->SetTitleSize(0.05);  
  
  //int max=  test->GetMaximum() + test->GetMaximum()*0.2; 
  //int max=  2000.;
  //h_pt->SetMaximum(max);
  //cout << max << endl;
 
  //h_pt->Draw();
  
  // int histLineColor = kOrange+7;
  //int histFillColor = kOrange-2;
  //float markerSize  = 1.0;

  TLatex latex;
  
  int n_ = 2;
  
  float x1_l = 0.75;
  float y1_l = 0.60;
  
  float dx_l = 0.60;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  
 TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);


  TFile *f1 = new TFile("SumHistsQCD.rot");
  TH1F *B_cnt = static_cast<TH1F*>(f1->Get("kcutscan")->Clone());
  int nbin = B_cnt->GetNbinsX();
  std::vector<double> BCK(nbin);
  for(int i=1;i<nbin+1;i++) BCK[i-1]=B_cnt->GetBinContent(i);
  f1->Close();



  TFile *f2 = new TFile("SumHistsModelB.rot");
  TH1F *S_cnt = static_cast<TH1F*>(f2->Get("kcutscan")->Clone());
  std::vector<double> SGL(nbin);
  for(int i=1;i<nbin+1;i++) SGL[i-1]=S_cnt->GetBinContent(i);
  f2->Close();

  std::vector<std::pair<int,double>> forSort(nbin);

  vector<double> srootb(nbin);
  for(int i=0;i<nbin;i++) {
    if(BCK[i]>0) {
      srootb[i]=SGL[i]/sqrt(BCK[i]);
    } else {
      srootb[i]=0.;
    }
    std::cout<<"cut "<<i<<" "<<SGL[i]<<" "<<BCK[i]<<" "<<srootb[i]<<std::endl;
    forSort[i]=std::make_pair (i,srootb[i]);
  }


  std::sort(forSort.begin(), forSort.end(), pairCompare);
  for(int i=0;i<forSort.size();i++) {
    std::cout<<i<<" "<<forSort[i].first<<" "<<forSort[i].second<<" "<<SGL[forSort[i].first]<<" "<<BCK[forSort[i].first]<<std::endl;
  }

  // find bin with best srootb
  int ipnt=0;
  float amax=0.;
  for(int i=0;i<nbin;i++) {
    if(srootb[i]>amax) {
      ipnt=i;
      amax=srootb[i];
    }
  }
  std::cout<<" optimum cut set is "<<ipnt<<" with srootb of "<<amax<<std::endl;
    
  // put into hist
  
  TH1F* A_pt = new TH1F("A_pt","s root b",nbin,0,nbin);
  for(int i=0;i<nbin;i++) A_pt->AddBinContent(i,srootb[i]);

  A_pt->GetYaxis()->SetTitle("");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  


  A_pt->SetDirectory(0);
  A_pt->SetLineColor(3);
  A_pt->SetLineWidth(3);
  A_pt->SetStats(0);
  A_pt->Draw("");


  
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  //  lgd->Draw();

 
  if (dolog) {
    canv->Print(canvName+"_log.pdf",".pdf");
    canv->Print(canvName+"_log.png",".png");}
  else{ 
    canv->Print(canvName+".pdf",".pdf");
    canv->Print(canvName+".png",".png");}
  return;
}


