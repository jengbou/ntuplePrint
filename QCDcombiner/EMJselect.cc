#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
using std::vector;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
Int_t           fCurrent; //!current Tree number in a TChain                       

// Fixed size dimensions of array or collections stored in the TTree if any.       





void EMJselect(const char* inputfilename,const char* outputfilename,
float HTcut, float alphaMaxcut, float NemfracCut,float CemfracCut,int NemergingCut) {
  // "ntuple.root", "histos.root"
  // suggest cuts 1000., 0.2,0.9,0.9,1


  // read the Tree generated by tree1w and fill two histograms
  // note that we use "new" to create the TFile and TTree objects,
  // to keep them alive after leaving this function.


  TFile *f = new TFile(inputfilename);

  // get histogram of events before trigger
  TH1F* eventCountPreTrigger = static_cast<TH1F*>(f->Get("eventCountPreTrigger/eventCountPreTrigger")->Clone());

  TTree *tt = (TTree*)f->Get("emJetAnalyzer/emJetTree");

  Int_t nVtx, event;
  Float_t met_pt, met_phi;


  vector<int> *jet_index=new vector<int>;
  vector<int> *jet_source=new vector<int>;
  vector<float> *jet_pt = new vector<float>;
  vector<float> *jet_eta = new vector<float>;
  vector<float> *jet_phi = new vector<float>;
  vector<float> *jet_alphaMax = new vector<float>;
  vector<float> *jet_cef = new vector<float>;
  vector<float> *jet_nef = new vector<float>;
  vector<float> *jet_chf = new vector<float>;
  vector<float> *jet_nhf = new vector<float>;
  vector<float> *jet_phf = new vector<float>;
  vector<vector<float> > *track_pt = 0;
  vector<vector<float> > *track_eta = 0;
  vector<vector<int> > *track_source = 0;
  vector<vector<int> > *track_index = 0;
  vector<vector<int> > *track_jet_index = 0;
  vector<vector<int> > *track_vertex_index = 0;
  vector<vector<int> > *track_algo = 0;
  vector<vector<float> > *track_vertex_weight =0;
  vector<vector<float> > *track_ipZ =0;

  //get event count pre trigger



  //for ntuple
  tt->SetBranchAddress("nVtx",&nVtx);
  tt->SetBranchAddress("event",&event);
  tt->SetBranchAddress("met_pt",&met_pt);
  tt->SetBranchAddress("met_phi",&met_phi);
  tt->SetBranchAddress("jet_index",&jet_index);
  tt->SetBranchAddress("jet_source",&jet_source);
  tt->SetBranchAddress("jet_pt",&jet_pt);
  tt->SetBranchAddress("jet_eta",&jet_eta);
  tt->SetBranchAddress("jet_phi",&jet_phi);
  tt->SetBranchAddress("jet_cef",&jet_cef);
  tt->SetBranchAddress("jet_nef",&jet_nef);
  tt->SetBranchAddress("jet_chf",&jet_chf);
  tt->SetBranchAddress("jet_nhf",&jet_nhf);
  tt->SetBranchAddress("jet_phf",&jet_phf);
  tt->SetBranchAddress("jet_alphaMax",&jet_alphaMax);
  tt->SetBranchAddress("track_pt",&track_pt);
  tt->SetBranchAddress("track_eta",&track_eta);
  tt->SetBranchAddress("track_source",&track_source);
  tt->SetBranchAddress("track_index",&track_index);
  tt->SetBranchAddress("track_jet_index",&track_jet_index);
  tt->SetBranchAddress("track_algo",&track_algo);
  tt->SetBranchAddress("track_vertex_index",&track_vertex_index);
  tt->SetBranchAddress("track_vertex_weight",&track_vertex_weight);
  tt->SetBranchAddress("track_ipZ",&track_ipZ);

  // create a histograms
  TH1F *acount = new TH1F("acount","counts",20,0.,20.);
  TH1F *count = new TH1F("count","counts",3,0,3);
  count->SetStats(0);
  count->SetCanExtend(TH1::kAllAxes);

  TH1F *hjetcut = new TH1F("hjetcut","jetcut counts",20,0.,20.);
  TH1F *hjetchf = new TH1F("hjetchf","jet charged hadron fr",20,0.,1.2);
  TH1F *h_nemg = new TH1F("h_nemg","number of emerging jets",20,0.,20.);
  TH1F *hnjet = new TH1F("hnjet","number of jets",20,0.,20.);
  TH1F *hpt = new TH1F("hpt","jet pt distribution",200,0.,1000.);
  TH1F *heta   = new TH1F("heta","jet eta distribution",100,-4.,4.);
  TH1F *heta2   = new TH1F("heta2","jet eta distribution first 4 jets",100,-4.,4.);
  TH1F *halpha   = new TH1F("halpha","jet alphaMax distribution",100,0.,1.5);
  TH1F *H_T      = new TH1F("H_T"," HT distribution before cut", 100,0.,5000.);
  TH1F *H_T2      = new TH1F("H_T2"," HT distribution after cut", 100,0.,5000.);
  TH1F *hbcut_ntrkpt1 = new TH1F("hbcut_ntrkpt1","number tracks pt>1 before cut",10,0.,20.);
  TH1F *hacut_ntrkpt1 = new TH1F("hacut_ntrkpt1","number tracks pt>1 after cut",10,0.,20.);
  TH1F *hbcut_nef = new TH1F("hbcut_nef","neutral em fraction before cut",10,0.,1.2);
  TH1F *hacut_nef = new TH1F("hacut_nef","neutral em fraction after cut",10,0.,1.2);
  TH1F *hbcut_cef = new TH1F("hbcut_cef","charged em fraction before cut",50,0.,1.2);
  TH1F *hacut_cef = new TH1F("hacut_cef","charged em fraction after cut",50,0.,1.2);
  TH1F *hbcut_alphamax = new TH1F("hbcut_alphamax","alphamax before cut",50,0.,1.5);
  TH1F *hacut_alphamax = new TH1F("hacut_alphamax","alphamax after cut",50,0.,1.5);


  //read all entries and fill the histograms
  Int_t nentries = (Int_t)tt->GetEntries();


  // loop over events
  for (Int_t i=0; i<nentries; i++) {
    std::cout<<"event "<<i<<std::endl;
    count->Fill("All",1);  // count number of events
    acount->Fill(0.5);
    tt->GetEntry(i);
    std::cout<<"event number is "<<event<<" number of vertex is "<<nVtx<<std::endl;

    // make some basic plots on all events before any selections
    // jets
    vector<int> jet_ntrkpt1((*jet_index).size());
    hnjet->Fill((*jet_index).size()+0.5);
    for(Int_t j=0; j<(*jet_index).size(); j++) {
      hpt->Fill((*jet_pt)[j]);
      heta->Fill((*jet_eta)[j]);
      hjetchf->Fill((*jet_chf)[j]);
      if(j<4) heta2->Fill((*jet_eta)[j]);
      halpha->Fill((*jet_alphaMax)[j]);
      //      calculate  number of tracks with pt > 1
      jet_ntrkpt1[j]=0;
      vector<float> track_pts = track_pt->at(j);
      vector<int> track_sources = track_source->at(j);
      for (unsigned itrack=0; itrack<track_pts.size(); itrack++) {
	if(track_sources[itrack]==0) {
	  if(track_pts[itrack]>1) jet_ntrkpt1[j]+=1;
	}
      }
     }  // end of loop over jets

    // now start the event selections

    // require at least 4 jets
    if((*jet_index).size()<3) continue;
    count->Fill("4 jets",1);
    acount->Fill(1.5);

    // calculate HT and require it greater than some cut value
    double HT = (*jet_pt)[1]+(*jet_pt)[2]+(*jet_pt)[3]+(*jet_pt)[4];
    H_T->Fill(HT);
    if(HT<HTcut) continue;
    count->Fill("HT",1);
    acount->Fill(2.5);
    H_T2->Fill(HT);

    // do pT cuts on jets  
    bool sel=false;
    if(((*jet_pt)[1]>400)&&((*jet_pt)[2]>200)&&((*jet_pt)[3]>125)&&((*jet_pt)[4]>50)) {
      sel=true;
      count->Fill("jet pt cuts",1);
    acount->Fill(3.5);
    }

      if(!sel) continue;

      //now look and see if any of the jets are emerging

      bool emerging[4];
      emerging[0]=false;emerging[1]=false;emerging[2]=false;emerging[3]=false;
      int nemerging=0;
      for(int ij=0;ij<4;ij++) {
	hjetcut->Fill(0.5);
	hbcut_alphamax->Fill((*jet_alphaMax)[ij]);
	if((*jet_alphaMax)[ij]<alphaMaxcut) {
	  hacut_alphamax->Fill((*jet_alphaMax)[ij]);
	  hjetcut->Fill(1.5);
	  hbcut_nef->Fill((*jet_nef)[ij]);
	  if((*jet_nef)[ij]<NemfracCut) {
	    hacut_nef->Fill((*jet_nef)[ij]);
	    hjetcut->Fill(2.5);
	    hbcut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
	    if(jet_ntrkpt1[ij]>0) {
	      hacut_ntrkpt1->Fill(jet_ntrkpt1[ij]);
	      hjetcut->Fill(3.5);
	      hbcut_cef->Fill((*jet_cef)[ij]);
	      if((*jet_cef)[ij]<CemfracCut) {
	        hacut_cef->Fill((*jet_cef)[ij]);
	        emerging[ij]=true;
	        nemerging+=1.;
		std::cout<<" an emerging jet"<<std::endl;
	      }
	    }
	  }
        }
      }
      h_nemg->Fill(nemerging);

      // require at least N emerging jets
      if(nemerging<NemergingCut) continue;


      count->Fill("emerging",1);
    acount->Fill(4.5);





  }  // end of loop over events
  // We do not close the file. We want to keep the generated
  // histograms we open a browser and the TreeViewer
  if (gROOT->IsBatch()) return;
  //  TCanvas *c1 = new TCanvas("c1","Vertex Plots",200,10,700,500);
  //  hpt->Draw();
  //c1->Modified();
  //c1->Update();
  //  t1->StartViewer();

  TFile myfile(outputfilename,"RECREATE");
  count->LabelsDeflate();
  count->LabelsOption("v");
  //  count->LabelsOption("a");

  eventCountPreTrigger->Write();
  acount->Write();
  count->Write();
  hjetcut->Write();
  hpt->Write();
  hnjet->Write();
  heta->Write();
  heta2->Write();
  halpha->Write();
  H_T->Write();
  H_T2->Write();
  h_nemg->Write();
  hjetchf->Write();
  hbcut_ntrkpt1->Write();
  hacut_ntrkpt1->Write();
  hbcut_nef->Write();
  hacut_nef->Write();
  hbcut_cef->Write();
  hacut_cef->Write();
  hbcut_alphamax->Write();
  hacut_alphamax->Write();


  tt->ResetBranchAddresses();
  delete jet_index;
  delete jet_source;
  delete jet_pt;
  delete jet_eta;
  delete jet_phi;
  delete jet_alphaMax;
  delete track_pt;
  delete track_eta;
  delete track_source;
  delete track_index;
  delete track_jet_index;
  delete track_algo;
  delete track_vertex_index;
  delete track_vertex_weight;
  delete track_ipZ;


  //In the browser, click on "ROOT Files", then on "tree1.root"
  //You can click on the histogram icons in the right panel to draw
  //them in the TreeViewer, follow the instructions in the Help.
}
