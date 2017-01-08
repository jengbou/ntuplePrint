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

//TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
//Int_t           fCurrent; //!current Tree number in a TChain                       

bool EMJscanFirst=true;





vector<int> EMJscan(const char* inputfilename,
	     float HTcutmin,int NHTcut, float HTcutSS,
	     float pt1cutmin, int Npt1cut, float pt1cutSS,
               float pt2cutmin,  int Npt2cut,float pt2cutSS,
               float pt3cutmin,  int Npt3cut,float pt3cutSS,
               float pt4cutmin, int Npt4cut,float pt4cutSS,
	     int NemergingCutmin, int NNemergingCut, 
		    float alphaMaxcut, float NemfracCut,float CemfracCut,int ntrk1cut) {

 
  int iicut = NHTcut*Npt1cut*Npt2cut*Npt3cut*Npt4cut*NNemergingCut;
  vector<int> npass(iicut);
  for(int i=0;i<sizeof(npass);i++) npass[i]=0;

  TFile *f = new TFile(inputfilename);

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


  //read all entries and fill the histograms
  Int_t nentries = (Int_t)tt->GetEntries();


  // loop over events
  for (Int_t i=0; i<nentries; i++) {
    //    if(i%100 == 0) std::cout<<"event "<<i<<std::endl;
    tt->GetEntry(i);
    //    std::cout<<"event number is "<<event<<" number of vertex is "<<nVtx<<std::endl;

    // jets
    vector<int> jet_ntrkpt1((*jet_index).size());
    for(Int_t j=0; j<(*jet_index).size(); j++) {
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


    double HT = (*jet_pt)[1]+(*jet_pt)[2]+(*jet_pt)[3]+(*jet_pt)[4];
    // now start the event selections

      //now look and see if any of the jets are emerging

      bool emerging[4];
      emerging[0]=false;emerging[1]=false;emerging[2]=false;emerging[3]=false;
      int nemerging=0;
      for(int ij=0;ij<4;ij++) {
	if((*jet_alphaMax)[ij]<alphaMaxcut) {
	  if((*jet_nef)[ij]<NemfracCut) {
	    if(jet_ntrkpt1[ij]>ntrk1cut) {
	      if((*jet_cef)[ij]<CemfracCut) {
	        emerging[ij]=true;
	        nemerging+=1.;
		//		std::cout<<" an emerging jet"<<std::endl;
	      }
	    }
	  }
        }
      }



    // require at least 4 jets
    if((*jet_index).size()<3) continue;


    int icut=0;
    float HTcut,pt1cut,pt2cut,pt3cut,pt4cut;
    int NemergingCut;
    for(int iht=0;iht<NHTcut;iht++) {
      HTcut = HTcutmin+iht*HTcutSS;
      for(int ipt1=0;ipt1<Npt1cut;ipt1++) {
	pt1cut=pt1cutmin+ipt1*pt1cutSS;
        for(int ipt2=0;ipt2<Npt1cut;ipt2++) {
	  pt2cut=pt2cutmin+ipt2*pt2cutSS;
          for(int ipt3=0;ipt3<Npt1cut;ipt3++) {
	    pt3cut=pt3cutmin+ipt3*pt3cutSS;
            for(int ipt4=0;ipt4<Npt1cut;ipt4++) {
	      pt4cut=pt4cutmin+ipt4*pt4cutSS;
	      for(int inem=0;inem<NNemergingCut;inem++) {
	        NemergingCut=NemergingCutmin+inem;
        
		icut = icut+1;          
		if(EMJscanFirst) {
		  std::cout<<"icut "<<icut<<" corresponds to "<<
		  "HT cut of "<<HTcut<<
		  "pt1 cut of "<<pt1cut<<
		  "pt2 cut of "<<pt2cut<<
		  "pt3 cut of "<<pt3cut<<
		  "pt4 cut of "<<pt4cut<<
		  "nemerging cut of "<<NemergingCut<<
		  std::endl;
		  if(icut==iicut) EMJscanFirst= false;
		}
	          if(HT>HTcut) {
	            if((*jet_pt)[1]>pt1cut) {
	              if((*jet_pt)[2]>pt2cut) {
	                if((*jet_pt)[3]>pt3cut) {
	                  if((*jet_pt)[4]>pt4cut) {
	                    if(nemerging>NemergingCut) {
                              npass[icut]+=1;
	                    }
	                  }
	                }
	              }
	            }
	          }
	        }
              }
            }
          }
        }
      }


  }  // end of loop over events


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

  f->Close();
  

  return npass;
}
