#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <map>

#include "vector"
using std::vector;
#include "algorithm"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "EMJselect.h"
#include "EMJ16003.h"
#include "QCDhists.h"
#include "EMJbkg.h"

float pilecut=2.5;//cm
//float pilecut2=pilecut;//cm
float pvztrkcut=0.01;//cm
float ipXYsigcut=2.0;
//float METcut=0.0;

double minJetPt = 100.0;
const int varType = 3;// 1: alpha; 2: alpha2Dsig; 3 (default): alpha3Dsig

//btagging related
float bTagWP=0.5426;//0.90095;//0.627516656;//0.5426;//0.8484;//0.9535;//CSVv2T defined here:  https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
int bUnfType=1;//1: loose (0.5426); 2: medium (0.8484); 3: tight (0.9535)
int ptType=1;//1: low pT; 2: high pT
// btagging related

float DeltaR(float eta1, float phi1, float eta2, float phi2) {

    float dR=0.;
    float deta = std::fabs(eta1-eta2);
    float dphi = std::fabs(phi1-phi2);
    if(dphi>3.14159) dphi = 2.*3.14159-dphi;
    dR=std::sqrt(deta*deta+dphi*dphi);

    return dR;
}

int EMJbkg(bool otfile, bool hasPre, std::string cutsetname, const char* inputfilename, const char* outputfilename,
           float HTcut, float pt1cut, float pt2cut, float pt3cut, float pt4cut, float jetacut, float alphaMaxcut,
           float maxIPcut, float NemfracCut, float CemfracCut,int ntrk1cut, int ntagType, float PUdzCut,
           float sigzCut, float METcut, bool blind, bool isData) {
    //const int ntagType = 22; //0:>=1 tag; 1: ==1 tag; 2: == 2 tag (flav dep); 21: == 2 tag (T&P); 22: == 2 tag (flav dep; data unfolded); 3: == 3 tag; else: >= NemergingCut tag

//     std::cout << "Cuts: HT pt1 pt2 pt3 pt4 MET ntagType PUdz 3Dsig medIP alpha3Dcut"
//               << std::endl
//               << std::setw(12) << " "
//               << HTcut << " "
//               << pt1cut << " "
//               << pt2cut << " "
//               << pt3cut << " "
//               << pt4cut << " "
//               << METcut << " "
//               << ntagType << " "
//               << PUdzCut << " "
//               << sigzCut << " "
//               << maxIPcut << " "
//               << alphaMaxcut << std::endl;
    const std::string cutname=cutsetname;
    const std::string histcutname="Cutset"+cutsetname;

    int NemergingCut = 0;
    pilecut=PUdzCut;

    if (minJetPt>pt4cut) minJetPt = pt4cut;
    QCDhists* qcdtools=0;

    int npass=0;

    TFile *f = new TFile(inputfilename);
    if(f->IsZombie()) {
        std::cout << "File: " << inputfilename << " does not exist." << std::endl;
        return 0;
    }

    //histogram map
    std::map<TString,TH1*> hs1D; // map of usual 1D histogram
    std::map<TString,TH2*> hs2D; // map of usual 2D histogram

    // get histogram of events before trigger
    if(hasPre) {
        if(otfile) hs1D[TString("eventCountPreTrigger_"+histcutname)] = static_cast<TH1F*>(f->Get("eventCountPreTrigger/eventCountPreTrigger")->Clone(TString("eventCountPreTrigger_"+histcutname)));
    } else {
        if(otfile) hs1D[TString("eventCountPreTrigger_"+histcutname)] = new TH1F(TString("eventCountPreTrigger_"+histcutname),"eventCountPreTrigger",2,0.,2.);
    }


//     char *hnames[9]={(char*)"hjptaFR",(char*)"hetaaFR",(char*)"hntrkFR",
//                      (char*)"hmedipXYSigFR",(char*)"hlogmedipXYSigFR",(char*)"hmedtheta2DFR",
//                      (char*)"hlogmedtheta2DFR",(char*)"hmassFR",(char*)"htheta2DvipXYSigFR"};

    std::cout << "CSV batgger cut used = " << bTagWP << std::endl;

    TTree *tt = (TTree*)f->Get("emJetAnalyzer/emJetTree");

    Int_t nVtx, event, lumi, run, nTrueInt, nTracks;//, pv_indexInColl;//analysis_20170523_v0
    Float_t met_pt, met_phi;//, pv_index;
    bool hltTrig1n, hltTrig1d, hltTrig2n, hltTrig2d, hltTrig3n, hltTrig3d;
    //float pv_z;//analysis_20170523_v0

    vector<int> *pv_index=0;//analysis_20170609_v0
    vector<float> *pv_zv=0;//analysis_20170609_v0

    vector<int> *jet_index=0;
    vector<int> *jet_source=0;
    vector<float> *jet_pt = 0;
    vector<float> *jet_eta = 0;
    vector<float> *jet_phi = 0;
    vector<float> *jet_alphaMax = 0;
    vector<float> *jet_cef = 0;
    vector<float> *jet_nef = 0;
    vector<float> *jet_chf = 0;
    vector<float> *jet_nhf = 0;
    vector<float> *jet_theta2D = 0;
    //  vector<float> *jet_phf = 0;
    vector<float> *jet_csv = 0;
    vector<vector<float> > *track_pt = 0;
    vector<vector<float> > *track_eta = 0;
    vector<vector<int> > *track_source = 0;
    vector<vector<int> > *track_index = 0;
    vector<vector<int> > *track_jet_index = 0;
    vector<vector<int> > *track_vertex_index = 0;
    vector<vector<int> > *track_algo = 0;
    vector<vector<int> > *track_quality = 0;
    vector<vector<float> > *track_vertex_weight =0;
    vector<vector<float> > *track_pvWeight = 0;
    vector<vector<float> > *track_ipZ =0;
    vector<vector<float> > *track_ipXY = 0;
    vector<vector<float> > *track_ipXYSig = 0;
    vector<vector<float> > *track_ref_z =0;//analysis_20170523_v0
    vector<vector<int> > *track_nMissInnerHits = 0;
    vector<vector<int> > *track_nMissInnerPxlLayers = 0;
    vector<vector<int> > *track_nPxlLayers = 0;
    vector<vector<int> > *track_nHits = 0;

    /*
      vector<int> *vertex_index=0;
      vector<int> *vertex_source=0;
      vector<float> *vertex_z =0;
    */

    // gen particles
    vector<int> *gp_index = new vector<int>;
    vector<int> *gp_pdgId = new vector<int>;
    vector<float> *gp_pt = new vector<float>;
    vector<float> *gp_eta = new vector<float>;
    vector<float> *gp_phi = new vector<float>;
    vector<int> *gp_charge = new vector<int>;
    vector<int> *gp_status = new vector<int>;
    vector<float> *gp_vx = new vector<float>;
    vector<float> *gp_vy = new vector<float>;
    vector<float> *gp_vz = new vector<float>;
    //get event count pre trigger


    // gen particles
    tt->SetBranchAddress("gp_index",&gp_index);
    tt->SetBranchAddress("gp_pdgId",&gp_pdgId);
    tt->SetBranchAddress("gp_pt",&gp_pt);
    tt->SetBranchAddress("gp_eta",&gp_eta);
    tt->SetBranchAddress("gp_phi",&gp_phi);
    tt->SetBranchAddress("gp_charge",&gp_charge);
    tt->SetBranchAddress("gp_status",&gp_status);
    tt->SetBranchAddress("gp_vx",&gp_vx);
    tt->SetBranchAddress("gp_vy",&gp_vy);
    tt->SetBranchAddress("gp_vz",&gp_vz);

    //for ntuple
    //tt->SetBranchAddress("pv_indexInColl",&pv_indexInColl);//analysis_20170523_v0
    tt->SetBranchAddress("pv_z",&pv_zv);
    tt->SetBranchAddress("pv_index",&pv_index);//analysis_20170609_v0
    tt->SetBranchAddress("nVtx",&nVtx);
    tt->SetBranchAddress("nTrueInt",&nTrueInt);
    tt->SetBranchAddress("nTracks",&nTracks);
    tt->SetBranchAddress("event",&event);
    tt->SetBranchAddress("lumi",&lumi);
    tt->SetBranchAddress("run",&run);
    tt->SetBranchAddress("met_pt",&met_pt);
    tt->SetBranchAddress("met_phi",&met_phi);
    tt->SetBranchAddress("jet_index",&jet_index);
    tt->SetBranchAddress("jet_source",&jet_source);
    //tt->SetBranchAddress("jet_ptRaw",&jet_pt);//analysis_20170523_v0
    tt->SetBranchAddress("jet_pt",&jet_pt);
    tt->SetBranchAddress("jet_eta",&jet_eta);
    tt->SetBranchAddress("jet_phi",&jet_phi);
    tt->SetBranchAddress("jet_cef",&jet_cef);
    tt->SetBranchAddress("jet_nef",&jet_nef);
    tt->SetBranchAddress("jet_chf",&jet_chf);
    tt->SetBranchAddress("jet_nhf",&jet_nhf);
    tt->SetBranchAddress("jet_theta2D",&jet_theta2D);//median of assoc. trk. theta2Ds.
    //  tt->SetBranchAddress("jet_phf",&jet_phf);
    tt->SetBranchAddress("jet_csv",&jet_csv);
    tt->SetBranchAddress("jet_alphaMax",&jet_alphaMax);
    tt->SetBranchAddress("track_pt",&track_pt);
    tt->SetBranchAddress("track_eta",&track_eta);
    tt->SetBranchAddress("track_source",&track_source);
    tt->SetBranchAddress("track_index",&track_index);
    tt->SetBranchAddress("track_jet_index",&track_jet_index);
    tt->SetBranchAddress("track_algo",&track_algo);
    tt->SetBranchAddress("track_quality",&track_quality);
    tt->SetBranchAddress("track_vertex_index",&track_vertex_index);
    tt->SetBranchAddress("track_vertex_weight",&track_vertex_weight);
    tt->SetBranchAddress("track_pvWeight",&track_pvWeight);
    tt->SetBranchAddress("track_ipXY",&track_ipXY);
    tt->SetBranchAddress("track_ipXYSig",&track_ipXYSig);
    tt->SetBranchAddress("track_ref_z",&track_ref_z);//analysis_20170523_v0
    tt->SetBranchAddress("track_nMissInnerHits",&track_nMissInnerHits);
    tt->SetBranchAddress("track_nMissInnerPxlLayers",&track_nMissInnerPxlLayers);
    tt->SetBranchAddress("track_nPxlLayers",&track_nPxlLayers);
    tt->SetBranchAddress("track_nHits",&track_nHits);
    tt->SetBranchAddress("track_ipZ",&track_ipZ);
    tt->SetBranchAddress("HLT_HT400",&hltTrig1d);
    tt->SetBranchAddress("HLT_HT500",&hltTrig1n);
    tt->SetBranchAddress("HLT_HT250",&hltTrig2d);
    tt->SetBranchAddress("HLT_HT350",&hltTrig2n);
    tt->SetBranchAddress("HLT_PFHT600",&hltTrig3d);
    //tt->SetBranchAddress("HLT_PFHT800",&hltTrig3n);//76, 80
    tt->SetBranchAddress("HLT_PFHT900",&hltTrig3n);//G,H

    /*
      tt->SetBranchAddress("vertex_index",&vertex_index);
      tt->SetBranchAddress("vertex_source",&vertex_source);
      tt->SetBranchAddress("vertex_z",&vertex_z);
    */

    if(otfile) {
        TH1::SetDefaultSumw2();

        hs1D[TString("acount_"+histcutname)] = new TH1F(TString("acount_"+histcutname),"counts",20,0.,20.);
        hs1D[TString("count_"+histcutname)] = new TH1F(TString("count_"+histcutname),"counts",3,0.,3);

        hs1D[TString("count_"+histcutname)]->SetStats(0);
        hs1D[TString("count_"+histcutname)]->SetCanExtend(TH1::kAllAxes);
        hs1D[TString("count_"+histcutname)]->Fill("All",0);
        hs1D[TString("count_"+histcutname)]->Fill("filter",0);
        hs1D[TString("count_"+histcutname)]->Fill("4 jets",0);
        hs1D[TString("count_"+histcutname)]->Fill("HT",0);
        hs1D[TString("count_"+histcutname)]->Fill("jet pt1",0);
        hs1D[TString("count_"+histcutname)]->Fill("jet pt2",0);
        hs1D[TString("count_"+histcutname)]->Fill("jet pt3",0);
        hs1D[TString("count_"+histcutname)]->Fill("jet pt4",0);
        hs1D[TString("count_"+histcutname)]->Fill("MET",0);
        hs1D[TString("count_"+histcutname)]->Fill("emerging",0);
        hs1D[TString("count_"+histcutname)]->Fill("almostemerging",0);
        hs1D[TString("count_"+histcutname)]->Fill("FakeRateBKG",0);


        //
        hs1D[TString("hjptaSR_"+histcutname)] = new TH1F(TString("hjptaSR_"+histcutname)," pT of emergng jets (SR)",100,0.,1000.);

        hs1D[TString("hjptaFR_"+histcutname)] = new TH1F(TString("hjptaFR_"+histcutname)," pT of emergng jets (FR)",100,0.,1000.);

        hs1D[TString("hetaaSR_"+histcutname)] = new TH1F(TString("hetaaSR_"+histcutname)," eta of emergng jets (SR)",100,-4.,4.);

        hs1D[TString("hetaaFR_"+histcutname)] = new TH1F(TString("hetaaFR_"+histcutname)," eta of emergng jets (FR)",100,-4.,4.);

        hs1D[TString("hntrkSR_"+histcutname)] = new TH1F(TString("hntrkSR_"+histcutname),"number tracks pt>1 (SR)",50,0.,50.);

        hs1D[TString("hntrkFR_"+histcutname)] = new TH1F(TString("hntrkFR_"+histcutname),"number tracks pt>1 (FR)",50,0.,50.);

        hs1D[TString("hmedipXYSigSR_"+histcutname)] = new TH1F(TString("hmedipXYSigSR_"+histcutname),"median ip_{sig} emerging jets (SR)",1000,0.,100.);

        hs1D[TString("hmedipXYSigFR_"+histcutname)] = new TH1F(TString("hmedipXYSigFR_"+histcutname),"median ip_{sig} emerging jets (FR)",1000,0.,100.);

        hs1D[TString("hlogmedipXYSigSR_"+histcutname)] = new TH1F(TString("hlogmedipXYSigSR_"+histcutname),
                                                                  "median log_{10} ip_{sig} emerging jets (SR)",1000,-1.,4.);
        hs1D[TString("hlogmedipXYSigFR_"+histcutname)] = new TH1F(TString("hlogmedipXYSigFR_"+histcutname),
                                                                  "median log_{10} ip_{sig} emerging jets (FR)",1000,-1.,4.);

        hs1D[TString("hmedtheta2DSR_"+histcutname)] = new TH1F(TString("hmedtheta2DSR_"+histcutname),
                                                               "median #hat{#Theta_{2D}} emerging jets (SR)",200,0.,0.4);
        hs1D[TString("hmedtheta2DFR_"+histcutname)] = new TH1F(TString("hmedtheta2DFR_"+histcutname),
                                                               "median #hat{#Theta_{2D}} emerging jets (FR)",200,0.,0.4);

        hs1D[TString("hlogmedtheta2DSR_"+histcutname)] = new TH1F(TString("hlogmedtheta2DSR_"+histcutname),
                                                                  "median log_{10} #hat{#Theta_{2D}} emerging jets (SR)",1000,-3.5,0.5);
        hs1D[TString("hlogmedtheta2DFR_"+histcutname)] = new TH1F(TString("hlogmedtheta2DFR_"+histcutname),
                                                                  "median log_{10} #hat{#Theta_{2D}} emerging jets (FR)",1000,-3.5,0.5);

        hs1D[TString("hmassSR_"+histcutname)] = new TH1F(TString("hmassSR_"+histcutname),"mass emerging and non pairs (SR)",500,0.,5000.);

        hs1D[TString("hmassFR_"+histcutname)] = new TH1F(TString("hmassFR_"+histcutname),"mass emerging and non pairs (FR)",500,0.,5000.);

        hs2D[TString("htheta2DvipXYSigSR_"+histcutname)] = new TH2F(TString("htheta2DvipXYSigSR_"+histcutname),
                                                                    " #hat{#Theta}_{2D} vs. #hat{IP}^{2D}_{Sig} plot (SR)",100,0.,0.4,100,0.,10.0);
        hs2D[TString("htheta2DvipXYSigFR_"+histcutname)] = new TH2F(TString("htheta2DvipXYSigFR_"+histcutname),
                                                                    " #hat{#Theta}_{2D} vs. #hat{IP}^{2D}_{Sig} plot (FR)",100,0.,0.4,100,0.,10.0);


        hs1D[TString("h_nemg_"+histcutname)] = new TH1F(TString("h_nemg_"+histcutname),"number of emerging jets",20,0.,20.);

        hs1D[TString("h_nemgSR_"+histcutname)] = new TH1F(TString("h_nemgSR_"+histcutname),"number of emerging jets (SR)",20,0.,20.);

        hs1D[TString("h_nemgFR_"+histcutname)] = new TH1F(TString("h_nemgFR_"+histcutname),"number of emerging jets (FR)",20,0.,20.);

        hs1D[TString("h_nalemg_"+histcutname)] = new TH1F(TString("h_alnemg_"+histcutname),"number of almostemerging jets",20,0.,20.);

        hs1D[TString("hnjet_"+histcutname)] = new TH1F(TString("hnjet_"+histcutname),"number of jets",20,0.,20.);

        hs1D[TString("hnjetSR_"+histcutname)] = new TH1F(TString("hnjetSR_"+histcutname),"number of jets (SR)",20,0.,20.);

        hs1D[TString("hnjetFR_"+histcutname)] = new TH1F(TString("hnjetFR_"+histcutname),"number of jets (FR)",20,0.,20.);

        hs1D[TString("halpha_"+histcutname)] = new TH1F(TString("halpha_"+histcutname),"jet alpha distribution",1000,0.,1.0);

        hs1D[TString("halphaZero_"+histcutname)] = new TH1F(TString("halphaZero_"+histcutname),"jet alpha==0 distribution",1000,-0.1,0.1);

        hs1D[TString("H_TSR_"+histcutname)] = new TH1F(TString("H_TSR_"+histcutname)," HT distribution at end (SR)", 100,0.,5000.);

        hs1D[TString("H_TFR_"+histcutname)] = new TH1F(TString("H_TFR_"+histcutname)," HT distribution at end (FR)", 100,0.,5000.);

    }

    //read all entries and fill the histograms
    Int_t nentries = (Int_t)tt->GetEntries();


    double np0 = 0;
    double np1 = 0;
    double np2 = 0;
    double np3 = 0;
    int nbTaggedTotAll = 0;//all jets
    int nNotbTagTotAll = 0;//all jets
    int nbTaggedTot = 0;//leading four jets
    int nNotbTagTot = 0;//leading four jets


    switch (ntagType) {
    case 0 : {NemergingCut=0; break;}
    case 1 : {NemergingCut=1; break;}
    case 12: {NemergingCut=1; break;}
    case 2 : {NemergingCut=2; break;}
    case 21: {NemergingCut=2; break;}
    case 22: {NemergingCut=2; break;}
    case 3 : {NemergingCut=3; break;}
    }


    //======================
    // Loop over events
    //======================
    for (Int_t i=0; i<nentries; i++) {
        int nbTaggedAll = 0;//# of all b-tagged jets
        int nbTagged = 0;//# of b-tagged jets of leading four jets
        int nZeroAlpha = 0;
        //std::cout<<"***event "<<event<<std::endl;
 
        if(!hasPre && otfile) hs1D[TString("eventCountPreTrigger_"+histcutname)]->Fill(1);
    
        if(otfile) hs1D[TString("count_"+histcutname)]->Fill("All",1);  // count number of events
        if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(0); // events passing JetFilter
        tt->GetEntry(i);

        // PV failure removal
        if (pv_index->at(0) != 0) continue;

        // PVZ cut
        float pv_z = pv_zv->at(0);
        if(fabs(pv_z)>15) continue;//analysis_20170523_v0 onward

        // Counting variables for PVTrackFraction
        int nTotGoodTrk = 0;
        int nGoodPVzTrk = 0;

        const int NNNjet = jet_index->size();
        // Note: we're only interested in events with at least four jets
        if(NNNjet<4) continue;

        // jets
        vector<double> jet_fnonPU(NNNjet);
        vector<double> jet_fmaxtrkpt(NNNjet);
        vector<double> jet_fsmallIPsig(NNNjet);
        vector<int> jet_ntrkpt1(NNNjet);
        vector<float> jet_alpha(NNNjet);
        vector<float> jet_meanip(NNNjet);
        vector<float> jet_meanipsig(NNNjet);
        vector<float> jet_logmeanipsig(NNNjet);
        vector<float> jet_medip(NNNjet);
        vector<float> jet_medipsig(NNNjet);
        vector<float> jet_logmedipsig(NNNjet);
        vector<float> jet_medtheta2D(NNNjet);
        vector<float> jet_logmedtheta2D(NNNjet);
        vector<float> r0(NNNjet);
        vector<float> r0sig(NNNjet);
        vector<float> r1(NNNjet);
        vector<int> jntrack(NNNjet);
        vector<float> jet_e(NNNjet);
        vector<float> jet_theta(NNNjet);
        vector<float> jet_px(NNNjet);
        vector<float> jet_py(NNNjet);
        vector<float> jet_pz(NNNjet);

        vector<int> jet_pid_maxEt(NNNjet);
        vector<float> jet_maxET_part(NNNjet);

        if(otfile) hs1D[TString("hnjet_"+histcutname)]->Fill(NNNjet);

        //======================
        // First loop over jets
        //======================
        for(Int_t j=0; j<NNNjet; j++) {
            //      std::cout<<"jet j = "<<j<<std::endl;
            jet_theta[j]=2.*atan(exp(-jet_eta->at(j)));
            jet_e[j]=jet_pt->at(j)/sin(jet_theta[j]);
            jet_px[j]=jet_pt->at(j)*cos(jet_phi->at(j));
            jet_py[j]=jet_pt->at(j)*sin(jet_phi->at(j));
            jet_pz[j]=jet_pt->at(j)/tan(jet_theta[j]);

            jet_ntrkpt1[j]=0;
            jet_meanip[j]=-1.;
            jet_logmeanipsig[j]=-1.;
            jet_medip[j]=-1.;
            jet_medipsig[j]=-1.;
            jet_logmedipsig[j]=-1.;
            jet_medtheta2D[j]=-1.;
            jet_logmedtheta2D[j]=-3.5;
            jet_alpha[j]=-1.;
            jet_fnonPU[j]=0.;
            jet_fmaxtrkpt[j]=0.;
            jet_fsmallIPsig[j]=0.;

            if(r0.size()>0) {r0[j]=-1.;r0sig[j]=-1.;}
            if(r1.size()>0) r1[j]=-1.;
            vector<float> track_pts = track_pt->at(j);
            vector<int> track_sources = track_source->at(j);
            vector<int> track_qualitys = track_quality->at(j);
            //vector<float> track_vertex_weights = track_vertex_weight->at(j);
            vector<float> track_pvWeights = track_pvWeight->at(j);
            vector<float> track_ipXYs = track_ipXY->at(j);
            vector<float> track_ipXYSigs = track_ipXYSig->at(j);
            vector<float> sort_ip(track_pts.size());
            vector<float> sort_ipsig(track_pts.size());
            for(uint it=0;it<track_pts.size();it++) sort_ip[it]=0;
            for(uint it=0;it<track_pts.size();it++) sort_ipsig[it]=0;
            vector<float> jet_trkip;
            vector<float> jet_trkipsig;
            vector<float> track_ref_zs = track_ref_z->at(j);//analysis_20170523_v0

            // jet_alpha
            if (varType==1) jet_alpha[j] = qcdtools->GetAlpha(track_pts,track_sources,track_qualitys,track_pvWeights,track_ref_zs,pv_z,PUdzCut);
            else if (varType==2) {
                jet_alpha[j] = qcdtools->GetAlpha2Dsig(track_pts,track_sources,track_qualitys,track_ipXYSigs,track_ref_zs,pv_z,PUdzCut);
            }
            else if (varType==3) {
                jet_alpha[j] = qcdtools->GetAlpha3Dsig(track_pts,track_sources,track_qualitys,track_ipXYSigs,track_ref_zs,pv_z,PUdzCut,sigzCut);
            }
            else {
                std::cout << "[Warning] varType not defined! Default alpha3Dsig is used!" << std::endl;
                jet_alpha[j] = qcdtools->GetAlpha3Dsig(track_pts,track_sources,track_qualitys,track_ipXYSigs,track_ref_zs,pv_z);
            }
            //std::cout << "alphaMax = " << jet_alphaMax->at(j) << "; alpha new = " << jet_alpha[j] << " [here]" << std::endl;

            if (jet_alpha[j]==0) nZeroAlpha++;
            jntrack[j]=0;
            double sumptallnz=0.;
            double sumnonPU=0.;
            double sumsmallIPsig=0.;
            double ptmaxtrk=0.;

            //===================================
            // Loop over tracks of first jet loop
            //===================================
            for (unsigned itrack=0; itrack<track_pts.size(); itrack++) {
                if(track_sources[itrack]==0 && ((track_qualitys[itrack] & 4) > 0)) {
                    sumptallnz+=track_pts[itrack];

                    // PVTrackFraction use tracks without pilecut
                    nTotGoodTrk++;//Count good tracks for PVTrackFraction
                    if (fabs(pv_z-track_ref_zs[itrack])<pvztrkcut) nGoodPVzTrk++;//Count good tracks passing pvztkcut for PVTrackFraction

                    // pileup cut
                    if (fabs(pv_z-track_ref_zs[itrack])>pilecut) continue;// remove tracks with exceedingly large z
                    if (fabs(pv_z-track_ref_zs[itrack])<PUdzCut) sumnonPU+=track_pts[itrack];
                    if (fabs(track_ipXYSigs[itrack])<ipXYsigcut) sumsmallIPsig+=track_pts[itrack];

                    if(track_pts[itrack]>ptmaxtrk) {
                        ptmaxtrk=track_pts[itrack];
                    }
                    sort_ip[jntrack[j]]=fabs(track_ipXYs[itrack]);
                    sort_ipsig[jntrack[j]]=fabs(track_ipXYSigs[itrack]);
                    //  std::cout<<"track vertex weight is "<<track_vertex_weights[itrack]<<std::endl;
                    if(track_pts[itrack]>1) jet_ntrkpt1[j]+=1;
                    //  std::cout<<" track "<<itrack<<" ip "<<track_ipXYs[itrack]<<" mean ip "<<jet_meanip[j]<<std::endl;
                    jet_meanip[j]=jet_meanip[j]+fabs(track_ipXYs[itrack]);
                    jet_meanipsig[j]=jet_meanipsig[j]+fabs(track_ipXYSigs[itrack]);
                    jet_logmeanipsig[j]=jet_logmeanipsig[j]+fabs(track_ipXYSigs[itrack]);
                    jet_trkip.push_back(fabs(track_ipXYs[itrack]));
                    jet_trkipsig.push_back(fabs(track_ipXYSigs[itrack]));
                    jntrack[j]++;
                }
            }
            //=======================================
            // End loop over tracks of first jet loop
            //=======================================

            float atmp = jntrack[j];
            if(jntrack[j]>0) {
                jet_meanip[j]=jet_meanip[j]/atmp;
                jet_logmeanipsig[j]=jet_logmeanipsig[j]/atmp;
                jet_logmeanipsig[j]=log10(jet_logmeanipsig[j]);

                // median
                jet_medip[j] = CalcMedian(jet_trkip);
                jet_medipsig[j] = CalcMedian(jet_trkipsig);
                jet_logmedipsig[j] = log10(jet_medipsig[j]);
                jet_medtheta2D[j] = jet_theta2D->at(j);
                jet_logmedtheta2D[j] = (jet_medtheta2D[j]==-1 ? -3.5 : log10(jet_theta2D->at(j)));
            }
            if(sumptallnz>0) jet_fnonPU[j]=sumnonPU/sumptallnz;//analysis_20170523_v0
            if(sumnonPU>0) jet_fsmallIPsig[j]=sumsmallIPsig/sumnonPU;//analysis_20170523_v0
            jet_fmaxtrkpt[j]=ptmaxtrk/jet_pt->at(j);//max trk pt/jet pt//analysis_20170523_v0

            std::sort(sort_ip.begin(), sort_ip.end(), std::greater<float>());
            if(sort_ip.size()>0) r0[j]=sort_ip[0];
            if(sort_ip.size()>1) r1[j]=sort_ip[1];
            std::sort(sort_ipsig.begin(), sort_ipsig.end(), std::greater<float>());
            if(sort_ipsig.size()>0) r0sig[j]=sort_ipsig[0];
            //      std::cout<<"mean max are "<<jet_meanip[j]<<" "<<r0[j]<<std::endl;

            //===================================
            // Jet flavor id
            //===================================
            jet_pid_maxEt[j]=0;
            jet_maxET_part[j]=0;

            if (!isData) {
                // calculate some gen particle information for jet
                int NNNgp = gp_index->size();
                int igenmax=-1;
                float etgenmax=0.;
                for(Int_t igen=1; igen<NNNgp; igen++) {
                    if((abs(gp_pdgId->at(igen))<6)||(abs(gp_pdgId->at(igen))==21)) {  // quark or gluon
                        if(DeltaR(jet_eta->at(j),jet_phi->at(j),gp_eta->at(igen),gp_phi->at(igen))<0.4) {
                            if(gp_pt->at(igen)>etgenmax) {
                                igenmax=igen;
                                etgenmax=gp_pt->at(igen);
                            }
                        }
                    }
                }
                // fix glue to bbbar
                float igenmax2=-1;
                float etgenmax2=0.;
                if(igenmax>0) {
                    if(abs(gp_pdgId->at(igenmax))==21) {
                        for(Int_t igen=1; igen<NNNgp; igen++) {
                            if((abs(gp_pdgId->at(igen))==5)&&(gp_pt->at(igen)>10.)) {  // b
                                if(DeltaR(jet_eta->at(j),jet_phi->at(j),gp_eta->at(igen),gp_phi->at(igen))<0.4) {
                                    if(gp_pt->at(igen)>etgenmax2) {
                                        igenmax2=1;
                                        etgenmax2=gp_pt->at(igen);
                                    }
                                }
                            }
                        }
                    }
                }

                if(igenmax>-1) {
                    int ipid = gp_pdgId->at(igenmax);
                    if(abs(ipid)<6) {
                        jet_pid_maxEt[j]=gp_pdgId->at(igenmax);
                        jet_maxET_part[j] = etgenmax;
                    } else {
                        if(igenmax2==-1) {
                            jet_pid_maxEt[j]=7;
                            jet_maxET_part[j] = etgenmax;
                        } else {//g->bb
                            jet_pid_maxEt[j]=8;
                            jet_maxET_part[j] = etgenmax;

                        }
                    }
                } // end calculate some gen particle information for jet
            }//end of MC truth only block
            if(jet_csv->at(j) >= bTagWP) {
                if (isData) jet_pid_maxEt[j]=5;
                nbTaggedTotAll++;
                nbTaggedAll++;
                //std::cout << "Jet b-tagged." << std::endl;
            }
            else {
                if (isData) jet_pid_maxEt[j]=0;
                nNotbTagTotAll++;
                //std::cout << "Jet not b-tagged." << std::endl;
            }
        }
        //===================================
        // End of first loop over jets
        //===================================


        float PVTrackFraction = (float)nGoodPVzTrk/(float)nTotGoodTrk;
        //std::cout <<"PVTrackFraction = " << PVTrackFraction << std::endl;

        //now see which jets are emerging
        //    std::cout<<" in event "<<event<<" number of jets is "<<NNNjet<<std::endl;
        vector<bool> emerging(NNNjet);
        vector<bool> almostemerging(NNNjet);
        vector<bool> basicjet(NNNjet);
        for( int i=0;i<NNNjet;i++) {
            emerging[i]=false;
            almostemerging[i]=false;
            basicjet[i]=false;
        }
        int nemerging=0;
        int nalmostemerging=0;


        //===================================
        // Second loop over jets
        // Jet selection
        //===================================
        for(int ij=0;ij<NNNjet;ij++) {
            vector<int> track_qualitys = track_quality->at(ij);
            vector<float> track_ipXYs = track_ipXY->at(ij);
            vector<float> track_ipXYSigs = track_ipXYSig->at(ij);
            vector<int> track_sources = track_source->at(ij);
            //vector<float> track_vertex_weights = track_vertex_weight->at(ij);
            vector<float> track_pvWeights = track_pvWeight->at(ij);
            vector<float> track_ref_zs = track_ref_z->at(ij);//analysis_20170523_v0

            if(fabs(jet_eta->at(ij))<jetacut) { // jet eta cut
                if(jet_nef->at(ij)<NemfracCut) {  // neutral fraction
                    if(jet_ntrkpt1[ij]>ntrk1cut) {  // tracks pt>1
                        if(jet_cef->at(ij)<CemfracCut) {  //charged fraction
                            if(jet_fmaxtrkpt[ij]>0.6) continue;//analysis_20170523_v0
//                             if(jet_fnonPU[ij]<0.4) continue;// analysis_20170523_v0
//                             if(jet_fsmallIPsig[ij]<0.3) continue;// analysis_20170523_v0
                            //basicjet[ij]=true;
                            if(ij<4) {
                                basicjet[ij]=true;
                                // # of btagged for leading 4 jets
                                if(jet_csv->at(ij) >= bTagWP) {
                                    nbTaggedTot++;
                                    nbTagged++;
                                    //std::cout << "Jet b-tagged." << std::endl;
                                }
                                else {
                                    nNotbTagTot++;
                                    //std::cout << "Jet not b-tagged." << std::endl;
                                }
                            }
                            //if(jet_alpha[ij]<alphaMaxcut) { // alpha max
                            if(jet_alpha[ij]<alphaMaxcut && jet_alpha[ij]>-1) { // alpha max
                                almostemerging[ij]=true;

                                // uncomment if count only leading jets for event selection later
                                if(ij<4)
                                    nalmostemerging+=1;
                                //if(r0[ij]>maxIPcut) { // max IP cut
                                //if(jet_medip[ij]>maxIPcut) { // med IP cut
                                if(jet_medip[ij]>maxIPcut && jet_medtheta2D[ij]>0) { // med IP cut
                                    emerging[ij]=true;

                                    // uncomment if count only leading jets for event selection later
                                    if(ij<4)
                                        nemerging+=1;
                                    // look at tracks in the emerging jets
//                                     if(otfile){
//                                         hmaxipXYEJ->Fill(r0[ij]);
//                                         hmeanipXYEJ->Fill(jet_meanip[ij]);
//                                         hmedipXYEJ->Fill(jet_medip[ij]);
//                                         hmeanipXYSigEJ->Fill(jet_meanipsig[ij]);
//                                         hmedipXYSigEJ->Fill(jet_medipsig[ij]);
//                                         hlogmeanipXYSigEJ->Fill(jet_logmeanipsig[ij]);
//                                         hlogmedipXYSigEJ->Fill(jet_logmedipsig[ij]);
//                                         hmedtheta2DEJ->Fill(jet_medtheta2D[ij]);
//                                         hlogmedtheta2DEJ->Fill(jet_logmedtheta2D[ij]);
//                                     }
                                    if(jet_meanip[ij]>r0[ij]) std::cout<<"DANGER DANGER"<<std::endl;
//                                     for (unsigned itrack=0; itrack<track_ipXYs.size(); itrack++) {
//                                         if(track_sources[itrack]==0 && ((track_qualitys[itrack] & 4) > 0)) {
//                                             if (fabs(pv_z-track_ref_zs[itrack])>pilecut) continue;//analysis_20170523_v0
//                                             if(otfile) hipXYEJ->Fill(track_ipXYs[itrack]);
//                                             if(otfile) hipXYSigEJ->Fill(track_ipXYSigs[itrack]);
//                                             //if(otfile) htvwEJ->Fill(track_vertex_weights[itrack]);
//                                             if(otfile) htvwEJ->Fill(track_pvWeights[itrack]);
//                                         }
//                                     }
                                }//maxIPcut
                            }//alphaMaxCut
                        }//CemfractCut
                    }//ntrk1cut
                }//NemfracCut
            }//jet eta cut
//             if(!emerging[ij]) {
//                 if(otfile) {
//                     hmaxipXYnEJ->Fill(r0[ij]);
//                     hmeanipXYnEJ->Fill(jet_meanip[ij]);
//                     hmedipXYnEJ->Fill(jet_medip[ij]);
//                     hmeanipXYSignEJ->Fill(jet_meanipsig[ij]);
//                     hmedipXYSignEJ->Fill(jet_medipsig[ij]);
//                     hlogmeanipXYSignEJ->Fill(jet_logmeanipsig[ij]);
//                     hlogmedipXYSignEJ->Fill(jet_logmedipsig[ij]);
//                     hmedtheta2DnEJ->Fill(jet_medtheta2D[ij]);
//                     hlogmedtheta2DnEJ->Fill(jet_logmedtheta2D[ij]);
//                 }
//                 for (unsigned itrack=0; itrack<track_ipXYs.size(); itrack++) {
//                     if(track_sources[itrack]==0 && ((track_qualitys[itrack] & 4) > 0)) {
//                         if (fabs(pv_z-track_ref_zs[itrack])>pilecut) continue;//analysis_20170523_v0
//                         if(otfile) hipXYnEJ->Fill(track_ipXYs[itrack]);
//                         if(otfile) hipXYSignEJ->Fill(track_ipXYSigs[itrack]);
//                     }
//                 }

//             }
            //  std::cout<<"event pt alphaM cef nef ntrkpt1 r0 emerging  almost "<<event<<" "<<jet_pt->at(ij)<<" "<<jet_alpha[ij]<<" "<<jet_cef->at(ij)<<" "<<jet_nef->at(ij)<<" "<<jet_ntrkpt1[ij]<<" "<<r0[ij]<<" "<<emerging[ij]<<" "<<almostemerging[ij]<<std::endl;
        }
        //===================================
        // End of second loop over jets
        //===================================

        if(otfile) hs1D[TString("h_nalemg_"+histcutname)]->Fill(nalmostemerging);
        if(otfile) hs1D[TString("h_nemg_"+histcutname)]->Fill(nemerging);

//         if(otfile) {//basicjets
//             for(int i=0;i<NNNjet;i++) {
//                 if(!basicjet[i]) continue;
//                 if(jet_pt->at(i)>minJetPt && fabs(jet_eta->at(i))<jetacut) {
//                     halpha->Fill(jet_alpha[i]);
//                     if (jet_alpha[i]==0) halphaZero->Fill(jet_alpha[i]);
//                 }
//             }
//         }


        // *************************************************************
        // now start the event selection preparation
        // *************************************************************

        if(otfile) hs1D[TString("count_"+histcutname)]->Fill("filter",1);
        if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(1);

        // PVTrackFraction cut introduced after AN-16-146_v2
        if (PVTrackFraction<=0.1) continue;

        //if (nZeroAlpha>2) continue;//AN-16-146_v2
        // require at least 4 jets
        bool C4jet=true;
        int nemgGoodjet = 0;
        int nalemgGoodjet = 0;
        vector<int> goodjetIdx;
        //for(int nj=0;nj<NNNjet;nj++) {
        for(int nj=0;nj<4;nj++) {
            if (!basicjet[nj] || jet_pt->at(nj)<minJetPt) continue;//comment out this line for leading four jets
            goodjetIdx.push_back(nj);
            //std::cout << "Idx[" << nj <<"] is good jet" << std::endl;
            if (emerging[nj]) nemgGoodjet+=1;
            if (almostemerging[nj]) nalemgGoodjet+=1;
            if(otfile) {
                hs1D[TString("halpha_"+histcutname)]->Fill(jet_alpha[nj]);
                if (jet_alpha[nj]==0) hs1D[TString("halphaZero_"+histcutname)]->Fill(jet_alpha[nj]);
            }
        }

//         bool PVreco = false;
//         if (nZeroAlpha<3) PVreco=true;

        //if(NNNjet<4) continue;
        if(goodjetIdx.size()!=4) {
            //if(goodjetIdx.size()<4) {
            C4jet=false;
            //double HT0 = 0.;
            //for (int i = 0; i<NNNjet; i++) HT0 += jet_pt->at(i);
            //std::cout << "***event " << event << std::endl;
            //std::cout << "NNNjet = " << NNNjet << " =? jet_pt->size() = " << jet_pt->size() << "; HT = " << HT0 << "; jet_pt[" << NNNjet-1 << "] = " << jet_pt->at(NNNjet-1) << std::endl;
            continue;
        }
//         if(goodjetIdx[3]>3) {
//             C4jet=false;
//             continue;
//         }

//         double HT = 0;
//         for(uint i=0;i<4;i++) {// for 4 (basic) leading jets
//             HT += jet_pt->at(goodjetIdx[i]);
//         }
        // HT
        double HT = jet_pt->at(0)+jet_pt->at(1)+jet_pt->at(2)+jet_pt->at(3);

        // HLT efficiency plots:
        bool HLT=false;
        if (hltTrig3n) HLT=true;

        bool CHT=true;
        if(HT<HTcut) CHT=false;
        // jet pt
        bool Cpt1=false;
        bool Cpt2=false;
        bool Cpt3=false;
        bool Cpt4=false;
        if((jet_pt->at(0)>pt1cut)&&(fabs(jet_eta->at(0))<jetacut)) Cpt1=true;
        if((jet_pt->at(1)>pt2cut)&&(fabs(jet_eta->at(1))<jetacut)) Cpt2=true;
        if((jet_pt->at(2)>pt3cut)&&(fabs(jet_eta->at(2))<jetacut)) Cpt3=true;
        if((jet_pt->at(3)>pt4cut)&&(fabs(jet_eta->at(3))<jetacut)) Cpt4=true;
        // basicjet includes eta cut already
//         if(jet_pt->at(goodjetIdx[0])>pt1cut) Cpt1=true;
//         if(jet_pt->at(goodjetIdx[1])>pt2cut) Cpt2=true;
//         if(jet_pt->at(goodjetIdx[2])>pt3cut) Cpt3=true;
//         if(jet_pt->at(goodjetIdx[3])>pt4cut) Cpt4=true;

        // number emerging jets
        bool Cnem = true;
        if(nemerging<NemergingCut) Cnem=false;
        //if(nemerging!=NemergingCut) Cnem=false;
        //if(nemgGoodjet<NemergingCut) Cnem=false;
//         switch (ntagType) {
//         case 0 : {if(nemerging<1) Cnem=false; break;}
//         case 1 : {if(nemerging!=1) Cnem=false; break;}
//         case 12: {if(nemerging!=1) Cnem=false; break;}
//         case 2 : {if(nemerging!=2) Cnem=false; break;}
//         case 21: {if(nemerging!=2) Cnem=false; break;}
//         case 22: {if(nemerging!=2) Cnem=false; break;}
//         case 3 : {if(nemerging!=3) Cnem=false; break;}
//         default: {if(nemerging<NemergingCut) Cnem=false; break;}
//         }
//         switch (ntagType) {
//         case 0 : {if(nemgGoodjet<1) Cnem=false; break;}
//         case 1 : {if(nemgGoodjet!=1) Cnem=false; break;}
//         case 12: {if(nemgGoodjet!=1) Cnem=false; break;}
//         case 2 : {if(nemgGoodjet!=2) Cnem=false; break;}
//         case 21: {if(nemgGoodjet!=2) Cnem=false; break;}
//         case 22: {if(nemgGoodjet!=2) Cnem=false; break;}
//         case 3 : {if(nemgGoodjet!=3) Cnem=false; break;}
//         default: {if(nemgGoodjet<NemergingCut) Cnem=false; break;}
//         }

        bool Canem =true;
        //if(nalmostemerging>=4) Canem=false;
        //if(nalemgGoodjet>=4) Canem=false;

        bool Cmet = false;
        if(met_pt>METcut) Cmet = true;

        //blind
        if(blind) {
            Cnem=false;
            //Canem=false;
        }
        // *************************************************************
        // End of the event selection preparation
        // *************************************************************


        // *************************************************************
        // apply event selection cuts sequentially
        // *************************************************************

        //    std::cout<<"c4jet cht cpt1 cpt2 cpt3 cpt4 cnem "<<C4jet<<" "<<CHT<<" "<<Cpt1<<" "<<Cpt2<<" "<<Cpt3<<" "<<Cpt4<<" "<<Cnem<<std::endl;

        if(C4jet) {
            if(otfile) hs1D[TString("count_"+histcutname)]->Fill("4 jets",1);
            if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(2);

            // calculate HT and require it greater than some cut value
            //if(HLT && CHT && PVreco) {
            if(HLT && CHT) {
                if(otfile) hs1D[TString("count_"+histcutname)]->Fill("HT",1);
                if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(3);

                // do pT cuts on jets  
                if(Cpt1) {
                    if(otfile) hs1D[TString("count_"+histcutname)]->Fill("jet pt1",1);
                    if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(4);

                    if(Cpt2) {
                        if(otfile) hs1D[TString("count_"+histcutname)]->Fill("jet pt2",1);
                        if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(5);

                        if(Cpt3) {
                            if(otfile) hs1D[TString("count_"+histcutname)]->Fill("jet pt3",1);
                            if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(6);

                            if(Cpt4) {
                                int njetsFR = 4;
                                //int njetsFR = goodjetIdx.size();
                                if(otfile) hs1D[TString("count_"+histcutname)]->Fill("jet pt4",1);
                                if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(7);
                                //std::cout << "# btagged = " << nbTagged << " for evt#: " << event << std::endl;

                                // MET cut
                                if (!Cmet) continue;
                                if(otfile) hs1D[TString("count_"+histcutname)]->Fill("MET",1);
                                if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(8);

                                if(otfile) {//Data-driven fake background
                                    //if( ( Canem || (isData && nalmostemerging<4) ) ) {if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemerging);}
                                    if( Canem ) {if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemerging);}
                                    //if( (nemgGoodjet < 1) && ( Canem || (isData && nalemgGoodjet<4) ) ) {
                                    //if( (nemerging < 1) && ( Canem || (isData && nalmostemerging<4) ) ) {
                                    //if( (nemerging < 1) && Canem ) {
                                    if( Canem ) {
                                    //if( Canem || (isData && blind && nalemgGoodjet<4) ) {
                                        //if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemerging);
                                        //if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemgGoodjet);

                                        //std::cout << "# btagged = " << nbTagged << " for evt#: " << event << std::endl;
                                        // ********************************
                                        // Data-driven fake background
                                        // ********************************
                                        int ngoodjet = 0;
                                        for(int nj=0;nj<njetsFR;nj++) {
                                            if (!basicjet[nj] || jet_pt->at(nj)<minJetPt) continue;
                                            ngoodjet+=1;
                                        }
                                        // ngoodjet should be = goodjetIdx.size()
                                        //std::cout << "ngoodjet = " << ngoodjet << " ?= " << "goodjetIdx.size() = " << goodjetIdx.size() << std::endl;
                                        //double nGJwgt = nGJrewgt(goodjetIdx.size());
                                        double nGJwgt = 1.0;

                                        // Default:
                                        double frwgt0 = qcdtools->frWeightFT0(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                              jet_pid_maxEt,njetsFR,minJetPt,varType,cutname);
                                        double frwgt1 = qcdtools->frWeightFT1(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                              jet_pid_maxEt,njetsFR,minJetPt,varType,cutname);
                                        double frwgt2 = qcdtools->frWeightFT2(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                              jet_pid_maxEt,njetsFR,minJetPt,varType,cutname);
                                        double frwgt3 = qcdtools->frWeightFT3(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                              jet_pid_maxEt,njetsFR,minJetPt,varType,cutname);

                                        // Alternatives:
                                        double frwgt12 = qcdtools->frWeightFT12(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                minJetPt,bUnfType,nbTagged,ptType,varType,isData,cutname);
                                        double frwgt21 = qcdtools->frWeightT21(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                               njetsFR,minJetPt,varType);
                                        double frwgt22 = qcdtools->frWeightFT22(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                minJetPt,bUnfType,nbTagged,ptType,varType,isData,cutname);
                                        double frwgtge2 = qcdtools->frWeight1(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                              njetsFR,minJetPt,varType);

                                        // Choose event-wise weight
                                        double frwgt = 1.0;
                                        switch (ntagType) {
                                        case 0 : frwgt = nGJwgt*(1.0 - frwgt0); break;//>=1 tag
                                        case 1 : frwgt = nGJwgt*frwgt1; break;//==1 tag
                                        case 2 : frwgt = nGJwgt*frwgt2; break;//==2 tag
                                        case 3 : frwgt = nGJwgt*frwgt3; break;//==3 tag (not implemented yet)
                                        case 12: frwgt = nGJwgt*frwgt12; break;//==1 tag
                                        case 21: frwgt = nGJwgt*frwgt21; break;//==2 tag; tag-and-probe
                                        case 22: frwgt = nGJwgt*frwgt22; break;//==2 tag; ture b unfolded
                                        default: frwgt = nGJwgt*frwgtge2; break;//>=2tag (need to check calculation)
                                        }

                                        // Fill histograms
                                        hs1D[TString("hnjetFR_"+histcutname)]->Fill(ngoodjet,frwgt);
                                        hs1D[TString("count_"+histcutname)]->Fill("FakeRateBKG",frwgt);
                                        hs1D[TString("acount_"+histcutname)]->Fill(11,frwgt);
                                        hs1D[TString("H_TFR_"+histcutname)]->Fill(HT,frwgt);
                                        np0 += frwgt0;

                                        if (ntagType==1) np1 += frwgt1;
                                        else np1 += frwgt12;

                                        if (ntagType==21) np2 += frwgt21;
                                        else if (ntagType==22 || ntagType==12) np2 += frwgt22;
                                        else np2 += frwgt2;

                                        np3 += frwgt3;

                                        hs1D[TString("h_nemgFR_"+histcutname)]->Fill(0.0,qcdtools->frWeightFT0(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                               jet_pid_maxEt,njetsFR,minJetPt,varType,cutname));

                                        // ntag==1 cases:
                                        if (ntagType==1)
                                            hs1D[TString("h_nemgFR_"+histcutname)]->Fill(1.0,qcdtools->frWeightFT1(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                                   jet_pid_maxEt,njetsFR,minJetPt,varType,cutname));
                                        else
                                            hs1D[TString("h_nemgFR_"+histcutname)]->Fill(1.0,qcdtools->frWeightFT12(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                                    minJetPt,bUnfType,nbTagged,ptType,varType,isData,cutname));

                                        // ntag==2 cases:
                                        if (ntagType==21)
                                            hs1D[TString("h_nemgFR_"+histcutname)]->Fill(2.0,qcdtools->frWeightT21(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                                   njetsFR,minJetPt,varType));
                                        else if (ntagType==22)
                                            hs1D[TString("h_nemgFR_"+histcutname)]->Fill(2.0,qcdtools->frWeightFT22(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                                    minJetPt,bUnfType,nbTagged,ptType,varType,isData,cutname));
                                        else
                                            hs1D[TString("h_nemgFR_"+histcutname)]->Fill(2.0,qcdtools->frWeightFT2(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                                   jet_pid_maxEt,njetsFR,minJetPt,varType,cutname));
                                        // end ntag==2 cases

                                        hs1D[TString("h_nemgFR_"+histcutname)]->Fill(3.0,qcdtools->frWeightFT3(jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                                               jet_pid_maxEt,njetsFR,minJetPt,varType,cutname));
                                        //hs1D[TString("h_nemgFR_"+histcutname)]->Fill(4.0,1-frwgt0-frwgt1-frwgt2);//x-check w.r.t. 3.0 ==> verified
                                        switch (ntagType) {
                                        case 1: {
                                            // For Ntag==1
                                            for(int i3=0;i3<njetsFR;i3++) {
                                                int idx3 = goodjetIdx[i3];
                                                double jfr = qcdtools->fakerateF(jet_pt->at(idx3),jet_eta->at(idx3),jntrack[idx3],
                                                                                 varType,jet_pid_maxEt[idx3],cutname);
                                                for(int i31=0;i31<njetsFR;i31++) {
                                                    int idx31 = goodjetIdx[i31];
                                                    if (i31 != i3) jfr *= (1.0-qcdtools->fakerateF(jet_pt->at(idx31),jet_eta->at(idx31),
                                                                                         jntrack[idx31],varType,jet_pid_maxEt[idx31],cutname));
                                                }
                                                if (!basicjet[idx3] || jet_pt->at(idx3)<minJetPt) continue;
                                                hs1D[TString("hjptaFR_"+histcutname)]->Fill(jet_pt->at(idx3),jfr);
                                                hs1D[TString("hetaaFR_"+histcutname)]->Fill(jet_eta->at(idx3),jfr);
                                                hs1D[TString("hntrkFR_"+histcutname)]->Fill(jntrack[idx3],jfr);
                                                hs1D[TString("hmedipXYSigFR_"+histcutname)]->Fill(jet_medipsig[idx3],jfr);
                                                hs1D[TString("hlogmedipXYSigFR_"+histcutname)]->Fill(jet_logmedipsig[idx3],jfr);
                                                hs1D[TString("hmedtheta2DFR_"+histcutname)]->Fill(jet_medtheta2D[idx3],jfr);
                                                hs1D[TString("hlogmedtheta2DFR_"+histcutname)]->Fill(jet_logmedtheta2D[idx3],jfr);
                                                hs2D[TString("htheta2DvipXYSigFR_"+histcutname)]->Fill(jet_medtheta2D[idx3],jet_medipsig[idx3],jfr);

                                                for(int i4=0;i4<njetsFR;i4++) {
                                                    if (i4 == i3) continue;
                                                    int idx4 = goodjetIdx[i4];
                                                    if (!basicjet[idx4] || jet_pt->at(idx4)<minJetPt) continue;
                                                    double mass = sqrt(
                                                                       pow((jet_e[idx3]+jet_e[idx4]),2) -
                                                                       pow((jet_px[idx3]+jet_px[idx4]),2) -
                                                                       pow((jet_py[idx3]+jet_py[idx4]),2) -
                                                                       pow((jet_pz[idx3]+jet_pz[idx4]),2)
                                                                       );
                                                    hs1D[TString("hmassFR_"+histcutname)]->Fill(mass,jfr);
                                                }
                                            }
                                            break;
                                        }//end case 1 (Ntag==1)
                                        case 12: {
                                            char *hnames[9]={(char*)"hjptaFR",(char*)"hetaaFR",(char*)"hntrkFR",
                                                             (char*)"hmedipXYSigFR",(char*)"hlogmedipXYSigFR",(char*)"hmedtheta2DFR",
                                                             (char*)"hlogmedtheta2DFR",(char*)"hmassFR",(char*)"htheta2DvipXYSigFR"};
                                            // For Ntag==1; flav dep; data unfolded; only for 4 jet events
                                            double frwgttmp=0.;
                                            double frwgts[5];

                                            TMatrixD Mwgt(5,5);
                                            bool goodWgt = qcdtools->UnfoldWgtPtDep(Mwgt, bUnfType, jet_pt, isData);
                                            if (!goodWgt) {
                                                std::cout << "Bad unfolding matrix inversion!" << std::endl;
                                            }

                                            // loop # unfolded "true" b jets
                                            for (int nbT=0;nbT<5;nbT++){
                                                double evtWgt = 0.0;
                                                //if (isData) evtWgt = qcdtools->UnfoldWgtD(bUnfType, nbT, nbTagged, ptType);
                                                //else evtWgt = qcdtools->UnfoldWgt(bUnfType, nbT, nbTagged, ptType);
                                                //evtWgt=qcdtools->UnfoldWgtPtDep(bUnfType, nbT, nbTagged, jet_pt, isData);
                                                evtWgt=Mwgt(nbT, nbTagged);
                                                int ncomb=1;
                                                switch (nbT) {
                                                case 1: ncomb=4; break;
                                                case 2: ncomb=6; break;
                                                case 3: ncomb=4; break;
                                                default: break;
                                                }
                                                evtWgt/=ncomb;
                                                switch (nbT) {
                                                case 0: {
                                                    int flavors[] = {0,0,0,0};
                                                    qcdtools->frWeightUFT1(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,flavors,varType,isData,cutname);
                                                    frwgttmp+=(evtWgt*frwgts[4]);
                                                    qcdtools->fillFRPlots(histcutname,hnames,
                                                                          jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                          jet_medtheta2D, jet_logmedtheta2D,
                                                                          jet_e,jet_px,jet_py,jet_pz,
                                                                          frwgts, evtWgt,1.0);
                                                    break;}
                                                case 1: {
                                                    for (int fidx=0;fidx<4;fidx++) {
                                                        int flavors[] = {0,0,0,0};
                                                        flavors[fidx]=5;
                                                        qcdtools->frWeightUFT1(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,flavors,varType,isData,cutname);
                                                        frwgttmp+=(evtWgt*frwgts[4]);
                                                        qcdtools->fillFRPlots(histcutname,hnames,
                                                                              jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                              jet_medtheta2D, jet_logmedtheta2D,
                                                                              jet_e,jet_px,jet_py,jet_pz,
                                                                              frwgts, evtWgt,1.0);
                                                    }
                                                    break;}
                                                case 2: {
                                                    //flavor double counted due to looping
                                                    for (int fidx0=0;fidx0<4;fidx0++) {
                                                        for (int fidx1=fidx0+1;fidx1<4;fidx1++) {
                                                            int flavors[] = {0,0,0,0};
                                                            flavors[fidx0]=5;
                                                            flavors[fidx1]=5;
                                                            qcdtools->frWeightUFT1(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,flavors,varType,isData,cutname);
                                                            frwgttmp+=(evtWgt*frwgts[4]);
                                                            qcdtools->fillFRPlots(histcutname,hnames,
                                                                                  jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                                  jet_medtheta2D, jet_logmedtheta2D,
                                                                                  jet_e,jet_px,jet_py,jet_pz,
                                                                                  frwgts, evtWgt,1.0);
                                                        }
                                                    }
                                                    break;}
                                                case 3: {
                                                    for (int fidx=0;fidx<4;fidx++) {
                                                        int flavors[] = {5,5,5,5};
                                                        flavors[fidx]=0;
                                                        qcdtools->frWeightUFT1(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,flavors,varType,isData,cutname);
                                                        frwgttmp+=(evtWgt*frwgts[4]);
                                                        qcdtools->fillFRPlots(histcutname,hnames,
                                                                              jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                              jet_medtheta2D, jet_logmedtheta2D,
                                                                              jet_e,jet_px,jet_py,jet_pz,
                                                                              frwgts, evtWgt,1.0);
                                                    }
                                                    break;}
                                                case 4: {
                                                    int flavors[] = {5,5,5,5};
                                                    qcdtools->frWeightUFT1(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,flavors,varType,isData,cutname);
                                                    frwgttmp+=(evtWgt*frwgts[4]);
                                                    qcdtools->fillFRPlots(histcutname,hnames,
                                                                          jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                          jet_medtheta2D, jet_logmedtheta2D,
                                                                          jet_e,jet_px,jet_py,jet_pz,
                                                                          frwgts, evtWgt,1.0);
                                                    break;}
                                                }
                                            }
                                            if (fabs(frwgt-frwgttmp)>1.e-7) {
                                                std::cout << "[WARNING] frwgt = " << std::setprecision(9)
                                                          << frwgt << " != frwgttmp = " << frwgttmp << "!!!!" << std::endl;
                                            }
                                            break;
                                        }//end case 12 (Ntag==1; flav dep; data unfolded)
                                        case 2: {
                                            // For Ntag==2; flav dep
                                            double frwgttmp=0.;
                                            for(int i1=0;i1<njetsFR;i1++) {
                                                int idx1 = goodjetIdx[i1];
                                                if (!basicjet[idx1] || jet_pt->at(idx1)<minJetPt) continue;
                                                double jfr = qcdtools->fakerateF(jet_pt->at(idx1),jet_eta->at(idx1),jntrack[idx1],
                                                                       varType,jet_pid_maxEt[idx1],cutname);
                                                for(int i2=0;i2<njetsFR;i2++) {
                                                    int idx2 = goodjetIdx[i2];
                                                    if (!basicjet[idx2] || jet_pt->at(idx2)<minJetPt) continue;
                                                    if (i2==i1) continue;
                                                    double kfr = jfr*(1.0-qcdtools->fakerateF(jet_pt->at(idx2),jet_eta->at(idx2),jntrack[idx2],
                                                                                    varType,jet_pid_maxEt[idx2],cutname));
                                                    for(int i3=0;i3<njetsFR;i3++) {
                                                        if (i3==i1 || i3==i2) continue;
                                                        int idx3 = goodjetIdx[i3];
                                                        if (!basicjet[idx3] || jet_pt->at(idx3)<minJetPt) continue;
                                                        double lfr = kfr*qcdtools->fakerateF(jet_pt->at(idx3),jet_eta->at(idx3),jntrack[idx3],
                                                                                   varType,jet_pid_maxEt[idx3],cutname);
                                                        for(int i4=0;i4<njetsFR;i4++) {
                                                            if (i4==i1 || i4==i2 || i4==i3) continue;
                                                            int idx4 = goodjetIdx[i4];
                                                            if (!basicjet[idx4] || jet_pt->at(idx4)<minJetPt) continue;
                                                            lfr *= (1.0-qcdtools->fakerateF(jet_pt->at(idx4),jet_eta->at(idx4),jntrack[idx4],
                                                                                  varType,jet_pid_maxEt[idx4],cutname));
                                                        }
                                                        //take care of combinatorics: e.g., var(1) double counted by (1_tag,2)(3_tag,4) and (1_tag,4)(3_tag,2)
                                                        double varWgt = lfr/(ngoodjet-2);
                                                        frwgttmp += varWgt;

                                                        hs1D[TString("hjptaFR_"+histcutname)]->Fill(jet_pt->at(idx1),varWgt);
                                                        hs1D[TString("hetaaFR_"+histcutname)]->Fill(jet_eta->at(idx1),varWgt);
                                                        hs1D[TString("hntrkFR_"+histcutname)]->Fill(jntrack[idx1],varWgt);
                                                        hs1D[TString("hmedipXYSigFR_"+histcutname)]->Fill(jet_medipsig[idx1],varWgt);
                                                        hs1D[TString("hlogmedipXYSigFR_"+histcutname)]->Fill(jet_logmedipsig[idx1],varWgt);
                                                        hs1D[TString("hmedtheta2DFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],varWgt);
                                                        hs1D[TString("hlogmedtheta2DFR_"+histcutname)]->Fill(jet_logmedtheta2D[idx1],varWgt);
                                                        hs2D[TString("htheta2DvipXYSigFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],jet_medipsig[idx1],varWgt);

                                                        double mass = sqrt(
                                                                           pow((jet_e[idx1]+jet_e[idx2]),2) -
                                                                           pow((jet_px[idx1]+jet_px[idx2]),2) -
                                                                           pow((jet_py[idx1]+jet_py[idx2]),2) -
                                                                           pow((jet_pz[idx1]+jet_pz[idx2]),2)
                                                                           );
                                                        hs1D[TString("hmassFR_"+histcutname)]->Fill(mass,lfr);
                                                    }
                                                }
                                            }
                                            if (fabs(frwgt*2-frwgttmp)>1.e-7) {
                                                std::cout << "[WARNING] frwgtx2 = " << std::setprecision(9)
                                                          << frwgt*2 << " != frwgttmp = " << frwgttmp << "!!!!" << std::endl;
                                            }
                                            break;
                                        }//end case 2 (Ntag==2; flav dep)
                                        case 21: {
                                            // For Ntag==2; tag-and-probe
                                            double frwgttmp=0.;
                                            for(int i1=0;i1<njetsFR;i1++) {
                                                int idx1 = goodjetIdx[i1];
                                                if (!basicjet[idx1] || jet_pt->at(idx1)<minJetPt) continue;
                                                double jfr = qcdtools->fakerate(jet_pt->at(idx1),jet_eta->at(idx1),jntrack[idx1],varType);
                                                for(int i2=0;i2<njetsFR;i2++) {
                                                    int idx2 = goodjetIdx[i2];
                                                    if (!basicjet[idx2] || jet_pt->at(idx2)<minJetPt) continue;
                                                    if (i2==i1) continue;
                                                    double kfr = jfr*(1.0-qcdtools->fakerateTP(jet_pt->at(idx2),jet_eta->at(idx2),jntrack[idx2],varType));
                                                    //kfr *= (1.0-qcdtools->fakerate(jet_pt->at(idx2),jet_eta->at(idx2),jntrack[idx2],varType));
                                                    for(int i3=0;i3<njetsFR;i3++) {
                                                        if (i3==i1 || i3==i2) continue;
                                                        int idx3 = goodjetIdx[i3];
                                                        if (!basicjet[idx3] || jet_pt->at(idx3)<minJetPt) continue;
                                                        double lfr = kfr*qcdtools->fakerateTP(jet_pt->at(idx3),jet_eta->at(idx3),jntrack[idx3],varType);
                                                        for(int i4=0;i4<njetsFR;i4++) {
                                                            if (i4==i1 || i4==i2 || i4==i3) continue;
                                                            int idx4 = goodjetIdx[i4];
                                                            if (!basicjet[idx4] || jet_pt->at(idx4)<minJetPt) continue;
                                                            lfr *= (1.0-qcdtools->fakerateTP(jet_pt->at(idx4),jet_eta->at(idx4),jntrack[idx4],varType));
                                                            //lfr *= (1.0-qcdtools->fakerate(jet_pt->at(idx4),jet_eta->at(idx4),jntrack[idx4],varType));
                                                        }
                                                        //take care of combinatorics: e.g., var(1) double counted by (1_tag,2)(3_tag,4) and (1_tag,4)(3_tag,2)
                                                        double varWgt = lfr/(ngoodjet-2);
                                                        frwgttmp += varWgt;

                                                        hs1D[TString("hjptaFR_"+histcutname)]->Fill(jet_pt->at(idx1),varWgt);
                                                        hs1D[TString("hetaaFR_"+histcutname)]->Fill(jet_eta->at(idx1),varWgt);
                                                        hs1D[TString("hntrkFR_"+histcutname)]->Fill(jntrack[idx1],varWgt);
                                                        hs1D[TString("hmedipXYSigFR_"+histcutname)]->Fill(jet_medipsig[idx1],varWgt);
                                                        hs1D[TString("hlogmedipXYSigFR_"+histcutname)]->Fill(jet_logmedipsig[idx1],varWgt);
                                                        hs1D[TString("hmedtheta2DFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],varWgt);
                                                        hs1D[TString("hlogmedtheta2DFR_"+histcutname)]->Fill(jet_logmedtheta2D[idx1],varWgt);
                                                        hs2D[TString("htheta2DvipXYSigFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],jet_medipsig[idx1],varWgt);

                                                        double mass = sqrt(
                                                                           pow((jet_e[idx1]+jet_e[idx2]),2) -
                                                                           pow((jet_px[idx1]+jet_px[idx2]),2) -
                                                                           pow((jet_py[idx1]+jet_py[idx2]),2) -
                                                                           pow((jet_pz[idx1]+jet_pz[idx2]),2)
                                                                           );
                                                        hs1D[TString("hmassFR_"+histcutname)]->Fill(mass,lfr);
                                                    }
                                                }
                                            }
                                            if (fabs(frwgt*2-frwgttmp)>1.e-7) {
                                                std::cout << "[WARNING] frwgtx2 = " << std::setprecision(9)
                                                          << frwgt*2 << " != frwgttmp = " << frwgttmp << "!!!!" << std::endl;
                                            }
                                            break;
                                        }//end case 21 (Ntag==2; tag-and-probe)
                                        case 22: {
                                            char *hnames[9]={(char*)"hjptaFR",(char*)"hetaaFR",(char*)"hntrkFR",
                                                             (char*)"hmedipXYSigFR",(char*)"hlogmedipXYSigFR",(char*)"hmedtheta2DFR",
                                                             (char*)"hlogmedtheta2DFR",(char*)"hmassFR",(char*)"htheta2DvipXYSigFR"};
                                            // For Ntag==1; flav dep; data unfolded; only for 4 jet events
                                            double frwgttmp=0.;
                                            double frwgts[7];

                                            TMatrixD Mwgt(5,5);
                                            bool goodWgt = qcdtools->UnfoldWgtPtDep(Mwgt, bUnfType, jet_pt, isData);

//                                             for (int i=0;i<5;i++){
//                                                 for (int j=0;j<5;j++){
//                                                     std::cout << "Mwgt(" << i << "," << j <<")=" << Mwgt(i,j) << std::endl;
//                                                 }
//                                             }

                                            if (!goodWgt) {
                                                std::cout << "Bad unfolding matrix inversion!" << std::endl;
                                            }

                                            // loop # unfolded "true" b jets
                                            for (int nbT=0;nbT<5;nbT++){
                                                double evtWgt = 0.0;
                                                //if (isData) evtWgt = qcdtools->UnfoldWgtD(bUnfType, nbT, nbTagged, ptType);
                                                //else evtWgt = qcdtools->UnfoldWgt(bUnfType, nbT, nbTagged, ptType);
                                                //evtWgt=qcdtools->UnfoldWgtPtDep(bUnfType, nbT, nbTagged, jet_pt, isData);
                                                evtWgt=Mwgt(nbT, nbTagged);
                                                double ncomb=1;
                                                switch (nbT) {
                                                case 1: ncomb=4; break;
                                                case 2: ncomb=6; break;
                                                case 3: ncomb=4; break;
                                                default: break;
                                                }
                                                evtWgt/=ncomb;
                                                switch (nbT) {
                                                case 0: {
                                                    int flavors[] = {0,0,0,0};
                                                    qcdtools->frWeightUFT23(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                            flavors,varType,isData,cutname);
                                                    frwgttmp+=(evtWgt*frwgts[6]);
                                                    qcdtools->fillFRPlots22(histcutname,hnames,
                                                                            jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                            jet_medtheta2D, jet_logmedtheta2D,
                                                                            jet_e,jet_px,jet_py,jet_pz,
                                                                            frwgts, evtWgt, 1.0);
                                                    break;}
                                                case 1: {
                                                    for (int fidx=0;fidx<4;fidx++) {
                                                        int flavors[] = {0,0,0,0};
                                                        flavors[fidx]=5;
                                                        qcdtools->frWeightUFT23(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                flavors,varType,isData,cutname);
                                                        frwgttmp+=(evtWgt*frwgts[6]);
                                                        qcdtools->fillFRPlots22(histcutname,hnames,
                                                                                jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                                jet_medtheta2D, jet_logmedtheta2D,
                                                                                jet_e,jet_px,jet_py,jet_pz,
                                                                                frwgts, evtWgt, 1.0);
                                                    }
                                                    break;}
                                                case 2: {
                                                    //flavor double counted due to looping
                                                    for (int fidx0=0;fidx0<4;fidx0++) {
                                                        for (int fidx1=fidx0+1;fidx1<4;fidx1++) {
                                                            int flavors[] = {0,0,0,0};
                                                            flavors[fidx0]=5;
                                                            flavors[fidx1]=5;
                                                            qcdtools->frWeightUFT23(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                    flavors,varType,isData,cutname);
                                                            frwgttmp+=(evtWgt*frwgts[6]);
                                                            qcdtools->fillFRPlots22(histcutname,hnames,
                                                                                    jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                                    jet_medtheta2D, jet_logmedtheta2D,
                                                                                    jet_e,jet_px,jet_py,jet_pz,
                                                                                    frwgts, evtWgt, 1.0);
                                                        }
                                                    }
                                                    break;}
                                                case 3: {
                                                    for (int fidx=0;fidx<4;fidx++) {
                                                        int flavors[] = {5,5,5,5};
                                                        flavors[fidx]=0;
                                                        qcdtools->frWeightUFT23(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                                flavors,varType,isData,cutname);
                                                        frwgttmp+=(evtWgt*frwgts[6]);
                                                        qcdtools->fillFRPlots22(histcutname,hnames,
                                                                                jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                                jet_medtheta2D, jet_logmedtheta2D,
                                                                                jet_e,jet_px,jet_py,jet_pz,
                                                                                frwgts, evtWgt, 1.0);
                                                    }
                                                    break;}
                                                case 4: {
                                                    int flavors[] = {5,5,5,5};
                                                    qcdtools->frWeightUFT23(frwgts,jet_pt,jet_eta,goodjetIdx,jntrack,
                                                                            flavors,varType,isData,cutname);
                                                    frwgttmp+=(evtWgt*frwgts[6]);
                                                    qcdtools->fillFRPlots22(histcutname,hnames,
                                                                            jet_pt, jet_eta, goodjetIdx,jntrack,jet_medipsig,jet_logmedipsig,
                                                                            jet_medtheta2D, jet_logmedtheta2D,
                                                                            jet_e,jet_px,jet_py,jet_pz,
                                                                            frwgts, evtWgt, 1.0);
                                                    break;}
                                                }
                                            }
                                            if (fabs(frwgt-frwgttmp)>1.e-7) {
                                                std::cout << "[WARNING] frwgt = " << std::setprecision(9)
                                                          << frwgt << " != frwgttmp = " << frwgttmp << "!!!!" << std::endl;
                                            }
                                            break;
                                        }//end case 22 (Ntag==2; flav dep; data unfolded)
                                        case 3: {
                                            // For Ntag==3; flav dep
                                            double frwgttmp=0.;
                                            for(int i1=0;i1<njetsFR;i1++) {
                                                int idx1 = goodjetIdx[i1];
                                                if (!basicjet[idx1] || jet_pt->at(idx1)<minJetPt) continue;
                                                double jfr = qcdtools->fakerateF(jet_pt->at(idx1),jet_eta->at(idx1),jntrack[idx1],
                                                                       varType,jet_pid_maxEt[idx1],cutname);
                                                for(int i2=0;i2<njetsFR;i2++) {
                                                    int idx2 = goodjetIdx[i2];
                                                    if (!basicjet[idx2] || jet_pt->at(idx2)<minJetPt) continue;
                                                    if (i2==i1) continue;
                                                    double kfr = jfr*qcdtools->fakerateF(jet_pt->at(idx2),jet_eta->at(idx2),jntrack[idx2],
                                                                               varType,jet_pid_maxEt[idx2],cutname);
                                                    for(int i3=0;i3<njetsFR;i3++) {
                                                        if (i3==i1 || i3==i2) continue;
                                                        int idx3 = goodjetIdx[i3];
                                                        if (!basicjet[idx3] || jet_pt->at(idx3)<minJetPt) continue;
                                                        double lfr = kfr*qcdtools->fakerateF(jet_pt->at(idx3),jet_eta->at(idx3),jntrack[idx3],
                                                                                   varType,jet_pid_maxEt[idx3],cutname);
                                                        for(int i4=0;i4<njetsFR;i4++) {
                                                            if (i4==i1 || i4==i2 || i4==i3) continue;
                                                            int idx4 = goodjetIdx[i4];
                                                            if (!basicjet[idx4] || jet_pt->at(idx4)<minJetPt) continue;
                                                            double mfr = lfr*(1.0-qcdtools->fakerateF(jet_pt->at(idx4),jet_eta->at(idx4),jntrack[idx4],
                                                                                            varType,jet_pid_maxEt[idx4],cutname));
                                                            if (njetsFR>4) {
                                                                for(int i5=0;i5<njetsFR;i5++) {
                                                                    if (i5==i1 || i5==i2 || i5==i3 || i5==i4) continue;
                                                                    int idx5 = goodjetIdx[i5];
                                                                    if (!basicjet[idx5] || jet_pt->at(idx5)<minJetPt) continue;
                                                                    mfr *= (1.0-qcdtools->fakerateF(jet_pt->at(idx5),jet_eta->at(idx5),jntrack[idx5],
                                                                                          varType,jet_pid_maxEt[idx5],cutname));
                                                                }
                                                            }
                                                            double varWgt = mfr/2;
                                                            frwgttmp += varWgt;

                                                            hs1D[TString("hjptaFR_"+histcutname)]->Fill(jet_pt->at(idx1),varWgt);
                                                            hs1D[TString("hetaaFR_"+histcutname)]->Fill(jet_eta->at(idx1),varWgt);
                                                            hs1D[TString("hntrkFR_"+histcutname)]->Fill(jntrack[idx1],varWgt);
                                                            hs1D[TString("hmedipXYSigFR_"+histcutname)]->Fill(jet_medipsig[idx1],varWgt);
                                                            hs1D[TString("hlogmedipXYSigFR_"+histcutname)]->Fill(jet_logmedipsig[idx1],varWgt);
                                                            hs1D[TString("hmedtheta2DFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],varWgt);
                                                            hs1D[TString("hlogmedtheta2DFR_"+histcutname)]->Fill(jet_logmedtheta2D[idx1],varWgt);
                                                            hs2D[TString("htheta2DvipXYSigFR_"+histcutname)]->Fill(jet_medtheta2D[idx1],jet_medipsig[idx1],varWgt);

                                                            double mass = sqrt(
                                                                               pow((jet_e[idx1]+jet_e[idx4]),2) -
                                                                               pow((jet_px[idx1]+jet_px[idx4]),2) -
                                                                               pow((jet_py[idx1]+jet_py[idx4]),2) -
                                                                               pow((jet_pz[idx1]+jet_pz[idx4]),2)
                                                                               );
                                                            //hs1D[TString("hmassFR_"+histcutname)]->Fill(mass,mfr);
                                                            hs1D[TString("hmassFR_"+histcutname)]->Fill(mass,varWgt);
                                                        }
                                                    }
                                                }
                                            }
                                            if (fabs(frwgt*3-frwgttmp)>1.e-7) {
                                                std::cout << "[WARNING] frwgt = " << std::setprecision(9)
                                                          << frwgt*3 << " != frwgttmp = " << frwgttmp << "!!!!" << std::endl;
                                            }
                                            break;
                                        }//end case 3 (Ntag==3; flav dep)

                                        default: break;
                                        }//end switch
                                    }// Canem Data-driven fake background
                                }
                                // ********************************
                                // end otfile Data-driven fake background
                                // ********************************


                                // ********************************
                                // SR plots
                                // ********************************
                                // require at least N emerging jets
                                if(Cnem) {//NOTE: revert Canem and Cnem for the time being
                                    if(otfile) hs1D[TString("count_"+histcutname)]->Fill("emerging",1);
                                    if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(9);

                                    if(Canem) {
//                                         if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemerging);
//                                         //if(otfile) hs1D[TString("h_nemgSR_"+histcutname)]->Fill(nemgGoodjet);

                                        if(otfile) hs1D[TString("count_"+histcutname)]->Fill("almostemerging",1);
                                        if(otfile) hs1D[TString("acount_"+histcutname)]->Fill(10);

                                        npass+=1;
                                        std::cout<<"passing run lumi event filename is "<<run<<" "<<lumi<<" "<<event<<" "<<inputfilename<<std::endl;
                                        for(int i=0;i<4;i++) {
                                            int idx = goodjetIdx[i];
                                            std::cout<<"  for jet "<<i<<" pt eta nef cfe ntrkpt1 alphamax r0 medIP"<<std::endl;
                                            std::cout<<"     "<<jet_pt->at(idx)<<" "<<jet_eta->at(idx)<<" "<< jet_nef->at(idx)<<" "<< jet_cef->at(idx)<<" "<<jet_ntrkpt1[idx]<<" "<<jet_alpha[idx]<<" "<<r0[idx]<<" " << jet_medip[idx] << " " <<std::endl;
                                        }
                                        if(otfile) {
                                            hs1D[TString("H_TSR_"+histcutname)]->Fill(HT);

                                            int ngoodjetSR = 0;
                                            for(int nj=0;nj<4;nj++) {
                                                if (!basicjet[nj] || jet_pt->at(nj)<minJetPt) continue;
                                                ngoodjetSR+=1;
                                            }
                                            hs1D[TString("hnjetSR_"+histcutname)]->Fill(ngoodjetSR);

                                            for(int i5=0;i5<4;i5++) {
                                                int idx5 = goodjetIdx[i5];
                                                if (jet_pt->at(idx5)<minJetPt) continue;
                                                if (emerging[idx5] && otfile){
                                                    hs1D[TString("hjptaSR_"+histcutname)]->Fill(jet_pt->at(idx5));
                                                    hs1D[TString("hetaaSR_"+histcutname)]->Fill(jet_eta->at(idx5));
                                                    hs1D[TString("hntrkSR_"+histcutname)]->Fill(jntrack[idx5]);
                                                    hs1D[TString("hmedipXYSigSR_"+histcutname)]->Fill(jet_medipsig[idx5]);
                                                    hs1D[TString("hlogmedipXYSigSR_"+histcutname)]->Fill(jet_logmedipsig[idx5]);
                                                    hs1D[TString("hmedtheta2DSR_"+histcutname)]->Fill(jet_medtheta2D[idx5]);
                                                    hs1D[TString("hlogmedtheta2DSR_"+histcutname)]->Fill(jet_logmedtheta2D[idx5]);
                                                    hs2D[TString("htheta2DvipXYSigSR_"+histcutname)]->Fill(jet_medtheta2D[idx5],jet_medipsig[idx5]);
                                                }
                                                for(int i6=i5+1;i6<4;i6++) {
                                                    int idx6 = goodjetIdx[i6];
                                                    if (jet_pt->at(idx6)<minJetPt) continue;
                                                    if ((emerging[idx5]&&!emerging[idx6])||(!emerging[idx5]&&emerging[idx6])) {
                                                        double mass = sqrt(
                                                                           pow((jet_e[idx5]+jet_e[idx6]),2) -
                                                                           pow((jet_px[idx5]+jet_px[idx6]),2) -
                                                                           pow((jet_py[idx5]+jet_py[idx6]),2) -
                                                                           pow((jet_pz[idx5]+jet_pz[idx6]),2)
                                                                           );
                                                        hs1D[TString("hmassSR_"+histcutname)]->Fill(mass);
                                                    }
                                                }
                                            }
                                        }
                                        std::cout<<"npass  event is "<<npass<<" "<<event<<std::endl;
                                        std::cout<<"nemerging nalmostemerging "<<nemerging<<" "<<nalmostemerging<<std::endl;
                                        std::cout<<"nemgGoodjet nalemgGoodjet "<<nemgGoodjet<<" "<<nalemgGoodjet<<std::endl;
                                    }//Canem (SR)
                                }//Cnem (SR)
                                // ********************************
                                // End of SR plots
                                // ********************************


                            }//Cpt4
                        }//Cpt3
                    }//Cpt2
                }//Cpt1
            }//HLT&&CHT
        }//C4jet
        // *************************************************************
        // End of apply event selection cuts sequentially
        // *************************************************************

    }// end of loop over events
    // *************************************************************
    // End of loop over events
    // *************************************************************

    std::cout << "np0 = " << np0 << std::endl;
    std::cout << "np1 = " << np1 << std::endl;
    std::cout << "np2 = " << np2 << std::endl;
    std::cout << "np3 = " << np3 << std::endl;

    std::cout << "Total # btagged = " << nbTaggedTot << std::endl;
    std::cout << "Total # not btagged = " << nNotbTagTot << std::endl;
    std::cout << "Total # btagged (ALL) = " << nbTaggedTotAll << std::endl;
    std::cout << "Total # not btagged (ALL) = " << nNotbTagTotAll << std::endl;

    if(otfile) {
        TFile myfile(outputfilename,"RECREATE");
        hs1D[TString("count_"+histcutname)]->LabelsDeflate();
        hs1D[TString("count_"+histcutname)]->LabelsOption("v");
        //  hs1D[TString("count_"+histcutname)]->LabelsOption("a");
        hs1D[TString("eventCountPreTrigger_"+histcutname)]->Write();
        hs1D[TString("acount_"+histcutname)]->Write();
        hs1D[TString("count_"+histcutname)]->Write();
        hs1D[TString("hnjet_"+histcutname)]->Write();
        hs1D[TString("hnjetSR_"+histcutname)]->Write();
        hs1D[TString("hnjetFR_"+histcutname)]->Write();
        hs1D[TString("hetaaSR_"+histcutname)]->Write();
        hs1D[TString("hetaaFR_"+histcutname)]->Write();
        hs1D[TString("halpha_"+histcutname)]->Write();
        hs1D[TString("halphaZero_"+histcutname)]->Write();
        hs1D[TString("H_TSR_"+histcutname)]->Write();
        hs1D[TString("H_TFR_"+histcutname)]->Write();
        hs1D[TString("h_nemg_"+histcutname)]->Write();
        hs1D[TString("h_nemgSR_"+histcutname)]->Write();
        hs1D[TString("h_nemgFR_"+histcutname)]->Write();
        hs1D[TString("h_nalemg_"+histcutname)]->Write();
        hs1D[TString("hntrkSR_"+histcutname)]->Write();
        hs1D[TString("hntrkFR_"+histcutname)]->Write();
        hs1D[TString("hjptaSR_"+histcutname)]->Write();
        hs1D[TString("hjptaFR_"+histcutname)]->Write();
        hs1D[TString("hmassSR_"+histcutname)]->Write();
        hs1D[TString("hmassFR_"+histcutname)]->Write();
        hs1D[TString("hmedipXYSigSR_"+histcutname)]->Write();
        hs1D[TString("hmedtheta2DSR_"+histcutname)]->Write();
        hs1D[TString("hlogmedipXYSigSR_"+histcutname)]->Write();
        hs1D[TString("hlogmedtheta2DSR_"+histcutname)]->Write();
        hs1D[TString("hmedipXYSigFR_"+histcutname)]->Write();
        hs1D[TString("hmedtheta2DFR_"+histcutname)]->Write();
        hs1D[TString("hlogmedipXYSigFR_"+histcutname)]->Write();
        hs1D[TString("hlogmedtheta2DFR_"+histcutname)]->Write();

        //2d
        hs2D[TString("htheta2DvipXYSigSR_"+histcutname)]->Write();
        hs2D[TString("htheta2DvipXYSigFR_"+histcutname)]->Write();

        myfile.Close();
    }

    tt->ResetBranchAddresses();
  
    delete jet_index;
    delete jet_source;
    delete jet_pt;
    delete jet_eta;
    delete jet_phi;
    delete jet_alphaMax;
    delete jet_cef;
    delete jet_nef;
    delete jet_chf;
    delete jet_theta2D;
    //  delete jet_phf;
    delete jet_csv;
    delete track_pt;
    delete track_eta;
    delete track_source;
    delete track_index;
    delete track_jet_index;
    delete track_vertex_index;
    delete track_algo;
    delete track_quality;
    delete track_vertex_weight;
    delete track_pvWeight;
    delete track_ipZ;
    delete track_ipXY;
    delete track_ipXYSig;
    delete qcdtools;


    f->Close();
  


    return npass;
}
