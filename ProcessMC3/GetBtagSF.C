// without CMSSW:
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"


void GetBtagSF(){
    // setup calibration readers
    // Before CMSSW_8_0_12 and CMSSW_8_1_0, one reader was needed 
    // for every OperatingPoint/MeasurementType/SysType combination.
    BTagCalibrationStandalone calib("CSVv2", "BtagFiles/CSVv2_Moriond17_mistag_G_H.csv");
    BTagCalibrationStandaloneReader readerB(BTagEntryStandalone::OP_LOOSE, "central");  // nominal
    BTagCalibrationStandaloneReader readerB_up(BTagEntryStandalone::OP_LOOSE, "up");    // sys up
    BTagCalibrationStandaloneReader readerB_do(BTagEntryStandalone::OP_LOOSE, "down");  // sys down

    BTagCalibrationStandaloneReader readerL(BTagEntryStandalone::OP_LOOSE, "central");  // nominal
    BTagCalibrationStandaloneReader readerL_up(BTagEntryStandalone::OP_LOOSE, "up");    // sys up
    BTagCalibrationStandaloneReader readerL_do(BTagEntryStandalone::OP_LOOSE, "down");  // sys down

    readerB.load(calib, BTagEntryStandalone::FLAV_B,"comb");
    readerB_up.load(calib, BTagEntryStandalone::FLAV_B,"comb");
    readerB_do.load(calib, BTagEntryStandalone::FLAV_B,"comb");
    readerL.load(calib, BTagEntryStandalone::FLAV_UDSG,"incl");
    readerL_up.load(calib, BTagEntryStandalone::FLAV_UDSG,"incl");
    readerL_do.load(calib, BTagEntryStandalone::FLAV_UDSG,"incl");

//     double jet_scalefactor0 =  readerB.eval_auto_bounds("central", BTagEntryStandalone::FLAV_B, 1.2, 120.0, 0.5426);
//     std::cout<< "jet_scalefactor B = " << jet_scalefactor0 << std::endl;

    double jet_scalefactor =  readerB.eval(BTagEntryStandalone::FLAV_B, 1.0, 1000.0);
    std::cout<< "jet_scalefactor B = " << jet_scalefactor << std::endl;

    double jet_scalefactor_up =  readerB_up.eval(BTagEntryStandalone::FLAV_B, 1.0, 1000.0);
    std::cout<< "jet_scalefactor_up B = " << jet_scalefactor_up << std::endl;

    double jet_scalefactor_do =  readerB_do.eval(BTagEntryStandalone::FLAV_B, 1.0, 1000.0);
    std::cout<< "jet_scalefactor_do B = " << jet_scalefactor_do << std::endl;

    double jet_scalefactorL =  readerL.eval(BTagEntryStandalone::FLAV_UDSG, 1.0, 1000.0);
    std::cout<< "jet_scalefactor L = " << jet_scalefactorL << std::endl;

    double jet_scalefactorL_up =  readerL_up.eval(BTagEntryStandalone::FLAV_UDSG, 1.0, 1000.0);
    std::cout<< "jet_scalefactor_up L = " << jet_scalefactorL_up << std::endl;

    double jet_scalefactorL_do =  readerL_do.eval(BTagEntryStandalone::FLAV_UDSG, 1.0, 1000.0);
    std::cout<< "jet_scalefactor_do L = " << jet_scalefactorL_do << std::endl;

}

// in your event loop
// float MaxBJetPt = 669.9, MaxLJetPt = 999.9;  // value must be below the boundary
// for (auto event : events) {
//     for (auto b_jet : event.b_jets) {
//         float JetPt = b_jet.pt(); bool DoubleUncertainty = false;
//         if (JetPt>MaxBJetPt)  { // use MaxLJetPt for  light jets
//             JetPt = MaxBJetPt; 
//             DoubleUncertainty = true;
//         }

//         // Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
//         double jet_scalefactor = readerB.eval(BTagEntryStandalone::FLAV_B, b_jet.eta(), JetPt); 
//         double jet_scalefactor_up =  readerB_up.eval(BTagEntryStandalone::FLAV_B, b_jet.eta(), JetPt); 
//         double jet_scalefactor_do =  readerB_do.eval(BTagEntryStandalone::FLAV_B, b_jet.eta(), JetPt); 

//         if (DoubleUncertainty) {
//             jet_scalefactor_up = 2*(jet_scalefactor_up - jet_scalefactor) + jet_scalefactor; 
//             jet_scalefactor_do = 2*(jet_scalefactor_do - jet_scalefactor) + jet_scalefactor; 
//         }
//     }
//  }
