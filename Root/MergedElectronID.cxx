#include "HGamGamStar/MergedElectronID.h"

//______________________________________________________________________________
HG::MergedElectronID::MergedElectronID()
{

}

//______________________________________________________________________________
HG::MergedElectronID::~MergedElectronID()
{

}

//______________________________________________________________________________

EL::StatusCode HG::MergedElectronID::initialize(Config &config)
{
  m_PreselNPassBlayer = config.getInt("MergedElectrons.Preselection.NtracksPassingBlayer",1);
  m_PreselRhad        = config.getNum("MergedElectrons.Preselection.RhadMin",0.10);
  m_mergedElePtCut    = config.getNum("MergedElectrons.Selection.PtPreCutGeV",20.) * GeV;
  m_mergedEleEtaCut   = config.getNum("MergedElectrons.Selection.MaxAbsEta"  ,2.37);

  return EL::StatusCode::SUCCESS;
}

//______________________________________________________________________________
void HG::MergedElectronID::decorateMergedVariables(xAOD::Electron &ele,xAOD::TrackParticle &trk1,xAOD::TrackParticle &trk2){
  // Note: YOU MUST INITIALIZE ANY ACCESSOR in
  // HiggsGamGamStarCutflowAndMxAOD::AddElectronDecorations(); otherwise you may cause
  // crashes in the code when files are merged.
  // (E.g. every *electron* should have a decorator, even if it's not a merged electron!)

  // Get the index!
  unsigned int index1 = HG::MapHelpers::getMatchingTrackIndex(&ele,&trk1);
  unsigned int index2 = HG::MapHelpers::getMatchingTrackIndex(&ele,&trk2);

  TrkAcc::mergedTrackParticleIndex(trk1) = index1;
  TrkAcc::mergedTrackParticleIndex(trk2) = index2;

  EleAcc::delta_z0_tracks(ele) = TrkAcc::z0pv(trk1) - TrkAcc::z0pv(trk2);
  EleAcc::delta_z0sinTheta_tracks(ele) = TrkAcc::z0sinTheta(trk1) - TrkAcc::z0sinTheta(trk2);

  double pTrk1 = trk1.p4().Rho(); // same as p4().P()
  double pTrk2 = trk2.p4().Rho(); // same as p4().P()
  EleAcc::EOverP0P1(ele) = ele.e() / (pTrk1 + pTrk2);

  // Local calculation of the extrapolated track position (from Perigee)
  AngularPosition trk1AngPos = getExtrapolatedTrackPosition(&trk1, extrapolationStartPositionEnum::Perigee, false, false, false);
  AngularPosition trk2AngPos = getExtrapolatedTrackPosition(&trk2, extrapolationStartPositionEnum::Perigee, false, false, false);
  EleAcc::dRExtrapTrk12(ele) = trk1AngPos.deltaR(trk2AngPos);

  // Local calculation of the extrapolated track position (from LastMeasurement)
  AngularPosition trk1AngPos_LM = getExtrapolatedTrackPosition(&trk1, extrapolationStartPositionEnum::LastMeasurement, false, false, false);
  AngularPosition trk2AngPos_LM = getExtrapolatedTrackPosition(&trk2, extrapolationStartPositionEnum::LastMeasurement, false, false, false);
  EleAcc::dRExtrapTrk12_LM(ele) = trk1AngPos_LM.deltaR(trk2AngPos_LM);

  // Note - these tracks are already pt-sorted!
  TLorentzVector tlv1 = trk1.p4();
  TLorentzVector tlv2 = trk2.p4();
  tlv1.SetPtEtaPhiM( tlv1.Pt(), tlv1.Eta(), tlv1.Phi(), 0.510998 ); // ele.m == 0.510998
  tlv2.SetPtEtaPhiM( tlv2.Pt(), tlv2.Eta(), tlv2.Phi(), 0.510998 );

  // New DeltaPhi_trktrk_IP decoration
  float dPhiIP = trk1.charge() * tlv1.DeltaPhi( tlv2 );
  float cutValue = 0.010*TMath::Exp(-0.374538*(tlv2.Pt()/1000.-0.5)) + 0.001;
  EleAcc::deltaPhiTrksIP(ele) = dPhiIP;
  EleAcc::passDeltaPhiIPCut(ele) = (-0.02 < dPhiIP && dPhiIP < cutValue);

  return;
}

//______________________________________________________________________________
bool HG::MergedElectronID::passPIDCut(const xAOD::Electron &ele,const xAOD::TrackParticle &trk1,const xAOD::TrackParticle &trk2){

    // calculate shower shapes and other discriminating variables
    float deltaEta1 = ele.trackCaloMatchValue(xAOD::EgammaParameters::TrackCaloMatchType::deltaEta1);

    double f1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f1);
//     double fSide = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::fracs1);
    double Eratio = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Eratio);
    double wTotS1 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::wtots1);

    double rEta = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Reta);
    double rPhi = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::Rphi);
    double wEta2 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::weta2);

    double f3 = ele.showerShapeValue(xAOD::EgammaParameters::ShowerShapeType::f3);

    // get pt and eta bins
    unsigned iPt = getPtBin(&ele);
    unsigned iEta = getEtaBin(&ele);


    // now apply the cuts

    // first cut: RHad
    // const std::vector<std::string> cutRHad({"<0.03", "<0.03", "<0.03", "<0.03", "<0.03", "<0.03", "<0.03", "<0.03", "<0.03", ""});
    static const std::vector<std::vector<float>> cutRHad({
      {999999, 999999, 999999, 999999,  999999, 999999, 999999, 999999, 999999, 999999}, // pT < 20 GeV
      {0.0490, 0.0470, 0.0450, 0.0400, 0.01980, 0.0500, 0.0500, 0.0500, 0.0450, 999999}, // pT [20, 25] GeV
      {0.0490, 0.0470, 0.0450, 0.0400, 0.01980, 0.0500, 0.0500, 0.0500, 0.0450, 999999}, // pT [25, 30] GeV
      {0.0440, 0.0420, 0.0400, 0.0350, 0.01820, 0.0450, 0.0450, 0.0450, 0.0400, 999999}, // pT [30, 35] GeV
      {0.0440, 0.0420, 0.0400, 0.0350, 0.01820, 0.0450, 0.0450, 0.0450, 0.0400, 999999}, // pT [35, 40] GeV
      {0.0365, 0.0345, 0.0325, 0.0300, 0.01680, 0.0375, 0.0375, 0.0375, 0.0325, 999999}, // pT [40, 45] GeV
      {0.0365, 0.0345, 0.0325, 0.0300, 0.01680, 0.0375, 0.0375, 0.0375, 0.0325, 999999}, // pT [45, 50] GeV
      {0.0340, 0.0320, 0.0300, 0.0275, 0.01680, 0.0350, 0.0350, 0.0350, 0.0300, 999999}, // pT [50, 60] GeV
      {0.0290, 0.0270, 0.0250, 0.0240, 0.01680, 0.0300, 0.0300, 0.0300, 0.0260, 999999}, // pT [60, 80] GeV
      {0.0270, 0.0250, 0.0230, 0.0220, 0.01680, 0.0280, 0.0280, 0.0280, 0.0250, 999999}, // pT > 80 GeV
    });
    if ( !(EleAcc::RhadForPID(ele) < cutRHad[iPt][iEta]) ) return(false);

    // second cut: f3
    // const std::vector<std::string> cutF3({"", "", "", "", "", ">0.0005", ">0.0005", ">0.0005", ">0.0005", ">0.0005", ">0.0005"});
    // if (!passCut(f3, cutF3[iPt])) return(false);
    static const std::vector<std::vector<float>> cutF3({
      { 999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999, 999999}, // pT < 20 GeV
      {0.01800, 0.01750, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [20, 25] GeV
      {0.01800, 0.01750, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [25, 30] GeV
      {0.01800, 0.01680, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [30, 35] GeV
      {0.01800, 0.01680, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [35, 40] GeV
      {0.01800, 0.01680, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [40, 45] GeV
      {0.01800, 0.01680, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [45, 50] GeV
      {0.01800, 0.01680, 0.01800, 0.02000,  999999, 0.01560, 0.02850, 0.02680, 0.03830, 999999}, // pT [50, 60] GeV
      { 999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999, 999999}, // pT [60, 80] GeV
      { 999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999,  999999, 999999}, // pT > 80 GeV
    });
    if ( !(f3 < cutF3[iPt][iEta]) ) return(false);

    // // seventh cut: f1
    // // const std::vector<std::string> cutF1({"", "", "", "", "", ">0.1", ">0.1", ">0.1", "", "", ""});
    // const std::string cutF1(">0.005");
    // if (!passCut(f1, cutF1)) return(false);

    // third cut: REta / eta
    static const std::vector<float> cutREtaInEta({0.920, 0.915, 0.910, 0.905, 0.85, 0.90, 0.90, 0.90, 0.90, -999999});
    if ( !(rEta > cutREtaInEta[iEta]) ) return(false);

    // // fourth cut: REta / pT
    // // const std::vector<std::string> cutREtaInPt({"", ">0.70", ">0.75", ">0.80", ">0.85", ">0.85", ">0.85", ">0.90", ">0.90", ">0.90"});
    // const std::vector<std::string> cutREtaInPt({"", ">0.90", ">0.90", ">0.90", ">0.90", ">0.90", ">0.90", ">0.90", ">0.90", ">0.90"});
    // if (!passCut(rEta, cutREtaInPt[iPt])){
    //   return(false);
    // }

    // fourth cut: ERatio / eta
    static const std::vector<float> cutERatioInEta({0.90, 0.90, 0.85, 0.80, -999999, 0.90, 0.90, 0.90, -999999, -999999});
    if (f1 > 0.005 && !(Eratio > cutERatioInEta[iEta]) ){
      if (EleAcc::dRExtrapTrk12(ele) > 0.1) return(false);
    }

    // fifth cut: wTotS1 / eta
    static const std::vector<float> cutWTotS1InEta({3.0, 3.1, 3.2, 3.3, 3.5, 3.5, 2.5, 1.8, 1.8, 999999});
    if (f1 > 0.005 && !(wTotS1 < cutWTotS1InEta[iEta]) ){
      return(false);
    }

    // sixth cut: wEtaS2 / eta
    static const std::vector<float> cutWEtaS2InEta({0.0115, 0.0120, 0.0125, 0.0130, 999999, 0.0130, 0.0130, 0.0140, 0.0140, 999999});
    if ( !(wEta2 < cutWEtaS2InEta[iEta]) ) return(false);

    // seventh cut: deltaEta1 / pT
    // const std::vector<std::string> cutDeltaEta1InPt({"", "", "", "", "", ">-0.003&<0.003", ">-0.003&<0.003", ">-0.003&<0.003", ">-0.003&<0.003", ">-0.003&<0.003", ">-0.003&<0.003"});
    // const std::vector<std::string> cutDeltaEta1InPt({">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002", ">-0.002&<0.002"});
    static const std::vector<float> cutDeltaEta1InPt({0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025});
    if (f1 > 0.005 && !(fabs(deltaEta1) < cutDeltaEta1InPt[iPt]) ) return(false);

    // eigth cut: RHad / eta
    static const std::vector<float> cutRHadInEta({0.02, 0.02, 0.02, 0.02, 999999, 0.02, 0.02, 0.02, 0.02, 999999});
    if ( !( EleAcc::RhadForPID(ele) < cutRHadInEta[iEta]) ) return(false);

    // ninth cut: FHT2 / eta
    static const std::vector<float> cutFHT2InEta({-0.16, -0.16, -0.16, -0.16, -0.16, -999999, -0.3, -999999, -999999, -999999});
    if ( !(TrkAcc::TRT_PID_trans(trk2) > cutFHT2InEta[iEta]) ) return(false);

    // tenth cut: FHT1 / eta
    static const std::vector<float> cutFHT1InEta({-0.16, -0.16, -0.16, -0.16, -0.16, -999999, -0.3, -999999, -999999, -999999});
    if ( !(TrkAcc::TRT_PID_trans(trk1) > cutFHT1InEta[iEta]) ) return(false);

    // eleventh cut: EOverP / eta
    // const std::vector<std::string> cutEOverPInPt({"", ">0.2", ">0.2", ">0.2", ">0.2", ">0.2", ">0.25", "", "", ""});
    static const std::vector<float> cutEOverPInEta({0.8, 0.8, 0.8, 0.8, -999999, -999999, -999999, -999999, -999999, -999999});
    if ( !(EleAcc::EOverP0P1(ele) > cutEOverPInEta[iEta]) ) return(false);

    // twelvth cut: rPhi / eta
    static const std::vector<float> cutRPhiInEta({0.84, 0.84, -999999, -999999, -999999, -999999, 0.80, 0.84, 0.84, -999999});
    if ( !(rPhi > cutRPhiInEta[iEta]) ) return(false);

    // // thirteenth cut: fSide / eta
    // const std::vector<std::string> cutFSideInEta({">0.15", ">0.175", ">0.2", ">0.25", "", ">0.2", "", "<0.3", "", ""});
    // if (!passCut(fSide, cutFSideInEta[iEta])) return(false);

    // thirteenth cut: d0SigmaTrk1
    // This cut is covered by passIPCut in TrackHandler.cxx
    // const std::vector<std::string> cutD0SigmaTrk1({"<5.0"});
    // if (!passCut(TrkAcc::d0significance(trk1), cutD0SigmaTrk1[0])) return(false);

    // fourteenth cut: d0SigmaTrk2
    // This cut is covered by passIPCut in TrackHandler.cxx
    // const std::vector<std::string> cutD0SigmaTrk2({"<5.0"});
    // if (!passCut(TrkAcc::d0significance(trk2), cutD0SigmaTrk2[0])) return(false);


    // fSide: if (f1 > 0.005 && )

    // congrats - this event survived all the cuts
    return(true);
}

unsigned HG::MergedElectronID::getPtBin(const xAOD::Electron * const el) const {

    double pt = el->pt() / 1000.0;

    if (pt < 20.0){return(0);}
    else if (pt < 25.0){return(1);}
    else if (pt < 30.0){return(2);}
    else if (pt < 35.0){return(3);}
    else if (pt < 40.0){return(4);}
    else if (pt < 45.0){return(5);}
    else if (pt < 50.0){return(6);}
    else if (pt < 60.0){return(7);}
    else if (pt < 80.0){return(8);}
    else {return(9);}

}

unsigned HG::MergedElectronID::getEtaBin(const xAOD::Electron * const el) const {

    float eta = fabs(el->caloCluster()->etaBE(2));

    if (eta < 0.6){return(0);}
    else if (eta < 0.8){return(1);}
    else if (eta < 1.15){return(2);}
    else if (eta < 1.37){return(3);}
    else if (eta < 1.52){return(4);}
    else if (eta < 1.81){return(5);}
    else if (eta < 2.01){return(6);}
    else if (eta < 2.37){return(7);}
    else if (eta < 2.47){return(8);}
    else {return(9);}

}
AngularPosition HG::MergedElectronID::getExtrapolatedTrackPosition(
    const xAOD::TrackParticle * track,
    const HG::MergedElectronID::extrapolationStartPositionEnum extrapolationStartPosition,
    const bool kalmanUpdate,
    const bool verbose = false,
    const bool printTrajectory = false){

    TrackModel::HelixParameters startParameters;
    // std::cout << "-------------------------------\n";

    if (extrapolationStartPosition == extrapolationStartPositionEnum::FirstMeasurement
        || extrapolationStartPosition == extrapolationStartPositionEnum::LastMeasurement){

      xAOD::ParameterPosition pos;

      // parameters at
      // - xAOD::BeamLine
      // - xAOD::CalorimeterEntrance
      // - xAOD::CalorimeterExit
      // are not available
      if (extrapolationStartPosition == extrapolationStartPositionEnum::FirstMeasurement)
        pos = xAOD::FirstMeasurement; // innermost ID layer with hit belonging to the track -> index = 0
      else if (extrapolationStartPosition == extrapolationStartPositionEnum::LastMeasurement)
        pos = xAOD::LastMeasurement;  // outermost ID layer with hit belonging to the track -> index = 1

      unsigned int index = 0;
      bool found = track->indexOfParameterAtPosition(index, pos);
      if (!found){
        Error("getTrackPosition()", "Could not find the track parameters");
        assert(false);
      }

      xAOD::CurvilinearParameters_t params = track->trackParameters(index);

      double x = params(0) / 1000.0; double px = params(3) / 1000.0;
      double y = params(1) / 1000.0; double py = params(4) / 1000.0;
      double z = params(2) / 1000.0; double pz = params(5) / 1000.0;

      if (verbose){
        std::cout << "track position: (x,  y,  z)  = " << Form("(%f, %f, %f)", x, y, z) << std::endl;
        std::cout << "track momentum: (px, py, pz) = " << Form("(%f, %f, %f)", px, py, pz) << std::endl;
      }

      double pT = sqrt(px*px + py*py);
      double qOverPt = track->charge() / pT;
      double cotTheta = pz / pT;
      double phi = atan2(py, px);
      // std::cout << "qOverPt = " << qOverPt << std::endl;
      // std::cout << "phi     = " << phi << std::endl;
      // std::cout << "cotTheta = " << cotTheta << std::endl;

      startParameters = TrackModel::HelixParameters(x, y, z, phi, cotTheta, qOverPt);


    } else if (extrapolationStartPosition == extrapolationStartPositionEnum::Perigee){

      xAOD::DefiningParameters_t trkPar = track->definingParameters();

      // old parameters
      double d0 = trkPar(0) / 1000.0;      // std::cout << "d0      = " << d0 << std::endl;
      double z0 = trkPar(1) / 1000.0;      // std::cout << "z0      = " << z0 << std::endl;
      double phi = trkPar(2);              // std::cout << "phi     = " << phi << std::endl;
      double theta = trkPar(3);            // std::cout << "theta   = " << theta << std::endl;
      double qOverPt = trkPar(4) * 1000.0; // std::cout << "qOverPt = " << qOverPt << std::endl;


      if (kalmanUpdate){

        assert(false);

        // std::cout << "beamPosX = " << eventInfo->beamPosX() << std::endl;
        // std::cout << "beamPosY = " << eventInfo->beamPosY() << std::endl;
        // std::cout << "beamPosZ = " << eventInfo->beamPosZ() << std::endl;
        // std::cout << "beamPosSigmaX = " << eventInfo->beamPosSigmaX() << std::endl;
        // std::cout << "beamPosSigmaY = " << eventInfo->beamPosSigmaY() << std::endl;
        // std::cout << "beamPosSigmaZ = " << eventInfo->beamPosSigmaZ() << std::endl;

//lines below commented out since kalmanUpdate is not used -- to enable it, should also pass (a pointer to) xAOD::eventInfo here
//         double sx = eventInfo->beamPosSigmaX();
//         double sy = eventInfo->beamPosSigmaY();
//
//         const double measPar = 0.0;             // Constrain d0 to 0
//         // const double measCov = 1e-12;        // Uncertainty on d0 1 nm -- something small
//         const double measCov = 0.5 * (sx + sy);
//
//         const int mk = 0; // d0 is the first entry in the track parameter matrix
//         double r = measPar - trkPar(mk);
//         double R = measCov + trkCov(mk,mk); R = 1./R;
//
//         // compute updated state, here = TP + K * r = TP + TCov*H.T*R * r
//         AmgVector(5) trkParNew = trkPar + trkCov.col(mk) * R * r;
//
//         d0 = trkParNew(0) / 1000.0;      // std::cout << "d0      = " << d0 << std::endl;
//         z0 = trkParNew(1) / 1000.0;      // std::cout << "z0      = " << z0 << std::endl;
//         phi = trkParNew(2);              // std::cout << "phi     = " << phi << std::endl;
//         theta = trkParNew(3);            // std::cout << "theta   = " << theta << std::endl;
//         qOverPt = trkParNew(4) * 1000.0; // std::cout << "qOverPt = " << qOverPt << std::endl;
//         pt = track->charge() / qOverPt;

      }

      double x =   d0 * cos(phi);
      double y = - d0 * sin(phi);
      double cotTheta = cos(theta) / sin(theta);

      if (verbose){
        std::cout << "x        = " << x << std::endl;
        std::cout << "y        = " << y << std::endl;
        std::cout << "cotTheta = " << cotTheta << std::endl;
      }

      startParameters = TrackModel::HelixParameters(x, y, z0, phi, cotTheta, qOverPt);

    }

    TrackModel model(startParameters);
    TrackModel::HelixState startState = model.determineInitialHelixState(startParameters);

    // set extrapolation targets
    double targetR, targetZ;
    // targetR  = 1.150; //Solenoid R (in m)
    targetR  = 1.500; //ECalStripLayer R (in m)
    targetZ  = 3.512; //ID endplate (in m)

    TrackModel::HelixState endState = model.stateAtBoundaryRZ(targetR, targetZ, printTrajectory);

    if (verbose){
      std::cout << Form("Position   (%2.2f,%2.2f,%2.2f) -> (%2.2f,%2.2f,%2.2f)",
                        startState.position.x, startState.position.y, startState.position.z,
                        endState.position.x, endState.position.y, endState.position.z ) << std::endl;
      std::cout << Form("Momentum   (%2.2f,%2.2f,%2.2f) -> (%2.2f,%2.2f,%2.2f)",
                        startState.momentum.x, startState.momentum.y, startState.momentum.z,
                        endState.momentum.x, endState.momentum.y, endState.momentum.z ) << std::endl;
    }

    return(AngularPosition(endState.position.eta, endState.position.phi));
}

//______________________________________________________________________________
bool HG::MergedElectronID::passPreselection(const xAOD::Electron &ele,
                                            const xAOD::TrackParticle &trk1,
                                            const xAOD::TrackParticle &trk2){

  if (ele.pt() < m_mergedElePtCut) return false;
  if (fabs(ele.caloCluster()->etaBE(2)) > m_mergedEleEtaCut) return false;

  if (HG::EleAcc::RhadForPID(ele) > m_PreselRhad) return false;


  int nPassBlayer = 0;
  nPassBlayer += HG::TrkAcc::passBLayerRequirement(trk1) ? 1 : 0;
  nPassBlayer += HG::TrkAcc::passBLayerRequirement(trk2) ? 1 : 0;
  if (nPassBlayer < m_PreselNPassBlayer) return false;

  return true;
}
