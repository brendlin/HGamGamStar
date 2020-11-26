#ifndef HGamGamStar_HggStarVariables_H
#define HGamGamStar_HggStarVariables_H

#include "HGamAnalysisFramework/VarHandler.h"
#include "HGamAnalysisFramework/HGamCommon.h"
#include "HGamAnalysisFramework/TruthUtils.h"

#include "HGamGamStar/HggStarCommon.h"
#include "HGamGamStar/ExtraHggStarObjects.h"

#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"
#include "FourMomUtils/xAODP4Helpers.h"

namespace HG {

  TruthPtcls getHyyStarSignalDecayProducts(const xAOD::TruthParticle *ptcl);
  TruthPtcls FilterLeptons(const TruthPtcls& stableHiggsDecayProducts);
  TruthPtcls FilterDirectPhotons(const TruthPtcls& stableHiggsDecayProducts);
  void SetMergedFourMomentum(xAOD::Electron& ele,const float calibrated_e);

  //____________________________________________________________________________
  class m_lly : public VarBase<float> {
  public:
  m_lly() : VarBase("m_lly") { m_default = -99; }
    ~m_lly() { }

    float calculateValue(bool truth)
    {
      // For Reco:
      // getElectrons and getMuons only return elecs / muons selected as the candidate y*.
      // getPhotons only returns the leading photon candidate.
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      // >= 2 muons and >= 1 photon
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).M();

      // >= 2 electrons and >= 1 photon
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).M();

      // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
      if (!truth && eles->size() == 1)
        return ((*gams)[0]->p4() + (*eles)[0]->p4()).M();

      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_lly_gev : public VarBase<float> {
  public:
  m_lly_gev() : VarBase("m_lly_gev") { m_default = -99; }
    ~m_lly_gev() { }

    float calculateValue(bool truth);
    // Implemented in the cxx file.
  };

  //____________________________________________________________________________
  class m_ll : public VarBase<float> {
  public:
  m_ll() : VarBase("m_ll") { m_default = -99; }
    ~m_ll() { }

    float calculateValue(bool truth)
    {
      // For Reco:
      // getElectrons and getMuons only return elecs / muons selected as the candidate y*.
      // For Truth: best to use "m_yStar_undressed_h1"
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4()).M();
      if (eles->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4()).M();

      // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
      if (!truth && eles->size() == 1)
        return (*eles)[0]->p4().M();

      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_l1y : public VarBase<float> {
  public:
  m_l1y() : VarBase("m_l1y") { m_default = -99; }
    ~m_l1y() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*gams)[0]->p4()).M();
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*gams)[0]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_l2y : public VarBase<float> {
  public:
  m_l2y() : VarBase("m_l2y") { m_default = -99; }
    ~m_l2y() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[1]->p4() + (*gams)[0]->p4()).M();
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[1]->p4() + (*gams)[0]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class deltaR_ll : public VarBase<float> {
  public:
  deltaR_ll() : VarBase("deltaR_ll") { m_default = -99; }
    ~deltaR_ll() { }

    float calculateValue(bool truth)
    {
      // For Reco:
      // getElectrons and getMuons only return elecs / muons selected as the candidate y*.
      // For Truth: best to use "m_yStar_undressed_h1"
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return (*mus)[0]->p4().DeltaR((*mus)[1]->p4());
      if (eles->size() >= 2)
        return (*eles)[0]->p4().DeltaR((*eles)[1]->p4());
      return m_default;
    }
  };

  //____________________________________________________________________________
  class Resolved_deltaPhiRescaled2 : public VarBase<float> {
  public:
  Resolved_deltaPhiRescaled2() : VarBase("Resolved_deltaPhiRescaled2") { m_default = -99; m_recoOnly = true; }
    ~Resolved_deltaPhiRescaled2() { }

    float calculateValue(bool truth)
    {
      // This is the dR of the tracks at last measurement extrapolated to the EM calorimeter.
      // This is dR between the cluster barycenter from 3rd sampling (etaBE2, phiBE2)
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);
      if (eles->size() >= 2) {

        // This is to fix how the shower shape variables are saved.
        int flipSign1 = ( (*eles)[0]->trackParticle()->charge() > 0) ? 1 : -1; // -1 = flip
        int flipSign2 = ( (*eles)[1]->trackParticle()->charge() > 0) ? 1 : -1; // -1 = flip

        float dphi_trk1,dphi_trk2;
        // Christos says that track-matching uses DeltaPhiRescaled, since this matches barycenter better
        // and that Layer 2 is best.
        (*eles)[0]->trackCaloMatchValue(dphi_trk1, xAOD::EgammaParameters::deltaPhiRescaled2);
        (*eles)[1]->trackCaloMatchValue(dphi_trk2, xAOD::EgammaParameters::deltaPhiRescaled2);
        float dphi_e1e2 = xAOD::P4Helpers::deltaPhi((*eles)[0]->caloCluster()->phiBE(2) + dphi_trk1*flipSign1,
                                                    (*eles)[1]->caloCluster()->phiBE(2) + dphi_trk2*flipSign2);

        // This is to get the higher-pt tracks all bending in the same direction.
        int flipSignOverall = ((*eles)[0]->trackParticle()->pt() > (*eles)[1]->trackParticle()->pt()) ? flipSign1 : flipSign2;
        return flipSignOverall * dphi_e1e2;
      }
      return m_default;
    }
  };

  //____________________________________________________________________________
  class Resolved_deltaEta2 : public VarBase<float> {
  public:
  Resolved_deltaEta2() : VarBase("Resolved_deltaEta2") { m_default = -99; m_recoOnly = true; }
    ~Resolved_deltaEta2() { }

    float calculateValue(bool truth)
    {
      // This is the dR of the tracks at last measurement extrapolated to the EM calorimeter.
      // This is dR between the cluster barycenter from 3rd sampling (etaBE2, phiBE2)
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);
      if (eles->size() >= 2) {
        float deta_trk1,deta_trk2;
        (*eles)[0]->trackCaloMatchValue(deta_trk1, xAOD::EgammaParameters::deltaEta2);
        (*eles)[1]->trackCaloMatchValue(deta_trk2, xAOD::EgammaParameters::deltaEta2);
        float deta_e1e2 = ( ((*eles)[0]->caloCluster()->etaBE(2) - deta_trk1) -
                            ((*eles)[1]->caloCluster()->etaBE(2) - deta_trk2) );
        return deta_e1e2;
      }
      return m_default;
    }
  };

  //____________________________________________________________________________
  class Resolved_dRExtrapTrk12 : public VarBase<float> {
  public:
  Resolved_dRExtrapTrk12() : VarBase("Resolved_dRExtrapTrk12") { m_default = -99; m_recoOnly = true; }
    ~Resolved_dRExtrapTrk12() { }

    float calculateValue(bool truth); // See cxx file
  };

  //____________________________________________________________________________
  class deltaPhi_trktrk_IP : public VarBase<float> {
  public:
  deltaPhi_trktrk_IP() : VarBase("deltaPhi_trktrk_IP") { m_default = -99; }
    ~deltaPhi_trktrk_IP() { }

    float calculateValue(bool truth)
    {
      if (truth) {
        const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
        if (childleps->size() != 2) return m_default;

        if ( (*childleps)[0]->pt() > (*childleps)[1]->pt())
          return (*childleps)[0]->charge() * xAOD::P4Helpers::deltaPhi((*childleps)[0],(*childleps)[1]);
        return   (*childleps)[1]->charge() * xAOD::P4Helpers::deltaPhi((*childleps)[1],(*childleps)[0]);
      }

      const xAOD::MuonContainer *mus = (xAOD::MuonContainer*)HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return (*mus)[0]->charge() * xAOD::P4Helpers::deltaPhi((*mus)[0],(*mus)[1]);

      const xAOD::TrackParticleContainer *eleTrks = (xAOD::TrackParticleContainer*)HG::ExtraHggStarObjects::getInstance()->getElectronTracks(truth);
      if (eleTrks->size() >= 2) {

        TLorentzVector tlv1 = (*eleTrks)[0]->p4();
        TLorentzVector tlv2 = (*eleTrks)[1]->p4();
        tlv1.SetPtEtaPhiM( tlv1.Pt(), tlv1.Eta(), tlv1.Phi(), 0.510998 ); // ele.m == 0.510998
        tlv2.SetPtEtaPhiM( tlv2.Pt(), tlv2.Eta(), tlv2.Phi(), 0.510998 );

        return (*eleTrks)[0]->charge() * tlv1.DeltaPhi( tlv2 );
      }

      return m_default;
    }
  };

  //____________________________________________________________________________
  class deltaEta_trktrk_IP : public VarBase<float> {
  public:
  deltaEta_trktrk_IP() : VarBase("deltaEta_trktrk_IP") { m_default = -99; m_recoOnly = true; }
    ~deltaEta_trktrk_IP() { }

    float calculateValue(bool truth)
    {
      const xAOD::MuonContainer *mus = (xAOD::MuonContainer*)HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return xAOD::P4Helpers::deltaEta((*mus)[0],(*mus)[1]);

      const xAOD::TrackParticleContainer *eleTrks = (xAOD::TrackParticleContainer*)HG::ExtraHggStarObjects::getInstance()->getElectronTracks(truth);
      if (eleTrks->size() >= 2) {

        TLorentzVector tlv1 = (*eleTrks)[0]->p4();
        TLorentzVector tlv2 = (*eleTrks)[1]->p4();
        tlv1.SetPtEtaPhiM( tlv1.Pt(), tlv1.Eta(), tlv1.Phi(), 0.510998 ); // ele.m == 0.510998
        tlv2.SetPtEtaPhiM( tlv2.Pt(), tlv2.Eta(), tlv2.Phi(), 0.510998 );

        return tlv1.Eta() - tlv2.Eta();
      }

      return m_default;
    }
  };

  //____________________________________________________________________________
  class deltaR_track4mom : public VarBase<float> {
  public:
  deltaR_track4mom() : VarBase("deltaR_track4mom") { m_default = -99; m_recoOnly = true; }
    ~deltaR_track4mom() { }

    float calculateValue(bool truth); // See cxx file
  };

  //____________________________________________________________________________
  class deltaPhi2_trktrk_perigee : public VarBase<float> {
  public:
  deltaPhi2_trktrk_perigee() : VarBase("deltaPhi2_trktrk_perigee") { m_default = -99; m_recoOnly = true; }
    ~deltaPhi2_trktrk_perigee() { }

    float calculateValue(bool truth)
    {
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);

      const xAOD::Electron* ele1 = NULL;
      const xAOD::Electron* ele2 = NULL;
      int ele1_trkIndex = -1;
      int ele2_trkIndex = -1;
      if (eles->size() >= 2) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[1];
        ele1_trkIndex = 0;
        ele2_trkIndex = 0;
      }
      else if (eles->size() == 1) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[0];
        ele1_trkIndex = EleAcc::vtxTrkIndex1(*ele1);
        ele2_trkIndex = EleAcc::vtxTrkIndex2(*ele2);
      }
      else {
        return m_default;
      }

      const xAOD::TrackParticle* ele_tp1 = ele1->trackParticle(ele1_trkIndex);
      const xAOD::TrackParticle* ele_tp2 = ele2->trackParticle(ele2_trkIndex);

      if (!ele_tp1 || !ele_tp2) return m_default;

      // This is to fix how the shower shape variables are saved.
      int flipSign1 = ( ele_tp1->charge() > 0) ? 1 : -1; // -1 = flip
      int flipSign2 = ( ele_tp2->charge() > 0) ? 1 : -1; // -1 = flip

      // This is to get the higher-pt tracks all bending in the same direction.
      float dphi_trk1 = HG::EleAcc::TrackMatchingP_dPhi2( *ele1 )[ele1_trkIndex];
      float dphi_trk2 = HG::EleAcc::TrackMatchingP_dPhi2( *ele2 )[ele2_trkIndex];
      float dphi_e1e2 = xAOD::P4Helpers::deltaPhi(ele1->caloCluster()->phiBE(2) + dphi_trk1*flipSign1,
                                                  ele2->caloCluster()->phiBE(2) + dphi_trk2*flipSign2);

      int flipSignOverall = (ele_tp1->pt() > ele_tp2->pt()) ? flipSign1 : flipSign2;
      return flipSignOverall * dphi_e1e2;
    }
  };

  //____________________________________________________________________________
  class deltaEta2_trktrk_perigee : public VarBase<float> {
  public:
  deltaEta2_trktrk_perigee() : VarBase("deltaEta2_trktrk_perigee") { m_default = -99; m_recoOnly = true; }
    ~deltaEta2_trktrk_perigee() { }

    float calculateValue(bool truth)
    {
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);

      const xAOD::Electron* ele1 = NULL;
      const xAOD::Electron* ele2 = NULL;
      int ele1_trkIndex = -1;
      int ele2_trkIndex = -1;
      if (eles->size() >= 2) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[1];
        ele1_trkIndex = 0;
        ele2_trkIndex = 0;
      }
      else if (eles->size() == 1) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[0];
        ele1_trkIndex = EleAcc::vtxTrkIndex1(*ele1);
        ele2_trkIndex = EleAcc::vtxTrkIndex2(*ele2);
      }
      else {
        return m_default;
      }

      // This is to get the higher-pt tracks all bending in the same direction.
      float deta_trk1 = HG::EleAcc::TrackMatchingP_dEta2( *ele1 )[ele1_trkIndex];
      float deta_trk2 = HG::EleAcc::TrackMatchingP_dEta2( *ele2 )[ele2_trkIndex];
      float deta_e1e2 = ( (ele1->caloCluster()->etaBE(2) - deta_trk1) -
                          (ele2->caloCluster()->etaBE(2) - deta_trk2) );

      return deta_e1e2;
    }
  };

  //____________________________________________________________________________
  class deltaRL2_trktrk_perigee : public VarBase<float> {
  public:
  deltaRL2_trktrk_perigee() : VarBase("deltaRL2_trktrk_perigee") { m_default = -99; m_recoOnly = true; }
    ~deltaRL2_trktrk_perigee() { }

    float calculateValue(bool truth); // See cxx file
  };
  
  //____________________________________________________________________________
  class passVBFpresel : public VarBase<bool> {
  public:
  passVBFpresel() : VarBase("passVBFpresel") { m_default = false; m_recoOnly = true; }
    ~passVBFpresel() { }

    bool calculateValue(bool truth); // See cxx file
  };

  //____________________________________________________________________________
  class deltaPhi2_trktrk_LM : public VarBase<float> {
  public:
  deltaPhi2_trktrk_LM() : VarBase("deltaPhi2_trktrk_LM") { m_default = -99; m_recoOnly = true; }
    ~deltaPhi2_trktrk_LM() { }

    float calculateValue(bool truth)
    {
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);

      const xAOD::Electron* ele1 = NULL;
      const xAOD::Electron* ele2 = NULL;
      int ele1_trkIndex = -1;
      int ele2_trkIndex = -1;
      if (eles->size() >= 2) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[1];
        ele1_trkIndex = 0;
        ele2_trkIndex = 0;
      }
      else if (eles->size() == 1) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[0];
        ele1_trkIndex = EleAcc::vtxTrkIndex1(*ele1);
        ele2_trkIndex = EleAcc::vtxTrkIndex2(*ele2);
      }
      else {
        return m_default;
      }

      const xAOD::TrackParticle* ele_tp1 = ele1->trackParticle(ele1_trkIndex);
      const xAOD::TrackParticle* ele_tp2 = ele2->trackParticle(ele2_trkIndex);

      if (!ele_tp1 || !ele_tp2) return m_default;

      // This is to fix how the shower shape variables are saved.
      int flipSign1 = ( ele_tp1->charge() > 0) ? 1 : -1; // -1 = flip
      int flipSign2 = ( ele_tp2->charge() > 0) ? 1 : -1; // -1 = flip

      // This is to get the higher-pt tracks all bending in the same direction.
      float dphi_trk1 = HG::EleAcc::TrackMatchingLM_dPhi2( *ele1 )[ele1_trkIndex];
      float dphi_trk2 = HG::EleAcc::TrackMatchingLM_dPhi2( *ele2 )[ele2_trkIndex];
      float dphi_e1e2 = xAOD::P4Helpers::deltaPhi(ele1->caloCluster()->phiBE(2) + dphi_trk1*flipSign1,
                                                  ele2->caloCluster()->phiBE(2) + dphi_trk2*flipSign2);

      int flipSignOverall = (ele_tp1->pt() > ele_tp2->pt()) ? flipSign1 : flipSign2;
      return flipSignOverall * dphi_e1e2;
    }
  };

  //____________________________________________________________________________
  class deltaEta2_trktrk_LM : public VarBase<float> {
  public:
  deltaEta2_trktrk_LM() : VarBase("deltaEta2_trktrk_LM") { m_default = -99; m_recoOnly = true; }
    ~deltaEta2_trktrk_LM() { }

    float calculateValue(bool truth)
    {
      const xAOD::ElectronContainer *eles = (xAOD::ElectronContainer*)HG::VarHandler::getInstance()->getElectrons(truth);

      const xAOD::Electron* ele1 = NULL;
      const xAOD::Electron* ele2 = NULL;
      int ele1_trkIndex = -1;
      int ele2_trkIndex = -1;
      if (eles->size() >= 2) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[1];
        ele1_trkIndex = 0;
        ele2_trkIndex = 0;
      }
      else if (eles->size() == 1) {
        ele1 = (*eles)[0];
        ele2 = (*eles)[0];
        ele1_trkIndex = EleAcc::vtxTrkIndex1(*ele1);
        ele2_trkIndex = EleAcc::vtxTrkIndex2(*ele2);
      }
      else {
        return m_default;
      }

      // This is to get the higher-pt tracks all bending in the same direction.
      float deta_trk1 = HG::EleAcc::TrackMatchingLM_dEta2( *ele1 )[ele1_trkIndex];
      float deta_trk2 = HG::EleAcc::TrackMatchingLM_dEta2( *ele2 )[ele2_trkIndex];
      float deta_e1e2 = ( (ele1->caloCluster()->etaBE(2) - deta_trk1) -
                          (ele2->caloCluster()->etaBE(2) - deta_trk2) );

      return deta_e1e2;
    }
  };

  //____________________________________________________________________________
  class deltaRL2_trktrk_LM : public VarBase<float> {
  public:
  deltaRL2_trktrk_LM() : VarBase("deltaRL2_trktrk_LM") { m_default = -99; m_recoOnly = true; }
    ~deltaRL2_trktrk_LM() { }

    float calculateValue(bool truth); // See cxx file
  };

  //____________________________________________________________________________
  class deltaPhi_naiveExtrap : public VarBase<float> {
  public:
  deltaPhi_naiveExtrap() : VarBase("deltaPhi_naiveExtrap") { m_default = -99; m_recoOnly = true;}
    ~deltaPhi_naiveExtrap() { }

    float calculateValue(bool truth)
    {
      const xAOD::TrackParticleContainer *eleTrks = (xAOD::TrackParticleContainer*)HG::ExtraHggStarObjects::getInstance()->getElectronTracks(truth);

      if (eleTrks->size() < 2) return m_default;

      const xAOD::TrackParticle* trk1 = (*eleTrks)[0];
      const xAOD::TrackParticle* trk2 = (*eleTrks)[1];

      float loopRadius_trk1 = 1.668*trk1->pt()/1000.;
      //float etaToAngle1 = M_PI/2. - 2*atan( exp( -fabs(trk1->eta()) ) );
      float phiMag2T_trk1 = ( (fabs(trk1->eta()) < 1.52) ?
                              // 1.500m is the radius of the barrel calorimeter
                              // 1.210m is where the magnet coil is. Not perfect... but close enough?
                              asin(0.750/loopRadius_trk1) :
                              // 3.512m is the start of the ID endplate (in m)
                              // 3.680m is my guess based on another source...
                              // NOT 4668.5mm is the front face (cold) of FCAL1 (endcap EM)
                              // The solenoid stops at 2.642
                              0.5*(2.642 / loopRadius_trk1) * (trk1->pt() / fabs(trk1->p4().Pz())) );
                              /* 0.5*(3.512 / loopRadius_trk1) * tan( etaToAngle1 ) ); */

      float loopRadius_trk2 = 1.668*trk2->pt()/1000.;
      //float etaToAngle2 = M_PI/2. - 2*atan( exp( -fabs(trk2->eta()) ) );
      float phiMag2T_trk2 = ( (fabs(trk2->eta()) < 1.52) ?
                              asin(0.750/loopRadius_trk2) :
                              0.5*(2.642 / loopRadius_trk2) * (trk2->pt() / fabs(trk2->p4().Pz())) );

      float deltaphi = xAOD::P4Helpers::deltaPhi( trk1->phi() - trk1->charge() * phiMag2T_trk1,
                                                  trk2->phi() - trk2->charge() * phiMag2T_trk2 );
      int flipSignOverall = (trk1->pt() > trk2->pt()) ? trk1->charge() : trk2->charge();
      return flipSignOverall * deltaphi;
    }
  };

  //____________________________________________________________________________
  class deltaPhi_overScaled_naiveExtrap : public VarBase<float> {
  public:
  deltaPhi_overScaled_naiveExtrap() : VarBase("deltaPhi_overScaled_naiveExtrap") { m_default = -99; m_recoOnly = true;}
    ~deltaPhi_overScaled_naiveExtrap() { }

    float calculateValue(bool truth)
    {
      const xAOD::TrackParticleContainer *eleTrks = (xAOD::TrackParticleContainer*)HG::ExtraHggStarObjects::getInstance()->getElectronTracks(truth);

      if (eleTrks->size() < 2) return m_default;

      const xAOD::TrackParticle* trk1 = (*eleTrks)[0];
      const xAOD::TrackParticle* trk2 = (*eleTrks)[1];

      float largerPt = (trk1->pt() > trk2->pt()) ? trk1->pt() : trk2->pt();

      float loopRadius_trk1 = 1.668*largerPt/1000.;
      float phiMag2T_trk1 = ( (fabs(trk1->eta()) < 1.52) ?
                              asin(0.750/loopRadius_trk1) :
                              0.5*(2.642 / loopRadius_trk1) * (largerPt / fabs(trk1->p4().Pz())) );

      float loopRadius_trk2 = 1.668*largerPt/1000.;
      float phiMag2T_trk2 = ( (fabs(trk2->eta()) < 1.52) ?
                              asin(0.750/loopRadius_trk2) :
                              0.5*(2.642 / loopRadius_trk2) * (largerPt / fabs(trk2->p4().Pz())) );

      float deltaphi = xAOD::P4Helpers::deltaPhi( trk1->phi() - trk1->charge() * phiMag2T_trk1,
                                                  trk2->phi() - trk2->charge() * phiMag2T_trk2 );
      int flipSignOverall = (trk1->pt() > trk2->pt()) ? trk1->charge() : trk2->charge();
      return flipSignOverall * deltaphi;
    }
  };

  //____________________________________________________________________________
  class pt_lly : public VarBase<float> {
  public:
  pt_lly() : VarBase("pt_lly") { m_default = -99; }
    ~pt_lly() { }

    float calculateValue(bool truth)
    {
      // For Reco:
      // getElectrons and getMuons only return elecs / muons selected as the candidate y*.
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).Pt();
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).Pt();

      // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
      if (!truth && eles->size() == 1 && gams->size() >= 1)
        return ((*gams)[0]->p4() + (*eles)[0]->p4()).Pt();

      return m_default;
    }
  };

  //____________________________________________________________________________
  class pt_ll : public VarBase<float> {
  public:
  pt_ll() : VarBase("pt_ll") { m_default = -99; }
    ~pt_ll() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4()).Pt();
      if (eles->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4()).Pt();

      // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
      if (!truth && eles->size() == 1)
        return (*eles)[0]->p4().Pt();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_lly_track4mom : public VarBase<float> {
  public:
  m_lly_track4mom() : VarBase("m_lly_track4mom") { m_default = -99; m_recoOnly = true; }
    ~m_lly_track4mom() { }

    float calculateValue(bool truth)
    {
      if (truth)
      { return m_default; }

      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (gams->size() < 1) return m_default;

      const xAOD::IParticleContainer *trks = ExtraHggStarObjects::getInstance()->getElectronTracks();
      if (trks->size() < 2) return m_default;

      TLorentzVector tlv1 = (*trks)[0]->p4();
      TLorentzVector tlv2 = (*trks)[1]->p4();
      tlv1.SetPtEtaPhiM( tlv1.Pt(), tlv1.Eta(), tlv1.Phi(), 0.510998 ); // ele.m == 0.510998
      tlv2.SetPtEtaPhiM( tlv2.Pt(), tlv2.Eta(), tlv2.Phi(), 0.510998 );

      return (tlv1 + tlv2 + (*gams)[0]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class m_llyy : public VarBase<float> {
  public:
  m_llyy() : VarBase("m_llyy") { m_default = -99; }
    ~m_llyy() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4() + (*gams)[1]->p4()).M();
      if (eles->size() >= 2 && gams->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4() + (*gams)[1]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_ll_track4mom : public VarBase<float> {
  public:
  m_ll_track4mom() : VarBase("m_ll_track4mom") { m_default = -99; m_recoOnly = true; }
    ~m_ll_track4mom() { }

    float calculateValue(bool truth)
    {
      if (truth)
      { return m_default; }

      const xAOD::IParticleContainer *trks = ExtraHggStarObjects::getInstance()->getElectronTracks();
      if (trks->size() < 2) return m_default;

      TLorentzVector tlv1 = (*trks)[0]->p4();
      TLorentzVector tlv2 = (*trks)[1]->p4();
      tlv1.SetPtEtaPhiM( tlv1.Pt(), tlv1.Eta(), tlv1.Phi(), 0.510998 ); // ele.m == 0.510998
      tlv2.SetPtEtaPhiM( tlv2.Pt(), tlv2.Eta(), tlv2.Phi(), 0.510998 );

      return (tlv1 + tlv2).M();
    }
  };
  
  //____________________________________________________________________________
  class trk_lead_pt : public VarBase<float> {
  public:
    trk_lead_pt() : VarBase("trk_lead_pt") { m_default = -99; m_recoOnly = true; }
    ~trk_lead_pt() { }
    
    float calculateValue(bool truth)
    {
      if (truth)
      { return m_default; }

      const xAOD::IParticleContainer *trks = ExtraHggStarObjects::getInstance()->getElectronTracks();
      if (!trks->size()) return m_default;

      return (*trks)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class pt_llyy : public VarBase<float> {
  public:
  pt_llyy() : VarBase("pt_llyy") { m_default = -99; }
    ~pt_llyy() { }

    float calculateValue(bool truth)
    {
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4() + (*gams)[1]->p4()).Pt();
      if (eles->size() >= 2 && gams->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4() + (*gams)[1]->p4()).Pt();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class deltaPhi_ll_y : public VarBase<float> {
  public:
  deltaPhi_ll_y() : VarBase("deltaPhi_ll_y") { m_default = -99; }
    ~deltaPhi_ll_y() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*mus)[1]->p4()).DeltaPhi( (*gams)[0]->p4() );
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*eles)[1]->p4()).DeltaPhi( (*gams)[0]->p4() );
      return m_default;
    }
  };

  //____________________________________________________________________________
  class eta_y1 : public VarBase<float> {
  public:
    eta_y1() : VarBase("eta_y1") { m_default = -99; }
    ~eta_y1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1)
      { return m_default; }

      return (*gams)[0]->eta();
    }
  };

  //____________________________________________________________________________
  class m_emu : public VarBase<float> {
  public:
  m_emu() : VarBase("m_emu") { m_default = -99; }
    ~m_emu() { }

    float calculateValue(bool truth)
    {
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 1 && eles->size() >=1)
        return ((*eles)[0]->p4() + (*mus)[0]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_emuy : public VarBase<float> {
  public:
  m_emuy() : VarBase("m_emuy") { m_default = -99; }
    ~m_emuy() { }

    float calculateValue(bool truth)
    {
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 1 && eles->size() >=1 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*mus)[0]->p4() + (*gams)[0]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class Dy_j_j_50 : public VarBase<float> {
  public:
    Dy_j_j_50() : VarBase("Dy_j_j_50") { m_default = -99; }
    ~Dy_j_j_50() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);

      if (jets->size() < 2)
      { return m_default; }

      return (*jets)[1]->pt() < 50 * HG::GeV ? m_default : fabs((*jets)[0]->rapidity() - (*jets)[1]->rapidity());
    }
  };

  //____________________________________________________________________________
  class pT_l1 : public VarBase<float> {
  public:
  pT_l1() : VarBase("pT_l1") { m_default = -99; }
    ~pT_l1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return ((*mus)[0]->pt());
      if (eles->size() >= 2)
        return ((*eles)[0]->pt());
      return m_default;
    }
  };

  //____________________________________________________________________________
  class eta_j1 : public VarBase<float> {
  public:
    eta_j1() : VarBase("eta_j1") { m_default = -99; }
    ~eta_j1() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);

      if (jets->size() < 1)
      { return m_default; }

      return (*jets)[0]->eta();
    }
  };

  //____________________________________________________________________________
  class N_j_gap : public VarBase<int> {
  public:
    N_j_gap() : VarBase("N_j_gap") { m_default = -99; }
    ~N_j_gap() { }

    int calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      int ngapjets = 0;
      int njets = jets->size();
      if (njets >= 2) {
        float j1_eta = (*jets)[0]->eta();
        float j2_eta = (*jets)[1]->eta(); 
        for (int i = 2; i < njets; i++) {
          if ( ((*jets)[i]->eta() < j1_eta && (*jets)[i]->eta() > j2_eta) || ((*jets)[i]->eta() > j1_eta && (*jets)[i]->eta() < j2_eta) ) { ngapjets++; }
        }
      }

      return ngapjets;
    }
  };

  //____________________________________________________________________________
  class Zy_centrality: public VarBase<float> {
  public:
  Zy_centrality(): VarBase("Zy_centrality") {m_default = -99; }
  ~Zy_centrality() { }
    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1 && jets->size() >= 2)
        return abs((((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).Eta() - ((*jets)[0]->rapidity() + (*jets)[1]->rapidity())/2)/((*jets)[0]->rapidity() - (*jets)[1]->rapidity()));
      if (eles->size() >= 2 && gams->size() >= 1 && jets->size() >= 2)
        return abs((((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).Eta() - ((*jets)[0]->rapidity() + (*jets)[1]->rapidity())/2)/((*jets)[0]->rapidity() - (*jets)[1]->rapidity()));
      return m_default;
    }
  }; 

  //____________________________________________________________________________
  class DR_Zy_jj : public VarBase<float> {
  public:
    DR_Zy_jj() : VarBase("DR_Zy_jj") { m_default = -99; }
    ~DR_Zy_jj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);

      if (jets->size() < 2)
      { return m_default; }

      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1) return m_default;

      if (mus->size() >= 2){
        return fabs(((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).DeltaR((*jets)[0]->p4() + (*jets)[1]->p4()));
      }
      if (eles->size() >= 2){
        return fabs(((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).DeltaR((*jets)[0]->p4() + (*jets)[1]->p4()));
      }
      return m_default;
    }

  };


  //____________________________________________________________________________
  class pT_l1_h1 : public VarBase<float> {
  public:
  pT_l1_h1() : VarBase("pT_l1_h1") { m_default = -99; m_truthOnly = true; }
    ~pT_l1_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
      if (childleps->size() != 2) return m_default;

      return (*childleps)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class pT_l2_h1 : public VarBase<float> {
  public:
  pT_l2_h1() : VarBase("pT_l2_h1") { m_default = -99; m_truthOnly = true; }
    ~pT_l2_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
      if (childleps->size() != 2) return m_default;

      return (*childleps)[1]->pt();
    }
  };

  //____________________________________________________________________________
  class deltaR_l1l2_h1 : public VarBase<float> {
  public:
  deltaR_l1l2_h1() : VarBase("deltaR_l1l2_h1") { m_default = -99; m_truthOnly = true; }
    ~deltaR_l1l2_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
      if (childleps->size() != 2) return m_default;

      return (*childleps)[0]->p4().DeltaR((*childleps)[1]->p4());
    }
  };

  //____________________________________________________________________________
  class ystar_pdg_flavor : public VarBase<int> {
  public:
  ystar_pdg_flavor() : VarBase("ystar_pdg_flavor") { m_default = -99; m_truthOnly = true; }
    ~ystar_pdg_flavor() { }

    int calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
      if (childleps->size() != 2) return m_default;

      return (*childleps)[0]->absPdgId();
    }
  };

  //____________________________________________________________________________
  class isNonHyyStarHiggs : public VarBase<char> {
  public:
  isNonHyyStarHiggs() : VarBase("isNonHyyStarHiggs") { m_default = false; m_truthOnly = true; }
    ~isNonHyyStarHiggs() { }

    // Note that this variable is intended for the GamGamStar cutflow, so background and
    // Leptonic dalitz should pass. Set by hand in CutflowAndMxAOD
  };

  //____________________________________________________________________________
  class pT_yDirect_h1 : public VarBase<float> {
  public:
  pT_yDirect_h1() : VarBase("pT_yDirect_h1") { m_default = -99; m_truthOnly = true; }
    ~pT_yDirect_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childphot = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsPhotons();
      if (childphot->size() != 1) return m_default;

      return (*childphot)[0]->pt();
    }
  };

  //____________________________________________________________________________
  class m_yStar_undressed_h1 : public VarBase<float> {
  public:
  m_yStar_undressed_h1() : VarBase("m_yStar_undressed_h1") { m_default = -99; m_truthOnly = true; }
    ~m_yStar_undressed_h1() { }

    // "undressed" here means no final-state radiation off of the leptons included.

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *childleps = HG::ExtraHggStarObjects::getInstance()->getTruthHiggsLeptons();
      if (childleps->size() != 2) return m_default;

      return ((*childleps)[0]->p4() + (*childleps)[1]->p4()).M();
    }
  };

  //____________________________________________________________________________
  class yyStarChannel : public VarBase<int> {
  public:
  yyStarChannel() : VarBase("yyStarChannel") { m_default = -99; }
    ~yyStarChannel() { }

    // Set by hand in CutflowAndMxAOD
    // Set the reco  value by specifying var::yyStarChannel.setValue(val)
    // Set the truth value by specifying var::yyStarChannel.setTruthValue(val)
  };
  
  //____________________________________________________________________________
  class yyStarChannelSimple : public VarBase<int> {
  public:
  yyStarChannelSimple() : VarBase("yyStarChannelSimple") { m_default = -99; m_truthOnly = true; }
    ~yyStarChannelSimple() { }

    // Set by hand in CutflowAndMxAOD
    // Set the truth value by specifying var::yyStarChannelSimple.setTruthValue(val)
  };

  //____________________________________________________________________________
  class ZyChannel : public VarBase<int> {
  public:
  ZyChannel() : VarBase("ZyChannel") { m_default = -99; }
    ~ZyChannel() { }

    // Set by hand in CutflowAndMxAOD
    // Set the reco  value by specifying var::ZyChannel.setValue(val)
    // Set the truth value by specifying var::ZyChannel.setTruthValue(val)
  };

  //____________________________________________________________________________
  class yyStarCategory : public VarBase<int> {
  public:
  yyStarCategory() : VarBase("yyStarCategory") { m_default = -99; }
    ~yyStarCategory() { }
    
    int calculateValue(bool truth); // See cxx file
  };
  
  //____________________________________________________________________________
  class yyStarCategory_electronOnly : public VarBase<int> {
  public:
  yyStarCategory_electronOnly() : VarBase("yyStarCategory_electronOnly") { m_default = -99; }
    ~yyStarCategory_electronOnly() { }

    int calculateValue(bool truth); // See cxx file
  };

  //____________________________________________________________________________
  class vertexTruthFitRadius : public VarBase<float> {
    public:
     vertexTruthFitRadius() : VarBase("vertexTruthFitRadius") { m_default = -99; }
     ~vertexTruthFitRadius() { }
  //
  };


  //____________________________________________________________________________
  class Dphi_lly_jj : public VarBase<float> {
  public:
    Dphi_lly_jj() : VarBase("Dphi_lly_jj") { m_default = -99; }
    ~Dphi_lly_jj() { }

    float calculateValue(bool truth)
    {
      // Note - for truth, this is calculated from the raw containers, and not
      // the true Higgs

      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);

      if (js->size() < 2)
      { return m_default; }

      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1) return m_default;

      if (mus->size() >= 2){
        return fabs(((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
      }
      if (eles->size() >= 2){
        return fabs(((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
      }
      if (!truth && eles->size() == 1){
        // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
        return fabs(((*gams)[0]->p4() + (*eles)[0]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
      }
      return m_default;
    }

  };

  //____________________________________________________________________________
  class Zepp_lly : public VarBase<float> {
  public:
    Zepp_lly() : VarBase("Zepp_lly") { m_default = -99; }
    ~Zepp_lly() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);

      if (js->size() < 2)
      { return m_default; }

      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1) return m_default;

      if (mus->size() >= 2){
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).Eta() - ((*js)[0]->eta() + (*js)[1]->eta()) / 2.0;
      }
      if (eles->size() >= 2){
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).Eta() - ((*js)[0]->eta() + (*js)[1]->eta()) / 2.0;
      }
      if (!truth && eles->size() == 1){
        // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
        return ((*eles)[0]->p4() + (*gams)[0]->p4()).Eta() - ((*js)[0]->eta() + (*js)[1]->eta()) / 2.0;
      }
      return m_default;
    }

  };

  //____________________________________________________________________________
  class pTt_lly : public VarBase<float> {
  public:
    pTt_lly() : VarBase("pTt_lly") { m_default = -99; }
    ~pTt_lly() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1) return m_default;

      TLorentzVector g1 = (*gams)[0]->p4();

      if (mus->size() >= 2){
        TLorentzVector g2 = (*mus)[0]->p4() + (*mus)[1]->p4();
        return fabs(g1.Px() * g2.Py() - g2.Px() * g1.Py()) / (g1 - g2).Pt() * 2.0;
      }
      if (eles->size() >= 2){
        TLorentzVector g2 = (*eles)[0]->p4() + (*eles)[1]->p4();
        return fabs(g1.Px() * g2.Py() - g2.Px() * g1.Py()) / (g1 - g2).Pt() * 2.0;
      }
      if (!truth && eles->size() == 1){
        // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
        TLorentzVector g2 = (*eles)[0]->p4();
        return fabs(g1.Px() * g2.Py() - g2.Px() * g1.Py()) / (g1 - g2).Pt() * 2.0;
      }
      return m_default;
    }

  };

  //____________________________________________________________________________
  class pT_llyjj : public VarBase<float> {
  public:
    pT_llyjj() : VarBase("pT_llyjj") { m_default = -99; }
    ~pT_llyjj() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *js = HG::VarHandler::getInstance()->getJets(truth);

      if (js->size() < 2)
      { return m_default; }

      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      if (gams->size() < 1) return m_default;

      if (mus->size() >= 2){
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).Pt();
      }
      if (eles->size() >= 2){
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).Pt();
      }
      if (!truth && eles->size() == 1){
        // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
        return ((*eles)[0]->p4() + (*gams)[0]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).Pt();
      }
      return m_default;
    }

  };

  //____________________________________________________________________________
  class DRmin_y_ystar_2jets : public VarBase<float> {
  public:
    DRmin_y_ystar_2jets() : VarBase("DRmin_y_ystar_2jets") { m_default = -99; }
    ~DRmin_y_ystar_2jets() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      if (gams->size() < 1) { return m_default; }
      if (jets->size() < 2) { return m_default; }

      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);

      double dR2min = 99.0, dR2 = 0.0;
      TLorentzVector gamStar;

      for (int j = 0; j < 2; j++) {
        auto jet = jets->at(j) ;

        // Check photon
        dR2 = xAOD::P4Helpers::deltaR2(*jet, (*gams)[0]->eta(), (*gams)[0]->phi(), false);
        if (dR2 < dR2min) { dR2min = dR2; }

        // Check gammaStar
        if (mus->size() >= 2){
          gamStar = (*mus)[0]->p4() + (*mus)[1]->p4();
          dR2 = xAOD::P4Helpers::deltaR2(*jet, gamStar.Eta(), gamStar.Phi(), false);
          if (dR2 < dR2min) { dR2min = dR2; }
        }
        if (eles->size() >= 2){
          gamStar = (*eles)[0]->p4() + (*eles)[1]->p4();
          dR2 = xAOD::P4Helpers::deltaR2(*jet, gamStar.Eta(), gamStar.Phi(), false);
          if (dR2 < dR2min) { dR2min = dR2; }
        }
        if (!truth && eles->size() == 1){
          // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
          gamStar = (*eles)[0]->p4();
          dR2 = xAOD::P4Helpers::deltaR2(*jet, gamStar.Eta(), gamStar.Phi(), false);
          if (dR2 < dR2min) { dR2min = dR2; }
        }
      }

      if (dR2min == 99) { return m_default; }

      return sqrt(dR2min);
    }
  };

  //____________________________________________________________________________
  class DRmin_y_leps_2jets : public VarBase<float> {
  public:
    DRmin_y_leps_2jets() : VarBase("DRmin_y_leps_2jets") { m_default = -99; }
    ~DRmin_y_leps_2jets() { }

    float calculateValue(bool truth)
    {
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);

      if (gams->size() < 1) { return m_default; }
      if (jets->size() < 2) { return m_default; }

      double dR2min = 99.0, dR2 = 0.0;

      for (int j = 0; j < 2; j++) {
        auto jet = jets->at(j) ;

        // Check photon
        dR2 = xAOD::P4Helpers::deltaR2(*jet, (*gams)[0]->eta(), (*gams)[0]->phi(), false);
        if (dR2 < dR2min) { dR2min = dR2; }

        // Check gammaStar
        if (mus->size() >= 2){
          for (int k = 0; k < 2; k++) {
            dR2 = xAOD::P4Helpers::deltaR2(*jet, (*mus)[k]->eta(), (*mus)[k]->phi(), false);
            if (dR2 < dR2min) { dR2min = dR2; }
          }
        }
        if (eles->size() >= 2){
          for (int k = 0; k < 2; k++) {
            dR2 = xAOD::P4Helpers::deltaR2(*jet, (*eles)[k]->eta(), (*eles)[k]->phi(), false);
            if (dR2 < dR2min) { dR2min = dR2; }
          }
        }
        if (!truth && eles->size() == 1){
          // If the electron container size is 1, it is merged (p4 set in SetMergedFourMomentum)
          dR2 = xAOD::P4Helpers::deltaR2(*jet, (*eles)[0]->eta(), (*eles)[0]->phi(), false);
          if (dR2 < dR2min) { dR2min = dR2; }
        }
      }

      if (dR2min == 99) { return m_default; }

      return sqrt(dR2min);
    }
  };

  //____________________________________________________________________________
  class m_lly2 : public VarBase<float> {
  public:
  m_lly2() : VarBase("m_lly2") { m_default = -99; }
    ~m_lly2() { }

    float calculateValue(bool truth)
    {
      // For Reco:
      // getElectrons and getMuons only return elecs / muons selected as the candidate y*.
      // getPhotons only returns the leading photon candidate.
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);

      // >= 2 muons and >= 2 photons
      if (mus->size() >= 2 && gams->size() >= 2)
	return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[1]->p4()).M();

      // >= 2 electrons and >= 2 photon2
      if (eles->size() >= 2 && gams->size() >= 2)
	return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[1]->p4()).M();

      return m_default;
    }
  };


  //____________________________________________________________________________

  void AssignZbosonIndices(const xAOD::IParticleContainer& leps,int& return_lep1i,int& return_lep2i,
                           double& return_mll,bool sortby_pt,
                           double closest_to_mev,float lead_pt_cut_gev=0); // Z = 91188

  bool eventIsNonHyyStarHiggs(const xAOD::TruthParticleContainer* allTruthParticles);
  bool eventIsHyyHiggs(const xAOD::TruthParticleContainer* allTruthParticles);
  bool isDirectlyFromHiggs(const xAOD::TruthParticle *ptcl);

  void DecorateLeptonDressing(const xAOD::IParticleContainer& leps, const xAOD::TruthParticleContainer& truthLeps);

}

namespace var {
  extern HG::m_lly m_lly;
  extern HG::m_lly_gev m_lly_gev;
  extern HG::m_ll m_ll;
  extern HG::m_l1y m_l1y;
  extern HG::m_l2y m_l2y;
  extern HG::deltaR_ll deltaR_ll;
  extern HG::Resolved_deltaPhiRescaled2 Resolved_deltaPhiRescaled2;
  extern HG::Resolved_deltaEta2 Resolved_deltaEta2;
  extern HG::Resolved_dRExtrapTrk12 Resolved_dRExtrapTrk12;
  extern HG::deltaPhi_trktrk_IP deltaPhi_trktrk_IP;
  extern HG::deltaEta_trktrk_IP deltaEta_trktrk_IP;
  extern HG::deltaR_track4mom deltaR_track4mom;
  extern HG::deltaPhi2_trktrk_perigee deltaPhi2_trktrk_perigee;
  extern HG::deltaEta2_trktrk_perigee deltaEta2_trktrk_perigee;
  extern HG::deltaRL2_trktrk_perigee deltaRL2_trktrk_perigee;
  extern HG::deltaPhi2_trktrk_LM deltaPhi2_trktrk_LM;
  extern HG::deltaEta2_trktrk_LM deltaEta2_trktrk_LM;
  extern HG::deltaRL2_trktrk_LM deltaRL2_trktrk_LM;
  extern HG::deltaPhi_naiveExtrap deltaPhi_naiveExtrap;
  extern HG::deltaPhi_overScaled_naiveExtrap deltaPhi_overScaled_naiveExtrap;
  extern HG::pt_lly pt_lly;
  extern HG::pt_ll pt_ll;
  extern HG::m_lly_track4mom m_lly_track4mom;
  extern HG::m_ll_track4mom m_ll_track4mom;
  extern HG::deltaPhi_ll_y deltaPhi_ll_y;
  extern HG::eta_y1 eta_y1;
  extern HG::pt_llyy pt_llyy;
  extern HG::m_llyy m_llyy;
  extern HG::m_emu m_emu;
  extern HG::m_emuy m_emuy;
  extern HG::Dy_j_j_50 Dy_j_j_50;
  extern HG::pT_l1 pT_l1;
  extern HG::eta_j1 eta_j1;
  extern HG::N_j_gap N_j_gap;
  extern HG::Zy_centrality Zy_centrality;
  extern HG::DR_Zy_jj DR_Zy_jj;
  extern HG::pT_l1_h1 pT_l1_h1;
  extern HG::pT_l2_h1 pT_l2_h1;
  extern HG::deltaR_l1l2_h1 deltaR_l1l2_h1;
  extern HG::ystar_pdg_flavor ystar_pdg_flavor;
  extern HG::isNonHyyStarHiggs isNonHyyStarHiggs;
  extern HG::pT_yDirect_h1 pT_yDirect_h1;
  extern HG::m_yStar_undressed_h1 m_yStar_undressed_h1;
  extern HG::yyStarChannel yyStarChannel;
  extern HG::yyStarChannelSimple yyStarChannelSimple;
  extern HG::ZyChannel ZyChannel;
  extern HG::vertexTruthFitRadius vertexTruthFitRadius;
  extern HG::trk_lead_pt trk_lead_pt;
  extern HG::yyStarCategory yyStarCategory;
  extern HG::yyStarCategory_electronOnly yyStarCategory_electronOnly;
  extern HG::Dphi_lly_jj Dphi_lly_jj;
  extern HG::Zepp_lly Zepp_lly;
  extern HG::pTt_lly pTt_lly;
  extern HG::pT_llyjj pT_llyjj;
  extern HG::DRmin_y_ystar_2jets DRmin_y_ystar_2jets;
  extern HG::DRmin_y_leps_2jets DRmin_y_leps_2jets;
  extern HG::m_lly2 m_lly2;
  extern HG::passVBFpresel passVBFpresel;
}


#endif // HGamGamStar_HggStarVariables_H
