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
  TLorentzVector MergedEleTLV(const xAOD::TrackParticle& trk1, const xAOD::TrackParticle& trk2, const xAOD::Electron& ele);

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

      // If the electron container size is 1, then take the (cluster) e-gamma mass (yystar)
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail() && gams->size() >= 1)
        return ((*gams)[0]->p4() + *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV()).M();

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

      // For merged electrons, take the TLV set in ExtraHggStarObjects
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail())
        return HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV()->M();

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
  class Resolved_dRExtrapTrk12 : public VarBase<float> {
  public:
  Resolved_dRExtrapTrk12() : VarBase("Resolved_dRExtrapTrk12") { m_default = -99; m_recoOnly = true; }
    ~Resolved_dRExtrapTrk12() { }

    float calculateValue(bool truth); // See cxx file
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
        int flipSign1 = ( (*eles)[0]->trackParticle()->charge() > 0) ? 1 : -1; // -1 = flip
        int flipSign2 = ( (*eles)[1]->trackParticle()->charge() > 0) ? 1 : -1; // -1 = flip
        float dphi_trk1,dphi_trk2;
        // Christos says that track-matching uses DeltaPhiRescaled, since this matches barycenter better
        // and that Layer 2 is best.
        (*eles)[0]->trackCaloMatchValue(dphi_trk1, xAOD::EgammaParameters::deltaPhiRescaled2);
        (*eles)[1]->trackCaloMatchValue(dphi_trk2, xAOD::EgammaParameters::deltaPhiRescaled2);
        float dphi_e1e2 = xAOD::P4Helpers::deltaPhi((*eles)[0]->caloCluster()->phiBE(2) + dphi_trk1*flipSign1,
                                                    (*eles)[1]->caloCluster()->phiBE(2) + dphi_trk2*flipSign2);
        return dphi_e1e2;
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
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail() && gams->size() >= 1)
        return ((*gams)[0]->p4() + *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV()).Pt();

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
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail())
        return HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV()->Pt();
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

      return ((*trks)[0]->p4() + (*trks)[1]->p4() + (*gams)[0]->p4()).M();
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

      return ((*trks)[0]->p4() + (*trks)[1]->p4()).M();
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
        return fabs(((*mus)[0]->p4() + (*mus)[1]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
      }
      if (eles->size() >= 2){
        return fabs(((*eles)[0]->p4() + (*eles)[1]->p4()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
      }
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
        return fabs(((*gams)[0]->p4() + *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV()).DeltaPhi((*js)[0]->p4() + (*js)[1]->p4()));
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
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
        return (*HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV() + (*gams)[0]->p4()).Eta() - ((*js)[0]->eta() + (*js)[1]->eta()) / 2.0;
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
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
        TLorentzVector g2 = *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV();
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
      if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
        return (*HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV() + (*gams)[0]->p4() + (*js)[0]->p4() + (*js)[1]->p4()).Pt();
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
        if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
          gamStar = *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV();
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
        if (!truth && HG::ExtraHggStarObjects::getInstance()->mergedElectronTLVAvail()){
          auto gamStar = *HG::ExtraHggStarObjects::getInstance()->getMergedElectronTLV();
          dR2 = xAOD::P4Helpers::deltaR2(*jet, gamStar.Eta(), gamStar.Phi(), false);
          if (dR2 < dR2min) { dR2min = dR2; }
        }
      }

      if (dR2min == 99) { return m_default; }

      return sqrt(dR2min);
    }
  };

  //____________________________________________________________________________

  void AssignZbosonIndices(const xAOD::IParticleContainer& leps,int& return_lep1i,int& return_lep2i,
                           double& return_mll,bool sortby_pt,double closest_to); // Z = 91188

  bool eventIsNonHyyStarHiggs(const xAOD::TruthParticleContainer* allTruthParticles);
  bool isDirectlyFromHiggs(const xAOD::TruthParticle *ptcl);

  void DecorateLeptonDressing(const xAOD::IParticleContainer& leps, const xAOD::TruthParticleContainer& truthLeps);

}

namespace var {
  extern HG::m_lly m_lly;
  extern HG::m_lly_gev m_lly_gev;
  extern HG::m_ll m_ll;
  extern HG::deltaR_ll deltaR_ll;
  extern HG::Resolved_dRExtrapTrk12 Resolved_dRExtrapTrk12;
  extern HG::Resolved_deltaPhiRescaled2 Resolved_deltaPhiRescaled2;
  extern HG::Resolved_deltaEta2 Resolved_deltaEta2;
  extern HG::pt_lly pt_lly;
  extern HG::pt_ll pt_ll;
  extern HG::m_lly_track4mom m_lly_track4mom;
  extern HG::m_ll_track4mom m_ll_track4mom;
  extern HG::deltaPhi_ll_y deltaPhi_ll_y;
  extern HG::eta_y1 eta_y1;
  extern HG::pt_llyy pt_llyy;
  extern HG::m_llyy m_llyy;
  extern HG::pT_l1_h1 pT_l1_h1;
  extern HG::pT_l2_h1 pT_l2_h1;
  extern HG::deltaR_l1l2_h1 deltaR_l1l2_h1;
  extern HG::ystar_pdg_flavor ystar_pdg_flavor;
  extern HG::isNonHyyStarHiggs isNonHyyStarHiggs;
  extern HG::pT_yDirect_h1 pT_yDirect_h1;
  extern HG::m_yStar_undressed_h1 m_yStar_undressed_h1;
  extern HG::yyStarChannel yyStarChannel;
  extern HG::Dphi_lly_jj Dphi_lly_jj;
  extern HG::Zepp_lly Zepp_lly;
  extern HG::pTt_lly pTt_lly;
  extern HG::pT_llyjj pT_llyjj;
  extern HG::DRmin_y_ystar_2jets DRmin_y_ystar_2jets;
  extern HG::DRmin_y_leps_2jets DRmin_y_leps_2jets;
}


#endif // HGamGamStar_HggStarVariables_H
