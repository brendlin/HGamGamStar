#include "HGamGamStar/HggStarVariables.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "HGamAnalysisFramework/TruthHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

#include "MCUtils/PIDUtils.h"

namespace var {
  HG::m_lly m_lly;
  HG::m_lly_gev m_lly_gev;
  HG::m_ll m_ll;
  HG::m_l1y m_l1y;
  HG::m_l2y m_l2y;
  HG::deltaR_ll deltaR_ll;
  HG::Resolved_dRExtrapTrk12 Resolved_dRExtrapTrk12;
  HG::Resolved_deltaPhiRescaled2 Resolved_deltaPhiRescaled2;
  HG::Resolved_deltaEta2 Resolved_deltaEta2;
  HG::pt_lly pt_lly;
  HG::pt_ll pt_ll;
  HG::m_lly_track4mom m_lly_track4mom;
  HG::m_ll_track4mom m_ll_track4mom;
  HG::deltaR_track4mom deltaR_track4mom;
  HG::deltaPhi_ll_y deltaPhi_ll_y;
  HG::eta_y1 eta_y1;
  HG::pt_llyy pt_llyy;
  HG::m_llyy m_llyy;
  HG::m_emu m_emu;
  HG::m_emuy m_emuy;
  HG::N_j_btag N_j_btag;
  HG::m_jj_50 m_jj_50;
  HG::Dy_j_j_50 Dy_j_j_50;
  HG::Zy_centrality Zy_centrality;
  HG::DR_Zy_jj DR_Zy_jj;
  HG::pT_l1_h1 pT_l1_h1;
  HG::pT_l2_h1 pT_l2_h1;
  HG::deltaR_l1l2_h1 deltaR_l1l2_h1;
  HG::ystar_pdg_flavor ystar_pdg_flavor;
  HG::isNonHyyStarHiggs isNonHyyStarHiggs;
  HG::pT_yDirect_h1 pT_yDirect_h1;
  HG::m_yStar_undressed_h1 m_yStar_undressed_h1;
  HG::yyStarChannel yyStarChannel;
  HG::ZyChannel ZyChannel;
  HG::vertexTruthFitRadius vertexTruthFitRadius;
  HG::trk_lead_pt trk_lead_pt;
  HG::yyStarCategory yyStarCategory;
  HG::Dphi_lly_jj Dphi_lly_jj;
  HG::Zepp_lly Zepp_lly;
  HG::pTt_lly pTt_lly;
  HG::pT_llyjj pT_llyjj;
  HG::DRmin_y_ystar_2jets DRmin_y_ystar_2jets;
  HG::DRmin_y_leps_2jets DRmin_y_leps_2jets;
}

// A special implementation of calculateValue that references another "var"
float HG::m_lly_gev::calculateValue(bool truth)
{
  if (truth) return var::m_lly.truth()/1000.;
  return var::m_lly()/1000.;
}

int HG::yyStarCategory::calculateValue(bool truth)
{
  if (truth) return CategoryEnum::CATEGORYUNKNOWN;
  bool passVBF = var::m_jj()/1000>400 && var::Deta_j_j() > 2.5;
  if (var::yyStarChannel()==ChannelEnum::DIMUON){
      if(passVBF) return CategoryEnum::VBF_DIMUON;
      else return CategoryEnum::GGF_DIMUON;
  }
  else if (var::yyStarChannel()==ChannelEnum::RESOLVED_DIELECTRON){
      if(passVBF) return CategoryEnum::VBF_RESOLVED_DIELECTRON;
      else return CategoryEnum::GGF_RESOLVED_DIELECTRON;
  }
  else if(var::yyStarChannel()==ChannelEnum::MERGED_DIELECTRON){
      if(passVBF) return CategoryEnum::VBF_MERGED_DIELECTRON;
      else return CategoryEnum::GGF_MERGED_DIELECTRON;
  }
  return CategoryEnum::CATEGORYUNKNOWN;
}

float HG::Resolved_dRExtrapTrk12::calculateValue(bool /* truth*/)
{
  if (var::yyStarChannel() != ChannelEnum::RESOLVED_DIELECTRON) return m_default;
  float deta_e1e2 = var::Resolved_deltaEta2();
  float dphi_e1e2 = var::Resolved_deltaPhiRescaled2();
  return sqrt(dphi_e1e2*dphi_e1e2 + deta_e1e2*deta_e1e2);
}

void HG::AssignZbosonIndices(const xAOD::IParticleContainer& leps,int& return_lep1i,int& return_lep2i,
                             double& return_mll,bool sortby_pt,
                             double closest_to_mev, float lead_pt_cut_gev){

  double min_delta = DBL_MAX;

  for (unsigned int i=0;i<leps.size();++i) {
    const xAOD::IParticle* lepi = leps[i];
    for (unsigned int j=0;j<leps.size();++j) {
      if (j == i) continue;

      const xAOD::IParticle* lepj = leps[j];

      if (lepi->pt() < lepj->pt()) continue;
      
      if(lepi->pt() < lead_pt_cut_gev*HG::GeV) continue;

      if (lepi->type() == xAOD::Type::TrackParticle &&
          ((xAOD::TrackParticle*)lepi)->charge() == ((xAOD::TrackParticle*)lepj)->charge()) continue;

      if (lepi->type() == xAOD::Type::Electron &&
          ((xAOD::Electron*)lepi)->charge() == ((xAOD::Electron*)lepj)->charge()) continue;

      if (lepi->type() == xAOD::Type::Muon &&
          ((xAOD::Muon*)lepi)->charge() == ((xAOD::Muon*)lepj)->charge()) continue;

      TLorentzVector tmp = lepi->p4() + lepj->p4();

      double metric = (sortby_pt ? tmp.Pt() : tmp.M() );

      if ( fabs( metric - closest_to_mev ) < min_delta) {
        min_delta = fabs( metric - closest_to_mev );
        return_lep1i = i;
        return_lep2i = j;
        return_mll = tmp.M();
      }

    }
  }
  return;
}

TLorentzVector HG::MergedEleTLV(const xAOD::TrackParticle& trk1, const xAOD::TrackParticle& trk2, const xAOD::Electron& ele)
{
  // New way: use the photon calibration, but the track TLVs.
  // The di-track "mass" is usually very close to the true one!
  float calibrated_e = HG::EleAcc::calibratedPhotonEnergy(ele);
  if (calibrated_e < 0) HG::fatal("Something went wrong - the energy of the photon is not calibrated correctly.");
  float scale_pt = calibrated_e/(trk1.e() + trk2.e());
  TLorentzVector tlv1;
  TLorentzVector tlv2;
  tlv1.SetPtEtaPhiM( trk1.pt() * scale_pt, trk1.eta(), trk1.phi(), ele.m() ); // ele.m == 0.510998
  tlv2.SetPtEtaPhiM( trk2.pt() * scale_pt, trk2.eta(), trk2.phi(), ele.m() );
  return (tlv1 + tlv2);

  // Old way: electron pt
  //
  // float scale_pt = ele.pt()/(trk1.pt() + trk2.pt());
  // TLorentzVector tlv1;
  // TLorentzVector tlv2;
  // tlv1.SetPtEtaPhiM( trk1.pt() * scale_pt, trk1.eta(), trk1.phi(), ele.m() ); // ele.m == 0.510998
  // tlv2.SetPtEtaPhiM( trk2.pt() * scale_pt, trk2.eta(), trk2.phi(), ele.m() );
  // return (tlv1 + tlv2);
}

HG::TruthPtcls HG::getHyyStarSignalDecayProducts(const xAOD::TruthParticle *ptcl)
{
  // Recursive, starting from the Higgs
  // STABLE particles returned.

  if (ptcl == nullptr) { HG::fatal("getHyyStarSignalDecayProducts FATAL: particle is NULL"); }

  TruthPtcls decay(SG::VIEW_ELEMENTS);

  // skip anything that is not an electron, muon or photon
  if (not MCUtils::PID::isPhoton( ptcl->pdgId() ) &&
      not MCUtils::PID::isElectron( ptcl->pdgId() ) &&
      not MCUtils::PID::isMuon( ptcl->pdgId() ) &&
      not MCUtils::PID::isHiggs( ptcl->pdgId() )) return decay;

  if (HG::isStable(ptcl))
  {
    decay.push_back(ptcl);
    return decay;
  }

  for (size_t ichild = 0; ichild < ptcl->nChildren(); ++ichild)
  {
   if (ptcl->child(ichild))
    {
      for (auto p : getHyyStarSignalDecayProducts(ptcl->child(ichild)))
      { decay.push_back(p); }
    }
  }

  return decay;
}

HG::TruthPtcls HG::FilterLeptons(const TruthPtcls& stableHiggsDecayProducts) {
  // Assuming leptons from taus or hadron decays already excluded

  TruthPtcls childleps(SG::VIEW_ELEMENTS);

  for (auto child : stableHiggsDecayProducts) {
    if (MCUtils::PID::isElectron(child->pdgId()) || MCUtils::PID::isMuon(child->pdgId()))
    {
      childleps.push_back(child);
    }
  }

  childleps.sort(TruthHandler::comparePt);

  return childleps;
}

HG::TruthPtcls HG::FilterDirectPhotons(const TruthPtcls& stableHiggsDecayProducts) {
  // Assuming leptons from taus or hadron decays already excluded

  TruthPtcls directphots(SG::VIEW_ELEMENTS);

  for (auto child : stableHiggsDecayProducts) {
    if (MCUtils::PID::isPhoton(child->pdgId()) && isDirectlyFromHiggs(child) )
    {
      directphots.push_back(child);
    }
  }

  return directphots;
}

bool HG::isDirectlyFromHiggs(const xAOD::TruthParticle *ptcl)
{
  // (Recursive)

  if (ptcl == nullptr) { HG::fatal("isFromHiggs FATAL: particle is NULL"); }

  if (MCUtils::PID::isHiggs(ptcl->pdgId())) { return true; }

  if (ptcl->parent() == nullptr) { return false; }

  if ( MCUtils::PID::isHiggs( ptcl->parent()->pdgId() ) ) { return true; }
  if (ptcl->parent()->pdgId() != ptcl->pdgId()) { return false; }

  return isFromHiggs(ptcl->parent());
}


bool HG::eventIsNonHyyStarHiggs(const xAOD::TruthParticleContainer * allParticles) {

  TruthPtcls higgses = getFinalHiggsBosons(allParticles);

  if (higgses.size() == 0) return false; // E.g. Background events

  TruthPtcls decayProds = getHyyStarSignalDecayProducts(higgses[0]);

  TruthPtcls childleps = HG::FilterLeptons(decayProds);
  if (childleps.size() != 2) return true;

  TruthPtcls directphots = HG::FilterDirectPhotons(decayProds);
  if (directphots.size() != 1) return true;

  return false;
}

void HG::DecorateLeptonDressing(const xAOD::IParticleContainer& leps, const xAOD::TruthParticleContainer& truthLeps){

  //leps from truthHandler(), truthLeps from TruthMuons or TruthElectrons Container which contains dressed quantities
  static SG::AuxElement::Accessor<float> pt_dressed("pt_dressed");
  static SG::AuxElement::Accessor<float> eta_dressed("eta_dressed");
  static SG::AuxElement::Accessor<float> phi_dressed("phi_dressed");
  static SG::AuxElement::Accessor<float> e_dressed("e_dressed");
  static SG::AuxElement::Accessor<int> nPhotons_dressed("nPhotons_dressed");

  if (leps.size() >0){
    for (auto lep: leps){
      for(auto truthLep: truthLeps){
        if (lep->auxdata<int>("barcode") == truthLep->barcode()) {
          pt_dressed(*lep)=pt_dressed(*truthLep);
          eta_dressed(*lep)=eta_dressed(*truthLep);
          phi_dressed(*lep)=phi_dressed(*truthLep);
          e_dressed(*lep)=e_dressed(*truthLep);
          nPhotons_dressed(*lep)=nPhotons_dressed(*truthLep);
          break;
        }
      }
    }
  }
  return;
}
