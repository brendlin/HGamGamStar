#include "HGamGamStar/HggStarVariables.h"

#include "HGamAnalysisFramework/TruthHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

namespace var {
  HG::m_lly m_lly;
  HG::m_ll m_ll;
  HG::pt_lly pt_lly;
  HG::pt_ll pt_ll;
  HG::m_lly_track4mom m_lly_track4mom;
  HG::m_ll_track4mom m_ll_track4mom;
  HG::deltaPhi_ll_y deltaPhi_ll_y;
  HG::eta_y1 eta_y1;
  HG::pt_llyy pt_llyy;
  HG::m_llyy m_llyy;
  HG::pT_l1_h1 pT_l1_h1;
  HG::pT_l2_h1 pT_l2_h1;
  HG::deltaR_l1l2_h1 deltaR_l1l2_h1;
  HG::ystar_pdg_flavor ystar_pdg_flavor;
  HG::isNonHyyStarHiggs isNonHyyStarHiggs;
  HG::pT_yDirect_h1 pT_yDirect_h1;
  HG::m_yStar_undressed_h1 m_yStar_undressed_h1;
}

void HG::AssignZbosonIndices(const xAOD::IParticleContainer& leps,int& return_lep1i,int& return_lep2i,
                             double& return_mll,double closest_to){

  double min_delta = DBL_MAX;
  bool sortby_pt = true;

  for (unsigned int i=0;i<leps.size();++i) {
    const xAOD::IParticle* lepi = leps[i];
    for (unsigned int j=0;j<leps.size();++j) {
      if (j == i) continue;

      const xAOD::IParticle* lepj = leps[j];

      if (lepi->pt() < lepj->pt()) continue;

      if (lepi->type() == xAOD::Type::TrackParticle &&
          ((xAOD::TrackParticle*)lepi)->charge() == ((xAOD::TrackParticle*)lepj)->charge()) continue;

      if (lepi->type() == xAOD::Type::Electron &&
          ((xAOD::Electron*)lepi)->charge() == ((xAOD::Electron*)lepj)->charge()) continue;

      if (lepi->type() == xAOD::Type::Muon &&
          ((xAOD::Muon*)lepi)->charge() == ((xAOD::Muon*)lepj)->charge()) continue;

      TLorentzVector tmp = lepi->p4() + lepj->p4();

      double metric = (sortby_pt ? tmp.Pt() : tmp.M() );

      if ( fabs( metric - closest_to ) < min_delta) {
        min_delta = fabs( metric - closest_to );
        return_lep1i = i;
        return_lep2i = j;
        return_mll = tmp.M();
      }

    }
  }
  return;
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

HG::TruthPtcls HG::FilterLeptons(TruthPtcls stableHiggsDecayProducts) {
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

HG::TruthPtcls HG::FilterDirectPhotons(TruthPtcls stableHiggsDecayProducts) {
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
        }
      }
    }
  }
  return;
}
