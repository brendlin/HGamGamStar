#include "HGamGamStar/HggStarVariables.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "HGamAnalysisFramework/TruthHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

#include "MCUtils/PIDUtils.h"

namespace var {
  HG::m_lly m_lly;
  HG::m_lly_gev m_lly_gev;
  HG::m_llyj m_llyj;
  HG::m_ll m_ll;
  HG::m_l1y m_l1y;
  HG::m_l2y m_l2y;
  HG::deltaR_ll deltaR_ll;
  HG::Resolved_deltaPhiRescaled2 Resolved_deltaPhiRescaled2;
  HG::Resolved_deltaEta2 Resolved_deltaEta2;
  HG::Resolved_dRExtrapTrk12 Resolved_dRExtrapTrk12;
  HG::deltaPhi_trktrk_IP deltaPhi_trktrk_IP;
  HG::deltaEta_trktrk_IP deltaEta_trktrk_IP;
  HG::deltaR_track4mom deltaR_track4mom;
  HG::deltaPhi2_trktrk_perigee deltaPhi2_trktrk_perigee;
  HG::deltaEta2_trktrk_perigee deltaEta2_trktrk_perigee;
  HG::deltaRL2_trktrk_perigee deltaRL2_trktrk_perigee;
  HG::deltaPhi2_trktrk_LM deltaPhi2_trktrk_LM;
  HG::deltaEta2_trktrk_LM deltaEta2_trktrk_LM;
  HG::deltaRL2_trktrk_LM deltaRL2_trktrk_LM;
  HG::deltaPhi_naiveExtrap deltaPhi_naiveExtrap;
  HG::deltaPhi_overScaled_naiveExtrap deltaPhi_overScaled_naiveExtrap;
  HG::pt_lly pt_lly;
  HG::pt_ll pt_ll;
  HG::m_lly_track4mom m_lly_track4mom;
  HG::m_ll_track4mom m_ll_track4mom;
  HG::deltaPhi_ll_y deltaPhi_ll_y;
  HG::eta_y1 eta_y1;
  HG::pt_llyy pt_llyy;
  HG::m_llyy m_llyy;
  HG::m_emu m_emu;
  HG::m_emuy m_emuy;
  HG::Dy_j_j_50 Dy_j_j_50;
  HG::pT_l1 pT_l1;
  HG::eta_j1 eta_j1;
  HG::N_j_gap N_j_gap;
  HG::N_j_Mix N_j_Mix;
  HG::N_j_Mix30 N_j_Mix30;
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
  HG::yyStarChannelSimple yyStarChannelSimple;
  HG::ZyChannel ZyChannel;
  HG::vertexTruthFitRadius vertexTruthFitRadius;
  HG::trk_lead_pt trk_lead_pt;
  HG::yyStarCategory yyStarCategory;
  HG::yyStarCategory_electronOnly yyStarCategory_electronOnly;
  HG::Dphi_lly_jj Dphi_lly_jj;
  HG::Dphi_yj Dphi_yj;
  HG::Zepp_lly Zepp_lly;
  HG::pTt_lly pTt_lly;
  HG::pT_llyjj pT_llyjj;
  HG::pt_llyj pt_llyj;
  HG::DRmin_y_ystar_2jets DRmin_y_ystar_2jets;
  HG::DRmin_y_leps_2jets DRmin_y_leps_2jets;
  HG::m_lly2 m_lly2;
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
  const xAOD::IParticleContainer *jets = HG::VarHandler::getInstance()->getJets(truth);
  bool passVBF = false;
  if (jets->size() > 1){ //we have two or more jets - prerequisite to pass VBF
      passVBF =  
      (*jets)[0]->pt() > 25 * HG::GeV   &&
      (*jets)[1]->pt() > 25 * HG::GeV   &&
      fabs(var::Zepp_lly()) < 2.0       &&
      var::DRmin_y_leps_2jets() > 1.5   &&
      var::Dphi_lly_jj() > 2.8          &&
      var::m_jj() > 500.0 * HG::GeV     &&
      var::Deta_j_j() > 2.7;
      //if all cuts pass finally check if we have forward jets and cut on min pT
      if (passVBF){
        float fw_jet0_pt = fabs((*jets)[0]->eta()) > 2.5 ? (*jets)[0]->pt() : 999 * HG::GeV;
        float fw_jet1_pt = fabs((*jets)[1]->eta()) > 2.5 ? (*jets)[1]->pt() : 999 * HG::GeV;
        passVBF = (fw_jet0_pt > 30 * HG::GeV && fw_jet1_pt > 30 * HG::GeV);
      }
  }
  bool passPtT = (var::pTt_lly() > 100.0 * HG::GeV);

  //prioritize first VBF, then high pTt, then fall back to Inclusive
  if (var::yyStarChannel()==ChannelEnum::DIMUON){
      if(passVBF) return CategoryEnum::VBF_DIMUON;
      if(passPtT) return CategoryEnum::HIPTT_DIMUON;
      return CategoryEnum::GGF_DIMUON;
  }
  else if (var::yyStarChannel()==ChannelEnum::RESOLVED_DIELECTRON){
      if(passVBF) return CategoryEnum::VBF_RESOLVED_DIELECTRON;
      if(passPtT) return CategoryEnum::HIPTT_RESOLVED_DIELECTRON;
      return CategoryEnum::GGF_RESOLVED_DIELECTRON;
  }
  else if(var::yyStarChannel()==ChannelEnum::MERGED_DIELECTRON){
      if(passVBF) return CategoryEnum::VBF_MERGED_DIELECTRON;
      if(passPtT) return CategoryEnum::HIPTT_MERGED_DIELECTRON;
      return CategoryEnum::GGF_MERGED_DIELECTRON;
  }
  return CategoryEnum::CATEGORYUNKNOWN;
}

int HG::yyStarCategory_electronOnly::calculateValue(bool truth)
{
  int chan = var::yyStarChannel();
  int cat = var::yyStarCategory();

  // If it is a muon channel, set to 0 (we do not want it here)
  if (chan == ChannelEnum::DIMUON) return CategoryEnum::CATEGORYUNKNOWN;

  // Assuming the muon category is always 1, 4, 7, ...
  // here is how to translate from 2, 3, 5, 6, 8, 9, ... to 1, 2, 3, 4, 5, 6:
  cat = cat - int( (cat+1)/3 );

  return cat;
}

float HG::Resolved_dRExtrapTrk12::calculateValue(bool /* truth*/)
{
  if (var::yyStarChannel() != ChannelEnum::RESOLVED_DIELECTRON) return m_default;
  float deta_e1e2 = var::Resolved_deltaEta2();
  float dphi_e1e2 = var::Resolved_deltaPhiRescaled2();
  return sqrt(dphi_e1e2*dphi_e1e2 + deta_e1e2*deta_e1e2);
}

float HG::deltaR_track4mom::calculateValue(bool /* truth */)
{
  float deta_e1e2 = var::deltaEta_trktrk_IP();
  float dphi_e1e2 = var::deltaPhi_trktrk_IP();
  return sqrt(dphi_e1e2*dphi_e1e2 + deta_e1e2*deta_e1e2);
}

float HG::deltaRL2_trktrk_perigee::calculateValue(bool /* truth */)
{
  float deta_e1e2 = var::deltaEta2_trktrk_perigee();
  float dphi_e1e2 = var::deltaPhi2_trktrk_perigee();
  return sqrt(dphi_e1e2*dphi_e1e2 + deta_e1e2*deta_e1e2);
}

float HG::deltaRL2_trktrk_LM::calculateValue(bool /* truth */)
{
  float deta_e1e2 = var::deltaEta2_trktrk_LM();
  float dphi_e1e2 = var::deltaPhi2_trktrk_LM();
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

void HG::SetMergedFourMomentum(xAOD::Electron& ele,const float calibrated_e) {
  // Using the example from:
  // acode-browser.usatlas.bnl.gov/lxr/source/athena/Reconstruction/egamma/egammaTools/src/EMFourMomBuilder.cxx?v=21.2
  // calibrated_e should come from calibration of the merged electron as a converted photon.

  const float E = calibrated_e;
  const float eta = HG::EleAcc::vtxEta(ele);
  const float phi = HG::EleAcc::vtxPhi(ele);
  const float mass = HG::EleAcc::vtxM(ele);
  const double pt = (E > mass) ? sqrt(E*E - mass*mass)/cosh(eta) : 0;

  ele.setP4(pt, eta, phi, mass);

  // std::cout << "New pt: " << ele.pt() << " eta: " << ele.eta()
  //           << " phi: " << ele.phi() << " Mass: " << ele.m() << std::endl;

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

bool HG::eventIsHyyHiggs(const xAOD::TruthParticleContainer * allParticles) {

  TruthPtcls higgses = getFinalHiggsBosons(allParticles);

  if (higgses.size() == 0) return false;

  TruthPtcls decayProds = getHyyStarSignalDecayProducts(higgses[0]);

  TruthPtcls childleps = HG::FilterLeptons(decayProds);
  if (childleps.size() != 0) return false;

  TruthPtcls directphots = HG::FilterDirectPhotons(decayProds);
  if (directphots.size() != 2) return false;

  return true;
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
