#ifndef HGamGamStar_HggStarVariables_H
#define HGamGamStar_HggStarVariables_H

#include "HGamAnalysisFramework/VarHandler.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

#include "xAODBase/IParticleContainer.h"
#include "xAODBase/IParticle.h"

namespace HG {

  TruthPtcls getHyyStarSignalDecayProducts(const xAOD::TruthParticle *ptcl);
  TruthPtcls FilterLeptons(TruthPtcls stableHiggsDecayProducts);
  TruthPtcls FilterDirectPhotons(TruthPtcls stableHiggsDecayProducts);

  //____________________________________________________________________________
  class m_lly : public VarBase<float> {
  public:
  m_lly() : VarBase("m_lly") { m_default = -99; }
    ~m_lly() { }

    float calculateValue(bool truth)
    {
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).M();
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).M();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_ll : public VarBase<float> {
  public:
  m_ll() : VarBase("m_ll") { m_default = -99; }
    ~m_ll() { }

    float calculateValue(bool truth)
    {
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4()).M();
      if (eles->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4()).M();
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
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      const xAOD::IParticleContainer *gams = HG::VarHandler::getInstance()->getPhotons(truth);
      if (mus->size() >= 2 && gams->size() >= 1)
        return ((*mus)[0]->p4() + (*mus)[1]->p4() + (*gams)[0]->p4()).Pt();
      if (eles->size() >= 2 && gams->size() >= 1)
        return ((*eles)[0]->p4() + (*eles)[1]->p4() + (*gams)[0]->p4()).Pt();
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
      (void)truth;
      const xAOD::IParticleContainer *eles = HG::VarHandler::getInstance()->getElectrons(truth);
      const xAOD::IParticleContainer *mus = HG::VarHandler::getInstance()->getMuons(truth);
      if (mus->size() >= 2)
        return ((*mus)[0]->p4() + (*mus)[1]->p4()).Pt();
      if (eles->size() >= 2)
        return ((*eles)[0]->p4() + (*eles)[1]->p4()).Pt();
      return m_default;
    }
  };

  //____________________________________________________________________________
  class m_llyy : public VarBase<float> {
  public:
  m_llyy() : VarBase("m_llyy") { m_default = -99; }
    ~m_llyy() { }

    float calculateValue(bool truth)
    {
      (void)truth;
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
  class pT_l1_h1 : public VarBase<float> {
  public:
  pT_l1_h1() : VarBase("pT_l1_h1") { m_default = -99; m_truthOnly = true; }
    ~pT_l1_h1() { }

    float calculateValue(bool truth)
    {
      if (not truth)
      { return m_default; }

      const xAOD::TruthParticleContainer *higgses = (xAOD::TruthParticleContainer*)HG::VarHandler::getInstance()->getHiggsBosons();

      if (higgses->size() == 0) return m_default;

      TruthPtcls decayProds = getHyyStarSignalDecayProducts((*higgses)[0]);
      TruthPtcls childleps = FilterLeptons(decayProds);

      if (childleps.size() != 2) return m_default;
      if (childleps[0]->absPdgId() != childleps[1]->absPdgId()) return m_default;

      return childleps[0]->pt();
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

      const xAOD::TruthParticleContainer *higgses = (xAOD::TruthParticleContainer*)HG::VarHandler::getInstance()->getHiggsBosons();

      if (higgses->size() == 0) return m_default;

      TruthPtcls decayProds = getHyyStarSignalDecayProducts((*higgses)[0]);
      TruthPtcls childleps = FilterLeptons(decayProds);

      if (childleps.size() != 2) return m_default;
      if (childleps[0]->absPdgId() != childleps[1]->absPdgId()) return m_default;

      return childleps[1]->pt();
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

      const xAOD::TruthParticleContainer *higgses = (xAOD::TruthParticleContainer*)HG::VarHandler::getInstance()->getHiggsBosons();

      if (higgses->size() == 0) return m_default;

      TruthPtcls decayProds = getHyyStarSignalDecayProducts((*higgses)[0]);
      TruthPtcls childleps = FilterLeptons(decayProds);

      if (childleps.size() != 2) return m_default;
      if (childleps[0]->absPdgId() != childleps[1]->absPdgId()) return m_default;

      return childleps[0]->p4().DeltaR(childleps[1]->p4());
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

      const xAOD::TruthParticleContainer *higgses = (xAOD::TruthParticleContainer*)HG::VarHandler::getInstance()->getHiggsBosons();

      if (higgses->size() == 0) return m_default;

      TruthPtcls decayProds = getHyyStarSignalDecayProducts((*higgses)[0]);
      TruthPtcls childleps = FilterLeptons(decayProds);

      if (childleps.size() != 2) return m_default;
      if (childleps[0]->absPdgId() != childleps[1]->absPdgId()) return m_default;

      return childleps[0]->absPdgId();
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

      const xAOD::TruthParticleContainer *higgses = (xAOD::TruthParticleContainer*)HG::VarHandler::getInstance()->getHiggsBosons();

      if (higgses->size() == 0) return m_default;

      TruthPtcls decayProds = getHyyStarSignalDecayProducts((*higgses)[0]);
      TruthPtcls childphot = FilterDirectPhotons(decayProds);

      if (childphot.size() != 1) return m_default;
      return childphot[0]->pt();
    }
  };

  //____________________________________________________________________________

  void AssignZbosonIndices(const xAOD::IParticleContainer& leps,int& return_lep1i,int& return_lep2i,
                           double& return_mll,double closest_to=91188.);

  bool eventIsNonHyyStarHiggs(const xAOD::TruthParticleContainer* allTruthParticles);
  bool isDirectlyFromHiggs(const xAOD::TruthParticle *ptcl);

}

namespace var {
  extern HG::m_lly m_lly;
  extern HG::m_ll m_ll;
  extern HG::pt_lly pt_lly;
  extern HG::pt_ll pt_ll;
  extern HG::pt_llyy pt_llyy;
  extern HG::m_llyy m_llyy;
  extern HG::pT_l1_h1 pT_l1_h1;
  extern HG::pT_l2_h1 pT_l2_h1;
  extern HG::deltaR_l1l2_h1 deltaR_l1l2_h1;
  extern HG::ystar_pdg_flavor ystar_pdg_flavor;
  extern HG::isNonHyyStarHiggs isNonHyyStarHiggs;
  extern HG::pT_yDirect_h1 pT_yDirect_h1;
}


#endif // HGamGamStar_HggStarVariables_H
