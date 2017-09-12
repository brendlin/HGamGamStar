#ifndef HGamGamStar_HggStarVariables_H
#define HGamGamStar_HggStarVariables_H

#include "HGamAnalysisFramework/VarHandler.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

namespace HG {

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

}

namespace var {
  extern HG::m_lly m_lly;
  extern HG::m_ll m_ll;
}

#endif // HGamGamStar_HggStarVariables_H
