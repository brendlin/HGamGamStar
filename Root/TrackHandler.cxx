#include "HGamGamStar/TrackHandler.h"

//______________________________________________________________________________
HG::TrackHandler::TrackHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
{

}

//______________________________________________________________________________
HG::TrackHandler::~TrackHandler()
{

}

//______________________________________________________________________________
EL::StatusCode HG::TrackHandler::initialize(Config &config)
{
  HgammaHandler::initialize(config);

  // Read in configuration information
  m_containerName = config.getStr(m_name + ".ContainerName", "GSFTrackParticles");

  //electron selection
  m_doTrqCuts  = config.getBool(m_name + ".Selection.ApplyTRQCuts", true);
  m_nSiMin     = config.getInt (m_name + ".Selection.nSiMin",7);
  m_nPixMin    = config.getInt (m_name + ".Selection.nPixMin",2);

  m_etaCut     = config.getNum(m_name + ".Selection.MaxAbsEta", 2.47);
  m_ptCut      = config.getNum(m_name + ".Selection.PtPreCutGeV", 0.3) * GeV;

  return EL::StatusCode::SUCCESS;
}

//______________________________________________________________________________
xAOD::TrackParticleContainer HG::TrackHandler::getCorrectedContainer()
{
  bool calib = false;
  xAOD::TrackParticleContainer shallowContainer = getShallowContainer(calib);

  // sort the tracks
  shallowContainer.sort(comparePt);

  return shallowContainer;
}

//______________________________________________________________________________
xAOD::TrackParticleContainer HG::TrackHandler::applySelection(xAOD::TrackParticleContainer &container)
{
  xAOD::TrackParticleContainer selected(SG::VIEW_ELEMENTS);

  for (auto track : container) {
    selected.push_back(track);
  }

  return selected;
}

//______________________________________________________________________________
CP::SystematicCode HG::TrackHandler::applySystematicVariation(const CP::SystematicSet &sys)
{
  return CP::SystematicCode::Ok;
}
