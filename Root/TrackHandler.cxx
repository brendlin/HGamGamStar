#include "HGamGamStar/TrackHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

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

  m_eleEtaCut  = config.getNum("ElectronHandler.Selection.MaxAbsEta", 2.47);
  m_elePtCut   = config.getNum("ElectronHandler.Selection.PtPreCutGeV", 25.0) * GeV;

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

  // Reserve this for ID variables.
  // Silicon, Pix hit cuts are applied in findTracksFromElectrons()

  for (auto track : container) {
    selected.push_back(track);
  }

  return selected;
}

//______________________________________________________________________________
CP::SystematicCode HG::TrackHandler::applySystematicVariation(const CP::SystematicSet &/*sys*/)
{
  return CP::SystematicCode::Ok;
}

//______________________________________________________________________________
xAOD::TrackParticleContainer HG::TrackHandler::findTracksFromElectrons(xAOD::TrackParticleContainer& container,const xAOD::ElectronContainer& elecs)
{

  xAOD::TrackParticleContainer selected(SG::VIEW_ELEMENTS);

  int nGoodQuality = 0;

  // std::cout << "~~~~~~" << std::endl;
  for (auto electron : elecs) {

    if ( std::abs( electron->caloCluster()->etaBE(2) ) > m_eleEtaCut ) continue;
    if ( electron->pt() < m_elePtCut ) continue;

    // std::cout << Form("Electron pt: %.0f eta: %.3f has %lu tracks",
    //                   electron->pt(),
    //                   electron->caloCluster()->etaBE(2),
    //                   electron->nTrackParticles()) << std::endl;

    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {

      const xAOD::TrackParticle* ele_tp = electron->trackParticle(i);

      int nSiHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfSiliconHitsAndDeadSensors(ele_tp);
      int nPixHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(ele_tp);
      // int passBLayerRequirement = ElectronSelectorHelpers::passBLayerRequirement(ele_tp);

      if ( std::abs(ele_tp->eta()) > m_etaCut ) continue;
      if ( ele_tp->pt() < m_ptCut ) continue;

      if (nSiHitsPlusDeadSensors  < 7) continue;
      if (nPixHitsPlusDeadSensors < 2) continue;
      // if (!passBLayerRequirement) continue;

      // std::cout << Form("Electron tp; pt: %.0f eta: %.3f",ele_tp->pt(),ele_tp->eta()) << std::endl;

      nGoodQuality++;

      bool found = false;

      // Check if it exists already
      for (auto track : selected) {
        if (ele_tp->p4() == track->p4())
        {
          // std::cout << "Found a duplicate track!" << std::endl;
          found = true;
          break;
        }
      }

      if (found) continue;

      for (auto track : container)
      {
        // std::cout << Form(" - GSF tp; pt: %.0f eta: %.3f ",track->pt(),track->eta())
        //           << (ele_tp->p4() == track->p4()) << std::endl;
        if (ele_tp->p4() == track->p4())
        {
          selected.push_back(track);
          found = true;
          break;
        }
      }

      if (!found) HG::fatal("Could not find TrackParticle associated to an electron!");
    }

  }

  // std::cout << "Number of elecss in Total: " << elecs.size() << std::endl;
  // std::cout << "Number of OQ trk in Total: " << nGoodQuality << std::endl;
  // std::cout << "Final size of container: " << container.size() << std::endl;
  // std::cout << "Final size of selected: " << selected.size() << std::endl;

  return selected;

}

xAOD::ElectronContainer
HG::TrackHandler::GetElecsAssociatedToTracks(const xAOD::TrackParticle& trk1,
                                             const xAOD::TrackParticle& trk2,
                                             xAOD::ElectronContainer& preSelElecs)
{

  xAOD::ElectronContainer selected(SG::VIEW_ELEMENTS);

  for (auto electron : preSelElecs) {

    if ( std::abs( electron->caloCluster()->etaBE(2) ) > m_eleEtaCut ) continue;
    if ( electron->pt() < m_elePtCut ) continue;

    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {

      const xAOD::TrackParticle* ele_tp = electron->trackParticle(i);

      if (ele_tp->p4() == trk1.p4() || ele_tp->p4() == trk2.p4())
      {
        selected.push_back(electron);
      }
    }
  }

  return selected;
}
