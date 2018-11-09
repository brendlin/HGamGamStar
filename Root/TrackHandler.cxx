#include "HGamGamStar/TrackHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

SG::AuxElement::Accessor< std::vector<int> > HG::TrackHandler::MatchedElectrons("MatchedElectrons");
SG::AuxElement::Accessor<char>  HG::TrackHandler::passIPCut("passIPCut");
SG::AuxElement::Accessor<float>  HG::TrackHandler::d0significance("d0significance");
SG::AuxElement::Accessor<float>  HG::TrackHandler::z0sinTheta("z0sinTheta");

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
  
  m_d0BySigd0Cut = config.getNum("ElectronHandler.Selection.d0BySigd0Max", 5.0);
  m_z0Cut = config.getNum("ElectronHandler.Selection.z0Max", 0.5);

  return EL::StatusCode::SUCCESS;
}

//______________________________________________________________________________
xAOD::TrackParticleContainer HG::TrackHandler::getCorrectedContainer()
{
  bool calib = false;
  xAOD::TrackParticleContainer shallowContainer = getShallowContainer(calib);

  // sort the tracks
  shallowContainer.sort(comparePt);
  
  for (auto trk : shallowContainer){
    decorateIPCut(*trk);
  }

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
          // Push the track into the "selected" container
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

//______________________________________________________________________________
xAOD::ElectronContainer
HG::TrackHandler::GetElecsAssociatedToTracks(xAOD::TrackParticle& trk1,
                                             xAOD::TrackParticle& trk2,
                                             xAOD::ElectronContainer& preSelElecs)
{
  // The electron container outputted here should be
  // the _final selected electrons_ that you will write out, otherwise
  // the indices linking the electrons and the tracks will not be correct.

  xAOD::ElectronContainer selected(SG::VIEW_ELEMENTS);

  MatchedElectrons(trk1).clear();
  MatchedElectrons(trk2).clear();

  unsigned int index_selected = 0;

  for (auto electron : preSelElecs) {

    bool matches = false;

    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {
      const xAOD::TrackParticle* ele_tp = electron->trackParticle(i);
      bool matches_trk1 = ele_tp->p4() == trk1.p4();
      bool matches_trk2 = ele_tp->p4() == trk2.p4();

      if (matches_trk1) MatchedElectrons(trk1).push_back(index_selected);
      if (matches_trk2) MatchedElectrons(trk2).push_back(index_selected);

      matches = matches || matches_trk1 || matches_trk2;
    }

    if (matches){
      selected.push_back(electron);
      index_selected++;
    }

  }

  if (MatchedElectrons(trk1).size() == 0 || MatchedElectrons(trk2).size() == 0)
    HG::fatal("Something went wrong - did not find the matching electrons that we should have.");

  return selected;
}

//______________________________________________________________________________
size_t HG::TrackHandler::nMatchedElectrons(const xAOD::TrackParticle& trk) const
{
  return MatchedElectrons(trk).size();
}

//______________________________________________________________________________
bool HG::TrackHandler::passIPCuts(xAOD::TrackParticle& trk)
{
    if (passIPCut.isAvailable(trk) && !passIPCut(trk)) { return false; }

    return true;
}

//______________________________________________________________________________
void HG::TrackHandler::decorateIPCut(xAOD::TrackParticle& trk)
{
  passIPCut(trk) = true;
  const xAOD::EventInfo *eventInfo = 0;

  if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
    fatal("Cannot access EventInfo");
  }

  double d0sig = xAOD::TrackingHelpers::d0significance(&trk, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
  
  d0significance(trk) = fabs(d0sig);

  if (fabs(d0sig) > m_d0BySigd0Cut) { passIPCut(trk) = false; }

  const xAOD::VertexContainer *vertexCont = 0;

  if (m_event->retrieve(vertexCont, "PrimaryVertices").isFailure()) { passIPCut(trk) = false; return; }

  const xAOD::Vertex *pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);

  if (pvx == nullptr) { passIPCut(trk) = false; return; }

  double z0 = trk.z0() + trk.vz() - pvx->z();
  z0 = z0 * sin(trk.theta());
  
  z0sinTheta(trk) = z0;

  if (fabs(z0) > m_z0Cut) { passIPCut(trk) = false; }
}
