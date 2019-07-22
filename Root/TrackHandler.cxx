#include "HGamGamStar/TrackHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include "xAODEgamma/ElectronxAODHelpers.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "PhotonVertexSelection/PhotonVertexHelpers.h"

#include "HGamAnalysisFramework/TruthUtils.h"

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

  //track selection
  // NOTE: if using track-based selection, these are applied on top of that, after
  // the indices are chosen.
  m_nSiMin     = config.getInt (m_name + ".Selection.nSiMin",7);
  m_nPixMin    = config.getInt (m_name + ".Selection.nPixMin",0);

  // do index-based track selection (see related function)
  // BE CAREFUL - this selection has hard-coded values to match the DAOD!
  m_doIndexBasedTrackSelection = config.getBool (m_name + ".Selection.IndexBasedTrackSelection",true);

  m_truth_nSiMin = config.getInt (m_name + ".TruthSelection.nSiMin",3);

  m_etaCut     = config.getNum(m_name + ".Selection.MaxAbsEta", 2.47);
  m_ptCut      = config.getNum(m_name + ".Selection.PtPreCutGeV", 0.5) * GeV;

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

    // If we already decorated everything, then do not decorate again!
    if (TrkAcc::TRT_PID_trans.isAvailable(*trk)) break;

    decorateIPCut(*trk);
    decorateAdditionalCuts(*trk);
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
bool HG::TrackHandler::passTrackPreselection(const xAOD::TrackParticle* trk,bool doTruthClassify) const
{

  if (doTruthClassify) {
    // (Truth-level) track preselection for truth classification

    int nSi = xAOD::EgammaHelpers::numberOfSiHits( trk );
    if( nSi < m_truth_nSiMin ) return false;
  }
  else {

    // (Reco-level) track preselection

    int nSiHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfSiliconHitsAndDeadSensors(trk);
    int nPixHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(trk);

    if ( std::abs(trk->eta()) > m_etaCut ) return false;
    if ( trk->pt() < m_ptCut ) return false;

    if (nSiHitsPlusDeadSensors  < m_nSiMin) return false;
    if (nPixHitsPlusDeadSensors < m_nPixMin) return false;
  }

    return true;
}


//______________________________________________________________________________
bool HG::TrackHandler::passIndexBasedTrackSelection(xAOD::Electron* ele,int i,
                                                    int& tmp_vtxTrkIndex1,int& tmp_vtxTrkIndex2) const
{
  // This function checks the indices of the objects as well as the tracking properties.
  // It is intended to match what is in the DOAD.
  // Start by giving it (by reference) tmp_vtxTrkIndex1 and tmp_vtxTrkIndex2 as -999, and it will
  // fill them in for you. Consistency between the DAOD result and this function is done later.

  const xAOD::TrackParticle* ele_tp = ele->trackParticle(i);

  int nSiHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfSiliconHitsAndDeadSensors(ele_tp);
  bool passBL = ElectronSelectorHelpers::passBLayerRequirement(ele_tp);

  bool passIndexRequirement = false;

  // Merged-ID track selection, intended to mirror the DAOD
  if (ele->caloCluster()->pt() > 5000)
  {
    // must pass preselection
    if (nSiHitsPlusDeadSensors >= 7 && passBL)
    {

      if (tmp_vtxTrkIndex1 < 0) {
        // First good track is assigned to the first "vertex" index
        tmp_vtxTrkIndex1 = i;
        passIndexRequirement = true;
      }
      else if (tmp_vtxTrkIndex2 < 0 &&
               ele_tp->charge() != ele->trackParticle(tmp_vtxTrkIndex1)->charge())
      {
        // If index2 has not been assigned yet, and this track is OS compared to index1, assign it.
        tmp_vtxTrkIndex2 = i;
        passIndexRequirement = true;
      }
    }
    // Two tracks were already assigned. Track is not eligible for Merged object.
  }

  // In any case, if it is the best-matched electron, we consider it for the merged ID
  if (i == 0) passIndexRequirement = true;

  if (!passIndexRequirement) return false;

  int nPixHitsPlusDeadSensors = ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(ele_tp);

  if ( std::abs(ele_tp->eta()) > m_etaCut ) return false;
  if ( ele_tp->pt() < m_ptCut ) return false;

  if (nSiHitsPlusDeadSensors  < m_nSiMin) return false;
  if (nPixHitsPlusDeadSensors < m_nPixMin) return false;

  return true;

  // For the resolved case, the 0th index is taken.
  if (i == 0) return true;
  return false;
}

//______________________________________________________________________________
xAOD::TrackParticleContainer HG::TrackHandler::findTracksFromElectrons(xAOD::TrackParticleContainer& container,
                                                                       xAOD::ElectronContainer& elecs,
                                                                       TrackElectronMap& trkEleMap)
{

  xAOD::TrackParticleContainer selected(SG::VIEW_ELEMENTS);

  int nGoodQuality = 0;

  // std::cout << "~~~~~~" << std::endl;
  for (auto electron : elecs) {

    // std::cout << Form("Electron pt: %.0f eta: %.3f has %lu tracks",
    //                   electron->pt(),
    //                   electron->caloCluster()->etaBE(2),
    //                   electron->nTrackParticles()) << std::endl;

    // For the new code wherein we select only two OS tracks
    int tmp_vtxTrkIndex1 = -999;
    int tmp_vtxTrkIndex2 = -999;

    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {

      const xAOD::TrackParticle* ele_tp = electron->trackParticle(i);

      if (!ele_tp) continue;

      // std::cout << Form("Electron tp; pt: %.0f eta: %.3f",ele_tp->pt(),ele_tp->eta()) << std::endl;

      if (m_doIndexBasedTrackSelection) {
        // Index-based track selection
        if (!passIndexBasedTrackSelection(electron,i,tmp_vtxTrkIndex1,tmp_vtxTrkIndex2)) continue;
      }
      else {
        // Non-index-based track selection
        if (!passTrackPreselection(ele_tp,false)) continue;
      }

      // Add track and electron to map
      MapHelpers::AddTrackElectronMapEntry(ele_tp,electron,trkEleMap);

      nGoodQuality++;

      bool found = false;

      // Check if it exists in the output TrackParticleContainer already
      for (auto track : selected) {
        xAOD::TrackParticle* track_p = (xAOD::TrackParticle*)HG::MapHelpers::getTheOriginalPointer(*track);
        if (track_p == ele_tp)
        {
          // std::cout << "Found a duplicate track!" << std::endl;
          found = true;
          break;
        }
      }

      // If it does not exist in the output TrackParticleContainer, add it.
      if (found) continue;

      xAOD::TrackParticle* container_tp = HG::MapHelpers::FindTrackParticle(&container,ele_tp);
      selected.push_back(container_tp);

    }

    // Electron: Check consistency between DAOD code and this code
    if (m_doIndexBasedTrackSelection) {
      if (tmp_vtxTrkIndex2 < 0) tmp_vtxTrkIndex1 = -999; // to maintain consistency with DAOD version

      if (EleAcc::vtxTrkIndex1.isAvailable(*electron)) {
        if (tmp_vtxTrkIndex1 != EleAcc::vtxTrkIndex1(*electron))
          HG::fatal(Form("Local vtxTrkIndex1 does not match DAOD! %d vs %d",
                         tmp_vtxTrkIndex1,EleAcc::vtxTrkIndex1(*electron)));
        if (tmp_vtxTrkIndex2 != EleAcc::vtxTrkIndex2(*electron))
          HG::fatal(Form("Local vtxTrkIndex2 does not match DAOD! %d vs %d",
                         tmp_vtxTrkIndex2,EleAcc::vtxTrkIndex2(*electron)));
      }
      else {
        EleAcc::vtxTrkIndex1(*electron) = tmp_vtxTrkIndex1;
        EleAcc::vtxTrkIndex2(*electron) = tmp_vtxTrkIndex2;
      }
    }

  }

  // std::cout << "Number of elecss in Total: " << elecs.size() << std::endl;
  // std::cout << "Number of OQ trk in Total: " << nGoodQuality << std::endl;
  // std::cout << "Final size of container: " << container.size() << std::endl;
  // std::cout << "Final size of selected: " << selected.size() << std::endl;

  return selected;

}

//______________________________________________________________________________
HG::TruthTrackMap HG::TrackHandler::MakeTruthTrackMapFromGSFContainer(xAOD::TrackParticleContainer& tracks)
{
  if ( !HG::isMC() ) fatal("Should not call MakeTruthTrackMap on data!");

  TruthTrackMap truthTrkMap;

  for (auto trkParticle : tracks )
  {

    int nSi = xAOD::EgammaHelpers::numberOfSiHits( trkParticle );
    if( nSi < m_truth_nSiMin ) continue;

    MapHelpers::AddTruthTrackMapEntry(trkParticle,truthTrkMap);
  }

  return truthTrkMap;
}

//______________________________________________________________________________
HG::TruthTrackMap HG::TrackHandler::MakeTruthTrackMapFromElectronContainer(const xAOD::ElectronContainer& elecs)
{
  if ( !HG::isMC() ) fatal("Should not call MakeTruthTrackMap on data!");

  TruthTrackMap truthTrkMap;

  for (auto electron : elecs) {
    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {
      const xAOD::TrackParticle* trkParticle = electron->trackParticle(i);

      int nSi = xAOD::EgammaHelpers::numberOfSiHits( trkParticle );
      if( nSi < m_truth_nSiMin ) continue;

      MapHelpers::AddTruthTrackMapEntry(trkParticle,truthTrkMap);
    }
  }

  return truthTrkMap;
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

  TrkAcc::MatchedElectrons(trk1).clear();
  TrkAcc::MatchedElectrons(trk2).clear();

  unsigned int index_selected = 0;

  for (auto electron : preSelElecs) {

    bool matches = false;

    for (unsigned int i=0; i<electron->nTrackParticles(); ++i) {
      const xAOD::TrackParticle* ele_tp = electron->trackParticle(i);
      bool matches_trk1 = ele_tp->p4() == trk1.p4();
      bool matches_trk2 = ele_tp->p4() == trk2.p4();

      if (matches_trk1) TrkAcc::MatchedElectrons(trk1).push_back(index_selected);
      if (matches_trk2) TrkAcc::MatchedElectrons(trk2).push_back(index_selected);

      matches = matches || matches_trk1 || matches_trk2;
    }

    if (matches){
      selected.push_back(electron);
      index_selected++;
    }

  }

  if (TrkAcc::MatchedElectrons(trk1).size() == 0 || TrkAcc::MatchedElectrons(trk2).size() == 0)
    HG::fatal("Something went wrong - did not find the matching electrons that we should have.");

  return selected;
}

//______________________________________________________________________________
size_t HG::TrackHandler::nMatchedElectrons(const xAOD::TrackParticle& trk) const
{
  return TrkAcc::MatchedElectrons(trk).size();
}

//______________________________________________________________________________
bool HG::TrackHandler::passIPCuts(xAOD::TrackParticle& trk)
{
  if (TrkAcc::passIPCut.isAvailable(trk) && !TrkAcc::passIPCut(trk)) { return false; }

    return true;
}

//______________________________________________________________________________
float HG::TrackHandler::calculateIPSig(const xAOD::TrackParticle& trk) const
{
  const xAOD::EventInfo *eventInfo = 0;

  if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
    fatal("Cannot access EventInfo");
  }

  float d0sig = xAOD::TrackingHelpers::d0significance(&trk, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
  return d0sig;
}

//______________________________________________________________________________
void HG::TrackHandler::decorateIPCut(xAOD::TrackParticle& trk)
{
  TrkAcc::passIPCut(trk) = true;

  float d0sig  = calculateIPSig(trk);
   
  TrkAcc::d0significance(trk) = d0sig;

  if (fabs(d0sig) > m_d0BySigd0Cut) { TrkAcc::passIPCut(trk) = false; }

  const xAOD::VertexContainer *vertexCont = 0;

  if (m_event->retrieve(vertexCont, "PrimaryVertices").isFailure()) { TrkAcc::passIPCut(trk) = false; return; }

  const xAOD::Vertex *pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);

  if (pvx == nullptr) { TrkAcc::passIPCut(trk) = false; return; }

  float z0 = trk.z0() + trk.vz() - pvx->z();
  TrkAcc::z0pv(trk) = z0; // delta z0 with respect to pvx

  float z0sinTheta = z0 * sin(trk.theta());
  TrkAcc::z0sinTheta(trk) = z0sinTheta;

  if (fabs(z0sinTheta) > m_z0Cut) { TrkAcc::passIPCut(trk) = false; }
}

float HG::TrackHandler::calculateTRT_PID(const xAOD::TrackParticle& trk) const
{
  float t_TRT_PID(0.0);
  trk.summaryValue(t_TRT_PID, xAOD::eProbabilityHT);
  const double tau = 15.0;
  if (t_TRT_PID >= 1.0) t_TRT_PID = 1.0 - 1.0e-15;
  if (t_TRT_PID < 1.0e-30) t_TRT_PID = 1.0e-30;
  return  - log(1.0/t_TRT_PID - 1.0) / tau;
}

void HG::TrackHandler::decorateAdditionalCuts(xAOD::TrackParticle& trk)
{
  // Add reco-level decorators
  TrkAcc::TRT_PID_trans(trk) = calculateTRT_PID(trk);
  TrkAcc::passBLayerRequirement(trk) = ElectronSelectorHelpers::passBLayerRequirement(&trk);
  TrkAcc::pt(trk) = trk.pt();
  TrkAcc::p(trk) = 1./(trk.qOverP());

  // Add default decorators
  TrkAcc::mergedTrackParticleIndex(trk) = -1;

  TrkAcc::nPixHitsAndDeadSens(trk) = ElectronSelectorHelpers::numberOfPixelHitsAndDeadSensors(&trk);
  TrkAcc::nSCTHitsAndDeadSens(trk) = ElectronSelectorHelpers::numberOfSCTHitsAndDeadSensors(&trk);

  uint8_t arg;
  TrkAcc::SharedIBL(trk) = trk.summaryValue(arg, xAOD::numberOfInnermostPixelLayerSharedHits      ) ? (int)arg : -999;
  TrkAcc::SharedBL(trk)  = trk.summaryValue(arg, xAOD::numberOfNextToInnermostPixelLayerSharedHits) ? (int)arg : -999;
  TrkAcc::SplitIBL(trk)  = trk.summaryValue(arg, xAOD::numberOfInnermostPixelLayerSplitHits       ) ? (int)arg : -999;
  TrkAcc::SplitBL(trk)   = trk.summaryValue(arg, xAOD::numberOfNextToInnermostPixelLayerSplitHits ) ? (int)arg : -999;

  TrkAcc::NIBL(trk) = trk.summaryValue(arg, xAOD::numberOfInnermostPixelLayerHits ) ? (int)arg : -999;
  // If NIBL is 0 and expectInnermostPixelLayerHit variable is defined,
  // fill NIBL with -(expectIBLHit) (0 = do not expect a hit or -1 = expected a hit)
  if (TrkAcc::NIBL(trk) == 0 && trk.summaryValue(arg, xAOD::expectInnermostPixelLayerHit)) {
    TrkAcc::NIBL(trk) = -(int)arg;
  }

  TrkAcc::NBL(trk) = trk.summaryValue(arg, xAOD::numberOfNextToInnermostPixelLayerHits) ? (int)arg : -999;
  if ( TrkAcc::NBL(trk) == 0 && trk.summaryValue(arg, xAOD::expectNextToInnermostPixelLayerHit)) {
    TrkAcc::NBL(trk) = -(int)arg;
  }

  // Make passTTVA branch (default is -1)
  TrkAcc::PassTTVA(trk) = 0;

  // Decorate MC particles with some truth information:
  if (HG::isMC()) {
    const xAOD::TruthParticle* truthPart = xAOD::TruthHelpers::getTruthParticle(trk);
    TrkAcc::isTrueHiggsElectron(trk) = truthPart && HG::isFromHiggs(truthPart) && HG::isGoodTruthElectron(truthPart);

    TrkAcc::PdgID(trk)   = truthPart ? truthPart->pdgId()   : -999;
    TrkAcc::Barcode(trk) = truthPart ? truthPart->barcode() : -999;
    TrkAcc::TruthE(trk)  = truthPart ? truthPart->p4().E()  : -999;
  }

}
