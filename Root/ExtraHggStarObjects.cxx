#include "HGamGamStar/ExtraHggStarObjects.h"

#include "xAODCore/AuxContainerBase.h"
#include "xAODCore/AuxInfoBase.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODTruth/TruthParticle.h"

#include "HGamAnalysisFramework/HgammaUtils.h"
#include "HGamGamStar/HggStarVariables.h"

//____________________________________________________________________________
namespace HG {
  ExtraHggStarObjects *ExtraHggStarObjects::m_ptr = nullptr;

  ExtraHggStarObjects *ExtraHggStarObjects::getInstance()
  {
    if (m_ptr == nullptr)
    { m_ptr = new ExtraHggStarObjects(); }

    return m_ptr;
  }

}

//____________________________________________________________________________
HG::ExtraHggStarObjects::ExtraHggStarObjects()
  : m_event(nullptr)
  , m_store(nullptr)
  , m_tracksAvail(false)
  , m_higgsLepsAvail(false)
  , m_higgsPhotonsAvail(false)
{

}

//____________________________________________________________________________
HG::ExtraHggStarObjects::~ExtraHggStarObjects()
{

}

//______________________________________________________________________________
void HG::ExtraHggStarObjects::setElectronTrackContainer(const xAOD::IParticleContainer *tracks)
{
  if (tracks) {
    m_tracks = *tracks;
    m_tracksAvail = true;
  }
  return;
}

//____________________________________________________________________________
const xAOD::IParticleContainer *HG::ExtraHggStarObjects::getElectronTracks(bool truth) const
{
  if (truth) {
    throw std::runtime_error("ExtraHggStarObjects::getElectronTrack should not be attempted in truth!");
  }

  if (!m_tracksAvail)
  { throw std::runtime_error("Track container requested but not set in ExtraHggStarObjects, throwing exception"); }

  return &m_tracks;
}

//______________________________________________________________________________
void HG::ExtraHggStarObjects::setTruthHiggsDecayProducts(const xAOD::TruthParticleContainer* all_particles)
{

  // // We are going to make a "deep" copy of these, in order to save them in the store:
  xAOD::TruthParticleContainer* photonCopy = new xAOD::TruthParticleContainer();
  xAOD::AuxContainerBase* photonCopyAux = new xAOD::AuxContainerBase();
  photonCopy->setStore(photonCopyAux);

  xAOD::TruthParticleContainer* leptonsCopy = new xAOD::TruthParticleContainer();
  xAOD::AuxContainerBase* leptonsCopyAux = new xAOD::AuxContainerBase();
  leptonsCopy->setStore(leptonsCopyAux);

  TruthPtcls higgses = getFinalHiggsBosons(all_particles);
  if (higgses.size() > 0)
  {

    TruthPtcls decayProds = getHyyStarSignalDecayProducts(higgses[0]);
    TruthPtcls childphot = FilterDirectPhotons(decayProds);

    if (childphot.size() == 1)
    {
      // copy to output container
      xAOD::TruthParticle *element = new xAOD::TruthParticle();
      photonCopy->push_back(element);
      *element = *(childphot[0]);
      setOriginalObjectLink(*(childphot[0]), *element);

      if (photonCopy->size() != 1) HG::fatal("Could not find the photon in original container - this should not happen.");
    }

    TruthPtcls childleps = FilterLeptons(decayProds);

    if (childleps.size() == 2 && (childleps[0]->absPdgId() == childleps[1]->absPdgId()))
    {
      for (auto lep : childleps) {

        // Copy to output container
        xAOD::TruthParticle *element = new xAOD::TruthParticle();
        leptonsCopy->push_back(element);
        *element = *lep;
        setOriginalObjectLink(*lep, *element);

      }
      if (leptonsCopy->size() != 2) HG::fatal("Could not find the leptons in original container - this should not happen.");
    }
  }

  if (m_store->record(photonCopy,"HiggsDirectPhotons").isFailure())
  { fatal("Cannot store deep copy of HiggsDirectPhotons to TStore, exiting."); }

  if (m_store->record(photonCopyAux,"HiggsDirectPhotonsAux").isFailure())
  { fatal("Cannot store deep copy of HiggsDirectPhotonsAux to TStore, exiting."); }

  if (m_store->record(leptonsCopy,"HiggsDecayLeptons").isFailure())
  { fatal("Cannot store deep copy of HiggsDecayLeptons to TStore, exiting."); }

  if (m_store->record(leptonsCopyAux,"HiggsDecayLeptonsAux").isFailure())
  { fatal("Cannot store deep copy of HiggsDecayLeptonsAux to TStore, exiting."); }

  m_higgsPhotons = *photonCopy;
  m_higgsPhotonsAvail = true;

  m_higgsLeps = *leptonsCopy;
  m_higgsLepsAvail = true;
  return;
}

//____________________________________________________________________________
const xAOD::TruthParticleContainer *HG::ExtraHggStarObjects::getTruthHiggsLeptons() const
{
  if (!m_higgsLepsAvail)
  { throw std::runtime_error("Higgs leptons requested but not set in ExtraHggStarObjects, throwing exception"); }

  return &m_higgsLeps;
}

//____________________________________________________________________________
const xAOD::TruthParticleContainer *HG::ExtraHggStarObjects::getTruthHiggsPhotons() const
{
  if (!m_higgsPhotonsAvail)
  { throw std::runtime_error("Higgs photon requested but not set in ExtraHggStarObjects, throwing exception"); }

  return &m_higgsPhotons;
}

//____________________________________________________________________________
void HG::ExtraHggStarObjects::setEventAndStore(xAOD::TEvent *event, xAOD::TStore *store)
{
  m_event = event;
  m_store = store;
}

//______________________________________________________________________________
void HG::ExtraHggStarObjects::clearContainers()
{
  m_tracks.clear();
  m_tracksAvail = false;

  m_higgsLeps.clear();
  m_higgsLepsAvail = false;

  m_higgsPhotons.clear();
  m_higgsPhotonsAvail = false;
}
