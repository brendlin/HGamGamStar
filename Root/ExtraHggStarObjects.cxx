#include "HGamGamStar/ExtraHggStarObjects.h"

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
  : m_tracksAvail(false)
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
void HG::ExtraHggStarObjects::clearContainers()
{
  m_tracks.clear();
  m_tracksAvail = false;
}
