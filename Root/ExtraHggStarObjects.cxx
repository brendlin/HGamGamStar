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

//______________________________________________________________________________
void HG::ExtraHggStarObjects::clearContainers()
{
  m_tracks.clear();
  m_tracksAvail = false;
}
