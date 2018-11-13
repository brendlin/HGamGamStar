#ifndef TrackElectronMap_h
#define TrackElectronMap_h

#include <map>
#include <vector>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODEgamma/ElectronContainer.h>

namespace HG
{

class TrackElectronMap
{

  public:
   TrackElectronMap(const xAOD::ElectronContainer&, bool isMC);
   ~TrackElectronMap();


   const xAOD::TrackParticle* getTrackMatchingTruth( const xAOD::TruthParticle* ) const;

   std::vector<const xAOD::Electron*> getElectronsMatchingTrack( const xAOD::TrackParticle* ) const;

   int getMatchingTrackIndex(const xAOD::Electron*, const xAOD::TrackParticle* ) const;

   std::vector<int> getMatchingTrackIndex(std::vector<const xAOD::Electron*>&, const xAOD::TrackParticle* ) const;


   float getTruthMatchProbability(const xAOD::TrackParticle* trackParticle) const;

  private:
   const bool m_isMC;
   std::map<const xAOD::TruthParticle*, const xAOD::TrackParticle*> m_truthTrackMap;
   std::multimap<const xAOD::TrackParticle*, const xAOD::Electron*> m_trackElectronMap;

};

}// End naespace
#endif
