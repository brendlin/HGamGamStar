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
   /**Constructor requires an ElectronContainer and to be told if its MC or not */
   TrackElectronMap(const xAOD::ElectronContainer&, bool isMC);
   ~TrackElectronMap();

   /** Find track best matched to the truth particle*/
   const xAOD::TrackParticle* getTrackMatchingTruth( const xAOD::TruthParticle* ) const;

   /** Find all electrons that have loosely matched a certain track*/
   std::vector<const xAOD::Electron*> getElectronsMatchingTrack( const xAOD::TrackParticle* ) const;

   /** Get the index of the track associated to the electron */
   int getMatchingTrackIndex(const xAOD::Electron*, const xAOD::TrackParticle* ) const;

   /** Get track index for  all electrons for a particular track */
   std::vector<int> getMatchingTrackIndex(std::vector<const xAOD::Electron*>&, const xAOD::TrackParticle* ) const;

   /** Get the truth match probability for a particular track particle*/
   float getTruthMatchProbability(const xAOD::TrackParticle* trackParticle) const;

  private:
   const bool m_isMC;
   std::map<const xAOD::TruthParticle*, const xAOD::TrackParticle*> m_truthTrackMap;
   std::multimap<const xAOD::TrackParticle*, const xAOD::Electron*> m_trackElectronMap;

};

}// End naespace
#endif
