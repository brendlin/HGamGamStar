#ifndef TrackElectronMap_h
#define TrackElectronMap_h

#include <map>
#include <vector>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODEgamma/ElectronContainer.h>

namespace HG
{

  typedef std::map<const xAOD::TruthParticle*, const xAOD::TrackParticle*> TruthTrackMap;
  typedef std::multimap<const xAOD::TrackParticle*, const xAOD::Electron*> TrackElectronMap;

  namespace MapHelpers
  {

   /** Find track best matched to the truth particle*/
   const xAOD::TrackParticle* getTrackMatchingTruth( const xAOD::TruthParticle* truth, const TruthTrackMap& trkTruthMap );

   /** Find all electrons that have loosely matched a certain track*/
   std::vector<const xAOD::Electron*> getElectronsMatchingTrack( const xAOD::TrackParticle* track, const TrackElectronMap& trkEleMap );

   /** Get the index of the track associated to the electron */
   int getMatchingTrackIndex(const xAOD::Electron*, const xAOD::TrackParticle* );

   /** Get track index for  all electrons for a particular track */
   std::vector<int> getMatchingTrackIndices(std::vector<const xAOD::Electron*>&, const xAOD::TrackParticle* );

   /** Get the truth match probability for a particular track particle*/
   float getTruthMatchProbability(const xAOD::TrackParticle* trackParticle);

   const xAOD::IParticle* getTheOriginalPointer(const xAOD::IParticle& part);
   void AddTruthTrackMapEntry(const xAOD::TrackParticle* trkParticle, TruthTrackMap& trkTruthMap);
   void AddTrackElectronMapEntry(const xAOD::TrackParticle* trkParticle, const xAOD::Electron* electron, TrackElectronMap& trkEleMap);


   /** Find the electron in a particular container */
   xAOD::Electron* FindElectron(xAOD::ElectronContainer* cont, const xAOD::Electron* toFind);

  }// End namespace

}// End namespace
#endif
