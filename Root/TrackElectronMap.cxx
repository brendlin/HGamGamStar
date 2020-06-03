#include "HGamGamStar/TrackElectronMap.h"
#include "xAODEgamma/ElectronxAODHelpers.h"
#include "xAODTruth/xAODTruthHelpers.h"

#include "xAODBase/IParticleHelpers.h"
#include "HGamAnalysisFramework/HGamCommon.h"

#include <iostream>

void HG::MapHelpers::DumpMapInfo(  const TrackElectronMap& trkEleMap )
{
  std::cout << "Dumping map : " << std::endl;
  for(auto pair :  trkEleMap )
    std::cout << pair.first << "  " <<  pair.first->pt() << "  " <<  pair.first->eta() << "  "<<  pair.first->phi() << " El:  "<< pair.second << "  " <<  pair.second->pt() << "  " <<  pair.second->eta() << "  "<<  pair.second->phi() << std::endl;

  std::cout << "End dumping map  "<< std::endl;

}
//______________________________________________________________________________
float HG::MapHelpers::getTruthMatchProbability(const xAOD::TrackParticle* trackParticle)
{
  float truthProb = 0.;
  if (trackParticle->isAvailable<float>("truthMatchProbability")) {
    truthProb = trackParticle->auxdata<float>("truthMatchProbability");
  }
  return truthProb;

}

/////______________________________________________________________________________
const xAOD::TrackParticle* HG::MapHelpers::getTrackMatchingTruth( const xAOD::TruthParticle* truth,const TruthTrackMap& truthTrkMap)
{
  // Find the original pointer (if any) -- in order to work with deep copies
  const xAOD::TruthParticle* origPointer = (xAOD::TruthParticle*)getTheOriginalPointer(*truth);

  auto truthPair =  truthTrkMap.find( origPointer );
  if( truthPair != truthTrkMap.end()){
    return truthPair->second ;
  }

  //std::cout << "Lepton was not reconstructed" << std::endl;
  return nullptr;
}

//______________________________________________________________________________
std::vector<const xAOD::Electron*> HG::MapHelpers::getElectronsMatchingTrack( const xAOD::TrackParticle* track, const TrackElectronMap& trkEleMap )
{
  // Find the original pointer (if any) -- in order to work with deep copies
  const xAOD::TrackParticle* origPointer = (xAOD::TrackParticle*)getTheOriginalPointer(*track);

  std::vector<const xAOD::Electron*> electrons;
  auto Trks_Electrons = trkEleMap.equal_range( origPointer );
  for( auto mapIt = Trks_Electrons.first; mapIt != Trks_Electrons.second; ++mapIt){
    electrons.push_back( mapIt->second );
  }
  return electrons;
}

int HG::MapHelpers::getMatchingTrackIndex(const xAOD::Electron* electron, const xAOD::TrackParticle* track)
{
  // Find the original pointer (if any) -- in order to work with deep copies
  const xAOD::TrackParticle* origTrackPointer = (xAOD::TrackParticle*)getTheOriginalPointer(*track);

  int index = -1;
  // Loop over all tracks in the electron
  for( unsigned int trk_i(0); trk_i < electron->nTrackParticles(); ++trk_i){
    // Check to see if the track in the electron is the track matching to the truth
    if( electron->trackParticle(trk_i) == origTrackPointer )
    {
      //Save the index
      index = trk_i;
      break;
    }
  }
  return index;
}

//______________________________________________________________________________
std::vector<int> HG::MapHelpers::getMatchingTrackIndices(std::vector<const xAOD::Electron*>& electrons,
                                                                      const xAOD::TrackParticle* track )
{
  std::vector<int> indices;
  for( const auto& electron : electrons){
    indices.push_back( getMatchingTrackIndex(electron, track) );
  }
  return indices;
}

//______________________________________________________________________________
void HG::MapHelpers::AddTrackElectronMapEntry(const xAOD::TrackParticle* trkParticle,
                                              const xAOD::Electron* electron,
                                              TrackElectronMap& trkEleMap)
{
  const xAOD::TrackParticle* trkParticle_p = (xAOD::TrackParticle*)getTheOriginalPointer(*trkParticle);
  const xAOD::Electron* electron_p = (xAOD::Electron*)getTheOriginalPointer(*electron);
  trkEleMap.insert( std::pair<const xAOD::TrackParticle*,const xAOD::Electron*>( trkParticle_p, electron_p) );

  return;
}

//______________________________________________________________________________
void HG::MapHelpers::RemoveTrackElectronMapElectron(const xAOD::Electron* electron,
                                                    TrackElectronMap& trkEleMap)
{
  const xAOD::Electron* electron_p = (xAOD::Electron*)getTheOriginalPointer(*electron);

  for(auto it = trkEleMap.begin(); it != trkEleMap.end(); ) {
    if(it->second == electron_p)
      it = trkEleMap.erase(it);
    else
      ++it;
  }

  return;
}

//______________________________________________________________________________
void HG::MapHelpers::AddTruthTrackMapEntry(const xAOD::TrackParticle* trkParticle,TruthTrackMap& truthTrkMap) {
  if ( !HG::isMC() ) fatal("Should not call MakeTruthTrackMap on data!");

  const xAOD::TrackParticle* trkParticle_p = (xAOD::TrackParticle*)getTheOriginalPointer(*trkParticle);

  const xAOD::TruthParticle* truthPart = xAOD::TruthHelpers::getTruthParticle(*trkParticle_p);
  if(!truthPart) return;

  const xAOD::TruthParticle* truthPart_p = (xAOD::TruthParticle*)getTheOriginalPointer(*truthPart);

  // Add truth and track to map
  // If truth particle is already in the map choose the track with the higher match probability
  auto truthPair =  truthTrkMap.find( truthPart_p );
  if( truthPair != truthTrkMap.end()){
    if( getTruthMatchProbability(trkParticle_p) > getTruthMatchProbability(truthPair->second) )
      truthTrkMap[truthPart_p] = trkParticle_p;
  } else {
    truthTrkMap.insert( std::pair<const xAOD::TruthParticle*, const xAOD::TrackParticle*>(truthPart_p,trkParticle_p) );
  }
  return;
}

//______________________________________________________________________________
const xAOD::IParticle* HG::MapHelpers::getTheOriginalPointer(const xAOD::IParticle& part)
{
  // Find the original pointer (if any) -- in order to work with deep copies
  const xAOD::IParticle* compPointer = &part;
  if (xAOD::getOriginalObject(part)) {
    compPointer = xAOD::getOriginalObject(part);
  }
  return compPointer;
}

//______________________________________________________________________________
xAOD::Electron* HG::MapHelpers::FindElectron(xAOD::ElectronContainer* cont,
                                             const xAOD::Electron* toFind)
{
  // Find a container electron corresponding to the specified electron pointer.

  if (!cont) HG::fatal("FindElectron: Passed the function a null container.");

  for (auto el : *cont)
  {
    xAOD::Electron* el_p = (xAOD::Electron*)getTheOriginalPointer(*el);
    if (el_p == toFind) return el;
  }

  HG::fatal("FindElectron: Could not find a pointer in the container.");
  return nullptr;
}

//______________________________________________________________________________
xAOD::TrackParticle* HG::MapHelpers::FindTrackParticle(xAOD::TrackParticleContainer* cont,
                                                       const xAOD::TrackParticle* toFind)
{
  // Find a container TrackParticle corresponding to the specified TrackParticle pointer.

  if (!cont) HG::fatal("FindTrackParticle: Passed the function a null container.");

  for (auto tp : *cont)
  {
    xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*tp);
    if (tp_p == toFind) return tp;
  }

  HG::fatal("FindTrackParticle: Could not find a pointer in the container.");
  return nullptr;
}

//______________________________________________________________________________
const xAOD::TrackParticle* HG::MapHelpers::FindTrackParticle(const xAOD::TrackParticleContainer* cont,
                                                             const xAOD::TrackParticle* toFind)
{
  // Find a container TrackParticle corresponding to the specified TrackParticle pointer.

  if (!cont) HG::fatal("FindTrackParticle: Passed the function a null container.");

  for (auto tp : *cont)
  {
    const xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*tp);
    if (tp_p == toFind) return tp;
  }

  HG::fatal("FindTrackParticle: Could not find a (const) pointer in the container.");
  return nullptr;
}

//______________________________________________________________________________
int HG::MapHelpers::FindTrackParticleIndex(const xAOD::TrackParticleContainer* cont,
                                           const xAOD::TrackParticle* toFind)
{
  // Find a container TrackParticle corresponding to the specified TrackParticle pointer.

  if (!cont) HG::fatal("FindTrackParticle: Passed the function a null container.");

  for (unsigned int i=0;i<cont->size();i++)
  {
    const xAOD::TrackParticle* tp = (*cont)[i];
    const xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*tp);
    if (tp_p == toFind) return (int)i;
  }

  HG::fatal("FindTrackParticleIndex: Could not find a pointer in the container.");
  return -1;
}
