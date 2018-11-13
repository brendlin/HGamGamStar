#include "HGamGamStar/TrackElectronMap.h"
#include "xAODEgamma/ElectronxAODHelpers.h"
#include "xAODTruth/xAODTruthHelpers.h"


////______________________________________________________________________________
HG::TrackElectronMap::TrackElectronMap( const xAOD::ElectronContainer& electrons, bool isMC )
: m_isMC( isMC )
{
  //Loop over all electrons and fill maps
  for( const xAOD::Electron* electron: electrons ){
    if(!electron) continue;
    for( unsigned int trk_i(0); trk_i < electron->nTrackParticles(); ++trk_i){
      auto trkParticle = electron->trackParticle(trk_i);
      if(!trkParticle)
        continue;
      // Ignore TRT only tracks
      int nSi  =  xAOD::EgammaHelpers::numberOfSiHits( trkParticle );
      if( nSi < 3)
        continue;

      //Get Truth
      const xAOD::TruthParticle* truthPart = nullptr;

      if(m_isMC)
      {
        truthPart = xAOD::TruthHelpers::getTruthParticle(*trkParticle);
        if(!truthPart)
          continue;
      }

      // Add track and electron to map
      m_trackElectronMap.insert( std::pair<const xAOD::TrackParticle*,const xAOD::Electron*>( trkParticle, electron) );

      if(!m_isMC)
        continue;
      // Add truth and track to map
      // If truth particle is already in the map choose the track with the higher match probability
      auto truthPair =  m_truthTrackMap.find( truthPart );
      if( truthPair != m_truthTrackMap.end()){
        if( getTruthMatchProbability(trkParticle) >  getTruthMatchProbability(truthPair->second) )
          m_truthTrackMap[truthPart] = trkParticle;
      } else {
        m_truthTrackMap.insert( std::pair<const xAOD::TruthParticle*, const xAOD::TrackParticle*>(truthPart,trkParticle) );
      }
    }
  }
}

////______________________________________________________________________________
HG::TrackElectronMap::~TrackElectronMap()
{

}
/////______________________________________________________________________________


float HG::TrackElectronMap::getTruthMatchProbability(const xAOD::TrackParticle* trackParticle) const
{
  float truthProb = 0.;
  if (trackParticle->isAvailable<float>("truthMatchProbability")) {
    truthProb = trackParticle->auxdata<float>("truthMatchProbability");
  }
  return truthProb;

}


const xAOD::TrackParticle* HG::TrackElectronMap::getTrackMatchingTruth( const xAOD::TruthParticle* truth) const
{
  auto truthPair =  m_truthTrackMap.find( truth );
  if( truthPair != m_truthTrackMap.end()){
    return truthPair->second ;
  } else {
    //std::cout << "Lepton was not reconstructed" << std::endl;
    return nullptr;
  }
}

std::vector<const xAOD::Electron*> HG::TrackElectronMap::getElectronsMatchingTrack( const xAOD::TrackParticle* track ) const
{
  std::vector<const xAOD::Electron*> electrons;
  auto Trks_Electrons = m_trackElectronMap.equal_range( track );
  for( auto mapIt = Trks_Electrons.first; mapIt != Trks_Electrons.second; ++mapIt){
    electrons.push_back( mapIt->second );
  }
  return electrons;
}

int HG::TrackElectronMap::getMatchingTrackIndex(const xAOD::Electron* electron, const xAOD::TrackParticle* track) const
{
  int index = -1;
  // Loop over all tracks in the electron
  for( unsigned int trk_i(0); trk_i < electron->nTrackParticles(); ++trk_i){
    // Check to see if the track in the electron is the track matching to the truth
    if( electron->trackParticle(trk_i) == track )
    {
      //Save the index
      index = trk_i;
      break;
    }
  }
  return index;
}


std::vector<int> HG::TrackElectronMap::getMatchingTrackIndex(std::vector<const xAOD::Electron*>& electrons,
                                                             const xAOD::TrackParticle* track ) const
{
  std::vector<int> indices;
  for( const auto& electron : electrons){
    indices.push_back( getMatchingTrackIndex(electron, track) );
  }
  return indices;
}
