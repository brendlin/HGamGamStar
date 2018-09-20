#include "HGamGamStar/MergedElectronID.h"
//______________________________________________________________________________
HG::MergedElectronID::MergedElectronID()
{

}

//______________________________________________________________________________
HG::MergedElectronID::~MergedElectronID()
{

}

//______________________________________________________________________________

// EL::StatusCode HG::MergedElectronID::initialize(Config &config)
// {
//   //TODO: add initialization with cuts
//   return EL::StatusCode::SUCCESS;
// }

//______________________________________________________________________________
bool HG::MergedElectronID::passPIDCut(xAOD::Electron *ele,xAOD::TrackParticle *trk1,xAOD::TrackParticle *trk2){
  //electron ID calculations to be added here
  std::cout<<"MRGDELID call"<<std::endl;
  return true;
}
