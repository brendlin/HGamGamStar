#ifndef HGamGamStar_MergedElectronID
#define HGamGamStar_MergedElectronID

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"


namespace HG {

  class MergedElectronID {

  private:

  public:

    /// constructor
    MergedElectronID();

    /// destructor
    virtual ~MergedElectronID();

    virtual EL::StatusCode initialize(Config &config);
    
    bool passPIDCut(xAOD::Electron &ele,xAOD::TrackParticle &trk1,xAOD::TrackParticle &trk2);

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID
