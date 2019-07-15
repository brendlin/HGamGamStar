#ifndef HGamGamStar_MergedElectronID_v2
#define HGamGamStar_MergedElectronID_v2

#include "HGamGamStar/HggStarCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamGamStar/HggStarVariables.h"
#include "HGamGamStar/AngularPosition.h"
#include "HGamGamStar/TrackModel.h"
#include "HGamGamStar/TrackHandler.h"
#include "ElectronPhotonSelectorTools/ElectronSelectorHelpers.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"


namespace HG {

  class MergedElectronID_v2 {

  private:



    // Preselection cuts
    int m_PreselNPassBlayer;
    double m_PreselRhad;
    double m_mergedElePtCut;





  public:

    /// constructor
    MergedElectronID_v2();

    /// destructor
    virtual ~MergedElectronID_v2();

    virtual EL::StatusCode initialize(Config &config);

    bool passPIDCut(const xAOD::Electron &ele) const ;

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID_v2
