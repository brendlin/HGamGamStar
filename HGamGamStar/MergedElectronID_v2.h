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
#include <string>
#include <TH2.h>


namespace HG {

  class MergedElectronID_v2 {

  private:

    // Preselection cuts
    // These are applied via the original MergedElectronID passPreselection function!
    /* int m_PreselNPassBlayer; */
    /* float m_PreselRhad; */
    /* flaot m_mergedElePtCut; */

    bool m_isMC;
    std::string m_deltaEtaFileName;
    std::string m_deltaEtaHistName;
    TH2* m_sdetaCorr;

    float correctDeta1(const xAOD::Electron &ele, bool isMC) const;
  

  public:

    /// constructor
    MergedElectronID_v2();

    /// destructor
    virtual ~MergedElectronID_v2();

    virtual EL::StatusCode initialize(Config &config);

    bool passPIDCut(const xAOD::Electron &ele, bool isMC) const ;

  };

} // namespace HG

#endif // HGamGamStar_MergedElectronID_v2
